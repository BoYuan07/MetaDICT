#' @title Microbiome data integration method via shared dictionary learning.
#'
#' @description A method for microbiome data integration. This method is designed to 
#' remove batch effects and preserve biological variation while integrating heterogeneous datasets. 
#' MetaDICT can better avoid overcorrection when unobserved confounding variables are present.
#'
#' @details MetaDICT is a two-step approach. It initially estimates the batch effects by covariate balancing,
#' then refines the estimation via shared dictionary learning.
#'
#' @importFrom stats optim
#' @importFrom edgeR calcNormFactors
#' @importFrom stats varimax dist median glm predict
#'
#' @param count The integrated count table (taxa-by-sample matrix).  
#'   Should be provided as either a \code{matrix} or a \code{data.frame}.
#' @param meta The integrated meta table containing sample information  
#'   and batch IDs. The data must include a column named 'batch' containing all batch IDs. The row names of the meta should match the sample names in the count table. 
#' @param covariates The covariates used in data integration. Default is `"all"`.
#' @param tree The phylogenetic tree (optional if a distance matrix or taxonomy is provided).
#' @param taxonomy The taxonomy table (optional if a distance matrix or phylogenetic tree is provided). The row names of the taxonomy table should match the taxa names in the count table.
#' @param distance_matrix A \code{matrix} measuring the dissimilarity of taxa.  
#'   Default is \code{NULL}, in which case MetaDICT generates a distance matrix  
#'   based on phylogenetic and taxonomic information.
#' @param tax_level The taxonomic level of the count table.
#' @param customize_parameter A logical variable. Set to \code{TRUE} if the  
#'   \code{alpha} and \code{beta} parameters are customized. If \code{FALSE},  
#'   MetaDICT determines these parameters based on the number of covariates.
#' @param alpha A parameter controlling the rank of the final corrected count table.  
#'   A larger \code{alpha} leads to a lower-rank shared dictionary.
#' @param beta A parameter controlling the smoothness of the estimated measurement efficiency.  
#'   A larger \code{beta} results in more similar measurement efficiencies across taxa.
#' @param normalization The normalization method. Options are `"Upper quantile"`, `"RSim"` or `"TSS"`.  
#'   Set to \code{NULL} if normalization is not needed.
#' @param max_iter The maximum number of iterations for the optimization process. Default is \code{10000}.
#' @param imputation A logical variable. Whether to allow MetaDICT to perform imputation  
#'   based on dictionary learning results. Default is \code{FALSE}.
#' @param verbose A logical variable. Whether to generate verbose output. Default is \code{TRUE}.
#' @param optim_trace A logical variable. Whether to print optimization steps. Default is \code{FALSE}.
#'
#' @returns A \code{list} with the following components:
#'   \item{\code{count}}{ (\code{data.frame}) – The corrected count table.  
#'     Rows represent taxa, and columns represent samples.}
#'   \item{\code{D}}{ (\code{matrix}) – The estimated shared dictionary.}
#'   \item{\code{R}}{ (\code{matrix}) – The estimated sample representation.}
#'   \item{\code{w}}{ (\code{matrix}) – The estimated measurement efficiency.  
#'     Rows represent datasets, and columns represent taxa.}
#'   \item{\code{meta}}{ (\code{data.frame}) – The meta table used in the covariate balancing step.}
#'   \item{\code{dist_mat}}{ (\code{matrix}) – The distance matrix measuring taxa dissimilarity.}
#'
#' @examples 
#'  data(exampleData)
#'  O = exampleData$O
#'  meta = exampleData$meta
#'  dist_mat = exampleData$dist_mat
#'  metadict_res = MetaDICT(O, meta, distance_matrix = dist_mat)
#' 
#' @export
#' 
MetaDICT <- function(count, meta, covariates = "all", tree = NULL, taxonomy = NULL, distance_matrix = NULL,
tax_level = NULL, customize_parameter = FALSE, alpha = 0.1, beta = 0.01,  
normalization = "uq", max_iter = 10000, imputation = FALSE, verbose = TRUE, optim_trace = FALSE){

    # convert input to required format of MetaDICT
    metadict_input <- data_check(count = count, meta = meta, covariates = covariates, tree = tree, 
                        distance_matrix = distance_matrix, taxonomy = taxonomy, tax_level = tax_level, verbose = verbose)

    O.list <- metadict_input$count_list
    meta.list <- metadict_input$meta_list
    dist_mat <- metadict_input$dist_mat
    controls <- metadict_input$controls
    neighbor <- controls$neighbor
    sigma <- controls$sigma
    r <- controls$r
    m <- length(O.list)
    d <- nrow(dist_mat)
    gamma <- 1
    O <- do.call(cbind, O.list)
    sample_num <- sapply(O.list,ncol)

    if (!customize_parameter){
        alpha <- controls$alpha
        beta <- controls$beta
    }
    if (verbose){
        message(paste("Paremeters are set to be", "alpha = ", alpha, "beta = ", beta))
    }

    # normalization
    if (verbose){
        if (is.null(normalization)){
            message("Normalization is skipped.")
        }else if (! normalization %in% c("uq", "rsim")){
            stop(paste("Normalization method is not supported by MetaDICT!",
            "Please run normalization first then used the normalized counts as input while
            set Normalization = FALSE.", sep = "\n"))
        }else{
            message(paste("Normalization starts with method", normalization))  
        }
    }
    if(is.null(normalization)){
        O_norm <- O.list
    }else if(normalization == "uq"){
        O_norm <- lapply(O.list,function(x)uq(x)$P)
    }else if(normalization == "rsim"){
        O_norm <- lapply(O.list,function(x)rsim(x)$P)
    }else if(normalization == "tss"){
        O_norm <- lapply(O.list,function(x)tss(x))
    }
    message("Normalization completed. Beginning model fitting...") 
    
    scale <- max(unlist(O_norm))
    O_list_scaled <- lapply(O_norm,function(x)x/scale)

    # Laplacian matrix of sequencing graph
    adj_mat <- matrix(0,d,d)
    for(i in 1:d){
        idx <- order(dist_mat[i,],decreasing = F)[2:(neighbor+1)]
        adj_mat[i,idx] <- exp(-dist_mat[i,idx]/sigma)
        adj_mat[idx,i] <- exp(-dist_mat[i,idx]/sigma)
    }
    degree_mat <- diag(rowSums(adj_mat))
    L <- degree_mat-adj_mat

    # step 1: covariates balancing
    initial <- init_algorithm(m,d,O_list_scaled,r,meta.list)

    # step 2: shared dictionary learning
    ## prepare the input vector
    x0 <- convert_to_vec(m,initial$w,initial$D,initial$R)

    ## set the (0,1) constraint of measurement efficiency
    lower <- c(rep(0,m*d),rep(-Inf,d*r+r*sum(sample_num)))
    upper <- c(rep(1,m*d),rep(Inf,d*r+r*sum(sample_num)))

    if (optim_trace){
        trace = 3
    }else{
        trace = 0
    }

    optim.res <- optim(x0, fn = (function(x) target_func(x,O_list_scaled,alpha,beta,gamma,m,d,r,sample_num,L)), 
                        gr = (function(x) gradient_func(x,O_list_scaled,alpha,beta,gamma,m,d,r,sample_num,L)), 
                        method = "L-BFGS-B", lower = lower, upper = upper, 
                        control = list(maxit = max_iter, trace = trace))
    

    x <- optim.res$par
    para <- convert_from_vec(x,m,r,d,sample_num)
    w_list <- para$w
    D <- para$D
    R_list <- para$R
    
    # apply varimax to make D interpretable
    var_res <- varimax(D, normalize = FALSE)
    D_rot <- var_res$loadings
    T_mat <- var_res$rotmat

    X_list <- list()
    R_list_rot <- list()
    for(i in 1:m){
        R_list_rot[[i]] <- t(T_mat)%*%R_list[[i]]*scale
        X_list[[i]] <- D_rot%*%R_list_rot[[i]]
    }

    # corrected count table
    X <- do.call(cbind,X_list)
    X[X<0] <- 0

    if (!imputation){
        X[O==0] <- 0
    }

    effective_r <- effective_rank(D)

    error_each <- sapply(1:length(O_norm), function(i) norm(O_norm[[i]]-diag(w_list[i,])%*%X_list[[i]], "F")^2/norm(O_norm[[i]],"F")^2)


    if(optim.res$convergence==0){
        message("Successful convergence.")
    }
    if(optim.res$convergence==1){
        message("The iteration limit max_iter has been reached. Please consider increasing max_iter.")
    }
    message(paste("Maximum relative error:", max(error_each), "\n", "Effective rank of D", effective_r))

    if (verbose){
        message("Finished.")
    }
    meta_filtered <- do.call(rbind,meta.list)
    res.metadict <- as.data.frame(X)
    colnames(res.metadict) <- colnames(count)
    rownames(res.metadict) <- rownames(count)
    return(list(count = res.metadict, D = D_rot, R = R_list_rot, w = w_list, meta = meta_filtered, dist_mat = dist_mat))
}

#=================covariate balancing step===================#
init_algorithm <- function(m,d,O_list,r,meta.list){
    # initialize the measurement efficiency matrix
    w_list <- matrix(0, nrow = m, ncol = d)

    # the first dataset is used as the reference
    w_list[1,] <- 1

    for(i in 2:m){
        O1 <- O_list[[1]]
        O2 <- O_list[[i]]
        batchnum <- as.factor(c(rep(1,ncol(O1)),rep(i,ncol(O2))))

        # add library size in the meta table
        if(is.null(meta.list)){
            meta <- data.frame("lib" = c(colSums(O1),colSums(O2)),"batch" = batchnum)
        }else{
            meta <- rbind(meta.list[[1]],meta.list[[i]])
            meta$batch <- batchnum
            meta$lib <- c(colSums(O1),colSums(O2))
        }

        # propensity score estimation
        mylogit <- glm(batch ~ ., data = meta, family = "binomial")
        meta$psvalue <- predict(mylogit, type="response")

        # weighting
        meta$weight <- ifelse(meta$batch==i,1/meta$psvalue,1/(1-meta$psvalue))
        O_adj_1 <- t(t(O1)*meta$weight[which(meta$batch==1)])
        O_adj_2 <- t(t(O2)*meta$weight[which(meta$batch==i)])

        # estimated measurement efficiency for dataset i
        w_list[i,] <- (rowMeans(O_adj_2)+1e-6)/(rowMeans(O_adj_1)+1e-6)
    }

    # rescale the matrix
    w_list <- t(t(w_list)/colSums(w_list))

    # adjust the count list
    O_list_adj <- lapply(1:m,function(i)diag(1/w_list[i,])%*%O_list[[i]])

    # estimate D and R
    O <- do.call(cbind,O_list_adj)
    svd.res <- svd(O)
    D <- svd.res$u[,1:r]
    R_list <- lapply(1:m,function(i)t(D)%*%O_list_adj[[i]])

    # output
    return(list(w = w_list, D = D, R = R_list))
}

#=================optimization gradient and target functions===================#

target_func <- function(x, O_list, alpha, beta, gamma, m, d, r, sample_num, L){
  res <- convert_from_vec(x, m, r, d, sample_num)
  w_list <- res$w
  D <- res$D
  R_list <- res$R
  target <- 0
  for(i in 1:m){
    w <- w_list[i,]
    W <- diag(w)
    O_diff <- O_list[[i]]-W%*%D%*%R_list[[i]]
    W_diff <- t(w)%*%L%*%w
    target <- target+gamma*norm(O_diff,"F")**2+beta*W_diff/(d*d)+norm(R_list[[i]],"F")**2*alpha/(2*r*sample_num[i])
    
  }
  target <- target+norm(D,"F")**2*alpha/(2*d*r)
  return(target)
}

gradient_func = function(x, O_list, alpha, beta, gamma, m, d, r, sample_num, L){
  res <- convert_from_vec(x,m,r,d,sample_num)
  w_list <- res$w
  D <- res$D
  R_list <- res$R
  gradw <- matrix(0,m,d)
  gradD <- matrix(0,d,r)
  gradR <- list()
  for(i in 1:m){
    w <- w_list[i,]
    W <- diag(w)
    O_diff <- O_list[[i]]-W%*%D%*%R_list[[i]]
    W_diff <- t(w)%*%L%*%w
    gradD <- gradD+gamma*(-2*W%*%O_diff%*%t(R_list[[i]]))
    gradR[[i]] <- gamma*(-2)*t(D)%*%W%*%O_diff+R_list[[i]]*alpha/(r*sample_num[i])
    B <- gamma*(-2)*(O_diff)%*%t(R_list[[i]])%*%t(D)
    gradw[i,] <- diag(B)+2*t(w)%*%L*beta/(d*d)
  }
  gradD <- gradD+alpha*D/(d*r)
  grad <- convert_to_vec(m,gradw,gradD,gradR)
  return(grad)
}

#=================helper function of optimization===================#
convert_from_vec <- function(x, m, r, d, sample_num){
  w.vec <- x[1:(d*m)]
  D.vec <- x[(d*m+1):(d*m+d*r)]
  R_list <- list()
  n0 <- d*m+d*r
  for(i in 1:length(sample_num)){
    num <- sample_num[[i]]
    R_list[[i]] <- matrix(x[(n0+1):(n0+num*r)],r,num)
    n0 <- n0+num*r
  }
  D <- matrix(D.vec,d,r)
  w_list <- matrix(w.vec,m,d)
  return(list(w = w_list, D = D, R = R_list))
}

convert_to_vec <- function(m, w_list, D, R_list){
  x <- c(c(w_list),c(D),unlist(R_list))
  return(x)
}

#====================Convergence Diagnosis===========================#

effective_rank <- function(A) {
  svd_vals <- svd(A)$d  # Compute singular values
  # Normalize singular values to get probabilities
  p <- svd_vals / sum(svd_vals)
  # Compute entropy
  entropy <- -sum(p * log(p), na.rm = TRUE)
  # Compute effective rank
  return(exp(entropy))
}

#====================Normalization helper===========================#

# upper quantile normalization
uq <- function(X){
    dds <- edgeR::calcNormFactors(as.matrix(X+1), method = "upperquartile")
    upperQ <- dds*colSums(X)
    upq.res <- scale(X,center=FALSE,scale=upperQ)
  return(list('P' = upq.res, 'sf' = upperQ))
}

# rsim
## calculate statistic
CStat <- function(X){
  d <- nrow(X)
  R <- X
  S1 <- apply(R,1,order)
  S1 <- S1 - colMeans(S1);
  S1 <- S1 / sqrt(colSums(S1^2));
  corr_s <- crossprod(S1)
  med <- as.data.frame(matrixStats::colMedians(corr_s))
  return(as.numeric(med[,1]))

}

# main function for rsim
rsim <- function(X,eta=0.1){
    d <- nrow(X)
    v <- CStat(X)
    for (gamma in seq(0.8, 0, by = -0.1)){
        I0.1 <- which(v>gamma)
        if (length(I0.1) < 0.1*d){
            break
        }
    }
    X0 <- X[I0.1,]
    v0 <- replicate(3,CStat(X0[sample(1:nrow(X0),0.5*nrow(X0)),]))
    w <- v[v>gamma]
    f1 <- sapply(w,function(x)mean(v>x))
    f0 <- sapply(w,function(x)mean(v0>x))
    pi_val <- sum(f1*f0)/sum(f0^2)
    vord <- order(v,decreasing = T)
    res <- sapply(1:length(vord),function(x)(1-pi_val*length(vord)*mean(v0>v[x])/(which(vord==x))))
    lowerx <- max(which(res[vord]<eta))
    ref <- vord[1:lowerx]
    tc.cn <- apply(X,2,function(x)sum(x[ref]))
    f.cn <- tc.cn/(mean(tc.cn))
    f.cn <- ifelse(f.cn==0,1,f.cn)
    cn.res <- scale(X,center=FALSE,scale=f.cn)
  return(list('P' = cn.res, 'I0' = ref, 'pi0'= pi, 'sf'=f.cn))
}

# TSS
tss <- function(X){
    return(t(t(X)/colSums(X)))
}

