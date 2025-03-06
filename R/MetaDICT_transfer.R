#' @title Batch correction for new datasets using existing dictionary.
#'
#' @description This function adds new studies to an integrated dataset using a pre-learned dictionary.  
#' The corrected data can be directly used with machine learning models trained on the  
#' previously integrated dataset, enabling seamless application without retraining.  
#' 
#' @details This function estimates measurement efficiency and debiased representations for new studies  
#' while keeping the dictionary unchanged.  
#'
#' @importFrom stats optim
#' @importFrom edgeR calcNormFactors
#'
#' @param newdata The integrated count table of new studies.  
#'   Rows represent taxa, and columns represent samples.  
#'   Should be provided as either a \code{matrix} or a \code{data.frame}.
#' @param newmeta The integrated meta table (\code{meta}) for the new studies,  
#'   containing sample information and batch IDs.
#' @param integrated_result The output list from a previous MetaDICT integration task.
#' @param customize_parameter A logical variable.  
#'   Set to \code{TRUE} if the \code{beta} parameter is customized.  
#'   If \code{FALSE}, MetaDICT determines \code{beta} based on the number of covariates.
#' @param beta A parameter controlling the smoothness of the estimated measurement efficiency.  
#'   A larger \code{beta} results in more similar measurement efficiency across taxa.
#' @param normalization The normalization method. Options are `"Upper quantile"`, `"RSim"` or `"TSS"`.  
#'   Set to \code{NULL} if normalization is not needed.
#'   This should be the same as in the previous integration task.
#' @param max_iter The maximum number of iterations for the optimization process.  
#'   Default is \code{10000}.
#' @param imputation A logical variable.  
#'   Whether to allow MetaDICT to perform imputation based on dictionary learning results.  
#'   Default is \code{FALSE}.
#' @param verbose A logical variable.  
#'   Whether to generate verbose output. Default is \code{TRUE}.
#' @param optim_trace A logical variable.  
#'   Whether to print optimization steps. Default is \code{FALSE}.
#' 
#' @returns A \code{list} with the following components:
#' \itemize{
#'   \item{\code{count}}{ (\code{data.frame}) – The corrected count table.  
#'     Rows represent taxa, and columns represent samples.}
#'   \item{\code{D}}{ (\code{matrix}) – The estimated shared dictionary.}
#'   \item{\code{R}}{ (\code{matrix}) – The estimated sample representation.}
#'   \item{\code{w}}{ (\code{matrix}) – The estimated measurement efficiency.  
#'     Rows represent datasets, and columns represent taxa.}
#'   \item{\code{meta}}{ (\code{data.frame}) – The meta table used in the covariate balancing step.}
#'   \item{\code{dist_mat}}{ (\code{matrix}) – The distance matrix measuring taxa dissimilarity.}
#' }
#' 
#' @export
#' 

metadict_add_new_data <- function(newdata, newmeta, integrated_result, customize_parameter = FALSE, beta = 0.01,  
normalization = "uq", max_iter = 10000, imputation = FALSE, verbose = TRUE, optim_trace = FALSE){

    D <- integrated_result$D 
    O_old <- integrated_result$count
    meta_old <- integrated_result$meta
    dist_mat <- integrated_result$dist_mat

    r <- ncol(D)

    covariates <- intersect(colnames(meta_old), colnames(newmeta))

   # convert input to required format of MetaDICT
    metadict_input <- data_check(count = newdata, meta = newmeta, covariates = covariates,  
                        distance_matrix = dist_mat, verbose = verbose)

    O.list <- metadict_input$count_list
    meta.list <- metadict_input$meta_list
    dist_mat <- metadict_input$dist_mat
    controls <- metadict_input$controls
    neighbor <- controls$neighbor
    sigma <- controls$sigma
    r <- ncol(D)
    m <- length(O.list)
    d <- nrow(dist_mat)
    gamma <- 1
    O <- do.call(cbind, O.list)
    sample_num <- sapply(O.list,ncol)

    if (!customize_parameter){
        alpha <- 0.001
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
    if(normalization == "uq"){
        O_norm <- lapply(O.list,function(x)uq(x)$P)
    }else if(normalization == "rsim"){
        O_norm <- lapply(O.list,function(x)rsim(x)$P)
    }else if(normalization == "tss"){
        O_norm <- lapply(O.list,function(x)tss(x))
    }else if(!normalization){
        O_norm <- O.list
    }

    O_all.list <- append(O_norm, list(O_old), after = 0)
    meta_all.list <- append(meta.list, list(meta_old), after =  0)


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
    initial <- init_algorithm_transfer(O_all.list, D, meta_all.list)
    w_list <- initial$w
    R_list <- initial$R

    sample_num <- sapply(O_norm,ncol)
    x0 <- convert_to_vec_transfer(m, w, R_list)

    lower <- c(rep(0,m*d),rep(-Inf,r*sum(sample_num)))
    upper <- c(rep(1,m*d),rep(Inf,r*sum(sample_num)))

     if (optim_trace){
        trace = 3
    }else{
        trace = 0
    }
    
    optim.res <- optim(x0, fn = (function(x) target_func_transfer(x,O_norm,D,alpha,beta,gamma,m,d,r,sample_num,L)), 
    gr = (function(x) gradient_func_transfer(x,O_norm,D,alpha,beta,gamma,m,d,r,sample_num,L)), 
    method = "L-BFGS-B", lower = lower, upper = upper, control = list(maxit = max_iter, trace = trace))

    x <- optim.res$par
    para <- convert_from_vec_transfer(x,m,r,d,sample_num)
    w_list <- para$w
    R_list <- para$R

    X_list <- list()
    for(i in 1:m){
        X_list[[i]] <- D%*%R_list[[i]]
    }

    # corrected count table
    X <- do.call(cbind,X_list)
    X[X<0] <- 0

    if (!imputation){
        X[O==0] <- 0
    }


    error_each <- sapply(1:length(O_norm), function(i) norm(O_norm[[i]]-diag(w_list[i,])%*%X_list[[i]], "F")^2/norm(O_norm[[i]],"F")^2)


    if(optim.res$convergence==0){
        message("Successful convergence.")
    }
    if(optim.res$convergence==1){
        message("The iteration limit max_iter has been reached. Please consider increasing max_iter.")
    }
    message(paste("Maximum relative error:", max(error_each), "\n"))

    if (verbose){
        message("Finished.")
    }
    meta_filtered <- do.call(rbind,meta.list)
    res.metadict <- as.data.frame(X)
    colnames(res.metadict) <- colnames(count)
    rownames(res.metadict) <- rownames(count)

    return(list(count = res.metadict, D = D, R = R_list, w = w_list, meta = meta_filtered, dist_mat = dist_mat))
}


#=================optimization gradient and target functions===================#

target_func_transfer <- function(x, O_list, D, alpha, beta, gamma, m, d, r, sample_num, L){
  res <- convert_from_vec_transfer(x,m,r,d,sample_num)
  w_list <- res$w
  R_list <- res$R
  target <- 0
  for(i in 1:m){
    w <- w_list[i,]
    W <- diag(w)
    O_diff <- O_list[[i]]-W%*%D%*%R_list[[i]]
    W_diff <- t(w)%*%L%*%w
    target <- target+gamma*norm(O_diff,"F")**2+beta*W_diff/(d*d)+norm(R_list[[i]],"F")**2*alpha/(2*r*sample_num[i])
  }
  return(target)
}

gradient_func_transfer <- function(x, O_list, D, alpha, beta, gamma, m, d, r, sample_num, L){
  res <- convert_from_vec_transfer(x,m,r,d,sample_num)
  w_list <- res$w
  R_list <- res$R
  gradw <- matrix(0,m,d)
  gradR <- list()
  for(i in 1:m){
    w <- w_list[i,]
    W <- diag(w)
    O_diff <- O_list[[i]]-W%*%D%*%R_list[[i]]
    W_diff <- t(w)%*%L%*%w
    gradR[[i]] <- gamma*(-2)*t(D)%*%W%*%O_diff+R_list[[i]]*alpha/(r*sample_num[i])
    B <- gamma*(-2)*(O_diff)%*%t(R_list[[i]])%*%t(D)
    gradw[i,] <- diag(B)+2*t(w)%*%L*beta/(d*d)
  }
  grad <- convert_to_vec_transfer(m,gradw,gradR)
  return(grad)
}





#=================covariate balancing step===================#

init_algorithm_transfer <- function(O_list,D,meta.list){
    m <- length(O_list)
    d <- nrow(D)
    
    # initialize measurement efficiency only for new datasets
    w_list <- matrix(0,nrow = (m-1),ncol = d)

    for(i in 2:m){
        O1 <- O_list[[1]]
        O2 <- O_list[[i]]
        batchnum <- as.factor(c(rep(1,ncol(O1)),rep(i,ncol(O2))))
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
        meta$weight <- ifelse(meta$batch==i,1/meta$psvalue,1/(1-meta$psvalue))

        # weighting
        O_adj_1 <- t(t(O1)*meta$weight[which(meta$batch==1)])
        O_adj_2 <- t(t(O2)*meta$weight[which(meta$batch==i)])

        # estimated measurement efficiency for dataset
        w_list[(i-1),] <- (rowMeans(O_adj_2)+1e-6)/(rowMeans(O_adj_1)+1e-6)
    }

    # adjust the count list
    O_list_adj <- lapply(2:m,function(i)diag(1/w_list[(i-1),])%*%O_list[[i]])

    # update representation for R using existing dictionary
    R_list <- lapply(1:(m-1),function(i)t(D)%*%O_list_adj[[i]])

    # output
    return(list(w = w_list, R = R_list))
}

#=================helper function of optimization===================#
convert_from_vec_transfer <- function(x, m, r, d, sample_num){
  w.vec <- x[1:(d*m)]
  R_list <- list()
  n0 <- d*m
  for(i in 1:length(sample_num)){
    num <- sample_num[[i]]
    R_list[[i]] <- matrix(x[(n0+1):(n0+num*r)],r,num)
    n0 <- n0+num*r
  }
  w_list <- matrix(w.vec,m,d)
  return(list(w = w_list, R = R_list))
}


convert_to_vec_transfer <- function(m, w_list, R_list){
  x <- c(c(w_list),unlist(R_list))
  return(x)
}
