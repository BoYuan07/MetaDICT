#' @title Microbiome data integration method via shared dictionary learning.
#'
#' @description A method for the integration of microbiome data. This method is designed to remove batch effects and preserve biological variation while integrating heterogeneous datasets. MetaDICT can better avoid overcorrection when unobserved confounding variables are present.
#'
#' @details MetaDICT is a two-step approach. It initially estimates the batch effects by weighting methods in causal inference literature then refines the estimation via a shared dictionary learning.
#'
#' @importFrom stats optim
#'
#' @examples
#'
#'  data(exampleData)
#'  alpha = 0.01
#'  beta = 0.01
#'  gamma = 1
#'  metadict_res = metadict(O_list,dist,meta.list = meta.list)
#'
#' @param O_list A \code{list} object contains count tables for different datasets. The rows represent taxa, and the columns represent samples. All count tables must have the same taxa.
#' @param alpha The parameter is related to the rank of the final corrected count table. A larger alpha leads to a lower rank in the corrected count table.
#' @param beta The parameter is related to the smoothness of the estimated measurement efficiency. A larger beta results in more similar measurement efficiency across taxa.
#' @param gamma The parameter is related to how well the model fits the data.
#' @param dist Sequence dissimilarity matrix. Dissimilarity can be measured by phylogenetic distance or taxonomic dissimilarity.
#' @param meta.list The \code{list} object contains metadata tables for different datasets.
#' @param normalization Normalization method. Upper quantile or RSim.
#' @param max_iter Maximum number of iterations for the optimization process.
#' @param neighbor The number of nearest neighbors used to construct the taxa neighbor graph. The default is 5.
#' @param sigma Hyperparameter for the Gaussian kernel used to determine graph edge weights.
#'
#' @returns a \code{list} with components:
#' \itemize{
#' \item{\code{X}, a matrix. Integrated count table. Rows represent taxa. Columns represent samples.}
#' \item{\code{D}, a matrix. Estimated shared dictionary.}
#' \item{\code{R}, a matrix. Estimated sample representation.}
#' \item{\code{w}, a matrix. Estimated measurement efficiency. Rows represent dataset, columns represent taxa.}
#' \item{\code{meta}, a data frame. Meta table without batch ID. Sample features must be present in all datasets. Default is NULL.}
#' }
#'
#' @export

metadict = function(O_list,dist,meta.list = NULL, alpha = 1,beta = 0.01, gamma = 1, normalization = "uq", max_iter=10000,neighbor=5,sigma=10){
  m = length(O_list)
  d = nrow(dist)

  O = do.call("cbind",O_list)
  r = sum(svd(O)$d>1e-3)
  
  sample_num = sapply(O_list,ncol)
  if(normalization == "uq"){
      O_norm = lapply(O_list,function(x)uq(x)$P)
  }else if(normalization == "rsim"){
      O_norm = lapply(O_list,function(x)rsim(x)$P)
  }
  else if(!normalization){
      O_norm = O_list
  }
  
  scale = max(unlist(O_norm))
  O_list_test = lapply(O_norm,function(x)x/scale)
  A = matrix(0,d,d)
  for(i in 1:d){
    idx = order(dist[i,],decreasing = F)[2:(neighbor+1)]
    A[i,idx] = exp(-dist[i,idx]/sigma)
    A[idx,i] = exp(-dist[i,idx]/sigma)
  }
  D1 = diag(rowSums(A))
  L = D1-A
  initial = init_algo(m,d,O_list_test,r,meta.list)

  x0 = convert_to_vec(m,initial[[1]],initial[[2]],initial[[3]])
  lower = c(rep(0,m*d),rep(-Inf,d*r+r*sum(sample_num)))
  upper = c(rep(1,m*d),rep(Inf,d*r+r*sum(sample_num)))
  optim.res = optim(x0, fn = (function(x) target_func(x,O_list_test,alpha,beta,gamma,m,d,r,sample_num,L)), gr = (function(x) gradient_func(x,O_list_test,alpha,beta,gamma,m,d,r,sample_num,L)), method = "L-BFGS-B",
      lower = lower, upper = upper, control = list(maxit = max_iter))
  if(optim.res$convergence==0){
      print("Successful convergence.")
  }
  if(optim.res$convergence==1){
      print("The iteration limit maxit had been reached.")
  }
  x = optim.res$par
  para = convert_from_vec(x,m,r,d,sample_num)
  w_list = para[[1]]
  D = para[[2]]
  R_list = para[[3]]
  X_list = list()
  for(i in 1:m){
    X_list[[i]] = D%*%R_list[[i]]*scale
  }
  X = do.call(cbind,X_list)
  X[X<0] = 0
  res = list()
  res[["X"]] = X
  
  var_res = varimax(D, normalize = FALSE)
  D_rot = var_res$loadings
  T_mat = var_res$rotmat
  res[["D"]] = D_rot
  res[["R"]] = lapply(R_list,function(x)t(T_mat)%*%x)
  res[["w"]] = w_list
  res[["meta"]] = do.call(rbind,meta.list)
  return(res)
}

target_func = function(x,O_list,alpha,beta,gamma,m,d,r,sample_num,L){
  res = convert_from_vec(x,m,r,d,sample_num)
  w_list = res[[1]]
  D = res[[2]]
  R_list = res[[3]]
  target = 0
  for(i in 1:m){
    w = w_list[i,]
    W = diag(w)
    O_diff = O_list[[i]]-W%*%D%*%R_list[[i]]
    W_diff = t(w)%*%L%*%w
    target = target+gamma*norm(O_diff,"F")**2+beta*W_diff/(d*d)+norm(R_list[[i]],"F")**2*alpha/(2*r*sample_num[i])
    
  }
  target = target+norm(D,"F")**2*alpha/(2*d*r)
  return(target)
}

gradient_func = function(x,O_list,alpha,beta,gamma,m,d,r,sample_num,L){
  res = convert_from_vec(x,m,r,d,sample_num)
  w_list = res[[1]]
  D = res[[2]]
  R_list = res[[3]]
  gradw = matrix(0,m,d)
  gradD = matrix(0,d,r)
  gradR = list()
  for(i in 1:m){
    w = w_list[i,]
    W = diag(w)
    O_diff = O_list[[i]]-W%*%D%*%R_list[[i]]
    W_diff = t(w)%*%L%*%w
    gradD = gradD+gamma*(-2*W%*%O_diff%*%t(R_list[[i]]))
    gradR[[i]] = gamma*(-2)*t(D)%*%W%*%O_diff+R_list[[i]]*alpha/(r*sample_num[i])
    B = gamma*(-2)*(O_diff)%*%t(R_list[[i]])%*%t(D)
    gradw[i,] = diag(B)+2*t(w)%*%L*beta/(d*d)
  }
  gradD = gradD+alpha*D/(d*r)
  grad = convert_to_vec(m,gradw,gradD,gradR)
  return(grad)
}

init_algo = function(m,d,O_list,r,meta.list){
  w_list = matrix(0,nrow = m,ncol = d)
  w_list[1,] = 1
  for(i in 2:m){
    O1 = O_list[[1]]
    O2 = O_list[[i]]
    batchnum = as.factor(c(rep(1,ncol(O1)),rep(i,ncol(O2))))
    if(is.null(meta.list)){
      meta = data.frame("lib" = c(colSums(O1),colSums(O2)),"batch" = batchnum)
    }else{
      meta = rbind(meta.list[[1]],meta.list[[i]])
      meta$batch = batchnum
      meta$lib = c(colSums(O1),colSums(O2))
    }
    mylogit <- glm(batch ~ ., data = meta, family = "binomial")
    meta$psvalue = predict(mylogit, type="response")
    meta$weight = ifelse(meta$batch==i,1/meta$psvalue,1/(1-meta$psvalue))
    O_adj_1 = t(t(O1)*meta$weight[which(meta$batch==1)])
    O_adj_2 = t(t(O2)*meta$weight[which(meta$batch==i)])
    w_list[i,] = (rowMeans(O_adj_2)+1e-6)/(rowMeans(O_adj_1)+1e-6)
  }
  w_list = t(t(w_list)/colSums(w_list))
  O_list1 = lapply(1:m,function(i)diag(1/w_list[i,])%*%O_list[[i]])
  O = do.call(cbind,O_list1)
  svd.res = svd(O)
  D = svd.res$u[,1:r]
  R_list = lapply(1:m,function(i)t(D)%*%O_list1[[i]])
  initial_val = list()
  initial_val[[1]] = w_list
  initial_val[[2]] = D
  initial_val[[3]] = R_list
  return(initial_val)
}

convert_from_vec = function(x,m,r,d,sample_num){
  x.w_list = x[1:(d*m)]
  x.D = x[(d*m+1):(d*m+d*r)]
  R_list = list()
  n0 = d*m+d*r
  for(i in 1:length(sample_num)){
    num = sample_num[[i]]
    R_list[[i]] = matrix(x[(n0+1):(n0+num*r)],r,num)
    n0 = n0+num*r
  }
  
  D = matrix(x.D,d,r)
  w_list = matrix(x.w_list,m,d)
  return(list(w_list,D,R_list))
}

convert_to_vec = function(m,w_list,D,R_list){
  x = c(c(w_list),c(D),unlist(R_list))
  return(x)
}

convert_to_vec = function(m,w_list,D,R_list){
  x = c(c(w_list),c(D),unlist(R_list))
  return(x)
}

uq<-function(X){
    dds = edgeR::calcNormFactors(as.matrix(X+1), method = "upperquartile")
    upperQ = dds*colSums(X)
    upq.res <- scale(X,center=FALSE,scale=upperQ)
  return(list('P' = upq.res, 'sf' = upperQ))
}

CStat = function(X){
  d <- nrow(X)
  R <- X
  S1 <- apply(R,1,order)
  S1 <- S1 - colMeans(S1);
  S1 <- S1 / sqrt(colSums(S1^2));
  corr_s <- crossprod(S1)
  med <- as.data.frame(matrixStats::colMedians(corr_s))
  return(as.numeric(med[,1]))

}

rsim = function(X,eta=0.01){
    d = nrow(X)
    v = CStat(X)
    I0.1 = which(v>0.8)
    X0 = X[I0.1,]
    v0 = replicate(3,CStat(X0[sample(1:nrow(X0),0.5*nrow(X0)),]))
    w = v[v>0.8]
    f1 = sapply(w,function(x)mean(v>x))
    f0 = sapply(w,function(x)mean(v0>x))
    pi = sum(f1*f0)/sum(f0^2)
    vord = order(v,decreasing = T)
    res = sapply(1:length(vord),function(x)(1-pi*length(vord)*mean(v0>v[x])/(which(vord==x))))
    lowerx = max(which(res[vord]<eta))
    ref = vord[1:lowerx]
    tc.cn <- apply(X,2,function(x)sum(x[ref]))
    f.cn <- tc.cn/(mean(tc.cn))
    f.cn <- ifelse(f.cn==0,1,f.cn)
    cn.res <- scale(X,center=FALSE,scale=f.cn)
  return(list('P' = cn.res, 'I0' = ref, 'pi0'= pi, 'sf'=f.cn))
}

