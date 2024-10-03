#' @title Community detection.
#'
#' @description A \code{k}-nearest neighbor graph based on Euclidean distance. Then community detection method is applied to identify communities. Various values of \code{k} within a specific range are tried and the one that yields the highest average Silhouette score is selected.
#'
#' @importFrom RANN nn2
#' @importFrom igraph cluster_louvain
#' @importFrom igraph cluster_walktrap
#'
#' @param X Representation. The rows represent samples, the columns represent features.
#' @param max_k The largest number of connected neighbors.
#' @param min_k The smallest number of connected neighbors.
#' @param method Community detection methods such as Louvain and Walktrap.
#' @param resolution The resolution parameter of Louvain algorithm.
#'
#' @returns a \code{list} with components:
#' \itemize{
#' \item{\code{cluster}, estimated cluster labels.}
#' \item{\code{graph}, the \code{k}-nearest neighbor graph.}
#' }
#'
#' @export

community_detection = function(X, max_k = 50, resolution=1, method = "Louvain", min_k = 5){
  avg_silwidth = numeric()
  k.list = c()
  for(k in min_k:max_k){
   membership = cluster_core(X,k,resolution,method)$membership
   if(length(unique(membership))!=1){
        avg_silwidth= c(avg_silwidth,mean(cluster::silhouette(membership, dist(X, method = "euclidean"))[,3], na.rm = TRUE))
        k.list = c(k.list,k)
   }
  }
  best_k = k.list[which.max(avg_silwidth)]
  best_res = cluster_core(X,best_k,resolution,method)
  return(list("cluster" = best_res$membership, "graph" = best_res$graph))
}

cluster_core = function(X,k,resolution,method = "Louvain"){
    knn.info <- RANN::nn2(X, k=k)
    knn <- knn.info$nn.idx
    adj <- matrix(0, nrow(X), nrow(X))
    for(i in seq_len(nrow(X))) {
        adj[i,knn[i,2:k]] <- 1
    }
    g <- igraph::graph.adjacency(adj, mode="undirected")
    g <- igraph::simplify(g)
    
    if(method == "Louvain"){
        km <- igraph::cluster_louvain(g,resolution = resolution)
    }else if(method == "Walktrap"){
        km <- igraph::cluster_walktrap(g)
    }
    return(list("membership"=km$membership,"graph"=g))
}


#' @title PCoA plots for discrete variable with bray-curtis distance.
#'
#' @importFrom vegan adonis2
#' @importFrom ecodist bcdist
#' @import ggplot2
#'
#' @param A Abundance matrix. The rows represent taxa, the columns represent samples.
#' @param covariate A discrete sample covariate.
#' @param main Graph title.
#' @param colorset Color set. Default is \code{Set1}.
#'
#' @returns a PCoA plot.
#'
#' @export
pcoa.plot.discrete = function(A,covariate,main,colorset = "Set1"){
  dist_matrix = bcdist(t(A))
  mds.stuff = cmdscale(dist_matrix, eig=T, x.ret=T)
  mds.var.per = round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
  mds.values = mds.stuff$points
  mds.data = data.frame( X=mds.values[,1],
                         Y=mds.values[,2],
                         Sample = covariate)

  r2 = permanova_pcoa(dist_matrix,covariate)
  ggplot(data=mds.data, aes(x=X, y=Y,color=Sample))+
    geom_point(size=1)+
    xlab(paste("PCoA1 -", mds.var.per[1], '%', sep=""))+
    ylab(paste("PCoA2 -", mds.var.per[2], '%', sep=""))+
    labs(title = main,
              subtitle = paste("R2 =",round(r2, digits = 4)))+
    theme_bw()+
    theme(legend.title=element_blank())+
    scale_color_brewer(palette=colorset)+
    theme(legend.key.height=unit(0.5,"cm"))
}

permanova_pcoa <- function(distP, Y) {
  df.Y = as.data.frame(Y)
  Re = adonis2(distP~Y, data = df.Y)
  return(Re$R2[1])
}

#' @title PCoA plots for continuous variable with bray-curtis distance.
#'
#' @importFrom vegan adonis2
#' @importFrom ecodist bcdist
#' @import ggplot2
#' @import viridis
#'
#' @param A Abundance matrix. The rows represent taxa, the columns represent samples.
#' @param covariate A continuous sample covariate.
#' @param main Graph title.
#'
#' @returns a PCoA plot.
#'
#' @export
pcoa.plot.continuous = function(A,covariate,main){
  dist_matrix = bcdist(t(A))
  mds.stuff = cmdscale(dist_matrix, eig=T, x.ret=T)
  mds.var.per = round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
  mds.values = mds.stuff$points
  mds.data = data.frame(X=mds.values[,1],
                         Y=mds.values[,2],
                         Signal = covariate)
  r2 = permanova_pcoa(dist_matrix,covariate)
  ggplot(data=mds.data, aes(x=X, y=Y,col=Signal))+
    geom_point(size=1)+
    scale_colour_gradientn(colors = viridis(10))+
    theme(legend.title=element_blank()) +
    xlab(paste("PCoA1 -", mds.var.per[1], '%', sep=""))+
    ylab(paste("PCoA2 -", mds.var.per[2], '%', sep=""))+
    labs(title = main,
              subtitle = paste("R2 =",round(r2, digits = 4)))+
    theme_bw()+
    theme(legend.key.height=unit(0.5,"cm"))
}


