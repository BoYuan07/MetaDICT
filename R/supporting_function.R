#' @title Taxa/Sample Community detection.
#'
#' @description A \code{k}-nearest neighbor graph is constructed based on Euclidean distance. Then community detection method is applied to identify communities. 
#' Various values of \code{k} within a specific range are tried and the one that yields the highest average Silhouette score is selected.
#'
#' @importFrom RANN nn2
#' @importFrom igraph cluster_louvain
#' @importFrom igraph cluster_walktrap
#'
#' @param X Input data.  
#'   Rows represent clustering objects, and columns represent features.
#' @param max_k The largest number of connected neighbors.
#' @param min_k The smallest number of connected neighbors.
#' @param method The community detection method to use. Options include `"Louvain"` and `"Walktrap"`.
#' @param resolution The resolution parameter for the Louvain algorithm.
#'
#' @returns A \code{list} with the following components:
#' \itemize{
#'   \item{\code{cluster}}{ – The estimated cluster labels.}
#'   \item{\code{graph}}{ – The \code{k}-nearest neighbor graph.}
#' }
#' 
#' @export
community_detection <- function(X, max_k = 10, method = "Louvain", resolution=1, min_k = 2){
  avg_silwidth <- numeric()
  k.list <- c()
  for(k in min_k:max_k){
   membership <- cluster_core(X,k,resolution,method)$membership
   if(length(unique(membership))!=1){
        avg_silwidth <- c(avg_silwidth,mean(cluster::silhouette(membership, dist(X, method = "euclidean"))[,3], na.rm = TRUE))
        k.list <- c(k.list,k)
   }
  }
  best_k <- k.list[which.max(avg_silwidth)]
  best_res <- cluster_core(X,best_k,resolution,method)
  clustring <- paste("Cluster", best_res$membership)
  return(list("cluster" = clustring, "graph" = best_res$graph))
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
    return(list("membership" = km$membership, "graph" = g))
}


#' @title PCoA plots for discrete variables.
#'
#' @importFrom vegan adonis2
#' @importFrom ecodist bcdist
#' @import ggplot2
#'
#' @param X Abundance matrix.  
#'   Rows represent taxa, and columns represent samples.
#' @param covariate A discrete sample covariate.
#' @param title The title of the graph.
#' @param R2 A logical variable.  
#'   Whether to display the R² statistic in the subtitle. Default is \code{TRUE}.
#' @param dissimilarity The dissimilarity type.  
#'   Options include:  
#'   - `"Bray-Curtis"` for Bray-Curtis dissimilarity.  
#'   - `"Euclidean"` for generalized UniFrac dissimilarity.
#' @param colorset The color set for visualization. Default is \code{"Set1"}.
#' @param point_size The size of the points in the plot. Default is \code{1}.
#' 
#' @returns a PCoA plot.
#'
#' @export
pcoa.plot.discrete = function(X, covariate, title, R2 = TRUE, dissimilarity = "Bray-Curtis", colorset = "Set1",point_size = 1){
  if(dissimilarity == "Bray-Curtis"){
     dist_matrix <- bcdist(t(X))
  }else if(dissimilarity == "Euclidean"){
    dist_matrix <- dist(t(X),method = "euclidean")
  }
 
  mds.stuff <- cmdscale(dist_matrix, eig=T, x.ret=T)
  mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
  mds.values <- mds.stuff$points
  mds.data <- data.frame(X=mds.values[,1],
                         Y=mds.values[,2],
                         Sample = covariate)
  p <- ggplot(data=mds.data, aes(x=X, y=Y,color=Sample))+
      geom_point(size=point_size)+
      xlab(paste("PCoA1 -", mds.var.per[1], '%', sep=""))+
      ylab(paste("PCoA2 -", mds.var.per[2], '%', sep=""))+
      labs(title = title)+
      theme_bw()+
      theme(legend.title=element_blank())+
      scale_color_brewer(palette=colorset)+
      theme(legend.key.height=unit(0.5,"cm"),
          legend.title = element_text(size = 14),  
          legend.text  = element_text(size = 12),
          plot.title = element_text(size = 20))

  if(R2){
    r2 <- permanova_pcoa(dist_matrix,covariate)
    p <- p + labs(subtitle = paste("R2 =",round(r2, digits = 4)))+
    theme(plot.subtitle = element_text(size = 16))
  }
  return(p)
}

permanova_pcoa <- function(distP, Y) {
  df.Y = as.data.frame(Y)
  Re = adonis2(distP~Y, data = df.Y)
  return(Re$R2[1])
}

#' @title PCoA plots for continuous variables.
#'
#' @importFrom vegan adonis2
#' @importFrom ecodist bcdist
#' @import ggplot2
#' @import viridis
#'
#' @param X Abundance matrix.  
#'   Rows represent taxa, and columns represent samples.
#' @param covariate A discrete sample covariate.
#' @param title The title of the graph.
#' @param R2 A logical variable.  
#'   Whether to display the R² statistic in the subtitle. Default is \code{TRUE}.
#' @param dissimilarity The dissimilarity type to use.  
#'   Options include:  
#'   - `"Bray-Curtis"` for Bray-Curtis dissimilarity.  
#'   - `"Euclidean"` for generalized UniFrac dissimilarity.
#' @param point_size The size of the points in the plot. Default is \code{1}.
#' 
#' @returns a PCoA plot.
#'
#' @export
pcoa.plot.continuous = function(X, covariate, title, R2 = TRUE, dissimilarity = "Bray-Curtis", point_size = 1){
  if(dissimilarity == "Bray-Curtis"){
     dist_matrix <- bcdist(t(X))
  }else if(dissimilarity == "Euclidean"){
    dist_matrix <- dist(t(X),method = "euclidean")
  }
  mds.stuff <- cmdscale(dist_matrix, eig=T, x.ret=T)
  mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
  mds.values <- mds.stuff$points
  mds.data <- data.frame(X=mds.values[,1],
                         Y=mds.values[,2],
                         Signal = covariate)

  p <- ggplot(data=mds.data, aes(x=X, y=Y,color=Sample))+
      geom_point(size=point_size)+
      scale_colour_gradientn(colors = viridis(10))+
      xlab(paste("PCoA1 -", mds.var.per[1], '%', sep=""))+
      ylab(paste("PCoA2 -", mds.var.per[2], '%', sep=""))+
      labs(title = title)+
      theme_bw()+
      scale_color_brewer(palette=colorset)+
      theme(legend.key.height=unit(0.5,"cm"),
          legend.title = element_text(size = 14),  
          legend.text  = element_text(size = 12),
          plot.title = element_text(size = 20))

  if(R2){
    r2 <- permanova_pcoa(dist_matrix,covariate)
    p <- p + labs(subtitle = paste("R2 =",round(r2, digits = 4)))+
    theme(plot.subtitle = element_text(size = 16))
  }
  return(p)
}


