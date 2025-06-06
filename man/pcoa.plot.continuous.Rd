% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supporting_function.R
\name{pcoa.plot.continuous}
\alias{pcoa.plot.continuous}
\title{PCoA plots for continuous variables.}
\usage{
pcoa.plot.continuous(
  X,
  covariate,
  title,
  R2 = TRUE,
  dissimilarity = "Bray-Curtis",
  point_size = 1
)
}
\arguments{
\item{X}{Abundance matrix.
Rows represent taxa, and columns represent samples.}

\item{covariate}{A discrete sample covariate.}

\item{title}{The title of the graph.}

\item{R2}{A logical variable.
Whether to display the R² statistic in the subtitle. Default is \code{TRUE}.}

\item{dissimilarity}{The dissimilarity type to use.
Options include:
\itemize{
\item \code{"Bray-Curtis"} for Bray-Curtis dissimilarity.
\item \code{"Euclidean"} for generalized UniFrac dissimilarity.
}}

\item{point_size}{The size of the points in the plot. Default is \code{1}.}
}
\value{
a PCoA plot.
}
\description{
PCoA plots for continuous variables.
}
