% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supporting_function.R
\name{pcoa.plot.discrete}
\alias{pcoa.plot.discrete}
\title{PCoA plots for discrete variable with bray-curtis distance.}
\usage{
pcoa.plot.discrete(A, covariate, main, colorset = "Set1")
}
\arguments{
\item{A}{Abundance matrix. The rows represent taxa, the columns represent samples.}

\item{covariate}{A discrete sample covariate.}

\item{main}{Graph title.}

\item{colorset}{Color set. Default is \code{Set1}.}
}
\value{
a PCoA plot.
}
\description{
PCoA plots for discrete variable with bray-curtis distance.
}
