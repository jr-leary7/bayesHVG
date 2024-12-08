% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeNaiveGeneStatistics.R
\name{computeNaiveGeneStatistics}
\alias{computeNaiveGeneStatistics}
\title{Compute the naive mean, variance, and dispersion of gene expression.}
\usage{
computeNaiveGeneStatistics(sc.obj = NULL)
}
\arguments{
\item{sc.obj}{An object of class \code{Seurat} or \code{SingleCellExperiment}. Defaults to NULL.}
}
\value{
A \code{data.frame} containing the per-gene naive mean, variance, and dispersion.
}
\description{
Useful as a comparison to the Bayesian results generated by \code{\link{findVariableFeaturesBayes}}, this function computes the Frequentist counterparts to the statistics estimated by that function.
}
\author{
Jack R. Leary
}
