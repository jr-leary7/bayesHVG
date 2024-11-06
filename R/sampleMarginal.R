#' Sample from the marginal distribution of a fitted model.
#'
#' @name sampleMarginal
#' @author Jack R. Leary
#' @description A simple wrapper function to generate a given number of samples from the marginal distribution of a fitted model.
#' @param marginal.dist The marginal distribution to be sampled from. Defaults to NULL.
#' @param n An integer specifying the number of samples. Defaults to 1000.
#' @import INLA
#' @return A numeric vector of samples.
#' @seealso \code{\link[INLA]{inla.rmarginal}}
#' @seealso \code{\link[INLA]{inla.is.marginal}}
#' @export

sampleMarginal <- function(marginal.dist = NULL, n = 1000L) {
  # check inputs
  if (is.null(marginal.dist) || !INLA::inla.is.marginal(marginal.dist)) { stop("Please provide a valid marginal distribution to sampleMarginal().") }
  # generate samples
  marginal_samples <- INLA::inla.rmarginal(n, marginal = marginal.dist)
  return(marginal_samples)
}
