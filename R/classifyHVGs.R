#' Classify genes as highly variable.
#'
#' @name classifyHVGs
#' @author Jack R. Leary
#' @description After estimating per-gene statistics, classify genes as highly variable or not in a variety of ways.
#' @param sc.obj An object of class \code{Seurat} or \code{SingleCellExperiment}. Defaults to NULL.
#' @import magrittr
#' @importFrom SingleCellExperiment rowData
#' @importFrom dplyr
#' @importFrom stats quantile
#' @importFrom S4Vectors DataFrame
#' @return Depending on the input, either an object of class \code{Seurat} or \code{SingleCellExperiment} with HVG metadata added.
#' @seealso \code{\link{findVariableFeaturesBayes}}
#' @export

classifyHVGs <- function(sc.obj = NULL) {

}
