#' Identify highly variable genes in a Bayesian manner.
#'
#' @name findVariableFeaturesBayes
#' @author Jack R. Leary
#' @description This function implements HVG identification using Bayesian inference to model the posterior distribution of the mean and variance of each gene.
#' @param expr.mat Either an object of class \code{Seurat} or \code{SingleCellExperiment}. Defaults to NULL. 
#' @param subject.id A string specifying the metadata column that contains subject IDs. Defaults to NULL.
#' @param n.cores An integer specifying the number of cores to be used when fitting the Bayesian hierarchical model. Defaults to 4.
#' @import INLA 
#' @import magrittr
#' @importFrom SingleCellExperiment colData 
#' @importFrom BiocGenerics counts 
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom dplyr mutate select 
#' @importFrom tidyr pivot_longer 
#' @importFrom stats as.formula quantile 
#' @importFrom withr with_output_sink
#' @importFrom tidyselect starts_with 
#' @return Depending on the input, either an object of class \code{Seurat} or \code{SingleCellExperiment} with HVG metadata added. 
#' @seealso \code{\link[Seurat]{FindVariableFeatures}}
#' @export 

findVariableFeaturesBayes <- function(expr.mat = NULL, 
                                      subject.id = NULL, 
                                      n.cores = 4L) {
  # check inputs 
  if (is.null(expr.mat)) { stop("Please provide all inputs to findVariableFeaturesBayes().") }
  # extract (sparse) counts matrix 
  if (inherits(expr.mat, "SingleCellExperiment")) {
    if (!is.null(subject.id)) {
      subject_vec <- SingleCellExperiment::colData(expr.mat)[[subject.id]]
    }
    expr.mat <- BiocGenerics::counts(expr.mat)
  } else if (inherits(expr.mat, "Seurat")) {
    if (!is.null(subject.id)) {
      subject_vec <- expr.mat@meta.data[[subject.id]]
    }
    expr.mat <- Seurat::GetAssayData(expr.mat,
                                     layer = "counts",
                                     assay = Seurat::DefaultAssay(expr.mat))
  }
  # convert to data.frame for modeling 
  expr_df <- as.data.frame(expr.mat) %>% 
             dplyr::mutate(gene = rownames(.), .before = 1)
  if (!is.null(subject.id)) {
    expr_df <- dplyr::mutate(subject = subject_vec, .before = 2)
    expr_df <- tidyr::pivot_longer(expr_df, 
                                   cols = !c(gene, subject),  
                                   names_to = "cell", 
                                   values_to = "count")
  } else {
    expr_df <- tidyr::pivot_longer(expr_df, 
                                   cols = !gene,  
                                   names_to = "cell", 
                                   values_to = "count")
  }
  # create model formula 
  if (!is.null(subject.id)) {
    model_formula <- stats::as.formula("count ~ 1 + f(gene, model = 'iid') + f(subject, model = 'iid')")
  } else {
    model_formula <- stats::as.formula("count ~ 1 + f(gene, model = 'iid')")
  }
  # fit hierarchical bayesian model via integrated nested laplace approximation
  withr::with_output_sink(tempfile(), {
    bayes_fit <- INLA::inla(model_formula, 
                            data = expr_df2,
                            family = "nbinomial",
                            control.compute = list(dic = TRUE, cpo = TRUE),
                            control.predictor = list(compute = TRUE), 
                            control.inla = list(strategy = "gaussian"), 
                            num.threads = 4L)
  })
  # extract estimates and generate credible intervals 
  gene_effects <- bayes_fit$summary.random$gene
  gene_summary <- data.frame(gene = gene_effects$ID,
                             mean = exp(gene_effects$mean),
                             mean_ci_lower = exp(gene_effects$`0.025quant`),
                             mean_ci_upper = exp(gene_effects$`0.975quant`),
                             var = exp(gene_effects$sd^2 + 2 * gene_effects$mean) - exp(2 * gene_effects$mean), 
                             var_ci_lower = exp((gene_effects$`0.025quant`)^2) - exp((gene_effects$`0.025quant`) * 2),
                             var_ci_upper = exp((gene_effects$`0.975quant`)^2) - exp((gene_effects$`0.975quant`) * 2)) %>% 
                  dplyr::mutate(dispersion = var / mean) %>% 
                  magrittr::set_rownames(.$gene)
}
