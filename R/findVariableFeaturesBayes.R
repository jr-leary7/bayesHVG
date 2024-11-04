#' Identify highly variable genes in a Bayesian manner.
#'
#' @name findVariableFeaturesBayes
#' @author Jack R. Leary
#' @description This function implements HVG identification using Bayesian inference to model the posterior distribution of the mean and variance of each gene.
#' @param expr.mat An object of class \code{Seurat} or \code{SingleCellExperiment}. Defaults to NULL.
#' @param subject.id A string specifying the metadata column in \code{expr.mat} that contains subject IDs. Defaults to NULL.
#' @param min.cell.depth
#' @param min.genes.per.cell
#' @param n.cores An integer specifying the number of cores to be used when fitting the Bayesian hierarchical model. Defaults to 4.
#' @import INLA
#' @import magrittr
#' @importFrom SingleCellExperiment colData
#' @importFrom BiocGenerics counts
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom Matrix colSums rowSums
#' @importFrom dplyr mutate select filter rowwise ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom stats as.formula quantile
#' @importFrom withr with_output_sink
#' @importFrom stringr str_detect
#' @return Depending on the input, either an object of class \code{Seurat} or \code{SingleCellExperiment} with HVG metadata added.
#' @seealso \code{\link[Seurat]{FindVariableFeatures}}
#' @export

findVariableFeaturesBayes <- function(expr.mat = NULL,
                                      subject.id = NULL,
                                      min.cell.depth = 1000L,
                                      min.genes.per.cell = 10L,
                                      n.samples = 5000L,
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
  # filter out low quality cells and genes
  keep_cells <- which(Matrix::colSums(expr.mat) >= min.cell.depth)
  keep_genes <- which(Matrix::rowSums(expr.mat > 0) >= min.genes.per.cell)
  expr.mat <- expr.mat[keep_genes, keep_cells]
  # convert to data.frame for modeling
  expr_df <- as.data.frame(expr.mat) %>%
             dplyr::mutate(gene = rownames(.), .before = 1)
  if (!is.null(subject.id)) {
    expr_df <- dplyr::mutate(subject = subject_vec, .before = 2)
    expr_df <- tidyr::pivot_longer(expr_df,
                                   cols = !c(gene, subject),
                                   names_to = "cell",
                                   values_to = "count") %>%
               dplyr::mutate(gene = factor(gene, levels = unique(gene)),
                             subject = factor(subject, levels = unique(subject)))
  } else {
    expr_df <- tidyr::pivot_longer(expr_df,
                                   cols = !gene,
                                   names_to = "cell",
                                   values_to = "count") %>%
               dplyr::with_groups(gene,
                                  dplyr::slice_sample,
                                  n = 500L) %>%
               #dplyr::filter(gene %in% c("CD8A", "NKG7", "MS4A1", "MALAT1", "CD3G", "CD14")) %>%
               dplyr::mutate(gene = factor(gene, levels = unique(gene)))
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
                            data = expr_df,
                            family = "nbinomial",
                            control.compute = list(dic = TRUE, cpo = TRUE),
                            control.predictor = list(compute = TRUE),
                            control.inla = list(strategy = "gaussian"),
                            num.threads = n.cores,
                            verbose = TRUE)
  })
  # extract estimates and generate credible intervals for mean & variance
  phi_row <- dplyr::mutate(bayes_fit$summary.hyperpar,
                           name = rownames(bayes_fit$summary.hyperpar),
                           .before = 1) %>%
             dplyr::filter(stringr::str_detect(name, "phi") | stringr::str_detect(name, "size"))
  phi_mean <- phi_row$mean
  phi_ci_lower <- phi_row$`0.025quant`
  phi_ci_upper <- phi_row$`0.975quant`
  gene_effects <- bayes_fit$marginals.random$gene
  gene_names_indices <- names(gene_effects)
  gene_indices <- as.numeric(stringr::str_extract(gene_names_indices, "\\d+"))
  gene_names <- levels(expr_df$gene)[gene_indices]
  sampleMarginal <- function(marginal, n = 1000L) {
    sample_res <- INLA::inla.rmarginal(n, marginal = marginal)
    return(sample_res)
  }
  mu_samples <- sapply(gene_effects, \(x) {
    samples_log_mu <- sampleMarginal(x, n = n.samples)
    samples_mu <- exp(samples_log_mu)
    return(samples_mu)
  })
  phi_marginal <- bayes_fit$marginals.hyperpar[[phi_row$name]]
  phi_samples <- sampleMarginal(phi_marginal, n_samples)
  phi_matrix <- matrix(phi_samples,
                       nrow = n_samples,
                       ncol = length(gene_names),
                       byrow = FALSE)
  var_samples <- mu_samples + (mu_samples^2) / phi_matrix
  dispersion_samples <- var_samples / mu_samples
  gene_summary <- data.frame(gene = gene_names,
                             mu_mean = apply(mu_samples, 2, mean),
                             mu_ci_lower = apply(mu_samples, 2, stats::quantile, probs = 0.025),
                             mu_ci_upper = apply(mu_samples, 2, stats::quantile, probs = 0.975),
                             var_mean = apply(var_samples, 2, mean),
                             var_ci_lower = apply(var_samples, 2, stats::quantile, probs = 0.025),
                             var_ci_upper = apply(var_samples, 2, stats::quantile, probs = 0.975),
                             dispersion_mean = apply(dispersion_samples, 2, mean),
                             dispersion_ci_lower = apply(dispersion_samples, 2, stats::quantile, probs = 0.025),
                             dispersion_ci_upper = apply(dispersion_samples, 2, stats::quantile, probs = 0.975)) %>% View()


  compute_mu <- function(effect.marginal = NULL){
    mu_mean_log <- INLA::inla.emarginal(\(x) x, effect.marginal)
    mu_ci_lower_log <- INLA::inla.qmarginal(0.025, effect.marginal)
    mu_ci_upper_log <- INLA::inla.qmarginal(0.975, effect.marginal)
    mu_mean <- exp(mu_mean_log)
    mu_ci_lower <- exp(mu_ci_lower_log)
    mu_ci_upper <- exp(mu_ci_upper_log)
    return(c(mu_mean, mu_ci_lower, mu_ci_upper))
  }
  mu_values <- t(sapply(gene_effects, compute_mu))
  colnames(mu_values) <- c("mu_mean", "mu_ci_lower", "mu_ci_upper")
  gene_summary <- data.frame(gene = gene_names,
                             mu_mean = mu_values[, "mu_mean"],
                             mu_ci_lower = mu_values[, "mu_ci_lower"],
                             mu_ci_upper = mu_values[, "mu_ci_upper"]) %>%
                  dplyr::rowwise() %>%
                  dplyr::mutate(var_mean = mu_mean + (mu_mean^2) / phi_mean,
                                var_ci_lower = mu_ci_lower + (mu_ci_lower^2) / phi_ci_upper,
                                var_ci_upper = mu_ci_upper + (mu_ci_upper^2) / phi_ci_lower) %>%
                  dplyr::ungroup() %>%
                  dplyr::mutate(dispersion = var_mean / mu_mean) %>%
                  as.data.frame()

}
