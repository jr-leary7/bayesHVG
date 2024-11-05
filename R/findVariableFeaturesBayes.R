#' Identify highly variable genes in a Bayesian manner.
#'
#' @name findVariableFeaturesBayes
#' @author Jack R. Leary
#' @description This function implements HVG estimation using Bayesian inference to model the posterior distribution of the mean and variance of each gene.
#' @param sc.obj An object of class \code{Seurat} or \code{SingleCellExperiment}. Defaults to NULL.
#' @param subject.id A string specifying the metadata column in \code{expr.mat} that contains subject IDs. Defaults to NULL.
#' @param n.marginal.samples An integer specifying the number of samples to take from the marginal distribution when computing estimates and credible intervals. Defults to 5000.
#' @param n.cells.subsample An integer specifying the number of cells per-gene to subsample to when performing estimation. Defaults to 500.
#' @param n.cores An integer specifying the number of cores to be used when fitting the Bayesian hierarchical model. Defaults to 4.
#' @import INLA
#' @import magrittr
#' @importFrom SingleCellExperiment colData
#' @importFrom BiocGenerics counts
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom dplyr mutate select with_groups slice_sample filter arrange desc
#' @importFrom tidyr pivot_longer
#' @importFrom stats as.formula quantile
#' @importFrom withr with_output_sink
#' @importFrom stringr str_detect
#' @return Depending on the input, either an object of class \code{Seurat} or \code{SingleCellExperiment} with HVG metadata added.
#' @seealso \code{\link[Seurat]{FindVariableFeatures}}
#' @seealso \code{\link[scran]{modelGeneVar}}
#' @export

findVariableFeaturesBayes <- function(sc.obj = NULL,
                                      subject.id = NULL,
                                      n.marginal.samples = 5000L,
                                      n.cells.subsample = 500L,
                                      n.cores = 4L) {
  # check inputs
  if (is.null(sc.obj)) { stop("Please provide all inputs to findVariableFeaturesBayes().") }
  # extract (sparse) counts matrix
  if (inherits(sc.obj, "SingleCellExperiment")) {
    if (!is.null(subject.id)) {
      subject_vec <- SingleCellExperiment::colData(sc.obj)[[subject.id]]
    }
    expr_mat <- BiocGenerics::counts(sc.obj)
  } else if (inherits(sc.obj, "Seurat")) {
    if (!is.null(subject.id)) {
      subject_vec <- sc.obj@meta.data[[subject.id]]
    }
    expr_mat <- Seurat::GetAssayData(sc.obj,
                                     layer = "counts",
                                     assay = Seurat::DefaultAssay(sc.obj))
  }
  # convert to data.frame for modeling
  expr_df <- as.data.frame(expr_mat) %>%
             dplyr::mutate(gene = rownames(.), .before = 1)
  if (!is.null(subject.id)) {
    expr_df <- dplyr::mutate(expr_df,
                             subject = subject_vec,
                             .before = 2) %>%
               tidyr::pivot_longer(cols = !c(gene, subject),
                                   names_to = "cell",
                                   values_to = "count") %>%
               dplyr::with_groups(c(gene, subject),
                                  dplyr::slice_sample,
                                  n = n.cells.subsample) %>%
               dplyr::mutate(gene = factor(gene, levels = unique(gene)),
                             subject = factor(subject, levels = unique(subject)))
  } else {
    expr_df <- tidyr::pivot_longer(expr_df,
                                   cols = !gene,
                                   names_to = "cell",
                                   values_to = "count") %>%
               dplyr::with_groups(gene,
                                  dplyr::slice_sample,
                                  n = n.cells.subsample) %>%
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
                            control.compute = list(dic = TRUE, cpo = FALSE),
                            control.predictor = list(compute = TRUE, link = 1),
                            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                            control.family = list(hyper = list(theta = list(prior = "loggamma", param = c(1, 1)))),
                            num.threads = n.cores,
                            verbose = TRUE,
                            debug = TRUE)
  })
  # extract estimates and generate credible intervals for mean & variance
  phi_row <- dplyr::mutate(bayes_fit$summary.hyperpar,
                           name = rownames(bayes_fit$summary.hyperpar),
                           .before = 1) %>%
             dplyr::filter(stringr::str_detect(name, "size"))
  gene_effects <- bayes_fit$marginals.random$gene
  gene_names_indices <- names(gene_effects)
  gene_indices <- as.numeric(stringr::str_extract(gene_names_indices, "\\d+"))
  gene_names <- levels(expr_df$gene)[gene_indices]
  sampleMarginal <- function(marginal, n = 1000L) {
    sample_res <- INLA::inla.rmarginal(n, marginal = marginal)
    return(sample_res)
  }
  intercept_marginal <- bayes_fit$marginals.fixed[["(Intercept)"]]
  intercept_samples <- sampleMarginal(intercept_marginal, n = n.marginal.samples)
  mu_samples <- sapply(gene_effects, \(x) {
    samples_log_re <- sampleMarginal(x, n = n.marginal.samples)
    samples_log_mu <- intercept_samples + samples_log_re
    samples_mu <- exp(samples_log_mu)
    return(samples_mu)
  })
  phi_marginal <- bayes_fit$marginals.hyperpar[[phi_row$name]]
  phi_samples <- sampleMarginal(phi_marginal, n.marginal.samples)
  phi_matrix <- matrix(phi_samples,
                       nrow = n.marginal.samples,
                       ncol = length(gene_names),
                       byrow = FALSE)
  var_samples <- mu_samples + (mu_samples^2 / phi_matrix)
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
                             dispersion_ci_upper = apply(dispersion_samples, 2, stats::quantile, probs = 0.975)) %>%
                  magrittr::set_rownames(.$gene) %>%
                  dplyr::arrange(dplyr::desc(dispersion_mean))
  # add gene-level estimates to object metadata -- TODO
}
