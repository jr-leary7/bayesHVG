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
#' @importFrom SingleCellExperiment colData rowData
#' @importFrom BiocGenerics counts
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom dplyr mutate select with_groups slice_sample filter arrange desc left_join
#' @importFrom tidyr pivot_longer
#' @importFrom stats as.formula quantile
#' @importFrom withr with_output_sink
#' @importFrom stringr str_detect
#' @importFrom purrr map_dbl
#' @importFrom S4Vectors DataFrame
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
                            control.compute = list(dic = FALSE, cpo = FALSE),
                            control.predictor = list(compute = TRUE, link = 1),
                            control.inla = list(strategy = "simplified.laplace", int.strategy = "eb"),
                            control.family = list(hyper = list(theta = list(prior = "loggamma", param = c(1, 1)))),
                            num.threads = n.cores,
                            verbose = TRUE,
                            debug = TRUE)
  })
  # sample from marginal distribution for gene mean 
  intercept_marginal <- bayes_fit$marginals.fixed[["(Intercept)"]]
  intercept_samples <- sampleMarginal(intercept_marginal, n = n.marginal.samples)
  gene_effects <- bayes_fit$marginals.random$gene
  gene_names_indices <- names(gene_effects)
  gene_indices <- as.numeric(stringr::str_extract(gene_names_indices, "\\d+"))
  gene_names <- levels(expr_df$gene)[gene_indices]
  mu_samples <- sapply(gene_effects, \(x) {
    samples_log_re <- sampleMarginal(x, n = n.marginal.samples)
    samples_log_mu <- intercept_samples + samples_log_re
    samples_mu <- exp(samples_log_mu)
    return(samples_mu)
  })
  # sample from marginal distribution for overdispersion parameter 
  phi_row <- dplyr::mutate(bayes_fit$summary.hyperpar,
                           name = rownames(bayes_fit$summary.hyperpar),
                           .before = 1) %>%
             dplyr::filter(stringr::str_detect(name, "size"))
  phi_marginal <- bayes_fit$marginals.hyperpar[[phi_row$name]]
  phi_samples <- sampleMarginal(phi_marginal, n.marginal.samples)
  phi_matrix <- matrix(phi_samples,
                       nrow = n.marginal.samples,
                       ncol = length(gene_names),
                       byrow = FALSE)
  # estimate variance & dispersion samples based on formula for negative-binomial variance 
  var_samples <- mu_samples + (mu_samples^2 / phi_matrix)
  dispersion_samples <- var_samples / mu_samples
  mu_samples <- as.data.frame(mu_samples)
  var_samples <- as.data.frame(var_samples)
  dispersion_samples <- as.data.frame(dispersion_samples)
  # generate central tendency estimates and credible intervals for each parameter 
  gene_summary <- data.frame(gene = gene_names, 
                             mu = purrr::map_dbl(mu_samples, mean), 
                             mu_ci_lower = purrr::map_dbl(mu_samples, \(x) stats::quantile(x, probs = 0.025)), 
                             mu_ci_upper = purrr::map_dbl(mu_samples, \(x) stats::quantile(x, probs = 0.975)), 
                             var = purrr::map_dbl(var_samples, mean), 
                             var_ci_lower = purrr::map_dbl(var_samples, \(x) stats::quantile(x, probs = 0.025)), 
                             var_ci_upper = purrr::map_dbl(var_samples, \(x) stats::quantile(x, probs = 0.975)), 
                             dispersion = purrr::map_dbl(dispersion_samples, mean), 
                             dispersion_ci_lower = purrr::map_dbl(dispersion_samples, \(x) stats::quantile(x, probs = 0.025)), 
                             dispersion_ci_upper = purrr::map_dbl(dispersion_samples, \(x) stats::quantile(x, probs = 0.975))) %>%
                  magrittr::set_rownames(.$gene) %>%
                  dplyr::arrange(dplyr::desc(dispersion))
  # add gene-level estimates to object metadata
  if (inherits(sc.obj, "SingleCellExperiment")) {
    gene_summary_s4 <- SingleCellExperiment::rowData(sc.obj) %>%
                       as.data.frame() %>%
                       dplyr::mutate(gene = rownames(.), .before = 1) %>%
                       dplyr::left_join(gene_summary, by = "gene") %>%
                       S4Vectors::DataFrame()
    SingleCellExperiment::rowData(sc.obj) <- gene_summary_s4
  } else if (inherits(sc.obj, "Seurat")) {
    orig_metadata <- sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data
    if (ncol(orig_metadata) > 0) {
      new_metadata <- dplyr::mutate(orig_metadata,
                                    gene = rownames(orig_metadata),
                                    .before = 1) %>%
                      dplyr::left_join(gene_summary, by = "gene")
      sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data <- new_metadata
    } else {
      gene_summary <- gene_summary[rownames(sc.obj), ]
      sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data <- gene_summary
    }
  }
}
