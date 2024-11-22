#' Identify highly variable genes in a Bayesian manner.
#'
#' @name findVariableFeaturesBayes
#' @author Jack R. Leary
#' @description This function implements HVG estimation using Bayesian variational inference to approximate the posterior distribution of the mean and dispersion of each gene.
#' @param sc.obj An object of class \code{Seurat} or \code{SingleCellExperiment}. Defaults to NULL.
#' @param subject.id A string specifying the metadata column in \code{expr.mat} that contains subject IDs. Defaults to NULL.
#' @param n.cells.subsample An integer specifying the number of cells per-gene to subsample to when performing estimation. Defaults to 500.
#' @param n.chains (Optional) An integer specifying the number of chains used when performing variational inference. Defaults to 4.
#' @param thin.rate (Optional) An integer specifying the thinning rate of the VI algorithm. Defaults to 5. 
#' @param n.cores An integer specifying the number of cores to be used when fitting the Bayesian hierarchical model. Defaults to 4.
#' @param random.seed A double specifying the random seed to be used when fitting the model. Defaults to 312.
#' @param verbose (Optional) A Boolean specifying whether or not verbose model output should be printed to the console. Defaults to FALSE. 
#' @import cmdstanr
#' @import magrittr
#' @importFrom SingleCellExperiment colData rowData
#' @importFrom BiocGenerics counts
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom dplyr mutate select with_groups slice_sample filter summarise arrange desc left_join pull row_number
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect matches all_of everything
#' @importFrom stats quantile
#' @importFrom withr with_output_sink
#' @importFrom brms set_prior brm bf negbinomial
#' @importFrom posterior as_draws_df
#' @importFrom S4Vectors DataFrame
#' @importFrom SeuratObject Version
#' @return Depending on the input, either an object of class \code{Seurat} or \code{SingleCellExperiment} with gene-level statistics added to the appropriate metadata slot.
#' @seealso \code{\link[Seurat]{FindVariableFeatures}}
#' @seealso \code{\link[scran]{modelGeneVar}}
#' @export

findVariableFeaturesBayes <- function(sc.obj = NULL,
                                      subject.id = NULL,
                                      n.cells.subsample = 500L,
                                      n.chains = 4L,
                                      thin.rate = 5L, 
                                      n.cores = 4L,
                                      random.seed = 312, 
                                      verbose = FALSE) {
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
  # convert counts matrix to long data.frame for modeling
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
  # convert from tibble to data.frame & convert count to integer to save space
  expr_df <- as.data.frame(expr_df) %>% 
             dplyr::mutate(count = as.integer(count)) %>% 
             dplyr::select(-cell)
  # create model formula
  if (!is.null(subject.id)) {
    model_formula <- brms::bf(count ~ 1 + (1 | gene) + (1 | subject), shape ~ 1 + (1 | gene))
  } else {
    model_formula <- brms::bf(count ~ 1 + (1 | gene), shape ~ 1 + (1 | gene))
  }
  # set up priors
  priors <- c(brms::set_prior("normal(0, 10)", class = "Intercept", resp = "mu"),
              brms::set_prior("student_t(3, 0, 10)", class = "sd", resp = "mu"),
              brms::set_prior("normal(0, 5)", class = "Intercept", resp = "shape"),
              brms::set_prior("student_t(3, 0, 10)", class = "sd", resp = "shape"))
  # fit negative-binomial hierarchical bayesian model via variational inference
  if (verbose) {
    brms_fit <- brms::brm(model_formula,
                          data = expr_df,
                          family = brms::negbinomial(link = "log", link_shape = "log"),
                          chains = n.chains,
                          iter = 1000,
                          warmup = 250,
                          thin = thin.rate, 
                          cores = n.cores,
                          silent = 2,
                          backend = "cmdstanr",
                          algorithm = "meanfield",
                          seed = random.seed)
  } else {
    withr::with_output_sink(tempfile(), {
      brms_fit <- brms::brm(model_formula,
                            data = expr_df,
                            family = brms::negbinomial(link = "log", link_shape = "log"),
                            chains = n.chains,
                            iter = 1000,
                            warmup = 250,
                            thin = thin.rate, 
                            cores = n.cores,
                            silent = 2,
                            backend = "cmdstanr",
                            algorithm = "meanfield",
                            seed = random.seed)
    })
  }
  # draw samples from approximate posterior
  posterior_samples <- as.data.frame(posterior::as_draws_df(brms_fit))
  # estimate posterior gene means
  mu_intercept <- dplyr::pull(posterior_samples, b_Intercept)
  mu_random_effects <- dplyr::select(posterior_samples, tidyselect::matches("r_gene\\[.*Intercept")) %>%
                       dplyr::mutate(intercept = mu_intercept)
  mu_samples_long <- tidyr::pivot_longer(mu_random_effects,
                                         cols = !intercept,
                                         names_to = "gene",
                                         values_to = "mu_re") %>%
                     as.data.frame() %>%
                     dplyr::mutate(gene = gsub(",Intercept\\]", "", gsub("r_gene\\[", "", gene)),
                                   mu = exp(intercept + mu_re), 
                                   sample = dplyr::row_number())
  mu_summary <- dplyr::with_groups(mu_samples_long, 
                                   gene,
                                   dplyr::summarise,
                                   mu_mean = mean(mu),
                                   mu_var = var(mu),
                                   mu_ci_ll = stats::quantile(mu, 0.025),
                                   mu_ci_ul = stats::quantile(mu, 0.975))
  # estimate posterior gene dispersions
  theta_intercept <- dplyr::pull(posterior_samples, b_shape_Intercept)
  theta_random_effects <- dplyr::select(posterior_samples, tidyselect::matches("r_gene__shape\\[.*Intercept")) %>%
                          dplyr::mutate(intercept = theta_intercept)
  theta_samples_long <- tidyr::pivot_longer(theta_random_effects,
                                            cols = !intercept,
                                            names_to = "gene",
                                            values_to = "theta_re") %>%
                        as.data.frame() %>% 
                        dplyr::mutate(gene = gsub(",Intercept\\]", "", gsub("r_gene__shape\\[", "", gene)),
                                      theta = 1 / exp(intercept + theta_re), 
                                      sample = dplyr::row_number())
  theta_summary <- dplyr::with_groups(theta_samples_long, 
                                      gene,
                                      dplyr::summarise,
                                      theta_mean = mean(theta),
                                      theta_var = var(theta),
                                      theta_ci_ll = stats::quantile(theta, 0.025),
                                      theta_ci_ul = stats::quantile(theta, 0.975))
  # estimate posterior variance 
  var_samples_long <- dplyr::inner_join(mu_samples_long, 
                                        theta_samples_long, 
                                        by = c("gene", "sample")) %>% 
                      dplyr::rowwise() %>% 
                      dplyr::mutate(sigma2 = mu * (1 + mu / theta)) %>% 
                      dplyr::ungroup()
  var_summary <- dplyr::with_groups(var_samples_long, 
                                    gene,
                                    dplyr::summarise,
                                    sigma2_mean = mean(sigma2),
                                    sigma2_var = var(sigma2),
                                    sigma2_ci_ll = stats::quantile(sigma2, 0.025),
                                    sigma2_ci_ul = stats::quantile(sigma2, 0.975))
  # coerce summaries to a single data.frame 
  gene_summary <- dplyr::inner_join(mu_summary, theta_summary, by = "gene") %>%
                  dplyr::inner_join(var_summary, by = "gene") %>% 
                  magrittr::set_rownames(.$gene)
  # add gene-level estimates to object metadata
  if (inherits(sc.obj, "SingleCellExperiment")) {
    gene_summary_s4 <- SingleCellExperiment::rowData(sc.obj) %>%
                       as.data.frame() %>%
                       dplyr::mutate(gene = rownames(.), .before = 1) %>%
                       dplyr::left_join(gene_summary, by = "gene") %>%
                       S4Vectors::DataFrame()
    SingleCellExperiment::rowData(sc.obj) <- gene_summary_s4
  } else if (inherits(sc.obj, "Seurat")) {
    if (substr(SeuratObject::Version(sc.obj), 1, 1) == "5") {
      orig_metadata <- sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data
    } else {
      orig_metadata <- sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.features
    }
    if (ncol(orig_metadata) > 0) {
      new_metadata <- dplyr::mutate(orig_metadata,
                                    gene = rownames(sc.obj),
                                    .before = 1) %>%
                      dplyr::left_join(gene_summary, by = "gene")
      if (substr(SeuratObject::Version(sc.obj), 1, 1) == "5") {
        sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data <- new_metadata
      } else {
        sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.features <- new_metadata
      }
    } else {
      gene_summary <- gene_summary[rownames(sc.obj), ]
      if (substr(SeuratObject::Version(sc.obj), 1, 1) == "5") {
        sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data <- gene_summary
      } else {
        sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.features <- gene_summary
      }
    }
  }
  return(sc.obj)
}
