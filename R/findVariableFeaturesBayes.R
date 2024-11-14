#' Identify highly variable genes in a Bayesian manner.
#'
#' @name findVariableFeaturesBayes
#' @author Jack R. Leary
#' @description This function implements HVG estimation using Bayesian variational inference to approximate the posterior distribution of the mean, variance, and dispersion of each gene.
#' @param sc.obj An object of class \code{Seurat} or \code{SingleCellExperiment}. Defaults to NULL.
#' @param subject.id A string specifying the metadata column in \code{expr.mat} that contains subject IDs. Defaults to NULL.
#' @param n.cells.subsample An integer specifying the number of cells per-gene to subsample to when performing estimation. Defaults to 500.
#' @param n.chains (Optional) An integer specifying the number of chains used when performing variational inference. Defaults to 4. 
#' @param n.cores An integer specifying the number of cores to be used when fitting the Bayesian hierarchical model. Defaults to 4.
#' @import cmdstanr
#' @import magrittr
#' @importFrom SingleCellExperiment colData rowData
#' @importFrom BiocGenerics counts
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom dplyr mutate select with_groups slice_sample filter arrange desc left_join
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect matches all_of everything  
#' @importFrom stats as.formula quantile
#' @importFrom withr with_output_sink
#' @importFrom brms set_prior brm bf negbinomial
#' @importFrom S4Vectors DataFrame
#' @return Depending on the input, either an object of class \code{Seurat} or \code{SingleCellExperiment} with HVG metadata added.
#' @seealso \code{\link[Seurat]{FindVariableFeatures}}
#' @seealso \code{\link[scran]{modelGeneVar}}
#' @export

findVariableFeaturesBayes <- function(sc.obj = NULL,
                                      subject.id = NULL,
                                      n.cells.subsample = 500L,
                                      n.chains = 4L, 
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
  # create model formula
  if (!is.null(subject.id)) {
    model_formula <- stats::as.formula("count ~ 1 + f(gene, model = 'iid') + f(subject, model = 'iid')")
  } else {
    model_formula <- stats::as.formula("count ~ 1 + f(gene, model = 'iid')")
  }
  # set up priors 
  priors <- c(brms::set_prior("normal(0, 5)", class = "Intercept", resp = "mu"),
              brms::set_prior("student_t(3, 0, 10)", class = "sd", resp = "mu"), 
              brms::set_prior("normal(0, 5)", class = "Intercept", resp = "shape"),
              brms::set_prior("student_t(3, 0, 10)", class = "sd", resp = "shape"))
  # fit negative-binomial hierarchical bayesian model via variational inference  
  withr::with_output_sink(tempfile(), {
    brms_fit <- brms::brm(brms::bf(count ~ 1 + (1 | gene), shape ~ gene),
                          data = expr_df2,
                          family = brms::negbinomial(link = "log", link_shape = "log"),
                          chains = n.chains,  
                          iter = 1000,
                          warmup = 250, 
                          cores = n.cores, 
                          silent = 0, 
                          backend = "cmdstanr", 
                          algorithm = "meanfield", 
                          seed = random.seed)
  })
  # sample from approximate posterior distribution for gene mean
  mu_random_effects <- dplyr::select(posterior_samples, tidyselect::matches("r_gene.*Intercept"))
  mu_intercept <- posterior_samples$b_Intercept
  mu_long <- tidyr::pivot_longer(mu_random_effects, 
                                 cols = tidyselect::everything(), 
                                 names_to = "gene", 
                                 values_to = "mu_re")
  # sample from approximate posterior distribution for overdispersion parameter
  shape_columns <- grep("^b_shape_gene.*", colnames(posterior_samples), value = TRUE)
  shape_summary <- dplyr::select(posterior_samples, tidyselect::all_of(shape_columns)) %>%
                   tidyr::pivot_longer(cols = tidyselect::everything(), 
                                       names_to = "gene", 
                                       values_to = "shape_log") %>%
                   dplyr::mutate(gene = gsub("b_shape_gene", "", gene),
                                 shape = exp(shape_log)) %>%
                   dplyr::with_groups(gene, 
                                      dplyr::summarise, 
                                      theta_mean = mean(shape),
                                      theta_sd = sd(shape),
                                      theta_lower = stats::quantile(shape, 0.025),
                                      theta_upper = stats::quantile(shape, 0.975)) %>%
                   dplyr::arrange(dplyr::desc(theta_mean))
  # estimate variance based on formula for negative-binomial variance
  
  # generate central tendency estimates and credible intervals for each parameter
  
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
                                    gene = rownames(sc.obj),
                                    .before = 1) %>%
                      dplyr::left_join(gene_summary, by = "gene")
      sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data <- new_metadata
    } else {
      gene_summary <- gene_summary[rownames(sc.obj), ]
      sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data <- gene_summary
    }
  }
}
