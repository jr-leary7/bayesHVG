# Changes in v0.0.3

+ Updated support for v4 vs. v5 `Seurat` objects. 
+ Exposed choice of variational inference algorithm to the user via argument `VI.algorithm` in `findvariableFeaturesBayes()`. 
+ Updated documentation thoroughly. 
+ Removed QR decomposition on fixed effects as it wasn't necessary. 

# Changes in v0.0.2

+ Added conditional support for legacy `Seurat` v4 objects in addition to the default `Seurat` v5 objects. 
+ Added per-gene estimated posterior variances and dispersions based on NB variance definition, along with credible intervals for each per-gene. 
+ Implemented ability to select HVGs by either estimated dispersion or estimated variance in `classifyHVGs()`. 
+ Changed model fitting process to support within-chain parallelism as long as enough cores are available. 
+ Added function `computeNaiveGeneStatistics()` to estimate (in a Frequentist manner) per-gene mean, variance, and dispersion.
+ Sped up `findvariableFeaturesBayes()` via:
  - Adding GPU acculeration support for OpenCL-compatible devices. 
  - Performing QR decomposition on covariates before fitting. 
  - Added compiler optimization flags for Stan to C++ code conversion. 
  - Set `normalize = FALSE` in call to `brms::brm()`, which increases efficiency of sampling from approximate posterior. 

# Changes in v0.0.1

+ Initial package skeleton. 
+ Added main function `findVariableFeaturesBayes()` (still in development right now). 
+ Added helper function `sampleMarginal()` to make sampling from posterior marginal distribution possible. 
+ Added support for adding gene statistics metadata to `SingleCellExperiment` or `Seurat` objects after estimation. 
+ Changed main model backend to `brms` via `cmdstanr` instead of `INLA`, since `INLA` doesn't appear to support per-group estimation of the Negative-binomial overdispersion parameter.
+ Added function `classifyHVGs()` to add a label to HVGs in `Seurat` or `SingleCellExperiment` object metadata based on several different methods. 
+ Added function `theme_bayesHVG()`, which implements a publication-ready theme for `ggplot2`. 
