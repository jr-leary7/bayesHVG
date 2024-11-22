# Changes ion v0.0.2

+ Added conditional support for `Seurat` v4 objects in addition to `Seurat` v5 objects. 
+ Added per-gene estimated variances based on NB variance definition, along with credible intervals per-gene. 
+ Implemented ability to select HVGs by either estimated dispersion or estimated variance in `classifyHVGs()`. 
+ Changed model fitting process to support within-chain parallelism as long as enough cores are available. 

# Changes in v0.0.1

+ Initial package skeleton. 
+ Added main function `findVariableFeaturesBayes()` (still in development right now). 
+ Added helper function `sampleMarginal()` to make sampling from posterior marginal distribution possible. 
+ Added support for adding gene statistics metadata to `SingleCellExperiment` or `Seurat` objects after estimation. 
+ Changed main model backend to `brms` via `cmdstanr` instead of `INLA`, since `INLA` doesn't appear to support per-group estimation of the negative-binomial dispersion parameter.
+ Added function `classifyHVGs()` to add a label to HVGs in `Seurat` or `SingleCellExperiment` object metadata based on several different methods. 
+ Added function `theme_bayesHVG()`, which implements a publication-ready theme for `ggplot2`. 
