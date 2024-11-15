# Changes in v0.0.1

+ Initial package skeleton. 
+ Added main function `findVariableFeaturesBayes()` (still in development right now). 
+ Added helper function `sampleMarginal()` to make sampling from posterior marginal distribution possible. 
+ Added support for adding gene statistics metadata to `SingleCellExperiment` or `Seurat` objects after estimation. 
+ Changed main model backend to `brms` via `cmdstanr` instead of `INLA`, since `INLA` doesn't appear to support per-group estimation of the negative-binomial dispersion parameter.
+ Added function `classifyHVGs()` to add a label to HVGs in `Seurat` or `SingleCellExperiment` object metadata based on several different methods. 
