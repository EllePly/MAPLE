Gaussian fitting is performed at the sample-matrix level using build_gaussians_maple(). The resulting per-sample Gaussian result lists are then passed to select_single_hit_gaussians_maple(), which selects proteins with one dominant Gaussian peak and identifies condition-level common single-hit proteins across replicates.

gaussian_fitting.R
    Gaussian mixture fitting adapted from PrInCE.

single_hit_selection.R
    Select proteins with a single dominant Gaussian peak and identify condition-level common proteins across biological replicates.
