# Pairwise protein-protein co-elution feature calculation for MAPLE

# This module calculates pairwise co-elution features between proteins within
# the same sample. These features are used for protein complex and interaction
# analysis.

# Input profiles are assumed to have already been processed by MAPLE
# (normalised, smoothed, calibrated and fraction-aligned). Therefore,
# no additional clean_profiles() step is applied and only a single
# Pearson-correlation-derived distance feature is calculated.

# This is distinct from inter-sample feature calculation:
#   pairwise_features.R      = Protein A vs Protein B within one sample
#   intersample_features.R   = Protein X in Sample A vs Protein X in Sample B

# Required packages:
#   purrr
#   Hmisc
#   lsa
#   progress
library(dplyr)
library(purrr)
library(Hmisc)
library(lsa)
library(progress)


# Calculate co-apex distance between protein pairs using Gaussian components.
#
# For each protein pair, the co-apex distance is defined as the minimum
# Euclidean distance between Gaussian components in the two-dimensional space
# defined by peak centre (mu) and peak width (sigma).
calculate_pairwise_co_apex_maple <- function(
    gaussians,
    proteins = NULL,
    show_progress = TRUE
) {
  
  if (is.null(proteins)) {
    proteins <- names(gaussians)
  }
  
  gaussian_names <- intersect(names(gaussians), proteins)
  n_proteins <- length(proteins)
  
  gaussian_centers <- purrr::map(gaussians, c("coefs", "mu"))
  gaussian_sigmas <- purrr::map(gaussians, c("coefs", "sigma"))
  
  valid_gaussians <- names(gaussians)[
    lengths(gaussian_centers) > 0 &
      lengths(gaussian_sigmas) > 0
  ]
  
  gaussian_centers <- gaussian_centers[valid_gaussians]
  gaussian_sigmas <- gaussian_sigmas[valid_gaussians]
  gaussian_names_valid <- names(gaussian_centers)
  
  if (length(gaussian_names_valid) == 0) {
    return(matrix(
      NA_real_,
      nrow = n_proteins,
      ncol = n_proteins,
      dimnames = list(proteins, proteins)
    ))
  }
  
  gaussian_matrix <- cbind(
    unlist(gaussian_centers),
    unlist(gaussian_sigmas)
  )
  
  CA <- as.matrix(dist(gaussian_matrix))
  
  gaussian_indices <- rep(
    gaussian_names_valid,
    lengths(gaussian_centers)
  )
  
  co_apex <- matrix(
    NA_real_,
    nrow = n_proteins,
    ncol = n_proteins,
    dimnames = list(proteins, proteins)
  )
  
  proteins_with_gaussians <- intersect(gaussian_names_valid, proteins)
  
  if (show_progress) {
    pb <- progress::progress_bar$new(
      format = "Calculating co-apex distances [:bar] :percent eta: :eta",
      total = length(proteins_with_gaussians),
      clear = FALSE,
      width = 60
    )
  }
  
  for (i in seq_along(proteins_with_gaussians)) {
    
    protein_A <- proteins_with_gaussians[i]
    idxs_A <- which(gaussian_indices == protein_A)
    
    for (j in seq(i, length(proteins_with_gaussians))) {
      
      protein_B <- proteins_with_gaussians[j]
      idxs_B <- which(gaussian_indices == protein_B)
      
      dist_value <- min(CA[idxs_A, idxs_B], na.rm = TRUE)
      
      co_apex[protein_A, protein_B] <- dist_value
      co_apex[protein_B, protein_A] <- dist_value
    }
    
    if (show_progress) {
      pb$tick()
    }
  }
  
  return(co_apex)
}


# Calculate pairwise protein-protein co-elution features within one sample.
calculate_pairwise_features_maple <- function(
    profile_matrix,
    gaussians,
    min_pairs = 0,
    pearson_R= TRUE,
    pearson_P = TRUE,
    cosine_correlation = TRUE,
    euclidean_distance = TRUE,
    co_peak = TRUE,
    co_apex = TRUE,
    n_pairs = FALSE,
    max_euclidean_quantile = 0.9,
    show_progress = TRUE
) {
  
  if (is(profile_matrix, "MSnSet")) {
    profile_matrix <- exprs(profile_matrix)
  }
  
  if (min_pairs < 3) {
    min_pairs <- 3
  }
  profile_data <- profile_matrix
  proteins <- rownames(profile_data)
  n_proteins <- length(proteins)
  
  if (is.null(proteins) || n_proteins < 2) {
    stop("profile_matrix must contain at least two proteins with rownames.")
  }
  
  feature_matrices <- list()
  
  pairs <- crossprod(t(!is.na(profile_matrix)))
  
  if (pearson_R) {
    cor_R <- suppressWarnings(
      1 - cor(
        t(profile_matrix),
        use = "pairwise.complete.obs"
      )
    )
    feature_matrices[["cor_R"]] <- cor_R
    message("pearson_R_finished")
  }
  
  if (pearson_P) {
    cor_P <- suppressWarnings(
      Hmisc::rcorr(t(profile_matrix))$P
    )
    feature_matrices[["cor_P"]] <- cor_P
    message("pearson_P finished")
  }
  
  if (cosine_correlation) {
    profile_flip <- data.frame(t(profile_matrix))
    profile_for_cosine <- as.matrix(profile_flip)
    cosine_sim <- lsa::cosine(profile_for_cosine)
    cosine_distance <- 1 - cosine_sim
    feature_matrices[["cosine_distance"]] <- cosine_distance
    message("cosine_distance finished")
  }
  
  if (euclidean_distance) {
    eucl <- as.matrix(
      dist(profile_data, method = "euclidean")
    )
    feature_matrices[["euclidean_distance"]] <- eucl
    message("euclidean_distance finished")
  }
  
  if (co_peak) {
    maxes <- apply(profile_data, 1, which.max)
    co_peak_mat <- as.matrix(dist(maxes))
    feature_matrices[["co_peak"]] <- co_peak_mat
    message("co_peak finished")
  }
  
  if (co_apex) {
    CA <- calculate_pairwise_co_apex_maple(
      gaussians = gaussians,
      proteins = proteins,
      show_progress = show_progress
    )
    feature_matrices[["co_apex"]] <- CA
    message("co_apex finished")
  }
  
  feature_matrices[["n_pairs"]] <- pairs
  
  if (
    !all(purrr::map_dbl(feature_matrices, ~ base::nrow(.x)) == n_proteins) ||
    !all(purrr::map_dbl(feature_matrices, ~ base::ncol(.x)) == n_proteins)
  ) {
    stop("At least one feature matrix did not have the correct dimensions.")
  }
  
  if (length(feature_matrices) == 0) {
    stop("No features were calculated.")
  }
  
  first_matrix <- feature_matrices[[1]]
  tri <- upper.tri(first_matrix)
  idxs <- which(tri, arr.ind = TRUE)
  
  dat <- data.frame(
    protein_A = rownames(first_matrix)[idxs[, 1]],
    protein_B = rownames(first_matrix)[idxs[, 2]],
    stringsAsFactors = FALSE
  )
  
  dat <- cbind(
    dat,
    purrr::map(feature_matrices, ~ .[tri])
  )
  
  if ("euclidean_distance" %in% colnames(dat)) {
    vec <- dat$euclidean_distance
    threshold <- quantile(
      vec,
      probs = max_euclidean_quantile,
      na.rm = TRUE
    )
    dat$euclidean_distance[vec >= threshold] <- threshold
  }
  
  dat <- dplyr::filter(dat, n_pairs >= min_pairs)
  
  if (!n_pairs) {
    dat$n_pairs <- NULL
  }
  
  return(dat)
}
