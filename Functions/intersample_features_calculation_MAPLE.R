# Inter-sample protein elution profile feature analysis for MAPLE
#
# This module compares protein elution profiles between samples.
#
# Core workflow:
#   1. Fill missing SEC fraction columns across samples.
#   2. Calculate inter-sample profile features for each matched protein.
#   3. Run all pairwise sample comparisons.
#   4. Convert similarity features into distance-like features.
#   5. Parse sample-pair metadata into within-condition and between-condition comparisons.
#   6. Calculate protein-specific z-scores using within-condition variability as background.
#
# Required packages:
#   dplyr

library(dplyr)


# Fill missing fraction columns in a list of protein profile matrices/data frames.
#
# Missing internal fractions are interpolated from available neighbouring
# fractions. Remaining missing values, including boundary values, can optionally
# be replaced with zero.
fill_missing_fractions_maple <- function(
    data_list,
    target_fractions = paste0("F", 37:75),
    fill_boundary_na_with_zero = TRUE
) {
  
  target_nums <- as.numeric(gsub("^F", "", target_fractions))
  filled_list <- list()
  
  for (sample in names(data_list)) {
    
    mat <- data_list[[sample]]
    
    original_cols <- colnames(mat)
    original_nums <- as.numeric(gsub("^F", "", original_cols))
    
    if (any(is.na(original_nums))) {
      stop(
        "All profile columns must be named as F followed by numbers, e.g. F38, F39. ",
        "Problem found in sample: ", sample
      )
    }
    
    new_mat <- matrix(
      NA_real_,
      nrow = nrow(mat),
      ncol = length(target_fractions)
    )
    
    rownames(new_mat) <- rownames(mat)
    colnames(new_mat) <- target_fractions
    
    for (gene in rownames(mat)) {
      
      y <- as.numeric(mat[gene, ])
      valid <- is.finite(original_nums) & is.finite(y)
      
      if (sum(valid) < 2) {
        next
      }
      
      interpolation_function <- approxfun(
        x = original_nums[valid],
        y = y[valid],
        rule = 1
      )
      
      new_mat[gene, ] <- interpolation_function(target_nums)
    }
    
    if (isTRUE(fill_boundary_na_with_zero)) {
      new_mat[is.na(new_mat)] <- 0
    }
    
    filled_list[[sample]] <- as.data.frame(new_mat)
    
    message("Filled missing fractions for: ", sample)
  }
  
  return(filled_list)
}


# Calculate cosine similarity between two protein elution profiles.
cosine_similarity_feature_maple <- function(x, y) {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(NA_real_)
  }
  
  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  
  if (!is.finite(denom) || denom == 0) {
    return(NA_real_)
  }
  
  sum(x * y) / denom
}


# Calculate Euclidean distance between two protein elution profiles.
euclidean_distance_feature_maple <- function(x, y) {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(NA_real_)
  }
  
  sqrt(sum((x - y)^2))
}


# Calculate co-peak distance between two protein elution profiles.
#
# Co-peak distance is the absolute difference between the maximum-intensity
# fraction positions in the two profiles.
co_peak_distance_feature_maple <- function(x, y) {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 3) {
    return(NA_real_)
  }
  
  abs(which.max(x) - which.max(y))
}


# Calculate co-apex shift using Gaussian peak centres and widths.
#
# For proteins with Gaussian fitting results in both samples, the co-apex shift
# is calculated as the minimum Euclidean distance between Gaussian components
# in the two-dimensional space defined by peak centre (mu) and width (sigma).
co_apex_shift_feature_maple <- function(
    protein,
    gaussians_A,
    gaussians_B
) {
  
  if (is.null(gaussians_A) || is.null(gaussians_B)) {
    return(NA_real_)
  }
  
  if (is.null(gaussians_A[[protein]]) || is.null(gaussians_B[[protein]])) {
    return(NA_real_)
  }
  
  entry_A <- gaussians_A[[protein]]
  entry_B <- gaussians_B[[protein]]
  
  if (
    is.null(entry_A$coefs) ||
    is.null(entry_B$coefs) ||
    is.null(entry_A$coefs$mu) ||
    is.null(entry_B$coefs$mu) ||
    is.null(entry_A$coefs$sigma) ||
    is.null(entry_B$coefs$sigma)
  ) {
    return(NA_real_)
  }
  
  muA <- entry_A$coefs$mu
  sigmaA <- entry_A$coefs$sigma
  
  muB <- entry_B$coefs$mu
  sigmaB <- entry_B$coefs$sigma
  
  matA <- cbind(muA, sigmaA)
  matB <- cbind(muB, sigmaB)
  
  if (nrow(matA) == 0 || nrow(matB) == 0) {
    return(NA_real_)
  }
  
  d <- as.matrix(dist(rbind(matA, matB)))
  
  nA <- nrow(matA)
  nB <- nrow(matB)
  
  dAB <- d[1:nA, (nA + 1):(nA + nB), drop = FALSE]
  
  min(dAB, na.rm = TRUE)
}


# Calculate inter-sample profile features for matched proteins.
#
# This function compares the elution profile of the same protein between two
# samples. It can optionally use Gaussian fitting results to calculate co-apex
# shifts.
calculate_features_intersample_maple <- function(
    profile_matrix_A,
    profile_matrix_B,
    gaussians_A = NULL,
    gaussians_B = NULL,
    pearson_R_raw = TRUE,
    pearson_R_cleaned = TRUE,
    pearson_P = TRUE,
    cosine_correlation = TRUE,
    euclidean_distance = TRUE,
    co_peak = TRUE,
    co_apex = TRUE
) {
  
  common_proteins <- intersect(
    rownames(profile_matrix_A),
    rownames(profile_matrix_B)
  )
  
  if (length(common_proteins) == 0) {
    stop("No common proteins found between profile_matrix_A and profile_matrix_B.")
  }
  
  common_fractions <- intersect(
    colnames(profile_matrix_A),
    colnames(profile_matrix_B)
  )
  
  common_fractions <- common_fractions[
    order(as.numeric(gsub("^F", "", common_fractions)))
  ]
  
  if (length(common_fractions) < 3) {
    stop("Fewer than 3 common fraction columns were found.")
  }
  
  dataA <- profile_matrix_A[common_proteins, common_fractions, drop = FALSE]
  dataB <- profile_matrix_B[common_proteins, common_fractions, drop = FALSE]
  
  result <- data.frame(
    protein = common_proteins,
    stringsAsFactors = FALSE
  )
  
  for (g in common_proteins) {
    
    x <- as.numeric(dataA[g, ])
    y <- as.numeric(dataB[g, ])
    
    ok <- is.finite(x) & is.finite(y)
    x_ok <- x[ok]
    y_ok <- y[ok]
    
    if (pearson_R_raw) {
      result[result$protein == g, "pearson_R_raw"] <- suppressWarnings(
        cor(x, y, use = "pairwise.complete.obs")
      )
    }
    
    if (pearson_R_cleaned) {
      if (length(x_ok) >= 3) {
        result[result$protein == g, "pearson_R_cleaned"] <- suppressWarnings(
          cor(x_ok, y_ok)
        )
      } else {
        result[result$protein == g, "pearson_R_cleaned"] <- NA_real_
      }
    }
    
    if (pearson_P) {
      if (length(x_ok) >= 3) {
        result[result$protein == g, "pearson_P"] <- suppressWarnings(
          cor.test(x_ok, y_ok)$p.value
        )
      } else {
        result[result$protein == g, "pearson_P"] <- NA_real_
      }
    }
    
    if (cosine_correlation) {
      result[result$protein == g, "cosine_sim"] <-
        cosine_similarity_feature_maple(x, y)
    }
    
    if (euclidean_distance) {
      result[result$protein == g, "euclidean"] <-
        euclidean_distance_feature_maple(x, y)
    }
    
    if (co_peak) {
      result[result$protein == g, "co_peak_diff"] <-
        co_peak_distance_feature_maple(x, y)
    }
    
    if (co_apex) {
      result[result$protein == g, "co_apex_shift"] <-
        co_apex_shift_feature_maple(
          protein = g,
          gaussians_A = gaussians_A,
          gaussians_B = gaussians_B
        )
    }
  }
  
  return(result)
}


# Run inter-sample profile feature calculation for all sample pairs.
run_all_intersample_feature_comparisons_maple <- function(
    data_list,
    gaussians_list = NULL,
    sample_names = NULL,
    pearson_R_raw = TRUE,
    pearson_R_cleaned = TRUE,
    pearson_P = TRUE,
    cosine_correlation = TRUE,
    euclidean_distance = TRUE,
    co_peak = TRUE,
    co_apex = TRUE
) {
  
  if (is.null(sample_names)) {
    sample_names <- names(data_list)
  }
  
  combos <- combn(sample_names, 2, simplify = FALSE)
  
  results <- list()
  
  for (combo in combos) {
    
    sampleA <- combo[1]
    sampleB <- combo[2]
    
    message("Running inter-sample comparison: ", sampleA, " vs ", sampleB)
    
    gaussians_A <- NULL
    gaussians_B <- NULL
    
    if (!is.null(gaussians_list)) {
      gaussians_A <- gaussians_list[[sampleA]]
      gaussians_B <- gaussians_list[[sampleB]]
    }
    
    res <- calculate_features_intersample_maple(
      profile_matrix_A = data_list[[sampleA]],
      profile_matrix_B = data_list[[sampleB]],
      gaussians_A = gaussians_A,
      gaussians_B = gaussians_B,
      pearson_R_raw = pearson_R_raw,
      pearson_R_cleaned = pearson_R_cleaned,
      pearson_P = pearson_P,
      cosine_correlation = cosine_correlation,
      euclidean_distance = euclidean_distance,
      co_peak = co_peak,
      co_apex = co_apex
    )
    
    results[[paste0(sampleA, "_vs_", sampleB)]] <- res
  }
  
  return(results)
}


# Convert similarity features into distance-like features.
#
# After conversion, larger values indicate larger profile differences.
# Pearson R and cosine similarity are transformed as 1 - similarity.
convert_similarity_to_distance_maple <- function(
    feature_results,
    pearson_R_raw = TRUE,
    pearson_R_cleaned = TRUE,
    cosine_sim = TRUE
) {
  
  converted_results <- lapply(feature_results, function(df) {
    
    if (isTRUE(pearson_R_raw) && "pearson_R_raw" %in% colnames(df)) {
      df$pearson_R_raw <- 1 - df$pearson_R_raw
    }
    
    if (isTRUE(pearson_R_cleaned) && "pearson_R_cleaned" %in% colnames(df)) {
      df$pearson_R_cleaned <- 1 - df$pearson_R_cleaned
    }
    
    if (isTRUE(cosine_sim) && "cosine_sim" %in% colnames(df)) {
      df$cosine_sim <- 1 - df$cosine_sim
    }
    
    return(df)
  })
  
  return(converted_results)
}


# Parse sample-pair names and classify comparisons as within or between.
#
# Pair names are expected to follow:
#   SampleA_vs_SampleB
#
# By default, biological condition names are inferred by removing trailing
# replicate numbers, e.g. WT1 -> WT.
parse_pair_info_maple <- function(
    pair_names,
    condition_fun = function(s) sub("\\d+$", "", s),
    sort_condition_labels = TRUE
) {
  
  pair_info <- do.call(rbind, lapply(pair_names, function(name) {
    
    parts <- strsplit(name, "_vs_")[[1]]
    
    if (length(parts) != 2) {
      stop("Pair name must contain exactly one '_vs_': ", name)
    }
    
    sampleA <- parts[1]
    sampleB <- parts[2]
    
    groupA <- condition_fun(sampleA)
    groupB <- condition_fun(sampleB)
    
    data.frame(
      pair = name,
      sampleA = sampleA,
      sampleB = sampleB,
      groupA = groupA,
      groupB = groupB,
      stringsAsFactors = FALSE
    )
  }))
  
  pair_info <- pair_info %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      g1 = if (isTRUE(sort_condition_labels)) sort(c(groupA, groupB))[1] else groupA,
      g2 = if (isTRUE(sort_condition_labels)) sort(c(groupA, groupB))[2] else groupB,
      comparison = paste(g1, g2, sep = "_vs_"),
      type = ifelse(groupA == groupB, "within", "between")
    ) %>%
    dplyr::ungroup()
  
  return(pair_info)
}


# Calculate protein-specific z-scores for one feature.
#
# Within-condition sample-pair comparisons are used to define a protein-specific
# background distribution. All within- and between-condition comparisons are
# then converted to z-scores relative to this background.
compute_z_feature_all_maple <- function(
    feature_name,
    dataset,
    pair_info
) {
  
  all_feature_df <- lapply(names(dataset), function(nm) {
    
    df <- dataset[[nm]]
    
    if (!feature_name %in% colnames(df)) {
      stop("Feature not found in dataset[[", nm, "]]: ", feature_name)
    }
    
    df %>%
      dplyr::select(protein, dplyr::all_of(feature_name)) %>%
      dplyr::rename(value = dplyr::all_of(feature_name)) %>%
      dplyr::mutate(pair = nm)
    
  }) %>%
    dplyr::bind_rows()
  
  all_feature_df <- dplyr::left_join(
    all_feature_df,
    pair_info,
    by = "pair"
  )
  
  within_stats <- all_feature_df %>%
    dplyr::filter(type == "within") %>%
    dplyr::group_by(protein) %>%
    dplyr::summarise(
      mean_within = mean(value, na.rm = TRUE),
      sd_within = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  all_z <- all_feature_df %>%
    dplyr::left_join(within_stats, by = "protein") %>%
    dplyr::mutate(
      z = ifelse(
        is.na(sd_within) | sd_within == 0,
        NA_real_,
        (value - mean_within) / sd_within
      )
    )
  
  return(all_z)
}
