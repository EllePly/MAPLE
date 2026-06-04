# Calibration functions for MAPLE
#
# These functions estimate and apply fraction-axis calibration between
# SEC protein elution profiles using shared single-hit Gaussian peak centres.
#
# Core logic:
#   1. Extract shared single-hit Gaussian peak centres between two samples.
#   2. Fit a linear calibration model:
#        reference_mu = slope * query_mu + intercept
#   3. Apply the calibration model to remap the query profile onto the
#      reference fraction axis using interpolation.
#
# Required packages:
#   dplyr
#   ggplot2

library(dplyr)
library(ggplot2)

# Estimate calibration parameters for one pair of samples.
#
# This low-level function compares shared single-hit Gaussian peak centres
# between one reference sample and one query sample, then fits a linear model:
#
#   reference_mu = slope * query_mu + intercept
#
# The resulting slope and intercept are used to map the query sample onto the
# reference fraction axis.
estimate_calibration_parameters_maple <- function(
    gaussians_singlehit_list,
    reference_sample,
    query_sample,
    max_mu_diff = 5,
    weighted_lm = FALSE,
    weight_bin_edges = seq(35, 75, by = 5),
    return_plot = FALSE
) {
  
  reference_single <- gaussians_singlehit_list[[reference_sample]]
  query_single <- gaussians_singlehit_list[[query_sample]]
  
  if (is.null(reference_single)) {
    stop("reference_sample not found in gaussians_singlehit_list: ", reference_sample)
  }
  
  if (is.null(query_single)) {
    stop("query_sample not found in gaussians_singlehit_list: ", query_sample)
  }
  
  common_proteins <- Reduce(
    intersect,
    list(names(reference_single), names(query_single))
  )
  
  if (length(common_proteins) == 0) {
    stop("No common proteins found between ", reference_sample, " and ", query_sample)
  }
  
  mu_df <- data.frame(
    GeneName = common_proteins,
    reference_mu = vapply(common_proteins, function(gene) {
      if (!is.null(reference_single[[gene]]$coefs$mu)) {
        reference_single[[gene]]$coefs$mu
      } else {
        NA_real_
      }
    }, numeric(1)),
    query_mu = vapply(common_proteins, function(gene) {
      if (!is.null(query_single[[gene]]$coefs$mu)) {
        query_single[[gene]]$coefs$mu
      } else {
        NA_real_
      }
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  
  mu_df <- mu_df %>%
    dplyr::filter(
      is.finite(reference_mu),
      is.finite(query_mu)
    ) %>%
    dplyr::mutate(
      mu_diff = reference_mu - query_mu
    ) %>%
    dplyr::filter(abs(mu_diff) <= max_mu_diff)
  
  if (nrow(mu_df) < 3) {
    stop(
      "Too few proteins after filtering to estimate calibration model: ",
      reference_sample, " vs ", query_sample
    )
  }
  
  if (isTRUE(weighted_lm)) {
    
    bin_labels <- paste(
      head(weight_bin_edges, -1),
      tail(weight_bin_edges, -1),
      sep = "-"
    )
    
    mu_df$mu_bin <- cut(
      mu_df$query_mu,
      breaks = weight_bin_edges,
      labels = bin_labels,
      include.lowest = TRUE
    )
    
    bin_counts <- table(mu_df$mu_bin)
    bin_weights <- 1 / (1 + bin_counts)
    
    mu_df$weights <- as.numeric(bin_weights[as.character(mu_df$mu_bin)])
    
    lm_model <- lm(
      reference_mu ~ query_mu,
      data = mu_df,
      weights = weights
    )
    
  } else {
    
    mu_df$weights <- NA_real_
    
    lm_model <- lm(
      reference_mu ~ query_mu,
      data = mu_df
    )
  }
  
  slope <- unname(coef(lm_model)[["query_mu"]])
  intercept <- unname(coef(lm_model)[["(Intercept)"]])
  r_squared <- summary(lm_model)$r.squared
  
  comparison_name <- paste0(reference_sample, "_", query_sample, "_lm")
  
  parameter_row <- data.frame(
    Comparison = comparison_name,
    Reference = reference_sample,
    Query = query_sample,
    Slope = slope,
    Intercept = intercept,
    R2 = r_squared,
    N_proteins = nrow(mu_df),
    Weighted = weighted_lm,
    stringsAsFactors = FALSE
  )
  
  p_scatter <- NULL
  
  if (isTRUE(return_plot)) {
    
    p_scatter <- ggplot(mu_df, aes(x = query_mu, y = reference_mu)) +
      geom_point(
        aes(size = if (weighted_lm) weights else NULL),
        alpha = 0.7
      ) +
      geom_abline(
        intercept = intercept,
        slope = slope,
        linewidth = 1
      ) +
      geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed"
      ) +
      theme_minimal() +
      labs(
        title = paste0(reference_sample, " vs ", query_sample),
        x = paste0(query_sample, " Gaussian peak centre"),
        y = paste0(reference_sample, " Gaussian peak centre")
      ) +
      guides(size = "none")
  }
  
  return(list(
    parameters = parameter_row,
    model = lm_model,
    mu_table = mu_df,
    plot = p_scatter
  ))
}

# Estimate calibration parameters from all query samples to one reference sample.
#
# This wrapper applies estimate_calibration_parameters_maple() to every sample
# in gaussians_singlehit_list except the reference sample. It returns a combined
# calibration parameter table, together with optional fitted models, diagnostic
# plots and Gaussian peak-centre tables.
estimate_all_to_reference_calibration_maple <- function(
    gaussians_singlehit_list,
    reference_sample,
    samples = NULL,
    max_mu_diff = 5,
    weighted_lm = FALSE,
    weight_bin_edges = seq(35, 75, by = 5),
    return_models = TRUE,
    return_plots = FALSE
) {
  
  if (is.null(samples)) {
    samples <- setdiff(names(gaussians_singlehit_list), reference_sample)
  } else {
    samples <- setdiff(samples, reference_sample)
  }
  
  results <- lapply(samples, function(query_sample) {
    
    message(
      "Estimating calibration parameters: ",
      reference_sample, " vs ", query_sample
    )
    
    estimate_calibration_parameters_maple(
      gaussians_singlehit_list = gaussians_singlehit_list,
      reference_sample = reference_sample,
      query_sample = query_sample,
      max_mu_diff = max_mu_diff,
      weighted_lm = weighted_lm,
      weight_bin_edges = weight_bin_edges,
      return_plot = return_plots
    )
  })
  
  names(results) <- samples
  
  calparam <- dplyr::bind_rows(
    lapply(results, function(x) x$parameters)
  )
  
  out <- list(
    calparam = calparam
  )
  
  if (isTRUE(return_models)) {
    out$models <- lapply(results, function(x) x$model)
  }
  
  if (isTRUE(return_plots)) {
    out$plots <- lapply(results, function(x) x$plot)
  }
  
  out$mu_tables <- lapply(results, function(x) x$mu_table)
  
  return(out)
}


calibrate_profile_list_maple <- function(
    data_list,
    calparam,
    reference_sample,
    group = NULL,
    samples = NULL
) {
  
  if (is.null(samples)) {
    
    if (is.null(group) || group == "" || group == "ALL") {
      samples <- names(data_list)
    } else {
      samples <- names(data_list)[grep(paste0("^", group), names(data_list))]
    }
  }
  
  samples_to_calibrate <- setdiff(samples, reference_sample)
  
  calibrated_data_list <- list()
  
  for (sample in samples_to_calibrate) {
    
    comparison_name <- paste0(reference_sample, "_", sample, "_lm")
    
    row_match <- calparam[calparam$Comparison == comparison_name, ]
    
    if (nrow(row_match) == 0) {
      warning(
        "Calibration parameter not found for ",
        comparison_name,
        ". Skipping this sample."
      )
      next
    }
    
    slope <- row_match$Slope[1]
    intercept <- row_match$Intercept[1]
    
    data <- data_list[[sample]]
    
    if (is.null(data)) {
      warning("Sample not found in data_list: ", sample)
      next
    }
    
    calibrated_data <- data.frame(
      matrix(nrow = nrow(data), ncol = ncol(data))
    )
    
    rownames(calibrated_data) <- rownames(data)
    colnames(calibrated_data) <- colnames(data)
    
    original_fractions <- as.numeric(gsub("F", "", colnames(data)))
    
    if (any(is.na(original_fractions))) {
      stop(
        "Fraction columns must be named as F followed by numbers, e.g. F38, F39."
      )
    }
    
    calibrated_fractions <- slope * original_fractions + intercept
    
    for (gene in rownames(data)) {
      
      interp_function <- approxfun(
        calibrated_fractions,
        as.numeric(data[gene, ]),
        rule = 2
      )
      
      calibrated_data[gene, ] <- interp_function(original_fractions)
    }
    
    calibrated_data_list[[sample]] <- calibrated_data
    
    message(
      "Calibration completed: ",
      sample,
      " mapped to reference sample ",
      reference_sample,
      " using model ",
      comparison_name
    )
  }
  
  return(calibrated_data_list)
}


add_reference_sample_maple <- function(
    calibrated_data_list,
    data_list,
    reference_sample
) {
  
  if (!reference_sample %in% names(data_list)) {
    stop("reference_sample not found in data_list: ", reference_sample)
  }
  
  calibrated_data_list[[reference_sample]] <- data_list[[reference_sample]]
  
  return(calibrated_data_list)
}
