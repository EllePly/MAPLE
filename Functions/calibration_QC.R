# Calibration QC functions for MAPLE
#
# These functions evaluate profile similarity before and after fraction-axis
# calibration using Pearson R2, RMSE and cosine similarity.
#
# Required packages:
#   dplyr
#   tidyr
#   ggplot2
#   gridExtra

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)


# Calculate RMSE between two protein elution profiles.
profile_rmse_maple <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 3) {
    return(NA_real_)
  }
  
  sqrt(mean((x[ok] - y[ok])^2))
}


# Calculate cosine similarity between two protein elution profiles.
cosine_similarity_maple <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(NA_real_)
  }
  
  if (sqrt(sum(x^2)) == 0 || sqrt(sum(y^2)) == 0) {
    return(NA_real_)
  }
  
  sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}


# Compare profile similarity before and after calibration for sample pairs.
#
# For each sample pair, common proteins are matched between the raw and
# calibrated profile lists. Pearson R2, RMSE and cosine similarity are then
# calculated before and after calibration.
compare_profiles_before_after_calibration_maple <- function(
    sample_pairs,
    rawdata_list,
    caldata_list
) {
  
  pairwise_qc_list <- list()
  
  for (pair in sample_pairs) {
    
    sampleA <- pair[1]
    sampleB <- pair[2]
    
    if (!sampleA %in% names(rawdata_list) || !sampleB %in% names(rawdata_list)) {
      warning("Sample pair not found in rawdata_list: ", sampleA, " vs ", sampleB)
      next
    }
    
    if (!sampleA %in% names(caldata_list) || !sampleB %in% names(caldata_list)) {
      warning("Sample pair not found in caldata_list: ", sampleA, " vs ", sampleB)
      next
    }
    
    common_genes <- Reduce(intersect, list(
      rownames(rawdata_list[[sampleA]]),
      rownames(rawdata_list[[sampleB]]),
      rownames(caldata_list[[sampleA]]),
      rownames(caldata_list[[sampleB]])
    ))
    
    if (length(common_genes) == 0) {
      warning("No common genes found for pair: ", sampleA, " vs ", sampleB)
      next
    }
    
    # Raw profiles
    mat_raw_A <- rawdata_list[[sampleA]][common_genes, , drop = FALSE] %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
      as.matrix()
    
    mat_raw_B <- rawdata_list[[sampleB]][common_genes, , drop = FALSE] %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
      as.matrix()
    
    common_raw_cols <- intersect(colnames(mat_raw_A), colnames(mat_raw_B))
    common_raw_cols <- common_raw_cols[
      order(as.numeric(gsub("F", "", common_raw_cols)))
    ]
    
    mat_raw_A <- mat_raw_A[, common_raw_cols, drop = FALSE]
    mat_raw_B <- mat_raw_B[, common_raw_cols, drop = FALSE]
    
    # Calibrated profiles
    mat_calib_A <- caldata_list[[sampleA]][common_genes, , drop = FALSE] %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
      as.matrix()
    
    mat_calib_B <- caldata_list[[sampleB]][common_genes, , drop = FALSE] %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
      as.matrix()
    
    common_calib_cols <- intersect(colnames(mat_calib_A), colnames(mat_calib_B))
    common_calib_cols <- common_calib_cols[
      order(as.numeric(gsub("F", "", common_calib_cols)))
    ]
    
    mat_calib_A <- mat_calib_A[, common_calib_cols, drop = FALSE]
    mat_calib_B <- mat_calib_B[, common_calib_cols, drop = FALSE]
    
    # Profile similarity metrics
    r_before <- sapply(seq_along(common_genes), function(i) {
      cor(
        as.numeric(mat_raw_A[i, ]),
        as.numeric(mat_raw_B[i, ]),
        use = "pairwise.complete.obs"
      )^2
    })
    
    r_after <- sapply(seq_along(common_genes), function(i) {
      cor(
        as.numeric(mat_calib_A[i, ]),
        as.numeric(mat_calib_B[i, ]),
        use = "pairwise.complete.obs"
      )^2
    })
    
    rmse_before <- sapply(seq_along(common_genes), function(i) {
      profile_rmse_maple(mat_raw_A[i, ], mat_raw_B[i, ])
    })
    
    rmse_after <- sapply(seq_along(common_genes), function(i) {
      profile_rmse_maple(mat_calib_A[i, ], mat_calib_B[i, ])
    })
    
    cosine_before <- sapply(seq_along(common_genes), function(i) {
      cosine_similarity_maple(mat_raw_A[i, ], mat_raw_B[i, ])
    })
    
    cosine_after <- sapply(seq_along(common_genes), function(i) {
      cosine_similarity_maple(mat_calib_A[i, ], mat_calib_B[i, ])
    })
    
    pairwise_df <- data.frame(
      Gene = common_genes,
      R2_before = r_before,
      R2_after = r_after,
      RMSE_before = rmse_before,
      RMSE_after = rmse_after,
      Cosine_before = cosine_before,
      Cosine_after = cosine_after,
      stringsAsFactors = FALSE
    )
    
    list_name <- paste0(sampleA, "_vs_", sampleB)
    pairwise_qc_list[[list_name]] <- pairwise_df
  }
  
  return(pairwise_qc_list)
}


# Generate calibration QC plots for a list returned by
# compare_profiles_before_after_calibration_maple().
plot_calibration_qc_maple <- function(
    pairwise_qc_list,
    plot_type = c("R2", "RMSE", "Cosine")
) {
  
  plot_type <- match.arg(plot_type)
  plot_list <- list()
  
  for (pair in names(pairwise_qc_list)) {
    
    df_pair <- pairwise_qc_list[[pair]]
    
    if (plot_type == "R2") {
      
      p <- ggplot(df_pair, aes(x = R2_before, y = R2_after)) +
        geom_point(alpha = 0.4) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
        coord_fixed() +
        xlim(0, 1) +
        ylim(0, 1) +
        labs(
          title = pair,
          x = "R2 before calibration",
          y = "R2 after calibration"
        ) +
        theme_minimal() +
        theme(
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
        )
      
    } else if (plot_type == "RMSE") {
      
      df_rmse <- df_pair %>%
        tidyr::pivot_longer(
          cols = c(RMSE_before, RMSE_after),
          names_to = "Stage",
          values_to = "RMSE"
        ) %>%
        dplyr::mutate(
          Stage = dplyr::recode(
            Stage,
            RMSE_before = "Before calibration",
            RMSE_after = "After calibration"
          )
        )
      
      p <- ggplot(df_rmse, aes(x = RMSE, fill = Stage)) +
        geom_density(alpha = 0.45) +
        labs(
          title = pair,
          x = "Profile RMSE",
          y = "Density"
        ) +
        theme_minimal() +
        theme(
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          legend.position = "bottom"
        )
      
    } else if (plot_type == "Cosine") {
      
      df_cosine <- df_pair %>%
        tidyr::pivot_longer(
          cols = c(Cosine_before, Cosine_after),
          names_to = "Stage",
          values_to = "Cosine"
        ) %>%
        dplyr::mutate(
          Stage = dplyr::recode(
            Stage,
            Cosine_before = "Before calibration",
            Cosine_after = "After calibration"
          )
        )
      
      p <- ggplot(df_cosine, aes(x = Cosine, fill = Stage)) +
        geom_density(alpha = 0.45) +
        coord_cartesian(xlim = c(0, 1)) +
        labs(
          title = pair,
          x = "Cosine similarity",
          y = "Density"
        ) +
        theme_minimal() +
        theme(
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          legend.position = "bottom"
        )
    }
    
    plot_list[[pair]] <- p
  }
  
  return(plot_list)
}


# Save a list of ggplot objects into a multipage PDF.
save_multipage_pdf_maple <- function(
    plot_list,
    output_file,
    width = 10,
    height = 10,
    plots_per_page = 9,
    ncol = 3
) {
  
  if (length(plot_list) == 0) {
    stop("plot_list is empty.")
  }
  
  n_pages <- ceiling(length(plot_list) / plots_per_page)
  
  pdf(output_file, width = width, height = height)
  
  for (page in seq_len(n_pages)) {
    
    idx_start <- (page - 1) * plots_per_page + 1
    idx_end <- min(page * plots_per_page, length(plot_list))
    
    page_plots <- plot_list[idx_start:idx_end]
    
    gridExtra::grid.arrange(
      grobs = page_plots,
      ncol = ncol
    )
  }
  
  dev.off()
  
  message("Multipage PDF saved to: ", output_file)
}
