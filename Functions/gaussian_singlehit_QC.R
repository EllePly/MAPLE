# Gaussian fitting QC functions for MAPLE
#
# These functions are used to evaluate Gaussian fitting results and support
# parameter optimisation. In particular, FWHM distributions can be inspected
# across biological replicates to assess whether selected Gaussian fitting
# parameters produce realistic and consistent peak widths.

library(dplyr)
library(ggplot2)
library(purrr)
library(gridExtra)


calculate_FWHM_maple <- function(single_result) {
  
  for (name in names(single_result)) {
    
    single_result[[name]] <- lapply(single_result[[name]], function(entry) {
      
      if (
        is.list(entry) &&
        "coefs" %in% names(entry) &&
        is.list(entry$coefs) &&
        "sigma" %in% names(entry$coefs)
      ) {
        
        sigma <- entry$coefs$sigma
        FWHM <- 2 * sqrt(2 * log(2)) * sigma
        
        entry$coefs[["FWHM"]] <- unname(FWHM)
        return(entry)
        
      } else {
        
        message("Skipping entry with unexpected structure in sample: ", name)
        return(entry)
      }
    })
  }
  
  return(single_result)
}


extract_FWHM_table_maple <- function(
    gaussians_with_fwhm,
    condition_fun = function(s) sub("\\d+$", "", s)
) {
  
  fwhm_df <- purrr::map_dfr(names(gaussians_with_fwhm), function(sample_name) {
    
    fwhm_values <- sapply(
      gaussians_with_fwhm[[sample_name]],
      function(entry) {
        if (
          is.list(entry) &&
          "coefs" %in% names(entry) &&
          "FWHM" %in% names(entry$coefs)
        ) {
          entry$coefs$FWHM
        } else {
          NA_real_
        }
      }
    )
    
    data.frame(
      Sample = sample_name,
      Condition = condition_fun(sample_name),
      FWHM = as.numeric(unname(fwhm_values))
    )
  }) %>%
    dplyr::filter(!is.na(FWHM), is.finite(FWHM))
  
  return(fwhm_df)
}


plot_FWHM_distribution_maple <- function(
    gaussians_with_fwhm,
    condition_fun = function(s) sub("\\d+$", "", s),
    density = TRUE,
    binwidth = 0.5,
    xlim = c(3, 10),
    alpha = 0.45,
    ncol = 3,
    nrow = NULL,
    save_plot = FALSE,
    output_file = "FWHM_distribution_by_condition.png",
    width = 8,
    height = 10
) {
  
  fwhm_df <- extract_FWHM_table_maple(
    gaussians_with_fwhm = gaussians_with_fwhm,
    condition_fun = condition_fun
  )
  
  if (nrow(fwhm_df) == 0) {
    stop("No valid FWHM values were found.")
  }
  
  plot_list <- list()
  
  for (cond in unique(fwhm_df$Condition)) {
    
    df_cond <- fwhm_df %>%
      dplyr::filter(Condition == cond)
    
    p <- ggplot(df_cond, aes(x = FWHM, fill = Sample))
    
    if (isTRUE(density)) {
      p <- p +
        geom_histogram(
          aes(y = after_stat(density)),
          binwidth = binwidth,
          alpha = alpha,
          position = "identity",
          color = "black",
          linewidth = 0.2
        ) +
        labs(y = "Density")
    } else {
      p <- p +
        geom_histogram(
          binwidth = binwidth,
          alpha = alpha,
          position = "identity",
          color = "black",
          linewidth = 0.2
        ) +
        labs(y = "Count")
    }
    
    p <- p +
      coord_cartesian(xlim = xlim) +
      labs(
        title = paste("FWHM distribution -", cond),
        x = "FWHM",
        fill = "Replicate"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 9, color = "black"),
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      )
    
    plot_list[[cond]] <- p
  }
  
  combined_plot <- gridExtra::arrangeGrob(
    grobs = plot_list,
    ncol = ncol,
    nrow = nrow
  )
  
  if (isTRUE(save_plot)) {
    ggsave(
      output_file,
      plot = combined_plot,
      width = width,
      height = height,
      bg = "white"
    )
    
    message("FWHM distribution plot saved to: ", output_file)
  }
  
  return(combined_plot)
}