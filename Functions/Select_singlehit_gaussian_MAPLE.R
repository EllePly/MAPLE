# Single-hit Gaussian peak selection for MAPLE
#
# This function selects proteins with one dominant Gaussian peak within a
# user-defined peak height and peak-centre range. It also identifies proteins
# that pass the single-hit criteria consistently across replicates within each
# biological condition.

select_single_hit_gaussians_maple <- function(
    gaussian,
    data_list,
    samples = NULL,
    A_threshold = 1,
    Mu_threshold_max = 37,
    Mu_threshold_min = 8,
    condition_fun = function(s) sub("\\d+$", "", s),
    export_global = FALSE
) {
  
  single_result <- list()
  
  # 1. Automatically identify samples shared between Gaussian results and data list.
  if (is.null(samples)) {
    samples <- intersect(names(gaussian), names(data_list))
  }
  
  if (length(samples) == 0) {
    stop("No overlapping samples between gaussian and data_list.")
  }
  
  # 2. Select single-hit Gaussian proteins for each sample.
  for (sample in samples) {
    
    if (!sample %in% names(data_list)) {
      message("Sample ", sample, " is not present in data_list. Skipping.")
      next
    }
    
    if (!sample %in% names(gaussian)) {
      message("Sample ", sample, " is not present in gaussian. Skipping.")
      next
    }
    
    gene_list <- rownames(data_list[[sample]])
    
    filtered_genes <- lapply(gene_list, function(gene) {
      
      entry <- gaussian[[sample]][[gene]]
      
      if (is.null(entry)) return(NULL)
      if (
        is.null(entry$coefs) ||
        is.null(entry$coefs$A) ||
        is.null(entry$coefs$mu) ||
        is.null(entry$coefs$sigma)
      ) {
        return(NULL)
      }
      
      A_values <- entry$coefs$A
      Mu_values <- entry$coefs$mu
      Sigma_values <- entry$coefs$sigma
      
      # Skip entries with inconsistent coefficient lengths.
      if (
        length(A_values) == 0 ||
        length(Mu_values) != length(A_values) ||
        length(Sigma_values) != length(A_values)
      ) {
        return(NULL)
      }
      
      # Retain proteins with exactly one Gaussian peak above the height threshold.
      if (sum(A_values > A_threshold) == 1) {
        
        max_A_index <- which.max(A_values)
        
        if (
          Mu_values[max_A_index] >= Mu_threshold_min &&
          Mu_values[max_A_index] <= Mu_threshold_max
        ) {
          
          filtered_entry <- list(
            n_gaussians = 1,
            R2 = entry$R2,
            iterations = entry$iterations,
            coefs = list(
              A = A_values[max_A_index],
              mu = Mu_values[max_A_index],
              sigma = Sigma_values[max_A_index]
            ),
            curveFit = entry$curveFit
          )
          
          return(setNames(list(filtered_entry), gene))
        }
      }
      
      return(NULL)
    })
    
    single_result[[sample]] <- do.call(
      c,
      Filter(Negate(is.null), filtered_genes)
    )
    
    message(
      "Sample ", sample, " finished. Retained ",
      length(single_result[[sample]]), " genes."
    )
  }
  
  # 3. Infer biological condition for each sample.
  sample_conditions <- vapply(
    names(single_result),
    condition_fun,
    character(1)
  )
  
  condition_samples <- split(names(single_result), sample_conditions)
  
  # 4. Identify common single-hit proteins across replicates within each condition.
  singlehit_gene_lists <- lapply(single_result, names)
  
  condition_common_genes <- lapply(names(condition_samples), function(cond) {
    
    reps <- condition_samples[[cond]]
    
    Reduce(intersect, singlehit_gene_lists[reps])
  })
  
  names(condition_common_genes) <- names(condition_samples)
  
  # 5. Return Gaussian results restricted to condition-level common genes.
  condition_common_gaussians <- lapply(names(condition_samples), function(cond) {
    
    reps <- condition_samples[[cond]]
    
    out <- lapply(reps, function(s) {
      single_result[[s]][condition_common_genes[[cond]]]
    })
    
    names(out) <- reps
    out
  })
  
  names(condition_common_gaussians) <- names(condition_samples)
  
  # 6. Print summary.
  message("Detected ", length(condition_samples), " biological conditions:")
  
  for (cond in names(condition_samples)) {
    message(
      " - ", cond, ": ",
      length(condition_samples[[cond]]), " replicates, common genes = ",
      length(condition_common_genes[[cond]])
    )
  }
  
  # 7. Optional export to the global environment.
  # This is disabled by default and is mainly retained for interactive analysis.
  if (isTRUE(export_global)) {
    assign("single_result", single_result, envir = .GlobalEnv)
    assign("condition_samples", condition_samples, envir = .GlobalEnv)
    assign("condition_common_genes", condition_common_genes, envir = .GlobalEnv)
    assign("condition_common_gaussians", condition_common_gaussians, envir = .GlobalEnv)
    
    message(
      "Objects exported to the global environment: ",
      "single_result, condition_samples, condition_common_genes, ",
      "condition_common_gaussians."
    )
  }
  
  invisible(list(
    single_result = single_result,
    condition_samples = condition_samples,
    condition_common_genes = condition_common_genes,
    condition_common_gaussians = condition_common_gaussians,
    condition_of_sample = sample_conditions
  ))
}