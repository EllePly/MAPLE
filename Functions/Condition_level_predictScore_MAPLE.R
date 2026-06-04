
# Condition-level profile remodelling analysis for MAPLE
#
# This module builds a condition-level protein profile remodelling score.
#
# Core workflow:
#   1. Average replicate protein profiles by condition.
#   2. Generate leave-one-out pseudo-average comparisons within each condition.
#   3. Use pseudo-average within-condition comparisons as the expected
#      background profile variability.
#   4. Train an Isolation Forest model on within-condition feature distributions.
#   5. Compare condition-averaged profiles against a reference condition.
#   6. Assign each protein a condition-level profile remodelling score.
#
# Required packages:
#   dplyr
#   purrr
#   tibble
#   stringr
#   isotree
#
# Required MAPLE functions:
#   calculate_features_intersample_maple()
#   convert_similarity_to_distance_maple()

library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(isotree)


# Average replicate protein profiles by biological condition.
average_by_condition_maple <- function(
    data_list,
    groups,
    condition_match = c("prefix", "exact")
) {
  
  condition_match <- match.arg(condition_match)
  condition_average_list <- list()
  
  for (group in groups) {
    
    if (condition_match == "prefix") {
      samples <- names(data_list)[grepl(paste0("^", group), names(data_list))]
    } else {
      samples <- names(data_list)[names(data_list) == group]
    }
    
    if (length(samples) == 0) {
      warning("No samples found for group: ", group)
      next
    }
    
    all_genes <- Reduce(union, lapply(samples, function(s) {
      rownames(data_list[[s]])
    }))
    
    common_cols <- Reduce(intersect, lapply(samples, function(s) {
      colnames(data_list[[s]])
    }))
    
    common_cols <- common_cols[
      order(as.numeric(gsub("^F", "", common_cols)))
    ]
    
    arr <- array(
      NA_real_,
      dim = c(length(all_genes), length(common_cols), length(samples)),
      dimnames = list(all_genes, common_cols, samples)
    )
    
    for (i in seq_along(samples)) {
      
      s <- samples[i]
      
      genes_use <- intersect(all_genes, rownames(data_list[[s]]))
      
      mat <- data_list[[s]][genes_use, common_cols, drop = FALSE] %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
        as.matrix()
      
      arr[genes_use, , i] <- mat
    }
    
    avg_mat <- apply(arr, c(1, 2), function(x) {
      if (all(is.na(x))) {
        NA_real_
      } else {
        mean(x, na.rm = TRUE)
      }
    })
    
    avg_mat <- as.data.frame(avg_mat)
    
    condition_average_list[[group]] <- avg_mat
    
    message(
      "Averaged condition: ", group,
      " | samples: ", paste(samples, collapse = ", "),
      " | genes: ", nrow(avg_mat),
      " | fractions: ", ncol(avg_mat)
    )
  }
  
  return(condition_average_list)
}


# Generate original and leave-one-out pseudo-average profiles for one condition.
make_leave_one_out_profile_list_maple <- function(
    data_list,
    group,
    min_replicates = 2
) {
  
  samples <- names(data_list)[grepl(paste0("^", group, "[0-9]+$"), names(data_list))]
  
  if (length(samples) < min_replicates) {
    stop(
      "Expected at least ", min_replicates, " samples for group ",
      group, ", but found: ", paste(samples, collapse = ", ")
    )
  }
  
  out <- list()
  
  for (s in samples) {
    out[[s]] <- as.matrix(data_list[[s]])
  }
  
  for (held_out in samples) {
    
    avg_samples <- setdiff(samples, held_out)
    
    common_genes <- Reduce(intersect, lapply(avg_samples, function(s) {
      rownames(data_list[[s]])
    }))
    
    common_cols <- Reduce(intersect, lapply(avg_samples, function(s) {
      colnames(data_list[[s]])
    }))
    
    common_cols <- common_cols[
      order(as.numeric(gsub("^F", "", common_cols)))
    ]
    
    avg_mat <- Reduce("+", lapply(avg_samples, function(s) {
      as.matrix(data_list[[s]][common_genes, common_cols, drop = FALSE])
    })) / length(avg_samples)
    
    avg_id <- paste(gsub(group, "", avg_samples), collapse = "")
    avg_name <- paste0(group, "_AVG", avg_id)
    
    out[[avg_name]] <- avg_mat
  }
  
  return(out)
}


# Generate leave-one-out pseudo-average comparison table for one condition.
make_leave_one_out_comparison_table_maple <- function(
    data_list,
    group,
    min_replicates = 2
) {
  
  samples <- names(data_list)[grepl(paste0("^", group, "[0-9]+$"), names(data_list))]
  
  if (length(samples) < min_replicates) {
    stop(
      "Expected at least ", min_replicates, " samples for group ",
      group, ", but found: ", paste(samples, collapse = ", ")
    )
  }
  
  comp_df <- purrr::map_dfr(samples, function(held_out) {
    
    avg_samples <- setdiff(samples, held_out)
    
    avg_id <- paste(gsub(group, "", avg_samples), collapse = "")
    avg_name <- paste0(group, "_AVG", avg_id)
    
    tibble::tibble(
      group = group,
      sampleA = avg_name,
      sampleB = held_out,
      comparison = paste0(avg_name, "_vs_", held_out),
      comparison_type = "within_condition_pseudo_average",
      avg_samples = paste(avg_samples, collapse = ";"),
      held_out = held_out
    )
  })
  
  return(comp_df)
}


# Build pseudo-average within-condition training set for Isolation Forest.
build_pseudo_average_training_set_maple <- function(
    data_list,
    groups,
    model_features = c(
      "pearson_R_raw",
      "pearson_R_cleaned",
      "cosine_sim",
      "euclidean",
      "co_peak_diff"
    ),
    min_replicates = 2
) {
  
  pseudo_profile_list <- list()
  pseudo_comparison_table <- list()
  
  for (g in groups) {
    
    pseudo_profiles_g <- make_leave_one_out_profile_list_maple(
      data_list = data_list,
      group = g,
      min_replicates = min_replicates
    )
    
    pseudo_comparisons_g <- make_leave_one_out_comparison_table_maple(
      data_list = data_list,
      group = g,
      min_replicates = min_replicates
    )
    
    pseudo_profile_list <- c(pseudo_profile_list, pseudo_profiles_g)
    pseudo_comparison_table[[g]] <- pseudo_comparisons_g
  }
  
  pseudo_comparison_table <- dplyr::bind_rows(pseudo_comparison_table)
  
  results_pseudo_average_within <- list()
  
  for (i in seq_len(nrow(pseudo_comparison_table))) {
    
    sampleA <- pseudo_comparison_table$sampleA[i]
    sampleB <- pseudo_comparison_table$sampleB[i]
    comp_name <- pseudo_comparison_table$comparison[i]
    
    message("Running pseudo-average comparison: ", comp_name)
    
    res <- calculate_features_intersample_maple(
      pseudo_profile_list[[sampleA]],
      pseudo_profile_list[[sampleB]],
      pearson_R_raw = TRUE,
      pearson_R_cleaned = TRUE,
      pearson_P = TRUE,
      cosine_correlation = TRUE,
      euclidean_distance = TRUE,
      co_peak = TRUE,
      co_apex = FALSE
    )
    
    results_pseudo_average_within[[comp_name]] <- res
  }
  
  results_pseudo_average_within_directional <-
    convert_similarity_to_distance_maple(results_pseudo_average_within)
  
  pseudo_train_long <- purrr::imap_dfr(
    results_pseudo_average_within_directional,
    function(df, nm) {
      
      if (!("protein" %in% colnames(df)) && !("GeneName" %in% colnames(df))) {
        df <- tibble::rownames_to_column(df, "protein")
      }
      
      if ("GeneName" %in% colnames(df)) {
        df <- df %>%
          dplyr::rename(protein = GeneName)
      }
      
      df %>%
        dplyr::mutate(
          comparison = nm,
          comparison_type = "within_condition_pseudo_average"
        )
    }
  )
  
  X_train <- pseudo_train_long %>%
    dplyr::select(dplyr::all_of(model_features)) %>%
    as.data.frame()
  
  return(list(
    pseudo_profile_list = pseudo_profile_list,
    pseudo_comparison_table = pseudo_comparison_table,
    pseudo_feature_results = results_pseudo_average_within_directional,
    pseudo_train_long = pseudo_train_long,
    X_train = X_train,
    model_features = model_features
  ))
}


# Train Isolation Forest model using pseudo-average within-condition background.
train_condition_remodelling_model_maple <- function(
    X_train,
    ntrees = 1000,
    sample_size = min(1024, nrow(X_train)),
    ndim = 1,
    nthreads = 4,
    seed = 100
) {
  
  if (anyNA(X_train)) {
    stop("X_train contains NA values. Please remove or impute missing values before training.")
  }
  
  set.seed(seed)
  
  model <- isotree::isolation.forest(
    X_train,
    ntrees = ntrees,
    sample_size = sample_size,
    ndim = ndim,
    nthreads = nthreads
  )
  
  return(model)
}


# Build condition-average comparison table against one reference condition.
make_condition_average_comparison_table_maple <- function(
    query_conditions,
    reference_condition
) {
  
  tibble::tibble(
    query_condition = query_conditions,
    reference_condition = reference_condition,
    sampleA = query_condition,
    sampleB = reference_condition,
    comparison = paste0(query_condition, "_vs_", reference_condition)
  )
}


# Score condition-level profile remodelling using a trained Isolation Forest model.
score_condition_remodelling_maple <- function(
    condition_average_profiles,
    model,
    reference_condition,
    query_conditions = NULL,
    model_features = c(
      "pearson_R_raw",
      "pearson_R_cleaned",
      "cosine_sim",
      "euclidean",
      "co_peak_diff"
    ),
    score_column = "IF_score_conditionLevel"
) {
  
  if (is.null(query_conditions)) {
    query_conditions <- setdiff(names(condition_average_profiles), reference_condition)
  }
  
  avg_comparison_table <- make_condition_average_comparison_table_maple(
    query_conditions = query_conditions,
    reference_condition = reference_condition
  )
  
  results_condition_average <- list()
  
  for (i in seq_len(nrow(avg_comparison_table))) {
    
    sampleA <- avg_comparison_table$sampleA[i]
    sampleB <- avg_comparison_table$sampleB[i]
    comp_name <- avg_comparison_table$comparison[i]
    
    message("Running condition-average comparison: ", comp_name)
    
    res <- calculate_features_intersample_maple(
      condition_average_profiles[[sampleA]],
      condition_average_profiles[[sampleB]],
      pearson_R_raw = TRUE,
      pearson_R_cleaned = TRUE,
      pearson_P = TRUE,
      cosine_correlation = TRUE,
      euclidean_distance = TRUE,
      co_peak = TRUE,
      co_apex = FALSE
    )
    
    results_condition_average[[comp_name]] <- res
  }
  
  results_condition_average_directional <-
    convert_similarity_to_distance_maple(results_condition_average)
  
  condition_avg_long <- purrr::imap_dfr(
    results_condition_average_directional,
    function(df, nm) {
      
      if (!("protein" %in% colnames(df)) && !("GeneName" %in% colnames(df))) {
        df <- tibble::rownames_to_column(df, "protein")
      }
      
      if ("GeneName" %in% colnames(df)) {
        df <- df %>%
          dplyr::rename(protein = GeneName)
      }
      
      parts <- stringr::str_split(nm, "_vs_", simplify = TRUE)
      
      df %>%
        dplyr::mutate(
          comparison = nm,
          query_condition = parts[1],
          reference_condition = parts[2],
          comparison_type = "between_condition_average"
        )
    }
  )
  
  X_condition <- condition_avg_long %>%
    dplyr::select(dplyr::all_of(model_features)) %>%
    as.data.frame()
  
  if (anyNA(X_condition)) {
    stop("X_condition contains NA values. Please remove or impute missing values before scoring.")
  }
  
  condition_avg_long[[score_column]] <- predict(
    model,
    X_condition,
    type = "score"
  )
  
  return(list(
    condition_comparison_table = avg_comparison_table,
    condition_feature_results = results_condition_average_directional,
    condition_score_long = condition_avg_long,
    X_condition = X_condition
  ))
}


# Summarise condition-level remodelling scores.
summarise_condition_remodelling_scores_maple <- function(
    condition_score_long,
    score_column = "IF_score_conditionLevel"
) {
  
  condition_score_long %>%
    dplyr::group_by(query_condition, reference_condition) %>%
    dplyr::summarise(
      n_proteins = dplyr::n_distinct(protein),
      median_IF = median(.data[[score_column]], na.rm = TRUE),
      mean_IF = mean(.data[[score_column]], na.rm = TRUE),
      q75_IF = quantile(.data[[score_column]], 0.75, na.rm = TRUE),
      max_IF = max(.data[[score_column]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(median_IF))
}