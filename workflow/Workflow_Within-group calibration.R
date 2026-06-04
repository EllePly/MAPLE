# -------------------------------------------------------------------------
# Within-group calibration
# -------------------------------------------------------------------------

sample_groups <- list(
  WTold = c("WTold1", "WTold2", "WTold3"),
  CRLS1KOold = c("CRLS1KOold1", "CRLS1KOold2", "CRLS1KOold3")
)

reference_by_group <- c(
  WTold = "WTold2",
  CRLS1KOold = "CRLS1KOold1"
)

# Estimate within-group calibration parameters.
within_group_calibration_results <- list()

for (group_name in names(sample_groups)) {
  
  message("Estimating within-group calibration for: ", group_name)
  
  within_group_calibration_results[[group_name]] <-
    estimate_all_to_reference_calibration_maple(
      gaussians_singlehit_list = single_hit_result$condition_common_gaussians[[group_name]],
      reference_sample = reference_by_group[[group_name]],
      samples = sample_groups[[group_name]],
      weighted_lm = TRUE,
      return_plots = TRUE
    )
}

calparam_within_group <- dplyr::bind_rows(
  lapply(within_group_calibration_results, function(x) x$calparam)
)

# Apply within-group calibration.
profile_list_within_group_calibrated <- list()

for (group_name in names(sample_groups)) {
  
  reference_sample_group <- reference_by_group[[group_name]]
  samples_group <- sample_groups[[group_name]]
  
  message("Applying within-group calibration for: ", group_name)
  
  calibrated_group <- calibrate_profile_list_maple(
    data_list = profile_list_calibrated,
    calparam = calparam_within_group,
    reference_sample = reference_sample_group,
    samples = samples_group
  )
  
  calibrated_group <- add_reference_sample_maple(
    calibrated_data_list = calibrated_group,
    data_list = profile_list_calibrated,
    reference_sample = reference_sample_group
  )
  
  profile_list_within_group_calibrated <- c(
    profile_list_within_group_calibrated,
    calibrated_group
  )
}

# Generate within-group sample pairs for QC.
sample_pairs_within_group <- unlist(
  lapply(sample_groups, function(samples) {
    combn(samples, 2, simplify = FALSE)
  }),
  recursive = FALSE
)

within_group_calibration_qc <- compare_profiles_before_after_calibration_maple(
  sample_pairs = sample_pairs_within_group,
  rawdata_list = profile_list_calibrated,
  caldata_list = profile_list_within_group_calibrated
)

within_group_r2_plots <- plot_calibration_qc_maple(
  pairwise_qc_list = within_group_calibration_qc,
  plot_type = "R2"
)