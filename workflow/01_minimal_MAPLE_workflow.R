# Minimal MAPLE workflow
#
# This script demonstrates the core MAPLE analysis workflow using a small
# example protein profile list.
#
# Workflow:
#   1. Load MAPLE functions
#   2. Load example protein elution profiles
#   3. Run Gaussian fitting
#   4. Select single-hit Gaussian proteins
#   5. Perform Gaussian FWHM QC
#   6. Estimate global fraction-axis calibration
#   7. Apply calibration
#   8. Evaluate calibration QC
#   9. Calculate inter-sample profile features
#   10. Calculate condition-level profile remodelling scores


# -------------------------------------------------------------------------
# 1. Load required MAPLE functions
# -------------------------------------------------------------------------

source("Functions/gaussian_fitting.R")
source("Functions/single_hit_selection.R")
source("Functions/gaussian_qc.R")
source("Functions/calibration_functions.R")
source("Functions/calibration_qc.R")
source("Functions/intersample_features.R")
source("Functions/condition_level_remodelling.R")


source("C://Users/62682/OneDrive - The Kids Research Institute Australia/Manuscripts/Co-fractionation lipid-protein interaction manuscript/Github/Select_singlehit_gaussian_MAPLE.R")
source("C://Users/62682/OneDrive - The Kids Research Institute Australia/Manuscripts/Co-fractionation lipid-protein interaction manuscript/Github/intersample_features_calculation_MAPLE.R")
source("C://Users/62682/OneDrive - The Kids Research Institute Australia/Manuscripts/Co-fractionation lipid-protein interaction manuscript/Github/gaussian_singlehit_QC.R")
source("C://Users/62682/OneDrive - The Kids Research Institute Australia/Manuscripts/Co-fractionation lipid-protein interaction manuscript/Github/gaussian_fitting_MAPLE.R")
source("C://Users/62682/OneDrive - The Kids Research Institute Australia/Manuscripts/Co-fractionation lipid-protein interaction manuscript/Github/Condition_level_predictScore_MAPLE.R")
source("C://Users/62682/OneDrive - The Kids Research Institute Australia/Manuscripts/Co-fractionation lipid-protein interaction manuscript/Github/calibration_QC.R")
source("C://Users/62682/OneDrive - The Kids Research Institute Australia/Manuscripts/Co-fractionation lipid-protein interaction manuscript/Github/calibration_functions_MAPLE.R")


# -------------------------------------------------------------------------
# 2. Load example data
# -------------------------------------------------------------------------

# The example dataset should be a named list of protein elution profile
# matrices or data frames. Each list element is one sample.
#
# Expected sample names in this minimal example:
#   WTold1, WTold2, WTold3
#   CRLS1KOold1, CRLS1KOold2, CRLS1KOold3

data_list <- readRDS(
  "Example_Data/test_profile_list_QC_WTold_CRLS1KOold.rds"
)

groups_use <- c("WTold", "CRLS1KOold")
reference_sample <- "WTold2"
reference_condition <- "WTold"


# -------------------------------------------------------------------------
# 3. Fill missing fractions
# -------------------------------------------------------------------------

profile_list_filled <- fill_missing_fractions_maple(
  data_list = data_list,
  target_fractions = paste0("F", 38:75),
  fill_boundary_na_with_zero = TRUE
)


# -------------------------------------------------------------------------
# 4. Gaussian fitting
# -------------------------------------------------------------------------

gaussian_results <- lapply(profile_list_filled, function(profile_matrix) {
  
  build_gaussians_maple(
    profile_matrix = profile_matrix,
    min_points = 1,
    min_consecutive = 5,
    impute_NA = FALSE,
    smooth = FALSE,
    smooth_width = 1,
    max_gaussians = 3,
    criterion = "AIC",
    max_iterations = 50,
    min_R_squared = 0.9,
    use_R_squared = TRUE,
    use_RMSE = TRUE,
    max_RMSE = 0.1,
    pass_rule = "both",
    optimize_by = "RMSE",
    method = "random",
    filter_gaussians_center = TRUE,
    filter_gaussians_height = 0.1,
    filter_gaussians_variance_min = 0.1,
    filter_gaussians_variance_max = 10,
    show_progress = TRUE
  )
})


# -------------------------------------------------------------------------
# 5. Single-hit Gaussian protein selection
# -------------------------------------------------------------------------

single_hit_result <- select_single_hit_gaussians_maple(
  gaussian = gaussian_results,
  data_list = profile_list_filled,
  A_threshold = 0.9,
  Mu_threshold_min = 38,
  Mu_threshold_max = 75,
  export_global = FALSE
)

single_hit_common <- unlist(
  single_hit_result$condition_common_gaussians,
  recursive = FALSE
)

names(single_hit_common) <- sub("^.*\\.", "", names(single_hit_common))


# -------------------------------------------------------------------------
# 6. Gaussian FWHM QC
# -------------------------------------------------------------------------

single_hit_common_fwhm <- calculate_FWHM_maple(
  single_hit_common
)

p_fwhm <- plot_FWHM_distribution_maple(
  gaussians_with_fwhm = single_hit_common_fwhm,
  density = TRUE,
  binwidth = 0.25,
  xlim = c(0, 15),
  ncol = 1,
  save_plot = FALSE
)

# To display in RStudio:
 grid::grid.newpage()
grid::grid.draw(p_fwhm)


# -------------------------------------------------------------------------
# 7. Estimate global calibration parameters
# -------------------------------------------------------------------------

calibration_result <- estimate_all_to_reference_calibration_maple(
  gaussians_singlehit_list = single_hit_common,
  reference_sample = reference_sample,
  weighted_lm = TRUE,
  return_plots = TRUE
)

calibration_parameters <- calibration_result$calparam


# -------------------------------------------------------------------------
# 8. Apply calibration
# -------------------------------------------------------------------------

profile_list_calibrated <- calibrate_profile_list_maple(
  data_list = profile_list_filled,
  calparam = calibration_parameters,
  reference_sample = reference_sample,
  group = "ALL"
)

profile_list_calibrated <- add_reference_sample_maple(
  calibrated_data_list = profile_list_calibrated,
  data_list = profile_list_filled,
  reference_sample = reference_sample
)


# -------------------------------------------------------------------------
# 9. Calibration QC
# -------------------------------------------------------------------------

sample_pairs_qc <- list(
  c("WTold2", "WTold1"),
  c("WTold2", "WTold3"),
  c("WTold2", "CRLS1KOold1")
)

calibration_qc <- compare_profiles_before_after_calibration_maple(
  sample_pairs = sample_pairs_qc,
  rawdata_list = profile_list_filled,
  caldata_list = profile_list_calibrated
)

r2_qc_plots <- plot_calibration_qc_maple(
  pairwise_qc_list = calibration_qc,
  plot_type = "R2"
)

# Example display:
# r2_qc_plots[[1]]


# -------------------------------------------------------------------------
# 10. Inter-sample feature calculation
# -------------------------------------------------------------------------

intersample_features <- run_all_intersample_feature_comparisons_maple(
  data_list = profile_list_calibrated,
  gaussians_list = gaussian_results,
  pearson_R_raw = TRUE,
  pearson_R_cleaned = TRUE,
  pearson_P = TRUE,
  cosine_correlation = TRUE,
  euclidean_distance = TRUE,
  co_peak = TRUE,
  co_apex = TRUE
)

intersample_features_directional <- convert_similarity_to_distance_maple(
  intersample_features
)

pair_info <- parse_pair_info_maple(
  names(intersample_features_directional)
)

z_euclidean <- compute_z_feature_all_maple(
  feature_name = "euclidean",
  dataset = intersample_features_directional,
  pair_info = pair_info
)


# -------------------------------------------------------------------------
# 11. Condition-level profile remodelling score
# -------------------------------------------------------------------------

pseudo_training_set <- build_pseudo_average_training_set_maple(
  data_list = profile_list_calibrated,
  groups = groups_use,
  min_replicates = 2
)

iso_model <- train_condition_remodelling_model_maple(
  X_train = pseudo_training_set$X_train,
  ntrees = 100,
  nthreads = 1
)

condition_average_profiles <- average_by_condition_maple(
  data_list = profile_list_calibrated,
  groups = groups_use
)

condition_scores <- score_condition_remodelling_maple(
  condition_average_profiles = condition_average_profiles,
  model = iso_model,
  reference_condition = reference_condition,
  query_conditions = c("CRLS1KOold"),
  model_features = pseudo_training_set$model_features
)

condition_score_summary <- summarise_condition_remodelling_scores_maple(
  condition_scores$condition_score_long
)

print(condition_score_summary)


# -------------------------------------------------------------------------
# End of minimal workflow
# -------------------------------------------------------------------------