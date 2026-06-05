# Pairwise protein-protein co-elution and complex prediction workflow
#
# This workflow demonstrates how to calculate pairwise protein co-elution
# features within one sample and use a reference protein complex annotation
# matrix, such as CORUM, to train/predict protein interactions using PrInCE.


# -------------------------------------------------------------------------
# 1. Select one sample for pairwise co-elution analysis
# -------------------------------------------------------------------------

# Select one sample from the calibrated profile list.
sample_use <- names(profile_list_calibrated)[1]

profile_test <- profile_list_calibrated[[sample_use]]
gaussian_test <- gaussian_results[[sample_use]]


# -------------------------------------------------------------------------
# 2. Subset profiles for a minimal test run
# -------------------------------------------------------------------------

# For a full analysis, use the complete profile matrix.
# For testing, only the first 100 proteins are used to reduce runtime.
profile_test_small <- profile_test[
  1:min(100, nrow(profile_test)),
  ,
  drop = FALSE
]

gaussian_test_small <- gaussian_test[
  rownames(profile_test_small)
]


# -------------------------------------------------------------------------
# 3. Calculate pairwise protein-protein co-elution features
# -------------------------------------------------------------------------

# This step compares Protein A vs Protein B within the same sample.
# The resulting feature table can be used for complex or interaction prediction.
#
# Returned features include:
#   cor_R              : Pearson-correlation-derived distance
#   cor_P              : Pearson correlation p value
#   cosine_distance    : cosine-distance-like feature
#   euclidean_distance : Euclidean distance between profiles
#   co_peak            : distance between maximum-intensity fraction positions
#   co_apex            : Gaussian component co-apex distance

pairwise_maple <- calculate_pairwise_features_maple(
  profile_matrix = profile_test_small,
  gaussians = gaussian_test_small,
  pearson_R = TRUE,
  pearson_P = TRUE,
  cosine_correlation = TRUE,
  euclidean_distance = TRUE,
  co_peak = TRUE,
  co_apex = TRUE,
  show_progress = TRUE
)


# -------------------------------------------------------------------------
# 4. Load reference network construction functions
# -------------------------------------------------------------------------

source("Functions/reference_networks.R")


# -------------------------------------------------------------------------
# 5. Load protein complex annotation table
# -------------------------------------------------------------------------

# The input complex annotation should be in wide format:
#   - Each column represents one protein complex.
#   - Each non-empty row value is one protein member.
#
# Example:
#   Complex_I     Complex_II
#   NDUFA1        SDHA
#   NDUFA2        SDHB
#   NDUFS1        SDHC

CORUM <- read.csv("CORUM_complexes_wide_format.csv")

corum_list <- wide_complex_table_to_list_maple(CORUM)


# -------------------------------------------------------------------------
# 6. Optional identifier mapping
# -------------------------------------------------------------------------

# This step is only needed if the protein identifiers in the complex annotation
# do not match the identifiers used in the MAPLE profile matrix.
#
# Example:
#   CORUM table uses UniProt IDs
#   MAPLE profiles use Gene Names
#
# In that case, provide a two-column mapping table:
#   UniprotID    GeneName

uniprot_to_gene <- read.delim("path/to/uniprot_to_gene_mapping.tsv")

colnames(uniprot_to_gene) <- c("UniprotID", "GeneName")

CORUM_genename <- map_complex_ids_maple(
  complex_list = corum_list,
  id_map = uniprot_to_gene,
  from_col = "UniprotID",
  to_col = "GeneName"
)


# If the complex annotation already uses the same identifiers as the MAPLE
# profile matrix, skip the mapping step and use:
#
# CORUM_genename <- corum_list


# -------------------------------------------------------------------------
# 7. Build reference adjacency matrix
# -------------------------------------------------------------------------

# The reference adjacency matrix labels protein pairs from the same annotated
# complex as positive pairs.
#
#   same complex     = 1
#   different/no annotation = 0

ref_CORUM <- build_reference_adjacency_from_complex_list_maple(
  complex_list = CORUM_genename
)


# -------------------------------------------------------------------------
# 8. Predict protein interactions using PrInCE
# -------------------------------------------------------------------------

# This step uses the PrInCE predict_interactions() function.
# MAPLE provides the pairwise feature table, while PrInCE performs the
# supervised interaction prediction using the reference adjacency matrix.

pairwise_prediction <- predict_interactions(
  pairwise_maple,
  ref_CORUM,
  classifier = "LR",
  verbose = FALSE,
  models = 1,
  cv_folds = 5,
  trees = 500
)