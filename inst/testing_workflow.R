# 1.) Load and configure markers -----------------------------------------------

markers <- load_markers(file_path = "~/Desktop/example_data/input_markers.xlsx")

markers$broad <- build_broad_marker_config(
  marker_list    = markers$broad,
  priority_order = c("vasculature", "immune"),
  default_threshold = 0.25
)

# 2.) Create SCE ---------------------------------------------------------------

sce <- create_sce(
  counts        = "~/Desktop/example_data/test_input_sc_matrix.rds",   # or a dgCMatrix directly
  cell_metadata = "~/Desktop/example_data/cell_metadata.rds"  # optional, also accepts .csv/.tsv
)

input_dir <- "~/Desktop/OneDrive - University of Leeds/analysis/Multiomics/Data/cellvoter_inputs/"

patient_samples <- list.files(input_dir, full.names = TRUE)

sce <- create_sce(counts  =  patient_samples[4])

# 3.) QC and removal of  low-quality cells (optional) --------------------------

sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)

# 4.) Normalize the data -------------------------------------------------------

sce <- normalize_counts(sce)

# 5.) Build analysis tracks ----------------------------------------------------

# This clusters both gene sets, attaches marker config and filtered fine markers
# to metadata, and stores the reduced altExp under "user_panel"

# Need to look into this where we have clusters with a small numbe of cells we skip
# assignment, but in cases where we have a cluster that has just one cell, we have an
# issue with a variance calculation
sce <- prepare_sce(sce, markers, k = 20, resolution = 1.2, n_hvgs = 2000, n_pcs = 50)

# 6.) Run ensemble annotation --------------------------------------------------
# Default run
results <- run_cellvoter(sce)

# default run with full output (scores, p-values, etc. for all methods and clusters)
results <- run_cellvoter(sce, return_full_output = TRUE)

# Or with custom parameters
results_diff_params <- run_cellvoter(
  sce,
  return_full_output = TRUE,
  annotation_args = list(
    broad_args   = list(test_type = "t", min_prop = 0.1),
    rank_args    = list(test_type = "t", min_prop = 0.1),
    extract_args = list(fdr_threshold = 0.01)
  )
)

# 7.) Assign ensemble labels ---------------------------------------------------
# you can tweak and re-run this step freely without re-running step 6
consensus <- resolve_consensus_labels(
  label_list        = results$labels,
  method_names      = results$method_names,
  tie_breaker_names = results$tie_breaker_names,
  unassigned_label  = "unknown",
  allow_even_split  = FALSE,
  ordered_tiebreak  = TRUE
)

# 8.) Attach final labels to SCE -----------------------------------------------

results$sce$cellVoteR_label  <- consensus$label
results$sce$cellVoteR_method <- consensus$method

# 9.) Inspect ------------------------------------------------------------------
table(results$sce$cellVoteR_label)
table(results$sce$cellVoteR_method)

# Per-method breakdown before consensus
table(results$labels$method_1)
table(results$labels$method_2)
table(results$labels$method_3)
table(results$labels$method_4)
table(results$labels$global_1)
table(results$labels$global_2)

# Full scores per cluster per method (if return_full_output = TRUE)
results$full_output$method_1$scores
results$full_output$global_1$scores

# Viewing the Vignette (Probably needs some updating for v2) -------------------

vignette("CellVoteR", package = "CellVoteR")

# END --------------------------------------------------------------------------
