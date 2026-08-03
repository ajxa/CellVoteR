# 1.) Load and configure markers -----------------------------------------------
markers <- load_markers(file_path = "~/Desktop/example_data/input_markers.xlsx")

markers$broad <- build_broad_marker_config(
  marker_list    = markers$broad,
  priority_order = c("vasculature", "immune"),
  default_threshold = 0.25
)

# 2.) Load data and create SCE -------------------------------------------------

input_dir <- "~/Desktop/OneDrive - University of Leeds/analysis/Multiomics/Data/cellvoter_inputs"

patient_samples <- list.files(input_dir, full.names = TRUE)
patient <-  gsub("_raw_counts\\.rds", "", basename(patient_samples))
names(patient_samples) <- patient

outlabs <- vector(mode = "list", length = length(patient_samples))
names(outlabs) <- names(patient_samples)

patient_samples <- patient_samples[5:12]


for (i in seq_along(patient_samples)) {

  patient_id <- names(patient_samples)[i]

  cli::cli_alert_info("\nrunning sample:\t {patient_id}")

  sce <- create_sce(counts = patient_samples[i])

  sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)

  sce <- normalize_counts(sce)

  sce <- prepare_sce(sce, markers)

  results <- run_cellvoter(sce)

  consensus <- resolve_consensus_labels(
    label_list        = results$labels,
    method_names      = results$method_names,
    tie_breaker_names = results$tie_breaker_names,
    unassigned_label  = "unknown",
    allow_even_split  = FALSE,
    ordered_tiebreak  = TRUE
  )

  results$sce$cellVoteR_label  <- consensus$label
  results$sce$cellVoteR_method <- consensus$method


  outlabs[[patient_id]] <- data.frame(
    sample_id = rownames(colData(results$sce)),
    colData(results$sce)[c("cellVoteR_label", "cellVoteR_method")],
    row.names = NULL
  )

}

outlabs[["Walton101_Recurrent"]] <- data.frame(
  sample_id = rownames(colData(results$sce)),
  colData(results$sce)[c("cellVoteR_label", "cellVoteR_method")],
  row.names = NULL
)


saveRDS(outlabs, "~/Desktop/OneDrive - University of Leeds/analysis/Multiomics/Data/cellvoter_labels.rds")
