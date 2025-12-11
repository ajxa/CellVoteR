#' Run the CellVoteR Labelling Method
#'
#' @param object Input Seurat object.
#' @param method Integer (1-6). 5 & 6 are Global (Default). 1-4 are Triage.
#' @param markers User markers (or NULL to use internal).
#' @param panel_name If markers is NULL, which internal panel to use? Default "idhwt_gbm_markers".
#' @param triage_markers Named list of broad markers. If NULL, uses internal defaults.
#' @param verbose Logical.
#' @param plot_heatmaps Logical.
#' @return List(object, heatmaps, scores)
#' @export
run_cellvoter <- function(object, method = 5,
                          markers = NULL,
                          panel_name = "idhwt_gbm_markers",
                          triage_markers = NULL,
                          verbose = TRUE, plot_heatmaps = FALSE) {

  if (!method %in% 1:6) stop("Method must be integer 1-6.")

  # --- 1. Process Markers ---
  marker_ref <- process_marker_input(markers, panel_name = panel_name)

  # --- 2. Handle Triage Defaults ---
  if (is.null(triage_markers)) {
    if (exists("cellvoter_data")) {
      triage_markers <- cellvoter_data$cell_groups
    } else {
      warning("Internal data missing. Using hardcoded fallback.")
      triage_markers <- list(Immune = "PTPRC", Endothelial = c("CDH5", "VWF"))
    }
  }

  # --- Step A: Data Prep ---
  if (method %in% c(3, 4, 6)) {
    if (verbose) message("Step 1: Filtering data to reference genes...")
    object <- restrict_to_reference(object, unique(marker_ref$gene))
  }

  # --- Step B: Broad Labelling ---
  if (method %in% 1:4) {
    if (method %in% c(1, 3)) {
      if (verbose) message("Step 2: Broad labelling via CLUSTERING (Triage)...")
      object <- run_broad_labelling(object, broad_markers = triage_markers)
    } else {
      if (verbose) message("Step 2: Broad labelling via MARKER THRESHOLDS (Triage)...")
      object$Broad_Type <- broad_label_by_cell(object, marker_list = triage_markers)
    }
  } else {
    if (verbose) message("Step 2: Skipping Broad Labelling (Global Method)...")
    object$Broad_Type <- "Global"
  }

  # --- Step C: Sub-clustering & Scoring ---
  final_labels <- rep("Unassigned", ncol(object))
  names(final_labels) <- colnames(object)
  heatmap_list <- list(); score_list <- list()
  background_genes <- rownames(object)

  for (b_type in unique(object$Broad_Type)) {
    if (verbose) message(sprintf("  - Processing group: %s", b_type))

    cells_in_group <- colnames(object)[object$Broad_Type == b_type]
    if (length(cells_in_group) == 0) next
    sub_obj <- subset(object, cells = cells_in_group)
    res <- process_subcluster(sub_obj)

    if (is.null(res)) {
      final_labels[cells_in_group] <- b_type
      next
    }

    # Filter Reference based on Broad Type
    if (method %in% 5:6) {
      relevant_ref <- marker_ref
    } else {
      valid_categories <- get_category_map(b_type)
      # Check if 'category' column exists (it should from process_marker_input)
      if ("category" %in% colnames(marker_ref)) {
        relevant_ref <- marker_ref[marker_ref$category %in% valid_categories, ]
      } else {
        relevant_ref <- marker_ref
      }
    }

    if (nrow(relevant_ref) == 0) {
      final_labels[cells_in_group] <- b_type; next
    }

    scores <- calculate_fisher_scores(res$markers, relevant_ref, background_genes)
    score_list[[b_type]] <- scores
    if (plot_heatmaps) heatmap_list[[b_type]] <- plot_fisher_heatmap(scores)

    cluster_map <- get_best_match(scores)
    new_labels <- cluster_map[as.character(Seurat::Idents(res$object))]
    new_labels[is.na(new_labels)] <- paste0(b_type, "_Unassigned")
    final_labels[names(new_labels)] <- new_labels
  }

  object$CellVoteR_Label <- final_labels
  return(list(object = object, heatmaps = heatmap_list, scores = score_list))
}

#' @noRd
get_category_map <- function(broad_type) {
  if (broad_type == "Immune") return("Immune")
  if (broad_type %in% c("Endothelial", "Vasculature")) return(c("Vasculature", "Endothelial"))
  # Everything else matches to Malignant/Normal/etc
  return(c("Cancer", "Normal", "Immune", "Vasculature", "Endothelial", "Malignant"))
}
