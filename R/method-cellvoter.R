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

  marker_ref <- process_marker_input(markers, panel_name = panel_name)

  triage_markers <- process_triage_input(triage_markers)

  if (method %in% c(3, 4, 6)) {
    if (verbose) message("Step 1: Filtering data to reference genes...")
    object <- restrict_to_reference(object, unique(marker_ref$gene))
  }

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

  # Sub-clustering & scoring
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

    # --- FILTERING LOGIC (Dynamic) ---
    if (method %in% 5:6) {
      # Global methods use everything
      relevant_ref <- marker_ref
    } else {
      # Triage methods (1-4)
      if (b_type == "Other") {
        # "Other" means it did NOT match any of the Triage keys.
        # So we include all reference categories EXCEPT the triage keys.
        exclude_cats <- names(triage_markers)

        # Alias Fix: If we exclude "Endothelial", also exclude "Vasculature"
        # (Handles the mismatch in the default GBM panel names)
        if ("Endothelial" %in% exclude_cats) exclude_cats <- c(exclude_cats, "Vasculature")

        relevant_ref <- marker_ref[!marker_ref$category %in% exclude_cats, ]
      } else {
        # Specific Match: Filter for the category matching the Broad Type
        relevant_ref <- marker_ref[marker_ref$category == b_type, ]

        # Alias Fix: If Broad Type is "Endothelial" but ref uses "Vasculature"
        if (nrow(relevant_ref) == 0 && b_type == "Endothelial") {
          relevant_ref <- marker_ref[marker_ref$category == "Vasculature", ]
        }
      }
    }

    if (nrow(relevant_ref) == 0) {
      # If no matching markers found (e.g. user provided "GroupA" but ref has no "GroupA"),
      # default to Broad Label (better than crashing or empty)
      final_labels[cells_in_group] <- b_type
      next
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

