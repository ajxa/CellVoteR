#' Run the CellVoteR Ensemble Pipeline
#'
#' @param object A Seurat object.
#' @param markers User markers (or NULL to use internal).
#' @param panel_name Name of internal panel to use if markers=NULL (e.g. "idhwt_gbm_markers").
#' @param triage_markers Named list for broad categories.
#' @param use_ensemble Logical. If TRUE, runs specific triage methods (1-4) in addition to globals.
#' @param verbose Logical.
#' @return A Seurat object with 'Ensemble_Resolved' labels.
#' @export
run_ensemble <- function(object,
                         markers = NULL,
                         panel_name = "idhwt_gbm_markers",
                         triage_markers = NULL,
                         use_ensemble = FALSE,
                         verbose = TRUE) {

  object <- ensure_seurat(object)

  # --- 1. Process Markers & Defaults ---
  valid_markers <- process_marker_input(markers, panel_name = panel_name)

  if (is.null(triage_markers)) {
    triage_markers <- if (exists("cellvoter_data")) cellvoter_data$cell_groups else list(Immune="PTPRC", Endothelial=c("CDH5","VWF"))
  }

  # --- 2. Define Strategy ---
  if (use_ensemble) {
    methods_to_use <- 1:4; breaker_methods <- 5:6
    if (verbose) message("Running Full Ensemble (1-4) with Global Breakers (5-6)...")
  } else {
    methods_to_use <- 5:6; breaker_methods <- NULL
    if (verbose) message("Running General Purpose Labelling (5-6)...")
  }

  if (!"data" %in% names(object@assays$RNA)) object <- Seurat::NormalizeData(object, verbose = FALSE)

  # --- Phase 1: Primary Methods ---
  vote_matrix <- matrix(NA, nrow = ncol(object), ncol = length(methods_to_use))
  colnames(vote_matrix) <- paste0("Method_", methods_to_use)

  for (i in seq_along(methods_to_use)) {
    m <- methods_to_use[i]
    if (verbose) message(sprintf("Running Method %d...", m))

    # Pass markers and triage down
    res <- run_cellvoter(object, method = m, markers = valid_markers,
                         triage_markers = triage_markers, verbose = FALSE)

    object[[paste0("Method_", m, "_Label")]] <- res$object$CellVoteR_Label
    vote_matrix[, i] <- res$object$CellVoteR_Label
  }

  # --- Phase 2: Breakers (if needed) ---
  breaker_votes <- list()
  if (use_ensemble) {
    for (m in breaker_methods) {
      if (verbose) message(sprintf("Running Breaker Method %d...", m))
      res <- run_cellvoter(object, method = m, markers = valid_markers,
                           triage_markers = triage_markers, verbose = FALSE)
      breaker_votes[[paste0("Method_", m)]] <- res$object$CellVoteR_Label
      object[[paste0("Method_", m, "_Label")]] <- res$object$CellVoteR_Label
    }
  }

  # --- Phase 3: Consensus ---
  raw_ensemble <- apply(vote_matrix, 1, function(row) {
    valid_votes <- row[!is.na(row) & !grepl("Unassigned", row)]
    if (length(valid_votes) == 0) return("Unknown")
    counts <- table(valid_votes)
    winners <- names(counts)[counts == max(counts)]
    if (length(winners) == 1) return(winners)
    return(paste(sort(winners), collapse = ":"))
  })

  object$Ensemble_Raw <- raw_ensemble

  # --- Phase 4: Resolution ---
  resolved_labels <- raw_ensemble
  clash_indices <- grep(":", resolved_labels)

  # A. Use Breakers
  if (use_ensemble && length(clash_indices) > 0) {
    for (i in clash_indices) {
      candidates <- strsplit(resolved_labels[i], ":")[[1]]
      b_votes <- c(breaker_votes$Method_5[i], breaker_votes$Method_6[i])
      b_votes <- b_votes[b_votes %in% candidates]
      if (length(b_votes) > 0) resolved_labels[i] <- names(sort(table(b_votes), decreasing = TRUE))[1]
    }
  }

  # B. Use Marker Intensity
  remaining_clashes <- grep(":", resolved_labels)
  if (length(remaining_clashes) > 0) {
    intens_resolved <- resolve_clashes_internal(
      object,
      colnames(object)[remaining_clashes],
      resolved_labels[remaining_clashes],
      valid_markers # Use valid_markers, NOT hardcoded marker_ref
    )
    resolved_labels[remaining_clashes] <- intens_resolved
  }

  object$Ensemble_Resolved <- resolved_labels
  return(object)
}

#' @noRd
resolve_clashes_internal <- function(object, cells, current_labels, ref) {
  expr_mat <- Seurat::GetAssayData(object, slot = "data")[, cells, drop = FALSE]
  resolved <- character(length(cells))
  names(resolved) <- cells

  for (cell in cells) {
    candidates <- strsplit(current_labels[cell], ":")[[1]]
    scores <- vapply(candidates, function(ct) {
      markers <- intersect(ref$gene[ref$cell_type == ct], rownames(expr_mat))
      if (length(markers) == 0) return(0)
      vals <- expr_mat[markers, cell]; expressed_vals <- vals[vals > 0]
      if (length(expressed_vals) == 0) return(0)
      return(mean(expressed_vals))
    }, FUN.VALUE = numeric(1))

    if (all(scores == 0)) resolved[cell] <- "Unknown" else resolved[cell] <- names(which.max(scores))
  }
  return(resolved)
}
