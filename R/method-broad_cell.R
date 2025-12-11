#' Broad Labelling by Cell (Marker Thresholding)
#'
#' @param object Seurat object or matrix.
#' @param marker_list Named list of markers.
#' @param expr_threshold Threshold for expression (0.1).
#' @param min_markers Min positive markers required (1).
#'
#' @return Vector of labels.
#' @export
#' @importFrom Matrix colSums
broad_label_by_cell <- function(object,
                                marker_list = NULL,
                                expr_threshold = 0.1,
                                min_markers = 1) {

  counts <- extract_counts(object)

  if (is.null(marker_list)) {
    marker_list <- list(Immune = c("PTPRC"), Endothelial = c("CDH5", "VWF"))
  }

  cell_names <- colnames(counts)
  results <- matrix(FALSE, nrow = length(cell_names), ncol = length(marker_list))
  colnames(results) <- names(marker_list)
  rownames(results) <- cell_names

  for (type in names(marker_list)) {
    markers <- intersect(marker_list[[type]], rownames(counts))
    if (length(markers) > 0) {
      sub_mat <- counts[markers, , drop = FALSE]
      detected <- Matrix::colSums(sub_mat > expr_threshold)
      results[, type] <- detected >= min_markers
    }
  }

  labels <- rep("Other", length(cell_names))
  if ("Endothelial" %in% colnames(results)) labels[results[, "Endothelial"]] <- "Endothelial"
  if ("Immune" %in% colnames(results)) labels[results[, "Immune"]] <- "Immune"

  return(labels)
}
