#' Filter cells based on QC metrics
#'
#' Removes cells that failed QC (QC_PASS == FALSE) and prints a summary of the operation.
#'
#' @param object The raw single-cell object (Seurat, SingleCellExperiment, or Matrix).
#' @param qc_metadata The dataframe output from \code{assess_cell_quality}.
#' @param verbose Logical. If TRUE, prints a summary of removed cells.
#' @param show_by_sample Logical. If TRUE and 'sample_id' is present in metadata, prints removal stats per sample.
#'
#' @return The filtered object (same class as input).
#' @export
#'
#' @examples
#' \dontrun{
#' qc <- assess_cell_quality(counts)
#' clean_data <- filter_cells(counts, qc)
#' }
filter_cells <- function(object, qc_metadata, verbose = TRUE, show_by_sample = FALSE) {

  # 1. Validation
  # Ensure QC metadata aligns with the object
  # We check if the number of rows in metadata matches columns in object
  obj_n_cells <- ncol(object)
  qc_n_cells <- nrow(qc_metadata)

  if (obj_n_cells != qc_n_cells) {
    stop(sprintf("Dimension mismatch: Object has %d cells but QC metadata has %d rows.",
                 obj_n_cells, qc_n_cells))
  }

  if (!"QC_PASS" %in% colnames(qc_metadata)) {
    stop("qc_metadata must contain a 'QC_PASS' column (output from assess_cell_quality).")
  }

  # 2. Identify Keepers
  # Use cell_id if strictly safer, otherwise rely on logical indexing
  keep_mask <- qc_metadata$QC_PASS
  n_start <- nrow(qc_metadata)
  n_keep <- sum(keep_mask)
  n_removed <- n_start - n_keep

  # 3. Subset the Object
  # Standard subsetting [ features, cells ] works for Matrix, Seurat, and SCE
  filtered_object <- object[, keep_mask]

  # 4. Reporting
  if (verbose) {
    pct_removed <- round((n_removed / n_start) * 100, 1)
    message(sprintf("Global Filter Summary:\n- Input cells: %d\n- Removed: %d (%s%%)\n- Remaining: %d",
                    n_start, n_removed, pct_removed, n_keep))

    if (show_by_sample && "sample_id" %in% colnames(qc_metadata)) {
      message("\nBreakdown by Sample:")

      # Create a summary table
      # We group by sample_id and count Total vs Removed
      samples <- unique(qc_metadata$sample_id)

      for (s in samples) {
        # Subset metadata for this sample
        sub_meta <- qc_metadata[qc_metadata$sample_id == s, ]
        n_s_total <- nrow(sub_meta)
        n_s_rem <- sum(!sub_meta$QC_PASS)
        pct_s <- round((n_s_rem / n_s_total) * 100, 1)

        message(sprintf("  - %s: Removed %d/%d (%s%%)", s, n_s_rem, n_s_total, pct_s))
      }
    }
  }

  return(filtered_object)
}
