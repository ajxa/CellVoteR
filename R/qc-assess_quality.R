#' Assess Cell Quality Metrics (Core)
#'
#' Calculates QC metrics for each cell within a single-cell matrix. It also
#' provides QC metrics for explicitly provided gene groups, e.g., mitochondrial
#' and ribosomal genes.
#'
#' @param object Input object (Seurat, SingleCellExperiment, or Matrix).
#' @param check_feature_groups Named list of gene sets. Each entry must contain:
#'   \itemize{
#'     \item \code{features}: Character vector of gene symbols present in the matrix.
#'     \item \code{max_pct}: Maximum allowed percentage (0-100).
#'   }
#'   Default is empty (no groups checked).
#' @param min_cells_per_sample Min cells per sample ID.
#' @param min_features Min detected genes per cell.
#' @param sample_ids Vector of sample IDs.
#'
#' @return A data.frame with QC metrics and pass/fail flags.
#' @export
assess_cell_quality <- function(object,
                                check_feature_groups = list(),
                                min_cells_per_sample = 100,
                                min_features = 200,
                                sample_ids = NULL) {
  # 1. Standardize Input
  counts <- extract_counts(object)
  total_counts <- Matrix::colSums(counts)
  n_features <- Matrix::colSums(counts > 0)

  # Initialise Output
  qc_df <- data.frame(
    cell_id = colnames(counts),
    n_features = n_features,
    n_counts = total_counts,
    pass_min_features = n_features >= min_features,
    stringsAsFactors = FALSE
  )

  vote_list <- list(features = qc_df$pass_min_features)
  missing_groups <- c()

  # 2. Iterate through explicit gene groups
  for (group_name in names(check_feature_groups)) {
    config <- check_feature_groups[[group_name]]

    # Strictly require 'features' vector
    target_genes <- config$features

    # Check if genes actually exist in the matrix (safety check)
    # We do NOT search for them, we just assume the user provided valid ones.
    # If the list is empty/NULL, we flag it.
    if (is.null(target_genes) || length(target_genes) == 0) {
      missing_groups <- c(missing_groups, group_name)
      qc_df[[paste0("pct_", group_name)]] <- 0
      qc_df[[paste0("pass_max_", group_name)]] <- TRUE
      next
    }

    # Calculate %
    # Subset matrix only on found genes
    # Note: If some genes in 'target_genes' are not in rownames(counts),
    # Matrix subsetting might fail or return NA depending on version.
    # Safe intersection:
    valid_genes <- intersect(target_genes, rownames(counts))

    if (length(valid_genes) == 0) {
      missing_groups <- c(missing_groups, group_name)
      qc_df[[paste0("pct_", group_name)]] <- 0
      qc_df[[paste0("pass_max_", group_name)]] <- TRUE
      next
    }

    group_counts <- Matrix::colSums(counts[valid_genes, , drop=FALSE])
    pct <- (group_counts / total_counts) * 100
    pct[is.nan(pct)] <- 0

    # Check Threshold
    pass_flag <- pct <= config$max_pct

    # Store
    qc_df[[paste0("pct_", group_name)]] <- pct
    qc_df[[paste0("pass_max_", group_name)]] <- pass_flag
    vote_list[[group_name]] <- pass_flag
  }

  # 3. Consolidated Warning
  if (length(missing_groups) > 0 && sum(total_counts) > 0) {
    warning("The following QC groups had no valid genes found in the matrix: ",
            paste(missing_groups, collapse = ", "))
  }

  # 4. Sample Size Check
  if (is.null(sample_ids)) sample_ids <- rep("Sample_1", ncol(counts))

  qc_df$sample_id <- sample_ids
  sample_counts <- table(sample_ids)
  valid_samples <- names(sample_counts)[sample_counts >= min_cells_per_sample]
  pass_sample <- sample_ids %in% valid_samples

  qc_df$pass_min_cells_sample <- pass_sample
  vote_list[["sample_size"]] <- pass_sample

  # 5. Final Vote
  qc_df$QC_PASS <- Reduce(`&`, vote_list)
  rownames(qc_df) <- colnames(counts)

  return(qc_df)
}
