#' Auto-Detect QC Feature Groups
#'
#' Scans the rownames of a Seurat object to identify genes matching specific patterns
#' (e.g., mitochondrial or ribosomal) and formats them for use in quality control checks.
#'
#' @param object A Seurat object.
#' @param group_configs A named list of configurations. Each element must be a list containing:
#' \itemize{
#'   \item \code{pattern}: A regex string to match gene names (e.g., "^MT-").
#'   \item \code{max_pct}: The maximum allowed percentage for this group (used later in QC).
#' }
#' Defaults to mitochondrial (`^MT-`, 20%) and ribosomal (`^RP[SL]`, 50%) genes.
#'
#' @return A named list where each element contains:
#' \itemize{
#'   \item \code{features}: A character vector of matching gene names found in the object.
#'   \item \code{max_pct}: The threshold percentage carried over from the config.
#' }
#' Returns an empty list structure for groups where no genes were found (with a CLI info message).
#'
#' @examples
#' \dontrun{
#' # Default usage
#' qc_feats <- find_qc_features(seu)
#'
#' # Custom patterns (e.g., hemoglobin)
#' custom_config <- list(
#'   hb = list(pattern = "^HB[A|B]", max_pct = 5)
#' )
#' qc_feats <- find_qc_features(seu, group_configs = custom_config)
#' }
#'
#' @export
find_qc_features <- function(
    object,
    group_configs = list(
      mito = list(pattern = "^MT-", max_pct = 20),
      ribo = list(pattern = "^RP[SL]", max_pct = 50)
    )
) {

  .valid_sObj_input(object)

  all_genes <- unique(rownames(object))

  final_list <- list()

  for (name in names(group_configs)) {
    config <- group_configs[[name]]

    found_genes <- grep(config$pattern, all_genes, value = TRUE, ignore.case = TRUE)

    if (length(found_genes) == 0) {
      cli::cli_alert_info("No genes found matching pattern {.val {config$pattern}} for group {.val {name}}.")
    }

    final_list[[name]] <- list(
      features = found_genes,
      max_pct = config$max_pct
    )
  }

  return(final_list)
}


#' Assess cell quality metrics for single-cell data
#'
#' Calculates comprehensive quality control (QC) metrics for a Seurat object,
#' including feature counts, library size, and specific gene group percentages
#' (e.g., mitochondrial or ribosomal genes). Supports filtering based on minimum
#' features, adaptive thresholds for library size, and user-defined feature groups.
#'
#' @param object A Seurat object containing single-cell data.
#' @param check_feature_groups A named list defining gene groups to check. Each element
#'   should be a list with `features` (vector of gene names) and `max_pct`
#'   (numeric threshold for maximum percentage). If `NULL` (default), attempts to
#'   auto-detect mitochondrial (`^MT-`) and ribosomal (`^RP[SL]`) genes using
#'   \code{\link{find_qc_features}}.
#' @param min_cells_per_sample Integer. Minimum number of cells required per sample
#'   for that sample's cells to be retained (default: 100). Requires `sample_col` to operate.
#' @param min_features Integer. Minimum number of detected features (genes) required
#'   for a cell to pass QC (default: 100).
#' @param nmads Numeric. Number of median absolute deviations (MADs) below the median
#'   log-transformed library size to set the adaptive low-count threshold (default: 3).
#'   Used to generate the `pass_adaptive_counts` metric but does not strictly filter
#'   unless manually enforced downstream.
#' @param sample_col Character vector or factor of the same length as the number of cells,
#'   indicating the sample ID for each cell. If `NULL`, all cells are treated as a single sample.
#' @param layer Character. The assay layer to use for calculating counts (default: "counts").
#' @param remove_failed_cells Logical. If `TRUE`, subsets the returned object to keep only
#'   cells that passed all QC checks (`QC_PASS == TRUE`). Default is `FALSE`.
#'
#' @return A Seurat object with updated metadata columns:
#' \itemize{
#'   \item \code{sc_n_features}: Number of features detected per cell.
#'   \item \code{sc_n_counts}: Total UMI counts per cell.
#'   \item \code{pass_min_features}: Logical, TRUE if cell exceeds `min_features`.
#'   \item \code{sc_adaptive_threshold}: The calculated low-count threshold for that cell.
#'   \item \code{pass_adaptive_counts}: Logical, TRUE if cell counts exceed adaptive threshold.
#'   \item \code{pct_[group]}: Percentage of counts belonging to specific feature groups (e.g., `pct_mito`).
#'   \item \code{pass_max_[group]}: Logical, TRUE if percentage is below `max_pct` for that group.
#'   \item \code{QC_PASS}: Logical, TRUE if the cell passed ALL active checks (features, groups, sample size).
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with auto-detection
#' seu <- assess_cell_quality(seu, remove_failed_cells = TRUE)
#'
#' # Custom usage with specific groups
#' groups <- list(
#'   mito = list(features = c("MT-ND1", "MT-CO1"), max_pct = 15),
#'   stress = list(features = c("JUN", "FOS"), max_pct = 5)
#' )
#' seu <- assess_cell_quality(seu, check_feature_groups = groups, min_features = 200)
#' }
#'
#' @importFrom Seurat GetAssayData AddMetaData
#' @importFrom Matrix colSums
#' @importFrom stats median mad
#' @export
assess_cell_quality <- function(
    object,
    check_feature_groups = NULL,
    min_cells_per_sample = 100,
    min_features = 100,
    nmads = 3,
    sample_col = NULL,
    layer = "counts",
    remove_failed_cells = FALSE
) {

  # Validation & Set-up ----

  .valid_sObj_input(object)

  cli::cli_h1("Assessing Cell Quality")

  if (is.null(check_feature_groups)) {
    cli::cli_alert_info("No feature groups provided. Auto-detecting Mito/Ribo patterns...")
    check_feature_groups <- .hush({find_qc_features(object)})
  }

  counts <- Seurat::GetAssayData(object, layer = layer)
  n_cells <- ncol(object)

  # Basic Metrics ----

  n_features <- Matrix::colSums(counts > 0)
  total_counts <- Matrix::colSums(counts)

  pass_mask <- n_features >= min_features

  log_counts <- log1p(total_counts)
  log_med <- median(log_counts)
  log_mad <- mad(log_counts)
  log_cutoff <- log_med - (nmads * log_mad)
  raw_cutoff <- round(expm1(log_cutoff), 2)
  pass_adaptive <- total_counts >= raw_cutoff


  meta_update <- data.frame(
    sc_n_features = n_features,
    sc_n_counts = total_counts,
    pass_min_features = pass_mask,
    # Adaptive columns (Not used in Vote)
    sc_adaptive_threshold = rep(raw_cutoff, n_cells),
    pass_adaptive_counts = pass_adaptive,
    row.names = colnames(object)
  )


  # Feature Group Checks ----

  for (group_name in names(check_feature_groups)) {

    config <- check_feature_groups[[group_name]]

    valid_genes <- intersect(config$features, rownames(counts))

    if (length(valid_genes) > 0) {

      group_sum <- Matrix::colSums(counts[valid_genes, , drop = FALSE])
      pct <- (group_sum / total_counts) * 100
      pct[is.nan(pct)] <- 0

    } else pct <- rep(0, n_cells)

    group_pass <- pct <= config$max_pct
    pass_mask <- pass_mask & group_pass

    meta_update[[paste0("pct_", group_name)]] <- pct
    meta_update[[paste0("pass_max_", group_name)]] <- group_pass

  }

  # Sample Size Check ----

  if (!is.null(sample_col) && length(sample_col) == n_cells) {

    sample_ids <- as.character(sample_col)

  } else {

    if (!is.null(sample_col)) {
      cli::cli_alert_warning("Invalid {.arg sample_col} length. Treating all cells as one sample.")
    }

    sample_ids <- rep("All_Cells", n_cells)
  }

  sample_counts <- table(sample_ids)
  valid_samples <- names(sample_counts)[sample_counts >= min_cells_per_sample]

  sample_pass <- sample_ids %in% valid_samples
  pass_mask <- pass_mask & sample_pass

  meta_update$sc_sample_id <- sample_ids
  meta_update$pass_min_cells_sample <- sample_pass
  meta_update$QC_PASS <- pass_mask

  # Return ----

  object <- Seurat::AddMetaData(
    object = object,
    metadata = meta_update,
    col.name = colnames(meta_update)
  )

  n_pass <- sum(pass_mask)
  cli::cli_alert_success("QC Complete: {n_pass}/{n_cells} cells passed ({round(n_pass/n_cells*100, 1)}%)")

  if(remove_failed_cells & n_pass < n_cells) {

    object <- object[, pass_mask]
    cli::cli_alert_warning("Removed failed cells from the Seurat object.")

  }

  return(object)
}
