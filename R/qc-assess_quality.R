#' Auto-detect QC feature groups
#'
#' Scans the row names of a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object to identify
#' genes matching specific patterns (e.g. mitochondrial or ribosomal) and
#' formats them for use in downstream quality control filtering.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @param group_configs A named list of configurations. Each element must be a
#'   list containing:
#'   \describe{
#'     \item{pattern}{A regex string to match gene names (e.g. \code{"^MT-"}).}
#'     \item{max_pct}{Numeric scalar. The maximum allowed percentage for this
#'       group, carried forward for use in QC filtering.}
#'   }
#'  Defaults to mitochondrial (`^MT-`, 20%) and ribosomal (`^RP[SL]`, 50%) genes.
#'
#' @return A named list where each element contains:
#'   \describe{
#'     \item{features}{Character vector of matching gene names found in the
#'       object. Empty character vector if no genes matched.}
#'     \item{max_pct}{Numeric threshold carried over from the config.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Default usage
#' qc_feats <- find_qc_features(sce)
#'
#' # Custom patterns (e.g. hemoglobin genes)
#' qc_feats <- find_qc_features(sce, group_configs = list(
#'   hb = list(pattern = "^HB[AB]", max_pct = 5)
#' ))
#' }
#'
#' @export
find_qc_features <- function(
    sce,
    group_configs = list(
      mito = list(pattern = "^MT-", max_pct = 20),
      ribo = list(pattern = "^RP[SL]", max_pct = 50)
      )
) {

  if (!inherits(sce, "SingleCellExperiment")) cli::cli_abort("{.arg sce} must be a {.cls SingleCellExperiment}. Received: {.cls {class(sce)}}")

  if (!is.list(group_configs) || is.null(names(group_configs))) cli::cli_abort("{.arg group_configs} must be a named list.")

  required_fields <- c("pattern", "max_pct")

  for (name in names(group_configs)) {

    config <- group_configs[[name]]

    if (!is.list(config)) cli::cli_abort("Entry {.val {name}} in {.arg group_configs} must be a list.")

    missing_fields <- setdiff(required_fields, names(config))

    if (length(missing_fields) > 0L) {

      cli::cli_abort(
        "Entry {.val {name}} in {.arg group_configs} is missing: {.val {missing_fields}}"
      )

    }

    if (!is.character(config$pattern) || length(config$pattern) != 1L) {

      cli::cli_abort("Field {.val pattern} in group {.val {name}} must be a single character string.")

    }

    if (!is.numeric(config$max_pct) || length(config$max_pct) != 1L || config$max_pct < 0) {

      cli::cli_abort("Field {.val max_pct} in group {.val {name}} must be a non-negative number.")

    }

  }

  all_genes  <- rownames(sce)
  final_list <- list()

  for (name in names(group_configs)) {

    config       <-   group_configs[[name]]
    found_genes  <-   grep(config$pattern, all_genes, value = TRUE, ignore.case = TRUE)

    if (length(found_genes) == 0L) {

      print_alert(
        "No genes matching pattern {.val {config$pattern}} for group {.val {name}}",
        type = "w"
      )

    } else {

      print_alert(
        "Found {.val {length(found_genes)}} gene{?s} for group {.val {name}}",
        type = "s"
      )
    }

    final_list[[name]] <- list(features = found_genes, max_pct  = config$max_pct)

    final_list[[name]] <- list(
      features = found_genes,
      max_pct  = config$max_pct
    )

  }

  return(final_list)

}


#' Resolve sample identifiers from colData or a supplied vector
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param sample_col A column name (character scalar) referencing
#'   \code{colData(sce)}, a character/factor vector of length \code{n_cells},
#'   or \code{NULL}.
#' @param n_cells Integer. Expected number of cells.
#' @return Character vector of sample identifiers.
#' @keywords internal
#' @noRd
resolve_sample_ids <- function(sce, sample_col, n_cells) {

  if (is.null(sample_col)) return(rep("all_cells", n_cells))

  if (is.character(sample_col) && length(sample_col) == 1L) {

    cd <- SummarizedExperiment::colData(sce)

    if (!sample_col %in% colnames(cd)) {

      cli::cli_abort(
        "Column {.val {sample_col}} not found in {.code colData(sce)}.
         Available: {.val {colnames(cd)}}"
      )

    }

    return(as.character(cd[[sample_col]]))

  }

  if (length(sample_col) != n_cells) {

    print_alert("{.code length(sample_col)} ({.val {length(sample_col)}}) does not match
                number of cells ({.val {n_cells}}).")

    print_alert("Treating all cells as a single sample.", type = "i")

    return(rep("all_cells", n_cells))
  }

  return(as.character(sample_col))

}



#' Assess Cell Quality Metrics
#'
#' Calculates comprehensive quality control (QC) metrics for a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object, including
#' feature counts, library size, and specific gene group percentages (e.g.
#' mitochondrial or ribosomal genes). Supports filtering based on minimum
#' features, adaptive thresholds for library size, and user-defined feature
#' groups.
#'
#' @section QC strategy:
#' Each cell is assessed against multiple independent criteria. A cell must
#' pass \strong{all} active checks to receive \code{QC_PASS = TRUE}:
#' \enumerate{
#'   \item Minimum detected features (\code{min_features}).
#'   \item Feature group percentage thresholds (e.g. mitochondrial < 20\%).
#'   \item Minimum cells per sample (\code{min_cells_per_sample}).
#' }
#' An adaptive library-size threshold is also computed and stored for
#' reference but is \strong{not} included in the pass/fail vote by default.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'   with a \code{counts} assay.
#' @param check_feature_groups A named list defining gene groups to check.
#'   Each element must be a list with \code{features} (character vector of
#'   gene names) and \code{max_pct} (numeric threshold for maximum
#'   percentage). If \code{NULL} (default), auto-detects mitochondrial
#'   (\code{^MT-}) and ribosomal (\code{^RP[SL]}) genes via
#'   \code{\link{find_qc_features}}.
#' @param min_cells_per_sample Integer scalar. Minimum number of cells
#'   required per sample for that sample's cells to be retained. Requires
#'   \code{sample_col} to operate. Defaults to \code{100}.
#' @param min_features Integer scalar. Minimum number of detected genes
#'   required for a cell to pass QC. Defaults to \code{100}.
#' @param nmads Numeric scalar. Number of median absolute deviations below
#'   the median log-transformed library size to set the adaptive low-count
#'   threshold. Stored for reference but not included in the pass/fail vote.
#'   Defaults to \code{3}.
#' @param sample_col Character scalar naming a column in \code{colData(sce)}
#'   containing sample identifiers, or a character/factor vector of the same
#'   length as the number of cells. If \code{NULL}, all cells are treated as
#'   a single sample.
#' @param assay_name Character scalar. Name of the assay to use. Defaults to
#'   \code{"counts"}.
#' @param remove_failed_cells Logical. If \code{TRUE}, subsets the returned
#'   object to retain only cells with \code{QC_PASS == TRUE}. Defaults to
#'   \code{FALSE}.
#'
#' @return The input SCE with additional columns in \code{colData()}:
#'   \describe{
#'     \item{\code{sc_n_features}}{Number of detected genes per cell.}
#'     \item{\code{sc_n_counts}}{Total UMI counts per cell.}
#'     \item{\code{pass_min_features}}{Logical. \code{TRUE} if cell exceeds
#'       \code{min_features}.}
#'     \item{\code{sc_adaptive_threshold}}{The computed low-count threshold.}
#'     \item{\code{pass_adaptive_counts}}{Logical. \code{TRUE} if cell counts
#'       exceed the adaptive threshold.}
#'     \item{\code{pct_[group]}}{Percentage of counts from each feature group.}
#'     \item{\code{pass_max_[group]}}{Logical. \code{TRUE} if percentage is
#'       below \code{max_pct} for that group.}
#'     \item{\code{sc_sample_id}}{Sample identity for each cell.}
#'     \item{\code{pass_min_cells_sample}}{Logical. \code{TRUE} if cell
#'       belongs to a sample meeting \code{min_cells_per_sample}.}
#'     \item{\code{QC_PASS}}{Logical. \code{TRUE} if the cell passed all
#'       active checks.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage with auto-detection
#' sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)
#'
#' # Custom feature groups
#' groups <- list(
#'   mito   = list(features = c("MT-ND1", "MT-CO1"), max_pct = 15),
#'   stress = list(features = c("JUN", "FOS"), max_pct = 5)
#' )
#' sce <- assess_cell_quality(sce,
#'   check_feature_groups = groups,
#'   min_features = 200
#' )
#'
#' # With per-sample filtering using a colData column
#' sce <- assess_cell_quality(sce, sample_col = "patient_id")
#' }
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay assayNames colData colData<-
#' @importFrom Matrix colSums
#' @importFrom stats median mad
#' @export
assess_cell_quality <- function(
    sce,
    check_feature_groups = NULL,
    min_cells_per_sample = 100,
    min_features = 100,
    nmads = 3,
    sample_col = NULL,
    assay_name = "counts",
    remove_failed_cells = FALSE
) {

  # -- Input validation --------------------------------------------------------

  if (!inherits(sce, "SingleCellExperiment")) cli::cli_abort("{.arg sce} must be a {.cls SingleCellExperiment}. Received: {.cls {class(sce)}}")

  if (!assay_name %in% SummarizedExperiment::assayNames(sce)) {
    cli::cli_abort(
      "Assay {.val {assay_name}} not found. Available: {.val {SummarizedExperiment::assayNames(sce)}}"
    )
  }

  checkmate::assert_count(min_cells_per_sample, positive = TRUE)
  checkmate::assert_count(min_features, positive = TRUE)
  checkmate::assert_number(nmads, lower = 0)
  checkmate::assert_flag(remove_failed_cells)

  if (is.null(check_feature_groups)) {
    print_alert("Auto-detecting mito/ribo patterns...")
    check_feature_groups <- find_qc_features(sce)
  }

  counts  <- SummarizedExperiment::assay(sce, assay_name)
  n_cells <- ncol(sce)

  n_features   <- Matrix::colSums(counts > 0)
  total_counts <- Matrix::colSums(counts)

  pass_mask <- n_features >= min_features

  # Adaptive library-size threshold (stored but not used for voting)
  log_counts <- log1p(total_counts)
  log_med    <- stats::median(log_counts)
  log_mad    <- stats::mad(log_counts)
  log_cutoff <- log_med - (nmads * log_mad)
  raw_cutoff <- round(expm1(log_cutoff), 2)

  pass_adaptive <- total_counts >= raw_cutoff

  meta_update <- S4Vectors::DataFrame(
    sc_n_features         = n_features,
    sc_n_counts           = total_counts,
    pass_min_features     = pass_mask,
    sc_adaptive_threshold = rep(raw_cutoff, n_cells),
    pass_adaptive_counts  = pass_adaptive,
    row.names             = colnames(sce)
  )

  for (group_name in names(check_feature_groups)) {

    config      <- check_feature_groups[[group_name]]
    valid_genes <- intersect(config$features, rownames(counts))

    if (length(valid_genes) > 0L) {

      group_sum <- Matrix::colSums(counts[valid_genes, , drop = FALSE])
      pct       <- (group_sum / total_counts) * 100
      pct[is.nan(pct)] <- 0

    } else pct <- rep(0, n_cells)

    group_pass <- pct <= config$max_pct
    pass_mask  <- pass_mask & group_pass

    meta_update[[paste0("pct_", group_name)]]      <- pct
    meta_update[[paste0("pass_max_", group_name)]]  <- group_pass
  }

  sample_ids <- resolve_sample_ids(sce, sample_col, n_cells)

  sample_counts  <- table(sample_ids)
  valid_samples  <- names(sample_counts)[sample_counts >= min_cells_per_sample]

  sample_pass <- sample_ids %in% valid_samples
  pass_mask   <- pass_mask & sample_pass

  meta_update$sc_sample_id           <- sample_ids
  meta_update$pass_min_cells_sample  <- sample_pass
  meta_update$QC_PASS                <- pass_mask

  # -- Attach QC metadata to colData -------------------------------------------

  existing_cd <- SummarizedExperiment::colData(sce)

  overlap <- intersect(colnames(existing_cd), colnames(meta_update))

  if (length(overlap) > 0L) {
    existing_cd <- existing_cd[, !colnames(existing_cd) %in% overlap, drop = FALSE]
  }

  SummarizedExperiment::colData(sce) <- cbind(existing_cd, meta_update)

  # -- (Optional) filtering ----------------------------------------------------

  n_pass <- sum(pass_mask)
  n_fail <- n_cells - n_pass

  if (remove_failed_cells && n_pass < n_cells) {
    sce <- sce[, pass_mask]
    print_alert("Removed {.val {n_fail}} cell{?s}", type = "w")
  }

  print_status_bar(
    n_pass       = round((n_pass / n_cells) * 100, 1),
    n_fail       = round((n_fail / n_cells) * 100, 1),
    label_text   = "passed",
    unknown_text = "failed",
    use_props    = TRUE
  )

  return(sce)
}
