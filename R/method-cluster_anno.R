#' Annotate clusters based on marker expression
#'
#' Assigns biological labels to existing numerical clusters by running differential expression
#' (DE) analysis and matching the resulting markers against a provided reference list.
#'
#' @details
#' \strong{Workflow:}
#' \enumerate{
#'   \item \strong{DE Analysis:} Runs \code{Seurat::FindAllMarkers} on the \code{initial_cluster_col}
#'   to identify defining features for each numeric cluster.
#'   \item \strong{Matching:} Filters the found markers against the \code{anno_markers} reference list.
#'   \item \strong{Scoring:} For each cluster, calculates a score for every potential annotation type.
#'   The score is the \strong{median avg_log2FC} of the matched markers.
#'   \item \strong{Assignment:} Assigns the annotation with the highest score to the cluster.
#'   Clusters with no matching markers in the reference are labelled as \code{unnassigned_lab}.
#' }
#'
#' @param sObj A Seurat object containing the data and existing cluster IDs.
#' @param anno_markers A named list of marker genes (e.g., \code{list(T_cells = c("CD3D", "CD3E"))})
#'   or a configuration object created by \code{build_broad_marker_config}.
#' @param initial_cluster_col Character. The metadata column name containing the current
#'   (usually numerical) cluster IDs.
#' @param label_out_col Character. The name of the new metadata column where the
#'   biological labels will be stored.
#' @param unnassigned_lab Character. Label to assign if a cluster does not match any
#'   reference markers. Default is "Other".
#' @param logfc.threshold Numeric. Minimum log2 fold-change required for a gene to be
#'   considered a marker. Passed to \code{FindAllMarkers}. Default is 0.1.
#' @param min.pct Numeric. Minimum percentage of cells expressing the gene. Passed to \code{FindAllMarkers}.
#' @param only.pos Logical. If \code{TRUE}, only returns positive markers. Default is \code{TRUE}.
#' @param p_val_adj_threshold Numeric. Maximum adjusted p-value for a marker to be kept
#'   before matching. Default is 0.1.
#' @param p_adjust_method Character. Method for p-value correction (e.g., "BH").
#' @param ... Additional arguments passed to \code{Seurat::FindAllMarkers}.
#'
#' @return A Seurat object with a new metadata column (\code{label_out_col}) containing
#' the biological annotations.
#'
#' @importFrom stats median setNames
#' @importFrom dplyr .data
#' @export
annotate_broad_clusters  <- function(
    sObj,
    anno_markers,
    initial_cluster_col,
    label_out_col,
    unnassigned_lab = "Other",
    logfc.threshold = 0.1,
    min.pct = 0.1,
    only.pos = TRUE,
    p_val_adj_threshold = 0.1,
    p_adjust_method = "BH",
    ...
) {

  .check_missing_args(
    required_args = c("sObj", "anno_markers", "initial_cluster_col", "label_out_col")
  )

  broad_marker_struct <- .valid_broad_marker_struct(anno_markers)
  if(broad_marker_struct) anno_markers <- lapply(anno_markers, `[[`, "markers")

  if(!(initial_cluster_col %in% colnames(sObj@meta.data))) {
    stop(paste0("The specified 'initial_cluster_col' (", initial_cluster_col, ") is not found in the Seurat object's metadata."))
  }

  cli::cli_progress_step("Finding cluster markers...", spinner = TRUE)

  cluster_markers <- Seurat::FindAllMarkers(
    object = sObj,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct,
    only.pos = only.pos,
    p.adjust.method = p_adjust_method,
    verbose = FALSE,
    ...
  ) %>%
    dplyr::filter(.data$p_val_adj < p_val_adj_threshold)

  marker_lookup <- .unpack_markers(anno_markers, val_col = "gene", ind_col = "type")

  cli::cli_progress_step("Annotating clusters...", spinner = TRUE)

  label_assignments <- cluster_markers %>%
    dplyr::filter(.data$gene %in% marker_lookup$gene) %>%
    dplyr::inner_join(marker_lookup, by = "gene", relationship = "many-to-many") %>%
    dplyr::group_by(.data$cluster, .data$type) %>%
    dplyr::summarise(
      score = median(.data$avg_log2FC),
      n_markers = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::slice_max(.data$score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()


  label_map <- stats::setNames(
    object = as.character(label_assignments$type),
    nm = as.character(label_assignments$cluster)
    )

  cluster_anno <- dplyr::recode(
    .x = as.character(sObj@meta.data[[initial_cluster_col]]),
    !!!label_map,
    .default = unnassigned_lab
  )

  sObj@meta.data[[label_out_col]] <- cluster_anno

  return(sObj)

}
