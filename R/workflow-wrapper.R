#' Run the complete CellVoteR ensemble pipeline
#'
#' Executes the full automated annotation workflow by running four distinct classification
#' strategies in parallel and deriving a final consensus label for every cell.
#'
#' @details
#' \strong{Ensemble Strategy:}
#' To maximize robustness, this function runs \code{\link{run_method}} four times with different configurations:
#' \enumerate{
#'   \item \strong{M1 (Base):} Clustering-based broad classification using \emph{all} features.
#'   \item \strong{M2 (Refined):} Clustering-based broad classification using \emph{only marker} features.
#'   \item \strong{M3 (Direct):} Marker-based broad classification using \emph{all} features.
#'   \item \strong{M4 (Focused):} Marker-based broad classification using \emph{only marker} features.
#' }
#'
#' \strong{Consensus & Tie-Breaking:}
#' 1. The metadata from all four runs is aggregated into a single Seurat object.
#' 2. \code{\link{tie_breakers}} is run on the clustering outputs (M1 and M2) to generate fallback labels.
#' 3. \code{\link{resolve_consensus_labels}} compares the four methods and assigns a final label based on majority voting.
#'
#' @param seu_obj A Seurat object (normalized).
#' @param broad_markers A configuration list for broad cell types (from \code{\link{build_broad_marker_config}}).
#' @param fine_markers A named list of character vectors containing high-resolution markers.
#'
#' @return A Seurat object containing the final consensus annotations in \code{meta.data$cellvoter_label}
#' and the method used in \code{meta.data$cellvoter_method}. Intermediate columns from all four runs
#' are also preserved.
#'
#' @export
#'
run_cellvoteR <- function(
    seu_obj,
    broad_markers,
    fine_markers
    ) {

  base_obj <- run_method(
    seu_obj = seu_obj,
    broad_markers = broad_markers,
    fine_markers = fine_markers,
    mode = "clustering"
  )

  res_m2 <- run_method(
    seu_obj = seu_obj,
    broad_markers = broad_markers,
    fine_markers = fine_markers,
    mode = "clustering",
    filter_markers = TRUE
  )

  res_m3 <- run_method(
    seu_obj = seu_obj,
    broad_markers = broad_markers,
    fine_markers = fine_markers,
    mode = "marker"
  )

  res_m4 <- run_method(
    seu_obj = seu_obj,
    broad_markers = broad_markers,
    fine_markers = fine_markers,
    mode = "marker",
    filter_markers = TRUE,
  )

  cols_to_add <- list(
    m2 = c("clustering_ref_broad_cluster", "clustering_ref_best_label"),
    m3 = c("marker_all_best_label"),
    m4 = c("marker_ref_best_label")
  )

  base_obj <- SeuratObject::AddMetaData(
    base_obj, metadata = res_m2@meta.data[, cols_to_add$m2], col.name = cols_to_add$m2
  )
  base_obj <- SeuratObject::AddMetaData(
    base_obj, metadata = res_m3@meta.data[, cols_to_add$m3], col.name = cols_to_add$m3
  )
  base_obj <- SeuratObject::AddMetaData(
    base_obj, metadata = res_m4@meta.data[, cols_to_add$m4], col.name = cols_to_add$m4
  )

  final_obj <- tie_breakers(
    seu_obj = base_obj,
    fine_markers = fine_markers,
    cluster_cols = c(
      "global_tie_1" = "clustering_all_broad_cluster",
      "global_tie_2" = "clustering_ref_broad_cluster"
    )
  )

  final_obj <- resolve_consensus_labels(final_obj)

  cli::cat_line()
  cli::cli_alert_success("Run completed!")
  cli::cli_alert_info("All labels are stored in the metadata of the returned Seurat object.")

  return(final_obj)

}
