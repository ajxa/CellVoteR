#' Supplement insufficient DE Markers
#'
#' Adds highly expressed genes to a cluster's marker list if the Differential Expression (DE)
#' analysis yields too few results. This ensures downstream enrichment tests do not fail due
#' to empty input sets.
#'
#' @details
#' The function identifies the top expressed genes in the specific cluster (based on \code{avg_expr})
#' that are not already present in \code{de_results}. It creates a "fallback" data frame where
#' statistical columns (p_val, p_val_adj) are set to 1, and \code{is_supplemental} is flagged as \code{TRUE}.
#'
#' @param de_results Data frame. The existing DE results (output of \code{FindAllMarkers}).
#' @param avg_expr Matrix. Average expression matrix (genes x clusters).
#' @param cl Character or Numeric. The cluster ID to supplement.
#' @param target_n Integer. The desired total number of markers (existing + supplemental).
#'
#' @return A data frame of supplemental markers formatted to match \code{de_results}.
#' @export
supplement_markers <- function(
    de_results,
    avg_expr,
    cl,
    target_n = 100
) {

  cli::cli_alert_warning("Supplementing cluster with insufficient DE markers: {.val {cl}}")

  existing <- de_results %>% dplyr::filter(cluster == cl)
  n_existing <- nrow(existing)

  top_expr_genes <- names(sort(avg_expr[, cl], decreasing = TRUE))
  needed_n <- target_n - n_existing

  new_genes <- setdiff(top_expr_genes, existing$gene)[1:needed_n]

  fallback_df <- data.frame(
    p_val = 1, p_val_adj = 1, avg_log2FC = 0.01,
    pct.1 = 1, pct.2 = 1,
    cluster = cl,
    gene = new_genes,
    is_supplemental = TRUE
  )

  return(fallback_df)

}




#' Perform sub-clustering and marker identification
#'
#' Re-clusters a subset of cells, identifies marker genes for the new clusters, and
#' ensures every cluster has a minimum number of features for downstream annotation.
#'
#' @details
#' \strong{Workflow:}
#' \enumerate{
#'   \item \strong{Re-clustering:} Calls \code{\link{cluster_cells}} on the provided \code{sub_sObj}.
#'   \item \strong{Marker Discovery:} Runs \code{Seurat::FindAllMarkers} on the new clusters.
#'   \item \strong{Quality Check:} If any cluster has fewer than \code{min_cluster_de_genes} markers,
#'   it calls \code{\link{supplement_markers}} to add top-expressed genes to the list.
#' }
#'
#' @param sub_sObj A Seurat object containing the subset of cells to re-cluster.
#' @param sub_marker_list A data frame or named list of reference markers.
#'   If a data frame is provided, it is converted to a list using \code{\link{distinct_cluster_markers}}.
#' @param subcluster_name Character. The prefix for the new cluster identities. Default is "subcluster".
#' @param min_cluster_de_genes Integer. If a cluster has fewer DE markers than this,
#'   supplemental markers are added. Default is 10.
#' @param supplement_target_n Integer. The total number of markers to aim for when supplementing.
#' @param min_logfc Numeric. Minimum log-fold change for DE analysis.
#' @param min_percent_genes Numeric. Minimum percentage expression for DE analysis.
#' @param min_adj_pval Numeric. Maximum adjusted p-value for DE markers.
#' @param p_adj_method Character. Correction method for p-values (default "BH").
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{clusters}: The re-clustered Seurat object.
#'   \item \code{de_results}: Data frame of all differential expression results (including supplemental).
#'   \item \code{de_markers}: A named list of markers per cluster (simplified format).
#'   \item \code{ref_markers}: The cleaned reference marker list used for context.
#'   \item \code{avg_expr}: The average expression matrix for the new clusters.
#'   \item \code{fisher_test_background}: Character vector of all genes in the object (background for enrichment tests).
#' }
#' @export
subcluster_cells <- function(
    sub_sObj,
    sub_marker_list,
    subcluster_name = NULL,
    min_cluster_de_genes = 10,
    supplement_target_n = 100,
    min_logfc = 0.1,
    min_percent_genes = 0.1,
    min_adj_pval = 0.1,
    p_adj_method = "BH"
){
  out <- list(
    clusters = NULL,
    ref_markers = NULL,
    fisher_test_background = NULL,
    de_results = NULL,
    de_markers = NULL,
    avg_expr = NULL
  )

  if (!exists("sub_marker_list")) cli::cli_abort("{.var sub_marker_list} is missing.")
  if (is.data.frame(sub_marker_list)){

    cli::cli_alert_info("Converting {.var sub_marker_list} from data.frame to list.")

    out[["ref_markers"]] <- distinct_cluster_markers(
      sub_marker_list, gene_col = "marker", cluster_col = "cell_type"
    )
  } else out[["ref_markers"]] <- sub_marker_list

  if(is.null(subcluster_name)) subcluster_name = "subcluster"

  out[["clusters"]] <- cluster_cells(
    sObj = sub_sObj,
    cluster_name = subcluster_name,
    set_new_ident = TRUE,
    subclustering = TRUE,
    verbose = FALSE
  )

  out[["fisher_test_background"]] <- unique(rownames(out$clusters))

  if(length(unique(Idents(out$clusters))) > 1){

    cli::cli_progress_step("Finding subcluster markers...", spinner = TRUE)

    out[["de_results"]] <- Seurat::FindAllMarkers(
      object = out$clusters,
      only.pos = TRUE,
      logfc.threshold = min_logfc,
      min.pct = min_percent_genes,
      test.use = "wilcox",
      assay = SeuratObject::DefaultAssay(out$clusters),
      p.adjust.method = p_adj_method,
      verbose = FALSE
    )

    out$de_results <- out$de_results %>%
      dplyr::filter(.data$p_val_adj < min_adj_pval) %>%
      dplyr::mutate(is_supplemental = FALSE)

    out[["de_markers"]] <- distinct_cluster_markers(out$de_results)
    cluster_gene_lens <- lengths(out$de_markers)

    clusters_to_supplement <- names(cluster_gene_lens[cluster_gene_lens < min_cluster_de_genes])

    cli::cli_progress_step("Calculating average expression per subcluster...", spinner = TRUE)

    suppressMessages({

      out[["avg_expr"]] <- Seurat::AverageExpression(out$clusters)[[SeuratObject::DefaultAssay(out$clusters)]]

      colnames(out$avg_expr) <- as.character(levels(Idents(out$clusters)))

    })

    if(length(clusters_to_supplement) > 0){

      supplement_dfs <- vector(mode = "list", length(clusters_to_supplement))

      # supplment clusters with < min_cluster_de_genes markers
      # with top (supplement_target_n) avg expr genes
      for (cluster in seq_along(clusters_to_supplement)) {
        cl <- clusters_to_supplement[cluster]
        supplement_dfs[[cluster]] <- supplement_markers(
          de_results = out$de_results,
          avg_expr = out$avg_expr,
          cl = cl,
          target_n = supplement_target_n
        )
      }

      supplement_dfs <- dplyr::bind_rows(supplement_dfs)

      out$de_results <- out$de_results %>%
        dplyr::filter(!(.data$cluster %in% clusters_to_supplement)) %>%
        dplyr::bind_rows(supplement_dfs)

      out$de_markers <- distinct_cluster_markers(dplyr::bind_rows(out$de_results))
    }

  }else{
    cli::cli_alert_warning("Only one cluster found - skipping marker identification.")
    out[["de_markers"]] <- list("1" = out$fisher_test_background)
  }

  return(out)
}
