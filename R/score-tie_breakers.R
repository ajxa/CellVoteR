#' Perform global clustering for tie-breaking
#'
#' Runs a complete differential expression and marker discovery workflow on a
#' globally defined cluster column. This serves as a "ground truth" or fallback
#' annotation layer for resolving conflicting cell type assignments.
#'
#' @details
#' \strong{Workflow:}
#' \enumerate{
#'   \item \strong{Setup:} Validates the input Seurat object and the specified global cluster column.
#'   \item \strong{DE Analysis:} Identifies markers for the global clusters using \code{Seurat::FindAllMarkers}.
#'   \item \strong{Supplementation:} If a cluster has too few DE markers (less than \code{min_cluster_de_genes}),
#'   it adds highly expressed genes to the marker list to ensure downstream scoring doesn't fail.
#' }
#'
#' @param sObj A Seurat object.
#' @param marker_list A list of reference markers (will be flattened for context).
#' @param global_cluster_col Character. The metadata column defining the global clusters to analyze.
#' @param min_cluster_de_genes Integer. Minimum number of markers required per cluster before supplementation kicks in.
#' @param supplement_target_n Integer. Target number of markers to reach if supplementation is needed.
#' @param min_logfc Numeric. Minimum log-fold change for DE.
#' @param min_percent_genes Numeric. Minimum percent expression for DE.
#' @param min_adj_pval Numeric. Maximum adjusted p-value for DE.
#' @param p_adj_method Character. Method for p-value correction (default "BH").
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{clusters}: The Seurat object (potentially cleaned/subsetted).
#'   \item \code{de_results}: The full table of differential expression results.
#'   \item \code{de_markers}: A simplified list of markers per cluster.
#'   \item \code{ref_markers}: The flattened reference list.
#'   \item \code{fisher_test_background}: The background gene universe.
#' }
#' @importFrom SeuratObject Idents `Idents<-`
#' @export
#'
global_consensus_clusters <- function(
    sObj,
    marker_list,
    global_cluster_col,
    min_cluster_de_genes = 10,
    supplement_target_n = 100,
    min_logfc = 0.1,
    min_percent_genes = 0.1,
    min_adj_pval = 0.1,
    p_adj_method = "BH"
) {

  # Input validation ----
  .valid_sObj_input(sObj)

  if (!global_cluster_col %in% colnames(sObj@meta.data)) {
    cli::cli_abort("Column {.val {global_cluster_col}} not found in metadata.")
  }

  # Setup output structure ----
  out <- list(
    clusters = NULL,
    ref_markers = purrr::flatten(marker_list),
    fisher_test_background = unique(rownames(sObj)),
    de_results = tibble::tibble(),
    de_markers = list(),
    avg_expr = NULL
  )

  clean_sObj <- sObj
  Idents(clean_sObj) <-  clean_sObj@meta.data[[global_cluster_col]]

  if (!is.factor(clean_sObj@meta.data[[global_cluster_col]])) {
    clean_sObj@meta.data[[global_cluster_col]] <- as.factor(clean_sObj@meta.data[[global_cluster_col]])
  }

  out[["clusters"]] <- clean_sObj

  # Check for single cluster case (cannot run DE)
  n_clusters <- length(unique(Idents(clean_sObj)))
  if (n_clusters < 2) {
    cli::cli_alert_warning("Only {n_clusters} cluster found. Skipping marker identification.")
    cluster_id <- as.character(unique(Idents(clean_sObj)))
    out[["de_markers"]] <- stats::setNames(list(out$fisher_test_background), cluster_id)
    return(out)
  }

  # Differential expression analysis ----
  cli::cli_progress_step("Finding cluster markers...", spinner = TRUE)

  tryCatch({
    de_res <- Seurat::FindAllMarkers(
      object = clean_sObj,
      only.pos = TRUE,
      logfc.threshold = min_logfc,
      min.pct = min_percent_genes,
      test.use = "wilcox",
      assay = Seurat::DefaultAssay(clean_sObj),
      p.adjust.method = p_adj_method,
      verbose = FALSE
    )
  }, error = function(e) {
    cli::cli_warn("FindAllMarkers failed: {e$message}")
    de_res <<- tibble::tibble() # Fallback to empty
  })

  if (nrow(de_res) > 0) {
    out[["de_results"]] <- de_res %>%
      dplyr::filter(.data$p_val_adj < min_adj_pval) %>%
      dplyr::mutate(is_supplemental = FALSE)
  }

  out[["de_markers"]] <- distinct_cluster_markers(out[["de_results"]])

  all_clusters <- levels(Idents(clean_sObj))
  missing_clusters <- setdiff(all_clusters, names(out[["de_markers"]]))

  if (length(missing_clusters) > 0) {
    empty_list <- rep(list(character(0)), length(missing_clusters))
    names(empty_list) <- missing_clusters
    out[["de_markers"]] <- c(out[["de_markers"]], empty_list)
  }

  # Supplementation Logic ----
  # Identifies clusters that didn't generate enough DE markers
  cluster_gene_lens <- lengths(out[["de_markers"]])
  clusters_to_supplement <- names(cluster_gene_lens[cluster_gene_lens < min_cluster_de_genes])

  if (length(clusters_to_supplement) == 0)  return(out)

  cli::cli_progress_step("Supplementing low-marker clusters...", spinner = TRUE)

  out[["avg_expr"]] <- Seurat::AverageExpression(clean_sObj, verbose = FALSE)[[Seurat::DefaultAssay(clean_sObj)]]
  colnames(out[["avg_expr"]]) <- as.character(levels(Idents(clean_sObj)))

  supplement_list <- lapply(clusters_to_supplement, function(cl) {
    supplement_markers(
      de_results = out[["de_results"]],
      avg_expr = out[["avg_expr"]],
      cl = cl,
      target_n = supplement_target_n
    )
  })

  supplement_dfs <- dplyr::bind_rows(supplement_list)

  if (nrow(supplement_dfs) > 0) {
    out[["de_results"]] <- out[["de_results"]] %>%
      dplyr::filter(!.data$cluster %in% clusters_to_supplement) %>%
      dplyr::bind_rows(supplement_dfs)

    out[["de_markers"]] <- distinct_cluster_markers(out[["de_results"]])
  }

  return(out)
}





#' Generate and assign tie-breaker annotations
#'
#' Iterates through specified clustering columns, calculates scores against reference markers,
#' and assigns "best fit" labels. These labels act as tie-breakers when primary annotation
#' methods disagree.
#'
#' @details
#' For every column pair in \code{cluster_cols}:
#' \enumerate{
#'   \item Runs \code{\link{global_consensus_clusters}} to find markers for the clusters in that column.
#'   \item Scores these clusters against \code{fine_markers} using \code{\link{score_subclusters}}.
#'   \item Determines the best label for each cluster using \code{\link{get_best_label}}.
#'   \item Projects these labels back onto the main Seurat object in the specified output column.
#' }
#'
#' @param seu_obj A Seurat object.
#' @param fine_markers A named list of high-resolution reference markers.
#' @param cluster_cols A named character vector.
#'   \itemize{
#'     \item \strong{Values:} The existing metadata columns to use as input (clustering source).
#'     \item \strong{Names:} The new metadata columns to create (where the tie-breaker labels will be stored).
#'   }
#'   Example: \code{c("global_tie_1" = "leiden_res_0.5", "global_tie_2" = "leiden_res_1.0")}
#'
#' @return The input Seurat object with new metadata columns containing the tie-breaker labels.
#'
#' @export
tie_breakers <- function(
    seu_obj,
    fine_markers,
    cluster_cols
){

  cli::cli_h1("Annotating global tie-breakers")

  for (i in seq_along(cluster_cols)) {

    out_col <- names(cluster_cols)[i]
    in_col  <- cluster_cols[i]

    res_list <- global_consensus_clusters(
      sObj = seu_obj,
      marker_list = fine_markers,
      global_cluster_col = in_col,
    )

    res_list$scores <- score_subclusters(
      sub_sObj = res_list$clusters,
      sub_de_markers = res_list$de_markers,
      sub_ref_markers = res_list$ref_markers,
      sub_fisher_background = res_list$fisher_test_background
    )

    res_list$best_label <- get_best_label(res_list$scores)

    res_list$clusters <- assign_labels(
      sub_sObj = res_list$clusters,
      sub_best_label_map = res_list$best_label,
      sub_sObj_col = in_col,
      sub_sObj_new_col = out_col
    )

    new_metadata <- res_list$clusters@meta.data[, out_col, drop = FALSE]

    seu_obj <- SeuratObject::AddMetaData(seu_obj, metadata = new_metadata, col.name = out_col)

  }

  return(seu_obj)

}




#' Resolve final cell labels via majority vote
#'
#' Combines multiple annotation inputs (methods) and tie-breakers to determine a single
#' final consensus label for each cell.
#'
#' @details
#' \strong{Decision Hierarchy:}
#' For each cell, the function evaluates the labels in \code{method_cols}:
#' \enumerate{
#'   \item \strong{Strong Majority:} If any label receives > 2 votes (out of 4+ inputs),
#'   it is selected immediately.
#'   \item \strong{Tie-Breaker 1:} If no majority exists, the label in the first
#'   \code{tie_breaker_cols} is selected.
#'   \item \strong{Tie-Breaker 2:} If the first tie-breaker is invalid (NA/empty),
#'   the second tie-breaker is used.
#'   \item \strong{Unresolved:} If all else fails, the cell is labeled as \code{unknown_label}.
#' }
#'
#' @param seu_obj A Seurat object containing all method and tie-breaker columns.
#' @param method_cols Character vector. The metadata columns representing primary votes
#'   (e.g., results from different clustering resolutions or marker scoring approaches).
#' @param tie_breaker_cols Character vector. The metadata columns to use as fallbacks.
#'   Order matters: the first column is checked before the second.
#' @param unknown_label Character. The label to assign if no consensus is reached. Default "Unknown".
#'
#' @return The Seurat object with two new metadata columns:
#' \itemize{
#'   \item \code{cellvoter_label}: The final resolved cell type.
#'   \item \code{cellvoter_method}: The logic used to reach that decision (e.g., "Majority Vote", "Tie-Breaker 1").
#' }
#' @export
resolve_consensus_labels <- function(
    seu_obj,
    method_cols = c("clustering_all_best_label",
                    "clustering_ref_best_label",
                    "marker_all_best_label",
                    "marker_ref_best_label"),
    tie_breaker_cols = c("global_tie_1", "global_tie_2"),
    unknown_label = "Unknown"
) {
  # Validation ----
  all_req_cols <- c(method_cols, tie_breaker_cols)
  missing_cols <- setdiff(all_req_cols, colnames(seu_obj@meta.data))

  if (length(missing_cols) > 0) {
    cli::cli_abort("The following required columns are missing from metadata: {.val {missing_cols}}")
  }

  cli::cli_h1("Resolving Consensus Labels")
  cli::cli_progress_bar("Voting...", total = ncol(seu_obj))

  meta_df <- seu_obj@meta.data[, all_req_cols]

  # Consensus Logic ----
  results <- apply(meta_df, 1, function(row) {

    votes <- row[method_cols]
    votes <- votes[!is.na(votes) & votes != ""]

    t1 <- row[tie_breaker_cols[1]]
    t2 <- row[tie_breaker_cols[2]]

    if (length(votes) > 0) {
      vote_counts <- table(votes)
      max_votes <- max(vote_counts)
      best_candidates <- names(vote_counts)[vote_counts == max_votes]
    } else max_votes <- 0

    # Decision Hierarchy ---

    if (max_votes > 2) return(c(label = best_candidates[1], method = "Majority Vote"))


    if (!is.na(t1) && t1 != "") {
      method_note <- if (!is.na(t2) && t1 == t2) "Tie-Breaker (Consensus)" else "Tie-Breaker 1"
      return(c(label = t1, method = method_note))
    }

    if (!is.na(t2) && t2 != "") return(c(label = t2, method = "Tie-Breaker 2"))

    return(c(label = unknown_label, method = "Unresolved"))
  })

  # Format output and return ----
  results_df <- as.data.frame(t(results))
  colnames(results_df) <- c("cellvoter_label", "cellvoter_method")


  seu_obj <- SeuratObject::AddMetaData(
    seu_obj,
    results_df$cellvoter_label,
    col.name = "cellvoter_label"
  )

  seu_obj <- SeuratObject::AddMetaData(
    seu_obj,
    results_df$cellvoter_method,
    col.name = "cellvoter_method"
  )

  cli::cli_alert_success("Consensus resolution complete.")

  return(seu_obj)
}
