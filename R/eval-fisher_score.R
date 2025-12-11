#' Calculate Fisher's Exact Test Scores
#'
#' @param cluster_markers Dataframe with columns 'gene' and 'cluster' (from FindAllMarkers).
#' @param reference_markers Dataframe with 'Gene' and 'Broad_Cell_Type' columns.
#' @param background_genes Character vector of all genes in the dataset.
#'
#' @export
calculate_fisher_scores <- function(cluster_markers, reference_markers, background_genes) {
  results_list <- list()
  clusters <- unique(cluster_markers$cluster)
  ref_types <- unique(reference_markers$Broad_Cell_Type)
  n_background <- length(background_genes)

  for (clust in clusters) {
    c_genes <- unique(cluster_markers$gene[cluster_markers$cluster == clust])
    for (rtype in ref_types) {
      r_genes <- unique(reference_markers$Gene[reference_markers$Broad_Cell_Type == rtype])

      a <- length(intersect(c_genes, r_genes))
      b <- length(c_genes) - a
      c <- length(r_genes) - a
      d <- n_background - (a + b + c)

      mat <- matrix(c(a, b, c, d), nrow = 2)
      # FIX: Explicit stats::
      pval <- stats::fisher.test(mat, alternative = "greater")$p.value

      results_list[[length(results_list) + 1]] <- data.frame(
        Cluster = clust, Reference_Type = rtype, n_overlap = a, Pval = pval,
        NegLog10P = -log10(pval + 1e-300), stringsAsFactors = FALSE
      )
    }
  }
  return(do.call(rbind, results_list))
}

#' Get best assignments from Fisher scores
#'
#' @param scores_df Dataframe output from calculate_fisher_scores.
#'
#' @export
get_best_match <- function(scores_df) {
  best_matches <- do.call(rbind, lapply(split(scores_df, scores_df$Cluster), function(df) {
    df[which.min(df$Pval), ]
  }))
  labels <- best_matches$Reference_Type
  names(labels) <- best_matches$Cluster
  return(labels)
}
