#' Calculate Fisher's Exact Test Scores
#'
#' @param cluster_markers Dataframe with columns 'gene' and 'cluster' (from FindAllMarkers).
#' @param reference_markers Dataframe with 'Gene' and 'Broad_Cell_Type' columns.
#' @param background_genes Character vector of all genes in the dataset.
#'
#' @export
#' Calculate Fisher's Exact Test Scores
#'
#' @param cluster_markers Dataframe with columns 'gene' and 'cluster' (from FindAllMarkers).
#' @param reference_markers Dataframe with 'gene' and 'cell_type' columns.
#' @param background_genes Character vector of all genes in the dataset.
#'
#' @export
calculate_fisher_scores <- function(cluster_markers, reference_markers, background_genes) {
  results_list <- list()
  clusters <- unique(cluster_markers$cluster)

  # --- FIX: Handle Column Name Mismatch ---
  # Check which column names exist (standardized 'cell_type' vs legacy 'Broad_Cell_Type')
  ref_cols <- colnames(reference_markers)
  type_col <- if ("cell_type" %in% ref_cols) "cell_type" else "Broad_Cell_Type"
  gene_col <- if ("gene" %in% ref_cols) "gene" else "Gene"

  # Validation
  if (!type_col %in% ref_cols || !gene_col %in% ref_cols) {
    warning("Reference markers missing required columns ('gene', 'cell_type'). Returning NULL.")
    return(NULL)
  }

  ref_types <- unique(reference_markers[[type_col]])
  n_background <- length(background_genes)

  for (clust in clusters) {
    c_genes <- unique(cluster_markers$gene[cluster_markers$cluster == clust])
    for (rtype in ref_types) {
      # Use dynamic column names
      r_genes <- unique(reference_markers[[gene_col]][reference_markers[[type_col]] == rtype])

      a <- length(intersect(c_genes, r_genes))
      b <- length(c_genes) - a
      c <- length(r_genes) - a
      d <- n_background - (a + b + c)

      mat <- matrix(c(a, b, c, d), nrow = 2)
      pval <- stats::fisher.test(mat, alternative = "greater")$p.value

      results_list[[length(results_list) + 1]] <- data.frame(
        Cluster = clust, Reference_Type = rtype, n_overlap = a, Pval = pval,
        NegLog10P = -log10(pval + 1e-300), stringsAsFactors = FALSE
      )
    }
  }

  # Safety: Return NULL if list is empty (avoids do.call error)
  if (length(results_list) == 0) return(NULL)

  return(do.call(rbind, results_list))
}

#' Get best assignments from Fisher scores
#'
#' @param scores_df Dataframe output from calculate_fisher_scores.
#'
#' @export
get_best_match <- function(scores_df) {
  # --- FIX: Safety Check for NULL/Empty Input ---
  if (is.null(scores_df) || nrow(scores_df) == 0) {
    return(NULL)
  }

  best_matches <- do.call(rbind, lapply(split(scores_df, scores_df$Cluster), function(df) {
    df[which.min(df$Pval), ]
  }))
  labels <- best_matches$Reference_Type
  names(labels) <- best_matches$Cluster
  return(labels)
}
