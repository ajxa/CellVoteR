#' Calculate Fisher's Exact Test for gene set enrichment
#'
#' Performs a one-tailed Fisher's Exact Test to determine if the overlap between
#' a query set (e.g., cluster markers) and a reference set (e.g., cell type markers)
#' is statistically significant given a background universe of genes.
#'
#' @details
#' The function constructs a contingency table:
#' \tabular{cc}{
#'   In Both Sets (a) \tab In Set 1 Only (b) \cr
#'   In Set 2 Only (c) \tab In Neither (d)
#' }
#' It tests the alternative hypothesis that the true odds ratio is greater than 1
#' (enrichment). If \code{length(set1)} is greater than or equal to the background,
#' it assumes \code{set1} represents the entire universe and returns p = 1 (or NA).
#'
#' @param set1 Character vector. Query genes (e.g., DE markers from a cluster).
#' @param set2 Character vector. Reference genes (e.g., known markers for a cell type).
#' @param background Character vector. The "universe" of all detected genes in the experiment.
#' @param return_pval Logical. If \code{TRUE}, returns only the p-value. If \code{FALSE},
#'   returns the full \code{htest} object. Default is \code{TRUE}.
#'
#' @return A numeric p-value or a \code{htest} object.
#' @importFrom stats fisher.test
#'
#' @export
fisher_score <- function(
    set1,
    set2,
    background,
    return_pval = TRUE
) {

  if(!is.character(set1)) cli::cli_abort("set1 must be a character vector of gene names.")
  if(!is.character(set2)) cli::cli_abort("set2 must be a character vector of gene names.")
  if(!is.character(background)) cli::cli_abort("background must be a character vector of gene names.")

  # case when cluster and background are the same
  if (length(set1) >= length(background)) {
    return(ifelse(return_pval, 1, NA))
  }

  a <- length(intersect(set1, set2)) # overlap
  b <- length(set1) - a # in set1 only
  c <- length(set2) - a # in set2 only
  d <- length(background) - (a+b+c) # neither
  mat <- matrix(c(a, b, c, d), nrow = 2)
  out <- fisher.test(mat, alternative = "greater")

  if(return_pval) return(out$p.value) else return(out)
}




#' Compute set similarity metrics
#'
#' Calculates the Jacquard index and overlap coefficient between two gene sets.
#'
#' @details
#' \strong{Jaccard Index:}
#' \deqn{J(A,B) = \frac{|A \cap B|}{|A \cup B|}}
#' Represents the intersection over the union.
#'
#' \strong{Overlap Coefficient:}
#' \deqn{O(A,B) = \frac{|A \cap B|}{\min(|A|, |B|)}}
#' Represents the intersection relative to the size of the smaller set. This is often
#' more useful when comparing gene sets of vastly different sizes (e.g., a small
#' specific reference list vs a large cluster DE list).
#'
#' @param cluster_genes Character vector. Genes found in the cluster.
#' @param ref_genes Character vector. Genes found in the reference signature.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{jaccard}: Numeric value (0-1).
#'   \item \code{overlap}: Numeric value (0-1).
#' }
#'
#' @export
calculate_similarity <- function(
    cluster_genes,
    ref_genes
) {
  intersection <- length(intersect(cluster_genes, ref_genes))
  union_size   <- length(union(cluster_genes, ref_genes))

  jaccard <- intersection / union_size

  overlap_coef <- intersection / min(length(cluster_genes), length(ref_genes))

  return(list(jaccard = jaccard, overlap = overlap_coef))
}




#' Score sub-clusters against reference signatures
#'
#' Evaluates the correspondence between sub-clusters and reference cell types using
#' multiple metrics: Fisher's exact test, set similarity indices, and average expression scores.
#'
#' @details
#' \strong{Scoring Logic:}
#' This function creates a grid of all combinations (Clusters x Reference Labels) and calculates:
#' \itemize{
#'   \item \strong{Enrichment:} Fisher's Exact Test p-values (via \code{\link{fisher_score}}).
#'   \item \strong{Similarity:} Jaccard and Overlap coefficients (via \code{\link{calculate_similarity}}).
#'   \item \strong{Module Score:} The average expression level of the reference markers across
#'   the entire sub-dataset (not per-cluster).
#' }
#'
#' \strong{Single Cluster Edge Case:}
#' If \code{sub_de_markers} contains only one cluster, Fisher's test cannot be reliably performed
#' (no background contrast). In this case, p-values are forced to 1 and \code{NegLog10P} to 0.
#'
#' @param sub_sObj A Seurat object containing the sub-clustered cells.
#' @param sub_de_markers A named list of DE markers for each cluster (from \code{subcluster_cells}).
#' @param sub_ref_markers A named list of reference markers to test against.
#' @param sub_fisher_background Character vector. The background gene universe (usually all genes in \code{sub_sObj}).
#' @param padj_method Character. The method used for p-value adjustment (passed to \code{p.adjust}).
#'   Default is "BH".
#'
#' @return A data frame (tibble) with one row per Cluster-Label pair, containing:
#' \itemize{
#'   \item \code{cluster}: The numeric/character ID of the subcluster.
#'   \item \code{cell_label}: The name of the reference cell type.
#'   \item \code{pval}: Raw Fisher's test p-value.
#'   \item \code{padj}: Adjusted p-value.
#'   \item \code{NegLog10P}: -log10(padj), useful for plotting heatmaps.
#'   \item \code{jaccard}: Jaccard similarity index.
#'   \item \code{overlap_coef}: Overlap coefficient.
#'   \item \code{avg_module_score}: Mean expression of the reference markers in the data.
#' }
#' @importFrom stats p.adjust
#'
#' @export
score_subclusters <- function(
    sub_sObj,
    sub_de_markers,
    sub_ref_markers,
    sub_fisher_background,
    padj_method = "BH"
) {
  # Input validation ----
  .valid_sObj_input(sub_sObj)
  if (is.null(names(sub_de_markers))) cli::cli_abort("{.arg sub_de_markers} must be a named list.")
  if (is.null(names(sub_ref_markers))) cli::cli_abort("{.arg sub_ref_markers} must be a named list.")

  if (length(sub_de_markers) == 0 || length(sub_ref_markers) == 0) return(data.frame())

  # Calculate average marker scores ----
  all_gene_means <- Matrix::rowMeans(Seurat::GetAssayData(sub_sObj, layer = "data"))

  ref_scores <- vapply(names(sub_ref_markers), function(ref_id) {
    valid_genes <- intersect(sub_ref_markers[[ref_id]], names(all_gene_means))
    if (length(valid_genes) == 0) return(0)
    return(mean(all_gene_means[valid_genes]))

  }, numeric(1))


  # Calculate fisher and similarity ----

  results <- tidyr::expand_grid(
    cluster = names(sub_de_markers),
    cell_label = names(sub_ref_markers)
  )

  is_single_cluster <- length(sub_de_markers) == 1

  stats_df <- purrr::pmap_dfr(results, function(cluster, cell_label) {

    de_genes  <- sub_de_markers[[cluster]]
    ref_genes <- sub_ref_markers[[cell_label]]

    sim <- calculate_similarity(de_genes, ref_genes)

    p_val <- 1

    if (!is_single_cluster) {
      p_val <- fisher_score(
        set1 = de_genes,
        set2 = ref_genes,
        background = sub_fisher_background
      )
    }

    list(
      pval = p_val,
      jaccard = sim$jaccard,
      overlap_coef = sim$overlap
    )
  })

  # Combine ----
  final_out_cols <- c(
    "cluster", "cell_label", "pval", "padj", "NegLog10P", "jaccard",
    "overlap_coef", "avg_module_score"
    )

  final_results <- dplyr::bind_cols(results, stats_df) %>%
    dplyr::mutate(
      avg_module_score = ref_scores[cell_label],
      padj = if (!is_single_cluster && !is.null(padj_method)) {
        p.adjust(pval, method = padj_method)
      } else {
        pval
      },
      NegLog10P = if (!is_single_cluster) -log10(.data$padj) else 0
    ) %>%
    dplyr::select(dplyr::all_of(final_out_cols))

  return(final_results)
}
