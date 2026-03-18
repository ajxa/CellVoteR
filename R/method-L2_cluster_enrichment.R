#' Rank markers for each cluster using scran::findMarkers()
#'
#' Runs \code{\link[scran]{findMarkers}} and stores the resulting ranked marker
#' tables in \code{metadata(sce)} for downstream annotation, or optionally
#' returns the ranked marker results directly as a list.
#'
#' This function is intended as a reusable preprocessing step for:
#' \itemize{
#'   \item broad label assignment from ranked broad markers, and
#'   \item fine label assignment via Fisher's Exact Test on top-ranked genes.
#' }
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param cluster_col Character scalar. Column in \code{colData(sce)} containing the cluster labels to test.
#' @param assay_type Character scalar. Assay to use. Defaults to \code{"logcounts"}.
#' @param test_type Character scalar. Differential expression test passed to
#'   \code{\link[scran]{findMarkers}}. One of \code{"wilcox"} or \code{"t"}.
#'   Defaults to \code{"wilcox"}.
#' @param direction Character scalar. Direction of testing. One of
#'   \code{"up"}, \code{"down"}, or \code{"any"}. Defaults to \code{"up"}.
#' @param pval_type Character scalar. A string specifying how p-values are to
#' be combined across pairwise comparisons for a given group/cluster. For more
#' information see \code{\link[scran]{combineMarkers}}. This can be one of
#' \code{"any"}, \code{"some"}, or \code{"all"}. Defaults to \code{"any"}.
#' @param min_prop Numeric scalar in \code{(0, 1]}. This sets the minimum
#' proportion of pairwise comparisons in which a gene must be
#' DE (in the expected direction) for it to contribute to that gene's
#' summary statistic. This is only relevant when \code{pval_type} is set to
#' \code{"any"}, and will ensure that a gene will only be high-ranked if it is
#' among the top-ranked genes in at least \code{min.prop} of the pairwise comparisons.
#' If the \code{pval_type} is set to \code{"some"} or \code{"all"}, then
#' this will be ignored. For more information see \code{\link[scran]{combineMarkers}}.
#' Defaults to 0.25, in order to be lenient and identify markers amongst
#' similar populations.
#' @param metadata_key Character scalar. Name of the metadata entry used to
#'   store the ranked marker results when \code{return_list = FALSE}. Defaults
#'   to \code{"broad_cluster_markers"}.
#' @param return_list Logical scalar. If \code{FALSE} (default), store the
#'   ranked marker tables in \code{metadata(sce)[[metadata_key]]} and return the
#'   modified \code{sce}. If \code{TRUE}, return the output list directly
#'   instead of modifying \code{sce}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object.
#'   Defaults to \code{\link[BiocParallel]{SerialParam}()}.
#'
#' @return
#' If \code{return_list = FALSE}, the input \code{SingleCellExperiment} with
#' ranked marker tables stored in \code{metadata(sce)[[metadata_key]]}.
#'
#' If \code{return_list = TRUE}, a named list with components:
#' \describe{
#'   \item{\code{marker_tables}}{The output of \code{\link[scran]{findMarkers}}.}
#'   \item{\code{params}}{A list of parameters used to generate the ranked
#'   marker tables.}
#' }
#'
#' @export
rank_cluster_markers <- function(sce,
                                 cluster_col,
                                 assay_type    =  "logcounts",
                                 test_type     =  c("wilcox", "t"),
                                 direction     =  c("up", "down", "any"),
                                 pval_type     =  c("any", "some", "all"),
                                 min_prop      =  0.25,
                                 metadata_key  =  "broad_cluster_markers",
                                 return_list   =  FALSE,
                                 BPPARAM       =  BiocParallel::SerialParam()
                                 ) {

  test_type <- match.arg(test_type)
  direction <- match.arg(direction)
  pval_type <- match.arg(pval_type)

  if (!inherits(sce, "SingleCellExperiment")) cli::cli_abort("{.arg sce} must be a {.cls SingleCellExperiment}.")

  if (!cluster_col %in% colnames(SummarizedExperiment::colData(sce))) cli::cli_abort("Column {.val {cluster_col}} was not found in {.fn colData}.")

  if (!assay_type %in% SummarizedExperiment::assayNames(sce)) {
    cli::cli_abort(c(
      "Assay {.val {assay_type}} is missing from {.arg sce}.",
      "i" = "Available assays: {.val {SummarizedExperiment::assayNames(sce)}}"
    ))
  }

  checkmate::assert_number(min_prop, lower = 0, upper = 1)
  checkmate::assert_string(metadata_key, min.chars = 1)

  cluster_ids <- as.character(SummarizedExperiment::colData(sce)[[cluster_col]])

  marker_tables <- scran::findMarkers(
    x = sce,
    groups = cluster_ids,
    assay.type = assay_type,
    test.type = test_type,
    direction = direction,
    pval.type = pval_type,
    min.prop = min_prop,
    BPPARAM = BPPARAM
  )

  out <- list(
    marker_tables = marker_tables,
    params = list(
      cluster_col = cluster_col,
      assay_type = assay_type,
      test_type = test_type,
      direction = direction,
      pval_type = pval_type,
      min_prop = min_prop
    )
  )

  if (isTRUE(return_list)) return(out)

  S4Vectors::metadata(sce)[[metadata_key]] <- out

  return(sce)
}


#' Assign broad labels based on the median rank of validated markers
#'
#' Assigns broad cluster labels using ranked marker tables produced by \code{\link{rank_cluster_markers}}.
#'
#' Ranked marker tables may be supplied either:
#' \itemize{
#'   \item directly via the \code{ranked_markers} argument, or
#'   \item indirectly via \code{metadata(sce)[[ranked_markers_key]]}.
#' }
#'
#' Broad marker definitions may be supplied either:
#' \itemize{
#'   \item directly via the \code{broad_config} argument, or
#'   \item indirectly via \code{metadata(sce)[[marker_config_key]][["broad"]]}.
#' }
#'
#' For each cluster, genes are ordered by ascending \code{FDR} then descending
#' effect size (\code{summary.AUC} for Wilcoxon, \code{summary.logFC} for
#' t-tests). Within each broad category, only markers passing both
#' \code{FDR <= fdr_threshold} and \code{effect size >= effect_threshold} are
#' retained. The category score is the median rank of passing markers; the
#' best category is chosen by lowest median rank with user-supplied
#' \code{priority} as a tie-breaker.
#'
#' @section Collapsed broad label edge case:
#' If all clusters receive the same broad label but more than one unsupervised
#' cluster exists, the original cluster labels are retained in \code{label_col}
#' rather than collapsing to a single label. This preserves cluster structure
#' for downstream subclustering.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param broad_config Named list or \code{NULL}. Validated broad marker
#'   definitions as produced by \code{\link{build_broad_marker_config}}. If
#'   \code{NULL}, extracted from \code{metadata(sce)[[marker_config_key]][["broad"]]}.
#' @param marker_config_key Character scalar. Metadata key containing the
#'   marker configuration. Only used when \code{broad_config = NULL}.
#'   Defaults to \code{"marker_config"}.
#' @param ranked_markers List or \code{NULL}. Ranked marker result from
#'   \code{\link{rank_cluster_markers}(return_list = TRUE)}. If \code{NULL},
#'   extracted from \code{metadata(sce)[[ranked_markers_key]]}.
#' @param ranked_markers_key Character scalar. Metadata key containing ranked
#'   marker tables. Only used when \code{ranked_markers = NULL}.
#'   Defaults to \code{"broad_cluster_markers"}.
#' @param cluster_col Character scalar. \code{colData} column containing
#'   cluster identifiers. Defaults to \code{"cluster_broad_hvg"}.
#' @param label_col Character scalar. Name of the output \code{colData} column.
#'   Defaults to \code{"broad_cluster"}.
#' @param fdr_threshold Numeric scalar. Maximum FDR for a marker to be
#'   considered significant. Defaults to \code{0.05}.
#' @param effect_threshold Numeric scalar. Minimum effect size (AUC for
#'   Wilcoxon, logFC for t-test). Defaults to \code{0.6}.
#' @param unassigned_label Character scalar. Label assigned to clusters with
#'   no passing markers. Defaults to \code{"other"}.
#'
#' @return The input \code{sce} with a new \code{colData} column
#'   \code{label_col} containing assigned broad labels.
#'
#' @export
label_broad_clusters <- function(sce,
                                 broad_config        = NULL,
                                 marker_config_key   = "marker_config",
                                 ranked_markers      = NULL,
                                 ranked_markers_key  = "broad_cluster_markers",
                                 cluster_col         = "cluster_broad_hvg",
                                 label_col           = "broad_cluster",
                                 fdr_threshold       = 0.05,
                                 effect_threshold    = 0.6,
                                 unassigned_label    = "other"
                                 ) {

  if (!inherits(sce, "SingleCellExperiment")) cli::cli_abort("{.arg sce} must be a {.cls SingleCellExperiment}.")

  checkmate::assert_string(cluster_col,        min.chars = 1L)
  checkmate::assert_string(label_col,          min.chars = 1L)
  checkmate::assert_string(unassigned_label,   min.chars = 1L)
  checkmate::assert_number(fdr_threshold,      lower = 0, upper = 1)
  checkmate::assert_number(effect_threshold,   lower = 0)

  if (!cluster_col %in% colnames(SummarizedExperiment::colData(sce))) cli::cli_abort("Column {.val {cluster_col}} not found in {.fn colData}.")

  broad_config <- .resolve_broad_config(sce, broad_config, marker_config_key)

  rank_res <- .resolve_ranked_markers(sce, ranked_markers, ranked_markers_key)

  if (!is.list(rank_res) || !all(c("marker_tables", "params") %in% names(rank_res))) {
    cli::cli_abort(
      "Ranked marker input must be a list with {.val marker_tables} and {.val params} entries."
    )
  }

  marker_tables <- rank_res$marker_tables
  test_type     <- rank_res$params$test_type
  eff_col       <- switch(test_type, wilcox = "summary.AUC", "summary.logFC")

  if (!eff_col %in% colnames(marker_tables[[1L]])) {
    cli::cli_abort(
      "Expected effect-size column {.val {eff_col}} not found in marker tables."
    )
  }

  new_labels <- vapply(
    X         = names(marker_tables),
    FUN       = .label_single_cluster,
    FUN.VALUE = character(1L),
    marker_tables    = marker_tables,
    broad_config     = broad_config,
    eff_col          = eff_col,
    fdr_threshold    = fdr_threshold,
    effect_threshold = effect_threshold,
    unassigned_label = unassigned_label
  )

  sce <- .assign_broad_labels(
    sce              = sce,
    new_labels       = new_labels,
    cluster_col      = cluster_col,
    label_col        = label_col,
    unassigned_label = unassigned_label
  )

  return(sce)
}


#' Resolve broad marker config from argument or SCE metadata
#' @noRd
.resolve_broad_config <- function(sce, broad_config, marker_config_key) {

  if (!is.null(broad_config)) {
    checkmate::assert_list(broad_config, names = "named", min.len = 1L)
    return(broad_config)
  }

  checkmate::assert_string(marker_config_key, min.chars = 1L)

  present_meta <- names(S4Vectors::metadata(sce))

  if (!marker_config_key %in% present_meta) {
    cli::cli_abort(
      "Metadata key {.val {marker_config_key}} not found in {.arg sce}. \\
       Available keys: {.val {present_meta}}."
    )
  }

  broad_config <- S4Vectors::metadata(sce)[[marker_config_key]][["broad"]]

  if (!is.list(broad_config) || length(broad_config) == 0L) {
    cli::cli_abort(
      "{.val {marker_config_key}$broad} must be a non-empty named list. \\
       Run {.fn build_broad_marker_config} first."
    )
  }

  return(broad_config)
}


#' Score a single cluster against all broad categories and return the best label
#' @noRd
.label_single_cluster <- function(cluster_id,
                                  marker_tables,
                                  broad_config,
                                  eff_col,
                                  fdr_threshold,
                                  effect_threshold,
                                  unassigned_label) {

  tab <- as.data.frame(marker_tables[[cluster_id]])
  tab <- tab[order(tab$FDR, -tab[[eff_col]]), , drop = FALSE]
  tab$global_rank <- seq_len(nrow(tab))

  # Score each broad category
  cat_scores <- lapply(names(broad_config), function(cat_name) {

    def   <- broad_config[[cat_name]]
    genes <- intersect(def$markers, rownames(tab))

    if (length(genes) == 0L) return(NULL)

    sub_tab <- tab[genes, , drop = FALSE]
    passing <- sub_tab$FDR <= fdr_threshold & sub_tab[[eff_col]] >= effect_threshold

    if (!any(passing)) return(NULL)

    out <- list( label = cat_name,
                 median_rank = median(sub_tab$global_rank[passing]),
                 priority = def$priority
    )
    return(out)
  })

  cat_scores <- Filter(Negate(is.null), cat_scores)

  if (length(cat_scores) == 0L) return(unassigned_label)

  score_df <- do.call(rbind, lapply(cat_scores, as.data.frame))
  score_df <- score_df[order(score_df$median_rank, score_df$priority), ]

  return(as.character(score_df$label[1L]))
}


#' Map per-cluster labels back to cells and assign to colData
#' @noRd
.assign_broad_labels <- function(sce,
                                 new_labels,
                                 cluster_col,
                                 label_col,
                                 unassigned_label) {

  cluster_ids        <- as.character(SummarizedExperiment::colData(sce)[[cluster_col]])
  n_original         <- length(unique(cluster_ids))
  unique_new_labels  <- unique(new_labels)

  # Edge case: all clusters collapsed to the same label
  if (length(unique_new_labels) == 1L && n_original > 1L) {


    print_alert("Broad annotation labelled all clusters as {.val {unique_new_labels}}", type = "w")
    print_alert("Retaining original cluster labels for subclustering.", type = "i")

    SummarizedExperiment::colData(sce)[[label_col]] <- SummarizedExperiment::colData(sce)[[cluster_col]]

    return(sce)
  }

  per_cell_labels <- unname(new_labels[cluster_ids])

  # Factor with unassigned level last
  non_other    <- sort(unique(per_cell_labels[per_cell_labels != unassigned_label]))
  label_levels <- c(non_other, unassigned_label)

  SummarizedExperiment::colData(sce)[[label_col]] <- factor(
    x = per_cell_labels,
    levels = label_levels
  )

  return(sce)
}


#' Rank cluster markers and assign broad cluster labels
#'
#' Convenience wrapper around \code{\link{rank_cluster_markers}} and
#' \code{\link{label_broad_clusters}} that computes ranked marker tables for an
#' existing clustering and then assigns broad labels using validated broad
#' marker definitions stored in \code{metadata(sce)} or supplied directly.
#'
#' This function is intended for workflows where:
#' \itemize{
#'   \item an unsupervised clustering has already been generated,
#'   \item a validated marker configuration has already been attached to the
#'   object (or is supplied directly via \code{broad_config}), and
#'   \item broad labels are assigned by testing whether curated broad markers
#'   are significantly and strongly ranked within each cluster.
#' }
#'
#' Marker ranking is performed with \code{\link[scran]{findMarkers}} via
#' \code{\link{rank_cluster_markers}}, and broad labels are then assigned by
#' \code{\link{label_broad_clusters}} using the median rank of markers that pass
#' the chosen FDR and effect-size thresholds.
#'
#' @section Edge case: collapsed broad labels:
#' In some datasets, all original unsupervised clusters may be assigned the same
#' broad label (for example, all clusters are labelled \code{"other"} or all
#' are labelled \code{"immune"}). In that situation, downstream subclustering
#' should still be able to operate on the original unsupervised clusters rather
#' than on a single collapsed broad label.
#'
#' To support this, \code{\link{label_broad_clusters}} checks whether broad
#' annotation has collapsed all clusters to one unique label while the original
#' clustering in \code{cluster_col} still contains more than one cluster. If
#' this occurs, the original cluster labels are retained in \code{label_col}
#' instead of the collapsed broad annotation. This preserves the original
#' cluster structure for downstream subclustering.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param broad_config Named list or \code{NULL}. Validated broad marker
#'   definitions as produced by \code{\link{build_broad_marker_config}}. If
#'   \code{NULL}, extracted from \code{metadata(sce)[[marker_config_key]][["broad"]]}.
#' @param cluster_col Character scalar. Column in \code{colData(sce)}
#'   containing the cluster labels to annotate. Defaults to
#'   \code{"cluster_broad_hvg"}.
#' @param marker_config_key Character scalar. Metadata entry containing the
#'   validated marker configuration list with a \code{$broad} element. Only
#'   used when \code{broad_config = NULL}. Defaults to \code{"marker_config"}.
#' @param ranked_markers_key Character scalar. Metadata entry used to store the
#'   ranked marker tables generated by \code{\link{rank_cluster_markers}}.
#'   Defaults to \code{"broad_cluster_markers"}.
#' @param ranked_markers List or \code{NULL}. Ranked marker result from
#'   \code{\link{rank_cluster_markers}(return_list = TRUE)}. If \code{NULL},
#'   marker ranking is performed internally and stored under
#'   \code{ranked_markers_key}.
#' @param label_col Character scalar. Name of the output broad label column in
#'   \code{colData(sce)}. Defaults to \code{"broad_cluster"}.
#' @param assay_type Character scalar. Assay to use for
#'   \code{\link{rank_cluster_markers}}. Defaults to \code{"logcounts"}.
#' @param test_type Character scalar. Differential expression test passed to
#'   \code{\link{rank_cluster_markers}}. One of \code{"wilcox"} or \code{"t"}.
#'   Defaults to \code{"wilcox"}.
#' @param direction Character scalar. Direction of testing passed to
#'   \code{\link{rank_cluster_markers}}. One of \code{"up"}, \code{"down"},
#'   or \code{"any"}. Defaults to \code{"up"}.
#' @param pval_type Character scalar. P-value combination mode passed to
#'   \code{\link{rank_cluster_markers}}. One of \code{"any"}, \code{"some"},
#'   or \code{"all"}. Defaults to \code{"any"}.
#' @param min_prop Numeric scalar in \code{(0, 1]}. Passed to
#'   \code{\link{rank_cluster_markers}}. Defaults to \code{0.25}.
#' @param fdr_threshold Numeric scalar. Maximum FDR to consider a gene
#'   significant in \code{\link{label_broad_clusters}}. Defaults to
#'   \code{0.05}.
#' @param effect_threshold Numeric scalar. Minimum effect size required for a
#'   marker to contribute to broad label assignment. Interpreted as AUC for
#'   Wilcoxon tests and log-fold change for t-tests. Defaults to \code{0.6}.
#' @param unassigned_label Character scalar. Label used when no broad category
#'   passes the assignment criteria. Defaults to \code{"other"}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object.
#'   Defaults to \code{\link[BiocParallel]{SerialParam}()}.
#'
#' @return The input \code{SingleCellExperiment} with:
#' \describe{
#'   \item{\code{metadata(sce)[[ranked_markers_key]]}}{Ranked marker tables
#'     produced by \code{\link{rank_cluster_markers}}.}
#'   \item{\code{colData(sce)[[label_col]]}}{Assigned broad labels. In the
#'     special case where all clusters collapse to a single broad label while
#'     \code{cluster_col} still contains multiple original clusters, this column
#'     will instead contain the original cluster labels.}
#' }
#'
#' @examples
#' \dontrun{
#' sce <- annotate_broad_clusters(
#'   sce        = sce,
#'   cluster_col = "cluster_broad_hvg",
#'   label_col   = "broad_cluster"
#' )
#'
#' # Supplying broad_config directly
#' sce <- annotate_broad_clusters(
#'   sce         = sce,
#'   broad_config = my_broad_config,
#'   cluster_col  = "cluster_broad_hvg",
#'   label_col    = "broad_cluster"
#' )
#' }
#'
#' @export
annotate_broad_clusters <- function(sce,
                                    broad_config        =  NULL,
                                    ranked_markers      =  NULL,
                                    cluster_col         =  "cluster_broad_hvg",
                                    marker_config_key   =  "marker_config",
                                    ranked_markers_key  =  "broad_cluster_markers",
                                    label_col           =  "broad_cluster",
                                    assay_type          =  "logcounts",
                                    test_type           =  c("wilcox", "t"),
                                    direction           =  c("up", "down", "any"),
                                    pval_type           =  c("any", "some", "all"),
                                    min_prop            =  0.25,
                                    fdr_threshold       =  0.05,
                                    effect_threshold    =  0.6,
                                    unassigned_label    =  "other",
                                    BPPARAM             =  BiocParallel::SerialParam()
                                    ) {

  test_type <- match.arg(test_type)
  direction <- match.arg(direction)
  pval_type <- match.arg(pval_type)

  if (is.null(ranked_markers)) {
    sce <- rank_cluster_markers(
      sce          = sce,
      cluster_col  = cluster_col,
      assay_type   = assay_type,
      test_type    = test_type,
      direction    = direction,
      pval_type    = pval_type,
      min_prop     = min_prop,
      metadata_key = ranked_markers_key,
      BPPARAM      = BPPARAM
    )
  }

  sce <- label_broad_clusters(
    sce                = sce,
    broad_config       = broad_config,
    marker_config_key  = marker_config_key,
    ranked_markers     = ranked_markers,
    ranked_markers_key = ranked_markers_key,
    cluster_col        = cluster_col,
    label_col          = label_col,
    fdr_threshold      = fdr_threshold,
    effect_threshold   = effect_threshold,
    unassigned_label   = unassigned_label
  )

  return(sce)
}
