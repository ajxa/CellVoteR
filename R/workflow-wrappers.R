# Valid arguments for each function - used to catch typos early
.valid_broad_args <- c("assay_type", "test_type", "direction",
                       "pval_type", "min_prop", "fdr_threshold",
                       "effect_threshold", "unassigned_label", "BPPARAM")

.valid_rank_args    <- c("assay_type", "test_type", "direction",
                         "pval_type", "min_prop", "BPPARAM")

.valid_extract_args <- c("fdr_threshold", "effect_threshold", "target_n")

.valid_score_args   <- character(0)   # no additional args currently

.valid_assign_args  <- character(0)   # no additional args currently


#' Validate a user-supplied args list against known valid parameter names
#'
#' @param args Named list supplied by the user.
#' @param valid_names Character vector of permitted names.
#' @param arg_label Character scalar. Name shown in error messages,
#'   e.g. \code{"rank_args"}.
#'
#' @noRd
.validate_fn_args <- function(args, valid_names, arg_label) {

  if (!is.list(args)) cli::cli_abort("{.arg {arg_label}} must be a named list, not a {.cls {class(args)}}.")

  if (length(args) == 0L) return(invisible(NULL))

  if (length(valid_names) == 0L) {
    cli::cli_abort(
      "{.arg {arg_label}} does not accept any additional arguments, \\
       but received: {.val {names(args)}}."
    )
  }

  unknown <- setdiff(names(args), valid_names)

  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "Unknown argument{?s} in {.arg {arg_label}}: {.val {unknown}}.",
        "i" = "Valid arguments are: {.val {valid_names}}."
      )
    )
  }

  invisible(NULL)
}


#' Shared rank -> extract -> score -> assign fine annotation pipeline
#'
#' @param sce A SingleCellExperiment with \code{cluster_col} in colData.
#' @param cluster_col Character scalar. colData column of cluster identifiers.
#' @param marker_panel Named nested list of fine marker gene sets.
#' @param background_genes Character vector of all genes in the dataset.
#' @param rank_args Named list of additional arguments passed to
#'   \code{\link{rank_cluster_markers}}. Valid entries:
#'   \code{assay_type}, \code{test_type}, \code{direction}, \code{pval_type},
#'   \code{min_prop}, \code{BPPARAM}.
#' @param extract_args Named list of additional arguments passed to
#'   \code{\link{extract_top_markers}}. Valid entries:
#'   \code{fdr_threshold}, \code{effect_threshold}, \code{target_n}.
#' @param score_args Named list of additional arguments passed to
#'   \code{\link{score_markers_against_panel}}. No additional arguments
#'   currently accepted beyond those supplied directly.
#' @param assign_args Named list of additional arguments passed to
#'   \code{\link{assign_fine_labels}}. No additional arguments currently
#'   accepted beyond those supplied directly.
#'
#' @return List with \code{$labels} (factor) and \code{$scores} (data.frame)
#'   as returned by \code{\link{assign_fine_labels}}.
#'
#' @noRd
.run_fine_annotation <- function(sce,
                                 cluster_col,
                                 marker_panel,
                                 background_genes,
                                 rank_args    = list(),
                                 extract_args = list(),
                                 score_args   = list(),
                                 assign_args  = list()
                                 ) {

  .validate_fn_args(rank_args,    .valid_rank_args,    "rank_args")
  .validate_fn_args(extract_args, .valid_extract_args, "extract_args")
  .validate_fn_args(score_args,   .valid_score_args,   "score_args")
  .validate_fn_args(assign_args,  .valid_assign_args,  "assign_args")

  marker_ranks <- do.call(
    rank_cluster_markers,
    c(list(sce         = sce,
           cluster_col = cluster_col,
           return_list = TRUE),
      rank_args)
  )

  top_markers <- do.call(
    extract_top_markers,
    c(list(ranked_markers = marker_ranks),
      extract_args)
  )

  label_scores <- do.call(
    score_markers_against_panel,
    c(list(top_markers      = top_markers,
           marker_panel     = marker_panel,
           background_genes = background_genes),
      score_args)
  )

  do.call(
    assign_fine_labels,
    c(list(cluster_col = SummarizedExperiment::colData(sce)[[cluster_col]],
           scores      = label_scores),
      assign_args)
  )
}


#' Run the CellVoteR ensemble annotation pipeline
#'
#' Orchestrates the four annotation methods and two global tie-breakers that
#' together form the CellVoteR ensemble. Each method runs broad annotation,
#' subclustering, marker ranking, panel scoring, and fine label assignment.
#' The tie-breakers skip the broad annotation step, operating directly on the
#' pre-existing unsupervised clusters from \code{\link{prepare_sce}}.
#'
#' The returned label list is designed to be passed directly to
#' \code{\link{resolve_consensus_labels}}, which the user calls independently
#' so that voting parameters can be adjusted and re-run without repeating the
#' annotation pipeline.
#'
#' @section Pipeline structure:
#' \preformatted{
#' prepare_sce()  <- must be run before this function
#'      |
#'      |- Method 1: annotate_broad_clusters() on full SCE
#'      |- Method 2: annotate_broad_clusters() on reduced altExp
#'      |- Method 3: annotate_broad_cells()    on full SCE
#'      |- Method 4: annotate_broad_cells()    on reduced altExp
#'      |       Each: broad label -> subcluster -> rank markers
#'      |             -> score panel -> assign fine labels
#'      |
#'      |- Tie-breaker 1: HVG clusters on full SCE  (no broad step)
#'      |_ Tie-breaker 2: panel clusters on reduced (no broad step)
#' }
#'
#' @section Ranking arguments — broad vs fine annotation:
#' Marker ranking via \code{\link{rank_cluster_markers}} is used at two
#' distinct points in the pipeline, and \code{annotation_args} exposes
#' independent control over each:
#'
#' \describe{
#'   \item{\code{broad_args}}{Controls ranking inside
#'     \code{\link{annotate_broad_clusters}} (methods 1 and 2 only). These
#'     ranks are used to assign broad cell lineage labels (e.g. immune,
#'     vasculature) based on the median rank of curated broad markers. Because
#'     broad marker sets are small and highly specific, a more lenient
#'     \code{min_prop} is often appropriate here.}
#'   \item{\code{rank_args}}{Controls ranking inside
#'     \code{\link{.run_fine_annotation}} for all six methods. These ranks are
#'     used to extract top marker genes per subcluster, which are then scored
#'     against the fine marker panel via Fisher's exact test. Methods 3 and 4
#'     use \code{rank_args} only, as \code{\link{annotate_broad_cells}} does
#'     not call \code{\link{rank_cluster_markers}}.}
#' }
#'
#' If you want both steps to behave identically, supply the same values to
#' both sublists. If you want to differentiate them — for example using a
#' lenient \code{min_prop} for broad assignment but a strict one for fine
#' annotation — supply different values:
#'
#' \preformatted{
#' results <- run_cellvoter(
#'   sce,
#'   annotation_args = list(
#'     broad_args   = list(test_type = "t", min_prop = 0.1),
#'     rank_args    = list(test_type = "t", min_prop = 0.25),
#'     extract_args = list(fdr_threshold = 0.01)
#'   )
#' )
#' }
#'
#' @section Usage:
#' \preformatted{
#' # Default run
#' results <- run_cellvoter(sce)
#'
#' # With custom annotation parameters
#' results <- run_cellvoter(
#'   sce,
#'   annotation_args = list(
#'     rank_args    = list(test_type = "t", min_prop = 0.1),
#'     extract_args = list(fdr_threshold = 0.01, target_n = 50L)
#'   )
#' )
#'
#' # Resolve consensus independently
#' consensus <- resolve_consensus_labels(
#'   label_list        = results$labels,
#'   method_names      = results$method_names,
#'   tie_breaker_names = results$tie_breaker_names,
#'   allow_even_split  = TRUE,
#'   unassigned_label  = "unknown"
#' )
#'
#' sce$cellVoteR_label  <- consensus$label
#' sce$cellVoteR_method <- consensus$method
#' }
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   processed by \code{\link{prepare_sce}}. Must have a \code{logcounts}
#'   assay, \code{metadata$marker_config},
#'   \code{metadata$filtered_fine_markers}, and a \code{"user_panel"}
#'   \code{altExp}.
#' @param return_full_output Logical scalar. When \code{FALSE} (default),
#'   only the per-cell label factors are returned under \code{$labels}. When
#'   \code{TRUE}, the full output from each \code{\link{.run_fine_annotation}}
#'   call - including per-cluster score tables - is also returned under
#'   \code{$full_output}.
#' @param annotation_args Named list of argument sublists passed through to
#'   the internal annotation pipeline. Valid sublists are:
#'   \describe{
#'     \item{\code{rank_args}}{Arguments for \code{\link{rank_cluster_markers}}:
#'       \code{assay_type}, \code{test_type}, \code{direction}, \code{pval_type},
#'       \code{min_prop}, \code{BPPARAM}.}
#'     \item{\code{extract_args}}{Arguments for \code{\link{extract_top_markers}}:
#'       \code{fdr_threshold}, \code{effect_threshold}, \code{target_n}.}
#'     \item{\code{score_args}}{Arguments for
#'       \code{\link{score_markers_against_panel}}. No additional arguments
#'       currently accepted.}
#'     \item{\code{assign_args}}{Arguments for \code{\link{assign_fine_labels}}.
#'       No additional arguments currently accepted.}
#'   }
#'   Only the sublists that differ from defaults need to be supplied.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{sce}}{The input SCE with broad and subcluster label columns
#'       added to \code{colData} for all four methods, and intermediate results
#'       from the reduced altExp methods stored in \code{metadata}:
#'       \describe{
#'         \item{\code{metadata(sce)$broad_cluster_markers_reduced}}{Ranked
#'           marker tables from the reduced cluster-based method (method 2).}
#'         \item{\code{metadata(sce)$broad_cell_enrichment_reduced}}{Cell
#'           enrichment scores and parameters from the reduced enrichment-based
#'           method (method 4).}
#'       }}
#'     \item{\code{labels}}{Named list of six per-cell label factors ready to
#'       pass to \code{\link{resolve_consensus_labels}} as \code{label_list}.
#'       Names: \code{method_1}, \code{method_2}, \code{method_3},
#'       \code{method_4}, \code{global_1}, \code{global_2}.}
#'     \item{\code{method_names}}{Character vector of the four primary method
#'       names.}
#'     \item{\code{tie_breaker_names}}{Character vector of the two tie-breaker
#'       names.}
#'     \item{\code{full_output}}{Only present when \code{return_full_output =
#'       TRUE}. Named list mirroring \code{$labels} but containing the complete
#'       \code{\link{assign_fine_labels}} output (labels + per-cluster score
#'       table) for each method.}
#'   }
#'
#' @seealso \code{\link{prepare_sce}}, \code{\link{resolve_consensus_labels}}
#'
#' @export
run_cellvoter <- function(sce,
                          return_full_output = FALSE,
                          annotation_args    = list()) {
  # --- 1.) validation ----
  if (!inherits(sce, "SingleCellExperiment")) cli::cli_abort("{.arg sce} must be a {.cls SingleCellExperiment}.")

  checkmate::assert_flag(return_full_output)
  checkmate::assert_list(annotation_args, names = "named")

  unknown_sublists <- setdiff(names(annotation_args),
                              c("broad_args", "rank_args", "extract_args", "score_args", "assign_args")
                              )

  if (length(unknown_sublists) > 0L) {
    cli::cli_abort(
      c(
        "Unknown sublist{?s} in {.arg annotation_args}: {.val {unknown_sublists}}.",
        "i" = "Valid sublists: {.val {c('broad_args', 'rank_args', 'extract_args', 'score_args', 'assign_args')}}."
      )
    )
  }

  required_meta <- c("marker_config", "filtered_fine_markers")
  missing_meta  <- setdiff(required_meta, names(S4Vectors::metadata(sce)))

  if (length(missing_meta) > 0L) {
    cli::cli_abort(c(
      "Missing metadata slot{?s}: {.val {missing_meta}}.",
      "i" = "Run {.fn prepare_sce} before {.fn run_cellvoter}."
    ))
  }

  if (!"user_panel" %in% SingleCellExperiment::altExpNames(sce)) {
    cli::cli_abort(c(
      "{.val user_panel} altExp not found.",
      "i" = "Run {.fn prepare_sce} before {.fn run_cellvoter}."
    ))
  }

  # --- 2.) Merge annotation_args with defaults ----
  args <- modifyList(
    list(
      broad_args   = list(),  # controls annotate_broad_clusters ranking
      rank_args    = list(),  # controls fine annotation ranking
      extract_args = list(),
      score_args   = list(),
      assign_args  = list()
    ),
    annotation_args
  )

  .validate_fn_args(args$broad_args,   .valid_broad_args,   "broad_args")
  .validate_fn_args(args$rank_args,    .valid_rank_args,    "rank_args")
  .validate_fn_args(args$extract_args, .valid_extract_args, "extract_args")
  .validate_fn_args(args$score_args,   .valid_score_args,   "score_args")
  .validate_fn_args(args$assign_args,  .valid_assign_args,  "assign_args")

  # --- 3.) Resolve shared inputs ----
  reduced      <- SingleCellExperiment::altExp(sce, "user_panel")
  fine_markers <- S4Vectors::metadata(sce)$filtered_fine_markers

  sce$cluster_broad_reduced <- reduced$cluster

  # Inner closure: captures args, avoids repeating all four sublist args
  .annotate <- function(sce, cluster_col, background_genes) {
    .run_fine_annotation(
      sce              = sce,
      cluster_col      = cluster_col,
      marker_panel     = fine_markers,
      background_genes = background_genes,
      rank_args        = args$rank_args,
      extract_args     = args$extract_args,
      score_args       = args$score_args,
      assign_args      = args$assign_args
    )
  }

  # Method 1: cluster-based, full gene set ----
  print_h2("Method 1: cluster-based annotation", info = "full gene set")

  sce <- do.call(
    annotate_broad_clusters,
    c(list(sce       = sce,
           label_col = "broad_cluster_m1"),
      args$broad_args)
  )

  # sce <- annotate_broad_clusters(sce, label_col = "broad_cluster_m1")

  sce <- subcluster_labels(sce,
                           group_col = "broad_cluster_m1",
                           out_col   = "broad_cluster_sub_m1")

  method_1 <- .annotate(sce, "broad_cluster_sub_m1", rownames(sce))

  # Method 2: cluster-based, reduced gene set ----

  cli::cat_line()
  print_h2("Method 2: cluster-based annotation", info = "reduced gene set")

  reduced <- do.call(
    annotate_broad_clusters,
    c(list(sce          = reduced,
           broad_config = S4Vectors::metadata(reduced)$marker_config$broad,
           cluster_col  = "cluster",
           label_col    = "broad_cluster_m2"),
      args$broad_args)
  )

  reduced <- subcluster_labels(reduced,
                               group_col    = "broad_cluster_m2",
                               feature_mode = "all",
                               out_col      = "broad_cluster_sub_m2")

  sce$broad_cluster_m2     <- reduced$broad_cluster_m2
  sce$broad_cluster_sub_m2 <- reduced$broad_cluster_sub_m2

  S4Vectors::metadata(sce)$broad_cluster_markers_reduced <- S4Vectors::metadata(reduced)$broad_cluster_markers

  method_2 <- .annotate(reduced, "broad_cluster_sub_m2", rownames(reduced))

  # Method 3: enrichment-based, full gene set ----

  cli::cat_line()
  print_h2("Method 3: enrichment-based annotation", info = "full gene set")

  sce <- annotate_broad_cells(sce, label_col = "broad_enrichment_m3")

  sce <- subcluster_labels(sce,
                           group_col = "broad_enrichment_m3",
                           out_col   = "broad_enrichment_sub_m3")

  method_3 <- .annotate(sce, "broad_enrichment_sub_m3", rownames(sce))

  # Method 4: enrichment-based, reduced gene set ----

  cli::cat_line()
  print_h2("Method 4: enrichment-based annotation", info = "reduced gene set")

  reduced <- annotate_broad_cells(reduced, label_col = "broad_enrichment_m4")

  reduced <- subcluster_labels(reduced,
                               group_col    = "broad_enrichment_m4",
                               feature_mode = "all",
                               out_col      = "broad_enrichment_sub_m4")

  method_4 <- .annotate(reduced, "broad_enrichment_sub_m4", rownames(reduced))

  sce$broad_enrichment_m4     <- reduced$broad_enrichment_m4
  sce$broad_enrichment_sub_m4 <- reduced$broad_enrichment_sub_m4

  S4Vectors::metadata(sce)$broad_cell_enrichment_reduced <- S4Vectors::metadata(reduced)$broad_cell_enrichment

  # Global tie-breakers ----

  cli::cat_line()
  print_h2("Global Tie-breaker 1", info = "full gene set")
  global_1 <- .annotate(sce, "cluster_broad_hvg", rownames(sce))

  cli::cat_line()
  print_h2("Global Tie-breaker 2", info = "reduced gene set")
  global_2 <- .annotate(reduced, "cluster", rownames(reduced))


  cli::cat_line()
  print_alert("Done: Pass the returned {.code $labels} to {.fn resolve_consensus_labels}.", type = "s")

  # --- 4.) Assemble output ----
  all_outputs <- list(
    method_1 = method_1,
    method_2 = method_2,
    method_3 = method_3,
    method_4 = method_4,
    global_1 = global_1,
    global_2 = global_2
  )

  out <- list(
    sce                =  sce,
    labels             =  lapply(all_outputs, `[[`, "labels"),
    method_names       =  c("method_1", "method_2", "method_3", "method_4"),
    tie_breaker_names  =  c("global_1", "global_2")
              )

  if (return_full_output) out$full_output <- all_outputs

  return(out)
}
