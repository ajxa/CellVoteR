#' Estimate clustering parameters using the cell count
#'
#' Computes sensible defaults for the number of principal components, SNN
#' graph neighbourhood size (K), and Leiden clustering resolution based on
#' the number of cells. All three parameters scale dynamically with dataset
#' size using bounded transformations and are intended to allow for slight over-clustering.
#'
#' @section Parameter scaling logic:
#' \describe{
#'   \item{\strong{PCs}}{Log2-scaled between \code{npc_min} and
#'     \code{npc_max}. Small datasets get fewer components to avoid
#'     overfitting noise; large datasets saturate at \code{npc_max}.}
#'   \item{\strong{K (SNN neighbours)}}{Square-root-scaled between
#'     \code{k_min} and \code{k_max}. Grows slowly because neighbourhood
#'     size affects graph topology - too large and distinct populations merge.}
#'   \item{\strong{Resolution (Leiden)}}{Square-root-scaled between
#'     \code{res_min} and \code{res_max}. Grows faster than K because
#'     larger datasets typically contain more distinct communities and
#'     require finer partitioning.}
#' }
#'
#' @param n_cells Integer scalar. Number of cells in the dataset.
#' @param min_cluster_cells Integer scalar. Minimum number of cells required
#'   to attempt clustering. If \code{n_cells} is below this threshold the
#'   function returns a parameter set with \code{skip = TRUE}. Defaults to
#'   \code{50}.
#' @param res_min,res_max Numeric scalars. Bounds for Leiden resolution.
#'   Defaults to \code{0.6} and \code{2.0}.
#' @param res_saturation_n Integer scalar. Cell count at which resolution
#'   reaches \code{res_max}. Defaults to \code{100000}.
#' @param k_min,k_max Integer scalars. Bounds for SNN nearest-neighbour K.
#'   Defaults to \code{10} and \code{50}.
#' @param k_saturation_n Integer scalar. Cell count at which K reaches
#'   \code{k_max}. Defaults to \code{200000}.
#' @param npc_min,npc_max Integer scalars. Bounds for number of principal
#'   components. Defaults to \code{20} and \code{50}.
#' @param npc_slope Numeric scalar. Controls the log2-based scaling rate
#'   for PCs. Defaults to \code{4}.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{n_pcs}{Integer. Number of principal components.}
#'     \item{dims}{Integer vector \code{seq_len(n_pcs)}.}
#'     \item{k}{Integer. Nearest-neighbour count for SNN graph.}
#'     \item{resolution}{Numeric. Leiden clustering resolution.}
#'     \item{skip}{Logical. \code{TRUE} if \code{n_cells} is below the
#'       minimum threshold.}
#'   }
#'
#' @examples
#' \dontrun{
#' params <- estimate_cluster_params(n_cells = ncol(sce))
#' if (!params$skip) {
#'   ctx <- prepare_context(sce, n_pcs = params$n_pcs, k = params$k)
#' }
#' }
#'
#' @export
estimate_cluster_params <- function(n_cells,
                                    min_cluster_cells = 50L,
                                    res_min = 0.6,
                                    res_max = 2.0,
                                    res_saturation_n = 100000L,
                                    k_min = 10L,
                                    k_max = 50L,
                                    k_saturation_n = 200000L,
                                    npc_min = 20L,
                                    npc_max = 50L,
                                    npc_slope = 4
                                    ) {
  # --- Validation ----

  checkmate::assert_count(n_cells, positive = TRUE)
  checkmate::assert_count(min_cluster_cells, positive = TRUE)
  checkmate::assert_number(res_min, lower = 0)
  checkmate::assert_number(res_max, lower = res_min)
  checkmate::assert_count(res_saturation_n, positive = TRUE)
  checkmate::assert_count(k_min, positive = TRUE)
  checkmate::assert_count(k_max, positive = TRUE)
  checkmate::assert_count(k_saturation_n, positive = TRUE)
  checkmate::assert_count(npc_min, positive = TRUE)
  checkmate::assert_count(npc_max, positive = TRUE)
  checkmate::assert_number(npc_slope, lower = 0)

  if (n_cells < min_cluster_cells) {

    print_alert(
      "Only {.val {n_cells}} cells - below minimum of {.val {min_cluster_cells}} for clustering",
      type = "w"
    )

    return(
      list(
        n_pcs      = 0L,
        dims       = integer(0),
        k          = 0L,
        resolution = 0,
        skip       = TRUE
        )
      )
  }

  # --- Dynamic PCs ----

  intercept  <- round((npc_slope * log2(min_cluster_cells)) - npc_min)
  calc_npcs  <- round(npc_slope * log2(n_cells) - intercept)
  final_npcs <- as.integer(max(npc_min, min(npc_max, calc_npcs)))

  # --- Dynamic K ----

  # K controls graph topology: too small = fragmented/noisy, too large =
  # over-smoothed. Grows slowly with sqrt because neighbourhood size has
  # diminishing returns - doubling cell count doesn't mean each cell needs
  # twice as many neighbours.

  k_ratio  <- sqrt(n_cells) / sqrt(k_saturation_n)
  calc_k   <- k_min + (k_ratio * (k_max - k_min))
  final_k  <- as.integer(round(max(k_min, min(k_max, calc_k))))

  final_k <- min(final_k, n_cells - 1L)

  # --- Dynamic clustering resolution ----

  # Resolution controls partition granularity within the graph. Scales faster
  # than K because larger datasets typically harbour more distinct communities.

  res_ratio  <- sqrt(n_cells) / sqrt(res_saturation_n)
  calc_res   <- res_min + (res_ratio * (res_max - res_min))
  final_res  <- round(max(res_min, min(res_max, calc_res)), 2)


  params <- list(
    n_pcs      = final_npcs,
    dims       = seq_len(final_npcs),
    k          = final_k,
    resolution = final_res,
    skip       = FALSE
  )

  return(params)
}


#' Select highly variable genes
#'
#' Models per-gene variance using
#' \code{\link[scran]{modelGeneVar}} and selects the top \code{n_hvgs}
#' genes by biological variance component.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{logcounts} assay.
#' @param n_hvgs Integer scalar. Number of HVGs to select.
#' @param BPPARAM A \code{BiocParallelParam} for parallelisation.
#'
#' @return Character vector of selected gene names.
#'
#' @keywords internal
#' @noRd
select_hvgs <- function(sce,
                        n_hvgs = 2000L,
                        BPPARAM = BiocParallel::SerialParam()) {

  gene_var <- scran::modelGeneVar(sce, BPPARAM = BPPARAM)

  n_available <- nrow(gene_var)

  if (n_hvgs > n_available) {

    print_alert("Requested {.val {n_hvgs}} HVGs but only {.val {n_available}} genes available. Using all.", "w")
    n_hvgs <- n_available

  }

  hvgs <- scran::getTopHVGs(gene_var, n = n_hvgs)

  return(hvgs)
}



#' Normalise counts via pooling-based size factors
#'
#' Computes size factors using \code{\link[scuttle]{computePooledFactors}} and
#' applies log-normalisation via \code{\link[scuttle]{logNormCounts}}. If a
#' \code{logcounts} assay already exists it is overwritten with a warning.
#'
#' @param sce A \code{SingleCellExperiment} with a \code{counts} assay.
#' @param min_cluster_size Integer scalar. Minimum cluster size for
#'   \code{\link[scran]{quickCluster}} used to seed the pooling. Defaults to \code{100}.
#' @param BPPARAM A \code{BiocParallelParam} for parallelisation.
#'
#' @return The input SCE with \code{sizeFactors()} set and a \code{logcounts} assay added.
#' @export
normalize_counts <- function(sce,
                             min_cluster_size = 100L,
                             BPPARAM = BiocParallel::SerialParam()) {

  if ("logcounts" %in% SummarizedExperiment::assayNames(sce)) print_alert("Overwriting existing {.val logcounts} assay.", type = "w")

  print_alert("Normalizing cells...")
  print_alert("Computing size factors via pooling...")

  min_cluster_size <- min(min_cluster_size, max(10L, floor(ncol(sce) / 2)))

  quick_clust <- scran::quickCluster(
    sce,
    min.size  = min_cluster_size,
    BPPARAM   = BPPARAM
  )

  # Lun et al. (2016) deconvolution strategy for scaling normalization of sparse count data.
  sce <- scuttle::computePooledFactors(
    x = sce,
    clusters = quick_clust,
    BPPARAM  = BPPARAM
  )

  sf <- BiocGenerics::sizeFactors(sce)

  pos_sf = sf[is.finite(sf) & sf > 0]

  if (length(pos_sf) == 0L) {
    cli::cli_abort(c(
      "Failed to compute valid positive size factors.",
      "i" = "Try a simpler normalization fallback for this dataset."
    ))
  }

  sf[!is.finite(sf) | sf <= 0] = min(pos_sf)
  BiocGenerics::sizeFactors(sce) = sf

  sce <- scuttle::logNormCounts(sce)

  print_alert("Done: size factors range [{.val {round(min(sf), 3)}} - {.val {round(max(sf), 3)}}]", type = "s")

  return(sce)
}


#' Validate broad marker configuration structure and genes
#'
#' Checks that the \code{broad} element of a marker configuration has been
#' processed by \code{\link{build_broad_marker_config}}, verifies that all
#' broad marker genes are present in the expression matrix, confirms
#' priority ranks are unique, and checks that no marker genes are shared
#' across broad categories.
#'
#' @section Validation steps:
#' \enumerate{
#'   \item \strong{Configuration structure}: Each entry must be a list
#'     containing \code{markers}, \code{expr_threshold}, \code{coexp_min},
#'     and \code{priority}. Raw character vectors (from
#'     \code{\link{load_markers}}) are detected and produce a specific error
#'     with the required fix.
#'   \item \strong{Gene presence (strict)}: Every marker gene in every broad
#'     category must exist in \code{available_genes}. Broad marker sets are
#'     typically small (2-5 genes) and cannot tolerate missing members.
#'   \item \strong{No overlapping markers}: Marker genes must be unique
#'     across broad categories. Shared markers would compromise the ability
#'     to discriminate between broad cell types.
#'   \item \strong{Priority uniqueness}: Duplicate priority ranks produce a
#'     warning, as they may cause non-deterministic tie-breaking during
#'     broad assignment.
#' }
#'
#' @param broad_config The \code{$broad} element of a marker configuration.
#' @param available_genes Character vector of gene names in the SCE
#'   (\code{rownames(sce)}).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{distinct_markers}{Character vector of all unique broad marker
#'       genes (guaranteed present in the expression matrix).}
#'     \item{n_markers}{Integer. Total number of distinct broad markers.}
#'   }
#'
#' @keywords internal
#' @noRd
validate_broad_markers <- function(broad_config, available_genes) {

  # --- Check `build_broad_marker_config` has been run ----

  if (!is.list(broad_config) || length(broad_config) == 0L) cli::cli_abort("{.arg marker_config$broad} must be a non-empty named list.")

  required_fields <- c("markers", "expr_threshold", "coexp_min", "priority")

  is_configured <- vapply(
    X = broad_config,
    FUN =  function(entry) {is.list(entry) && all(required_fields %in% names(entry))},
    FUN.VALUE =  logical(1)
  )

  if (!all(is_configured)) {

    unconfigured <- names(broad_config)[!is_configured]
    is_raw <- vapply(broad_config, is.character, logical(1))

    if (any(is_raw)) {
      cli::cli_abort(c(
        "{.arg marker_config$broad} has not been processed by {.fn build_broad_marker_config}.",
        "i" = "The raw marker list from {.fn load_markers} must be configured before
               use in the annotation pipeline.",
        "i" = "Run:",
        " " = "{.code marker_config$broad <- build_broad_marker_config(marker_config$broad, priority_order = c(...))}"
      ))
    }

    cli::cli_abort(c(
      "Malformed entries in {.arg marker_config$broad}.",
      "x" = "The following categories are missing required fields: {.val {unconfigured}}",
      "i" = "Each entry must be a list with: {.val {required_fields}}",
      "i" = "These are produced by {.fn build_broad_marker_config}."
    ))
  }

  # --- Check overlapping markers across categories ----

  all_markers_by_cat <- lapply(broad_config, `[[`, "markers")
  all_markers_flat   <- unlist(all_markers_by_cat, use.names = FALSE)

  if (any(duplicated(all_markers_flat))) {

    dup_genes <- unique(all_markers_flat[duplicated(all_markers_flat)])

    # Identify which categories share each duplicated gene
    overlap_detail <- vapply(
      X = dup_genes,
      FUN = function(gene) {
        cats <- names(all_markers_by_cat)[vapply(all_markers_by_cat, function(m) gene %in% m, logical(1))]
        paste(cats, collapse = ", ")
      },
      FUN.VALUE = character(1)
    )

    msg <- c("Broad marker genes must be unique across categories.")

    for (i in seq_along(dup_genes)) {
      msg <- c(msg, "x" = "{.val {dup_genes[i]}} shared by: {overlap_detail[i]}")
    }

    msg <- c(msg,
             "i" = "Overlapping markers compromise broad category discrimination.",
             "i" = "Please assign each marker gene to a single broad category."
    )

    cli::cli_abort(msg)
  }

  distinct_markers <- unique(all_markers_flat)

  # --- Check all broad marker genes are present ----

  all_missing <- list()

  for (cat in names(broad_config)) {

    markers <- broad_config[[cat]]$markers

    missing <- setdiff(markers, available_genes)

    if (length(missing) > 0L) all_missing[[cat]] <- missing

  }

  if (length(all_missing) > 0L) {

    msg <- c("Broad marker genes missing from expression matrix.")

    for (cat in names(all_missing)) {

      n_missing <- length(all_missing[[cat]])
      n_total   <- length(broad_config[[cat]]$markers)

      msg <- c(msg,
               "x" = "{.val {cat}}: {.val {n_missing}}/{.val {n_total}} missing - {.val {all_missing[[cat]]}}"
      )
    }

    msg <- c(msg,
             "i" = "All broad markers must be present. These small gene sets cannot tolerate missing members.",
             "i" = "Please check that marker gene names match the identifiers in your expression matrix."
    )

    cli::cli_abort(msg)
  }

  # --- Check priorities are unique ----

  priorities <- vapply(broad_config, `[[`, integer(1), "priority")

  if (any(duplicated(priorities))) {

    dup_ranks <- unique(priorities[duplicated(priorities)])
    dup_cats  <- names(priorities)[priorities %in% dup_ranks]

    cli::cli_warn(c(
      "!" = "Broad categories share the same priority rank: {.val {dup_cats}}",
      "i" = "This may cause non-deterministic tie-breaking during broad assignment.",
      "i" = "Consider re-running {.fn build_broad_marker_config} with a complete {.arg priority_order}."
    ))

  }

  return(
    list(
      distinct_markers = distinct_markers,
      n_markers        = length(distinct_markers)
    )
  )
}


#' Validate fine markers against available genes
#'
#' Checks fine marker overlap with the expression matrix. Aborts if no fine
#' markers are found at all and warns about partially missing markers with per-label detail,
#' ordered from the most to least affected label for easy triage.
#'
#' @section Return value:
#' The function returns its computed overlap statistics so that
#' \code{\link{prepare_context}} can reuse them for reduced feature space
#' selection without redundant calculation.
#'
#' @param fine_config Fine marker configuration from
#'   \code{\link{load_markers}} (\code{markers$fine}). A nested named list:
#'   top level keyed by broad category, second level keyed by fine label,
#'   values are character vectors of marker gene symbols.
#' @param available_genes Character vector of gene names in the SCE
#'   (\code{rownames(sce)}).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{total_distinct}{Character. The distinct fine marker genes (including any missing).}
#'     \item{total_present}{Character. The distinct fine markers present in \code{available_genes}.}
#'     \item{total_missing}{Character. The distinct fine markers missing from \code{available_genes}.}
#'     \item{filtered_fine}{Nested named list of fine markers with missing genes removed, structured the same as \code{fine_config}.}
#'     \item{overlap_pct}{Numeric. The total percentage overlap of \code{present_markers} within \code{available_genes}.}
#'     \item{missing_by_label}{Named list of missing markers keyed by \code{"category > label"}, ordered from highest to lowest percentage missing. Empty list if none missing.}
#'   }
#'
#' @keywords internal
#' @noRd
validate_fine_markers <- function(fine_config, available_genes) {

  # --- Marker statistics ----

  distinct_fine_markers <- unique(unlist(fine_config, use.names = FALSE))

  if (length(distinct_fine_markers) == 0L) {
    cli::cli_abort(c(
      "No fine markers found in {.arg marker_config$fine}.",
      "i" = "The fine marker configuration is empty or malformed.",
      "i" = "Please check the output of {.fn load_markers}."
    ))
  }

  distinct_present  <- intersect(distinct_fine_markers, available_genes)
  distinct_missing  <- setdiff(distinct_fine_markers, available_genes)

  total_found <- length(distinct_present)
  overlap_pct <- round(total_found / length(distinct_fine_markers) * 100, 2)

  if (total_found == 0L) {
    cli::cli_abort(c(
      "None of the {.val {length(distinct_fine_markers)}} fine marker genes were found in the expression matrix.",
      "i" = "Please check that marker gene names match the identifiers in your expression matrix."
    ))
  }

  # --- filter out missing markers from fine_config for downstream use ----

  if(length(distinct_missing > 1)){

    filtered_fine <- lapply(
      X = fine_config,
      FUN = function(x) lapply(x, function(y) y[!y %in% distinct_missing])
    )
  } else filtered_fine <- fine_config


  # --- Compute per label statistics ----

  compute_per_label_stats <- function(markers, available_genes) {

    missing <- setdiff(markers, available_genes)

    present <- intersect(markers, available_genes)

    out <- list(
        total        =  length(markers),
        found        =  length(present),
        missing      =  missing,
        pct_missing  =  if (length(markers) > 0L) round((length(missing) / length(markers)) * 100, 2) else 0
      )

    return(out)
  }

  label_stats <- purrr::list_flatten(fine_config, name_spec = "{outer} > {inner}")

  label_stats <- lapply(X = label_stats,
                        FUN = compute_per_label_stats,
                        available_genes = available_genes
                        )

  missing_labels <- vapply(X = label_stats,
                           FUN = function(x) length(x[["missing"]]) > 0,
                           FUN.VALUE = logical(1)
                           )

  n_labels_affected <- length(which(missing_labels))

  if(n_labels_affected > 0L) {

    missing_by_label <- label_stats[missing_labels]

    desc_ordering <- order(
      vapply(missing_by_label, `[[`, numeric(1), "pct_missing"),
      decreasing = TRUE
      )

    missing_by_label <- missing_by_label[desc_ordering]

    } else missing_by_label <- list()

  out <- list(
              total_distinct   =  distinct_fine_markers,
              total_present    =  distinct_present,
              total_missing    =  distinct_missing,
              filtered_fine    =  filtered_fine,
              missing_by_label =  missing_by_label,
              overlap_pct      =  overlap_pct
              )

  return(out)
}


#' @keywords internal
#' @noRd
validate_marker_config <- function(marker_config) {

  if (!is.list(marker_config)) cli::cli_abort("{.arg marker_config} must be a list.")

  missing <- setdiff(c("broad", "fine"), names(marker_config))

  if (length(missing) > 0L) {
    cli::cli_abort(c(
      "{.arg marker_config} is missing: {.val {missing}}",
      "i" = "Expected output of {.fn load_markers} with {.code $broad} processed by {.fn build_broad_marker_config}."
    ))
  }
}


#' Prepare Analysis Tracks on a SingleCellExperiment
#'
#' Runs the complete preprocessing pipeline for both feature spaces required
#' for CellVoteR directly on the SCE object using standard Bioconductor functions.
#'
#' @section Storage layout:
#'
#' \preformatted{
#' sce
#' |-- assays: counts, logcounts
#' |-- rowSubset("broad_hvgs")
#' |-- reducedDim("PCA_broad_hvg")
#' |-- colData$cluster_broad_hvg
#' |-- metadata$broad_hvg_params
#' |-- metadata$marker_config
#' |-- metadata$filterd_fine_markers
#' |-- metadata$missing_by_label
#' `-- altExp("user_panel")
#'     |-- assays: counts, logcounts
#'     |-- reducedDim("PCA")
#'     |-- colData$cluster
#'     |-- metadata$marker_config
#'     |-- metadata$filterd_fine_markers
#'     `-- metadata$params
#' }
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   with \code{logcounts} assay (from \code{\link{normalize_counts}}).
#' @param marker_config Named list with \code{$broad} (from
#'   \code{\link{build_broad_marker_config}}) and \code{$fine} (from
#'   \code{\link{load_markers}}).
#' @param n_hvgs Integer scalar. Number of HVGs. Defaults to \code{2000}.
#' @param overlap_feat_percent Numeric scalar (0-100). Minimum fine marker
#'   overlap. Defaults to \code{50}.
#' @param n_pcs,k,resolution Override automatic parameter estimation.
#'   \code{NULL} (default) uses \code{\link{estimate_cluster_params}}.
#' @param cluster_params_args Named list passed to
#'   \code{\link{estimate_cluster_params}}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}}.
#'   Defaults to \code{\link[BiocParallel]{SerialParam}()}.
#'
#' @return The input SCE with both analysis tracks populated. The \code{marker_config}
#' that is stored in the \code{metadata} of the \code{altExp} differs from that stored
#' in the main SCE in that the fine marker sets are filtered to only include
#' the genes present in the user panel and that overlap with the broad marker sets. The broad HVG track is built using the full \code{logcounts} assay and the full set of HVGs, while the user panel track is built using a subset of features based on the fine markers that overlap with the user panel and broad markers. The dimensionality reduction and clustering for each track are performed independently using their respective feature spaces, and their parameters are estimated separately using \code{\link{estimate_cluster_params}}.
#'
#' @examples
#' \dontrun{
#' sce <- normalize_counts(sce)
#' sce <- prepare_sce(sce, marker_config)
#'
#' # Broad HVG track
#' reducedDim(sce, "PCA_broad_hvg")
#' sce$cluster_broad_hvg
#'
#' # User panel track
#' reducedDim(altExp(sce, "user_panel"), "PCA")
#' altExp(sce, "user_panel")$cluster
#' }
#'
#' @export
prepare_sce <- function(sce,
                        marker_config,
                        n_hvgs = 2000L,
                        overlap_feat_percent = 50,
                        n_pcs = NULL,
                        k = NULL,
                        resolution = NULL,
                        cluster_params_args = list(),
                        BPPARAM = BiocParallel::SerialParam()) {

  if (!inherits(sce, "SingleCellExperiment")) cli::cli_abort("{.arg sce} must be a {.cls SingleCellExperiment}.")

  if (anyDuplicated(rownames(sce))) {
    cli::cli_abort(c(
      "Feature names in {.arg sce} must be unique.",
      "i" = "Duplicate row names make marker matching and feature tracking ambiguous."
    ))
  }

  if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
    cli::cli_abort(c(
      "SCE is missing a {.val logcounts} assay.",
      "i" = "Please run {.fn normalize_counts} before {.fn prepare_sce}."
    ))
  }

  checkmate::assert_number(overlap_feat_percent, lower = 0, upper = 100)

  #--- 2.) Marker universe ----

  validate_marker_config(marker_config)

  available_genes  <- rownames(sce)
  broad_validation <- validate_broad_markers(marker_config$broad, available_genes)
  fine_validation  <- validate_fine_markers(marker_config$fine, available_genes)

  if (fine_validation$overlap_pct < overlap_feat_percent) {
    cli::cli_abort(
      c(
        "x" = "Insufficient fine marker overlap ({.val {fine_validation$overlap_pct}}%).",
        "i" = "Minimum {.val {overlap_feat_percent}}% required."
      ),
      missing_feats    = fine_validation$total_missing,
      missing_by_label = fine_validation$missing_by_label
    )
  }

  if(length(fine_validation$missing_by_label) > 0){

    total_missing <- length(fine_validation$total_missing)
    n_labels_affected <- length(fine_validation$missing_by_label)

    missing_pct <- round(100 - fine_validation$overlap_pct, 0)

    print_alert(
      c(
        cli::style_bold(cli::col_yellow("Warning: ")),
        cli::style_italic("{.val {total_missing}} (~{missing_pct}%) fine markers across {.val {n_labels_affected}} labels are missing"),
        "\n",
        cli::col_blue(cli::symbol$info),
        " This may affect user panel clustering and wider annotation.",
        "\n - Will filter out missing markers for fisher annotation.",
        "\n - For more info, inspect:",
        cli::style_italic(cli::col_grey('{.code metadata(<sce>)$missing_by_label}'))
      ),
      type = "w",
      highlight_color = "white"
    )

    cli::cat_line()
  }

  S4Vectors::metadata(sce)$marker_config <- marker_config
  S4Vectors::metadata(sce)$filtered_fine_markers <- fine_validation$filtered_fine
  S4Vectors::metadata(sce)$missing_by_label <- fine_validation$missing_by_label

  #--- 3.) Setting params ----

  est_args <- c(list(n_cells = ncol(sce)), cluster_params_args)
  est      <- do.call(estimate_cluster_params, est_args)

  if (est$skip) cli::cli_abort("Too few cells ({.val {ncol(sce)}}) for clustering.")

  final_n_pcs      <- n_pcs      %||% est$n_pcs
  final_k          <- k          %||% est$k
  final_resolution <- resolution %||% est$resolution


  checkmate::assert_int(final_n_pcs, lower = 1)
  checkmate::assert_int(final_k, lower = 1)
  checkmate::assert_number(final_resolution, lower = 0)

  final_k <- min(final_k, ncol(sce) - 1L)

  #--- 4.) Building individual tracks ----

  print_alert("Preparing broad HVG track...")

  sce <- build_broad_hvg_track(sce,
                               n_hvgs = n_hvgs,
                               n_pcs = final_n_pcs,
                               k = final_k,
                               resolution = final_resolution,
                               broad_cluster_name = "cluster_broad_hvg",
                               BPPARAM = BPPARAM
                               )

  print_alert("Preparing user panel track...")

  user_panel <- build_user_panel_track(sce,
                                       broad_validation = broad_validation,
                                       fine_validation = fine_validation,
                                       panel_marker_config = marker_config,
                                       n_pcs = final_n_pcs,
                                       k = final_k,
                                       resolution = final_resolution,
                                       BPPARAM = BPPARAM
                                       )

  SingleCellExperiment::altExp(sce, "user_panel") <- user_panel

  print_alert("All tracks prepared", type = "s")

  return(sce)
}


#' @keywords internal
#' @noRd
run_pca_and_cluster = function(sce,
                               n_pcs,
                               k,
                               resolution,
                               dimred_name,
                               cluster_col,
                               subset_row = NULL,
                               BPPARAM = BiocParallel::SerialParam()
                               ) {

  if (is.null(subset_row)){

    n_features <- nrow(sce)

  } else if (is.logical(subset_row)) {

    n_features <- sum(subset_row, na.rm = TRUE)

  } else n_features <- length(subset_row)

  max_pcs <- min(n_pcs, n_features - 1L, ncol(sce) - 1L)

  if (max_pcs < 1L) {
    cli::cli_abort(c(
      "Cannot run PCA for.",
      "i" = "Need at least 2 cells and 2 selected features."
    ))
  }

  if (max_pcs < n_pcs) print_alert("Capped PCA at {.val {max_pcs}} components", type = "w")

  sce <- scran::fixedPCA(
    x = sce,
    subset.row = subset_row,
    rank = max_pcs,
    name = dimred_name,
    BPPARAM = BPPARAM
  )

  k_eff <- min(k, ncol(sce) - 1L)

  clusters <- scran::clusterCells(
    x = sce,
    use.dimred = dimred_name,
    BLUSPARAM = bluster::SNNGraphParam(
      k = k_eff,
      cluster.fun = "leiden",
      cluster.args = list(
        resolution_parameter = resolution,
        objective_function = "modularity"
      )
    )
  )

  sce[[cluster_col]] <- factor(clusters)

  # print_alert("   {.val {nlevels(sce[[cluster_col]])}} clusters found")

  out <- list(
    sce = sce,
    n_pcs = max_pcs,
    k = k_eff,
    resolution = resolution
  )

  return(out)
}


#' @keywords internal
#' @noRd
build_broad_hvg_track <- function(sce,
                                  n_hvgs,
                                  n_pcs,
                                  k,
                                  resolution,
                                  broad_cluster_name = "cluster_broad_hvg",
                                  BPPARAM
                                  ) {

  hvgs <- select_hvgs(sce, n_hvgs, BPPARAM)

  SingleCellExperiment::rowSubset(sce, "broad_hvgs") <- rownames(sce) %in% hvgs

  res <- run_pca_and_cluster(
    sce = sce,
    n_pcs = n_pcs,
    k = k,
    resolution = resolution,
    dimred_name = "PCA_broad_hvg",
    cluster_col = broad_cluster_name,
    subset_row = hvgs,
    BPPARAM = BPPARAM
  )

  sce <- res$sce

  S4Vectors::metadata(sce)$broad_hvg_params <- list(
    n_hvgs = length(hvgs),
    n_pcs = res$n_pcs,
    k = res$k,
    resolution = res$resolution
  )

  return(sce)
}


#' @keywords internal
#' @noRd
make_clean_panel_sce <- function(sce,
                                 markers_present,
                                 drop_coldata = c("cluster_broad_hvg"),
                                 drop_rowdata = c("broad_hvgs")
                                 ) {

  panel_sce <- sce[markers_present, ]

  S4Vectors::metadata(panel_sce) <- list()

  SingleCellExperiment::altExps(panel_sce) <- S4Vectors::SimpleList()
  SingleCellExperiment::reducedDims(panel_sce) <- S4Vectors::SimpleList()

  existing_col <- colnames(SummarizedExperiment::colData(panel_sce))
  keep_col <- setdiff(existing_col, intersect(existing_col, drop_coldata))

  SummarizedExperiment::colData(panel_sce) <-
    SummarizedExperiment::colData(panel_sce)[, keep_col, drop = FALSE]

  existing_row <- colnames(SummarizedExperiment::rowData(panel_sce))
  keep_row <- setdiff(existing_row, intersect(existing_row, drop_rowdata))

  SummarizedExperiment::rowData(panel_sce) <-
    SummarizedExperiment::rowData(panel_sce)[, keep_row, drop = FALSE]

  return(panel_sce)
}


#' @keywords internal
#' @noRd
build_user_panel_track <- function(sce,
                                   broad_validation,
                                   fine_validation,
                                   panel_marker_config,
                                   n_pcs,
                                   k,
                                   resolution,
                                   BPPARAM
                                   ) {

  all_distinct <- unique(c(broad_validation$distinct_markers,
                           fine_validation$total_distinct))

  markers_present <- unique(c(broad_validation$distinct_markers,
                              fine_validation$total_present)
                            )

  markers_missing  <- setdiff(all_distinct, markers_present)

  panel_sce <- make_clean_panel_sce(sce = sce,
                                    markers_present = markers_present,
                                    drop_coldata = c("cluster_broad_hvg"),
                                    drop_rowdata = c("broad_hvgs")
                                    )

  res <- run_pca_and_cluster(sce = panel_sce,
                             n_pcs = n_pcs,
                             k = k,
                             resolution = resolution,
                             dimred_name = "PCA",
                             cluster_col = "cluster",
                             subset_row = NULL,
                             BPPARAM = BPPARAM
                             )

  panel_sce <- res$sce

  S4Vectors::metadata(panel_sce)[["marker_config"]] <- panel_marker_config
  S4Vectors::metadata(panel_sce)[["filtered_fine"]] <- fine_validation$filtered_fine

  S4Vectors::metadata(panel_sce)[["params"]] <- list(
    input_markers = all_distinct,
    matched_markers = markers_present,
    missing_markers = markers_missing,
    n_pcs = res$n_pcs,
    k = res$k,
    resolution = res$resolution
    )

  return(panel_sce)
}
