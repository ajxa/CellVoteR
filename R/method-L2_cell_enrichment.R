#' Assign cell labels using broad-marker enrichment
#'
#' Assigns broad labels to individual cells using aggregated log-normalised
#' expression across small curated marker sets for each broad category.
#'
#' For each category, expression is aggregated across that category's marker
#' genes using one of \code{"sum"}, \code{"mean"}, or \code{"median"}.
#' A category is considered to pass for a cell if the aggregated expression
#' exceeds the category-specific \code{expr_threshold}. If no categories pass,
#' the cell is labelled as \code{other_label}. If more than one category passes,
#' the tie is resolved using the user-supplied \code{priority}, with lower
#' numeric values treated as higher priority.
#'
#' This function assumes that \code{broad_config} has already been validated in
#' upstream preprocessing, including checks that all broad markers are present,
#' marker sets do not overlap, and priorities are unique.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing a \code{logcounts} assay.
#' @param marker_config_key Character. Metadata entry containing a named list of
#' validated marker definitions. Defaults to \code{"marker_config"}.
#' @param label_col Character scalar. Name of the output column to create in
#'   \code{colData(sce)}. Defaults to \code{"broad_enrichment"}.
#' @param assay_name Character scalar. Assay to use for expression values.
#'   Defaults to \code{"logcounts"}.
#' @param aggregate_fun Character scalar specifying how marker expression should
#'   be aggregated within each category. One of \code{"sum"}, \code{"mean"}, or
#'   \code{"median"}. Defaults to \code{"sum"}.
#' @param other_label Character scalar. Label assigned when no category passes.
#'   Defaults to \code{"other"}.
#'
#' @return The input \code{SingleCellExperiment} with:
#' \describe{
#'   \item{\code{colData(sce)[[label_col]]}}{A factor of per-cell broad labels,
#'     with \code{other_label} as the last level.}
#'   \item{\code{metadata(sce)$broad_cell_enrichment}}{A list containing the
#'     per-cell category score matrix, logical pass matrix, and assignment
#'     parameters.}
#' }
#'
#' @export
annotate_broad_cells <- function(sce,
                                 marker_config_key = "marker_config",
                                 label_col = "broad_enrichment",
                                 assay_name = "logcounts",
                                 aggregate_fun = c("sum", "mean", "median"),
                                 other_label = "other"
                                 ) {

  aggregate_fun <- match.arg(aggregate_fun)

  if (!inherits(sce, "SingleCellExperiment")) cli::cli_abort("{.arg sce} must be a {.cls SingleCellExperiment}.")

  if (!assay_name %in% SummarizedExperiment::assayNames(sce)) {
    cli::cli_abort(c(
      "Assay {.val {assay_name}} is missing from {.arg sce}.",
      "i" = "Available assays: {.val {SummarizedExperiment::assayNames(sce)}}"
    ))
  }

  if (!marker_config_key %in% names(S4Vectors::metadata(sce))) {
    cli::cli_abort(c(
      "Metadata entry {.val {marker_config_key}} is missing from {.arg sce}.",
      "i" = "Available metadata entries: {.val {names(S4Vectors::metadata(sce))}}"
    ))
  }

  broad_config <- S4Vectors::metadata(sce)[[marker_config_key]][["broad"]]

  if (!is.list(broad_config) || length(broad_config) == 0L) cli::cli_abort("{.arg broad_config} must be a non-empty named list.")

  checkmate::assert_string(label_col, min.chars = 1)
  checkmate::assert_string(assay_name, min.chars = 1)
  checkmate::assert_string(other_label, min.chars = 1)

  cat_markers <- lapply(broad_config, `[[`, "markers")
  category_names <- names(cat_markers)

  if (is.null(category_names) || any(category_names == "")) cli::cli_abort("{.arg broad_config} must be a named list.")

  all_markers <- unlist(cat_markers, use.names = FALSE)
  marker_groups <- rep(category_names, lengths(cat_markers))

  priority_vec <- vapply(broad_config, function(x) x$priority, numeric(1))
  expr_thresholds <- vapply(broad_config, function(x) x$expr_threshold, numeric(1))

  sub_expr <- SummarizedExperiment::assay(sce, assay_name)[all_markers, , drop = FALSE]

  if (aggregate_fun == "sum") {

    cat_expr <- t(scuttle::sumCountsAcrossFeatures(x = sub_expr, id = marker_groups))

  } else {

    split_idx <- split(seq_along(all_markers), marker_groups)

    cat_expr_list <- lapply(category_names, function(cat_name) {

      idx <- split_idx[[cat_name]]

      mat <- sub_expr[idx, , drop = FALSE]

      if (aggregate_fun == "mean") {

        if (inherits(mat, "sparseMatrix")) vals = Matrix::colMeans(mat) else vals = colMeans(mat)

      } else vals <- apply(as.matrix(mat), 2, stats::median)

    return(as.numeric(vals))

    })

    cat_expr <- do.call(cbind, cat_expr_list)
    colnames(cat_expr) <- category_names
    rownames(cat_expr) <- colnames(sce)

  }

  passed <- sweep(cat_expr, 2, expr_thresholds, FUN = ">")

  labels <- rep(other_label, nrow(cat_expr))

  for (i in seq_len(nrow(cat_expr))) {

    passed_cats <- colnames(cat_expr)[passed[i, ]]

    if (length(passed_cats) == 0L) next

    if (length(passed_cats) == 1L) labels[i] <- passed_cats; next

    tied_priorities <- priority_vec[passed_cats]
    winner <- passed_cats[order(tied_priorities, passed_cats)][1]
    labels[i] <- winner

  }

  non_other <- sort(unique(labels[labels != other_label]))
  label_levels <- c(non_other, other_label)

  sce[[label_col]] <- factor(labels, levels = label_levels)

  S4Vectors::metadata(sce)$broad_cell_enrichment <- list(
    scores = cat_expr,
    passed = passed,
    matched_markers = cat_markers,
    params = list(
      label_col = label_col,
      assay_name = assay_name,
      aggregate_fun = aggregate_fun,
      expr_thresholds = expr_thresholds,
      priority = priority_vec,
      other_label = other_label
    )
  )

  return(sce)
}
