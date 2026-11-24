#' Resolve same-category ties left by the consensus vote
#'
#' A secondary step applied after \code{\link{resolve_consensus_labels}()}.
#' Cells that the consensus left unassigned are re-examined: if every method
#' voted for a label from the \emph{same} category (for example
#' \code{GPC, GPC, OPC, OPC}, all \code{"cancer"}), the cell is given a
#' category-level label such as \code{"cancer_unresolved"}. The category is
#' agreed on; only the fine label (state) is not, and the output keeps that
#' visible rather than guessing a state.
#'
#' @details
#' \strong{Decision rule (per cell):}
#' \enumerate{
#'   \item If the consensus label is not \code{unassigned_label}, it is kept.
#'   \item Otherwise the method votes are mapped to categories with
#'     \code{label_category}. If every vote is present, every vote maps to a
#'     known category, and all categories are identical, the cell is labelled
#'     \code{paste0(category, suffix)}.
#'   \item In all other cases (votes split across categories, a missing vote,
#'     or a label with no or an ambiguous category) the cell stays
#'     \code{unassigned_label}.
#' }
#'
#' Cross-category ties (e.g. \code{Astrocyte, Astrocyte, GPC, GPC}, normal
#' versus cancer) are deliberately left unassigned: in such cells the
#' malignant/normal status itself is in doubt.
#'
#' The step has no fitted parameters. The category for each label comes from
#' the \code{category} column of the marker table, via the
#' \code{label_category} element returned by \code{\link{load_markers}()}.
#'
#' @param consensus The list returned by \code{\link{resolve_consensus_labels}()}
#'   (elements \code{label} and \code{method}).
#' @param label_list Named list of per-method labels, as passed to
#'   \code{\link{resolve_consensus_labels}()} (e.g. \code{$labels} from
#'   \code{\link{run_cellvoter}()}).
#' @param method_names Character vector. Names within \code{label_list} of the
#'   voting methods (e.g. \code{$method_names} from \code{run_cellvoter()}).
#'   Tie-breakers are not used by this step.
#' @param label_category Named character vector mapping fine labels to their
#'   category, as returned in \code{$label_category} by
#'   \code{\link{load_markers}()}.
#' @param unassigned_label Character scalar. The label used for unresolved
#'   cells in \code{consensus}. Defaults to \code{"Unknown"}, matching
#'   \code{resolve_consensus_labels()}.
#' @param suffix Character scalar appended to the category name to form the
#'   output label. Defaults to \code{"_unresolved"}.
#'
#' @return A named list of three equal-length factors:
#'   \describe{
#'     \item{\code{label}}{The label after this step.}
#'     \item{\code{method}}{The decision rule, as in
#'       \code{resolve_consensus_labels()}, with \code{"Same-category tie"}
#'       for cells changed by this step.}
#'     \item{\code{label_before_tie_step}}{The consensus label before this
#'       step, so the step can be switched off downstream.}
#'   }
#'
#' @examples
#' \dontrun{
#' markers <- load_markers(marker_panels$GBM$GBM_RNA_reference_panel_v1)
#' # ... prepare_sce(), run_cellvoter() ...
#' res <- run_cellvoter(sce)
#'
#' consensus <- resolve_consensus_labels(
#'   label_list        = res$labels,
#'   method_names      = res$method_names,
#'   tie_breaker_names = res$tie_breaker_names
#' )
#'
#' final <- resolve_category_ties(
#'   consensus      = consensus,
#'   label_list     = res$labels,
#'   method_names   = res$method_names,
#'   label_category = markers$label_category
#' )
#'
#' table(final$method)
#' }
#'
#' @seealso \code{\link{resolve_consensus_labels}}, \code{\link{load_markers}}
#' @export
resolve_category_ties <- function(consensus,
                                  label_list,
                                  method_names,
                                  label_category,
                                  unassigned_label = "Unknown",
                                  suffix = "_unresolved") {

  # --- 1.) Input validation ----
  checkmate::assert_list(consensus, names = "named")
  if (!all(c("label", "method") %in% names(consensus))) {
    cli::cli_abort(
      "{.arg consensus} must contain {.val label} and {.val method}, as returned by {.fn resolve_consensus_labels}."
    )
  }
  checkmate::assert_list(label_list, names = "named", min.len = 1L)
  checkmate::assert_character(method_names, min.len = 1L, any.missing = FALSE)
  checkmate::assert_character(label_category, names = "named")
  checkmate::assert_string(unassigned_label, min.chars = 1L)
  checkmate::assert_string(suffix)

  missing_names <- setdiff(method_names, names(label_list))
  if (length(missing_names) > 0L) {
    cli::cli_abort(
      "The following names are absent from {.arg label_list}: {.val {missing_names}}."
    )
  }

  cons_label  <- as.character(consensus$label)
  cons_method <- as.character(consensus$method)
  n_cells     <- length(cons_label)

  lens <- lengths(label_list[method_names])
  if (any(lens != n_cells) || length(cons_method) != n_cells) {
    cli::cli_abort(
      "{.arg consensus} and every method in {.arg label_list} must have the same length. \\
       Found consensus: {.val {n_cells}}, methods: {.val {lens}}."
    )
  }

  # --- 2.) Map votes to categories ----
  vote_mat <- do.call(cbind, lapply(label_list[method_names], as.character))
  vote_mat[!is.na(vote_mat) & !nzchar(vote_mat)] <- NA_character_

  cat_mat <- matrix(unname(label_category[vote_mat]),
                    nrow = n_cells, ncol = length(method_names))

  unmapped <- setdiff(unique(vote_mat[!is.na(vote_mat)]), names(label_category))
  if (length(unmapped) > 0L) {
    print_alert(
      "{length(unmapped)} voted label{?s} not found in {.arg label_category} \\
       (cells with {?this/these} vote{?s} stay unassigned): {.val {unmapped}}",
      type = "w"
    )
  }

  # --- 3.) Same-category ties ----
  first_cat    <- cat_mat[, 1L]
  all_present  <- rowSums(is.na(cat_mat)) == 0L
  all_same     <- rowSums(cat_mat != first_cat, na.rm = TRUE) == 0L
  is_unassigned <- !is.na(cons_label) & cons_label == unassigned_label

  resolve <- is_unassigned & all_present & all_same

  out_label  <- cons_label
  out_method <- cons_method

  out_label[resolve]  <- paste0(first_cat[resolve], suffix)
  out_method[resolve] <- "Same-category tie"

  print_alert(
    "{sum(resolve)} of {sum(is_unassigned)} unassigned cell{?s} resolved to a category-level label.",
    type = "i"
  )

  list(
    label                 = factor(out_label),
    method                = factor(out_method),
    label_before_tie_step = factor(cons_label)
  )
}
