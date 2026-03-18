#' Resolve final cell labels via majority vote
#'
#' Combines multiple annotation inputs and two tie-breakers to determine a
#' single consensus label for each cell.
#'
#' @details
#' \strong{Decision hierarchy (evaluated per cell):}
#' \enumerate{
#'   \item \strong{Strong majority:} If any label receives strictly more votes
#'     than all others, and that count exceeds the majority threshold, it is
#'     selected immediately. The threshold behaviour is controlled by
#'     \code{allow_even_split}.
#'   \item \strong{Split — tie-breaker agreement:} If there is a genuine
#'     leading candidate (at least one label has more votes than the others),
#'     and both tie-breakers agree with each other and match that candidate,
#'     the label is assigned.
#'   \item \strong{Split — priority order} (when \code{ordered_tiebreak = TRUE},
#'     the default): Tie-breaker 1 is tried first against the leading candidate.
#'     If it does not match, tie-breaker 2 is tried.
#'   \item \strong{Split — either agrees} (when \code{ordered_tiebreak = FALSE}):
#'     Either tie-breaker agreeing with the leading candidate is sufficient,
#'     with no priority between them.
#'   \item \strong{Unresolved:} If every method disagrees (no leading
#'     candidate), or no tie-breaker can resolve the split,
#'     \code{unassigned_label} is assigned.
#' }
#'
#' @param label_list Named list of factors or character vectors, one entry per
#'   cell. Must contain all names in \code{method_names} and
#'   \code{tie_breaker_names}, and all elements must be the same length.
#' @param method_names Character vector. Names within \code{label_list}
#'   representing the primary voting methods.
#' @param tie_breaker_names Character vector of length 2. Names within
#'   \code{label_list} to use as ordered fallbacks.
#' @param unassigned_label Character scalar. Label assigned when no consensus
#'   is reached. Defaults to \code{"Unknown"}.
#' @param ordered_tiebreak Logical scalar. When \code{TRUE} (default),
#'   tie-breakers are consulted in order: tie-breaker 1 is tried before
#'   tie-breaker 2. When \code{FALSE}, either tie-breaker agreeing with the
#'   leading candidate is sufficient, with no priority between them.
#' @param allow_even_split Logical scalar. When \code{FALSE} (default), a
#'   majority requires strictly more than 50\% of method votes (e.g. 3 out of
#'   4). When \code{TRUE}, a label with more votes than any other is accepted
#'   as a majority even if it does not exceed 50\% (e.g. 2 out of 4 when all
#'   others have 1).
#'
#' @return A named list with two equal-length factors:
#'   \describe{
#'     \item{\code{label}}{The resolved cell type label.}
#'     \item{\code{method}}{The decision rule used to reach that label.}
#'   }
#'
#' @export
resolve_consensus_labels <- function(label_list,
                                     method_names,
                                     tie_breaker_names,
                                     unassigned_label = "Unknown",
                                     ordered_tiebreak = TRUE,
                                     allow_even_split = FALSE
                                     ) {
  # --- 1.) Input validation ----
  checkmate::assert_list(label_list,  names = "named", min.len = 1L)
  checkmate::assert_character(method_names,      min.len = 1L, any.missing = FALSE)
  checkmate::assert_character(tie_breaker_names, len = 2L,     any.missing = FALSE)
  checkmate::assert_string(unassigned_label, min.chars = 1L)
  checkmate::assert_flag(ordered_tiebreak)
  checkmate::assert_flag(allow_even_split)

  all_required  <- c(method_names, tie_breaker_names)
  missing_names <- setdiff(all_required, names(label_list))

  if (length(missing_names) > 0L) {
    cli::cli_abort(
      "The following names are absent from {.arg label_list}: {.val {missing_names}}."
    )
  }

  lengths <- lengths(label_list[all_required])
  if (length(unique(lengths)) != 1L) {
    cli::cli_abort(
      "All elements of {.arg label_list} must have the same length. \\
       Found lengths: {.val {lengths}}."
    )
  }

  # ---- 2.) format the labels ----
  method_mat <- do.call(cbind, lapply(label_list[method_names], as.character))
  tb1        <- as.character(label_list[[tie_breaker_names[1L]]])
  tb2        <- as.character(label_list[[tie_breaker_names[2L]]])

  n_cells   <- nrow(method_mat)
  n_methods <- length(method_names)

  out_label  <- character(n_cells)
  out_method <- character(n_cells)

  # --- 3.) voting logic ----
  for (i in seq_len(n_cells)) {

    votes <- method_mat[i, ]
    votes <- votes[!is.na(votes) & nzchar(votes)]

    t1 <- tb1[i];  t1_valid <- !is.na(t1) && nzchar(t1)
    t2 <- tb2[i];  t2_valid <- !is.na(t2) && nzchar(t2)

    # majority:
    vote_counts   <- table(votes)
    sorted_counts <- sort(vote_counts, decreasing = TRUE)
    max_votes     <- if (length(sorted_counts) > 0L) sorted_counts[[1L]] else 0L

    strict_plurality <- length(sorted_counts) <= 1L ||sorted_counts[[1L]] > sorted_counts[[2L]]

    is_majority <- strict_plurality && (
      if (allow_even_split) max_votes >= (n_methods / 2) else max_votes >  (n_methods / 2)
    )

    if (is_majority) {
      out_label[i]  <- names(sorted_counts)[[1L]]
      out_method[i] <- "Majority vote"
      next
    }

    # tie-breaking:

    has_leading_candidate <- strict_plurality && max_votes > 1L

    if (has_leading_candidate) {

      candidate <- names(sorted_counts)[[1L]]

      tbs_agree <- t1_valid && t2_valid && (t1 == t2)

      if (tbs_agree && t1 == candidate) {
        out_label[i]  <- candidate
        out_method[i] <- "Tie-breaker (both agree)"
        next
      }

      if (ordered_tiebreak) {

        if (t1_valid && t1 == candidate) {
          out_label[i]  <- candidate
          out_method[i] <- "Tie-breaker 1"
          next
        }

        if (t2_valid && t2 == candidate) {
          out_label[i]  <- candidate
          out_method[i] <- "Tie-breaker 2"
          next
        }

      } else {

        tb1_breaks <- t1_valid && t1 == candidate
        tb2_breaks <- t2_valid && t2 == candidate

        if (tb1_breaks || tb2_breaks) {
          out_label[i]  <- candidate
          out_method[i] <- "Tie-breaker (either agrees)"
          next
        }
      }
    }

    out_label[i]  <- unassigned_label
    out_method[i] <- "Unresolved"
  }

  final_labels <- list(label  = factor(out_label), method = factor(out_method))

  return(final_labels)
}
