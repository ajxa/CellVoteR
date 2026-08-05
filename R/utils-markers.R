#' Load and structure CellVoteR marker definitions
#'
#' Reads a marker definition file (csv, tab-separated txt, xlsx, rds) or receives an R
#' \code{\link[base]{data.frame}} directly and organises it into a hierarchical list suitable for two-tier (broad -> fine)
#' cell-type annotation pipelines. The returned broad markers are a simple
#' named list of character vectors intended to be passed to
#' \code{\link{build_broad_marker_config}()} as a subsequent step to attach
#' priority, threshold, and co-expression settings.
#'
#' @section Expected file layout:
#' The input file must contain at least four columns (names configurable):
#' \describe{
#'   \item{type}{Either \code{"broad"} or \code{"fine"}.}
#'   \item{category}{Broad cell-type category (e.g. \code{"immune"},
#'     \code{"vasculature"}). Every broad category must also appear as a
#'     category in the fine rows.}
#'   \item{label}{Fine-grained cell-type label (e.g. \code{"CD8_T"},
#'     \code{"endothelial"}).}
#'   \item{marker}{Gene symbol.}
#' }
#'
#' @section Supported input types:
#' \describe{
#'   \item{\code{.csv}}{Comma-separated values, read via \code{\link[utils]{read.csv}}.}
#'   \item{\code{.txt}}{Tab-separated values, read via \code{\link[utils]{read.delim}}.}
#'   \item{\code{.xlsx}}{Excel workbook, read via \code{\link[openxlsx]{read.xlsx}}.}
#'   \item{\code{.rds}}{A saved R object, read via \code{\link[base]{readRDS}}. This
#'     must contain a \code{data.frame} matching the expected layout.}
#'   \item{\code{data.frame}}{A \code{data.frame} passed directly via either the
#'     \code{file_path} or \code{markers} argument.}
#' }
#'
#' @section Category reconciliation:
#' Fine-type rows whose \code{category} does not match any broad category are
#' reassigned to \code{unnamed_broad_cat_label} (default \code{"other"}) with
#' an informational message. Broad categories that have \strong{no}
#' corresponding fine rows cause an error.
#'
#' @section Typical workflow when loading external marker files:
#' \preformatted{
#' markers <- load_markers("markers.csv")
#'
#' markers$broad <- build_broad_marker_config(
#'   marker_list            = markers$broad,
#'   priority_order         = c("vasculature", "immune"),
#'   per_category_overrides = list(
#'     immune = list(coexp_min = 2)
#'   )
#' )
#' }
#'
#' @section Using internal marker panels:
#' CellVoteR includes a set of curated marker panels for IDHwt glioblastoma (GBM) under the
#' \code{marker_panels$GBM} list entry. These panels can be passed directly as a
#' \code{data.frame} to \code{load_markers()}:
#' \preformatted{
#' markers <- load_markers(marker_panels$GBM$gbmap_neftel_full)
#'
#' # Or using the explicit markers parameter:
#' markers <- load_markers(markers = marker_panels$GBM$gbmap_neftel_full)
#'
#' }
#'
#' @param file_path Character scalar. Path to a \code{.csv}, \code{.txt}
#'   (tab-separated), or \code{.xlsx} file containing marker definitions.
#' @param markers Alternative argument for directly passing a \code{data.frame} containing
#'   marker definitions. If supplied, this overrides \code{file_path}.
#' @param unnamed_broad_cat_label Character scalar. Label assigned to fine-type
#'   categories that do not map to any broad category. Defaults to
#'   \code{"other"}.
#' @param type_col,cat_col,label_col,marker_col Character scalars giving the
#'   column names in the input file for type, category, label, and marker
#'   respectively.
#' @param unique_types Character vector of permitted values in \code{type_col}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{broad}{Named list of character vectors - one element per broad
#'     category, values are marker gene symbols. Pass this to
#'     \code{\link{build_broad_marker_config}()} to generate the full
#'     configuration.}
#'   \item{fine}{Named list of lists: top level keyed by category, second
#'     level keyed by label, values are character vectors of markers. The
#'     \code{unnamed_broad_cat_label} group (if any) is placed last.}
#' }
#'
#' @examples
#' \dontrun{
#' # Load from file
#' markers <- load_markers("markers.csv")

#' # Load from .rds file
#' markers <- load_markers("markers.rds")

#' # Load directly from data frame / internal package dataset
#' markers <- load_markers(marker_panels$GBM$gbmap_neftel_full)
#'
#' # Inspect raw broad markers before configuring
#' names(markers$broad)
#' markers$broad[["immune"]]
#'
#' # Then configure broad markers for annotation
#' markers$broad <- build_broad_marker_config(
#'   marker_list    = markers$broad,
#'   priority_order = c("vasculature", "immune")
#' )
#' }
#'
#' @seealso \code{\link{build_broad_marker_config}} for configuring the broad
#'   markers returned by this function.
#'
#' @export
#' @export
load_markers <- function(file_path = NULL,
                         markers = NULL,
                         unnamed_broad_cat_label = "other",
                         type_col = "type",
                         unique_types = c("broad", "fine"),
                         cat_col = "category",
                         label_col = "label",
                         marker_col = "marker") {

  checkmate::assert_string(unnamed_broad_cat_label)
  checkmate::assert_character(unique_types, min.len = 1L, any.missing = FALSE)

  # 1.) Resolve the input ----

  if (is.null(file_path) && is.null(markers)) {
    cli::cli_abort("Please supply either a file path or a data.frame to 'file_path' or 'markers'.")
  }

  if (!is.null(markers)) {

    checkmate::assert_data_frame(markers, min.rows = 1L)
    marker_list <- markers

  } else if (is.data.frame(file_path)) {

    checkmate::assert_data_frame(file_path, min.rows = 1L)
    marker_list <- file_path

  } else if (is.character(file_path)) {

    checkmate::assert_string(file_path)
    checkmate::assert_file_exists(file_path, access = "r")

    ext <- tolower(tools::file_ext(file_path))

    marker_list <- switch(
      EXPR = ext,
      csv  = utils::read.csv(file_path, stringsAsFactors = FALSE),
      txt  = utils::read.delim(file_path, stringsAsFactors = FALSE),
      xlsx = openxlsx::read.xlsx(file_path),
      rds  = readRDS(file_path),
      cli::cli_abort("Unsupported file extension '{ext}'. Expected csv, txt, xlsx, or rds.")
    )

    if (!is.data.frame(marker_list)) {
      cli::cli_abort("The object loaded from '{file_path}', is not a data.frame.")
    }
  } else cli::cli_abort("'file_path' must be a character string path or a data.frame object.")

  if (nrow(marker_list) == 0L) cli::cli_abort("Marker data is empty.")

  # 2.) Column validation ----

  required_cols <- c(type_col, cat_col, label_col, marker_col)
  missing_cols  <- setdiff(required_cols, colnames(marker_list))

  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         "\n  Found: ", paste(colnames(marker_list), collapse = ", "),
         call. = FALSE)
  }

  # 3.) Trim white-space & drop empty rows ----

  for (col in required_cols) {
    marker_list[[col]] <- trimws(marker_list[[col]])
  }

  marker_list <- marker_list[!is.na(marker_list[[marker_col]]) &
                               nchar(marker_list[[marker_col]]) > 0L, ]

  # 4.) Type validation ----

  observed_types <- unique(marker_list[[type_col]])
  unexpected_types <- setdiff(observed_types, unique_types)

  if (length(unexpected_types) > 0L) {
    stop("Unexpected values in '", type_col, "': ",
         paste(unexpected_types, collapse = ", "),
         "\n  Allowed: ", paste(unique_types, collapse = ", "),
         call. = FALSE)
  }

  for (req_type in unique_types) {
    if (!req_type %in% observed_types) {
      stop("No '", req_type, "' rows found in '", type_col, "' column.",
           call. = FALSE)
    }
  }

  # 5.) Check for duplicate marker entries ----

  dup_rows <- duplicated(marker_list[, required_cols])

  if (any(dup_rows)) {
    n_dup <- sum(dup_rows)

    dup_breakdown <- sort(table(marker_list[dup_rows, "label"]), decreasing = TRUE)

    dup_items <- stats::setNames(
      paste0(names(dup_breakdown), ": ", as.vector(dup_breakdown)),
      rep("*", length(dup_breakdown))
    )

    print_alert("Removed {n_dup} duplicate row{?s} from input:", type = "w")
    cli::cli_bullets(dup_items)

    marker_list <- marker_list[!dup_rows, ]
  }

  # 6.) Split by type ----

  cat_split  <- split(marker_list, marker_list[[type_col]])
  broad_cats <- unique(cat_split[["broad"]][[cat_col]])
  fine_cats  <- unique(cat_split[["fine"]][[cat_col]])

  missing_broad_cats <- broad_cats[!broad_cats %in% fine_cats]

  if (length(missing_broad_cats) > 0L) {
    cli::cli_abort(c(
      "Broad categories with no matching fine entries:",
      "x" = "{.val {missing_broad_cats}}",
      "i" = "Every broad category must also appear as a category in the fine rows."
    ))
  }

  other_fine_cats <- setdiff(fine_cats, broad_cats)

  if (length(other_fine_cats) > 0L) {
    reassign_idx <- cat_split[["fine"]][[cat_col]] %in% other_fine_cats
    cat_split[["fine"]][[cat_col]][reassign_idx] <- unnamed_broad_cat_label

    print_alert("{length(other_fine_cats)} fine categor{?y/ies} re-assigned to {.val {unnamed_broad_cat_label}}: {other_fine_cats}")
  }

  # 7.) Build marker list structure ----

  out <- list()

  out$broad <- split(cat_split[["broad"]][[marker_col]],
                     cat_split[["broad"]][[cat_col]])

  fine_split <- split(cat_split[["fine"]], cat_split[["fine"]][[cat_col]])

  fine_split <- lapply(fine_split, \(x) split(x[[marker_col]], x[[label_col]]))

  if (unnamed_broad_cat_label %in% names(fine_split)) {
    other_group <- fine_split[[unnamed_broad_cat_label]]
    fine_split[[unnamed_broad_cat_label]] <- NULL
    fine_split[[unnamed_broad_cat_label]] <- other_group
  }

  out$fine <- fine_split

  return(out)

}



#' Build broad marker configuration
#'
#' Converts a named list of broad-category marker vectors into a structured
#' configuration list suitable for priority-based cell-type assignment. Each
#' category is annotated with an expression threshold, minimum co-expression
#' count, and a numeric priority rank derived from \code{priority_order}.
#'
#' @section Priority order:
#' The \code{priority_order} vector determines how ties are resolved when a
#' cell or cluster matches multiple broad categories. Categories listed
#' earlier receive higher priority (lower rank number). At most one category
#' may be omitted from \code{priority_order}, in which case it is
#' automatically assigned the lowest priority. If two or more categories are
#' omitted they would share the same rank, making tie-breaking ambiguous -
#' this is treated as an error.
#'
#'
#' @param marker_list Named list of character vectors. Names are broad category
#'   labels (e.g. \code{"immune"}, \code{"vasculature"}), values are marker
#'   gene symbols.
#' @param priority_order Character vector defining the assignment priority.
#'   Categories listed first receive lower (higher-priority) rank values.
#'   All entries must correspond to names in \code{marker_list}. At most one
#'   category in \code{marker_list} may be absent from this vector (see
#'   Priority order section).
#' @param default_threshold Numeric scalar (>= 0). Default expression threshold
#'   applied to every category. Defaults to \code{0.1}.
#' @param default_coexp Positive integer scalar. Minimum number of co-expressed
#'   markers required for a category call. Defaults to \code{1}. Deprecated in this version,
#'   but retained for backward compatibility - will be removed in a future release.
#' @param per_category_overrides Optional named list of named lists, keyed by
#'   category name. Each inner list may contain \code{expr_threshold} and/or
#'   \code{coexp_min} to override the defaults for that category. Categories
#'   not listed use the defaults. Defaults to \code{NULL} (no overrides).
#'
#' @return A named list (one element per category, ordered by priority rank
#'   ascending) where each element is a list with components:
#'   \describe{
#'     \item{markers}{Character vector of marker gene symbols.}
#'     \item{expr_threshold}{Numeric expression threshold for the category.}
#'     \item{coexp_min}{Integer minimum co-expression count.}
#'     \item{priority}{Integer priority rank (1 = highest priority).}
#'   }
#'
#' @examples
#' \dontrun{
#'   markers <- list(
#'     immune      = c("CD45", "CD3", "CD8"),
#'     vasculature = c("PECAM1", "VWF"),
#'     stromal     = c("COL1A1", "VIM")
#'   )
#'
#'   # All categories ranked explicitly
#'   cfg <- build_broad_marker_config(
#'     marker_list    = markers,
#'     priority_order = c("vasculature", "immune", "stromal")
#'   )
#'
#'   # One category omitted - automatically gets lowest priority
#'   cfg <- build_broad_marker_config(
#'     marker_list    = markers,
#'     priority_order = c("vasculature", "immune")
#'   )
#'   # vasculature=1, immune=2, stromal=3
#'
#'   # Per-category overrides
#'   cfg <- build_broad_marker_config(
#'     marker_list    = markers,
#'     priority_order = c("vasculature", "immune", "stromal"),
#'     per_category_overrides = list(
#'       immune = list(coexp_min = 2, expr_threshold = 0.2)
#'     )
#'   )
#' }
#'
#' @export
build_broad_marker_config <- function(
    marker_list,
    priority_order,
    default_threshold = 0.1,
    default_coexp = 1,
    per_category_overrides = NULL
) {

  # -- Input validation ---

  if (!is.list(marker_list) || is.null(names(marker_list))){
    cli::cli_abort("{.arg marker_list} must be a named list.")
  }

  if (any(nchar(names(marker_list)) == 0L)) {
    cli::cli_abort("All elements of {.arg marker_list} must be named.")
  }

  if (anyDuplicated(names(marker_list))) {
    dups <- unique(names(marker_list)[duplicated(names(marker_list))])
    cli::cli_abort("Duplicate category names in {.arg marker_list}: {.val {dups}}")
  }

  checkmate::assert_character(priority_order, min.len = 1L, any.missing = FALSE)
  checkmate::assert_number(default_threshold, lower = 0)
  checkmate::assert_count(default_coexp, positive = TRUE)

  # -- Validate priority_order entries exist ---

  missing_prio <- setdiff(priority_order, names(marker_list))

  if (length(missing_prio) > 0L) {
    cli::cli_abort(
      "Priority categories not found in {.arg marker_list}: {.val {missing_prio}}"
    )
  }

  # -- Validate priority coverage ---

  unranked <- setdiff(names(marker_list), priority_order)

  if (length(unranked) > 0L) {

    n_ranked   <- length(priority_order)
    n_unranked <- length(unranked)

    if (n_ranked == 0L) {

      cli::cli_abort(c(
        "No categories have an assigned priority.",
        "x" = "All {.val {n_unranked}} categor{?y/ies} would receive the same default rank.",
        "i" = "Please supply a {.arg priority_order} covering at least the categories
               that should take precedence when a cell matches multiple groups.",
        "i" = "Unranked categories: {.val {unranked}}"
      ))

    }

    if (n_unranked > 1L) {

      cli::cli_abort(c(
        "Multiple categories missing from {.arg priority_order}.",
        "x" = "{.val {n_unranked}} categories would share the same lowest rank: {.val {unranked}}",
        "i" = "When multiple categories are unranked, ties cannot be resolved during
               broad assignment. Please include all categories in {.arg priority_order},
               or leave at most one category unranked."
      ))

    }

    print_alert(
      "Category {.val {unranked}} not in {.arg priority_order} - assigned lowest priority (rank {.val {n_ranked + 1L}})",
      type = "w"
    )
  }

  # -- Validate per-category overrides ---

  if (!is.null(per_category_overrides)) {

    if (!is.list(per_category_overrides) || is.null(names(per_category_overrides))) {
      cli::cli_abort("{.arg per_category_overrides} must be a named list.")
    }

    unknown_cats <- setdiff(names(per_category_overrides), names(marker_list))

    if (length(unknown_cats) > 0L) {
      print_alert(
        "Overrides for unknown categories will be ignored: {.val {unknown_cats}}", type = "w"
      )
    }

    allowed_fields <- c("expr_threshold", "coexp_min")

    for (cat in intersect(names(per_category_overrides), names(marker_list))) {

      bad_fields <- setdiff(names(per_category_overrides[[cat]]), allowed_fields)

      if (length(bad_fields) > 0L) {
        print_alert(
          "Unknown override fields for {.val {cat}} will be ignored: {.val {bad_fields}}",
          type = "w"
        )
      }

    }

  }

  # -- Validate marker vectors ---

  empty_cats <- names(marker_list)[vapply(marker_list, length, integer(1)) == 0L]

  if (length(empty_cats) > 0L) {
    print_alert("Categories with no markers: {.val {empty_cats}}", type = "w")
  }

  # -- Build config ---

  config_list <- lapply(names(marker_list), function(cat_name) {

    markers <- as.character(marker_list[[cat_name]])

    prio_rank <- match(cat_name, priority_order,
                       nomatch = length(priority_order) + 1L)

    cat_threshold <- default_threshold

    cat_coexp     <- default_coexp

    if (!is.null(per_category_overrides[[cat_name]])) {

      ovr <- per_category_overrides[[cat_name]]

      if (!is.null(ovr$expr_threshold)) cat_threshold <- ovr$expr_threshold

      if (!is.null(ovr$coexp_min))      cat_coexp     <- ovr$coexp_min

    }

    return(
      list(
            markers        =   markers,
            expr_threshold =   cat_threshold,
            coexp_min      =   cat_coexp,
            priority       =   prio_rank
            )
      )
  })

  names(config_list) <- names(marker_list)

  config_list <- config_list[order(vapply(config_list, `[[`, integer(1), "priority"))]

  return(config_list)

}




#' Extract distinct markers per cluster
#'
#' Converts a data frame of differential expression (DE) results into a named list,
#' where each element represents a cluster and contains a vector of unique gene symbols.
#'
#' @param de_results A data frame containing DE results (e.g., from Seurat's \code{FindAllMarkers}).
#' @param cluster_col Character. The name of the column containing cluster identities.
#'   Default is "cluster".
#' @param gene_col Character. The name of the column containing gene symbols.
#'   Default is "gene".
#'
#' @return A named list of character vectors. Names correspond to cluster IDs.
#' @keywords internal
distinct_cluster_markers <- function(
    de_results,
    cluster_col = "cluster",
    gene_col = "gene"
) {

  return(
    lapply(split(de_results[[gene_col]], de_results[[cluster_col]]), unique)
  )

}





#' Validate broad marker structure
#'
#' Checks if a broad marker list is named and contains the required structural fields.
#'
#' @param broad_marker_list A list of broad markers.
#' @param struct_fields A character vector of fields required inside each list element.
#' Defaults to `c("markers", "expr_threshold", "coexp_min", "priority")`.
#'
#' @return `TRUE` if the structure is valid, `FALSE` otherwise. Throws an error if elements are unnamed.
#' @keywords internal
.valid_broad_marker_struct <- function(
    broad_marker_list,
    struct_fields = c("markers", "expr_threshold", "coexp_min", "priority")
){

  if(length(names(broad_marker_list)) != length(broad_marker_list)){
    cli::cli_abort("{.arg broad_marker_list} elements are not named.")
  }

  return(all(sapply(broad_marker_list, \(x) all(struct_fields %in% names(x)))))
}




#' Check for missing features in a set
#'
#' Verifies that all features in `set1` are present in `set2`.
#'
#' @param set1 A character vector of query features (e.g., markers to look for).
#' @param set2 A character vector of target features (e.g., available genes in data).
#'
#' @return NULL. Throws an error if any features in `set1` are missing from `set2`.
#' @keywords internal
.check_missing_feats <- function(set1, set2){

  missing_genes <- setdiff(set1, set2)

  if (length(missing_genes) > 0) {
    cli::cli_abort(c(
      "x" = "The following features are not present in {.arg {set2}}:",
      "i" = "{.val {missing_genes}}",
      "!" = "Please check features/marker list to ensure they exist before re-running."
    ))
  }
}




#' Unpack marker list to Tibble
#'
#' Converts a named list of markers (either simple or nested) into a long-format tibble.
#'
#' @param marker_list A list of markers. Can be a simple named list of vectors or a nested list.
#' @param val_col Name of the output column for values (default: "gene").
#' @param ind_col Name of the output column for indices/names (default: "type").
#'
#' @return A tibble with columns defined by `val_col` and `ind_col`.
#' @importFrom utils stack
#' @importFrom magrittr `%>%`
#' @importFrom rlang `:=`
#' @importFrom dplyr .data
#' @keywords internal
.unpack_markers <- function(marker_list, val_col = "gene", ind_col = "type") {

  if (!is.list(marker_list)) cli::cli_abort("{.arg marker_list} is not a list ")

  if (length(marker_list) == 0) {
    return(tibble::tibble(!!ind_col := character(), !!val_col := character()))
  }

  is_nested <- is.list(marker_list[[1]])

  if (is_nested){
    out <- lapply(marker_list, stack) %>% dplyr::bind_rows()
  } else out <- stack(marker_list)

  out <- out %>%
    dplyr::rename(!!val_col := .data$values, !!ind_col := .data$ind) %>%
    tibble::as_tibble()

  return(out)
}




#' Define a complete feature-space
#'
#' Extracts a unique vector of all features present in both the broad and fine marker lists.
#' Handles both structured (complex) and simple broad marker lists.
#'
#' @param broad_marker_list A list of broad markers (simple or structured).
#' @param fine_marker_list A list of fine markers.
#'
#' @return A unique character vector of all features found in the inputs.
#' @keywords internal
.feature_space <- function(broad_marker_list, fine_marker_list){

  if(.valid_broad_marker_struct(broad_marker_list) ){
    all_broad_markers <- unlist(lapply(broad_marker_list, `[[`, "markers"), use.names = FALSE)
  } else {
    all_broad_markers <-  unlist(broad_marker_list, use.names = FALSE)
  }

  fine_markers <- unlist(fine_marker_list, use.names = FALSE)

  all_features <- unique(c(all_broad_markers, fine_markers))

  return(all_features)
}




#' Retrieve gene synonyms from Ensembl
#'
#' Connects to the Ensembl BioMart database to retrieve gene synonyms for a
#' provided list of gene names. It queries both the 'external_gene_name' and
#' 'external_synonym' filters to maximize match rates.
#'
#' @param lookup_genes Character vector. A list of gene names (or synonyms)
#'   to search for.
#'
#' @return A named list where:
#'   \itemize{
#'     \item \strong{Names} are the standard Ensembl external gene names.
#'     \item \strong{Values} are character vectors of unique synonyms for that gene.
#'   }
#'   Returns \code{NULL} if the connection to Ensembl fails or if no genes are found.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks if \code{biomaRt} is installed.
#'   \item safely connects to the \code{hsapiens_gene_ensembl} dataset.
#'   \item Queries Ensembl twice: once matching against gene names and once
#'     matching against existing synonyms.
#'   \item Cleans the results by converting empty strings \code{""} to \code{NA}
#'     and removing incomplete records.
#' }
#'
#' @examples
#' \dontrun{
#'   # Basic usage
#'   syns <- gene_synonyms(c("TP53", "BRCA1"))
#'
#'   # Access synonyms for TP53
#'   print(syns$TP53)
#' }
#' @export
#'
gene_synonyms <- function(lookup_genes) {

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    cli::cli_abort(c(
      "x" = "Package 'biomaRt' is not installed.",
      "i" = "Please install it and try again."
    ))
  }

  if (missing(lookup_genes)) {
    cli::cli_abort(c(
      "x" = "No lookup genes provided.",
      "i" = "Please provide a vector of gene names."
    ))
  }

  mart <- tryCatch({
    biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  }, error = function(e) {
    print_alert(glue::glue("Error connecting to biomaRt: {e$message}"), type = "w", face = "i")
    return(NULL)
  })

  if (is.null(mart)) return(NULL)

  hsap_version <- biomaRt::searchDatasets(mart, pattern = "hsapiens")$version[[1]]

  print_alert(
    glue::glue("Using Ensembl version: {hsap_version }"),
    type = "i", face = "i", color = "grey90"
    )

  res <- list()

  res$by_name <- biomaRt::getBM(
      attributes = c("external_gene_name", "external_synonym"),
      filters = c("external_gene_name"),
      values = lookup_genes,
      mart = mart
    )

    res$by_synonym <- biomaRt::getBM(
      attributes = c("external_gene_name", "external_synonym"),
      filters = c("external_synonym"),
      values = lookup_genes,
      mart = mart
    )

  res <- dplyr::bind_rows(res)

  if (nrow(res) == 0) {
    print_alert("Did not find any lookup genes in Ensembl", type = "w")
    return(NULL)
  }

  res <- res %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.character), ~dplyr::na_if(., ""))) %>%
    tidyr::drop_na()

  out <- split(res, as.factor(res$external_gene_name)) |>
    lapply(\(x) unique(x[["external_synonym"]]))

  return(out)
}
