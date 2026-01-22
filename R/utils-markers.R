#' Clean and split broad and fine markers
#'
#' This function validates that all broad categories exist within the fine marker dataframe.
#' It reassigns unmatched fine categories to a generic label, splits both datasets into
#' hierarchical lists, and compiles a comprehensive vector of all unique markers found
#' across both datasets.
#'
#' @param fine_marker_df A dataframe containing specific cell markers and their labels.
#' @param broad_marker_df A dataframe containing higher-level categories to validate against.
#' @param fine_marker_col String. The column name for markers in `fine_marker_df`. Defaults to "marker".
#' @param fine_label_col String. The column name for the specific cell type label (e.g., "B cell") in `fine_marker_df`. Defaults to "fine_label".
#' @param fine_category_col String. The column name for the broader category (e.g., "immune") in `fine_marker_df`. Defaults to "fine_category".
#' @param broad_category_col String. The column name for categories in `broad_marker_df`. Defaults to "broad_category".
#' @param broad_marker_col String. The column name for markers in `broad_marker_df`. Defaults to "marker".
#' @param unnamed_broad_cat_label String. The label to assign to fine categories that do not match any broad category. Defaults to "other".
#'
#' @return A named list containing three elements:
#' \itemize{
#'   \item \strong{fine}: A nested list where top-level keys are `fine_category` and sub-keys are `fine_label`. Values are character vectors of markers.
#'   \item \strong{broad}: A named list where keys are `broad_category` and values are character vectors of markers.
#'   \item \strong{all}: A single character vector containing every unique marker found in both the fine and broad datasets.
#' }

#' @examples
#' # res <- clean_markers(fine_df, broad_df)
#' # head(res$all) # View all unique markers
#'
#' @export
#'
clean_markers <- function(
    fine_marker_df,
    broad_marker_df,
    fine_marker_col = "marker",
    fine_label_col = "fine_label",
    fine_category_col = "fine_category",
    broad_category_col = "broad_category",
    broad_marker_col = "marker",
    unnamed_broad_cat_label = "other"
){

  broad_cats <- unique(broad_marker_df[[broad_category_col]])
  fine_cats <- unique(fine_marker_df[[fine_category_col]])
  broad_cats_present <- broad_cats %in% unique(fine_marker_df[[fine_category_col]])

  if(any(!broad_cats_present)){
    missing_cats <- broad_cats[!broad_cats_present]
    cli::cli_abort(c(
      "Some broad categories are missing from the fine marker dataframe:",
      "x" = "{missing_cats}",
      "i" = "please ensure ALL broad categories are represented in the fine_marker_df"
    ))
  }

  fine_cats_to_other <- fine_cats[!fine_cats %in% broad_cats]

  if(length(fine_cats_to_other) > 0){
    print_alert(
      text = glue::glue(
        "{length(fine_cats_to_other)} categories re-assigned to '{unnamed_broad_cat_label}':\t",
        "{paste(fine_cats_to_other, collapse = ', ')}"
      )
    )

    fine_marker_df[[fine_category_col]][fine_marker_df[[fine_category_col]] %in% fine_cats_to_other] <- unnamed_broad_cat_label

  }

  out <- list()

  fine_split <- split(fine_marker_df, fine_marker_df[[fine_category_col]])

  out$fine <-  lapply(fine_split, \(x) split(x[[fine_marker_col]], x[[fine_label_col]]))

  out$broad <- split(broad_marker_df[[broad_marker_col]], broad_marker_df[[broad_category_col]])

  out$all <- unique(unlist(out, use.names = F))

  return(out)
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


#' Build broad marker configuration
#'
#' specific list structure required for annotation or scoring functions.
#' It assigns priorities, default expression thresholds, and co-expression requirements
#' to each cell type category.
#'
#' @param marker_list A named list where names are cell types (e.g., "T_cells", "B_cells")
#'   and values are character vectors of gene markers.
#' @param priority_order A character vector defining the hierarchy of cell types.
#'   Categories appearing earlier in this vector are assigned a lower numeric priority score
#'   (1 = highest priority).
#' @param default_threshold Numeric. The default expression detection threshold. Default is 0.1.
#' @param default_coexp Numeric/Integer. The minimum number of markers required to be
#'   co-expressed. Default is 1.
#'
#' @return A named list of configuration lists. Each element contains:
#' \itemize{
#'   \item \code{markers}: Vector of genes.
#'   \item \code{expr_threshold}: The threshold set.
#'   \item \code{coexp_min}: The co-expression minimum.
#'   \item \code{priority}: Integer rank based on \code{priority_order}.
#' }
#' @keywords internal
build_broad_marker_config <- function(
    marker_list,
    priority_order,
    default_threshold = 0.1,
    default_coexp = 1
) {

  if (!all(priority_order %in% names(marker_list))) {
    missing <- setdiff(priority_order, names(marker_list))
    cli::cli_abort("Priority categories {.val {missing}} not found in the input marker list.")
  }

  config_list <- lapply(names(marker_list), function(cat_name) {

    prio_rank <- match(cat_name, priority_order, nomatch = length(priority_order) + 1)

    list(
      markers        = as.vector(marker_list[[cat_name]], mode = "character"),
      expr_threshold = default_threshold,
      coexp_min      = default_coexp,
      priority       = prio_rank
    )
  })

  names(config_list) <- names(marker_list)
  return(config_list)
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
