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




#' Expand a marker reference with gene aliases
#'
#' Augments a reference marker data frame by adding rows for gene aliases.
#' It effectively duplicates the metadata of an original marker and assigns it
#' to its alias, allowing the downstream scoring system to recognize both names.
#'
#' @param ref_df A data frame containing the reference markers and associated metadata.
#' @param alias_list A named list where names are the original markers (matching \code{ref_df})
#'   and values are vectors of alias names.
#' @param ref_marker_col Character. The column name in \code{ref_df} that holds the
#'   primary gene symbols. Default is "marker".
#'
#' @examples
#' \dontrun{
#'   # Define aliases: CD8A can also be detected as CD8
#'   aliases <- list(CD8A = "CD8")
#'
#'   # Expand the reference table
#'   new_ref <- expand_markers_with_aliases(reference_data, aliases)
#' }
#' @return A data frame with additional rows for the aliases. Duplicates are removed.
#' @importFrom dplyr .data
#' @keywords internal
expand_markers_with_aliases <- function(
    ref_df,
    alias_list,
    ref_marker_col = "marker"
) {


  alias_map <- stack(alias_list) %>%
    dplyr::rename(alias = .data$values, original_marker = .data$ind) %>%
    dplyr::mutate(original_marker = as.character(.data$original_marker),
                  alias = as.character(.data$alias))

  rem_cols <- c("original_marker", "alias")

  alias_rows <- alias_map %>%
    dplyr::inner_join(ref_df, by = c("original_marker" = ref_marker_col)) %>%
    dplyr::mutate(marker = .data$alias) %>%
    dplyr::select(-dplyr::all_of(rem_cols))

  updated_df <- dplyr::bind_rows(ref_df, alias_rows) %>%
    dplyr::distinct()

  n_added <- nrow(updated_df) - nrow(ref_df)
  cli::cli_alert_success("Added {n_added} alias-based markers to the reference set.")

  return(updated_df)
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




#' Expand gene symbols with synonyms from biomaRt
#'
#' Connects to Ensembl to find synonyms for a list of genes. If `valid_genes` is provided,
#' it filters the results to match genes present in your dataset and supports regex patterns.
#'
#' @param input_features A character vector of gene symbols OR a single regex pattern string.
#' @param valid_genes (Optional) A character vector of all gene names in your dataset.
#'   If provided, used to filter results and allow regex matching. Default NULL.
#' @param organism Character. "hsapiens", "mmusculus", etc.
#'
#' @return A character vector of unique gene symbols found in `valid_genes`.
#' @export
feature_synonyms <- function(input_features, valid_genes = NULL, organism = "hsapiens") {

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required for synonym lookup. Please install it.")
  }

  # 1. Determine if input is a regex or a list of genes
  if (!is.null(valid_genes) && length(input_features) == 1 && grepl("[^a-zA-Z0-9]", input_features)) {
    search_genes <- grep(input_features, valid_genes, value = TRUE, ignore.case = TRUE)
  } else {
    search_genes <- input_features
  }

  if (length(search_genes) == 0) return(character(0))

  # 2. Setup BiomaRt
  dataset_name <- switch(tolower(organism),
                         "human" = "hsapiens_gene_ensembl",
                         "hsapiens" = "hsapiens_gene_ensembl",
                         "mouse" = "mmusculus_gene_ensembl",
                         "mmusculus" = "mmusculus_gene_ensembl",
                         stop("Organism not supported automatically."))

  mart <- tryCatch({
    biomaRt::useMart("ensembl", dataset = dataset_name)
  }, error = function(e) {
    warning("Could not connect to biomaRt: ", e$message)
    return(NULL)
  })

  if (is.null(mart)) return(search_genes)

  # 3. Query
  results <- biomaRt::getBM(
    attributes = c("external_gene_name", "external_synonym"),
    filters = "external_gene_name",
    values = search_genes,
    mart = mart
  )

  all_variants <- unique(c(results$external_gene_name, results$external_synonym))

  # 4. Filter synonyms against what is actually in the data
  if (!is.null(valid_genes)) {
    return(intersect(all_variants, valid_genes))
  }

  return(all_variants)
}
