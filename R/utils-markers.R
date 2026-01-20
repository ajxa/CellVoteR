#' Process and Standardize User Marker Input
#'
#' Converts various input formats (Dataframe, List, CSV, Excel) into the
#' standardized 3-column format required by CellVoteR.
#'
#' @param input Input data. Can be:
#'   \itemize{
#'     \item \code{NULL}: Returns the internal default reference (idhwt_gbm_markers).
#'     \item \code{data.frame}: Must contain gene and cell type info.
#'     \item \code{list}: Named list (e.g. \code{list(T_cell = "CD3D")}).
#'     \item \code{character}: File path to .csv or .xlsx file.
#'   }
#' @param panel_name String. If \code{input} is NULL, selects this named panel from internal data.
#'                   Default: "idhwt_gbm_markers".
#' @param gene_col Name of the column containing gene symbols (default guesses: "gene", "symbol", "feature").
#' @param type_col Name of the column containing cell types (default guesses: "cell_type", "type", "cluster").
#' @param category_col Name of the column containing broad categories (default guesses: "category", "class", "broad_type").
#'
#' @return A data.frame with columns: \code{gene}, \code{cell_type}, \code{category}.
#' @export
process_marker_input <- function(input = NULL,
                                 panel_name = "idhwt_gbm_markers",
                                 gene_col = NULL,
                                 type_col = NULL,
                                 category_col = NULL) {

  # 1. Handle NULL (Default Internal Data) -------------------------------------
  if (is.null(input)) {
    if (!exists("cellvoter_data")) {
      stop("Internal data 'cellvoter_data' is missing. Please reinstall the package.")
    }

    # Check if the requested panel exists
    available_panels <- names(cellvoter_data$marker_panels)
    if (!panel_name %in% available_panels) {
      stop("Internal panel '", panel_name, "' not found. Available panels: ",
           paste(available_panels, collapse = ", "))
    }

    # Return the selected panel
    return(cellvoter_data$marker_panels[[panel_name]])
  }

  # 2. Handle File Paths (CSV / Excel) -----------------------------------------
  if (is.character(input) && length(input) == 1) {
    if (!file.exists(input)) stop("File not found: ", input)

    ext <- tolower(tools::file_ext(input))
    if (ext == "csv") {
      input <- utils::read.csv(input, stringsAsFactors = FALSE)
    } else if (ext %in% c("xls", "xlsx")) {
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        stop("Package 'openxslx' is required to read Excel files. Please install it")
      }
      input <- as.data.frame(openxlsx::read.xlsx(input))
    } else {
      stop("Unsupported file extension: .", ext)
    }
  }

  # 3. Handle Named List -------------------------------------------------------
  # Convert list(T_cell = "CD3D") -> DF(gene="CD3D", cell_type="T_cell")
  if (is.list(input) && !is.data.frame(input)) {
    if (is.null(names(input))) stop("Marker list must be named (e.g. list(TypeA = 'Gene1')).")

    df_list <- lapply(names(input), function(n) {
      data.frame(
        gene = input[[n]],
        cell_type = n,
        category = n, # Default category to the cell type itself if not specified
        stringsAsFactors = FALSE
      )
    })
    # Return immediately as we know the structure is correct
    return(do.call(rbind, df_list))
  }

  # 4. Handle Dataframe Standardization ----------------------------------------
  if (is.data.frame(input)) {
    cols <- colnames(input)

    # Helper to find column name if not provided
    find_col <- function(user_choice, candidates, mandatory = TRUE) {
      # If user specified a name, check it exists
      if (!is.null(user_choice)) {
        if (!user_choice %in% cols) stop("Column '", user_choice, "' not found in input.")
        return(user_choice)
      }
      # Otherwise search for candidates
      match <- grep(paste0("^", candidates, "$", collapse="|"), cols, ignore.case=TRUE, value=TRUE)
      if (length(match) > 0) return(match[1])

      if (mandatory) stop("Could not find a gene/type column. Please specify 'gene_col' or 'type_col'.")
      return(NULL)
    }

    # Identify Columns
    final_gene <- find_col(gene_col, c("gene", "genes", "symbol", "feature", "marker"))
    final_type <- find_col(type_col, c("cell_type", "type", "cluster", "annotation", "label", "broad_cell_type"))
    final_cat  <- find_col(category_col, c("category", "class", "broad_type", "group"), mandatory = FALSE)

    # Construct Standardized DF
    clean_df <- data.frame(
      gene = input[[final_gene]],
      cell_type = input[[final_type]],
      stringsAsFactors = FALSE
    )

    # Handle Category
    if (!is.null(final_cat)) {
      clean_df$category <- input[[final_cat]]
    } else {
      # Default: If no category provided, use the cell_type
      # This ensures downstream functions always have a 'category' column to filter on
      clean_df$category <- clean_df$cell_type
    }

    return(clean_df)
  }

  stop("Invalid input format. Must be Dataframe, Named List, or File Path.")
}


#' Process Triage Markers
#'
#' Helper to resolve triage markers from user input, internal defaults, or fallback.
#'
#' @param triage_markers Named list or NULL.
#' @return A named list of markers.
#' @export
process_triage_input <- function(triage_markers = NULL) {

  # 1. User supplied
  if (!is.null(triage_markers)) {
    if (!is.list(triage_markers)) stop("triage_markers must be a named list.")
    return(triage_markers)
  }

  # 2. Try Internal Data
  if (exists("cellvoter_data")) {
    return(cellvoter_data$cell_groups)
  }

  # 3. Fallback (Safety net)
  warning("Internal 'cellvoter_data' missing. Using hardcoded fallback for triage.")
  return(list(Immune = "PTPRC", Endothelial = c("CDH5", "VWF")))
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

#' Check for missing features in a Set
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
    dplyr::rename(!!val_col := values, !!ind_col := ind) %>%
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
