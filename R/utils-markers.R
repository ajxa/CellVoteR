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
