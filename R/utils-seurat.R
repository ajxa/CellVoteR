#' Validate Seurat object input
#'
#' Checks if the provided object inherits from the "Seurat" class.
#'
#' @param sObj The object to validate.
#'
#' @return NULL. Throws a an error if the object is not a Seurat object.
#' @keywords internal
.valid_sObj_input <- function(sObj){
  if(!inherits(sObj, "Seurat")) cli::cli_abort("{.arg sObj} must be a Seurat object.")
}

#' Check for the presence of the Normalised Data Layer in a Seurat object.
#'
#' Verifies that the Seurat object contains a "data" layer (indicating normalisation has occurred).
#'
#' @param sObj A Seurat object.
#'
#' @return NULL. Throws an error if the "data" layer is missing.
#' @keywords internal
.is_norm_sObj <- function(sObj){

  if(!"data" %in% SeuratObject::Layers(sObj)){
    cli::cli_abort("The {.val data} slot is empty. Please run {.fn NormalizeData} before classification.")
  }
}

#' Filter Seurat object features
#'
#' Subsets a Seurat object to retain only a specific list of features.
#' Performs validation to ensure the features actually exist in the object.
#'
#' @param sObj A Seurat object.
#' @param filter_features A character vector of features (genes) to keep.
#'
#' @return A subsetted Seurat object containing only the common features.
#' @keywords internal
.filter_sObj_features <- function(sObj, filter_features){

  .valid_sObj_input(sObj)

  common_features <- intersect(rownames(sObj), filter_features)

  if(length(common_features) == 0){
    cli::cli_abort("No common features found between {.arg sObj} and {.arg filter_features}.")
  } else {
    sObj <- sObj[common_features, ]
    cli::cli_alert_info("Using reduced ({length(common_features)}) features for labelling.")
  }

  return(sObj)
}
