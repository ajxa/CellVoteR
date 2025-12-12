#' Extract count matrix from various inputs
#'
#' @param object A Seurat object, SingleCellExperiment, or matrix.
#'
#' @return A sparse count matrix (dgCMatrix).
#' @noRd
#' @importFrom SummarizedExperiment assays assay assayNames
#' @importFrom SingleCellExperiment counts
#' @importFrom methods as is
extract_counts <- function(object) {

  counts <- NULL

  # 1. Handle Seurat Objects
  if (inherits(object, "Seurat")) {
    counts <- Seurat::GetAssayData(object, layer = "counts")

  }

  # 2. Handle SingleCellExperiment
  else if (inherits(object, "SingleCellExperiment")) {
    # Check if "counts" is a valid assay name
    # We use SummarizedExperiment generic to be safe
    nms <- SummarizedExperiment::assayNames(object)
    if ("counts" %in% nms) {
      counts <- SummarizedExperiment::assay(object, "counts")
    } else if (length(nms) > 0) {
      # Fallback to first assay if "counts" isn't named explicitly
      warning("No 'counts' assay found. Using the first assay: ", nms[1])
      counts <- SummarizedExperiment::assay(object, 1)
    } else {
      stop("Input is a SingleCellExperiment but no assays found.")
    }
  }

  # 3. Handle Matrix / Sparse Matrix
  else if (inherits(object, c("Matrix", "dgCMatrix", "matrix"))) {
    counts <- object
  }

  else {
    stop("Input format not supported. Please provide a Seurat object, SingleCellExperiment, or matrix.")
  }

  # Ensure it is strictly a sparse matrix
  if (!inherits(counts, "dgCMatrix")) {
    counts <- methods::as(counts, "dgCMatrix")
  }

  return(counts)
}

#' Ensure input is a Seurat object
#'
#' Converts SingleCellExperiment or Matrix objects to Seurat.
#'
#' @param object Input object (Seurat, SingleCellExperiment, or Matrix).
#'
#' @return A Seurat object.
#' @noRd
#' @importFrom Seurat CreateSeuratObject
ensure_seurat <- function(object) {

  if (inherits(object, "Seurat")) {
    return(object)
  }

  # Use the helper above to get the raw data safely
  counts <- extract_counts(object)

  # Create Seurat Object
  # min.cells/features = 0 ensures we don't accidentally drop data the user wanted to keep
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)

  return(seurat_obj)
}
