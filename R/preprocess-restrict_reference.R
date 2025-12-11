#' Restrict object to reference genes
#'
#' @param object Seurat object.
#' @param reference_genes Character vector of genes to keep.
#' @return Subsetted Seurat object.
#' @noRd
restrict_to_reference <- function(object, reference_genes) {
  valid_genes <- intersect(rownames(object), reference_genes)
  if (length(valid_genes) < 10) stop("Filtering to reference left <10 genes.")
  object <- object[valid_genes, ]
  return(object)
}
