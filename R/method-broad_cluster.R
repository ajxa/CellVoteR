#' Perform initial broad cell type classification via Clustering
#'
#' @param object Seurat object.
#' @param broad_markers Named list of markers.
#'
#' @return Seurat object with 'Broad_Type' metadata.
#' @export
run_broad_labelling <- function(object,
                                broad_markers = list(Immune = "PTPRC", Endothelial = c("CDH5", "VWF"))) {

  if (!"pca" %in% names(object@reductions)) {
    object <- Seurat::NormalizeData(object, verbose=FALSE)
    object <- Seurat::FindVariableFeatures(object, verbose=FALSE)
    object <- Seurat::ScaleData(object, verbose=FALSE)
    object <- Seurat::RunPCA(object, verbose=FALSE)
    object <- Seurat::FindNeighbors(object, dims = 1:30, verbose=FALSE)
    object <- Seurat::FindClusters(object, resolution = 1.0, verbose=FALSE)
  }

  avg_exp <- Seurat::AverageExpression(object, features = unlist(broad_markers), assays = "RNA")$RNA

  cluster_labels <- rep("Other", ncol(avg_exp))
  names(cluster_labels) <- colnames(avg_exp)

  for (cluster in colnames(avg_exp)) {
    scores <- numeric()
    for (type in names(broad_markers)) {
      genes <- broad_markers[[type]]
      valid_genes <- intersect(genes, rownames(avg_exp))
      if(length(valid_genes) > 0) {
        scores[type] <- mean(avg_exp[valid_genes, cluster])
      } else {
        scores[type] <- 0
      }
    }

    if (max(scores) > 0.1) {
      cluster_labels[cluster] <- names(which.max(scores))
    }
  }

  object$Broad_Type <- cluster_labels[as.character(Seurat::Idents(object))]
  return(object)
}
