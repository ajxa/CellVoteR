#' Lean Sub-clustering (Optimized)
#' @noRd
process_subcluster <- function(object, min_cells = 10) {
  n_cells <- ncol(object)
  if (n_cells < min_cells) return(NULL)

  if (n_cells < 30) {
    npcs <- 2; res <- 0.1
  } else if (n_cells < 100) {
    npcs <- 10; res <- 2.0
  } else {
    npcs <- 30; res <- 2.0
  }

  tryCatch({
    object <- Seurat::NormalizeData(object, verbose = FALSE)
    object <- Seurat::FindVariableFeatures(object, selection.method = "vst", verbose = FALSE)
    object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE)
    suppressWarnings({
      object <- Seurat::RunPCA(object, npcs = npcs, verbose = FALSE)
    })
    object <- Seurat::FindNeighbors(object, dims = 1:npcs, verbose = FALSE)
    object <- Seurat::FindClusters(object, resolution = res, verbose = FALSE)

    markers <- Seurat::FindAllMarkers(object, only.pos = TRUE,
                                      min.pct = 0.25, logfc.threshold = 0.25,
                                      verbose = FALSE)
    return(list(object = object, markers = markers))
  }, error = function(e) { return(NULL) })
}
