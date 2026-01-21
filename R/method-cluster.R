#' Calculate dynamic clustering parameters
#'
#' Determines the optimal resolution and number of principal components (PCs)
#' based on the number of cells in the dataset.
#'
#' @details
#' This function scales parameters dynamically to avoid over-clustering small datasets
#' or under-clustering large ones:
#' \itemize{
#'   \item \strong{Resolution:} Scales linearly between \code{res_min} and \code{res_max}
#'   based on the square root of \code{n_cells}.
#'   \item \strong{PCs:} Scales logarithmically based on \code{n_cells}.
#' }
#' If \code{n_cells} is less than \code{mins_cluster_cells}, clustering is flagged to be skipped.
#'
#' @param n_cells Integer. The number of cells in the dataset/subset.
#' @param mins_cluster_cells Integer. Minimum number of cells required to perform clustering.
#'   If \code{n_cells} is below this, returns \code{skip = TRUE}.
#' @param res_min Numeric. Minimum resolution value.
#' @param res_max Numeric. Maximum resolution value.
#' @param res_saturation_n Integer. The cell count at which resolution reaches \code{res_max}.
#' @param npc_min Integer. Minimum number of PCs to use.
#' @param npc_max Integer. Maximum number of PCs to use.
#' @param npc_slope Numeric. Slope factor for the logarithmic calculation of PCs.
#' @param subcluster Logical. Controls the verbosity label.
#' @param verbose Logical. If \code{TRUE}, prints the calculated parameters to the console.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{npcs}: Integer. Number of PCs to use.
#'   \item \code{dims}: Integer vector. Sequence of dimensions (e.g., \code{1:30}).
#'   \item \code{res}: Numeric. The calculated resolution.
#'   \item \code{skip}: Logical. Whether to skip clustering (TRUE if cells < \code{mins_cluster_cells}).
#' }
#'
#' @export

set_cluster_params <- function(
    n_cells,
    mins_cluster_cells = 500,
    res_min = 0.6,
    res_max = 2.0,
    res_saturation_n = 100000,
    npc_min = 20,
    npc_max = 50,
    npc_slope = 4,
    subcluster = FALSE,
    verbose = TRUE
){

  if (n_cells < mins_cluster_cells) {

    return(list(npcs = 0, dims = 0, res = 0, skip = TRUE))

  } else {

    # Dynamic resolution:
    ratio <- sqrt(n_cells) / sqrt(res_saturation_n )
    calc_res <- res_min + (ratio * (res_max - res_min))

    final_res <- max(res_min, min(calc_res, res_max))
    final_res <- round(final_res, 2)

    # Dynamic PCs:
    # Formula: y = slope * log2(x) - intercept
    intercept <- round((npc_slope * log2(mins_cluster_cells)) - npc_min)
    calc_npcs <- round(npc_slope * log2(n_cells) - intercept)
    final_npcs <- max(npc_min, min(npc_max, calc_npcs))

    params <- list(
      npcs = final_npcs,
      dims = seq_len(final_npcs),
      res  = final_res,
      skip = FALSE
    )
  }

  msg_type <- if(subcluster) "sub-clustering" else "clustering"

  if (verbose) {
    cli::cli_h2("{msg_type} ({n_cells} cells)")
    cli::cli_ul(c(
      "PCA components: {.val {params$npcs}}",
      "Resolution: {.val {params$res}}"
    ))
  }

  return(params)
}


#' Run general, automated clustering on a Seurat object
#'
#' Performs a complete clustering workflow including normalization, scaling,
#' PCA, finding neighbours, and graph-based clustering. It automatically adapts
#' clustering parameters (resolution and PCs) based on the number of cells using
#' \code{\link{set_cluster_params}}.
#'
#' @details
#' \strong{Workflow:}
#' \enumerate{
#'   \item Checks and removes existing \code{scale.data} and variable features to ensure consistency.
#'   \item (Optional) Subsets the object to specific \code{filter_features}.
#'   \item Normalizes and Scales data (\code{Seurat::NormalizeData}, \code{Seurat::ScaleData}).
#'   \item Runs PCA (\code{Seurat::RunPCA}).
#'   \item Finds Neighbors and Clusters (\code{Seurat::FindNeighbors}, \code{Seurat::FindClusters}).
#' }
#'
#' \strong{Note on Existing Data:} This function will suppress warnings and remove existing
#' \code{scale.data} and variable features from the active assay before recalculating them.
#'
#' @param sObj A Seurat object.
#' @param algorithm Integer. The modularity optimization algorithm to use (passed to \code{Seurat::FindClusters}).
#'   Default is \code{4} (Leiden).
#' @param randseed Integer. Seed for random number generation to ensure reproducibility.
#' @param filter_features Character vector or NULL. If provided, the Seurat object is subset to these
#'   features (genes) before clustering. Useful for sub-clustering on specific gene sets.
#' @param subclustering Logical. Affects the logging output labels (toggles between "clustering" and "sub-clustering").
#' @param cluster_name Character. The name of the metadata column where cluster IDs will be stored.
#' @param set_new_ident Logical. If \code{TRUE}, updates the Seurat object's active identity (\code{Idents(sObj)})
#'   to the new clusters.
#' @param verbose Logical. If \code{TRUE}, prints progress bars and alerts.
#'
#' @return A Seurat object with updated metadata containing the new cluster IDs.
#'   If the number of cells is below the threshold defined in \code{set_cluster_params},
#'   all cells are assigned to a single cluster named "1".
#'
#' @importFrom Seurat NormalizeData ScaleData FindVariableFeatures RunPCA FindNeighbors FindClusters
#' @importFrom SeuratObject Layers Idents `Idents<-` VariableFeatures `VariableFeatures<-`
#'
#' @export
cluster_cells <- function(
    sObj,
    algorithm = 4,
    randseed = 123,
    filter_features = NULL,
    subclustering = FALSE,
    cluster_name = "cellvoter_cluster",
    set_new_ident = TRUE,
    verbose = TRUE
){

  cli::cat_line()

  cluster_params <- set_cluster_params(
    n_cells = ncol(sObj),
    subcluster = subclustering,
    verbose = verbose
  )

  if(cluster_params$skip) {
    cli::cli_alert_info("Assigned all cells to a single cluster identity.")

    sObj[[cluster_name]] <- "1"

    if(set_new_ident) Seurat::Idents(sObj) <- sObj[[cluster_name, drop = TRUE]]

    return(sObj)
  }

  # detect seurat object assays ----
  active_assays <- SeuratObject::Layers(sObj[[sObj@active.assay]])

  if("scale.data" %in% active_assays){
    has_scaled_data <- nrow(sObj[[sObj@active.assay]]["scale.data"]) > 0
  }else has_scaled_data <- FALSE

  has_var_features <- length(VariableFeatures(sObj)) > 0

  if(has_scaled_data){

    suppressWarnings({
      sObj[[sObj@active.assay]]$scale.data = NULL
    })

    cli::cli_alert_warning("Removed existing scaled data for consistency - will be re-calculated.")
  }

  if(has_var_features){

    VariableFeatures(sObj[[sObj@active.assay]]) <- NULL

    cli::cli_alert_warning("Removed existing variable features for consistency - will be re-calculated.")
  }


  # (optional) filter the feature-space ----
  if(!is.null(filter_features)){

    common_features <- intersect(rownames(sObj), filter_features)

    sObj <- sObj[common_features, ]

    cli::cli_alert_info("Using reduced ({length(common_features)}) features for clustering.")

  } else {
    cli::cli_alert_info("Using all {nrow(sObj)} features for clustering.")
  }

  # normalise and scale ----
  cli::cli_progress_step("Normalizing, scaling, and finding variable features...", spinner = TRUE)
  sObj <- Seurat::NormalizeData(sObj, verbose = verbose)
  sObj <- Seurat::ScaleData(object = sObj, features = rownames(sObj), verbose = verbose)
  sObj <- Seurat::FindVariableFeatures(sObj, verbose = verbose)

  # dimension reduction --------------------------------------------------------
  cli::cli_progress_step("Running PCA...",  spinner = TRUE)
  .hush({
    sObj <- Seurat::RunPCA(
      object = sObj,
      npcs = cluster_params$npcs,
      features = VariableFeatures(sObj)
    )
  })

  # clustering -----------------------------------------------------------------
  cli::cli_progress_step("Finding neighbors and clusters...",  spinner = TRUE)

  suppressMessages({

    sObj <- Seurat::FindNeighbors(
      object = sObj, dims = cluster_params$dims,  reduction = "pca", verbose = verbose)

    sObj <- Seurat::FindClusters(
      object = sObj,
      resolution = cluster_params$res,
      random.seed = randseed,
      algorithm = algorithm,
      verbose = verbose,
      cluster.name = cluster_name
    )

  })

  # clean-up and out ------------------------------------------------------------
  if(set_new_ident){
    cli::cli_alert_warning("Set active Idents to {.emph {cluster_name}}")
    Idents(sObj) <- sObj@meta.data[[cluster_name]]
  }

  return(sObj)
}



