# ORIGINAL LOGIC GBMDECONVOLUTER MARKERS CELL TYPE PROPORTIONS -----------------
seu_obj_6k <- readRDS("../OneDrive - University of Leeds/adhoc/Ensemble_Cell_Typing/data/GBMDeconvoluteR_Ensemble.Rds")
markers <- readRDS("../OneDrive - University of Leeds/adhoc/Ensemble_Cell_Typing/data/markers.rds")

cells_names = unlist(markers$labels$main, use.names = FALSE)
cells_names = cells_names[!cells_names %in% c("PN", "MES3")]
cells_names = c(cells_names, "Unknown")

cells = prop.table(table(seu_obj_6k$Resolved_Label))
cells = cells[cells_names]


gbmdeconv_6k_props = data.frame(
  cell_type = names(cells),
  prop = round(as.numeric(cells), digits = 3)
  )


sum(gbmdeconv_6k_props$prop)

openxlsx::write.xlsx(
  x = gbmdeconv_6k_props,
  file =
  "../OneDrive - University of Leeds/adhoc/Ensemble_Cell_Typing/outputs/CosMx_6k_GBMDeconvoluetR_props.xlsx",
  colNames = TRUE, rowNames = FALSE)

rm(cells, cells_names, gbmdeconv_6k_props, markers)

# 6K COSMX DATA RUN ------------------------------------------------------------
# renv::install("ajxa/CellVoteR", rebuild = TRUE)
seu_obj_6k <- readRDS("../OneDrive - University of Leeds/adhoc/Ensemble_Cell_Typing/data/GBMDeconvoluteR_Ensemble.Rds")

seu_obj_6k <- Seurat::SplitObject(seu_obj_6k, split.by = "Patient_P_R")
seu_obj_6k <- seu_obj_6k[1:5]

devtools::load_all()


# --- 2. Prepare a Clean Input Object ---
counts <- Seurat::GetAssayData(seu_obj_6k[[1]], layer = "counts")
object <- ensure_seurat(counts)

valid_markers <- process_marker_input()
triage_markers <- process_triage_input()

methods_to_use <- 1:4; breaker_methods <- 5:6



default_assay <- Seurat::DefaultAssay(object)
assay_layers <- SeuratObject::Layers(object, assay = default_assay)

if (!"data" %in% assay_layers) {
  if (verbose) message("Normalizing data (required for clash resolution)...")
  object <- Seurat::NormalizeData(object, verbose = FALSE)
}

vote_matrix <- matrix(NA, nrow = ncol(object), ncol = length(methods_to_use))
colnames(vote_matrix) <- paste0("Method_", methods_to_use)


# Pass markers and triage down
res <- run_cellvoter(object,
                     method = 1,
                     markers = valid_markers,
                     triage_markers = triage_markers,
                     verbose = TRUE
                     )


foo <- run_broad_labelling(object, broad_markers = triage_markers)


run_broad_labelling <- function(object, broad_markers = NULL) {

  broad_markers <- process_triage_input(broad_markers)

  if (!"pca" %in% names(object@reductions)) {
    object <- Seurat::NormalizeData(object, verbose=FALSE)
    object <- Seurat::ScaleData(object, rownames(object), verbose=FALSE)
    object <- Seurat::FindVariableFeatures(object, verbose=FALSE)
    object <- Seurat::RunPCA(object, npcs = 40,  verbose=FALSE)
    object <- Seurat::FindNeighbors(object, dims = 1:40, verbose=FALSE)
    object <- Seurat::FindClusters(object, resolution = 2, algorithm = 4, verbose=FALSE)
  }

  # need to suppress the warning message
  avg_exp <- Seurat::AverageExpression(object, features = unlist(triage_markers), assays = "RNA")$RNA

  cluster_labels <- rep("Other", ncol(avg_exp))
  names(cluster_labels) <- colnames(avg_exp)

  for (cluster in colnames(avg_exp)) {
    scores <- numeric()
    for (type in names(triage_markers)) {
      genes <- triage_markers[[type]]
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




# --- 3. Run CellVoteR (The New Way) ---
# We set a seed because Seurat's RunPCA/FindClusters has slight randomness
set.seed(123)


# use_ensemble = TRUE activates Methods 1-4 (Triage) + 5-6 (Global)
# This matches your full original pipeline logic.
test_obj <- run_ensemble(test_obj, use_ensemble = TRUE, verbose = TRUE)
