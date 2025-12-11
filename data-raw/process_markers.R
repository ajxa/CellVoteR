# --- Triage Groups (Broad Classification) -------------------------------------
# Used by Methods 1-4 to split data before labelling
triage_groups <- list(
  Immune = c("PTPRC"),
  Endothelial = c("CDH5", "VWF", "CLDN5", "PECAM1")
)

# --- SAA GBMDeconvoluetR/Moreno/Neftel (The 22 types) -------------------------
path_to_SAA_markers <- "~/Desktop/OneDrive - University of Leeds/adhoc/Ensemble_Cell_Typing/data/markers.rds"

saa_markers <- readRDS(path_to_SAA_markers)

required_cols <- c("gene", "cell_type", "cell_category")

colnames(saa_markers$markers) <- required_cols

idhwt_gbm_markers <- saa_markers$markers

if (!all(required_cols %in% colnames(idhwt_gbm_markers))) {
  stop("idhwt_gbm_markers must have columns: gene, cell_type, category")
}

# --- Master List --------------------------------------------------------------
cellvoter_data <- list(
  marker_panels = list(
    idhwt_gbm_markers = idhwt_gbm_markers
    ),
  cell_groups = triage_groups
  )


usethis::use_data(cellvoter_data, internal = FALSE, overwrite = TRUE)

