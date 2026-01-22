# Inputs ------------------------------------------------------------------

data_path <- "~/Desktop/OneDrive - University of Leeds/adhoc/Ensemble_Cell_Typing/data"

sc_data <- list(
  cosmx = list.files(data_path, pattern = "sObj\\.rds", full.names = TRUE),
  wang = list.files(file.path(data_path, "Wang_patient_seurats"), full.names = TRUE),
  nomura =  list.files(file.path(data_path, "Nomura_patient_seurats"), full.names = TRUE)
)

input_markers <- list()

input_markers$fine <- openxlsx::read.xlsx(
  xlsxFile = list.files(file.path(data_path), pattern = "input_markers", full.names = TRUE),
  sheet = 1
)

input_markers$broad <- openxlsx::read.xlsx(
  xlsxFile = list.files(file.path(data_path), pattern = "input_markers", full.names = TRUE),
  sheet = 2
)

# Load Data  -------------------------------------------------------------------
seu_data <- readRDS(sc_data$cosmx)[[1]]
seu_data <- assess_cell_quality(seu_data, remove_failed_cells = T)

input_markers$synonyms <- gene_synonyms(unique(input_markers$fine$marker))

ref <- clean_markers(
  fine_marker_df = input_markers$fine,
  broad_marker_df = input_markers$broad
  )

ref$broad_config <- build_broad_marker_config(
  marker_list = ref$broad, priority_order = c("vasculature", "immune")
  )
