data_path <- "~/Desktop/OneDrive - University of Leeds/adhoc/Ensemble_Cell_Typing/data"

sc_data <- list(
  cosmx = list.files(data_path, pattern = "sObj\\.rds", full.names = TRUE),
  wang = list.files(file.path(data_path, "Wang_patient_seurats"), full.names = TRUE),
  nomura =  list.files(file.path(data_path, "Nomura_patient_seurats"), full.names = TRUE)
)


ref <- readRDS(file.path(data_path, "processed_markers.rds"))


seu_data <- readRDS(sc_data$cosmx)[[1]]
seu_data <- assess_cell_quality(seu_data, remove_failed_cells = T)



