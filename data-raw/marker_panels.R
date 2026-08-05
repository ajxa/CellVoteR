marker_panel_file <- "../OneDrive - University of Leeds/adhoc/Ensemble_Cell_Typing/data/CellVoteR_marker_panels.xlsx"


marker_panels_names <- openxlsx::getSheetNames(marker_panel_file)

marker_panels <- list(
  GBM = list()
  )

marker_panels$GBM <-  lapply(marker_panels_names, \(x) openxlsx::read.xlsx(marker_panel_file, sheet = x)) |>
  setNames(marker_panels_names)


usethis::use_data(marker_panels, overwrite = TRUE)
