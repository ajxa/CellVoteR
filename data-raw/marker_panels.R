marker_panel_file <- "../OneDrive - University of Leeds/analysis/CellVoteR/data/CellVoteR_marker_panels.xlsx"

marker_panels_names <- openxlsx::getSheetNames(marker_panel_file)

marker_panels <- list(
  GBM = list()
  )

marker_panels$GBM <-  lapply(marker_panels_names, \(x) openxlsx::read.xlsx(marker_panel_file, sheet = x)) |>
  setNames(marker_panels_names)

remove_dups <- function(panel_df, panel_name){

  duplicated_genes <- which(duplicated(panel_df))

  if(length(duplicated_genes) > 0) {

    cli::cli_alert_warning("Removed {length(duplicated_genes)} duplicated genes from {panel_name}")
    panel_df <- panel_df[-duplicated_genes, ]

  } else cli::cli_alert_success("No duplicated genes found in {panel_name}")

  return(panel_df)

}

marker_panels$GBM <- purrr::imap(marker_panels$GBM, remove_dups)

usethis::use_data(marker_panels, overwrite = TRUE)
