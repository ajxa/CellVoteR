# Adds GBM_RNA_reference_panel_v1 to the built-in marker panels.
#
# Run after data-raw/marker_panels.R (which rebuilds the original three panels
# from the shared xlsx). This script appends the reference panel from the CSV
# kept alongside it, without touching the other panels.

load("data/marker_panels.rda")

panel <- utils::read.csv("data-raw/GBM_RNA_reference_panel_v1.csv",
                         stringsAsFactors = FALSE)

# Match the layout of the other built-in panels: broad rows carry NA labels
panel$label[!nzchar(panel$label)] <- NA_character_

stopifnot(
  identical(colnames(panel), c("type", "category", "label", "marker")),
  nrow(panel) == 1031L,
  !anyDuplicated(panel)
)

marker_panels$GBM$GBM_RNA_reference_panel_v1 <- panel

usethis::use_data(marker_panels, overwrite = TRUE)
