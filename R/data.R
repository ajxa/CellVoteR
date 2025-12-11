#' CellVoteR Internal Data Reference
#'
#' The primary internal data structure containing marker panels and cell grouping definitions.
#'
#' @format A list with two main components:
#' \describe{
#'   \item{marker_panels}{A list of marker dataframes for specific tissues/contexts.
#'     \describe{
#'       \item{idhwt_gbm_markers}{The default panel for IDH-wildtype GBM. Contains 22 cell types
#'       across 'Malignant', 'Immune', and 'Normal Brain/Vasculature' categories.
#'       Derived from Neftel et al. (2019), Ruiz-Moreno et al. (2022), and Ajaib et al. (2022).}
#'     }
#'   }
#'   \item{cell_groups}{A list of markers used for broad cell type triage.
#'     \describe{
#'       \item{Immune}{Marker: PTPRC.}
#'       \item{Endothelial}{Markers: CDH5, VWF, CLDN5, PECAM1.}
#'     }
#'   }
#' }
"cellvoter_data"
