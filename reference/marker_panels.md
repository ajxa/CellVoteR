# Cell Type Marker Panels for Single-Cell Analysis

A nested list containing curated marker gene panels for cell-type and
cell-state annotation across different disease contexts. Currently only
includes IDHwt glioblastoma (GBM) marker sets under the `GBM` list
entry.

## Usage

``` r
data(marker_panels)
```

## Format

A named `list` of lists. The `GBM` sublist contains three `data.frame`
objects:

- GBM:

  A named list of marker panels for glioblastoma analysis:

  gbmap_neftel_full

  :   A data frame containing comprehensive marker genes across
      non-malignant and malignant glioblastoma cell populations.
      Non-malignant markers (normal brain, vasculature, and immune
      cells) originate from Ruiz-Moreno et al. (2025). Malignant cell
      state markers are derived from Neftel et al. (2019) and include
      all six fine-grained states: MES1, MES2, NPC1, NPC2, OPC, and AC.
      Includes rarer immune subtypes (monocytes, mast cells, dendritic
      cells) and normal brain progenitors (e.g., normal OPCs, radial
      glia).

  gbmap_neftel_reduced

  :   A streamlined version of `gbmap_neftel_full` optimised for smaller
      sample sizes or datasets where rare populations are absent. The
      Neftel et al. malignant states are collapsed into the four core
      meta-states (MES, NPC, OPC, AC). Non-malignant and rare cell types
      (e.g., normal OPCs) are excluded or simplified to minimise
      false-positive annotation caused by similarity between normal and
      neoplastic counterparts.

  nomura

  :   A data frame containing cell-type marker panels derived from the
      multi-layered transcriptional architecture of glioblastoma
      ecosystems defined by Nomura et al. (2025).

## Source

- **GBmap (Non-malignant, immune, and vascular states):** Ruiz-Moreno,
  C., Salas, S. M., Samuelsson, E., et al. (2025). Charting the
  Single-Cell and Spatial Landscape of IDH-Wild-Type Glioblastoma with
  GBmap. *Neuro-Oncology*, 27(9), 2281–2295.
  [doi:10.1093/neuonc/noaf113](https://doi.org/10.1093/neuonc/noaf113)

- **Neftel States (Malignant cellular states):** Neftel, C., Laffy, J.,
  Filbin, M. G., et al. (2019). An Integrative Model of Cellular States,
  Plasticity, and Genetics for Glioblastoma. *Cell*, 178(4),
  835–849.e21.
  [doi:10.1016/j.cell.2019.06.024](https://doi.org/10.1016/j.cell.2019.06.024)

- **Nomura Panel:** Nomura, M., Spitzer, A., Johnson, K. C., et al.
  (2025). The Multi-layered Transcriptional Architecture of Glioblastoma
  Ecosystems. *Nature Genetics*, 57(5), 1155–1167.
  [doi:10.1038/s41588-025-02167-5](https://doi.org/10.1038/s41588-025-02167-5)

## Details

Marker panels in `marker_panels$GBM$gbmap_neftel_full` and
`marker_panels$GBM$gbmap_neftel_reduced` integrate non-malignant,
immune, and vascular reference annotations from the GBmap atlas with
malignant state definitions from Neftel et al. (2019).
