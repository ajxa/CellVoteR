# CellVoteR <img src="man/figures/logo.png" align="right" height="130" alt="" />

CellVoteR is an ensemble-based pipeline for robust cell type annotation in single-cell RNA-seq data.

Rather than relying on a single best-match score, 
CellVoteR combines complementary annotation strategies across full and 
marker-reduced feature spaces, then resolves agreement and disagreement 
through a configurable consensus voting step.


## Workflow

<img src="man/figures/cellvoter-workflow.svg"
     alt="CellVoteR workflow"
     width="100%">

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("ajxa/CellVoteR")
```

## Quick Start

```r
markers <- load_markers(file_path = "path/to/input_markers.xlsx")

markers$broad <- build_broad_marker_config(
  marker_list       = markers$broad,
  priority_order    = c("vasculature", "immune"),
  default_threshold = 0.25
)

sce <- create_sce(
  counts        = "path/to/counts.rds",
  cell_metadata = "path/to/metadata.rds"
)

sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)
sce <- normalize_counts(sce)
sce <- prepare_sce(sce, markers)

results <- run_cellvoter(sce)

consensus <- resolve_consensus_labels(
  label_list        = results$labels,
  method_names      = results$method_names,
  tie_breaker_names = results$tie_breaker_names,
  unassigned_label  = "unknown"
)

sce$cellVoteR_label  <- consensus$label
sce$cellVoteR_method <- consensus$method

table(sce$cellVoteR_label)

```
## Learn More

- [Introduction to CellVoteR](articles/intro_to_CellVoteR.html)
- [User inputs and required formats](articles/user_inputs.html)
- [Preprocessing and analysis tracks](articles/preprocessing.html)
- [Annotation methods](articles/annotation_methods.html)
- [Consensus and results](articles/consensus_and_results.html)
- [Troubleshooting](articles/troubleshooting.html)
