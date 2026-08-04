# CellVoteR

CellVoteR is an ensemble-based pipeline for robust cell type annotation
in single-cell RNA-seq data.

Rather than relying on a single best-match score, CellVoteR combines
complementary annotation strategies across full and marker-reduced
feature spaces, then resolves agreement and disagreement through a
configurable consensus voting step.

![CellVoteR workflow](reference/figures/cellvoter-workflow.svg)

------------------------------------------------------------------------

## Installation

You can install the development version from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("ajxa/CellVoteR")
```

## Quick Start

``` r

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

## Getting started

[Introduction to
CellVoteR](https://ajxa.github.io/CellVoteR/articles/intro_to_CellVoteR.md)  
A short walkthrough of the complete CellVoteR workflow.

[User inputs and required
formats](https://ajxa.github.io/CellVoteR/articles/user_inputs.md)  
How to structure count matrices, metadata, and marker files.

[Preprocessing and analysis
tracks](https://ajxa.github.io/CellVoteR/articles/preprocessing.md)  
QC, normalisation, marker overlap checks, and feature-space
construction.

## Annotation workflow

[Annotation
methods](https://ajxa.github.io/CellVoteR/articles/annotation_methods.md)  
The four primary methods and two global tie-breakers.

[Consensus and
results](https://ajxa.github.io/CellVoteR/articles/consensus_and_results.md)  
How final labels are resolved and inspected.

## Practical guidance

[Troubleshooting](https://ajxa.github.io/CellVoteR/articles/troubleshooting.md)  
Common issues with marker overlap, unresolved labels, and large
datasets.

## Source

[Browse source code](https://github.com/ajxa/CellVoteR)  
[Report a bug](https://github.com/ajxa/CellVoteR/issues)

## License

[Full license](https://ajxa.github.io/CellVoteR/LICENSE.md)

## Citation

[How to cite
CellVoteR](https://ajxa.github.io/CellVoteR/authors.html#citation)
