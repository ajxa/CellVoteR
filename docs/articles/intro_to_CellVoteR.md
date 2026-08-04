# Introduction to CellVoteR

## Overview

CellVoteR is an ensemble-based pipeline for robust cell type annotation
in single-cell RNA-seq data.

Rather than assigning labels from a single marker score or clustering
strategy, CellVoteR runs complementary annotation methods across two
feature spaces and resolves the results through a configurable consensus
step. The workflow has four main stages:

1.  Prepare a SingleCellExperiment object from raw counts.
2.  Load and configure broad and fine marker sets.
3.  Build full and marker-reduced analysis tracks.
4.  Run ensemble annotation and resolve final consensus labels.

## Installation

``` r

# install.packages("devtools")
devtools::install_github("ajxa/CellVoteR")
```

``` r

library(CellVoteR)
```

------------------------------------------------------------------------

## Complete Workflow

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
table(sce$cellVoteR_method)
```

## What Happens Internally?

At a high level, CellVoteR compares four primary annotation methods and
two global tie-breakers. The four primary methods combine:

- two feature spaces: full highly variable gene space and reduced
  marker-panel space
- two broad annotation strategies: cluster-based and enrichment-based

The consensus step then assigns a final label when the methods agree
strongly, or uses the global tie-breakers to resolve partial
disagreement.

## Where to Go Next?

- For details on input file formats and marker structure, see [User
  inputs and required
  formats](https://ajxa.github.io/CellVoteR/articles/user_inputs.md).

- For information on preprocessing, normalisation and preparing the
  analysis tracks, see [Preprocessing and analysis
  tracks](https://ajxa.github.io/CellVoteR/articles/preprocessing.md).

- For information on the individual annotation methods, see [Annotation
  methods](https://ajxa.github.io/CellVoteR/articles/annotation_methods.md)
  .

- For details on the consensus resolution step, see [Consensus and
  results](https://ajxa.github.io/CellVoteR/articles/consensus_and_results.md).

- For any additional practical guidance, including troubleshooting tips
  and common issues, see
  [Troubleshooting](https://ajxa.github.io/CellVoteR/articles/troubleshooting.md).
