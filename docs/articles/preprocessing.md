# Preprocessing and Analysis Tracks

## Overview

This article describes the preprocessing steps used before running
CellVoteR annotation.

The main steps are:

1.  Create a `SingleCellExperiment` object.
2.  Assess cell quality.
3.  Normalise counts.
4.  Build full and reduced analysis tracks with
    [`prepare_sce()`](https://ajxa.github.io/CellVoteR/reference/prepare_sce.md).

## Create a SingleCellExperiment

``` r

sce <- create_sce(
  counts        = "path/to/counts.rds",
  cell_metadata = "path/to/metadata.rds"
)
```

## Quality Control

``` r

sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)
```

## Normalisation

``` r

sce <- normalize_counts(sce)
```

## Build Analysis Tracks

``` r

sce <- prepare_sce(sce, markers)
```

## Automatic Parameter Estimation

[`prepare_sce()`](https://ajxa.github.io/CellVoteR/reference/prepare_sce.md)
can estimate clustering and dimensionality-reduction parameters
automatically based on dataset size.

``` r

sce <- prepare_sce(
  sce,
  markers,
  n_hvgs     = 3000L,
  n_pcs      = 30L,
  k          = 20L,
  resolution = 0.8
)
```

## Marker Overlap Reporting

``` r

metadata(sce)$missing_by_label
```

## Resulting Object Structure

After preprocessing, the `SingleCellExperiment` contains the assays,
reduced dimensions, cluster labels, marker configuration, and reduced
marker-panel experiment needed by the annotation workflow.
