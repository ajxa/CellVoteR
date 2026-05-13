# Assign fine cell type labels from enrichment scores

Selects the best-matching fine cell type label for each cluster based on
minimum Fisher p-value, with ties broken by descending overlap
similarity.

## Usage

``` r
assign_fine_labels(cluster_col, scores)
```

## Arguments

- cluster_col:

  Factor or character vector. Cluster identifiers per cell, typically a
  `colData` column from a
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).
  Values must match `cluster_name` entries in `scores`.

- scores:

  A `data.frame` as returned by
  [`score_markers_against_panel`](score_markers_against_panel.md), with
  columns `cluster_name`, `marker_set_name`, `fisher_p`, and
  `similarity`.

## Value

A named list containing:

- `labels`:

  Factor of fine cell type labels, one per cell, in the same order as
  `cluster_col`.

- `scores`:

  A `data.frame` of the best label per cluster with columns
  `cluster_name`, `marker_set_name`, `fisher_p`, and `similarity`.
