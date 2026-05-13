# Prepare Analysis Tracks on a SingleCellExperiment

Runs the complete preprocessing pipeline for both feature spaces
required for CellVoteR directly on the SCE object using standard
Bioconductor functions.

## Usage

``` r
prepare_sce(
  sce,
  marker_config,
  n_hvgs = 2000L,
  overlap_feat_percent = 50,
  n_pcs = NULL,
  k = NULL,
  resolution = NULL,
  cluster_params_args = list(),
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- sce:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  with `logcounts` assay (from
  [`normalize_counts`](normalize_counts.md)).

- marker_config:

  Named list with `$broad` (from
  [`build_broad_marker_config`](build_broad_marker_config.md)) and
  `$fine` (from [`load_markers`](load_markers.md)).

- n_hvgs:

  Integer scalar. Number of HVGs. Defaults to `2000`.

- overlap_feat_percent:

  Numeric scalar (0-100). Minimum fine marker overlap. Defaults to `50`.

- n_pcs, k, resolution:

  Override automatic parameter estimation. `NULL` (default) uses
  [`estimate_cluster_params`](estimate_cluster_params.md).

- cluster_params_args:

  Named list passed to
  [`estimate_cluster_params`](estimate_cluster_params.md).

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html).
  Defaults to
  [`SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html).

## Value

The input SCE with both analysis tracks populated. The `marker_config`
that is stored in the `metadata` of the `altExp` differs from that
stored in the main SCE in that the fine marker sets are filtered to only
include the genes present in the user panel and that overlap with the
broad marker sets. The broad HVG track is built using the full
`logcounts` assay and the full set of HVGs, while the user panel track
is built using a subset of features based on the fine markers that
overlap with the user panel and broad markers. The dimensionality
reduction and clustering for each track are performed independently
using their respective feature spaces, and their parameters are
estimated separately using
[`estimate_cluster_params`](estimate_cluster_params.md).

## Storage layout


    sce
    |-- assays: counts, logcounts
    |-- rowSubset("broad_hvgs")
    |-- reducedDim("PCA_broad_hvg")
    |-- colData$cluster_broad_hvg
    |-- metadata$broad_hvg_params
    |-- metadata$marker_config
    |-- metadata$filterd_fine_markers
    |-- metadata$missing_by_label
    `-- altExp("user_panel")
        |-- assays: counts, logcounts
        |-- reducedDim("PCA")
        |-- colData$cluster
        |-- metadata$marker_config
        |-- metadata$filterd_fine_markers
        `-- metadata$params

## Examples

``` r
if (FALSE) { # \dontrun{
sce <- normalize_counts(sce)
sce <- prepare_sce(sce, marker_config)

# Broad HVG track
reducedDim(sce, "PCA_broad_hvg")
sce$cluster_broad_hvg

# User panel track
reducedDim(altExp(sce, "user_panel"), "PCA")
altExp(sce, "user_panel")$cluster
} # }
```
