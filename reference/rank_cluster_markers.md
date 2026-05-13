# Rank markers for each cluster using scran::findMarkers()

Runs [`findMarkers`](https://rdrr.io/pkg/scran/man/findMarkers.html) and
stores the resulting ranked marker tables in `metadata(sce)` for
downstream annotation, or optionally returns the ranked marker results
directly as a list.

## Usage

``` r
rank_cluster_markers(
  sce,
  cluster_col,
  assay_type = "logcounts",
  test_type = c("wilcox", "t"),
  direction = c("up", "down", "any"),
  pval_type = c("any", "some", "all"),
  min_prop = 0.25,
  metadata_key = "broad_cluster_markers",
  return_list = FALSE,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- sce:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

- cluster_col:

  Character scalar. Column in `colData(sce)` containing the cluster
  labels to test.

- assay_type:

  Character scalar. Assay to use. Defaults to `"logcounts"`.

- test_type:

  Character scalar. Differential expression test passed to
  [`findMarkers`](https://rdrr.io/pkg/scran/man/findMarkers.html). One
  of `"wilcox"` or `"t"`. Defaults to `"wilcox"`.

- direction:

  Character scalar. Direction of testing. One of `"up"`, `"down"`, or
  `"any"`. Defaults to `"up"`.

- pval_type:

  Character scalar. A string specifying how p-values are to be combined
  across pairwise comparisons for a given group/cluster. For more
  information see
  [`combineMarkers`](https://rdrr.io/pkg/scran/man/combineMarkers.html).
  This can be one of `"any"`, `"some"`, or `"all"`. Defaults to `"any"`.

- min_prop:

  Numeric scalar in `(0, 1]`. This sets the minimum proportion of
  pairwise comparisons in which a gene must be DE (in the expected
  direction) for it to contribute to that gene's summary statistic. This
  is only relevant when `pval_type` is set to `"any"`, and will ensure
  that a gene will only be high-ranked if it is among the top-ranked
  genes in at least `min.prop` of the pairwise comparisons. If the
  `pval_type` is set to `"some"` or `"all"`, then this will be ignored.
  For more information see
  [`combineMarkers`](https://rdrr.io/pkg/scran/man/combineMarkers.html).
  Defaults to 0.25, in order to be lenient and identify markers amongst
  similar populations.

- metadata_key:

  Character scalar. Name of the metadata entry used to store the ranked
  marker results when `return_list = FALSE`. Defaults to
  `"broad_cluster_markers"`.

- return_list:

  Logical scalar. If `FALSE` (default), store the ranked marker tables
  in `metadata(sce)[[metadata_key]]` and return the modified `sce`. If
  `TRUE`, return the output list directly instead of modifying `sce`.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object. Defaults to
  [`SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html).

## Value

If `return_list = FALSE`, the input `SingleCellExperiment` with ranked
marker tables stored in `metadata(sce)[[metadata_key]]`.

If `return_list = TRUE`, a named list with components:

- `marker_tables`:

  The output of
  [`findMarkers`](https://rdrr.io/pkg/scran/man/findMarkers.html).

- `params`:

  A list of parameters used to generate the ranked marker tables.

## Details

This function is intended as a reusable preprocessing step for:

- broad label assignment from ranked broad markers, and

- fine label assignment via Fisher's Exact Test on top-ranked genes.
