# Troubleshooting

This article collects common issues and practical checks for CellVoteR
workflows.

## Low marker overlap

If
[`prepare_sce()`](https://ajxa.github.io/CellVoteR/reference/prepare_sce.md)
warns about low fine marker overlap, you can inspect which labels are
missing markers with:

``` r

metadata(sce)$missing_by_label
```

Consider whether the missing genes are platform-specific (e.g. not
captured by your assay technology), or whether alternative gene symbols
should be used.

## High unresolved rate

If many cells are labelled `"unknown"` after consensus, try:

``` r

# 1. Allow even splits
consensus <- resolve_consensus_labels(
  ...,
  allow_even_split = TRUE
)

# 2. Disable ordered tie-breaking so either tie-breaker can resolve
consensus <- resolve_consensus_labels(
  ...,
  ordered_tiebreak = FALSE
)

# 3. Inspect which method combinations are causing disagreements
table(results$labels$method_1, results$labels$method_2)
```

## Collapsed broad labels

If all clusters receive the same broad label, CellVoteR retains the
original cluster structure automatically. This is expected behaviour for
highly homogeneous datasets (e.g. a sample consisting entirely of tumour
cells). In this case, the numeric cluster prefixes (e.g. `1_sc1`,
`2_sc1`) trigger testing against the full fine marker panel rather than
a lineage-specific subset.

## Large datasets

For datasets exceeding available RAM, convert the SCE to HDF5-backed
storage after
[`create_sce()`](https://ajxa.github.io/CellVoteR/reference/create_sce.md):

``` r

HDF5Array::saveHDF5SummarizedExperiment(sce, dir = "my_hdf5_sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment("my_hdf5_sce")
```

Whilst this has been implemented and tested, it is not yet a fully
supported workflow. Please report any issues which may arise.

### Parallelisation

Key functions accept a `BPPARAM` argument for parallelisation via
[`BiocParallel`](https://bioconductor.org/packages/BiocParallel/):

``` r

library(BiocParallel)

results <- run_cellvoter(
  sce,
  annotation_args = list(
    broad_args = list(BPPARAM = MulticoreParam(4)),
    rank_args  = list(BPPARAM = MulticoreParam(4))
  )
)
```
