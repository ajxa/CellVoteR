# Troubleshooting

## Overview

This article collects common issues and practical checks for CellVoteR
workflows.

## Low Marker Overlap

If
[`prepare_sce()`](https://ajxa.github.io/CellVoteR/reference/prepare_sce.md)
reports low marker overlap, inspect the missing marker report:

``` r

metadata(sce)$missing_by_label
```

Common causes include gene symbol mismatches, platform-specific gene
capture, or marker sets from a different species or annotation source.

## High Unresolved Rate

If many cells are assigned the unresolved label, inspect method-level
agreement:

``` r

table(results$labels$method_1, results$labels$method_2)
table(results$labels$method_3, results$labels$method_4)
```

You can also re-run consensus with less conservative settings:

``` r

consensus <- resolve_consensus_labels(
  label_list        = results$labels,
  method_names      = results$method_names,
  tie_breaker_names = results$tie_breaker_names,
  unassigned_label  = "unknown",
  allow_even_split  = TRUE,
  ordered_tiebreak  = FALSE
)
```

## Collapsed Broad Labels

If all clusters receive the same broad label, CellVoteR retains the
original cluster structure automatically. This can be expected in highly
homogeneous datasets.

## Large Datasets

For datasets that exceed available memory, consider storing the object
in HDF5-backed format after creating the initial SCE:

``` r

HDF5Array::saveHDF5SummarizedExperiment(sce, dir = "my_hdf5_sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment("my_hdf5_sce")
```

## Parallelisation

Some internal steps support parallel execution through `BiocParallel`:

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
