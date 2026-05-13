# Normalise counts via pooling-based size factors

Computes size factors using
[`computePooledFactors`](https://rdrr.io/pkg/scuttle/man/computePooledFactors.html)
and applies log-normalisation via
[`logNormCounts`](https://rdrr.io/pkg/scuttle/man/logNormCounts.html).
If a `logcounts` assay already exists it is overwritten with a warning.

## Usage

``` r
normalize_counts(
  sce,
  min_cluster_size = 100L,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `counts` assay.

- min_cluster_size:

  Integer scalar. Minimum cluster size for
  [`quickCluster`](https://rdrr.io/pkg/scran/man/quickCluster.html) used
  to seed the pooling. Defaults to `100`.

- BPPARAM:

  A `BiocParallelParam` for parallelisation.

## Value

The input SCE with `sizeFactors()` set and a `logcounts` assay added.
