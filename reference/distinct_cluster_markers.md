# Extract distinct markers per cluster

Converts a data frame of differential expression (DE) results into a
named list, where each element represents a cluster and contains a
vector of unique gene symbols.

## Usage

``` r
distinct_cluster_markers(
  de_results,
  cluster_col = "cluster",
  gene_col = "gene"
)
```

## Arguments

- de_results:

  A data frame containing DE results (e.g., from Seurat's
  `FindAllMarkers`).

- cluster_col:

  Character. The name of the column containing cluster identities.
  Default is "cluster".

- gene_col:

  Character. The name of the column containing gene symbols. Default is
  "gene".

## Value

A named list of character vectors. Names correspond to cluster IDs.
