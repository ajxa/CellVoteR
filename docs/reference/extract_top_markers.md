# Extract top marker genes from ranked marker tables

Generates a list of marker genes per cluster using ranked marker tables
produced by
[`rank_cluster_markers`](https://ajxa.github.io/CellVoteR/reference/rank_cluster_markers.md).

## Usage

``` r
extract_top_markers(
  sce = NULL,
  ranked_markers_key = "broad_cluster_markers",
  ranked_markers = NULL,
  fdr_threshold = 0.05,
  effect_threshold = 0.6,
  target_n = 100L
)
```

## Arguments

- sce:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).
  Required only when `ranked_markers` is not supplied.

- ranked_markers_key:

  Character scalar. Metadata entry containing ranked marker tables.
  Defaults to `"broad_cluster_markers"`.

- ranked_markers:

  List or `NULL`. Ranked marker result returned by
  [`rank_cluster_markers`](https://ajxa.github.io/CellVoteR/reference/rank_cluster_markers.md)`(return_list = TRUE)`.
  If supplied, this takes precedence over `sce` and
  `ranked_markers_key`.

- fdr_threshold:

  Numeric scalar. Maximum FDR threshold. Defaults to `0.05`.

- effect_threshold:

  Numeric scalar. Minimum effect size threshold. Defaults to `0.6`.

- target_n:

  Integer scalar. Number of genes to return per cluster. Defaults to
  `100L`.

## Value

A named list with one entry per cluster, each containing:

- `top_n`:

  Character vector of up to `target_n` genes, passing genes first, then
  supplemented backfill.

- `supplemented`:

  Character vector of backfill genes that did not meet the FDR/effect
  thresholds. Empty character vector if all `target_n` slots were filled
  by passing genes.

## Details

Ranked marker tables may be supplied either:

- directly via the `ranked_markers` argument, or

- indirectly via `metadata(sce)[[ranked_markers_key]]`.

If `ranked_markers` is supplied, it takes precedence and `sce` is only
used for class consistency. In that case, `ranked_markers_key` is
ignored.

Genes are first ordered by:

1.  ascending `FDR`, and

2.  descending effect size (AUC for Wilcoxon, logFC for t-tests).

Marker genes passing:

- `FDR <= fdr_threshold`, and

- `effect size >= effect_threshold`

are prioritised. The final output for each cluster consists of exactly
`target_n` genes: passing genes are included first, then remaining slots
are filled with top-ranked genes that did not pass thresholds.

This ensures a consistent gene set size across clusters while
prioritising high-confidence markers for downstream enrichment analysis.
