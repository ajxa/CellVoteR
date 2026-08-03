# Assign broad labels based on the median rank of validated markers

Assigns broad cluster labels using ranked marker tables produced by
[`rank_cluster_markers`](https://ajxa.github.io/CellVoteR/reference/rank_cluster_markers.md).

## Usage

``` r
label_broad_clusters(
  sce,
  broad_config = NULL,
  marker_config_key = "marker_config",
  ranked_markers = NULL,
  ranked_markers_key = "broad_cluster_markers",
  cluster_col = "cluster_broad_hvg",
  label_col = "broad_cluster",
  fdr_threshold = 0.05,
  effect_threshold = 0.6,
  unassigned_label = "other"
)
```

## Arguments

- sce:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

- broad_config:

  Named list or `NULL`. Validated broad marker definitions as produced
  by
  [`build_broad_marker_config`](https://ajxa.github.io/CellVoteR/reference/build_broad_marker_config.md).
  If `NULL`, extracted from
  `metadata(sce)[[marker_config_key]][["broad"]]`.

- marker_config_key:

  Character scalar. Metadata key containing the marker configuration.
  Only used when `broad_config = NULL`. Defaults to `"marker_config"`.

- ranked_markers:

  List or `NULL`. Ranked marker result from
  [`rank_cluster_markers`](https://ajxa.github.io/CellVoteR/reference/rank_cluster_markers.md)`(return_list = TRUE)`.
  If `NULL`, extracted from `metadata(sce)[[ranked_markers_key]]`.

- ranked_markers_key:

  Character scalar. Metadata key containing ranked marker tables. Only
  used when `ranked_markers = NULL`. Defaults to
  `"broad_cluster_markers"`.

- cluster_col:

  Character scalar. `colData` column containing cluster identifiers.
  Defaults to `"cluster_broad_hvg"`.

- label_col:

  Character scalar. Name of the output `colData` column. Defaults to
  `"broad_cluster"`.

- fdr_threshold:

  Numeric scalar. Maximum FDR for a marker to be considered significant.
  Defaults to `0.05`.

- effect_threshold:

  Numeric scalar. Minimum effect size (AUC for Wilcoxon, logFC for
  t-test). Defaults to `0.6`.

- unassigned_label:

  Character scalar. Label assigned to clusters with no passing markers.
  Defaults to `"other"`.

## Value

The input `sce` with a new `colData` column `label_col` containing
assigned broad labels.

## Details

Ranked marker tables may be supplied either:

- directly via the `ranked_markers` argument, or

- indirectly via `metadata(sce)[[ranked_markers_key]]`.

Broad marker definitions may be supplied either:

- directly via the `broad_config` argument, or

- indirectly via `metadata(sce)[[marker_config_key]][["broad"]]`.

For each cluster, genes are ordered by ascending `FDR` then descending
effect size (`summary.AUC` for Wilcoxon, `summary.logFC` for t-tests).
Within each broad category, only markers passing both
`FDR <= fdr_threshold` and `effect size >= effect_threshold` are
retained. The category score is the median rank of passing markers; the
best category is chosen by lowest median rank with user-supplied
`priority` as a tie-breaker.

## Collapsed broad label edge case

If all clusters receive the same broad label but more than one
unsupervised cluster exists, the original cluster labels are retained in
`label_col` rather than collapsing to a single label. This preserves
cluster structure for downstream subclustering.
