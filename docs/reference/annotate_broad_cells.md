# Assign cell labels using broad-marker enrichment

Assigns broad labels to individual cells using aggregated log-normalised
expression across small curated marker sets for each broad category.

## Usage

``` r
annotate_broad_cells(
  sce,
  marker_config_key = "marker_config",
  label_col = "broad_enrichment",
  assay_name = "logcounts",
  aggregate_fun = c("sum", "mean", "median"),
  other_label = "other"
)
```

## Arguments

- sce:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  containing a `logcounts` assay.

- marker_config_key:

  Character. Metadata entry containing a named list of validated marker
  definitions. Defaults to `"marker_config"`.

- label_col:

  Character scalar. Name of the output column to create in
  `colData(sce)`. Defaults to `"broad_enrichment"`.

- assay_name:

  Character scalar. Assay to use for expression values. Defaults to
  `"logcounts"`.

- aggregate_fun:

  Character scalar specifying how marker expression should be aggregated
  within each category. One of `"sum"`, `"mean"`, or `"median"`.
  Defaults to `"sum"`.

- other_label:

  Character scalar. Label assigned when no category passes. Defaults to
  `"other"`.

## Value

The input `SingleCellExperiment` with:

- `colData(sce)[[label_col]]`:

  A factor of per-cell broad labels, with `other_label` as the last
  level.

- `metadata(sce)$broad_cell_enrichment`:

  A list containing the per-cell category score matrix, logical pass
  matrix, and assignment parameters.

## Details

For each category, expression is aggregated across that category's
marker genes using one of `"sum"`, `"mean"`, or `"median"`. A category
is considered to pass for a cell if the aggregated expression exceeds
the category-specific `expr_threshold`. If no categories pass, the cell
is labelled as `other_label`. If more than one category passes, the tie
is resolved using the user-supplied `priority`, with lower numeric
values treated as higher priority.

This function assumes that `broad_config` has already been validated in
upstream preprocessing, including checks that all broad markers are
present, marker sets do not overlap, and priorities are unique.
