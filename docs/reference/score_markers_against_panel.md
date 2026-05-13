# Score cluster marker genes against fine cell type marker panels

For each cluster, identifies the relevant subset of the fine cell type
marker panel based on the broad label encoded in the cluster name, then
scores overlap using Fisher's exact test and overlap similarity.

## Usage

``` r
score_markers_against_panel(top_markers, marker_panel, background_genes)
```

## Arguments

- top_markers:

  Named list. Output of [`extract_top_markers`](extract_top_markers.md);
  each element must contain a `top_n` character vector.

- marker_panel:

  Named list of named lists. Top-level names are broad cell-type
  categories (e.g. `"immune"`, `"vasculature"`, `"other"`); inner names
  are fine cell types; values are character vectors of marker genes.

- background_genes:

  Character vector. Full set of genes in the dataset, used as the
  background population for Fisher's exact test.

## Value

A `data.frame` with columns:

- `cluster_name`:

  Cluster identifier from `top_markers`.

- `marker_set_name`:

  Fine cell type name from the panel.

- `fisher_p`:

  P-value from Fisher's exact test.

- `similarity`:

  Overlap similarity score.

Rows are ordered by `cluster_name` then ascending `fisher_p`.

## Details

Cluster names are expected to follow the convention `"<broad>_sc<N>"`,
e.g. `"immune_sc1"`, `"vasculature_sc2"`. The broad prefix is matched
against the top-level names of `marker_panel`.

**Edge case — collapsed broad labels:** When all clusters share the same
broad label (e.g. all assigned `"other"`), downstream subclustering
produces names like `"1_sc1"`, `"2_sc2"` whose prefix is numeric and
does not resolve to a panel entry. In this situation the function tests
against *all* marker sets in the panel.
