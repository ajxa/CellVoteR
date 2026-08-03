# Load and structure CellVoteR marker definitions

Reads a marker definition file (csv, tab-separated txt, or xlsx) and
organises it into a hierarchical list suitable for two-tier (broad -\>
fine) cell-type annotation pipelines. The returned broad markers are a
simple named list of character vectors intended to be passed to
[`build_broad_marker_config()`](https://ajxa.github.io/CellVoteR/reference/build_broad_marker_config.md)
as a subsequent step to attach priority, threshold, and co-expression
settings.

## Usage

``` r
load_markers(
  file_path,
  unnamed_broad_cat_label = "other",
  type_col = "type",
  unique_types = c("broad", "fine"),
  cat_col = "category",
  label_col = "label",
  marker_col = "marker"
)
```

## Arguments

- file_path:

  Character scalar. Path to a `.csv`, `.txt` (tab-separated), or `.xlsx`
  file containing marker definitions.

- unnamed_broad_cat_label:

  Character scalar. Label assigned to fine-type categories that do not
  map to any broad category. Defaults to `"other"`.

- type_col, cat_col, label_col, marker_col:

  Character scalars giving the column names in the input file for type,
  category, label, and marker respectively.

- unique_types:

  Character vector of permitted values in `type_col`.

## Value

A named list with components:

- broad:

  Named list of character vectors - one element per broad category,
  values are marker gene symbols. Pass this to
  [`build_broad_marker_config()`](https://ajxa.github.io/CellVoteR/reference/build_broad_marker_config.md)
  to generate the full configuration.

- fine:

  Named list of lists: top level keyed by category, second level keyed
  by label, values are character vectors of markers. The
  `unnamed_broad_cat_label` group (if any) is placed last.

## Expected file layout

The input file must contain at least four columns (names configurable):

- type:

  Either `"broad"` or `"fine"`.

- category:

  Broad cell-type category (e.g. `"immune"`, `"vasculature"`). Every
  broad category must also appear as a category in the fine rows.

- label:

  Fine-grained cell-type label (e.g. `"CD8_T"`, `"endothelial"`).

- marker:

  Gene symbol.

## Supported file formats

- `.csv`:

  Comma-separated values, read via
  [`read.csv`](https://rdrr.io/r/utils/read.table.html).

- `.txt`:

  Tab-separated values, read via
  [`read.delim`](https://rdrr.io/r/utils/read.table.html).

- `.xlsx`:

  Excel workbook, read via
  [`read.xlsx`](https://rdrr.io/pkg/openxlsx/man/read.xlsx.html).

## Category reconciliation

Fine-type rows whose `category` does not match any broad category are
reassigned to `unnamed_broad_cat_label` (default `"other"`) with an
informational message. Broad categories that have **no** corresponding
fine rows cause an error.

## Typical workflow


    markers <- load_markers("markers.csv")

    markers$broad <- build_broad_marker_config(
      marker_list    = markers$broad,
      priority_order = c("vasculature", "immune"),
      per_category_overrides = list(
        immune = list(coexp_min = 2)
      )
    )

## See also

[`build_broad_marker_config`](https://ajxa.github.io/CellVoteR/reference/build_broad_marker_config.md)
for configuring the broad markers returned by this function.

## Examples

``` r
if (FALSE) { # \dontrun{
markers <- load_markers("markers.csv")

# Inspect raw broad markers before configuring
names(markers$broad)
markers$broad[["immune"]]

# Then configure broad markers for annotation
markers$broad <- build_broad_marker_config(
  marker_list    = markers$broad,
  priority_order = c("vasculature", "immune")
)
} # }
```
