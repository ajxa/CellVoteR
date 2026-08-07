# Load and structure CellVoteR marker definitions

Reads a marker definition file (csv, tab-separated txt, xlsx, rds) or
receives an R [`data.frame`](https://rdrr.io/r/base/data.frame.html)
directly and organises it into a hierarchical list suitable for two-tier
(broad -\> fine) cell-type annotation pipelines. The returned broad
markers are a simple named list of character vectors intended to be
passed to
[`build_broad_marker_config()`](https://ajxa.github.io/CellVoteR/reference/build_broad_marker_config.md)
as a subsequent step to attach priority, threshold, and co-expression
settings.

## Usage

``` r
load_markers(
  file_path = NULL,
  markers = NULL,
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

- markers:

  Alternative argument for directly passing a `data.frame` containing
  marker definitions. If supplied, this overrides `file_path`.

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

## Supported input types

- `.csv`:

  Comma-separated values, read via
  [`read.csv`](https://rdrr.io/r/utils/read.table.html).

- `.txt`:

  Tab-separated values, read via
  [`read.delim`](https://rdrr.io/r/utils/read.table.html).

- `.xlsx`:

  Excel workbook, read via
  [`read.xlsx`](https://rdrr.io/pkg/openxlsx/man/read.xlsx.html).

- `.rds`:

  A saved R object, read via
  [`readRDS`](https://rdrr.io/r/base/readRDS.html). This must contain a
  `data.frame` matching the expected layout.

- `data.frame`:

  A `data.frame` passed directly via either the `file_path` or `markers`
  argument.

## Category reconciliation

Fine-type rows whose `category` does not match any broad category are
reassigned to `unnamed_broad_cat_label` (default `"other"`) with an
informational message. Broad categories that have **no** corresponding
fine rows cause an error.

## Typical workflow when loading external marker files


    markers <- load_markers("markers.csv")

    markers$broad <- build_broad_marker_config(
      marker_list            = markers$broad,
      priority_order         = c("vasculature", "immune"),
      per_category_overrides = list(
        immune = list(coexp_min = 2)
      )
    )

## Using internal marker panels

CellVoteR includes a set of curated marker panels for IDHwt glioblastoma
(GBM) under the `marker_panels$GBM` list entry. These panels can be
passed directly as a `data.frame` to `load_markers()`:


    markers <- load_markers(marker_panels$GBM$gbmap_neftel_full)

    # Or using the explicit markers parameter:
    markers <- load_markers(markers = marker_panels$GBM$gbmap_neftel_full)

## See also

[`build_broad_marker_config`](https://ajxa.github.io/CellVoteR/reference/build_broad_marker_config.md)
for configuring the broad markers returned by this function.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load from file
markers <- load_markers("markers.csv")
# Load from .rds file
markers <- load_markers("markers.rds")
# Load directly from data frame / internal package dataset
markers <- load_markers(marker_panels$GBM$gbmap_neftel_full)

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
