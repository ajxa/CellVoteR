# Build broad marker configuration

Converts a named list of broad-category marker vectors into a structured
configuration list suitable for priority-based cell-type assignment.
Each category is annotated with an expression threshold, minimum
co-expression count, and a numeric priority rank derived from
`priority_order`.

## Usage

``` r
build_broad_marker_config(
  marker_list,
  priority_order,
  default_threshold = 0.1,
  default_coexp = 1,
  per_category_overrides = NULL
)
```

## Arguments

- marker_list:

  Named list of character vectors. Names are broad category labels (e.g.
  `"immune"`, `"vasculature"`), values are marker gene symbols.

- priority_order:

  Character vector defining the assignment priority. Categories listed
  first receive lower (higher-priority) rank values. All entries must
  correspond to names in `marker_list`. At most one category in
  `marker_list` may be absent from this vector (see Priority order
  section).

- default_threshold:

  Numeric scalar (\>= 0). Default expression threshold applied to every
  category. Defaults to `0.1`.

- default_coexp:

  Positive integer scalar. Minimum number of co-expressed markers
  required for a category call. Defaults to `1`. Deprecated in this
  version, but retained for backward compatibility - will be removed in
  a future release.

- per_category_overrides:

  Optional named list of named lists, keyed by category name. Each inner
  list may contain `expr_threshold` and/or `coexp_min` to override the
  defaults for that category. Categories not listed use the defaults.
  Defaults to `NULL` (no overrides).

## Value

A named list (one element per category, ordered by priority rank
ascending) where each element is a list with components:

- markers:

  Character vector of marker gene symbols.

- expr_threshold:

  Numeric expression threshold for the category.

- coexp_min:

  Integer minimum co-expression count.

- priority:

  Integer priority rank (1 = highest priority).

## Priority order

The `priority_order` vector determines how ties are resolved when a cell
or cluster matches multiple broad categories. Categories listed earlier
receive higher priority (lower rank number). At most one category may be
omitted from `priority_order`, in which case it is automatically
assigned the lowest priority. If two or more categories are omitted they
would share the same rank, making tie-breaking ambiguous - this is
treated as an error.

## Examples

``` r
if (FALSE) { # \dontrun{
  markers <- list(
    immune      = c("CD45", "CD3", "CD8"),
    vasculature = c("PECAM1", "VWF"),
    stromal     = c("COL1A1", "VIM")
  )

  # All categories ranked explicitly
  cfg <- build_broad_marker_config(
    marker_list    = markers,
    priority_order = c("vasculature", "immune", "stromal")
  )

  # One category omitted - automatically gets lowest priority
  cfg <- build_broad_marker_config(
    marker_list    = markers,
    priority_order = c("vasculature", "immune")
  )
  # vasculature=1, immune=2, stromal=3

  # Per-category overrides
  cfg <- build_broad_marker_config(
    marker_list    = markers,
    priority_order = c("vasculature", "immune", "stromal"),
    per_category_overrides = list(
      immune = list(coexp_min = 2, expr_threshold = 0.2)
    )
  )
} # }
```
