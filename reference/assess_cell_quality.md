# Assess Cell Quality Metrics

Calculates comprehensive quality control (QC) metrics for a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object, including feature counts, library size, and specific gene group
percentages (e.g. mitochondrial or ribosomal genes). Supports filtering
based on minimum features, adaptive thresholds for library size, and
user-defined feature groups.

## Usage

``` r
assess_cell_quality(
  sce,
  check_feature_groups = NULL,
  min_cells_per_sample = 100,
  min_features = 100,
  nmads = 3,
  sample_col = NULL,
  assay_name = "counts",
  remove_failed_cells = FALSE
)
```

## Arguments

- sce:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object with a `counts` assay.

- check_feature_groups:

  A named list defining gene groups to check. Each element must be a
  list with `features` (character vector of gene names) and `max_pct`
  (numeric threshold for maximum percentage). If `NULL` (default),
  auto-detects mitochondrial (`^MT-`) and ribosomal (`^RP[SL]`) genes
  via
  [`find_qc_features`](https://ajxa.github.io/CellVoteR/reference/find_qc_features.md).

- min_cells_per_sample:

  Integer scalar. Minimum number of cells required per sample for that
  sample's cells to be retained. Requires `sample_col` to operate.
  Defaults to `100`.

- min_features:

  Integer scalar. Minimum number of detected genes required for a cell
  to pass QC. Defaults to `100`.

- nmads:

  Numeric scalar. Number of median absolute deviations below the median
  log-transformed library size to set the adaptive low-count threshold.
  Stored for reference but not included in the pass/fail vote. Defaults
  to `3`.

- sample_col:

  Character scalar naming a column in `colData(sce)` containing sample
  identifiers, or a character/factor vector of the same length as the
  number of cells. If `NULL`, all cells are treated as a single sample.

- assay_name:

  Character scalar. Name of the assay to use. Defaults to `"counts"`.

- remove_failed_cells:

  Logical. If `TRUE`, subsets the returned object to retain only cells
  with `QC_PASS == TRUE`. Defaults to `FALSE`.

## Value

The input SCE with additional columns in `colData()`:

- `sc_n_features`:

  Number of detected genes per cell.

- `sc_n_counts`:

  Total UMI counts per cell.

- `pass_min_features`:

  Logical. `TRUE` if cell exceeds `min_features`.

- `sc_adaptive_threshold`:

  The computed low-count threshold.

- `pass_adaptive_counts`:

  Logical. `TRUE` if cell counts exceed the adaptive threshold.

- `pct_[group]`:

  Percentage of counts from each feature group.

- `pass_max_[group]`:

  Logical. `TRUE` if percentage is below `max_pct` for that group.

- `sc_sample_id`:

  Sample identity for each cell.

- `pass_min_cells_sample`:

  Logical. `TRUE` if cell belongs to a sample meeting
  `min_cells_per_sample`.

- `QC_PASS`:

  Logical. `TRUE` if the cell passed all active checks.

## QC strategy

Each cell is assessed against multiple independent criteria. A cell must
pass **all** active checks to receive `QC_PASS = TRUE`:

1.  Minimum detected features (`min_features`).

2.  Feature group percentage thresholds (e.g. mitochondrial \< 20\\

3.  Minimum cells per sample (`min_cells_per_sample`).

An adaptive library-size threshold is also computed and stored for
reference but is **not** included in the pass/fail vote by default.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with auto-detection
sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)

# Custom feature groups
groups <- list(
  mito   = list(features = c("MT-ND1", "MT-CO1"), max_pct = 15),
  stress = list(features = c("JUN", "FOS"), max_pct = 5)
)
sce <- assess_cell_quality(sce,
  check_feature_groups = groups,
  min_features = 200
)

# With per-sample filtering using a colData column
sce <- assess_cell_quality(sce, sample_col = "patient_id")
} # }
```
