# Create/Load a SingleCellExperiment from raw counts

Entry point for the CellVoteR annotation workflow. Accepts a raw
(un-normalised) gene-by-cell counts matrix and constructs a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object ready for downstream QC, normalisation, and clustering.

## Usage

``` r
create_sce(
  counts = NULL,
  mtx_file = NULL,
  cells_file = NULL,
  genes_file = NULL,
  cell_metadata = NULL,
  gene_metadata = NULL
)
```

## Arguments

- counts:

  A gene-by-cell raw counts input. Either:

  - A `dgCMatrix` (in-memory sparse matrix with rownames and colnames).

  - A character scalar path to a `.rds` file containing a `dgCMatrix`.

  When provided, takes precedence over the MTX arguments.

- mtx_file:

  Character scalar. Path to a Matrix Market (`.mtx` or `.mtx.gz`) file.
  Must be accompanied by `cells_file` and `genes_file`.

- cells_file:

  Character scalar. Path to a text file containing one cell identifier
  per line (no header). Required when `mtx_file` is provided.

- genes_file:

  Character scalar. Path to a text file containing one gene identifier
  per line (no header). Required when `mtx_file` is provided.

- cell_metadata:

  Optional per-cell annotations. Either:

  - A `data.frame` with one row per cell.

  - A character scalar path to an `.rds` file containing a `data.frame`.

  - A character scalar path to a `.csv`, `.tsv`, or `.txt` file
    (tab-separated).

  Must have one row per column in the counts matrix. This will be stored
  in `colData()`.

- gene_metadata:

  Optional `data.frame` of per-gene annotations. Must have one row per
  row in the counts matrix. Stored in `rowData()`.

## Value

A
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
with a single `counts` assay stored as a `dgCMatrix`.

## Input precedence

The function accepts two mutually exclusive input modes. If `counts` is
provided it takes precedence and the MTX arguments are ignored (with a
warning if also supplied). Otherwise, all three MTX arguments must be
provided together.

## Supported input formats

- **counts (in-memory)**:

  A `dgCMatrix` with gene identifiers as rownames and cell barcodes as
  colnames.

- **counts (RDS file)**:

  A character path to a `.rds` file containing a `dgCMatrix` as
  described above.

- **MTX triplet files**:

  Three file paths supplied via `mtx_file`, `cells_file`, and
  `genes_file`. The MTX file contains the sparse matrix in Matrix Market
  format. The cells and genes files are text files with one identifier
  per line (or tab-separated, in which case the second column is used
  preferentially). All three must be provided together.

- **cell_metadata (RDS or CSV/TSV file)**:

  A character path to an `.rds` file containing a `data.frame`, or a
  path to a `.csv` or `.tsv`/`.txt` file. The file must contain one row
  per cell, in the same order as the columns of the counts matrix.

## Disk-backed extension

The returned SCE holds an in-memory `dgCMatrix` counts assay. For
datasets that exceed available RAM, the object can be converted to
disk-backed storage after creation:


    sce <- create_sce(counts = "large_counts.rds")
    HDF5Array::saveHDF5SummarizedExperiment(sce, dir = "my_hdf5_sce")
    sce <- HDF5Array::loadHDF5SummarizedExperiment("my_hdf5_sce")

## See also

[`load_markers`](https://ajxa.github.io/CellVoteR/reference/load_markers.md)
and
[`build_broad_marker_config`](https://ajxa.github.io/CellVoteR/reference/build_broad_marker_config.md)
for preparing the marker configuration to attach to the SCE via
`metadata()`.

## Examples

``` r
if (FALSE) { # \dontrun{
# From an in-memory sparse matrix
sce <- create_sce(counts = my_sparse_counts)

# From an RDS file with metadata from a CSV
sce <- create_sce(
  counts        = "raw_counts.rds",
  cell_metadata = "cell_metadata.csv"
)

# From MTX triplet files with metadata from an RDS
sce <- create_sce(
  mtx_file      = "data/matrix.mtx.gz",
  cells_file    = "data/barcodes.tsv",
  genes_file    = "data/features.tsv",
  cell_metadata = "data/cell_metadata.rds"
)
} # }
```
