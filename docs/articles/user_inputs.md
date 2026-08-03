# User Inputs and Required Formats

  

> In order to run the complete CellVoteR workflow, you must supply two
> inputs:
>
> 1.  Raw gene-by-cell counts matrix (sparse `dgCMatrix`, `RDS`, or
>     `MTX triplet` file).
>
> 2.  Marker configuration — a structured list of broad and fine cell
>     type marker genes.

## Preparing input marker

Markers are the backbone of CellVoteR’s annotation strategy. These are
organised into two tiers:

### Broad markers (lineage-specific)

Broad markers are used to define coarse cell lineages (e.g. *Immune*,
*Vasculature*, *Other*) and therefore, must be:

- **Small** — typically 2–5 genes per category.

- **Mutually exclusive** — no gene should appear in more than one broad
  category.

- **Biologically diagnostic** — genes that robustly delineate lineages
  even in heterogeneous datasets.

These broad markers are loaded and then configured with
\[build_broad_marker_config()\], which assigns expression thresholds and
priority rankings used for tie-breaking when a cell passes multiple
broad categories.

### Fine markers (cell-type specific)

Fine markers define sub-populations that correspond directly to each
broad lineage (e.g. *B cell*, *T cell*, *NK cell* within the *Immune*
broad marker category). They marker sets can be larger and more
redundant than broad markers, as they are only used for Fisher’s Exact
Test scoring during fine annotation. Moreover, they do not need to be
mutually exclusive, but should still sufficiently distinguish between to
distinct cell types that come from the same broad lineage, e.g, T cells
vs. B cells and Mural cells vs Endothelial cells.

### Loading markers

User-supplied markers can be loaded from either `Excel`, `CSV`, or `TXT`
files. Irrespective of format, the file must be structured to contain
both broad and fine markers, with the following four columns:

> - **type** — indicates whether the marker is a *broad* or *fine*
>
> - **category** — for broad markers, this is the coarse lineage
>   category; for fine markers, this is the same as the broad category
>   they belong to.
>
> - **label** - for broad markers, this field is left blank; for fine
>   markers, this is the specific cell type label that the marker
>   corresponds to.
>
> - **marker** — the gene symbol of the marker gene.

  

An example of the required format is as follows:

| type  | category    | label       | marker |
|-------|-------------|-------------|--------|
| broad | immune      |             | PTPRC  |
| broad | vasculature |             | CDH5   |
| broad | vasculature |             | VWF    |
| fine  | immune      | T cell      | CD2    |
| fine  | immune      | T cell      | CD3D   |
| fine  | immune      | T cell      | IL32   |
| fine  | immune      | B cell      | CD79A  |
| fine  | immune      | B cell      | CD79B  |
| fine  | vasculature | Mural cell  | IGFBP7 |
| fine  | vasculature | Mural cell  | FN1    |
| fine  | vasculature | Endothelial | A2M    |
| fine  | vasculature | Endothelial | IGFBP7 |

> When defining broad category markers, there is no need to assign a
> **label** as these are only used for coarse lineage assignment.
> However, for the fine markers, the **label** field is required as
> these are used to assign the final cell type labels to clusters and
> cells.

``` r

markers <- load_markers(file_path = "path/to/input_markers.xlsx")

# Inspect the structure
str(markers$broad)   # named list of character vectors
str(markers$fine)    # nested named list: broad category > fine cell type > genes
```

### Configuring broad markers

The \[build_broad_marker_config()\] function processes the raw broad
marker list, attaching expression thresholds and priority ranks used
during the enrichment-based annotation step.

``` r


markers$broad <- build_broad_marker_config(
  marker_list       = markers$broad,
  priority_order    = c("vasculature", "immune"),   # higher priority listed first
  default_threshold = 0.25                          # default logcounts threshold
)

str(markers$broad$immune)
```

> The **priority_order** argument is used to assign numeric priority
> ranks to each broad category, which is used for tie-breaking when a
> cell passes expression thresholds for more than one broad category.
> Categories listed earlier receive a lower (higher priority) numeric
> rank.

## Preparing expression data

CellVoteR is designed to work natively with
[*SingleCellExperiment*](https://bioconductor.org/packages/SingleCellExperiment/)
(SCE) objects. The \[create_sce()\] helper function can be used to
construct a basic SCE from a variety of input formats, which is then
processed and annotated by the downstream functions.

### From sparse matrix (in-memory)

The most direct way to create an SCE is from an in-memory sparse matrix
of class `dgCMatrix` (genes x cells), with optional cell metadata as a
`data.frame`:

``` r

sce <- create_sce(
  counts        = my_sparse_matrix,  # dgCMatrix, genes x cells
  cell_metadata = my_metadata_df     # data.frame, one row per cell (optional)
)
```

### From a file path

This function also accepts file paths, which is useful for large
datasets where the matrix is stored on disk: currently, disk-backed
matrices are not supported but this is a planned feature for a future
release:

``` r


sce <- create_sce(
  counts        = "path/to/counts.rds",
  cell_metadata = "path/to/metadata.rds"   # also accepts .csv or .tsv
)
```

## From a 10X Genomics output

The \[create_sce()\] function can also directly import from the standard
10X Genomics output directory, which contains the `matrix.mtx.gz`,
`barcodes.tsv`, and `features.tsv` files:

``` r


sce <- create_sce(
  mtx_file   = "path/to/matrix.mtx.gz",
  cells_file = "path/to/barcodes.tsv",
  genes_file = "path/to/features.tsv"
)
```

## Inputs quality control

The \[assess_cell_quality()\] helper function can be used to calculate
per-cell QC metrics and optionally remove low-quality cells before
downstream analysis. This is usually the first step after creating the
SCE, as it ensures that downstream processes are not affected by
low-quality cells that may have aberrant expression profiles.

``` r

sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)
```

> If the **remove_failed_cells** argument is set to **FALSE**, the QC
> metrics are still calculated and attached to the SCE, but no cells are
> removed. This allows you to inspect the QC metrics and decide on
> thresholds for filtering before actually removing any cells.

## Input normalisation

Following quality control, the raw counts matrix must be normalised to
account for differences in sequencing depth and other technical factors.
This is a critical step, as the annotation methods rely on accurate
expression values to score marker gene enrichment. CellVoteR uses a
pooling-based normalisation strategy ([Lun et
al. 2016](https://doi.org/10.1186/s13059-016-0947-7)) via
`scran::computePooledFactors`, followed by log-normalisation with
[`scuttle::logNormCounts`](https://rdrr.io/pkg/scuttle/man/logNormCounts.html).
This is implemented in the \[normalize_counts()\] helper function, which
takes the raw SCE and returns a new SCE with normalised logcounts.

``` r

sce <- normalize_counts(sce)
```

> The function adds a *logcounts* assay to the input SCE which contains
> the normalised expression values. The original raw counts are retained
> in the *counts* assay for reference and potential use in downstream
> processes that require raw counts (e.g. DE testing). The normalised
> logcounts are used for all downstream analyses, including HVG
> selection, PCA, clustering, and marker scoring.

## Building analysis tracks

The final preparatory step before annotation is to build the analysis
tracks - the feature spaces and cluster structures on which the
annotation methods depend. This is done by the \[prepare_sce()\]
function, which takes the normalised SCE and the marker configuration as
inputs and returns an SCE with all necessary metadata, reduced
dimensions, and cluster assignments attached for downstream annotation.

Internally, the \[prepare_sce()\] function applies a series of
processing steps:

1.  Validates broad and fine marker configurations against the
    expression matrix.

2.  Builds two independent feature spaces — the **full HVG space** and
    the **reduced marker-defined space**.

3.  Runs PCA and unsupervised clustering (Leiden via SNN graph) on each
    space.

4.  Attaches the marker configuration and filtered fine markers to the
    SCE metadata.

5.  Stores the reduced feature space as an **altExp** named
    *“user_panel”*.

``` r

sce <- prepare_sce(sce, markers)
```

The resulting SCE contains two parallel cluster structures: one based on
the full HVG space, and one based on the reduced marker-defined space.
Both are used for different annotation methods, as some methods rely on
broad clusters defined in the full HVG space, while others use clusters
defined in the reduced space to ensure they are driven by the supplied
markers.

The reduced marker-defined space is stored within the same SCE object as
an alternative experiment (altExp) named *“user_panel”*. This altExp
contains the log-normalised expression values for the subset of genes
that overlap with the supplied markers, and is used for PCA and
clustering in the reduced feature space.

The overall structure of the SCE after running
[`prepare_sce()`](https://ajxa.github.io/CellVoteR/reference/prepare_sce.md)
is as follows:

    sce
    ├── assays: counts, logcounts
    ├── rowSubset("broad_hvgs")           ←    HVGs used for broad clustering
    ├── reducedDim("PCA_broad_hvg")       ←    PCA on full HVG space
    ├── colData$cluster_broad_hvg         ←    Leiden clusters (full space)
    ├── colData$cluster_broad_reduced     ←    Leiden clusters (reduced space)
    ├── metadata$marker_config            ←    full marker configuration
    ├── metadata$filtered_fine_markers    ←    fine markers present in data
    ├── metadata$missing_by_label         ←    per-label missing marker report
    └── altExp("user_panel")              ←    reduced feature SCE
        ├── assays: counts, logcounts
        ├── reducedDim("PCA")
        ├── colData$cluster
        ├── metadata$marker_config
        ├── metadata$filtered_fine
        └── metadata$params

### Automatic parameter estimation

During the preparation step, key parameters for PCA and clustering are
estimated automatically based on the number of cells in the dataset,
using bounded $`log/sqrt`$ scaling via the \[estimate_cluster_params()\]
helper. This ensures that the feature spaces and cluster structures are
appropriately sized for the dataset at hand, without requiring manual
tuning of parameters. However, if you wish to override the automatically
estimated parameters, you can do so by passing custom values for the
number of HVGs, PCs, SNN neighbourhood size $`k`$, and Leiden resolution
directly to \[prepare_sce()\] as follows:

``` r

sce <- prepare_sce(
  sce,
  markers,
  n_hvgs     = 3000L,
  n_pcs      = 30L,
  k          = 20L,
  resolution = 0.8
)
```

### Marker overlap reporting

In order to ensure that the annotation methods have sufficient
information to work with, \[prepare_sce()\] checks the supplied marker
lists against the expression matrix and reports on the degree of
overlap. If a large proportion of fine markers are missing from the
dataset, this may indicate an issue with gene naming conventions or
platform-specific gene capture, which could affect annotation
performance. The function captures this information and details the
labels that are most affected, which can be inspected as follows:

``` r

metadata(sce)$missing_by_label
```

Markers that are missing from the dataset are automatically removed from
the fine marker sets used for scoring, so annotation still proceeds with
the genes that are present. This only occurs when the total number of
missing markers is \<50% (default )of the total distinct markers
supplied - this parameter can be relaxed or tightened accordingly, by
altering the `overlap_feat_percent` argument:

``` r


prepare_sce(
  sce,
  markers,
  overlap_feat_percent = 75       # more stringent
)
```
