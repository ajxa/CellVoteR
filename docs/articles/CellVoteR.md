# CellVoteR: Ensemble Cell Type Annotation for Single-Cell RNA-seq

## Overview

CellVoteR is an ensemble-based pipeline for robust cell type annotation
in single-cell RNA-seq (scRNA-seq) data. Rather than relying on a single
classification strategy, CellVoteR integrates four complementary
annotation methods across two feature spaces, then resolves
disagreements through a principled consensus voting step.

The core design philosophy is:

- **Divide and conquer** — broadly triage cells into lineages before
  applying fine-resolution annotation, preventing dominant populations
  from masking rare cell types

- **Redundancy** — four methods running in parallel reduces sensitivity
  to the failure modes of any single approach

- **Separation of concerns** — annotation (slow, compute-intensive) and
  consensus resolution (fast, parameter-sensitive) are decoupled, so the
  user can re-tune voting without repeating the pipeline

In order to run the complete workflow, CellVoteR requires two inputs to
be supplied:

1.  A **raw gene-by-cell counts matrix** (sparse `dgCMatrix`, `RDS`, or
    `MTX triplet` file).

2.  A **marker configuration** — a structured list of broad and fine
    cell type marker genes.

------------------------------------------------------------------------

## Installation

Currently, the package can be installed directly from Github:

``` r

# install.packages("devtools")
devtools::install_github("ajxa/CellVoteR")
```

``` r

library(CellVoteR)
```

------------------------------------------------------------------------

## Detailed Pipeline Steps

### Step 1: Preparing Marker Inputs

Markers are the backbone of CellVoteR’s annotation strategy. They are
organised into two tiers:

#### Broad markers (lineage-specific)

Broad markers define coarse cell lineages (e.g. *Immune*, *Vasculature*,
*Other*). Therefore, they must be:

- **Small** — typically 2–5 genes per category.

- **Mutually exclusive** — no gene should appear in more than one broad
  category.

- **Biologically diagnostic** — genes that robustly delineate lineages
  even in heterogeneous datasets.

These broad markers are loaded and then configured with
[`build_broad_marker_config()`](../reference/build_broad_marker_config.md),
which assigns expression thresholds and priority rankings used for
tie-breaking when a cell passes multiple broad categories.

#### Fine markers (cell-type specific)

Fine markers define sub-populations within each broad lineage (e.g. *B
cell*, *T cell*, *NK cell* within *Immune*). They can be larger gene
sets and are used for Fisher’s Exact Test scoring during fine
annotation. These marker do not need to be mutually exclusive, but
should sufficiently distinguish between to cell types from a common
lineage, e.g, T cells vs B cells and Mural cells vs Endothelial cells.

#### Loading markers

User-supplied markers can be loaded from either `Excel`, `CSV`, or `TXT`
files. The files must be structured to comprise four columns: `type`
(broad/fine), `category`, `label`, and `marker`:

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

> When defining broad category markers, leave the **label** field blank.
> For fine cell type markers within that broad category, assign a
> **label**.

``` r

markers <- load_markers(file_path = "path/to/input_markers.xlsx")

# Inspect the structure
str(markers$broad)   # named list of character vectors
str(markers$fine)    # nested named list: broad category > fine cell type > genes
```

#### Configuring broad markers

[`build_broad_marker_config()`](../reference/build_broad_marker_config.md)
processes the raw broad marker list, attaching expression thresholds and
priority ranks used during the enrichment-based annotation step.

``` r

markers$broad <- build_broad_marker_config(
  marker_list       = markers$broad,
  priority_order    = c("vasculature", "immune"),  # higher priority listed first
  default_threshold = 0.25                          # default logcounts threshold
)

# Each broad category now has markers, expr_threshold, coexp_min, and priority
str(markers$broad$immune)
```

> The **priority_order** argument controls tie-breaking when a cell
> passes expression thresholds for more than one broad category -
> categories listed earlier receive a lower (higher priority) numeric
> rank.

------------------------------------------------------------------------

### Step 2: Object Creation & QC

CellVoteR works natively with
[*SingleCellExperiment*](https://bioconductor.org/packages/SingleCellExperiment/)
objects. Use [`create_sce()`](../reference/create_sce.md) to construct
one from your raw data.

#### From (in-memory) sparse matrix

``` r

sce <- create_sce(
  counts        = my_sparse_matrix,  # dgCMatrix, genes x cells
  cell_metadata = my_metadata_df     # data.frame, one row per cell (optional)
)
```

#### From file path

[`create_sce()`](../reference/create_sce.md) also accepts file paths,
which is useful for large datasets where the matrix is stored on disk:

``` r

# From RDS files
sce <- create_sce(
  counts        = "path/to/counts.rds",
  cell_metadata = "path/to/metadata.rds"   # also accepts .csv or .tsv
)

# From MTX triplet files
sce <- create_sce(
  mtx_file   = "path/to/matrix.mtx.gz",
  cells_file = "path/to/barcodes.tsv",
  genes_file = "path/to/features.tsv"
)
```

------------------------------------------------------------------------

### Step 3: Quality Control

[`assess_cell_quality()`](../reference/assess_cell_quality.md)
calculates per-cell QC metrics and optionally removes low-quality cells
before downstream analysis.

``` r

sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)
```

> Cells that fail QC are flagged in **colData(sce)\$QC_PASS**. Setting
> **remove_failed_cells = TRUE** subsets the object to passing cells
> only.

------------------------------------------------------------------------

### Step 4: Normalisation

CellVoteR uses a pooling-based normalisation strategy ([Lun et
al. 2016](https://doi.org/10.1186/s13059-016-0947-7)) via
`scran::computePooledFactors`, followed by log-normalisation with
[`scuttle::logNormCounts`](https://rdrr.io/pkg/scuttle/man/logNormCounts.html).

``` r

sce <- normalize_counts(sce)
```

> This adds a *logcounts* assay to the SCE and sets `sizeFactors()`. The
> *logcounts* assay is required by all downstream steps.

------------------------------------------------------------------------

### Step 5: Building Analysis Tracks

The [`prepare_sce()`](../reference/prepare_sce.md) function co-ordinates
the key preprocessing step. This function performs/applies, the
following calculations and logic:

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

After this step, the initial SingleCellExepriment object structure is as
follows:

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

#### Automatic parameter estimation

Clustering parameters (number of PCs, SNN neighbourhood size $`k`$,
Leiden resolution) are estimated automatically from cell count using
bounded $`log/sqrt`$ scaling via
[`estimate_cluster_params()`](../reference/estimate_cluster_params.md).
However, these parameters can be overridden if required:

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

#### Marker overlap reporting

In some cases, you may have fine markers that are partially missing from
your dataset: the [`prepare_sce()`](../reference/prepare_sce.md)
function captures this information and details the labels that are most
affected. which can be inspected as follows:

``` r

# Inspect missing marker report after prepare_sce()
metadata(sce)$missing_by_label
```

Missing (in the data) marker are automatically removed from the fine
marker sets used for scoring, so annotation still proceeds with the
genes that are present. This only occurs when the total number of
missing markers is \<50% of the total distinct markers supplied - this
parameter can be relaxed or tightened accordingly, by altering the
`overlap_feat_percent` argument:

``` r


prepare_sce(
  sce,
  markers,
  overlap_feat_percent = 75       # more stringent
)
```

------------------------------------------------------------------------

### Step 6: Annotation Methods

The [`run_cellvoter()`](../reference/run_cellvoter.md) orchestrates all
of the individual annotation pipelines and returns a named list of
per-cell label factors, one per method:

``` r

results <- run_cellvoter(sce)
```

#### What this does internally

Four primary methods and two global tie-breakers are run:

| Method        | Feature space   | Broad strategy   | Subcluster feature mode |
|---------------|-----------------|------------------|-------------------------|
| Method 1      | Full (HVG)      | Cluster-based    | HVG                     |
| Method 2      | Reduced (panel) | Cluster-based    | All                     |
| Method 3      | Full (HVG)      | Enrichment-based | HVG                     |
| Method 4      | Reduced (panel) | Enrichment-based | All                     |
| Tie-breaker 1 | Full (HVG)      | None (global)    | —                       |
| Tie-breaker 2 | Reduced (panel) | None (global)    | —                       |

Each primary method follows the same pipeline:

    broad annotation
          ↓
    subcluster_labels()
          ↓
    rank_cluster_markers()          ←   DE testing per sub-cluster
          ↓
    extract_top_markers()           ←   select top N genes per cluster
          ↓
    score_markers_against_panel()   ←   Fisher's Exact Test + overlap similarity
          ↓
    assign_fine_labels()            ←   best label per cluster    →    mapped to cells

#### Broad annotation strategies

**Cluster-based** (`annotate_broad_clusters`): Runs DE testing on
pre-existing unsupervised clusters. Each cluster is assigned the broad
category whose curated markers have the lowest median rank among
significantly up-regulated genes (FDR ≤ 0.05, AUC ≥ 0.6 by default).

**Enrichment-based** (`annotate_broad_cells`): Assigns labels directly
to individual cells by aggregating expression across each broad
category’s marker genes and comparing against category-specific
thresholds. Does not depend on clustering.

> In some datasets (highly homogeneous tumour sample), all clusters may
> receive the same broad label. CellVoteR detects this and retains the
> original cluster structure for subclustering rather than collapsing to
> a single group, preserving the fine-resolution information for
> downstream processes.

#### Annotation Method Results

After running [`run_cellvoter()`](../reference/run_cellvoter.md), the
returned SingleCellExeriment object is populated with all of th
intermediate cluster labels, which can be accessed by querying the
colData() columns of the object:

| Column | Description |
|----|----|
| **cluster_broad_hvg** | Pre-existing HVG clusters from [`prepare_sce()`](../reference/prepare_sce.md) |
| **cluster_broad_reduced** | Pre-existing reduced clusters from [`prepare_sce()`](../reference/prepare_sce.md) |
| **broad_cluster_m1** | Broad labels, method 1 |
| **broad_cluster_sub_m1** | Sub-cluster labels, method 1 |
| **broad_cluster_m2** | Broad labels, method 2 |
| **broad_cluster_sub_m2** | Sub-cluster labels, method 2 |
| **broad_enrichment_m3** | Broad labels, method 3 |
| **broad_enrichment_sub_m3** | Sub-cluster labels, method 3 |
| **broad_enrichment_m4** | Broad labels, method 4 |
| **broad_enrichment_sub_m4** | Sub-cluster labels, method 4 |

#### Customising parameters

The main [`run_cellvoter()`](../reference/run_cellvoter.md) function
accepts an `annotation_args` list which can be used to customise various
parameters of the underlying internal functions, if required. There
following parameter lists that can be specified, which each control
specific aspects of the underlying logic:

- **rank_args** Controls parameters associated with
  [`rank_cluster_markers()`](../reference/rank_cluster_markers.md):

  - assay_type
  - test_type
  - direction
  - pval_type
  - min_prop
  - BPPARAM

- **broad_args** is identical to the `rank_args`, but only controls the
  ranking inside the `annotate_broad_clusters` function. This is useful
  for altering the parameters which control how the broad cell lineage
  labels (e.g, immune, vasculature) are defined - these broad markers
  are small and highly specific and so you may wish to use a more
  lenient `min_prop` compared to the fine labels and this allows for
  such fine control of the process.

- **extract_args** Controls parameters associated with
  [`extract_top_markers()`](../reference/extract_top_markers.md):

  - fdr_threshold
  - effect_threshold
  - target_n

An following is an example of how this can be used to independently
alter the broad and fine labelling logic:

``` r

results <- run_cellvoter(
  sce,
  return_full_output = TRUE,
  annotation_args = list(

    # Controls ranking inside annotate_broad_clusters (methods 1 and 2)
    # Lenient min_prop appropriate for small, specific broad marker sets
    broad_args = list(
      test_type = "wilcox",
      min_prop  = 0.1
    ),

    # Controls ranking inside .run_fine_annotation (all six methods)
    # Stricter settings for fine sub-cluster marker extraction
    rank_args = list(
      test_type = "wilcox",
      min_prop  = 0.25
    ),

    # Controls top marker extraction
    extract_args = list(
      fdr_threshold    = 0.05,
      effect_threshold = 0.6,
      target_n         = 100L
    )
  )
)
```

> **broad_args** vs **rank_args**: these both control
> [`rank_cluster_markers()`](../reference/rank_cluster_markers.md) but
> at different pipeline stages. *broad_args* affects how broad lineage
> labels are assigned (methods 1 and 2 only); *rank_args* affects how
> sub-cluster marker genes are extracted for Fisher scoring (all six
> methods). Methods 3 and 4 use *rank_args* only, as the
> enrichment-based broad step does not call *rank_cluster_markers()*.

#### Accessing full outputs

> Setting `return_full_output` = **TRUE**, returns the per-cluster
> Fisher scores and similarity values are for every method (by default
> this is set to *FALSE*).

``` r

# Per-cluster score table for method 1
results$full_output$method_1$scores

# Per-cluster score table for tie-breaker 2
results$full_output$global_2$scores
```

------------------------------------------------------------------------

### Step 7: Resolving Consensus Labels

Consensus resolution is intentionally a separate step. This means you
can adjust voting parameters and re-run
[`resolve_consensus_labels()`](../reference/resolve_consensus_labels.md)
as many times as needed without re-running the annotation pipeline.

``` r

consensus <- resolve_consensus_labels(
  label_list        = results$labels,
  method_names      = results$method_names,
  tie_breaker_names = results$tie_breaker_names,
  unassigned_label  = "unknown",
  allow_even_split  = FALSE,
  ordered_tiebreak  = TRUE
)
```

#### Decision hierarchy

For each cell, the following logic is applied in order:

1.  **Strong majority** — if one label receives strictly more than 50%
    of method votes (3 of 4 by default) it is assigned immediately.

2.  **Tie-breaker agreement** — if the two tie-breakers agree with each
    other *and* their shared label matches the leading candidate, that
    label is assigned.

3.  **Ordered tie-breaking** (`ordered_tiebreak = TRUE`, default) —
    tie-breaker 1 is tried first; if it matches the leading candidate
    the label is assigned. Otherwise tie-breaker 2 is tried.

4.  **Either tie-breaking** (`ordered_tiebreak = FALSE`) — either
    tie-breaker agreeing with the leading candidate is sufficient, with
    no priority between them.

5.  **Unresolved** — if every method disagrees (no leading candidate
    exists) or no tie-breaker resolves the split, `unassigned_label` is
    assigned.

#### Voting parameters

| Parameter | Default | Effect |
|----|----|----|
| **allow_even_split** | **FALSE** | When **TRUE**, a 2-of-4 plurality is accepted as majority |
| **ordered_tiebreak** | **TRUE** | When **FALSE**, either tie-breaker can resolve a split |
| **unassigned_label** | **“Unknown”** | Label for unresolved cells |

#### Attaching labels to the SCE

``` r

sce$cellVoteR_label  <- consensus$label
sce$cellVoteR_method <- consensus$method
```

------------------------------------------------------------------------

### Step 8: Inspecting Results

#### Summary tables

``` r

# Final label distribution
table(sce$cellVoteR_label)

# Decision method breakdown
table(sce$cellVoteR_method)

# Cross-tabulate label vs decision method
table(sce$cellVoteR_label, sce$cellVoteR_method)
```

#### Per-method label distributions

Comparing individual method outputs before consensus is useful for
understanding where methods agree and disagree:

``` r

table(results$labels$method_1)   # cluster-based, full
table(results$labels$method_2)   # cluster-based, reduced
table(results$labels$method_3)   # enrichment-based, full
table(results$labels$method_4)   # enrichment-based, reduced
table(results$labels$global_1)   # tie-breaker 1
table(results$labels$global_2)   # tie-breaker 2
```

#### Method agreement

A useful diagnostic is to look at how often all four primary methods
agree:

``` r

label_df <- as.data.frame(results$labels[1:4])

# Proportion of cells where all four methods agree
mean(apply(label_df, 1, function(x) length(unique(x)) == 1L))

# Per-cell method agreement count
label_df$n_agree <- apply(label_df, 1, function(x) max(table(x)))
table(label_df$n_agree)
```

#### Re-running consensus with different parameters

``` r

# allowing even splits
consensus_liberal <- resolve_consensus_labels(
  label_list        = results$labels,
  method_names      = results$method_names,
  tie_breaker_names = results$tie_breaker_names,
  unassigned_label  = "unknown",
  allow_even_split  = TRUE
)

table(consensus_liberal$label)

# Compare unresolved rate between settings
mean(consensus$label == "unknown")
mean(consensus_liberal$label == "unknown")
```

------------------------------------------------------------------------

## Complete Workflow

The full pipeline from raw data to annotated SCE:

``` r

library(CellVoteR)

# ── 1. Load and configure markers ──────────────────────────────────────────────
markers <- load_markers(file_path = "path/to/input_markers.xlsx")

markers$broad <- build_broad_marker_config(
  marker_list       = markers$broad,
  priority_order    = c("vasculature", "immune"),
  default_threshold = 0.25
)

# ── 2. Create SCE ──────────────────────────────────────────────────────────────
sce <- create_sce(
  counts        = "path/to/counts.rds",
  cell_metadata = "path/to/metadata.rds"
)

# ── 3. QC ──────────────────────────────────────────────────────────────────────
sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)

# ── 4. Normalise ───────────────────────────────────────────────────────────────
sce <- normalize_counts(sce)

# ── 5. Build analysis tracks ───────────────────────────────────────────────────
sce <- prepare_sce(sce, markers)

# ── 6. Run ensemble annotation ─────────────────────────────────────────────────
results <- run_cellvoter(sce)

# ── 7. Resolve consensus ───────────────────────────────────────────────────────
consensus <- resolve_consensus_labels(
  label_list        = results$labels,
  method_names      = results$method_names,
  tie_breaker_names = results$tie_breaker_names,
  unassigned_label  = "unknown"
)

# ── 8. Attach labels ───────────────────────────────────────────────────────────
sce$cellVoteR_label  <- consensus$label
sce$cellVoteR_method <- consensus$method

# ── 9. Inspect ─────────────────────────────────────────────────────────────────
table(sce$cellVoteR_label)
table(sce$cellVoteR_method)
```

------------------------------------------------------------------------

## Tips and Troubleshooting

### Low marker overlap

If [`prepare_sce()`](../reference/prepare_sce.md) warns about low fine
marker overlap, inspect which labels are most affected:

``` r

metadata(sce)$missing_by_label
```

Consider whether the missing genes are platform-specific (e.g. not
captured by your assay technology), or whether alternative gene symbols
should be used.

#### High unresolved rate

If many cells are labelled `"unknown"` after consensus, try:

``` r

# 1. Allow even splits
consensus <- resolve_consensus_labels(
  ...,
  allow_even_split = TRUE
)

# 2. Disable ordered tie-breaking so either tie-breaker can resolve
consensus <- resolve_consensus_labels(
  ...,
  ordered_tiebreak = FALSE
)

# 3. Inspect which method combinations are causing disagreements
table(results$labels$method_1, results$labels$method_2)
```

### Collapsed broad labels

If all clusters receive the same broad label, CellVoteR retains the
original cluster structure automatically. This is expected behaviour for
highly homogeneous datasets (e.g. a sample consisting entirely of tumour
cells). In this case, the numeric cluster prefixes (e.g. `1_sc1`,
`2_sc1`) trigger testing against the full fine marker panel rather than
a lineage-specific subset.

### Large datasets

For datasets exceeding available RAM, convert the SCE to HDF5-backed
storage after [`create_sce()`](../reference/create_sce.md):

``` r

HDF5Array::saveHDF5SummarizedExperiment(sce, dir = "my_hdf5_sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment("my_hdf5_sce")
```

### Parallelisation

Key functions accept a `BPPARAM` argument for parallelisation via
[`BiocParallel`](https://bioconductor.org/packages/BiocParallel/):

``` r

library(BiocParallel)

results <- run_cellvoter(
  sce,
  annotation_args = list(
    broad_args = list(BPPARAM = MulticoreParam(4)),
    rank_args  = list(BPPARAM = MulticoreParam(4))
  )
)
```

------------------------------------------------------------------------

## Session Info

    #> R version 4.4.1 (2024-06-14)
    #> Platform: aarch64-apple-darwin20
    #> Running under: macOS 26.3.1
    #> 
    #> Matrix products: default
    #> BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    #> LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    #> 
    #> locale:
    #> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    #> 
    #> time zone: Europe/London
    #> tzcode source: internal
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices datasets  utils     methods   base     
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] cli_3.6.6           knitr_1.51          rlang_1.2.0        
    #>  [4] xfun_0.54           otel_0.2.0          renv_1.0.7         
    #>  [7] textshaping_1.0.4   jsonlite_2.0.0      htmltools_0.5.9    
    #> [10] ragg_1.5.0          sass_0.4.10         rmarkdown_2.31     
    #> [13] evaluate_1.0.5      jquerylib_0.1.4     rmdformats_1.0.4   
    #> [16] fastmap_1.2.0       yaml_2.3.11         lifecycle_1.0.5    
    #> [19] bookdown_0.46       BiocManager_1.30.27 compiler_4.4.1     
    #> [22] fs_1.6.6            htmlwidgets_1.6.4   rstudioapi_0.17.1  
    #> [25] systemfonts_1.3.1   digest_0.6.39       R6_2.6.1           
    #> [28] bslib_0.9.0         tools_4.4.1         pkgdown_2.2.0      
    #> [31] cachem_1.1.0        desc_1.4.3
