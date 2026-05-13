# CellVoteR

An ensemble-based pipeline for robust cell type annotation in
single-cell RNA-seq data.

CellVoteR moves beyond simple “best-match” scoring by integrating four
complementary annotation strategies across two feature spaces, resolved
through a principled consensus voting step to generate a high-confidence
label for every cell.

## How It Works

CellVoteR determines cell identity through a three-stage process:

### 1. Broad Triage

Each cell is first assigned a broad lineage label (e.g. *Immune*,
*Vasculature*, *Other*) using one of two strategies:

- **Cluster-based:** Unsupervised clusters are labelled by testing
  whether curated broad marker genes are significantly and consistently
  top-ranked within each cluster.
- **Enrichment-based:** Individual cells are labelled directly by
  aggregating expression across small, distinct broad marker sets and
  applying category-specific thresholds.

### 2. Targeted Sub-clustering & Fine Annotation

Each broad lineage is sub-clustered independently. Sub-cluster identity
is then determined by scoring the top ranked marker genes against a
curated fine-resolution marker panel using Fisher’s Exact Test and
overlap similarity. This ensures, for example, that immune sub-clusters
are only compared against immune markers — reducing false positives and
improving resolution of rare populations.

Both strategies are applied to two feature spaces:

- **Full gene set** — the complete normalised expression matrix
- **Reduced gene set** — a targeted panel of user-supplied marker genes

This gives four primary annotation methods in total.

### 3. Global Consensus & Tie-breaking

Two global tie-breakers are computed in parallel without broad triage,
by clustering the full and reduced gene sets directly and scoring
against the fine marker panel. These serve as independent reference
votes.

All six method outputs are then passed to a configurable ensemble voting
step which applies a decision hierarchy to resolve disagreements and
assign a final label to each cell.

## Workflow Diagram

``` mermaid
graph TD
    classDef default fill:#ffffff,stroke:#333,stroke-width:1px,color:#000
    classDef start fill:#e1f5fe,stroke:#01579b,stroke-width:2px,color:#000
    classDef track fill:none,stroke:#999,stroke-width:2px,stroke-dasharray: 5 5
    classDef title fill:#eceff1,stroke:#455a64,stroke-width:2px,font-weight:bold,color:#000
    classDef decision fill:#fff9c4,stroke:#fbc02d,stroke-width:2px,color:#000
    classDef method fill:#e0f2f1,stroke:#00695c,stroke-width:1px,color:#000
    classDef tiebreak fill:#fce4ec,stroke:#880e4f,stroke-width:1px,color:#000
    classDef logic fill:#fff3e0,stroke:#e65100,stroke-width:1px,color:#000
    classDef endnode fill:#dcedc8,stroke:#33691e,stroke-width:2px,color:#000

    Input([Raw Counts + Markers]):::start
    Input --> QC[assess_cell_quality]
    QC --> Norm[normalize_counts]
    Norm --> Prep[prepare_sce\nBuilds full + reduced feature spaces\nand attaches marker config]

    Prep --> M1
    Prep --> M2
    Prep --> M3
    Prep --> M4
    Prep --> G1
    Prep --> G2

    subgraph PrimaryMethods [ ]
        direction TB
        TitleP[Primary Annotation Methods]:::title

        subgraph Full[Full Gene Set]
            M1[Method 1\nCluster-based\nBroad triage → subcluster\n→ Fisher score]:::method
            M3[Method 3\nEnrichment-based\nBroad triage → subcluster\n→ Fisher score]:::method
        end

        subgraph Reduced[Reduced Gene Set]
            M2[Method 2\nCluster-based\nBroad triage → subcluster\n→ Fisher score]:::method
            M4[Method 4\nEnrichment-based\nBroad triage → subcluster\n→ Fisher score]:::method
        end

        TitleP ~~~ Full
        TitleP ~~~ Reduced
    end

    subgraph TieBreakerTrack [ ]
        direction TB
        TitleT[Global Tie-breakers]:::title

        G1[Tie-breaker 1\nHVG clusters\nFull gene set\n→ Fisher score]:::tiebreak
        G2[Tie-breaker 2\nPanel clusters\nReduced gene set\n→ Fisher score]:::tiebreak

        TitleT ~~~ G1
        TitleT ~~~ G2
    end

    M1 & M2 & M3 & M4 --> Vote
    G1 & G2 --> Vote

    Vote[resolve_consensus_labels\nEnsemble voting]

    Vote --> Majority{Majority\nvote?}:::decision
    Majority -->|Yes| Final
    Majority -->|No, but\nleading candidate| TB{Tie-breakers\nagree?}:::decision
    TB -->|Both agree| Final
    TB -->|Priority order| Final
    TB -->|Neither resolves| Unresolved[Unassigned]:::logic
    Majority -->|All methods\ndisagree| Unresolved

    Final([Final Cell Label]):::endnode

    class PrimaryMethods,TieBreakerTrack track
```

## Installation

You can install the development version of CellVoteR with:

``` r

# install.packages("devtools")
devtools::install_github("ajxa/CellVoteR")
```

## Quick Start

### 1. Prepare inputs

CellVoteR requires two inputs: a raw counts matrix and a marker
configuration.

``` r

library(CellVoteR)

# Load and configure markers
markers <- load_markers(file_path = "path/to/input_markers.xlsx")

markers$broad <- build_broad_marker_config(
  marker_list    = markers$broad,
  priority_order = c("vasculature", "immune"),
  default_threshold = 0.25
)

# Create SCE from a sparse matrix or file path
sce <- create_sce(
  counts        = "path/to/counts.rds",
  cell_metadata = "path/to/metadata.rds"  # optional, also accepts .csv / .tsv
)
```

### 2. QC and preprocessing

``` r

sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)
sce <- normalize_counts(sce)

# Builds both feature spaces (full + reduced), clusters them, and
# attaches the marker configuration to the SCE
sce <- prepare_sce(sce, markers)
```

### 3. Run the ensemble

``` r

# Default run — uses all four methods and both tie-breakers
results <- run_cellvoter(sce)

# With custom annotation parameters
results <- run_cellvoter(
  sce,
  return_full_output = TRUE,
  annotation_args = list(
    broad_args   = list(test_type = "t", min_prop = 0.1),
    rank_args    = list(test_type = "t", min_prop = 0.25),
    extract_args = list(fdr_threshold = 0.01, target_n = 50L)
  )
)
```

### 4. Resolve consensus

The consensus step is intentionally separate so voting parameters can be
tweaked and re-run without repeating the annotation pipeline.

``` r

consensus <- resolve_consensus_labels(
  label_list        = results$labels,
  method_names      = results$method_names,
  tie_breaker_names = results$tie_breaker_names,
  unassigned_label  = "unknown",
  allow_even_split  = FALSE,
  ordered_tiebreak  = TRUE
)

# Attach final labels to the SCE
sce$cellVoteR_label  <- consensus$label
sce$cellVoteR_method <- consensus$method

# Inspect results
table(sce$cellVoteR_label)
table(sce$cellVoteR_method)

# Inspect per-method labels before consensus
table(results$labels$method_1)
table(results$labels$method_2)
table(results$labels$method_3)
table(results$labels$method_4)
table(results$labels$global_1)
table(results$labels$global_2)

# Inspect per-cluster scores for a method (requires return_full_output = TRUE)
results$full_output$method_1$scores
```

### 5. Custom markers

CellVoteR is tissue-agnostic. Broad markers should be small (ideally ≤ 3
genes), mutually exclusive across categories, and biologically
diagnostic of a lineage. Fine markers can be larger gene sets and are
used for Fisher scoring.

``` r

# Markers are loaded from an Excel workbook or nested YAML file
markers <- load_markers(file_path = "path/to/my_tissue_markers.xlsx")
```

See [`?load_markers`](reference/load_markers.md) and
[`?build_broad_marker_config`](reference/build_broad_marker_config.md)
for details on the expected input format.

## Key Design Decisions

| Decision | Rationale |
|----|----|
| Broad triage before fine annotation | Prevents dominant lineages from masking rare populations |
| Two feature spaces (full + reduced) | Full HVG space captures global structure; marker-defined space sharpens lineage boundaries |
| Four methods + two tie-breakers | Redundancy across strategies reduces sensitivity to any single method’s failure mode |
| Consensus step is separate | Users can re-tune voting parameters without re-running the full pipeline |
| Fisher’s Exact Test for scoring | Statistically principled enrichment test that accounts for background gene set size |
