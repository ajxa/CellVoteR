# CellVoteR <img src="https://github.com/ajxa/CellVoteR/blob/main/man/figures/main.png?raw=true" align="right" width="120" />

[ ](https://www.google.com/search?q=https://github.com/ajxa/CellVoteR/actions) [ ](https://opensource.org/licenses/MIT) **CellVoteR** is an ensemble-based pipeline for robust cell type classification in single-cell RNA-seq data.

It integrates multiple classification strategies - including clustering-based, marker-thresholding, and hierarchical triage—to generate a consensus label for every cell. By combining diverse methodologies ("voting"), CellVoteR resolves ambiguities that single methods often miss, providing a confidence-aware annotation workflow.

## Overview

CellVoteR is designed to be flexible and generalizable, offering:

-   **Ensemble Voting:** Runs up to 6 distinct classification methods and calculates a consensus vote.
-   **Clash Resolution:** Automatically resolves ties using raw marker intensity scores.
-   **Hierarchical Triage:** Optional logic to split cells into broad categories (e.g., Immune vs. Endothelial) before fine-grained sub-clustering.
-   **Quality Control:** Built-in functions to assess and filter cells based on feature counts and mitochondrial/ribosomal content.
-   **Extensible References:** Includes a curated IDHwt GBM marker panel but supports custom user-supplied marker lists for any tissue.

## Installation

You can install the development version of CellVoteR with:

``` r
# install.packages("devtools")
devtools::install_github("ajxa/CellVoteR")
```

## Quick Start

### 1. Quality Control

Before labelling, ensure your data is clean using the built-in QC engine.

``` r
library(CellVoteR)
library(Seurat)

# Load your data (Seurat, SingleCellExperiment, or Matrix)
counts <- Read10X(data.dir = "path/to/data")

# Assess quality (calculates mito/ribo % automatically)
qc_metrics <- assess_cell_quality(counts, 
                                  min_features = 200, 
                                  max_mito_pct = 20)

# Filter the object
clean_obj <- filter_cells(counts, qc_metrics)
```

### 2. Run the Ensemble

Run the full pipeline on your filtered object. You can use the default global methods or enable the full ensemble.

``` r
# Basic Run (General Purpose)
# Uses global clustering and marker matching
labelled_obj <- run_ensemble(clean_obj)

# Advanced Run (Full Ensemble)
# Runs triage methods (1-4) + global breakers (5-6) and resolves clashes
labelled_obj <- run_ensemble(clean_obj, use_ensemble = TRUE)

# View Results
table(labelled_obj$Ensemble_Resolved)
```

### 3. Custom Markers

CellVoteR works with any tissue. Simply provide a dataframe of markers.

``` r
# Define your own markers
my_markers <- data.frame(
  gene = c("EPCAM", "COL1A1", "CD3D"),
  cell_type = c("Epithelial", "Stromal", "T_Cell"),
  category = c("Epithelial", "Stromal", "Immune") 
)

# Run with custom panel
results <- run_ensemble(clean_obj, markers = my_markers)
```

## Methods

CellVoteR employs a "Mix and Match" architecture to generate votes:

| Method | Strategy                           | Scope                      |
|:-------|:-----------------------------------|:---------------------------|
| **1**  | Clustering + Fisher Test           | Triage (Immune/Endo split) |
| **2**  | Marker Thresholds                  | Triage (Immune/Endo split) |
| **3**  | Clustering (Ref Genes Only)        | Triage (Immune/Endo split) |
| **4**  | Marker Thresholds (Ref Genes Only) | Triage (Immune/Endo split) |
| **5**  | Clustering + Fisher Test           | Global (All cells)         |
| **6**  | Clustering (Ref Genes Only)        | Global (All cells)         |

*By default, `run_ensemble()` uses Methods 5 & 6. Setting `use_ensemble = TRUE` activates all 6 methods.*
