# CellVoteR <img src="man/figures/logo.png" align="right" height="130" alt="" />

An ensemble-based pipeline for robust cell type classification in single-cell RNAseq data. 

It moves beyond simple "best-match" scoring by integrating multiple classification 
strategies: 
lineage triage, targeted sub-clustering, and global consensus, to generate a 
high-confidence label for every cell.

## How It Works

CellVoteR determines cell identity through a three-stage voting process:

### 1. Broad Triage (divide & conquer)
The dataset is first split into broad biological categories (e.g., *Immune* vs. *Vascular* vs. *Tumour*) 
using defined marker thresholds or coarse clustering. 
This "Triage" step isolates distinct lineages, 
preventing dominant signals from obscuring rare cell types.

### 2. Targeted Sub-clustering
Each broad category is processed independently. Cells are sub-clustered to 
identify fine-grained states, and their identity is determined by 
performing Fisher's Exact Tests against a reference marker panel. 
This ensures that an immune cell is only compared against immune markers, 
reducing false positives.

### 3. Global Consensus (tie-breaker)
Simultaneously, the full dataset is clustered globally (without triage). 
These global labels serve as a "baseline vote." 
In the final ensemble step, if the triage methods disagree, 
the global result and raw marker intensities are used to break the tie and 
resolve the final identity.

## Workflow Diagram

```mermaid
graph TD
    %% --- Node Definitions ---
    Input([Input Data]):::start
    QC[Quality Control]:::proc
    
    %% Track 1: Triage
    subgraph T1 [Track 1: Broad Triage]
        direction TB
        Triage{Broad Triage}:::decision
        
        %% Paths
        subgraph P1 [ ]
            direction TB
            PathA[Immune Subset]:::proc
            SubA[Sub-Cluster &<br/>Fisher Score]:::proc
        end
        
        subgraph P2 [ ]
            direction TB
            PathB[Vascular Subset]:::proc
            SubB[Sub-Cluster &<br/>Fisher Score]:::proc
        end
        
        subgraph P3 [ ]
            direction TB
            PathC[Remaining Cells]:::proc
            SubC[Sub-Cluster &<br/>Fisher Score]:::proc
        end
    end

    %% Track 2: Global
    subgraph T2 [Track 2: Global Consensus]
        direction TB
        Global[Global Clustering]:::proc
        ScoreG[Global Fisher Score]:::proc
    end

    %% --- Logic Flow ---
    Input --> QC
    QC --> Triage
    QC --> Global
    
    Triage -- Immune --> PathA --> SubA
    Triage -- Endo --> PathB --> SubB
    Triage -- Other --> PathC --> SubC
    
    Global --> ScoreG
    
    %% --- Voting & Resolution ---
    SubA & SubB & SubC --> Vote[Ensemble Voting]:::proc
    ScoreG --> Vote
    
    Vote --> Resolve{Clash?}:::decision
    Resolve -- Yes --> TieBreak[Apply Global Breaker<br/>& Marker Intensity]:::proc
    Resolve -- No --> Final([Final Cell Label]):::endnode
    TieBreak --> Final

    %% --- Styling ---
    %% Theme-Neutral Colors: Light fills with black text ensure readability in Dark Mode
    classDef start fill:#e1f5fe,stroke:#01579b,stroke-width:2px,color:#000
    classDef proc fill:#f5f5f5,stroke:#333,stroke-width:1px,color:#000
    classDef decision fill:#fff9c4,stroke:#fbc02d,stroke-width:1px,stroke-dasharray: 5 5,color:#000
    classDef endnode fill:#dcedc8,stroke:#33691e,stroke-width:2px,color:#000
    
    linkStyle default stroke:#b0bec5,stroke-width:2px
```

## Installation

You can install the development version of CellVoteR with:

``` r
# install.packages("devtools")
devtools::install_github("ajxa/CellVoteR")
```

## Quick Start

### 1. Quality Control

Before labelling, ensure your data is clean using the built-in QC functions.

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
