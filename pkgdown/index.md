# CellVoteR <img src="man/figures/logo.png" align="right" height="130" alt="" />

CellVoteR is an ensemble-based pipeline for robust cell type annotation in single-cell RNA-seq data.

Rather than relying on a single best-match score, 
CellVoteR combines complementary annotation strategies across full and 
marker-reduced feature spaces, then resolves agreement and disagreement 
through a configurable consensus voting step.


<img src="man/figures/cellvoter-workflow.svg"
     alt="CellVoteR workflow"
     width="100%">

<hr>

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("ajxa/CellVoteR")
```


## Quick Start

```r
markers <- load_markers(file_path = "path/to/input_markers.xlsx")

markers$broad <- build_broad_marker_config(
  marker_list       = markers$broad,
  priority_order    = c("vasculature", "immune"),
  default_threshold = 0.25
)

sce <- create_sce(
  counts        = "path/to/counts.rds",
  cell_metadata = "path/to/metadata.rds"
)

sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)
sce <- normalize_counts(sce)
sce <- prepare_sce(sce, markers)

results <- run_cellvoter(sce)

consensus <- resolve_consensus_labels(
  label_list        = results$labels,
  method_names      = results$method_names,
  tie_breaker_names = results$tie_breaker_names,
  unassigned_label  = "unknown"
)

sce$cellVoteR_label  <- consensus$label
sce$cellVoteR_method <- consensus$method

table(sce$cellVoteR_label)

```

<div class="row">
<div class="col-md-6">

## Getting started

[Introduction to CellVoteR](articles/intro_to_CellVoteR.html)  
A short walkthrough of the complete CellVoteR workflow.

[User inputs and required formats](articles/user_inputs.html)  
How to structure count matrices, metadata, and marker files.

[Preprocessing and analysis tracks](articles/preprocessing.html)  
QC, normalisation, marker overlap checks, and feature-space construction.

</div>
<div class="col-md-6">

## Annotation workflow

[Annotation methods](articles/annotation_methods.html)  
The four primary methods and two global tie-breakers.

[Consensus and results](articles/consensus_and_results.html)  
How final labels are resolved and inspected.

</div>
</div>


<div class="row">
<div class="col-md-8">

## Practical guidance

[Troubleshooting](articles/troubleshooting.html)  
Common issues with marker overlap, unresolved labels, and large datasets.

</div>
</div>


<div class="row">
<div class="col-md-4">

## Source

[Browse source code](https://github.com/ajxa/CellVoteR)  
[Report a bug](https://github.com/ajxa/CellVoteR/issues)

</div>
<div class="col-md-4">

## License

[Full license](LICENSE.html)

</div>
<div class="col-md-4">

## Citation

[How to cite CellVoteR](authors.html#citation)

</div>
</div>
