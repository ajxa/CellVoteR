# Introduction to CellVoteR

## Overview

CellVoteR is an ensemble-based pipeline for robust cell type annotation
in single-cell RNA-seq data. Rather than assigning labels from a single
marker score or clustering strategy, CellVoteR runs complementary
annotation methods across two feature spaces and resolves the results
through a configurable consensus step. The workflow has four main
stages:

1.  Prepare a SingleCellExperiment object from raw counts.
2.  Load and configure broad and fine marker sets.
3.  Build full and marker-reduced analysis tracks.
4.  Run ensemble annotation and resolve final consensus labels.

## Installation

``` r

# install.packages("devtools")
devtools::install_github("ajxa/CellVoteR")
```

``` r

library(CellVoteR)
```

------------------------------------------------------------------------

## Complete Workflow

``` r

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
table(sce$cellVoteR_method)
```

## What Happens Internally?

At a high level, CellVoteR compares four primary annotation methods and
two global tie-breakers. The four primary methods combine:

- two feature spaces: full highly variable gene space and reduced
  marker-panel space
- two broad annotation strategies: cluster-based and enrichment-based

The consensus step then assigns a final label when the methods agree
strongly, or uses the global tie-breakers to resolve partial
disagreement.

## Where to Go Next?

For details on input file formats and marker structure, see
vignette(“user_inputs”). For preprocessing, analysis tracks, and
prepare_sce(), see the preprocessing article. For the annotation methods
and custom parameters, see the annotation methods article. For consensus
voting, output inspection, and troubleshooting, see the consensus and
results articles.

## Session Info

    #> R version 4.4.1 (2024-06-14)
    #> Platform: aarch64-apple-darwin20
    #> Running under: macOS 26.5.2
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
