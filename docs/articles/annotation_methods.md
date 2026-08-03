# Annotation Methods

## Overview

CellVoteR runs four primary annotation methods and two global
tie-breakers.

## Run CellVoteR

``` r

results <- run_cellvoter(sce)
```

## Primary Methods

| Method   | Feature space        | Broad strategy   | Fine annotation           |
|----------|----------------------|------------------|---------------------------|
| Method 1 | Full HVG space       | Cluster-based    | Subcluster marker scoring |
| Method 2 | Reduced marker space | Cluster-based    | Subcluster marker scoring |
| Method 3 | Full HVG space       | Enrichment-based | Subcluster marker scoring |
| Method 4 | Reduced marker space | Enrichment-based | Subcluster marker scoring |

## Global Tie-breakers

| Tie-breaker   | Feature space        | Strategy                             |
|---------------|----------------------|--------------------------------------|
| Tie-breaker 1 | Full HVG space       | Global clustering and marker scoring |
| Tie-breaker 2 | Reduced marker space | Global clustering and marker scoring |

## Fine Annotation Flow

``` text
broad annotation
↓
subcluster_labels()
↓
rank_cluster_markers()
↓
extract_top_markers()
↓
score_markers_against_panel()
↓
assign_fine_labels()
```

## Custom Parameters

``` r

results <- run_cellvoter(
  sce,
  return_full_output = TRUE,
  annotation_args = list(
    broad_args = list(
      test_type = "wilcox",
      min_prop  = 0.1
    ),
    rank_args = list(
      test_type = "wilcox",
      min_prop  = 0.25
    ),
    extract_args = list(
      fdr_threshold    = 0.05,
      effect_threshold = 0.6,
      target_n         = 100L
    )
  )
)
```

## Full Outputs

``` r

results$full_output$method_1$scores
results$full_output$global_2$scores
```
