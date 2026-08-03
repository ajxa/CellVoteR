# Consensus and Results

## Overview

This article describes how CellVoteR resolves method-level annotations
into final cell labels.

## Resolve Consensus Labels

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

## Attach Labels

``` r

sce$cellVoteR_label  <- consensus$label
sce$cellVoteR_method <- consensus$method
```

## Voting Parameters

| Parameter | Default | Effect |
|----|----|----|
| `allow_even_split` | `FALSE` | When `TRUE`, a 2-of-4 plurality can be accepted |
| `ordered_tiebreak` | `TRUE` | When `TRUE`, tie-breakers are checked in priority order |
| `unassigned_label` | `"Unknown"` | Label assigned when consensus cannot be resolved |

## Inspect Results

``` r

table(sce$cellVoteR_label)
table(sce$cellVoteR_method)
table(sce$cellVoteR_label, sce$cellVoteR_method)
```

## Compare Method Outputs

``` r

table(results$labels$method_1)
table(results$labels$method_2)
table(results$labels$method_3)
table(results$labels$method_4)
table(results$labels$global_1)
table(results$labels$global_2)
```

## Method Agreement

``` r

label_df <- as.data.frame(results$labels[1:4])

mean(apply(label_df, 1, function(x) length(unique(x)) == 1L))

label_df$n_agree <- apply(label_df, 1, function(x) max(table(x)))
table(label_df$n_agree)
```

## Re-run Consensus

``` r

consensus_liberal <- resolve_consensus_labels(
  label_list        = results$labels,
  method_names      = results$method_names,
  tie_breaker_names = results$tie_breaker_names,
  unassigned_label  = "unknown",
  allow_even_split  = TRUE
)
```
