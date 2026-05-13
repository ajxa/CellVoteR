# Run the CellVoteR ensemble annotation pipeline

Orchestrates the four annotation methods and two global tie-breakers
that together form the CellVoteR ensemble. Each method runs broad
annotation, subclustering, marker ranking, panel scoring, and fine label
assignment. The tie-breakers skip the broad annotation step, operating
directly on the pre-existing unsupervised clusters from
[`prepare_sce`](prepare_sce.md).

## Usage

``` r
run_cellvoter(sce, return_full_output = FALSE, annotation_args = list())
```

## Arguments

- sce:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  processed by [`prepare_sce`](prepare_sce.md). Must have a `logcounts`
  assay, `metadata$marker_config`, `metadata$filtered_fine_markers`, and
  a `"user_panel"` `altExp`.

- return_full_output:

  Logical scalar. When `FALSE` (default), only the per-cell label
  factors are returned under `$labels`. When `TRUE`, the full output
  from each `.run_fine_annotation` call - including per-cluster score
  tables - is also returned under `$full_output`.

- annotation_args:

  Named list of argument sublists passed through to the internal
  annotation pipeline. Valid sublists are:

  `rank_args`

  :   Arguments for [`rank_cluster_markers`](rank_cluster_markers.md):
      `assay_type`, `test_type`, `direction`, `pval_type`, `min_prop`,
      `BPPARAM`.

  `extract_args`

  :   Arguments for [`extract_top_markers`](extract_top_markers.md):
      `fdr_threshold`, `effect_threshold`, `target_n`.

  `score_args`

  :   Arguments for
      [`score_markers_against_panel`](score_markers_against_panel.md).
      No additional arguments currently accepted.

  `assign_args`

  :   Arguments for [`assign_fine_labels`](assign_fine_labels.md). No
      additional arguments currently accepted.

  Only the sublists that differ from defaults need to be supplied.

## Value

A named list with the following elements:

- `sce`:

  The input SCE with broad and subcluster label columns added to
  `colData` for all four methods, and intermediate results from the
  reduced altExp methods stored in `metadata`:

  `metadata(sce)$broad_cluster_markers_reduced`

  :   Ranked marker tables from the reduced cluster-based method (method
      2).

  `metadata(sce)$broad_cell_enrichment_reduced`

  :   Cell enrichment scores and parameters from the reduced
      enrichment-based method (method 4).

- `labels`:

  Named list of six per-cell label factors ready to pass to
  [`resolve_consensus_labels`](resolve_consensus_labels.md) as
  `label_list`. Names: `method_1`, `method_2`, `method_3`, `method_4`,
  `global_1`, `global_2`.

- `method_names`:

  Character vector of the four primary method names.

- `tie_breaker_names`:

  Character vector of the two tie-breaker names.

- `full_output`:

  Only present when `return_full_output = TRUE`. Named list mirroring
  `$labels` but containing the complete
  [`assign_fine_labels`](assign_fine_labels.md) output (labels +
  per-cluster score table) for each method.

## Details

The returned label list is designed to be passed directly to
[`resolve_consensus_labels`](resolve_consensus_labels.md), which the
user calls independently so that voting parameters can be adjusted and
re-run without repeating the annotation pipeline.

## Pipeline structure


    prepare_sce()  <- must be run before this function
         |
         |- Method 1: annotate_broad_clusters() on full SCE
         |- Method 2: annotate_broad_clusters() on reduced altExp
         |- Method 3: annotate_broad_cells()    on full SCE
         |- Method 4: annotate_broad_cells()    on reduced altExp
         |       Each: broad label -> subcluster -> rank markers
         |             -> score panel -> assign fine labels
         |
         |- Tie-breaker 1: HVG clusters on full SCE  (no broad step)
         |_ Tie-breaker 2: panel clusters on reduced (no broad step)

## Ranking arguments — broad vs fine annotation

Marker ranking via [`rank_cluster_markers`](rank_cluster_markers.md) is
used at two distinct points in the pipeline, and `annotation_args`
exposes independent control over each:

- `broad_args`:

  Controls ranking inside
  [`annotate_broad_clusters`](annotate_broad_clusters.md) (methods 1 and
  2 only). These ranks are used to assign broad cell lineage labels
  (e.g. immune, vasculature) based on the median rank of curated broad
  markers. Because broad marker sets are small and highly specific, a
  more lenient `min_prop` is often appropriate here.

- `rank_args`:

  Controls ranking inside `.run_fine_annotation` for all six methods.
  These ranks are used to extract top marker genes per subcluster, which
  are then scored against the fine marker panel via Fisher's exact test.
  Methods 3 and 4 use `rank_args` only, as
  [`annotate_broad_cells`](annotate_broad_cells.md) does not call
  [`rank_cluster_markers`](rank_cluster_markers.md).

If you want both steps to behave identically, supply the same values to
both sublists. If you want to differentiate them — for example using a
lenient `min_prop` for broad assignment but a strict one for fine
annotation — supply different values:


    results <- run_cellvoter(
      sce,
      annotation_args = list(
        broad_args   = list(test_type = "t", min_prop = 0.1),
        rank_args    = list(test_type = "t", min_prop = 0.25),
        extract_args = list(fdr_threshold = 0.01)
      )
    )

## Usage


    # Default run
    results <- run_cellvoter(sce)

    # With custom annotation parameters
    results <- run_cellvoter(
      sce,
      annotation_args = list(
        rank_args    = list(test_type = "t", min_prop = 0.1),
        extract_args = list(fdr_threshold = 0.01, target_n = 50L)
      )
    )

    # Resolve consensus independently
    consensus <- resolve_consensus_labels(
      label_list        = results$labels,
      method_names      = results$method_names,
      tie_breaker_names = results$tie_breaker_names,
      allow_even_split  = TRUE,
      unassigned_label  = "unknown"
    )

    sce$cellVoteR_label  <- consensus$label
    sce$cellVoteR_method <- consensus$method

## See also

[`prepare_sce`](prepare_sce.md),
[`resolve_consensus_labels`](resolve_consensus_labels.md)
