foo = subcluster_labels(
  sce = sce,
  group_col = "broad_enrichment",
  out_col = "broad_enrichment_sub"
  )



group_col = "broad_enrichment"
out_col = "broad_enrichment_sub"

feature_mode = "hvg"
subcluster_col_fmt = "%s_sc%s"
hvg_prop = 0.1
min_ncells = 50
seed = 1234
BPPARAM = BiocParallel::SerialParam()


rm(subcluster_col_fmt, hvg_prop, min_ncells, seed, BPPARAM, feature_mode,
   est_params, subset_row, vasc_sce, dec, hvgs, n_features, max_rank,
   group_col, out_col)

table(as.character(SummarizedExperiment::colData(sce)[[group_col]]))
unique(as.character(SummarizedExperiment::colData(sce)[[group_col]]))


vasc_sce <- subset(sce, select = broad_enrichment == "vasculature")
vasc_sce


est_params <- estimate_cluster_params(n_cells = ncol(vasc_sce))
subset_row <- NULL

dec <- scran::modelGeneVar(vasc_sce)
hvgs <- scran::getTopHVGs(dec, prop = hvg_prop)
subset_row <- hvgs
n_features <- length(hvgs)

max_rank <- min(est_params$n_pcs, n_features - 1L, ncol(vasc_sce) - 1L)


if (max_rank < 1L) {
  cli::cli_abort(
    c(
      "Cannot run PCA for a subcluster subset.",
      "i" = "Need at least 2 cells and 2 features."
    )
  )
}

cluster_input <- scran::fixedPCA(
  x = ,
  subset.row = subset_row,
  name = "sub_PCA",
  rank = max_rank,
  BPPARAM = BPPARAM
)
