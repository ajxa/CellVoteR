# 1.) Load and format the markers ----------------------------------------------

markers <- load_markers(file_path = "inst/input_markers.xlsx")

markers$broad <- build_broad_marker_config(
  marker_list = markers$broad,
  priority_order = c("vasculature", "immune"),
  default_threshold = 0.25
  )

# 2.) Load the single-cell expression data -------------------------------------
# data_path <- "~/Desktop/OneDrive - University of Leeds/adhoc/Ensemble_Cell_Typing/data"

# sc_data_path <- list(
#   cosmx = list.files(data_path, pattern = "sObj\\.rds", full.names = TRUE)
  # wang = list.files(file.path(data_path, "Wang_patient_seurats"), full.names = TRUE),
  # nomura =  list.files(file.path(data_path, "Nomura_patient_seurats"), full.names = TRUE)
# )

# raw data is in the count data slot typically
# This is stored as a sparse matrix, so can easily be saved to an rds file
# similarly, the metadata is stored as a data frame, and can also be saved as
# an rds file by extracting the metadata slot of the Seurat object and
# saving that separately.

# seu_data <- readRDS(sc_data_path$cosmx)
# rm(data_path, sc_data_path, seu_data)
#

# cell_meta <- readRDS("~/Desktop/example_data/cell_metadata.rds")
# data <- readRDS("~/Desktop/example_data/test_input_sc_matrix.rds")
#
# cell_meta[1:5,1:7]
# data[1:5,1:7]

# Loading a sparse matrix and cell_metadata saved as .rds objects
sce <- create_sce(
  counts = "~/Desktop/example_data/test_input_sc_matrix.rds",
  cell_metadata = "~/Desktop/example_data/cell_metadata.rds"
)

# Loading a sparse matrix saved as a .mtx file, with separate files for genes and cells
# mtx_data <- create_sce(
#   mtx_file = "inst/test_input_sc_matrix.mtx",
#   genes_file = "inst/genes.txt",
#   cells_file = "inst/cells.txt",
#   cell_metadata = cell_meta
#   )

# 3.) QC & Pre-process ----------------------------------------------------------
sce <- assess_cell_quality(sce, remove_failed_cells = TRUE)

sce <- normalize_counts(sce)

sce <- prepare_sce(sce, markers)
rm(markers)

reduced <- SingleCellExperiment::altExp(sce)

# The reduced feature set is stored in the `altExperiment` slot of the SingleCellExperiment class
# reduced <- SingleCellExperiment::altExp(sce)


clusters <- list(all = list(), reduced = list())
cells    <- list(all = list(), reduced = list())
ties     <- list(global_1 = list(), global_2 = list())

# (All) Cluster-based Labelling ------------------------------------------------

sce <- annotate_broad_clusters(sce, label_col = "broad_cluster")
sce <- subcluster_labels(sce, group_col = "broad_cluster", out_col = "broad_cluster_sub")

clusters$all$marker_ranks <- rank_cluster_markers(sce = sce,
                                                  cluster_col = "broad_cluster_sub",
                                                  return_list = TRUE
                                                  )

clusters$all$top_markers <- extract_top_markers(ranked_markers = clusters$all$marker_ranks)

clusters$all$label_scores <- score_markers_against_panel(
  top_markers = clusters$all$top_markers,
  marker_panel = sce@metadata$filtered_fine_markers,
  background_genes = rownames(sce)
)

clusters$all$final <- assign_fine_labels(cluster_col = sce$broad_cluster_sub,
                                         scores = clusters$all$label_scores)

# (Reduced) Cluster-based Labelling --------------------------------------------

reduced <- annotate_broad_clusters(sce = reduced,
                                   broad_config = reduced@metadata$marker_config$broad,
                                   cluster_col = "cluster",
                                   label_col = "broad_cluster"
                                   )

reduced <- subcluster_labels(sce = reduced,
                             group_col = "broad_cluster",
                             feature_mode = "all",
                             out_col = "broad_cluster_sub"
                             )


clusters$reduced$marker_ranks <- rank_cluster_markers(sce = reduced,
                                                      cluster_col = "broad_cluster_sub",
                                                      return_list = TRUE
                                                      )

clusters$reduced$top_markers <- extract_top_markers(ranked_markers = clusters$reduced$marker_ranks)

clusters$reduced$label_scores <- score_markers_against_panel(
  top_markers = clusters$reduced$top_markers,
  marker_panel = reduced@metadata$filtered_fine,
  background_genes = rownames(reduced)
)

clusters$reduced$final <- assign_fine_labels(cluster_col = reduced$broad_cluster_sub,
                                             scores = clusters$reduced$label_scores
                                             )

# (All) Cell Enrichment-based Labelling ----------------------------------------

sce <- annotate_broad_cells(sce, label_col = "broad_enrichment")
sce <- subcluster_labels(sce, group_col = "broad_enrichment", out_col = "broad_enrichment_sub")

cells$all$marker_ranks <- rank_cluster_markers(sce = sce,
                                               cluster_col = "broad_enrichment_sub",
                                               return_list = TRUE
                                               )

cells$all$top_markers <- extract_top_markers(ranked_markers = cells$all$marker_ranks)

cells$all$label_scores <- score_markers_against_panel(
  top_markers = cells$all$top_markers,
  marker_panel = sce@metadata$filtered_fine_markers,
  background_genes = rownames(sce)
)

cells$all$final <- assign_fine_labels(cluster_col = sce$broad_enrichment_sub,
                                      scores = cells$all$label_scores
                                      )

# (Reduced) Cell Enrichment-based Labelling ------------------------------------

reduced <- annotate_broad_cells(sce = reduced, label_col = "broad_enrichment")


reduced <- subcluster_labels(sce = reduced,
                             group_col = "broad_enrichment",
                             feature_mode = "all",
                             out_col = "broad_enrichment_sub"
                             )

cells$reduced$marker_ranks <- rank_cluster_markers(sce = reduced,
                                                   cluster_col = "broad_enrichment_sub",
                                                   return_list = TRUE
                                                   )

cells$reduced$top_markers <- extract_top_markers(ranked_markers = cells$reduced$marker_ranks)

cells$reduced$label_scores <- score_markers_against_panel(
  top_markers = cells$reduced$top_markers,
  marker_panel = reduced@metadata$filtered_fine,
  background_genes = rownames(reduced)
)

cells$reduced$final <- assign_fine_labels(cluster_col = reduced$broad_enrichment_sub,
                                          scores = cells$reduced$label_scores
                                          )

# Global 1 tie-breaker ----------------------------------------------------------

ties$global_1$marker_ranks <- rank_cluster_markers(sce = sce,
                                                   cluster_col = "cluster_broad_hvg",
                                                   return_list = TRUE
                                                   )

ties$global_1$top_markers <- extract_top_markers(ranked_markers = ties$global_1$marker_ranks)

ties$global_1$label_scores <- score_markers_against_panel(
  top_markers = ties$global_1$top_markers,
  marker_panel = sce@metadata$filtered_fine_markers,
  background_genes = rownames(sce)
)

ties$global_1$final <- assign_fine_labels(cluster_col = sce$cluster_broad_hvg,
                                          scores = ties$global_1$label_scores
                                          )

# Global 2 tie-breaker ---------------------------------------------------------

ties$global_2$marker_ranks <- rank_cluster_markers(sce = reduced,
                                                   cluster_col = "cluster",
                                                   return_list = TRUE
                                                   )

ties$global_2$top_markers <- extract_top_markers(ranked_markers = ties$global_2$marker_ranks)


ties$global_2$label_scores <- score_markers_against_panel(
  top_markers = ties$global_2$top_markers,
  marker_panel = reduced@metadata$filtered_fine,
  background_genes = rownames(reduced)
)

ties$global_2$final <- assign_fine_labels(cluster_col = reduced$cluster,
                                          scores = ties$global_2$label_scores
                                          )

# Label comparisons ------------------------------------------------------------
table(clusters$all$final$labels)
table(clusters$reduced$final$labels)

table(cells$all$final$labels)
table(cells$reduced$final$labels)

table(ties$global_1$final$labels)
table(ties$global_2$final$labels)
# Reconciling final labels -----------------------------------------------------

out_labs <- list(
  method_1 = clusters$all$final$labels,
  method_2 = clusters$reduced$final$labels,
  method_3 = cells$all$final$labels,
  method_4 = cells$reduced$final$labels,
  global_1 = ties$global_1$final$labels,
  global_2 = ties$global_2$final$labels
  )


foo = resolve_consensus_labels(
  label_list = out_labs,
  method_names = names(out_labs)[1:4],
  tie_breaker_names = names(out_labs)[5:6],
  allow_even_split = TRUE,
  unassigned_label = "unknown"
  )

sce$cellVoteR_label <- foo$label
sce$cellVoteR_method <- foo$method

# End --------------------------------------------------------------------------
