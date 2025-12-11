library(testthat)
library(Matrix)

# --- Helper ---
make_syn_data <- function() {
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4,   1, 3,   3,    2),
    j = c(1, 1, 1, 1,   2, 2,   3,    4),
    x = c(10, 10, 40, 40, 50, 5,  10,   50),
    dims = c(4, 4)
  )
  colnames(counts) <- paste0("Cell_", 1:4)
  rownames(counts) <- c("MT-1", "RPS-1", "GeneA", "GeneB")
  return(counts)
}

test_that("Core QC accepts explicit gene lists", {
  counts <- make_syn_data()

  # User (or helper) supplies EXPLICIT genes
  groups <- list(
    mito = list(features = c("MT-1"), max_pct = 20),
    ribo = list(features = c("RPS-1"), max_pct = 50)
  )

  res <- assess_cell_quality(counts, check_feature_groups = groups, min_features = 0, min_cells_per_sample = 1)

  # Check Metrics
  # Cell 1: MT-1 is 10/100 = 10% (Passes <20)
  expect_equal(res["Cell_1", "pct_mito"], 10)
  expect_true(res["Cell_1", "QC_PASS"])

  # Cell 2: MT-1 is 50/55 = 90% (Fails >20)
  expect_equal(round(res["Cell_2", "pct_mito"]), 91)
  expect_false(res["Cell_2", "QC_PASS"])
})

test_that("Core QC warns on missing features (non-empty matrix)", {
  counts <- make_syn_data()

  # Supply a gene that doesn't exist
  bad_groups <- list(
    missing_group = list(features = c("NonExistentGene"), max_pct = 10)
  )

  expect_warning(
    assess_cell_quality(counts, check_feature_groups = bad_groups, min_features = 0, min_cells_per_sample = 1),
    "missing_group"
  )
})

test_that("Core QC works with NO feature groups (basic mode)", {
  counts <- make_syn_data()

  # No groups provided -> should only check min_features
  res <- assess_cell_quality(counts, check_feature_groups = list(), min_features = 200, min_cells_per_sample = 1)

  # Cell 1 has 4 features, min is 200 -> Fail
  expect_false(res["Cell_1", "QC_PASS"])
  expect_true("n_features" %in% colnames(res))
  # Should NOT have pct_mito columns
  expect_false("pct_mito" %in% colnames(res))
})

test_that("Helper function resolves regex correctly", {
  counts <- make_syn_data() # Has "MT-1" and "RPS-1"

  # Define configs with regex
  regex_configs <- list(
    my_mito = list(pattern = "^MT-", max_pct = 15)
  )

  # Run the helper
  resolved <- find_qc_features(counts, group_configs = regex_configs)

  # Check output structure
  expect_type(resolved, "list")
  expect_equal(resolved$my_mito$features, "MT-1")
  expect_equal(resolved$my_mito$max_pct, 15)
})


test_that("Filter function removes failed cells correctly", {
  counts <- make_syn_data() # 4 cells
  # Let's say Cell_2 and Cell_3 fail
  # We manually construct a fake QC dataframe to control the test strictly
  qc_df <- data.frame(
    cell_id = colnames(counts),
    sample_id = c("S1", "S1", "S2", "S2"),
    QC_PASS = c(TRUE, FALSE, FALSE, TRUE), # Keep Cell 1 and 4
    stringsAsFactors = FALSE
  )

  # 1. Test Matrix Subsetting
  clean_mat <- filter_cells(counts, qc_df, verbose = FALSE)
  expect_equal(ncol(clean_mat), 2)
  expect_equal(colnames(clean_mat), c("Cell_1", "Cell_4"))

  # 2. Test Verbose Output (Capture the message)
  expect_message(
    filter_cells(counts, qc_df, verbose = TRUE),
    "Global Filter Summary"
  )

  # 3. Test Sample Breakdown output
  expect_message(
    filter_cells(counts, qc_df, verbose = TRUE, show_by_sample = TRUE),
    "Breakdown by Sample"
  )
})

test_that("Filter function catches mismatch errors", {
  counts <- make_syn_data()


  qc_df_short <- data.frame(
    cell_id = c("Cell_1", "Cell_2", "Cell_3"),
    QC_PASS = c(TRUE, TRUE, TRUE)
  )

  expect_error(
    filter_cells(counts, qc_df_short),
    "Dimension mismatch"
  )
})

test_that("Filter function handles Seurat objects", {
  skip_if_not_installed("Seurat")

  counts <- make_syn_data()
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts)

  qc_df <- data.frame(
    cell_id = colnames(counts),
    QC_PASS = c(TRUE, FALSE, TRUE, FALSE) # Keep 2
  )

  clean_seurat <- filter_cells(seurat_obj, qc_df, verbose = FALSE)

  expect_s4_class(clean_seurat, "Seurat")
  expect_equal(ncol(clean_seurat), 2)
})

test_that("Filter function handles SingleCellExperiment objects", {
  skip_if_not_installed("SingleCellExperiment")

  counts <- make_syn_data()
  sce_obj <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

  qc_df <- data.frame(
    cell_id = colnames(counts),
    QC_PASS = c(TRUE, FALSE, TRUE, FALSE) # Keep 2
  )

  clean_sce <- filter_cells(sce_obj, qc_df, verbose = FALSE)

  expect_s4_class(clean_sce, "SingleCellExperiment")
  expect_equal(ncol(clean_sce), 2)
})
