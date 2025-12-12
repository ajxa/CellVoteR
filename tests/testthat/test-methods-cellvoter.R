library(testthat)
library(Seurat)
library(Matrix)

# HELPER: SYNTHETHIC DATA ------------------------------------------------------
# We need enough cells (>10) to pass the sub-clustering minimums
make_test_seurat <- function() {
  # 20 genes, 30 cells
  # Cells 1-15: Group A (High GeneA)
  # Cells 16-30: Group B (High GeneB)
  counts <- Matrix::sparseMatrix(
    i = c(rep(1, 15), rep(2, 15)),
    j = 1:30,
    x = rep(10, 30),
    dims = c(20, 30)
  )
  rownames(counts) <- c("GeneA", "GeneB", paste0("Noise", 1:18))
  colnames(counts) <- paste0("Cell_", 1:30)

  obj <- CreateSeuratObject(counts)
  obj <- NormalizeData(obj, verbose = FALSE) # Required for some internal steps
  return(obj)
}

# TESTS FOR RUN_CELLVOTER ------------------------------------------------------

test_that("run_cellvoter runs Method 5 (Global) correctly", {
  obj <- make_test_seurat()

  # Define markers that match our data
  markers <- list(TypeA = "GeneA", TypeB = "GeneB")

  # Run Method 5 (Global Clustering + All Genes)
  # Suppress warnings (e.g. low resolution clustering on small data)
  suppressWarnings({
    res <- run_cellvoter(obj, method = 5, markers = markers, verbose = FALSE)
  })

  # 1. Check Output Structure
  expect_type(res, "list")
  expect_named(res, c("object", "heatmaps", "scores"))
  expect_s4_class(res$object, "Seurat")

  # 2. Check Metadata Creation
  expect_true("CellVoteR_Label" %in% colnames(res$object@meta.data))
  expect_true("Broad_Type" %in% colnames(res$object@meta.data))

  # Method 5 sets Broad_Type to "Global"
  expect_equal(unique(res$object$Broad_Type), "Global")
})

test_that("run_cellvoter Method 6 filters genes", {
  obj <- make_test_seurat() # Has 20 genes

  # Marker list has only 2 genes
  markers <- list(TypeA = "GeneA", TypeB = "GeneB")

  # Run Method 6 (Global + Filtered Genes)
  # Note: restrict_to_reference might fail if <10 genes remain.
  # We expect an error here because our marker list is too small (2 genes vs min 10)
  # This verifies the logic IS trying to filter the genes.

  expect_error(
    run_cellvoter(obj, method = 6, markers = markers, verbose = FALSE),
    "Filtering to reference left <10 genes"
  )
})

test_that("run_cellvoter Method 2 (Triage by Cell) assigns broad types", {
  obj <- make_test_seurat()

  # Setup Triage Markers so Cells 1-15 -> Immune, 16-30 -> Endothelial
  # GeneA matches first half, GeneB matches second half
  triage <- list(Immune = "GeneA", Endothelial = "GeneB")
  markers <- list(Subtype1 = "GeneA")

  suppressWarnings({
    res <- run_cellvoter(obj, method = 2,
                         markers = markers,
                         triage_markers = triage,
                         verbose = FALSE)
  })

  # Check Broad Labels
  meta <- res$object@meta.data

  # Cell 1 should be Immune (GeneA high)
  expect_equal(as.character(meta[1, "Broad_Type"]), "Immune")

  # Cell 16 should be Endothelial (GeneB high)
  expect_equal(as.character(meta[16, "Broad_Type"]), "Endothelial")
})

test_that("run_cellvoter handles inputs with no valid markers gracefully", {
  obj <- make_test_seurat()

  # Provide markers that don't exist in the data
  markers <- list(TypeX = "GeneZ")

  suppressWarnings({
    res <- run_cellvoter(obj, method = 5, markers = markers, verbose = FALSE)
  })

  # Should run but assign "Unassigned" or similar, not crash
  expect_s4_class(res$object, "Seurat")

  # If no markers match, scores list might be empty or valid_ref is empty
  # Our logic sets label to Broad_Type (e.g. "Global") if reference is empty
  expect_true(all(res$object$CellVoteR_Label == "Global"))
})

# END --------------------------------------------------------------------------
