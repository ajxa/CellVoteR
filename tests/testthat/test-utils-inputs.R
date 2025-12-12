library(testthat)
library(Matrix)

# HELPER: CREATE SYNTHETIC DATA ------------------------------------------------
make_sparse_matrix <- function() {
  Matrix::sparseMatrix(
    i = c(1, 2, 1), j = c(1, 2, 2), x = c(10, 20, 5),
    dims = c(3, 2),
    dimnames = list(c("GeneA", "GeneB", "GeneC"), c("Cell1", "Cell2"))
  )
}

# TESTS FOR EXTRACT_COUNTS -----------------------------------------------------
test_that("extract_counts handles Matrix objects (Sparse & Dense)", {
  # 1. Sparse Input
  sparse_mat <- make_sparse_matrix()
  res_sparse <- extract_counts(sparse_mat)

  expect_s4_class(res_sparse, "dgCMatrix")
  expect_equal(sum(res_sparse), 35)
  expect_equal(dim(res_sparse), c(3, 2))

  # 2. Dense Input
  dense_mat <- as.matrix(sparse_mat)
  res_dense <- extract_counts(dense_mat)

  # Should be converted to sparse automatically
  expect_s4_class(res_dense, "dgCMatrix")
  expect_equal(as.matrix(res_dense), dense_mat)
})

test_that("extract_counts handles Seurat objects", {
  skip_if_not_installed("Seurat")

  sparse_mat <- make_sparse_matrix()
  seurat_obj <- Seurat::CreateSeuratObject(counts = sparse_mat)

  res <- extract_counts(seurat_obj)

  expect_s4_class(res, "dgCMatrix")
  expect_equal(dim(res), c(3, 2))
  expect_equal(res["GeneA", "Cell1"], 10)
})

test_that("extract_counts handles SingleCellExperiment objects", {
  skip_if_not_installed("SingleCellExperiment")

  sparse_mat <- make_sparse_matrix()
  sce_obj <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = sparse_mat))

  res <- extract_counts(sce_obj)

  expect_s4_class(res, "dgCMatrix")
  expect_equal(res["GeneB", "Cell2"], 20)
})

test_that("extract_counts throws error on invalid input", {
  expect_error(extract_counts("Not a matrix"), "not supported")
  expect_error(extract_counts(list(a = 1)), "not supported")
})

# TEST FOR ENSURE_SEURAT -------------------------------------------------------
test_that("ensure_seurat converts Matrix to Seurat", {
  skip_if_not_installed("Seurat")

  mat <- make_sparse_matrix()
  res <- ensure_seurat(mat)

  expect_s4_class(res, "Seurat")

  expect_equal(Seurat::GetAssayData(res, layer = "counts")["GeneA", "Cell1"], 10)
})

test_that("ensure_seurat passes existing Seurat objects through unchanged", {
  skip_if_not_installed("Seurat")

  mat <- make_sparse_matrix()
  original <- Seurat::CreateSeuratObject(counts = mat)

  # Add some metadata to ensure we don't lose it
  original$test_meta <- "keep_me"

  res <- ensure_seurat(original)

  expect_identical(res, original)

  expect_equal(unname(res$test_meta[1]), "keep_me")
})

test_that("ensure_seurat converts SingleCellExperiment to Seurat", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleCellExperiment")

  mat <- make_sparse_matrix()
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat))

  res <- ensure_seurat(sce)

  expect_s4_class(res, "Seurat")
  expect_equal(ncol(res), 2)
})

# END --------------------------------------------------------------------------
