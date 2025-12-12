library(testthat)
library(Matrix)
library(Seurat)

# HELPER: MAKE A TINY SEURAT OBJECT --------------------------------------------
make_small_seurat <- function() {
  counts <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(1, 2), x = c(10, 10), dims = c(10, 2)
  )
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- c("Cell1", "Cell2")
  CreateSeuratObject(counts)
}

# TESTS FOR BROAD_LABEL_BY_CELL ------------------------------------------------
test_that("broad_label_by_cell works correctly", {
  obj <- make_small_seurat()

  # Gene1 is high in Cell1, Gene2 is high in Cell2
  markers <- list(GroupA = "Gene1", GroupB = "Gene2")

  labels <- broad_label_by_cell(obj, marker_list = markers, expr_threshold = 1)

  expect_equal(labels[1], "GroupA")
  expect_equal(labels[2], "GroupB")
})

test_that("broad_label_by_cell handles unmatched cells", {
  obj <- make_small_seurat()
  markers <- list(GroupA = "Gene99") # Gene doesn't exist/express

  labels <- broad_label_by_cell(obj, marker_list = markers)
  expect_true(all(labels == "Other"))
})

# END --------------------------------------------------------------------------
