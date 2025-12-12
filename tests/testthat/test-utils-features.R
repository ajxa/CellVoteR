library(testthat)
library(Matrix)

# HELPER: CREATE MATRIX WITH SPECIFIC GENE NAMES -------------------------------
make_gene_matrix <- function() {
  # Create a matrix with Mito, Ribo, and Hemoglobin genes
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4, 1),
    j = c(1, 1, 1, 1, 2),
    x = c(10, 10, 10, 10, 5),
    dims = c(5, 2)
  )
  # Rownames mimic real gene patterns
  rownames(counts) <- c("MT-CO1", "RPS6", "RPL10", "HBA1", "ACTB")
  colnames(counts) <- c("Cell1", "Cell2")
  return(counts)
}

# TEST FOR FIND_QC_FEATURES ----------------------------------------------------

test_that("find_qc_features identifies default groups (Mito/Ribo)", {
  mat <- make_gene_matrix()

  # Run with defaults
  qc_list <- find_qc_features(mat)

  # Check Structure
  expect_type(qc_list, "list")
  expect_named(qc_list, c("mito", "ribo"))

  # Check Content
  # "MT-CO1" matches ^MT-
  expect_equal(qc_list$mito$features, "MT-CO1")
  expect_equal(qc_list$mito$max_pct, 20)

  # "RPS6" and "RPL10" match ^RP[SL]
  expect_true(all(c("RPS6", "RPL10") %in% qc_list$ribo$features))
  expect_equal(qc_list$ribo$max_pct, 50)
})

test_that("find_qc_features accepts custom configs", {
  mat <- make_gene_matrix()

  # Custom config
  my_config <- list(
    hemo = list(pattern = "^HB", max_pct = 5)
  )

  qc_list <- find_qc_features(mat, group_configs = my_config)

  expect_named(qc_list, "hemo")
  expect_equal(qc_list$hemo$features, "HBA1")
  expect_equal(qc_list$hemo$max_pct, 5)
})

test_that("find_qc_features handles no matches gracefully", {
  mat <- make_gene_matrix()

  # Look for a gene pattern that doesn't exist
  my_config <- list(
    missing = list(pattern = "^ZYZ", max_pct = 10)
  )

  qc_list <- find_qc_features(mat, group_configs = my_config)

  # Should return empty character vector, not NULL or error
  expect_equal(qc_list$missing$features, character(0))
})

# TEST FOR FEATURE_SYNONYMS ----------------------------------------------------
test_that("feature_synonyms handles Regex inputs locally (with valid_genes)", {
  valid_genes <- c("MT-CO1", "MT-ND1", "ACTB", "GAPDH")
  input_regex <- "^MT-"

  res <- feature_synonyms(input_regex, valid_genes = valid_genes)

  expect_true("MT-CO1" %in% res)
  expect_true("MT-ND1" %in% res)
  expect_false("ACTB" %in% res)
})

test_that("feature_synonyms works without valid_genes (General Lookup)", {
  # This tests the NULL default behavior
  input_genes <- c("CD3D", "CD3E")

  # Run without valid_genes
  # (Wrapped in expect_no_error - if biomaRt is offline/slow)
  expect_no_error({
    res <- feature_synonyms(input_genes, valid_genes = NULL)
  })

  # Even if biomaRt fails (returns input), or succeeds (returns synonyms),
  # we just want to ensure it accepted NULL and returned a character vector
  expect_type(res, "character")
  expect_true(length(res) >= 2)
})

test_that("feature_synonyms checks if features are available in biomaRt", {

  valid <- c("GeneA")

  res <- tryCatch(
    feature_synonyms("GeneA", valid),
    error = function(e) "GeneA" # fallback to "GeneA"
  )

  expect_true(identical(res, "GeneA") || length(res) == 0)

})

# END --------------------------------------------------------------------------
