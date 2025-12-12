library(testthat)

# TESTS FOR PROCESS_MARKER_INPUT -----------------------------------------------

test_that("process_marker_input handles NULL (Defaults)", {
  # This relies on cellvoter_data being available (LazyData)
  skip_if_not(exists("cellvoter_data"))

  markers <- process_marker_input(NULL)
  expect_s3_class(markers, "data.frame")
  expect_true(all(c("gene", "cell_type", "cell_category") %in% colnames(markers)))
})

test_that("process_marker_input handles Named Lists", {
  input_list <- list(
    T_Cell = c("CD3D", "CD3E"),
    B_Cell = c("CD79A")
  )

  markers <- process_marker_input(input_list)

  expect_equal(nrow(markers), 3) # 2 T-cells + 1 B-cell
  expect_equal(markers$cell_type[3], "B_Cell")
  # Default category should match type
  expect_equal(markers$category[1], "T_Cell")
})

test_that("process_marker_input handles Dataframes with custom columns", {
  input_df <- data.frame(
    Symbol = c("GeneA", "GeneB"),
    Cluster = c("Type1", "Type2"),
    Broad_Type = c("Cat1", "Cat1"),
    stringsAsFactors = FALSE
  )

  markers <- process_marker_input(input_df, gene_col = "Symbol", type_col = "Cluster")

  expect_true("gene" %in% colnames(markers))
  expect_equal(markers$gene[1], "GeneA")
  expect_equal(markers$category[1], "Cat1") # Should auto-detect "Broad_Type"
})

test_that("process_marker_input throws errors for bad input", {
  expect_error(process_marker_input(list("GeneA")), "must be named")

  bad_df <- data.frame(A = 1, B = 2)
  expect_error(process_marker_input(bad_df), "Could not find a gene")
})

# TESTS FOR PROCESS_TRIAGE_MARKERS ---------------------------------------------
test_that("process_triage_input handles NULL (Defaults)", {
  # Should find internal data if available
  skip_if_not(exists("cellvoter_data"))

  res <- process_triage_input(NULL)
  expect_type(res, "list")
  # Verify it pulled the correct default keys
  expect_true("Immune" %in% names(res))
  expect_true("Endothelial" %in% names(res))
})

test_that("process_triage_input accepts user overrides", {
  # User supplies their own list
  my_triage <- list(GroupA = "Gene1")
  res <- process_triage_input(my_triage)

  expect_equal(res, my_triage)
  expect_equal(names(res), "GroupA")
})

test_that("process_triage_input validates input type", {
  expect_error(process_triage_input("Not a list"), "must be a named list")
})

# END --------------------------------------------------------------------------
