library(testthat)

# TESTS FOR CALCULATE_FISHER_SCORES --------------------------------------------
test_that("Fisher test correctly identifies perfect overlaps", {
  # Cluster 1 matches TypeX perfectly (3 genes overlap)
  cluster_markers <- data.frame(gene = c("A", "B", "C"), cluster = "1")
  ref_markers <- data.frame(Gene = c("A", "B", "C"), Broad_Cell_Type = "TypeX")
  background <- c("A", "B", "C", paste0("Noise", 1:97))

  scores <- calculate_fisher_scores(cluster_markers, ref_markers, background)

  # Expect highly significant p-value
  expect_lt(scores$Pval[1], 0.05)
  expect_equal(scores$n_overlap[1], 3)
})

test_that("Fisher test handles no overlap", {
  cluster_markers <- data.frame(gene = c("A", "B"), cluster = "1")
  ref_markers <- data.frame(Gene = c("Y", "Z"), Broad_Cell_Type = "TypeX")
  background <- c("A", "B", "Y", "Z", "Other")

  scores <- calculate_fisher_scores(cluster_markers, ref_markers, background)

  expect_equal(scores$n_overlap[1], 0)
  expect_gte(scores$Pval[1], 0.5)
})

# TESTS FOR GET BEST MATCH -----------------------------------------------------
test_that("get_best_match picks the winner", {
  scores <- data.frame(
    Cluster = c("1", "1"),
    Reference_Type = c("Winner", "Loser"),
    Pval = c(0.001, 0.5),
    stringsAsFactors = FALSE
  )

  res <- get_best_match(scores)
  expect_equal(as.character(res["1"]), "Winner")
})

# END --------------------------------------------------------------------------
