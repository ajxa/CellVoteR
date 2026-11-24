# Small marker table: two cancer labels, two normal labels, two immune labels
toy_markers <- data.frame(
  type     = c("broad", "fine", "fine", "fine", "fine", "fine", "fine"),
  category = c("immune", "cancer", "cancer", "normal", "normal", "immune", "immune"),
  label    = c(NA, "GPC", "OPC", "Astrocyte", "OPC normal", "Macrophage", "Microglia"),
  marker   = c("PTPRC", "EGFR", "PDGFRA", "GJA1", "CSPG4", "CD163", "P2RY12"),
  stringsAsFactors = FALSE
)

methods <- paste0("method_", 1:4)

make_votes <- function(...) {
  rows <- list(...)
  out <- lapply(seq_along(methods), function(j) vapply(rows, `[`, character(1), j))
  stats::setNames(out, methods)
}

run_step <- function(votes, consensus_label, label_category) {
  consensus <- list(
    label  = factor(consensus_label),
    method = factor(ifelse(consensus_label == "Unknown", "Unresolved", "Majority vote"))
  )
  suppressMessages(
    resolve_category_ties(consensus, votes, methods, label_category)
  )
}

test_that("load_markers keeps each label's original category", {
  m <- suppressMessages(load_markers(toy_markers))
  # cancer and normal are merged into "other" in $fine ...
  expect_true("other" %in% names(m$fine))
  expect_false(any(c("cancer", "normal") %in% names(m$fine)))
  # ... but stay distinguishable in $label_category
  expect_equal(unname(m$label_category[c("GPC", "OPC", "Astrocyte", "OPC normal", "Microglia")]),
               c("cancer", "cancer", "normal", "normal", "immune"))
})

test_that("a label listed under two categories is marked ambiguous", {
  dup <- rbind(toy_markers,
               data.frame(type = "fine", category = "normal", label = "OPC",
                          marker = "OLIG2"))
  m <- suppressMessages(load_markers(dup))
  expect_true(is.na(m$label_category[["OPC"]]))
  expect_equal(m$label_category[["GPC"]], "cancer")
})

test_that("same-category ties become category-level labels; others are untouched", {
  lc <- suppressMessages(load_markers(toy_markers))$label_category

  votes <- make_votes(
    c("GPC", "GPC", "OPC", "OPC"),                        # 1 cancer tie          -> cancer_unresolved
    c("Astrocyte", "Astrocyte", "GPC", "GPC"),            # 2 normal vs cancer    -> Unknown
    c("Macrophage", "Microglia", "Macrophage", "Microglia"), # 3 immune tie       -> immune_unresolved
    c("GPC", "GPC", "GPC", "OPC"),                        # 4 already resolved    -> GPC
    c("GPC", "OPC", "Astrocyte", "OPC normal"),           # 5 2 cancer + 2 normal -> Unknown
    c("OPC normal", "Astrocyte", "Astrocyte", "OPC normal") # 6 normal tie       -> normal_unresolved
  )
  cons <- c("Unknown", "Unknown", "Unknown", "GPC", "Unknown", "Unknown")

  res <- run_step(votes, cons, lc)

  expect_equal(as.character(res$label),
               c("cancer_unresolved", "Unknown", "immune_unresolved",
                 "GPC", "Unknown", "normal_unresolved"))
  expect_equal(as.character(res$method),
               c("Same-category tie", "Unresolved", "Same-category tie",
                 "Majority vote", "Unresolved", "Same-category tie"))
  expect_equal(as.character(res$label_before_tie_step), cons)
})

test_that("missing votes, unknown labels and ambiguous labels keep the cell unassigned", {
  lc <- suppressMessages(load_markers(toy_markers))$label_category
  lc_amb <- lc; lc_amb[["OPC"]] <- NA_character_

  votes <- make_votes(
    c("GPC", "GPC", "OPC", NA),          # missing vote
    c("GPC", "GPC", "OPC", ""),          # empty vote
    c("GPC", "GPC", "NotALabel", "OPC")  # label not in the lookup
  )
  res <- run_step(votes, rep("Unknown", 3), lc)
  expect_equal(as.character(res$label), rep("Unknown", 3))

  votes2 <- make_votes(c("GPC", "GPC", "OPC", "OPC"))
  expect_equal(as.character(run_step(votes2, "Unknown", lc_amb)$label), "Unknown")
})

test_that("custom unassigned label and suffix are respected", {
  lc <- suppressMessages(load_markers(toy_markers))$label_category
  votes <- make_votes(c("GPC", "OPC", "GPC", "OPC"))
  consensus <- list(label = factor("unassigned"), method = factor("Unresolved"))
  res <- suppressMessages(resolve_category_ties(consensus, votes, methods, lc,
                                                unassigned_label = "unassigned",
                                                suffix = "_state_unknown"))
  expect_equal(as.character(res$label), "cancer_state_unknown")
})

test_that("works end to end with resolve_consensus_labels()", {
  lc <- suppressMessages(load_markers(toy_markers))$label_category
  votes <- make_votes(
    c("GPC", "GPC", "OPC", "OPC"),
    c("Astrocyte", "Astrocyte", "GPC", "GPC"),
    c("GPC", "GPC", "GPC", "OPC")
  )
  # tie-breakers that do not break anything
  labs <- c(votes, list(global_1 = c("Microglia", "Microglia", "GPC"),
                        global_2 = c("Microglia", "Microglia", "GPC")))
  cons <- resolve_consensus_labels(labs, methods, c("global_1", "global_2"))
  res  <- suppressMessages(resolve_category_ties(cons, labs, methods, lc))
  expect_equal(as.character(res$label), c("cancer_unresolved", "Unknown", "GPC"))
})

test_that("input errors are caught", {
  lc <- suppressMessages(load_markers(toy_markers))$label_category
  votes <- make_votes(c("GPC", "GPC", "OPC", "OPC"))
  cons  <- list(label = factor(c("Unknown", "Unknown")), method = factor(c("Unresolved", "Unresolved")))
  expect_error(resolve_category_ties(cons, votes, methods, lc), "same length")
  expect_error(resolve_category_ties(list(label = factor("Unknown")), votes, methods, lc), "method")
  expect_error(resolve_category_ties(cons, votes, c(methods, "method_5"), lc), "absent")
})

test_that("the built-in reference panel loads with a category for every label", {
  m <- suppressMessages(load_markers(marker_panels$GBM$GBM_RNA_reference_panel_v1))
  expect_length(m$label_category, 24L)
  expect_false(anyNA(m$label_category))
  expect_equal(m$label_category[["OPC"]], "cancer")
  expect_equal(m$label_category[["OPC normal"]], "normal")
  expect_false("NEU" %in% names(m$label_category))
})

test_that("a missing consensus label is left as it is", {
  lc <- suppressMessages(load_markers(toy_markers))$label_category
  votes <- make_votes(c("GPC", "GPC", "OPC", "OPC"), c("GPC", "GPC", "OPC", "OPC"))
  consensus <- list(label = factor(c(NA, "Unknown")), method = factor(c(NA, "Unresolved")))
  res <- suppressMessages(resolve_category_ties(consensus, votes, methods, lc))
  expect_equal(as.character(res$label), c(NA, "cancer_unresolved"))
})
