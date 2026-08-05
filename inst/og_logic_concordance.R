markers <- load_markers(marker_panels$GBM$nomura)

markers$broad <- build_broad_marker_config(
  marker_list = markers$broad,
  priority_order = c("vasculature", "immune"),
  # default_threshold = 0.25                    # keep this a 0.1
)

sce <- list(
  Imp14_P = list.files("~/Desktop/example_data", pattern = "Imp14_P", full.names = TRUE),
  Imp14_R = list.files("~/Desktop/example_data", pattern = "Imp14_R", full.names = TRUE)
)

sce <- lapply(sce, \(x) create_sce(counts = x[[1]], cell_metadata = x[[2]]))

get_og_labs <- function(sce_data = sce, sample){

  lab_index <- seq(
    grep("(?i)method_1$", colnames(sce_data[[sample]]@colData)),
    grep("(?i)^resolved.+clash$", colnames(sce_data[[sample]]@colData))
  )

  return(as.data.frame(sce_data[[sample]]@colData[,lab_index]))
}

og_labs <- setNames(object = as.list(names(sce)), nm = names(sce))

og_labs <- lapply(og_labs, \(x) get_og_labs(sample = x))

rm(get_og_labs)


sce_clean <- lapply(sce, \(x){

  out <- assess_cell_quality(x, remove_failed_cells = TRUE)

  out <- normalize_counts(out)

  out <- prepare_sce(out, markers)

  return(out)

})


results <- setNames(object = vector("list", length(sce_clean)), nm = names(sce))
consensus <- setNames(object = vector("list", length(sce_clean)), nm = names(sce))

results$Imp14_P <- run_cellvoter(sce_clean$Imp14_P)
results$Imp14_R <- run_cellvoter(sce_clean$Imp14_R)

consensus$Imp14_P <- resolve_consensus_labels(
  label_list        = results$Imp14_P$labels,
  method_names      = results$Imp14_P$method_names,
  tie_breaker_names = results$Imp14_P$tie_breaker_names,
  unassigned_label  = "unknown",
  allow_even_split  = FALSE,
  ordered_tiebreak  = TRUE
)

consensus$Imp14_R <- resolve_consensus_labels(
  label_list        = results$Imp14_R$labels,
  method_names      = results$Imp14_R$method_names,
  tie_breaker_names = results$Imp14_R$tie_breaker_names,
  unassigned_label  = "unknown",
  allow_even_split  = FALSE,
  ordered_tiebreak  = TRUE
)


results$Imp14_P$sce$cellVoteR_label  <- consensus$Imp14_P$label
results$Imp14_P$sce$cellVoteR_method  <- consensus$Imp14_P$method

results$Imp14_R$sce$cellVoteR_label  <- consensus$Imp14_R$label
results$Imp14_R$sce$cellVoteR_method  <- consensus$Imp14_R$method



table(results$Imp14_P$labels$method_1)
table(og_labs$Imp14_P$Fisher_Method_1)
