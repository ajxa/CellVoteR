#' Plot Fisher Score Heatmap
#'
#' @param scores_df Dataframe output from calculate_fisher_scores.
#'
#' @export
plot_fisher_heatmap <- function(scores_df) {
  # FIX: Silence R CMD check notes
  Cluster <- Reference_Type <- NegLog10P <- NULL

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  # FIX: gtools::mixedsort
  scores_df$Cluster <- factor(scores_df$Cluster, levels = unique(gtools::mixedsort(as.character(scores_df$Cluster))))

  p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = Cluster, y = Reference_Type, fill = NegLog10P)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_viridis_c(option = "magma", name = "-log10(P)") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Cluster Identity Scores", x = "Cluster", y = "Reference Type") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  return(p)
}
