#' Execute a single annotation pipeline method
#'
#' Orchestrates the complete hierarchical annotation workflow for a specific
#' configuration (mode + feature filtering). It performs broad classification,
#' splits the data, sub-clusters/annotates each partition, and merges the results.
#'
#' @details
#' \strong{Naming Convention:}
#' Output columns are prefixed with the method ID, constructed as \code{mode_filter}.
#' For example, if \code{mode = "clustering"} and \code{filter_markers = TRUE}:
#' \itemize{
#'   \item \strong{ID:} \code{clustering_ref}
#'   \item \strong{Columns:} \code{clustering_ref_broad_label}, \code{clustering_ref_best_label}, etc.
#' }
#'
#' \strong{Workflow:}
#' \enumerate{
#'   \item \strong{Broad Classification:} Depending on \code{mode}, uses either
#'   \code{\link{annotate_broad_clusters}} (clustering-based) or
#'   \code{\link{annotate_broad_cells}} (marker-based) to assign initial categories.
#'   \item \strong{Splitting:} The object is split into a list of sub-objects based on these broad labels.
#'   \item \strong{Sub-processing:} For each broad category:
#'   \itemize{
#'     \item \code{\link{subcluster_cells}} is run to find high-resolution clusters.
#'     \item \code{\link{score_subclusters}} matches these against fine-grained markers.
#'     \item \code{\link{get_best_label}} determines the final identity.
#'   }
#'   \item \strong{Merging:} The sub-objects are merged back into a single Seurat object
#'   containing the new high-resolution metadata.
#' }
#'
#' @param seu_obj A Seurat object.
#' @param broad_markers A configuration list for broad cell types (see \code{\link{build_broad_marker_config}}).
#' @param fine_markers A named list of fine-grained markers for sub-clustering.
#' @param mode Character. The broad classification strategy: \code{"clustering"} or \code{"marker"}.
#' @param filter_markers Logical. If \code{TRUE}, restricts the feature space to the union of
#'   provided markers before running broad classification (optimizes speed/relevance).
#' @param return_single_object Logical. If \code{TRUE}, merges results back into a single Seurat object.
#'   If \code{FALSE}, returns a list of processed sub-objects (useful for debugging).
#' @param skip_subclustering Logical. If \code{TRUE}, stops after broad classification.
#' @param verbose Logical. Controls the amount of CLI output.
#'
#' @return A Seurat object with updated metadata columns for broad and fine labels.
#'   If \code{return_single_object = FALSE}, returns a list of Seurat objects and scoring tables.
#' @export
#'
run_method <- function(
    seu_obj,
    broad_markers,
    fine_markers,
    mode = c("clustering", "marker"),
    filter_markers = FALSE,
    return_single_object = TRUE,
    skip_subclustering = FALSE,
    verbose = FALSE
) {

  be_quiet <- function(expr) {
    if (!verbose) {
      trash <- utils::capture.output({
        res <- suppressMessages(expr)
      })
      return(res)
    } else expr
  }

  mode <- match.arg(mode)

  # Run configuration ----

  if (filter_markers) {
    filter_suffix <- "ref"
    feats <- .feature_space(broad_markers, fine_markers)
  } else {
    filter_suffix <- "all"
    feats <- NULL
  }

  method_args <- list(
    mode            =   mode,
    id              =   paste(mode, filter_suffix, sep = "_"),
    filter_features =   feats,
    filter_lab      =   paste(filter_suffix, "features")
  )

  cols <- list(
    broad_cluster  =  paste0(method_args$id, "_broad_cluster"),
    broad_label    =  paste0(method_args$id, "_broad_label"),
    subcluster     =  paste0(method_args$id, "_subcluster"),
    best_label     =  paste0(method_args$id, "_best_label")
  )

  cli::cli_h1("Method ({method_args$mode} - {method_args$filter_lab})")

  # Broad labelling ----
  bid <- NULL
  if (!verbose) {
    label_type <- ifelse(mode == "clustering", "clustering", "marker-based")

    bid <- cli::cli_progress_step(
      msg = paste0("Performing broad classification (", label_type, ")..."),
      spinner = TRUE
    )
  }

  if (mode == "clustering") {

    seu_obj <- be_quiet({
      cluster_cells(
        sObj = seu_obj,
        cluster_name = cols$broad_cluster,
        filter_features = method_args$filter_features,
        verbose = FALSE
      )
    })

    seu_obj <- be_quiet({
      annotate_broad_clusters(
        sObj = seu_obj,
        anno_markers = broad_markers,
        initial_cluster_col = cols$broad_cluster,
        label_out_col = cols$broad_label
      )
    })

  } else {

    seu_obj <- be_quiet({
      annotate_broad_cells(
        sObj = seu_obj,
        broad_marker_list = broad_markers,
        out_label_col = cols$broad_label,
        filter_features = method_args$filter_features
      )
    })
  }

  if (!verbose && !is.null(bid)) cli::cli_progress_done(id = bid)

  # Sub-clustering broad labels ----
  subclusters <- Seurat::SplitObject(seu_obj, split.by = cols$broad_label)

  if(skip_subclustering) return(subclusters)

  # ensures that all sub cluster labels have corresponding fine markers
  .check_missing_feats(names(subclusters), names(fine_markers))

  sid <- NULL
  if (!verbose) {
    n_subs <- length(subclusters)
    sid <- cli::cli_progress_step(
      msg = paste0("Sub-clustering ", n_subs, " broad label categories..."),
      spinner = TRUE
    )
  }

  for (sub_name in names(subclusters)){

    sub_res <- be_quiet({
      subcluster_cells(
        sub_sObj = subclusters[[sub_name]],
        sub_marker_list = fine_markers[[sub_name]],
        subcluster_name = cols$subcluster
      )
    })

    scores <- be_quiet({
      score_subclusters(
        sub_sObj = sub_res$clusters,
        sub_de_markers = sub_res$de_markers,
        sub_ref_markers = sub_res$ref_markers,
        sub_fisher_background = sub_res$fisher_test_background
      )
    })

    best_lbl <- get_best_label(sub_score_df = scores)

    sub_res$clusters <- be_quiet({
      assign_labels(
        sub_sObj = sub_res$clusters,
        sub_best_label_map = best_lbl,
        sub_sObj_col = cols$subcluster,
        sub_sObj_new_col = cols$best_label
      )
    })

    sub_res[["label_scores"]] <- scores
    sub_res[["best_label"]] <- best_lbl
    subclusters[[sub_name]] <- sub_res

  }

  if (!verbose && !is.null(sid)) cli::cli_progress_done(id = sid)

  # Merging sub-clustered objects ----

  if (return_single_object) {
    mid <- NULL
    if (!verbose) mid <- cli::cli_progress_step("Merging subclusters...")

    subclusters <- be_quiet({
      lapply(subclusters, `[[`, "clusters")
    })

    if (length(subclusters) > 1){
      subclusters <- merge(x = subclusters[[1]], y = subclusters[-1])
    } else subclusters <- subclusters[[1]]

    subclusters <- SeuratObject::JoinLayers(subclusters)

    if (!verbose && !is.null(mid)) {
      cli::cli_progress_done(id = mid)
      cli::cli_alert_success("Method complete.")
    }

    return(subclusters)
  }

  if (!verbose) cli::cli_alert_success("Method complete.")

  return(subclusters)

}
