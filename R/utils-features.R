#' Prepare QC feature groups from regex patterns
#'
#' Scans the dataset for genes matching the provided patterns and constructs
#' the list format required by `assess_cell_quality`.
#'
#' @param object Input object (Seurat, SingleCellExperiment, or Matrix).
#' @param group_configs Named list of definitions. Each entry should contain:
#'   \itemize{
#'     \item \code{pattern}: Regex string to find genes (e.g. "^MT-").
#'     \item \code{max_pct}: The threshold to preserve in the output.
#'   }
#'   Default includes Mito \code{("^MT-")} and Ribo \code{("^RP[SL]")}.
#'
#' @return A named list suitable for the `check_feature_groups` argument.
#' @export
find_qc_features <- function(object,
                             group_configs = list(
                               mito = list(pattern = "^MT-", max_pct = 20),
                               ribo = list(pattern = "^RP[SL]", max_pct = 50)
                             )) {

  counts <- extract_counts(object)
  all_genes <- rownames(counts)

  final_list <- list()

  for (name in names(group_configs)) {
    config <- group_configs[[name]]

    # Find the genes
    found_genes <- grep(config$pattern, all_genes, value = TRUE, ignore.case = TRUE)

    # Construct the explicit list element
    final_list[[name]] <- list(
      features = found_genes,
      max_pct = config$max_pct
    )
  }

  return(final_list)
}

#' Expand gene symbols with synonyms from biomaRt
#'
#' Connects to Ensembl to find synonyms for a list of genes. If `valid_genes` is provided,
#' it filters the results to match genes present in your dataset and supports regex patterns.
#'
#' @param input_features A character vector of gene symbols OR a single regex pattern string.
#' @param valid_genes (Optional) A character vector of all gene names in your dataset.
#'   If provided, used to filter results and allow regex matching. Default NULL.
#' @param organism Character. "hsapiens", "mmusculus", etc.
#'
#' @return A character vector of unique gene symbols found in `valid_genes`.
#' @export
feature_synonyms <- function(input_features, valid_genes = NULL, organism = "hsapiens") {

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required for synonym lookup. Please install it.")
  }

  # 1. Determine if input is a regex or a list of genes
  if (!is.null(valid_genes) && length(input_features) == 1 && grepl("[^a-zA-Z0-9]", input_features)) {
    search_genes <- grep(input_features, valid_genes, value = TRUE, ignore.case = TRUE)
  } else {
    search_genes <- input_features
  }

  if (length(search_genes) == 0) return(character(0))

  # 2. Setup BiomaRt
  dataset_name <- switch(tolower(organism),
                         "human" = "hsapiens_gene_ensembl",
                         "hsapiens" = "hsapiens_gene_ensembl",
                         "mouse" = "mmusculus_gene_ensembl",
                         "mmusculus" = "mmusculus_gene_ensembl",
                         stop("Organism not supported automatically."))

  mart <- tryCatch({
    biomaRt::useMart("ensembl", dataset = dataset_name)
  }, error = function(e) {
    warning("Could not connect to biomaRt: ", e$message)
    return(NULL)
  })

  if (is.null(mart)) return(search_genes)

  # 3. Query
  results <- biomaRt::getBM(
    attributes = c("external_gene_name", "external_synonym"),
    filters = "external_gene_name",
    values = search_genes,
    mart = mart
  )

  all_variants <- unique(c(results$external_gene_name, results$external_synonym))

  # 4. Filter synonyms against what is actually in the data
  if (!is.null(valid_genes)) {
    return(intersect(all_variants, valid_genes))
  }

  return(all_variants)
}
