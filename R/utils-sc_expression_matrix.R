#' Create/Load a SingleCellExperiment from raw counts
#'
#' Entry point for the CellVoteR annotation workflow. Accepts a raw (un-normalised)
#' gene-by-cell counts matrix and constructs a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object ready for
#' downstream QC, normalisation, and clustering.
#'
#' @section Input precedence:
#' The function accepts two mutually exclusive input modes. If \code{counts}
#' is provided it takes precedence and the MTX arguments are ignored (with a
#' warning if also supplied). Otherwise, all three MTX arguments must be
#' provided together.
#'
#' @section Supported input formats:
#' \describe{
#'   \item{\strong{counts (in-memory)}}{A \code{dgCMatrix} with gene
#'     identifiers as rownames and cell barcodes as colnames.}
#'   \item{\strong{counts (RDS file)}}{A character path to a \code{.rds} file
#'     containing a \code{dgCMatrix} as described above.}
#'   \item{\strong{MTX triplet files}}{Three file paths supplied via
#'     \code{mtx_file}, \code{cells_file}, and \code{genes_file}. The MTX
#'     file contains the sparse matrix in Matrix Market format. The cells and
#'     genes files are text files with one identifier per line (or
#'     tab-separated, in which case the second column is used preferentially).
#'     All three must be provided together.}
#'  \item{\strong{cell_metadata (RDS or CSV/TSV file)}}{A character path to
#'     an \code{.rds} file containing a \code{data.frame}, or a path to a
#'     \code{.csv} or \code{.tsv}/\code{.txt} file. The file must contain one
#'     row per cell, in the same order as the columns of the counts matrix.}
#' }
#'
#' @section Disk-backed extension:
#' The returned SCE holds an in-memory \code{dgCMatrix} counts assay. For
#' datasets that exceed available RAM, the object can be converted to
#' disk-backed storage after creation:
#' \preformatted{
#' sce <- create_sce(counts = "large_counts.rds")
#' HDF5Array::saveHDF5SummarizedExperiment(sce, dir = "my_hdf5_sce")
#' sce <- HDF5Array::loadHDF5SummarizedExperiment("my_hdf5_sce")
#' }
#'
#' @param counts A gene-by-cell raw counts input. Either:
#'   \itemize{
#'     \item A \code{dgCMatrix} (in-memory sparse matrix with rownames and
#'       colnames).
#'     \item A character scalar path to a \code{.rds} file containing a
#'       \code{dgCMatrix}.
#'   }
#'   When provided, takes precedence over the MTX arguments.
#' @param mtx_file Character scalar. Path to a Matrix Market (\code{.mtx} or
#'   \code{.mtx.gz}) file. Must be accompanied by \code{cells_file} and
#'   \code{genes_file}.
#' @param cells_file Character scalar. Path to a text file containing one cell
#'   identifier per line (no header). Required when \code{mtx_file} is
#'   provided.
#' @param genes_file Character scalar. Path to a text file containing one gene
#'   identifier per line (no header). Required when \code{mtx_file} is
#'   provided.
#' @param cell_metadata Optional per-cell annotations. Either:
#'   \itemize{
#'     \item A \code{data.frame} with one row per cell.
#'     \item A character scalar path to an \code{.rds} file containing a
#'       \code{data.frame}.
#'     \item A character scalar path to a \code{.csv}, \code{.tsv}, or
#'       \code{.txt} file (tab-separated).
#'   }
#'   Must have one row per column in the counts matrix. This will be stored in \code{colData()}.
#' @param gene_metadata Optional \code{data.frame} of per-gene annotations.
#'   Must have one row per row in the counts matrix. Stored in
#'   \code{rowData()}.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} with a
#'   single \code{counts} assay stored as a \code{dgCMatrix}.
#'
#' @examples
#' \dontrun{
#' # From an in-memory sparse matrix
#' sce <- create_sce(counts = my_sparse_counts)
#'
#' # From an RDS file with metadata from a CSV
#' sce <- create_sce(
#'   counts        = "raw_counts.rds",
#'   cell_metadata = "cell_metadata.csv"
#' )
#'
#' # From MTX triplet files with metadata from an RDS
#' sce <- create_sce(
#'   mtx_file      = "data/matrix.mtx.gz",
#'   cells_file    = "data/barcodes.tsv",
#'   genes_file    = "data/features.tsv",
#'   cell_metadata = "data/cell_metadata.rds"
#' )
#' }
#'
#' @seealso
#' \code{\link{load_markers}} and \code{\link{build_broad_marker_config}} for
#' preparing the marker configuration to attach to the SCE via
#' \code{metadata()}.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom Matrix readMM
#' @export
create_sce <- function(counts = NULL,
                       mtx_file = NULL,
                       cells_file = NULL,
                       genes_file = NULL,
                       cell_metadata = NULL,
                       gene_metadata = NULL) {

  # --- 1.) Resolve the counts ----

  mtx_args <- c(!is.null(mtx_file), !is.null(cells_file), !is.null(genes_file))

  if (!is.null(counts)) {

    if (any(mtx_args)) {
      print_alert(
        "Both counts and MTX file arguments provided - using 'counts'", type = "w"
      )
    }

    counts <- resolve_counts(counts)

  } else if (any(mtx_args)) {

    if (!all(mtx_args)) {

      provided <- c("mtx_file", "cells_file", "genes_file")[mtx_args]
      missing  <- c("mtx_file", "cells_file", "genes_file")[!mtx_args]

      cli::cli_abort(c(
        "Incomplete MTX input.",
        "x" = "Provided: {.arg {provided}}",
        "x" = "Missing:  {.arg {missing}}",
        "i" = "All three of {.arg mtx_file}, {.arg cells_file}, and
               {.arg genes_file} are required."
      ))
    }

    counts <- read_mtx_triplet(mtx_file, cells_file, genes_file)

  } else cli::cli_abort("Must provide either {.arg counts} or all three MTX file arguments.")

  validate_counts_matrix(counts)

  # --- 2.) Resolve the cell metadata ----

  if (!is.null(cell_metadata)) {
    cell_metadata <- resolve_cell_metadata(cell_metadata)

    validate_metadata(cell_metadata,
                      expected_n = ncol(counts),
                      dim_label  = "cells",
                      arg_name   = "cell_metadata")
  }

  # --- 3.) Resolve gene metadata ----

  if (!is.null(gene_metadata)) {

    validate_metadata(gene_metadata,
                      expected_n = nrow(counts),
                      dim_label = "genes",
                      arg_name = "gene_metadata"
                      )
  }

  # --- 4.) Construct the SCE object ----

  sce_args <- list(assays = list(counts = counts))

  if (!is.null(cell_metadata)) sce_args$colData <- S4Vectors::DataFrame(cell_metadata)
  if (!is.null(gene_metadata)) sce_args$rowData <- S4Vectors::DataFrame(gene_metadata)

  sce <- do.call(SingleCellExperiment::SingleCellExperiment, sce_args)

  print_alert(
    "loaded expression data: {.val {nrow(sce)}} genes x {.val {ncol(sce)}} cells", type = "s"
  )

  return(sce)
}

#' Resolve cell metadata from a data.frame or file path
#'
#' Accepts a \code{data.frame} directly, or a path to an \code{.rds},
#' \code{.csv}, \code{.tsv}, or \code{.txt} file.
#'
#' @param cell_metadata A \code{data.frame} or character scalar file path.
#' @return A \code{data.frame}.
#' @keywords internal
#' @noRd
resolve_cell_metadata <- function(cell_metadata) {

  if (is.data.frame(cell_metadata)) return(cell_metadata)

  if (!is.character(cell_metadata) || length(cell_metadata) != 1L) {
    cli::cli_abort(
      "{.arg cell_metadata} must be a {.cls data.frame} or a single file path."
    )
  }

  if (!file.exists(cell_metadata)) cli::cli_abort("File not found: {.path {cell_metadata}}.")

  ext <- tolower(tools::file_ext(cell_metadata))

  out <- switch(ext,
                rds = {
                  obj <- readRDS(cell_metadata)
                  if (!is.data.frame(obj)) {
                    cli::cli_abort(
                      "The RDS file {.path {cell_metadata}} must contain a {.cls data.frame}, \\
           not a {.cls {class(obj)}}."
                    )
                  }
                  obj
                },
                csv = utils::read.csv(cell_metadata,  stringsAsFactors = FALSE),
                tsv = ,
                txt = utils::read.delim(cell_metadata, stringsAsFactors = FALSE),
                cli::cli_abort(
                  "Unsupported file extension {.val {ext}} for {.arg cell_metadata}. \\
       Supported formats: {.val rds}, {.val csv}, {.val tsv}, {.val txt}."
                )
  )

  return(out)
}


#' Resolve the counts argument to a dgCMatrix
#'
#' Handles both in-memory \code{dgCMatrix} objects and RDS file paths.
#'
#' @param counts A \code{dgCMatrix} or character scalar RDS path.
#' @return A \code{dgCMatrix}.
#' @keywords internal
#' @noRd
resolve_counts <- function(counts) {

  if (is.character(counts) && length(counts) == 1L) {

    checkmate::assert_file_exists(counts, access = "r")

    ext <- tolower(tools::file_ext(counts))

    if (ext != "rds") {

      cli::cli_abort(
        "Unsupported file extension {.val {ext}} for {.arg counts}. Expected {.val rds}."
      )

    }

    print_alert("Reading file: {.path {counts}}", type = "p")

    counts <- readRDS(counts)

  }

  if (!inherits(counts, "dgCMatrix")) {

    cli::cli_abort(
      "{.arg counts} must be a {.cls dgCMatrix} or a path to an RDS file
       containing one. Received: {.cls {class(counts)}}"
    )

  }

  return(counts)
}


#' Read MTX triplet files into a dgCMatrix
#'
#' Reads a Matrix Market file and companion cell/gene identifier files,
#' validates dimensional consistency, and returns a named sparse matrix.
#'
#' @param mtx_file Character scalar. Path to the \code{.mtx}/\code{.mtx.gz}.
#' @param cells_file Character scalar. Path to the cell identifiers file.
#' @param genes_file Character scalar. Path to the gene identifiers file.
#' @return A \code{dgCMatrix} with gene rownames and cell colnames.
#' @keywords internal
#' @noRd
read_mtx_triplet <- function(mtx_file, cells_file, genes_file) {

  checkmate::assert_file_exists(mtx_file, access = "r")
  checkmate::assert_file_exists(cells_file, access = "r")
  checkmate::assert_file_exists(genes_file, access = "r")

  print_alert("Reading MTX triplet files:")
  print_alert("matrix: {.path {mtx_file}}", type = "p")
  print_alert("cells:  {.path {cells_file}}", type = "p")
  print_alert("genes:  {.path {genes_file}}", type = "p")

  mat   <- Matrix::readMM(mtx_file)
  mat   <- methods::as(mat, "CsparseMatrix")
  cells <- read_id_file(cells_file, label = "cell")
  genes <- read_id_file(genes_file, label = "gene")

  if (length(genes) != nrow(mat)) {
    cli::cli_abort(
      "Gene count mismatch: {.path {genes_file}} has {.val {length(genes)}}
       identifiers but matrix has {.val {nrow(mat)}} rows."
    )
  }

  if (length(cells) != ncol(mat)) {
    cli::cli_abort(
      "Cell count mismatch: {.path {cells_file}} has {.val {length(cells)}}
       identifiers but matrix has {.val {ncol(mat)}} columns."
    )
  }

  rownames(mat) <- genes
  colnames(mat) <- cells

  return(mat)
}


#' Read a single-column identifier file
#'
#' Reads a text file containing one identifier per line. If the file is
#' tab-separated with multiple columns (e.g. 10x \code{features.tsv} with
#' Ensembl ID, symbol, and type), the second column is used preferentially
#' as it typically contains gene symbols.
#'
#' @param path Character scalar. File path.
#' @param label Character scalar. Descriptor for error messages.
#' @return Character vector of identifiers.
#' @keywords internal
#' @noRd
read_id_file <- function(path, label) {

  lines <- readLines(path)
  lines <- lines[nchar(trimws(lines)) > 0L]

  if (length(lines) == 0L) {
    cli::cli_abort("{label} identifier file is empty: {.path {path}}")
  }

  first_line <- strsplit(lines[1L], "\t")[[1L]]

  if (length(first_line) > 1L) {
    parsed  <- strsplit(lines, "\t")
    col_idx <- if (length(first_line) >= 2L) 2L else 1L
    ids     <- vapply(parsed, function(x) x[col_idx], character(1))

    print_alert(
      "Detected multi-column {label} file - using column {.val {col_idx}}", type = "w"
    )

  } else ids <- trimws(lines)

  return(ids)
}


#' Validate a counts matrix
#' @param mat A \code{dgCMatrix}.
#' @keywords internal
#' @noRd
validate_counts_matrix <- function(mat) {

  if (nrow(mat) == 0L || ncol(mat) == 0L) cli::cli_abort("Counts matrix has zero-length dimension: {nrow(mat)} x {ncol(mat)}")

  if (is.null(rownames(mat))) cli::cli_abort("Counts matrix must have rownames (gene identifiers).")

  if (is.null(colnames(mat))) cli::cli_abort("Counts matrix must have colnames (cell barcodes/identifiers).")

  if (anyDuplicated(rownames(mat))) {

    n_dup <- sum(duplicated(rownames(mat)))

    cli::cli_abort(
      "Counts matrix has {.val {n_dup}} duplicate gene name{?s}.
       Please deduplicate before input."
    )
  }

  if (anyDuplicated(colnames(mat))) {

    n_dup <- sum(duplicated(colnames(mat)))

    cli::cli_abort(
      "Counts matrix has {.val {n_dup}} duplicate cell barcode{?s}.
       Please deduplicate before input."
    )
  }

  if (any(mat@x < 0)) cli::cli_abort("Counts matrix contains negative values. Expected raw, un-normalised counts.")

  sample_vals <- mat@x[seq_len(min(1000L, length(mat@x)))]

  if (length(sample_vals) > 0L && all(sample_vals != floor(sample_vals))) {

    print_alert(
      "Values do not appear to be integer counts...
        This workflow expects raw, un-normalised counts.",
      type = "w"
      )

  }

}


#' Validate optional metadata dimensions
#' @param df A \code{data.frame}.
#' @param expected_n Integer.
#' @param dim_label Character.
#' @param arg_name Character.
#' @keywords internal
#' @noRd
validate_metadata <- function(df, expected_n, dim_label, arg_name) {

  if (!is.data.frame(df)) cli::cli_abort("{.arg {arg_name}} must be a data.frame.")

  if (nrow(df) != expected_n) {

    cli::cli_abort(
      "{.arg {arg_name}} has {.val {nrow(df)}} rows but counts has
       {.val {expected_n}} {dim_label}."
    )

  }

}
