# Input conversion ------------------------------------------------------------

#' Convert supported inputs to SingleCellExperiment
#'
#' @param x Seurat object, SingleCellExperiment, or gene-by-cell matrix.
#' @param assay Seurat assay name if x is a Seurat object.
#' @param validate_counts If TRUE, reject negative, non-finite, or clearly
#'   normalized values in the count assay.
#' @return SingleCellExperiment.
#' @export
md_as_sce <- function(x, assay = "RNA", validate_counts = TRUE) {
  if (.md_is_sce(x)) {
    sce <- x
  } else if (.md_is_seurat(x)) {
    .md_require("Seurat", "converting Seurat to SingleCellExperiment")
    sce <- Seurat::as.SingleCellExperiment(x, assay = assay)
  } else if (.md_is_matrix(x)) {
    counts <- .md_sparse(x)
    sce <- SingleCellExperiment(assays = list(counts = counts))
  } else {
    .md_stop("Unsupported input. Provide a Seurat object, SingleCellExperiment, or gene-by-cell matrix.")
  }

  if (is.null(rownames(sce)) || any(rownames(sce) == "")) {
    .md_stop("Features must have rownames/gene identifiers.")
  }
  if (is.null(colnames(sce)) || any(colnames(sce) == "")) {
    colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))
  }
  if (anyDuplicated(rownames(sce))) {
    .md_stop("Feature identifiers must be unique. Store symbols in rowData and use unique feature IDs as rownames.")
  }
  if (anyDuplicated(colnames(sce))) .md_stop("Cell/barcode names must be unique.")
  if (isTRUE(validate_counts)) .md_validate_count_assay(sce)
  sce <- .md_add_logcounts(sce)
  sce
}

.md_validate_count_assay <- function(sce) {
  counts <- .md_sparse(.md_get_counts(sce))
  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    .md_stop("Count assay must contain at least one feature and one cell.")
  }
  vals <- counts@x
  if (length(vals) && any(!is.finite(vals))) .md_stop("Count assay contains non-finite values.")
  if (length(vals) && any(vals < 0)) .md_stop("Count assay contains negative values.")
  if (length(vals)) {
    take <- seq_len(min(length(vals), 100000L))
    non_integer <- mean(abs(vals[take] - round(vals[take])) > 1e-6)
    if (non_integer > 0.01) {
      .md_stop(
        "The selected count assay appears normalized (", round(100 * non_integer, 1),
        "% sampled non-zero values are non-integer). mitoDoubletR requires raw count-like data."
      )
    }
  }
  empty <- which(Matrix::colSums(counts) <= 0)
  if (length(empty)) {
    shown <- paste(head(colnames(counts)[empty], 5L), collapse = ", ")
    .md_stop(
      "Count assay contains ", length(empty), " empty cell(s)",
      if (nzchar(shown)) paste0(" (for example: ", shown, ")") else "",
      ". Remove empty barcodes before adjudication."
    )
  }
  if (nrow(counts) < 200L) .md_warn_once("few_input_genes", "Fewer than 200 genes are present; identity modelling may be unreliable.")
  if (ncol(counts) < 100L) .md_warn_once("few_input_cells", "Fewer than 100 cells are present; adaptive calibration will be unstable.")
  invisible(TRUE)
}
