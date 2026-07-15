#' mitoDoubletR: Mitosis-Aware Doublet Adjudication for scRNA-seq
#'
#' A second-stage adjudication framework for single-cell RNA-seq doublet calls.
#' Instead of replacing a doublet caller, mitoDoubletR asks whether a suspicious
#' barcode is better explained by one coherent identity carrying a
#' mitotic/cytokinetic or high-RNA/polyploid-like programme, or genuinely needs
#' two parent transcriptomes. See \code{\link{md_adjudicate}} for the one-command
#' workflow.
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim reducedDimNames
#' @importFrom SummarizedExperiment assay assay<- assayNames colData colData<- rowData rowData<-
#' @keywords internal
"_PACKAGE"

# Column names referenced inside ggplot2::aes() by the plotting helpers. Declared
# here so R CMD check does not flag them as undefined globals.
utils::globalVariables(c(
  "md_doublet_score_full", "md_doublet_score_no_cycle", "md_class",
  "md_mitosis_score", "md_lineage_conflict", "dim1", "dim2"
))
