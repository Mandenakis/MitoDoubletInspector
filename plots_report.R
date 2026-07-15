# Rescue candidate extraction --------------------------------------------------

#' Extract candidate rescued mitotic/polyploid-like singlets
#'
#' Returns cells originally flagged as doublet/multiplet by the imported or run
#' doublet caller, but adjudicated as `mitotic_singlet` or
#' `polyploid_like_singlet` by mitoDoubletR. These are the cells a naive filter
#' would have deleted.
#'
#' @param x SingleCellExperiment or Seurat object after `md_adjudicate()`.
#' @param rare_col Optional metadata column defining a rare population/cluster.
#' @param rare_value Optional value(s) inside `rare_col` to retain.
#' @param classes Rescue classes to return.
#' @param include_unresolved_high_rna If FALSE, omit polyploid-like candidates
#'   whose RNA profile is intrinsically indistinguishable from a homotypic
#'   doublet. They remain available in the full reason table.
#' @return Data.frame.
#' @export
md_rescue_candidates <- function(
  x,
  rare_col = NULL,
  rare_value = NULL,
  classes = c("mitotic_singlet", "polyploid_like_singlet"),
  include_unresolved_high_rna = TRUE
) {
  cd <- if (.md_is_seurat(x)) x[[]] else .md_coldata_df(x)
  if (!"md_class" %in% colnames(cd)) .md_stop("No md_class column found. Run md_adjudicate() first.")
  original_flag <- rep(FALSE, nrow(cd))
  if ("md_original_doublet_class" %in% colnames(cd)) {
    original_flag <- tolower(cd$md_original_doublet_class) %in% c("doublet", "multiplet")
  }
  keep <- original_flag & cd$md_class %in% classes
  if ("md_rescue_eligible" %in% colnames(cd)) keep <- keep & as.logical(cd$md_rescue_eligible)
  if (!isTRUE(include_unresolved_high_rna)) keep <- keep & cd$md_class != "polyploid_like_singlet"
  if (!is.null(rare_col)) {
    if (!rare_col %in% colnames(cd)) .md_stop("rare_col not found: ", rare_col)
    if (!is.null(rare_value)) keep <- keep & cd[[rare_col]] %in% rare_value
  }
  out <- cd[keep, grep("^md_|^seurat_clusters$|^cluster$|^celltype$|^cell_type$", colnames(cd), value = TRUE), drop = FALSE]
  out$cell <- rownames(out)
  out[, c("cell", setdiff(colnames(out), "cell")), drop = FALSE]
}
