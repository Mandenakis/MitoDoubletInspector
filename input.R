# Plots and reports ------------------------------------------------------------

.md_plot_df <- function(x) {
  if (.md_is_seurat(x)) x[[]] else .md_coldata_df(x)
}

#' Plot full versus no-cycle doublet scores
#'
#' @param x SingleCellExperiment or Seurat object after mitoDoubletR.
#' @return ggplot.
#' @export
md_plot_cycle_drop <- function(x) {
  df <- .md_plot_df(x)
  req <- c("md_doublet_score_full", "md_doublet_score_no_cycle", "md_class")
  if (!all(req %in% colnames(df))) .md_stop("Required columns missing: ", paste(setdiff(req, colnames(df)), collapse = ", "))
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = md_doublet_score_full, y = md_doublet_score_no_cycle, colour = md_class)
  ) +
    ggplot2::geom_point(alpha = 0.65, size = 1.1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::labs(
      x = "Doublet score, all genes",
      y = "Doublet score, cycle/cytokinesis genes removed",
      colour = "miDAS class",
      title = "Counterfactual doublet evidence"
    ) +
    ggplot2::theme_bw()
}

#' Plot mitotic signal versus lineage conflict
#'
#' @param x SingleCellExperiment or Seurat object after mitoDoubletR.
#' @return ggplot.
#' @export
md_plot_mitosis_conflict <- function(x) {
  df <- .md_plot_df(x)
  req <- c("md_mitosis_score", "md_lineage_conflict", "md_class")
  if (!all(req %in% colnames(df))) .md_stop("Required columns missing: ", paste(setdiff(req, colnames(df)), collapse = ", "))
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = md_mitosis_score, y = md_lineage_conflict, colour = md_class)
  ) +
    ggplot2::geom_point(alpha = 0.65, size = 1.1) +
    ggplot2::labs(
      x = "Mitotic programme score",
      y = "Lineage conflict score",
      colour = "miDAS class",
      title = "Mitotic programme versus transcriptomic hybridity"
    ) +
    ggplot2::theme_bw()
}

#' Plot mitoDoubletR classes on UMAP
#'
#' @param x SingleCellExperiment or Seurat object with UMAP coordinates.
#' @param reduction Reduction name. For Seurat, usually "umap". For SCE,
#'   usually "UMAP".
#' @return ggplot.
#' @export
md_plot_class_umap <- function(x, reduction = NULL) {
  if (.md_is_seurat(x)) {
    .md_require("SeuratObject", "extracting Seurat reductions")
    reduction <- reduction %||% "umap"
    emb <- SeuratObject::Embeddings(x, reduction = reduction)
    df <- cbind(as.data.frame(emb[, 1:2, drop = FALSE]), x[[]])
    colnames(df)[1:2] <- c("dim1", "dim2")
  } else {
    rds <- reducedDimNames(x)
    reduction <- reduction %||% if ("UMAP" %in% rds) "UMAP" else rds[1]
    if (length(rds) == 0L || is.na(reduction) || !reduction %in% rds) .md_stop("No reduced dimension found for plotting.")
    emb <- reducedDim(x, reduction)
    if (ncol(emb) < 2L) .md_stop("Reduced dimension needs at least two components.")
    df <- cbind(as.data.frame(emb[, 1:2, drop = FALSE]), .md_coldata_df(x))
    colnames(df)[1:2] <- c("dim1", "dim2")
  }
  if (!"md_class" %in% colnames(df)) .md_stop("No md_class column found.")
  ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, colour = md_class)) +
    ggplot2::geom_point(alpha = 0.75, size = 0.8) +
    ggplot2::labs(x = paste0(reduction, "_1"), y = paste0(reduction, "_2"), colour = "miDAS class") +
    ggplot2::theme_bw()
}

#' Write a mitoDoubletR report directory
#'
#' @param x SingleCellExperiment or Seurat object after mitoDoubletR.
#' @param outdir Output directory.
#' @param prefix File prefix.
#' @return Invisibly returns output paths.
#' @export
md_report <- function(x, outdir = "mitoDoubletR_report", prefix = "md") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  cd <- .md_plot_df(x)
  meta_path <- file.path(outdir, paste0(prefix, "_cell_metadata.csv"))
  summary_path <- file.path(outdir, paste0(prefix, "_class_summary.csv"))
  reason_path <- file.path(outdir, paste0(prefix, "_reason_table.csv"))

  utils::write.csv(cd, meta_path, row.names = TRUE)
  utils::write.csv(md_summary(x), summary_path, row.names = FALSE)
  if (!.md_is_seurat(x)) {
    utils::write.csv(md_reason_table(x), reason_path, row.names = FALSE)
  } else {
    wanted <- grep("^md_", colnames(cd), value = TRUE)
    tmp <- cd[, wanted, drop = FALSE]
    tmp$cell <- rownames(tmp)
    utils::write.csv(tmp, reason_path, row.names = FALSE)
  }

  plot_paths <- character()
  p1 <- tryCatch(md_plot_cycle_drop(x), error = function(e) NULL)
  if (!is.null(p1)) {
    pp <- file.path(outdir, paste0(prefix, "_cycle_drop.png"))
    ggplot2::ggsave(pp, p1, width = 7, height = 5, dpi = 300)
    plot_paths <- c(plot_paths, pp)
  }
  p2 <- tryCatch(md_plot_mitosis_conflict(x), error = function(e) NULL)
  if (!is.null(p2)) {
    pp <- file.path(outdir, paste0(prefix, "_mitosis_conflict.png"))
    ggplot2::ggsave(pp, p2, width = 7, height = 5, dpi = 300)
    plot_paths <- c(plot_paths, pp)
  }
  invisible(list(metadata = meta_path, summary = summary_path, reasons = reason_path, plots = plot_paths))
}

