# RNA residual scoring ---------------------------------------------------------

#' Score within-identity RNA-content residuals
#'
#' Fits a linear model for log10(nCount + 1) using nFeature, identity, sample,
#' mitochondrial fraction, and stress score when available. Positive residuals
#' indicate unexpectedly high RNA content *for the assigned identity* - the signal
#' expected from polyploid-like/high-RNA singlets - rather than a raw nCount
#' threshold that would penalise legitimately large cell states.
#'
#' @param sce SingleCellExperiment.
#' @param sample_col Optional sample/batch column.
#' @return SingleCellExperiment.
#' @export
md_score_rna_residuals <- function(sce, sample_col = NULL) {
  cd <- .md_coldata_df(sce)
  if (!all(c("md_nCount", "md_nFeature") %in% colnames(cd))) {
    sce <- md_add_qc_metrics(sce, sample_col = sample_col)
    cd <- .md_coldata_df(sce)
  }
  identity_col <- if ("md_input_cluster" %in% colnames(cd)) "md_input_cluster" else "md_identity_top"
  if (!identity_col %in% colnames(cd)) {
    colData(sce)$md_identity_top <- rep("unknown", ncol(sce))
    identity_col <- "md_identity_top"
    cd <- .md_coldata_df(sce)
  }

  identity <- as.character(cd[[identity_col]])
  identity[is.na(identity) | identity == ""] <- "unknown"
  df <- data.frame(
    y = log10(cd$md_nCount + 1),
    nfeature = log10(cd$md_nFeature + 1),
    identity = factor(identity),
    percent_mt = if ("md_percent_mt" %in% colnames(cd)) cd$md_percent_mt else 0,
    stress = if ("md_stress_score" %in% colnames(cd)) cd$md_stress_score else 0,
    stringsAsFactors = FALSE
  )
  if (!is.null(sample_col) && sample_col %in% colnames(cd)) {
    sample_value <- as.character(cd[[sample_col]])
    sample_value[is.na(sample_value) | sample_value == ""] <- "unknown"
    df$sample <- factor(sample_value)
    f <- y ~ nfeature + identity + sample + percent_mt + stress
  } else {
    f <- y ~ nfeature + identity + percent_mt + stress
  }

  # Drop factor terms that have a single level to keep lm() well posed.
  if (nlevels(df$identity) < 2L) f <- stats::update(f, . ~ . - identity)
  if (!is.null(df$sample) && nlevels(df$sample) < 2L) f <- stats::update(f, . ~ . - sample)

  fit <- tryCatch(stats::lm(f, data = df), error = function(e) NULL)
  if (is.null(fit)) {
    .md_warn_once("rna_resid_lm_failed", "RNA residual model failed; using centered log-counts.")
    resid <- df$y - mean(df$y, na.rm = TRUE)
  } else {
    resid <- stats::residuals(fit)
  }
  colData(sce)$md_rna_residual <- as.numeric(resid)
  group <- interaction(
    as.character(df$identity),
    if (!is.null(df$sample)) as.character(df$sample) else "all",
    drop = TRUE
  )
  colData(sce)$md_rna_residual_z <- .md_group_robust_z(resid, group, min_group_n = 20L)
  sce
}
