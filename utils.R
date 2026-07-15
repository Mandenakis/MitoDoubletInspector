# Main workflow ----------------------------------------------------------------

#' Mitosis-aware doublet adjudication
#'
#' Runs the complete mitoDoubletR workflow: conversion to SingleCellExperiment,
#' QC, background-corrected cell-cycle scoring, stress/ambient scoring, feature
#' space construction, doublet-score import or scDblFinder calls, de novo
#' identity reference building, identity coherence scoring, RNA residual
#' modelling, calibrated parent-mixture fitting, and final classification.
#'
#' @param x Seurat object, SingleCellExperiment, or gene-by-cell matrix.
#' @param assay Seurat assay name if input is Seurat.
#' @param sample_col Optional sample/batch column in metadata.
#' @param cluster_col Optional cluster column in metadata.
#' @param species "human" or "mouse".
#' @param doublet_score_col Existing doublet score metadata column.
#' @param doublet_class_col Existing doublet class metadata column.
#' @param run_scdblfinder If TRUE, try to run scDblFinder. If unavailable,
#'   imported scores are used.
#' @param protect_rare If TRUE, rare clusters are protected from overconfident
#'   doublet calls unless evidence is strong.
#' @param rare_cluster_min_n Cluster size below which rare protection is applied.
#' @param n_hvg Number of variable genes for feature-space construction.
#' @param bg_correct If TRUE (default), background-correct module scores.
#' @param mixture_space "proportion" (default) or "log" for the parent-mixture
#'   model.
#' @param max_mixture_candidates Runtime cap for parent-mixture modelling.
#' @param mixture_min_alpha Minimum minority-parent contribution.
#' @param mixture_null_n Number of high-confidence singlets used for empirical
#'   mixture calibration.
#' @param seed Random seed used by stochastic components without altering the
#'   caller's RNG state.
#' @param ncores Number of cores for scDblFinder/BiocParallel when available.
#' @param return_sce If TRUE, always return a SingleCellExperiment. Otherwise
#'   Seurat inputs receive md_* metadata and are returned as Seurat objects.
#' @param verbose Logical.
#' @return Object of same type as input where practical, otherwise SCE.
#' @export
md_adjudicate <- function(
  x,
  assay = "RNA",
  sample_col = NULL,
  cluster_col = NULL,
  species = c("human", "mouse"),
  doublet_score_col = NULL,
  doublet_class_col = NULL,
  run_scdblfinder = TRUE,
  protect_rare = TRUE,
  rare_cluster_min_n = 50L,
  n_hvg = 3000L,
  bg_correct = TRUE,
  mixture_space = c("proportion", "log"),
  max_mixture_candidates = 10000L,
  mixture_min_alpha = 0.15,
  mixture_null_n = 500L,
  seed = 1234L,
  ncores = 1L,
  return_sce = FALSE,
  verbose = TRUE
) {
  species <- match.arg(species)
  mixture_space <- match.arg(mixture_space)
  original <- x

  .md_msg("Converting input to SingleCellExperiment...", verbose = verbose)
  sce <- md_as_sce(x, assay = assay)

  if (!is.null(sample_col) && !.md_col_exists(sce, sample_col)) {
    .md_stop("sample_col '", sample_col, "' not found in metadata.")
  }
  if (!is.null(cluster_col) && !.md_col_exists(sce, cluster_col)) {
    .md_stop("cluster_col '", cluster_col, "' not found in metadata.")
  }

  .md_msg("Adding QC metrics...", verbose = verbose)
  sce <- md_add_qc_metrics(sce, sample_col = sample_col)

  .md_msg("Scoring cell-cycle, mitotic, cytokinesis, stress and ambient programmes...", verbose = verbose)
  sce <- md_score_cell_cycle(sce, species = species, bg_correct = bg_correct)
  sce <- md_score_stress_ambient(sce, species = species, bg_correct = bg_correct)

  .md_msg("Defining full/identity/cycle feature spaces...", verbose = verbose)
  sce <- md_define_feature_spaces(sce, species = species, n_hvg = n_hvg)

  sce <- md_run_or_import_doublet_scores(
    sce,
    run_scdblfinder = run_scdblfinder,
    doublet_score_col = doublet_score_col,
    doublet_class_col = doublet_class_col,
    sample_col = sample_col,
    cluster_col = cluster_col,
    ncores = ncores,
    seed = seed,
    verbose = verbose
  )

  .md_msg("Building identity reference profiles...", verbose = verbose)
  sce <- md_build_identity_reference(
    sce,
    cluster_col = cluster_col,
    sample_col = sample_col,
    protect_rare = protect_rare,
    rare_cluster_min_n = rare_cluster_min_n,
    verbose = verbose
  )

  .md_msg("Scoring identity coherence and lineage conflict...", verbose = verbose)
  sce <- md_score_identity_coherence(sce)

  .md_msg("Modelling within-identity RNA-content residuals...", verbose = verbose)
  sce <- md_score_rna_residuals(sce, sample_col = sample_col)

  sce <- md_fit_parent_mixture_models(
    sce,
    candidates = "flagged_or_extreme",
    mixture_space = mixture_space,
    min_alpha = mixture_min_alpha,
    null_n = mixture_null_n,
    max_candidates = max_mixture_candidates,
    seed = seed,
    verbose = verbose
  )

  .md_msg("Classifying cells under competing explanations...", verbose = verbose)
  sce <- md_classify_cells(sce, sample_col = sample_col)

  if (isTRUE(return_sce)) return(sce)
  .md_return_object(original, sce)
}
