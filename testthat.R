# Doublet score import and counterfactual rerun -------------------------------

#' Import or run doublet scores, with optional no-cycle counterfactual
#'
#' If scDblFinder is available and `run_scdblfinder = TRUE`, this runs
#' scDblFinder on all genes and again after removing the package's cycle feature
#' space. Otherwise it imports user-supplied score/class columns. The difference
#' between the two runs (`md_cycle_drop`) is a key piece of rescue evidence:
#' cells whose doublet score collapses once cell-cycle genes are removed were
#' suspicious largely because of proliferation, not because of two identities.
#'
#' @param sce SingleCellExperiment.
#' @param run_scdblfinder Logical.
#' @param doublet_score_col Existing score column to import.
#' @param doublet_class_col Existing class column to import.
#' @param sample_col Optional sample column.
#' @param cluster_col Optional cluster column passed to scDblFinder if present.
#' @param ncores Number of cores for BiocParallel if available.
#' @param seed Random seed shared by the full and no-cycle runs.
#' @param verbose Logical.
#' @return SingleCellExperiment.
#' @export
md_run_or_import_doublet_scores <- function(
  sce,
  run_scdblfinder = TRUE,
  doublet_score_col = NULL,
  doublet_class_col = NULL,
  sample_col = NULL,
  cluster_col = NULL,
  ncores = 1L,
  seed = 1234L,
  verbose = TRUE
) {
  cd <- .md_coldata_df(sce)

  if (!is.null(doublet_score_col)) {
    if (!doublet_score_col %in% colnames(cd)) .md_stop("doublet_score_col not found: ", doublet_score_col)
    imported <- as.numeric(cd[[doublet_score_col]])
    if (any(!is.finite(imported))) .md_stop("doublet_score_col contains non-finite values.")
    colData(sce)$md_doublet_score_full <- imported
  }
  if (!is.null(doublet_class_col)) {
    if (!doublet_class_col %in% colnames(cd)) .md_stop("doublet_class_col not found: ", doublet_class_col)
    colData(sce)$md_original_doublet_class <- as.character(cd[[doublet_class_col]])
  }

  can_run <- isTRUE(run_scdblfinder) && requireNamespace("scDblFinder", quietly = TRUE)
  if (can_run) {
    .md_msg("Running scDblFinder on full feature space...", verbose = verbose)
    full <- .md_try_scdblfinder(sce, sample_col, cluster_col, ncores, seed)
    if (!is.null(full)) {
      colData(sce)$md_doublet_score_full <- full$score
      colData(sce)$md_original_doublet_class <- full$class
    }

    fs <- metadata(sce)$md_feature_spaces
    no_cycle_features <- rownames(sce)
    if (!is.null(fs$cycle) && length(fs$cycle)) no_cycle_features <- setdiff(rownames(sce), fs$cycle)
    if (length(no_cycle_features) >= 200L) {
      .md_msg("Running scDblFinder after removing cycle/cytokinesis genes...", verbose = verbose)
      sce_no <- sce[no_cycle_features, ]
      no <- .md_try_scdblfinder(sce_no, sample_col, cluster_col, ncores, seed)
      if (!is.null(no)) {
        colData(sce)$md_doublet_score_no_cycle <- no$score
        colData(sce)$md_no_cycle_doublet_class <- no$class
      }
    } else {
      .md_warn_once("too_few_no_cycle", "Too few features for no-cycle scDblFinder rerun; skipping.")
    }
  } else if (isTRUE(run_scdblfinder)) {
    .md_warn_once(
      "scdblfinder_missing",
      "scDblFinder is not installed. Using imported doublet scores if supplied; otherwise doublet score columns are set to zero."
    )
  }

  if (!"md_doublet_score_full" %in% colnames(colData(sce))) {
    colData(sce)$md_doublet_score_full <- rep(0, ncol(sce))
  }
  if (!"md_original_doublet_class" %in% colnames(colData(sce))) {
    colData(sce)$md_original_doublet_class <- rep("unknown", ncol(sce))
  }
  if (!"md_doublet_score_no_cycle" %in% colnames(colData(sce))) {
    colData(sce)$md_doublet_score_no_cycle <- colData(sce)$md_doublet_score_full
  }
  if (!"md_no_cycle_doublet_class" %in% colnames(colData(sce))) {
    colData(sce)$md_no_cycle_doublet_class <- rep("unknown", ncol(sce))
  }

  full_score <- as.numeric(colData(sce)$md_doublet_score_full)
  no_score <- as.numeric(colData(sce)$md_doublet_score_no_cycle)
  group <- if (!is.null(sample_col) && sample_col %in% colnames(colData(sce))) colData(sce)[[sample_col]] else NULL
  full_rank <- .md_percentile(full_score, group)
  no_rank <- .md_percentile(no_score, group)
  colData(sce)$md_doublet_rank_full <- full_rank
  colData(sce)$md_doublet_rank_no_cycle <- no_rank
  colData(sce)$md_cycle_drop_raw <- full_score - no_score
  # scDblFinder scores from two separately trained models are not on a common
  # absolute scale. The primary counterfactual is therefore the within-sample
  # percentile drop; raw differences are retained only for audit.
  colData(sce)$md_cycle_drop <- pmax(0, full_rank - no_rank)
  counterfactual_available <- "md_no_cycle_doublet_class" %in% colnames(colData(sce)) &&
    any(as.character(colData(sce)$md_no_cycle_doublet_class) != "unknown")
  colData(sce)$md_counterfactual_available <- rep(counterfactual_available, ncol(sce))
  metadata(sce)$md_doublet_score_source <- if (can_run) "scDblFinder" else if (!is.null(doublet_score_col)) "imported" else "none"
  sce
}

.md_try_scdblfinder <- function(sce, sample_col, cluster_col, ncores, seed) {
  args <- list(sce = sce)
  if (!is.null(sample_col) && sample_col %in% colnames(colData(sce))) {
    args$samples <- sample_col
  }
  if (!is.null(cluster_col) && cluster_col %in% colnames(colData(sce))) {
    args$clusters <- cluster_col
  } else {
    args$clusters <- TRUE
  }
  if (requireNamespace("BiocParallel", quietly = TRUE) && ncores > 1L) {
    args$BPPARAM <- BiocParallel::MulticoreParam(workers = ncores)
  }
  out <- tryCatch(.md_with_seed(seed, {
    suppressMessages(do.call(scDblFinder::scDblFinder, args))
  }), error = function(e) {
    .md_warn_once(
      paste0("scdblfinder_error_", .md_key(conditionMessage(e))),
      "scDblFinder failed: ", conditionMessage(e), ". Continuing without this run."
    )
    NULL
  })
  if (is.null(out)) return(NULL)
  cd <- .md_coldata_df(out)
  if (!all(c("scDblFinder.score", "scDblFinder.class") %in% colnames(cd))) {
    .md_warn_once("scdblfinder_missing_columns", "scDblFinder output did not contain expected columns.")
    return(NULL)
  }
  list(score = as.numeric(cd$scDblFinder.score), class = as.character(cd$scDblFinder.class))
}
