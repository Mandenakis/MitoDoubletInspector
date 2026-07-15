# Identity reference and coherence scoring ------------------------------------

#' Build de novo identity reference profiles
#'
#' Builds cluster-level mean log-expression profiles for cosine identity scoring
#' and mean per-cell count-proportion profiles for additive mixture modelling.
#' High-confidence reference cells are preferred: original doublet calls,
#' top doublet scores, high stress, and high mitochondrial cells are excluded
#' where possible, so that suspicious cells are scored against a stable
#' reference rather than against themselves.
#'
#' @param sce SingleCellExperiment.
#' @param cluster_col Cluster column. If NULL, common names are inferred.
#' @param sample_col Optional sample column.
#' @param protect_rare If TRUE, rare clusters are retained but flagged.
#' @param rare_cluster_min_n Threshold below which a cluster is flagged as rare.
#' @param min_cells_per_profile Minimum reference cells per profile.
#' @param verbose Logical.
#' @return SingleCellExperiment.
#' @export
md_build_identity_reference <- function(
  sce,
  cluster_col = NULL,
  sample_col = NULL,
  protect_rare = TRUE,
  rare_cluster_min_n = 50L,
  min_cells_per_profile = 10L,
  verbose = TRUE
) {
  cluster_col <- cluster_col %||% .md_infer_col(sce, c("seurat_clusters", "cluster", "clusters", "label", "celltype", "cell_type"))
  if (is.null(cluster_col) || !.md_col_exists(sce, cluster_col)) {
    .md_warn_once("no_cluster_col", "No cluster column supplied/found. Using one global identity profile; lineage conflict will be weak.")
    clusters <- factor(rep("global", ncol(sce)))
  } else {
    cluster_values <- as.character(colData(sce)[[cluster_col]])
    cluster_values[is.na(cluster_values) | cluster_values == ""] <- "unknown"
    clusters <- factor(cluster_values)
  }

  fs <- metadata(sce)$md_feature_spaces
  if (is.null(fs$identity)) {
    .md_warn_once("no_feature_space", "No identity feature space found; running md_define_feature_spaces() with human defaults.")
    sce <- md_define_feature_spaces(sce, species = "human")
    fs <- metadata(sce)$md_feature_spaces
  }
  features <- .md_subset_features(sce, fs$identity, min_features = 20L)
  logcounts <- .md_sparse(assay(sce, "logcounts"))[features, , drop = FALSE]
  counts <- .md_sparse(.md_get_counts(sce))[features, , drop = FALSE]
  cd <- .md_coldata_df(sce)
  sample_group <- if (!is.null(sample_col) && sample_col %in% colnames(cd)) cd[[sample_col]] else NULL

  good <- rep(TRUE, ncol(sce))
  if ("md_original_doublet_class" %in% colnames(cd)) {
    good <- good & !tolower(cd$md_original_doublet_class) %in% c("doublet", "multiplet")
  }
  if ("md_doublet_score_full" %in% colnames(cd)) {
    good <- good & .md_percentile(cd$md_doublet_score_full, sample_group) < 0.95
  }
  if ("md_percent_mt" %in% colnames(cd)) {
    q <- .md_safe_quantile(cd$md_percent_mt, 0.95, default = Inf)
    good <- good & cd$md_percent_mt <= q
  }
  if ("md_stress_score" %in% colnames(cd)) {
    q <- .md_safe_quantile(cd$md_stress_score, 0.95, default = Inf)
    good <- good & cd$md_stress_score <= q
  }

  clv <- levels(clusters)
  profiles <- matrix(0, nrow = length(features), ncol = length(clv), dimnames = list(features, clv))
  prop_profiles <- matrix(0, nrow = length(features), ncol = length(clv), dimnames = list(features, clv))
  n_ref <- stats::setNames(integer(length(clv)), clv)
  n_total <- stats::setNames(integer(length(clv)), clv)

  for (cl in clv) {
    idx_all <- which(clusters == cl)
    idx <- idx_all[good[idx_all]]
    if (length(idx) < min_cells_per_profile) idx <- idx_all
    n_ref[cl] <- length(idx)
    n_total[cl] <- length(idx_all)
    profiles[, cl] <- as.numeric(.md_safe_rowmeans(logcounts[, idx, drop = FALSE]))
    cell_prop <- .md_as_proportions(counts[, idx, drop = FALSE])
    p <- as.numeric(.md_safe_rowmeans(cell_prop))
    prop_profiles[, cl] <- p / max(sum(p), .Machine$double.eps)
  }

  rare <- n_total < rare_cluster_min_n
  if (isTRUE(protect_rare) && any(rare)) {
    .md_msg("Rare cluster protection active for: ", paste(names(n_total)[rare], collapse = ", "), verbose = verbose)
  }

  mdmeta <- metadata(sce)
  mdmeta$md_identity_profiles <- profiles
  mdmeta$md_identity_proportions <- prop_profiles
  mdmeta$md_identity_features <- features
  mdmeta$md_identity_reference_n <- n_ref
  mdmeta$md_identity_total_n <- n_total
  mdmeta$md_identity_rare <- rare
  mdmeta$md_cluster_col <- cluster_col
  metadata(sce) <- mdmeta

  colData(sce)$md_input_cluster <- as.character(clusters)
  colData(sce)$md_input_cluster_is_rare <- as.logical(rare[as.character(clusters)])
  sce
}

#' Score identity coherence against reference profiles
#'
#' Computes cosine similarity to each identity reference and stores top identity,
#' second identity, coherence, and lineage conflict.
#'
#' @param sce SingleCellExperiment with identity reference profiles.
#' @param chunk_size Number of cells to score per chunk.
#' @return SingleCellExperiment.
#' @export
md_score_identity_coherence <- function(sce, chunk_size = 5000L) {
  profiles <- metadata(sce)$md_identity_profiles
  features <- metadata(sce)$md_identity_features
  if (is.null(profiles) || is.null(features)) .md_stop("No identity profiles found. Run md_build_identity_reference() first.")
  features <- intersect(features, rownames(sce))
  X <- .md_sparse(assay(sce, "logcounts"))[features, , drop = FALSE]
  P <- profiles[features, , drop = FALSE]
  Pn <- .md_profile_norm(P)

  top <- second <- rep(NA_character_, ncol(sce))
  top_score <- second_score <- rep(NA_real_, ncol(sce))
  all_scores <- matrix(NA_real_, nrow = ncol(Pn), ncol = ncol(sce), dimnames = list(colnames(Pn), colnames(sce)))

  if (ncol(Pn) == 1L) {
    cn <- .md_col_norms(X)
    cn[!is.finite(cn) | cn == 0] <- 1
    Xin <- X %*% Matrix::Diagonal(x = 1 / cn)
    scores <- as.numeric(crossprod(Pn, Xin))
    all_scores[1, ] <- scores
    top[] <- colnames(Pn)[1]
    second[] <- NA_character_
    top_score[] <- scores
    second_score[] <- 0
    coherence <- top_score
    conflict <- rep(0, ncol(sce))

    colData(sce)$md_identity_top <- top
    colData(sce)$md_identity_second <- second
    colData(sce)$md_identity_top_score <- top_score
    colData(sce)$md_identity_second_score <- second_score
    colData(sce)$md_identity_coherence <- coherence
    colData(sce)$md_lineage_conflict <- conflict
    colData(sce)$md_identity_reference_quality <- top_score

    mdmeta <- metadata(sce)
    mdmeta$md_identity_score_matrix <- all_scores
    metadata(sce) <- mdmeta
    return(sce)
  }

  starts <- seq(1L, ncol(sce), by = chunk_size)
  for (st in starts) {
    en <- min(ncol(sce), st + chunk_size - 1L)
    idx <- st:en
    Xi <- X[, idx, drop = FALSE]
    cn <- .md_col_norms(Xi)
    cn[!is.finite(cn) | cn == 0] <- 1
    Xin <- Xi %*% Matrix::Diagonal(x = 1 / cn)
    scores <- as.matrix(crossprod(Pn, Xin))
    all_scores[, idx] <- scores
    ord <- apply(scores, 2, order, decreasing = TRUE)
    if (is.matrix(ord)) {
      top_idx <- ord[1, ]
      second_idx <- if (nrow(ord) >= 2L) ord[2, ] else ord[1, ]
    } else {
      top_idx <- ord[1]
      second_idx <- ord[min(2, length(ord))]
    }
    top[idx] <- rownames(scores)[top_idx]
    second[idx] <- rownames(scores)[second_idx]
    top_score[idx] <- scores[cbind(top_idx, seq_along(idx))]
    second_score[idx] <- scores[cbind(second_idx, seq_along(idx))]
  }

  coherence <- pmax(0, top_score - second_score)
  # A ratio alone calls a cell with two equally terrible fits "conflicted".
  # Weight the small top/second margin by the absolute second-reference fit so
  # genuine two-lineage cells rank high while globally poor-quality cells do not.
  conflict <- second_score * (1 - coherence)
  conflict[!is.finite(conflict)] <- 0
  conflict <- .md_clip(conflict, 0, 1)

  colData(sce)$md_identity_top <- top
  colData(sce)$md_identity_second <- second
  colData(sce)$md_identity_top_score <- top_score
  colData(sce)$md_identity_second_score <- second_score
  colData(sce)$md_identity_coherence <- coherence
  colData(sce)$md_lineage_conflict <- conflict
  colData(sce)$md_identity_reference_quality <- top_score

  mdmeta <- metadata(sce)
  mdmeta$md_identity_score_matrix <- all_scores
  metadata(sce) <- mdmeta
  sce
}
