# QC and module scoring --------------------------------------------------------

#' Add basic QC metrics
#'
#' Adds `md_nCount`, `md_nFeature`, and `md_percent_mt`.
#'
#' @param sce SingleCellExperiment.
#' @param sample_col Optional sample/batch column.
#' @return SingleCellExperiment.
#' @export
md_add_qc_metrics <- function(sce, sample_col = NULL) {
  counts <- .md_sparse(.md_get_counts(sce))
  nCount <- Matrix::colSums(counts)
  nFeature <- Matrix::colSums(counts > 0)
  mt_features <- .md_gene_pattern(sce, "^(MT-|mt-)")
  pct_mt <- rep(0, ncol(sce))
  if (length(mt_features)) {
    mt_counts <- Matrix::colSums(counts[mt_features, , drop = FALSE])
    pct_mt <- 100 * mt_counts / pmax(nCount, 1)
  }
  colData(sce)$md_nCount <- as.numeric(nCount)
  colData(sce)$md_nFeature <- as.numeric(nFeature)
  colData(sce)$md_percent_mt <- as.numeric(pct_mt)

  if (!is.null(sample_col) && !.md_col_exists(sce, sample_col)) {
    .md_stop("sample_col '", sample_col, "' was not found in colData.")
  }
  sce
}

# Background-corrected module score (Tirosh 2016 / Seurat AddModuleScore style).
#
# Raw mean log-expression of a gene set is confounded by overall library size and
# by a cell's average expression level: high-RNA and polyploid-like cells score
# high on *every* module, including mitosis. That confound is fatal here because
# separating mitotic cells from high-RNA/polyploid-like cells is the package's
# central task. We therefore subtract, for each cell, the mean expression of a
# set of control genes drawn from the same average-expression bins as the
# signature genes. The difference isolates *specific* enrichment.
.md_module_score <- function(
  sce, genes, score_name,
  assay_name = "logcounts", min_genes = 3L,
  bg_correct = TRUE, ctrl = 100L, nbin = 24L, seed = 1234L
) {
  features <- .md_match_genes(sce, genes)
  mat <- .md_get_assay(sce, assay_name)
  ng <- nrow(mat)

  if (length(features) < min_genes) {
    .md_warn_once(
      paste0("few_genes_", score_name),
      "Only ", length(features), " genes matched for ", score_name,
      "; score will be set to zero."
    )
    score <- rep(0, ncol(sce))
  } else {
    sig_idx <- match(features, rownames(mat))
    sig_score <- as.numeric(.md_safe_colmeans(mat[sig_idx, , drop = FALSE]))

    do_bg <- isTRUE(bg_correct) && ng >= (min_genes + 10L)
    if (do_bg) {
      nbin_eff <- max(2L, min(as.integer(nbin), floor(ng / 5L)))
      avg <- as.numeric(.md_safe_rowmeans(mat))
      bins <- cut(
        rank(avg, ties.method = "first"),
        breaks = nbin_eff, labels = FALSE, include.lowest = TRUE
      )
      # Match controls to every signature gene's expression bin. This is closer
      # to AddModuleScore than sampling a fixed number per occupied bin and does
      # not let a densely populated bin dominate the background.
      ctrl_idx <- .md_with_seed(seed, {
        unlist(lapply(sig_idx, function(j) {
          pool <- setdiff(which(bins == bins[j]), sig_idx)
          if (!length(pool)) return(integer())
          take <- min(as.integer(ctrl), length(pool))
          pool[sample.int(length(pool), take)]
        }), use.names = FALSE)
      })
      ctrl_idx <- unique(ctrl_idx)
      if (length(ctrl_idx) >= min_genes) {
        ctrl_score <- as.numeric(.md_safe_colmeans(mat[ctrl_idx, , drop = FALSE]))
        score <- sig_score - ctrl_score
      } else {
        score <- sig_score
      }
    } else {
      score <- sig_score
    }
  }

  colData(sce)[[score_name]] <- score
  mdmeta <- metadata(sce)
  mdmeta$md_matched_genes <- mdmeta$md_matched_genes %||% list()
  mdmeta$md_matched_genes[[score_name]] <- features
  mdmeta$md_gene_set_coverage <- mdmeta$md_gene_set_coverage %||% list()
  mdmeta$md_gene_set_coverage[[score_name]] <- c(
    matched = length(features), requested = length(unique(genes)),
    fraction = length(features) / max(1L, length(unique(genes)))
  )
  metadata(sce) <- mdmeta
  sce
}

#' Score cell-cycle, mitotic, and cytokinesis programmes
#'
#' Uses background-corrected module scores so that overall RNA content does not
#' inflate the mitotic signal (see Details).
#'
#' @param sce SingleCellExperiment.
#' @param species "human" or "mouse".
#' @param gene_sets Optional custom gene-set list with names S, G2M, mitosis,
#'   cytokinesis.
#' @param bg_correct If TRUE (default), subtract a matched control-gene
#'   background from each module score.
#' @return SingleCellExperiment.
#' @export
md_score_cell_cycle <- function(sce, species = c("human", "mouse"), gene_sets = NULL, bg_correct = TRUE) {
  species <- match.arg(species)
  gene_sets <- gene_sets %||% md_default_gene_sets(species)
  sce <- .md_add_logcounts(sce)
  sce <- .md_module_score(sce, gene_sets$S, "md_s_score", bg_correct = bg_correct)
  sce <- .md_module_score(sce, gene_sets$G2M, "md_g2m_score", bg_correct = bg_correct)
  sce <- .md_module_score(sce, gene_sets$mitosis, "md_mitosis_score", bg_correct = bg_correct)
  sce <- .md_module_score(sce, gene_sets$cytokinesis, "md_cytokinesis_score", bg_correct = bg_correct)
  cd <- .md_coldata_df(sce)
  colData(sce)$md_cycle_score <- pmax(
    cd$md_s_score, cd$md_g2m_score, cd$md_mitosis_score, cd$md_cytokinesis_score,
    na.rm = TRUE
  )
  sce
}

#' Score stress and ambient-prone modules
#'
#' @param sce SingleCellExperiment.
#' @param species "human" or "mouse".
#' @param gene_sets Optional custom gene-set list.
#' @param bg_correct If TRUE (default), background-correct the stress score.
#' @return SingleCellExperiment.
#' @export
md_score_stress_ambient <- function(sce, species = c("human", "mouse"), gene_sets = NULL, bg_correct = TRUE) {
  species <- match.arg(species)
  gene_sets <- gene_sets %||% md_default_gene_sets(species)
  sce <- .md_add_logcounts(sce)
  sce <- .md_module_score(sce, gene_sets$stress, "md_stress_score", bg_correct = bg_correct)
  # Ambient is a contamination fraction, not a coherent programme: keep it as a
  # raw mean so it tracks absolute ambient-gene load rather than enrichment.
  sce <- .md_module_score(sce, gene_sets$ambient, "md_ambient_score", min_genes = 1L, bg_correct = FALSE)
  sce
}
