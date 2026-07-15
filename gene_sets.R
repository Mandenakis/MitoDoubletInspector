# Feature space construction ---------------------------------------------------

#' Define full, identity, and cycle feature spaces
#'
#' Creates three feature spaces stored in `metadata(sce)$md_feature_spaces`:
#' full HVGs, identity HVGs with cycle/stress/ambient features removed, and
#' cycle/cytokinesis features. This three-space design is the core of the
#' method: identity is judged in a cycle-neutral space, while suspiciousness that
#' is driven only by cell-cycle genes can be detected and discounted.
#'
#' @param sce SingleCellExperiment.
#' @param species "human" or "mouse".
#' @param n_hvg Number of variable features to retain for the full feature space.
#' @param gene_sets Optional custom gene-set list.
#' @return SingleCellExperiment.
#' @export
md_define_feature_spaces <- function(
  sce,
  species = c("human", "mouse"),
  n_hvg = 3000L,
  gene_sets = NULL
) {
  species <- match.arg(species)
  gene_sets <- gene_sets %||% md_default_gene_sets(species)
  sce <- .md_add_logcounts(sce)
  logcounts <- .md_sparse(assay(sce, "logcounts"))
  vars <- .md_row_variance_sparse(logcounts)
  names(vars) <- rownames(sce)
  vars[!is.finite(vars)] <- 0
  keep <- names(sort(vars, decreasing = TRUE))[seq_len(min(n_hvg, length(vars)))]

  cycle <- unique(c(
    .md_match_genes(sce, gene_sets$S),
    .md_match_genes(sce, gene_sets$G2M),
    .md_match_genes(sce, gene_sets$mitosis),
    .md_match_genes(sce, gene_sets$cytokinesis)
  ))
  stress <- unique(c(
    .md_match_genes(sce, gene_sets$stress),
    .md_match_genes(sce, gene_sets$ambient)
  ))
  ribo <- .md_gene_pattern(sce, "^(RPL|RPS|Rpl|Rps)")
  mito <- .md_gene_pattern(sce, "^(MT-|mt-)")
  hb <- .md_match_genes(sce, c("HBB", "HBA1", "HBA2", "HBD", "HBE1", "Hbb", "Hba-a1", "Hba-a2"))
  blacklist <- unique(c(cycle, stress, ribo, mito, hb))

  identity <- setdiff(keep, blacklist)
  if (length(identity) < 200L) {
    .md_warn_once(
      "few_identity_features",
      "Fewer than 200 identity features remain after filtering; relaxing to HVGs minus cycle genes only."
    )
    identity <- setdiff(keep, cycle)
  }

  mdmeta <- metadata(sce)
  mdmeta$md_feature_spaces <- list(
    full = keep,
    identity = identity,
    cycle = cycle,
    blacklist = blacklist
  )
  metadata(sce) <- mdmeta
  sce
}

