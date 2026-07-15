# Cell classification ----------------------------------------------------------

#' Classify cells by competing mitotic/doublet explanations
#'
#' Produces bounded evidence-support scores and final `md_class` calls. The
#' supports are transparent heuristic evidence summaries, not probabilities.
#' A doublet call requires concordance of at least two evidence families
#' (external caller, lineage conflict, empirical parent-mixture fit). This avoids
#' turning a single arbitrary threshold into a biological deletion decision.
#'
#' @param sce SingleCellExperiment.
#' @param sample_col Optional sample column for within-sample standardization.
#' @param ambiguity_margin If top support minus second support is below
#'   this value, label ambiguous.
#' @param min_top_support Minimum winning support before labelling
#'   ambiguous.
#' @return SingleCellExperiment.
#' @export
md_classify_cells <- function(
  sce,
  sample_col = NULL,
  ambiguity_margin = 0.08,
  min_top_support = 0.55
) {
  cd <- .md_coldata_df(sce)
  n <- nrow(cd)
  get <- function(nm, default = 0) {
    if (nm %in% colnames(cd)) return(cd[[nm]])
    if (length(default) == n) return(default)
    rep(default[1], n)
  }

  grp <- NULL
  if (!is.null(sample_col) && sample_col %in% colnames(cd)) grp <- cd[[sample_col]]

  Dfull <- as.numeric(get("md_doublet_score_full"))
  Dno <- as.numeric(get("md_doublet_score_no_cycle", Dfull))
  cycle_drop <- as.numeric(get("md_cycle_drop"))
  conflict <- as.numeric(get("md_lineage_conflict"))
  coherence <- as.numeric(get("md_identity_coherence"))
  mitosis <- as.numeric(get("md_mitosis_score"))
  cytokinesis <- as.numeric(get("md_cytokinesis_score"))
  rna_resid <- as.numeric(get("md_rna_residual"))
  stress <- as.numeric(get("md_stress_score"))
  ambient <- as.numeric(get("md_ambient_score"))
  mt <- as.numeric(get("md_percent_mt"))
  mix <- as.numeric(get("md_mixture_evidence", get("md_mixture_improvement")))
  original_flag <- tolower(as.character(get("md_original_doublet_class", "unknown"))) %in% c("doublet", "multiplet")

  pDfull <- if ("md_doublet_rank_full" %in% colnames(cd)) cd$md_doublet_rank_full else .md_percentile(Dfull, grp)
  pDno <- if ("md_doublet_rank_no_cycle" %in% colnames(cd)) cd$md_doublet_rank_no_cycle else .md_percentile(Dno, grp)
  has_cf <- as.logical(get("md_counterfactual_available", FALSE))
  pdrop <- ifelse(has_cf, .md_percentile(cycle_drop, grp), 0)
  pconf <- .md_percentile(conflict, grp)
  pcoh <- .md_percentile(coherence, grp)
  pmit <- .md_percentile(mitosis, grp)
  pcyto <- .md_percentile(cytokinesis, grp)
  prna <- .md_percentile(rna_resid, grp)
  pstress <- .md_percentile(stress, grp)
  pambient <- .md_percentile(ambient, grp)
  pmt <- .md_percentile(mt, grp)
  pmix <- if ("md_mixture_p_value" %in% colnames(cd)) {
    1 - .md_clip(as.numeric(cd$md_mixture_p_value), 0, 1)
  } else {
    .md_percentile(mix, grp)
  }

  supports <- cbind(
    singlet = 0.30 * pcoh + 0.25 * (1 - pconf) + 0.25 * (1 - pDno) + 0.20 * (1 - pmix),
    doublet = 0.30 * pDno + 0.25 * pconf + 0.35 * pmix + 0.10 * as.numeric(original_flag),
    mitotic_singlet = (0.30 * pmit + 0.20 * pcyto + 0.20 * pcoh + 0.15 * (1 - pconf) + 0.15 * pdrop) * (1 - 0.60 * pmix),
    polyploid_like_singlet = (0.40 * prna + 0.25 * pcoh + 0.20 * (1 - pconf) + 0.15 * (1 - pmix)) - 0.10 * pstress - 0.10 * pambient,
    stress_or_ambient = 0.45 * pstress + 0.25 * pambient + 0.20 * pmt + 0.10 * (1 - pcoh)
  )
  supports <- .md_clip(supports, 0, 1)
  ord <- t(apply(supports, 1, order, decreasing = TRUE))
  top_idx <- ord[, 1]
  second_idx <- ord[, 2]
  classes <- colnames(supports)[top_idx]
  top_support <- supports[cbind(seq_len(n), top_idx)]
  second_support <- supports[cbind(seq_len(n), second_idx)]
  confidence <- top_support - second_support

  caller_evidence <- pDno >= 0.90 | original_flag
  conflict_evidence <- pconf >= 0.80
  mixture_evidence <- pmix >= 0.90
  doublet_votes <- caller_evidence + conflict_evidence + mixture_evidence
  classes[classes == "doublet" & doublet_votes < 2L] <- "ambiguous"
  classes[classes == "mitotic_singlet" & !(pmit >= 0.80 | pcyto >= 0.80)] <- "ambiguous"
  classes[classes == "mitotic_singlet" & (pcoh < 0.50 | pmix >= 0.90)] <- "ambiguous"
  classes[classes == "polyploid_like_singlet" & (prna < 0.80 | pcoh < 0.50 | pmix >= 0.90)] <- "ambiguous"
  classes[classes == "stress_or_ambient" & pmax(pstress, pambient) < 0.80] <- "ambiguous"
  classes[confidence < ambiguity_margin | top_support < min_top_support] <- "ambiguous"

  # Rare clusters are not assumed to be singlets. Instead, a rare-cluster
  # doublet without strong mixture evidence is explicitly left ambiguous.
  if ("md_input_cluster_is_rare" %in% colnames(cd)) {
    rare <- as.logical(cd$md_input_cluster_is_rare)
    classes[rare & classes == "doublet" & !mixture_evidence] <- "ambiguous"
  }

  colData(sce)$md_support_singlet <- supports[, "singlet"]
  colData(sce)$md_support_doublet <- supports[, "doublet"]
  colData(sce)$md_support_mitotic <- supports[, "mitotic_singlet"]
  colData(sce)$md_support_polyploid_like <- supports[, "polyploid_like_singlet"]
  colData(sce)$md_support_stress_ambient <- supports[, "stress_or_ambient"]
  # Compatibility aliases. These are support scores, not calibrated probabilities.
  colData(sce)$md_score_singlet <- supports[, "singlet"]
  colData(sce)$md_score_doublet <- supports[, "doublet"]
  colData(sce)$md_score_mitotic <- supports[, "mitotic_singlet"]
  colData(sce)$md_score_polyploid_like <- supports[, "polyploid_like_singlet"]
  colData(sce)$md_score_stress_ambient <- supports[, "stress_or_ambient"]
  colData(sce)$md_class <- classes
  colData(sce)$md_confidence <- confidence
  colData(sce)$md_top_support <- top_support
  colData(sce)$md_doublet_evidence_votes <- as.integer(doublet_votes)
  tier <- ifelse(classes == "ambiguous", "ambiguous", ifelse(confidence >= 0.20, "high", "supported"))
  colData(sce)$md_decision_tier <- tier

  # Identifiability: RNA-only data cannot separate a high-RNA/polyploid-like
  # singlet from a homotypic (same-identity) doublet. Flag the cells where that
  # ambiguity is live so downstream users do not over-interpret the call.
  ident <- rep("resolved", n)
  ident[classes == "mitotic_singlet"] <- "mitotic_transcriptome_supported_not_proof_of_division"
  ident[classes == "polyploid_like_singlet" | (prna >= 0.80 & pmix < 0.90 & pcoh >= 0.50)] <- "high_rna_singlet_or_homotypic_doublet_unresolved"
  if (ncol(metadata(sce)$md_identity_profiles %||% matrix(nrow = 0, ncol = 0)) < 2L) {
    ident[] <- "insufficient_identity_reference"
  }
  colData(sce)$md_identifiability <- ident

  rescue_status <- rep("not_rescue", n)
  rescue_status[original_flag & classes == "mitotic_singlet"] <- "candidate_mitotic_rescue"
  rescue_status[original_flag & classes == "polyploid_like_singlet"] <- "candidate_high_rna_rescue_unresolved"
  colData(sce)$md_rescue_status <- rescue_status
  colData(sce)$md_rescue_eligible <- rescue_status != "not_rescue"

  colData(sce)$md_reason <- .md_reason_vector(sce)
  sce
}

.md_reason_vector <- function(sce) {
  cd <- .md_coldata_df(sce)
  n <- nrow(cd)
  q <- function(nm, p) .md_safe_quantile(cd[[nm]], p, default = Inf)
  hi_mit <- if ("md_mitosis_score" %in% colnames(cd)) cd$md_mitosis_score >= q("md_mitosis_score", 0.80) else rep(FALSE, n)
  hi_cyto <- if ("md_cytokinesis_score" %in% colnames(cd)) cd$md_cytokinesis_score >= q("md_cytokinesis_score", 0.80) else rep(FALSE, n)
  hi_conflict <- if ("md_lineage_conflict" %in% colnames(cd)) cd$md_lineage_conflict >= q("md_lineage_conflict", 0.80) else rep(FALSE, n)
  low_conflict <- if ("md_lineage_conflict" %in% colnames(cd)) cd$md_lineage_conflict <= q("md_lineage_conflict", 0.40) else rep(FALSE, n)
  hi_coh <- if ("md_identity_coherence" %in% colnames(cd)) cd$md_identity_coherence >= q("md_identity_coherence", 0.70) else rep(FALSE, n)
  hi_drop <- if ("md_cycle_drop" %in% colnames(cd)) cd$md_cycle_drop >= q("md_cycle_drop", 0.75) else rep(FALSE, n)
  hi_mix <- if ("md_mixture_improvement" %in% colnames(cd)) cd$md_mixture_improvement >= q("md_mixture_improvement", 0.80) else rep(FALSE, n)
  hi_rna <- if ("md_rna_residual" %in% colnames(cd)) cd$md_rna_residual >= q("md_rna_residual", 0.80) else rep(FALSE, n)
  hi_stress <- if ("md_stress_score" %in% colnames(cd)) cd$md_stress_score >= q("md_stress_score", 0.85) else rep(FALSE, n)
  hi_ambient <- if ("md_ambient_score" %in% colnames(cd)) cd$md_ambient_score >= q("md_ambient_score", 0.85) else rep(FALSE, n)
  ident <- if ("md_identifiability" %in% colnames(cd)) as.character(cd$md_identifiability) else rep("resolved", n)

  reasons <- character(n)
  for (i in seq_len(n)) {
    cls <- cd$md_class[i]
    bits <- character()
    if (isTRUE(hi_mit[i])) bits <- c(bits, "high mitotic programme")
    if (isTRUE(hi_cyto[i])) bits <- c(bits, "high cytokinesis programme")
    if (isTRUE(hi_drop[i])) bits <- c(bits, "doublet score drops after cycle-gene removal")
    if (isTRUE(hi_coh[i])) bits <- c(bits, "coherent dominant identity")
    if (isTRUE(low_conflict[i])) bits <- c(bits, "low lineage conflict")
    if (isTRUE(hi_conflict[i])) bits <- c(bits, "high lineage conflict")
    if (isTRUE(hi_mix[i])) bits <- c(bits, "two-parent mixture improves fit")
    if (isTRUE(hi_rna[i])) bits <- c(bits, "high RNA-content residual")
    if (isTRUE(hi_stress[i])) bits <- c(bits, "stress/dissociation programme")
    if (isTRUE(hi_ambient[i])) bits <- c(bits, "ambient-prone RNA signal")
    if (!length(bits)) bits <- "weak or balanced evidence"
    tail <- switch(
      ident[i],
      high_rna_singlet_or_homotypic_doublet_unresolved = " [RNA cannot exclude a homotypic doublet]",
      mitotic_transcriptome_supported_not_proof_of_division = " [RNA supports a mitotic programme but cannot prove cell division]",
      insufficient_identity_reference = " [fewer than two usable identity references]",
      ""
    )
    reasons[i] <- paste0(cls, ": ", paste(bits, collapse = "; "), ".", tail)
  }
  reasons
}

#' Return compact reason table
#'
#' @param sce SingleCellExperiment after classification.
#' @return Data.frame with key audit columns.
#' @export
md_reason_table <- function(sce) {
  cd <- .md_coldata_df(sce)
  wanted <- c(
    "md_class", "md_decision_tier", "md_confidence", "md_reason", "md_identifiability",
    "md_rescue_status", "md_rescue_eligible",
    "md_identity_top", "md_identity_second",
    "md_doublet_score_full", "md_doublet_score_no_cycle", "md_cycle_drop",
    "md_mitosis_score", "md_cytokinesis_score", "md_lineage_conflict",
    "md_identity_coherence", "md_rna_residual", "md_mixture_improvement",
    "md_mixture_p_value", "md_mixture_q_value", "md_doublet_evidence_votes",
    "md_parent_A", "md_parent_B", "md_parent_alpha"
  )
  out <- cd[, intersect(wanted, colnames(cd)), drop = FALSE]
  out$cell <- rownames(cd)
  out <- out[, c("cell", setdiff(colnames(out), "cell")), drop = FALSE]
  out
}

#' Summarise mitoDoubletR classes
#'
#' @param sce SingleCellExperiment or Seurat object carrying md_class metadata.
#' @return Data.frame.
#' @export
md_summary <- function(sce) {
  if (.md_is_seurat(sce)) {
    cd <- sce[[]]
  } else {
    cd <- .md_coldata_df(sce)
  }
  if (!"md_class" %in% colnames(cd)) .md_stop("No md_class column found. Run md_classify_cells() or md_adjudicate().")
  tab <- as.data.frame(table(md_class = cd$md_class), stringsAsFactors = FALSE)
  tab$percent <- 100 * tab$Freq / sum(tab$Freq)
  tab
}
