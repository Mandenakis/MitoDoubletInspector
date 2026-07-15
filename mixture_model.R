# Benchmarking harness ---------------------------------------------------------

#' Benchmark adjudication against simulated ground truth
#'
#' Simulates a labelled stress-test dataset (or uses one supplied via
#' `truth_col`), runs \code{\link{md_adjudicate}}, and scores the result. This
#' operationalizes the validation the method depends on: it reports the
#' confusion matrix and the quantities that actually matter for rare-cell rescue
#' - true-doublet sensitivity, and the false-removal rate of genuine
#' mitotic/polyploid-like biology.
#'
#' @param sce Input SingleCellExperiment to simulate from, or an object that
#'   already carries a truth column named by `truth_col`.
#' @param cluster_col Cluster column for simulation and adjudication.
#' @param sample_col Optional sample/batch column. Synthetic doublet parents are
#'   drawn within sample, and adjudication percentiles are sample-aware.
#' @param truth_col Name of an existing ground-truth column. If present, no new
#'   simulation is performed.
#' @param n_doublets,n_mitotic,n_polyploid_like,n_ambient Simulation sizes.
#' @param n_homotypic Number of homotypic synthetic doublets.
#' @param n_stress Number of stress-programme cells.
#' @param species "human" or "mouse".
#' @param run_scdblfinder Passed to md_adjudicate. Independent evaluation of a
#'   second-stage adjudicator should normally include the upstream caller.
#' @param seed Random seed.
#' @param return_object If TRUE, include the adjudicated SCE in the result.
#' @param ... Further arguments forwarded to \code{\link{md_adjudicate}}.
#' @return A list with `confusion`, `metrics`, `per_class`, and optionally
#'   `object`.
#' @export
md_benchmark_adjudication <- function(
  sce,
  cluster_col = NULL,
  sample_col = NULL,
  truth_col = "truth_class",
  n_doublets = 200L,
  n_homotypic = 100L,
  n_mitotic = 200L,
  n_polyploid_like = 200L,
  n_ambient = 200L,
  n_stress = 200L,
  species = c("human", "mouse"),
  run_scdblfinder = TRUE,
  seed = 1L,
  return_object = FALSE,
  ...
) {
  species <- match.arg(species)

  has_truth <- .md_is_sce(sce) && truth_col %in% colnames(colData(sce))
  if (has_truth) {
    sim <- sce
  } else {
    sim <- md_simulate_adjudication_dataset(
      sce, cluster_col = cluster_col, sample_col = sample_col,
      n_doublets = n_doublets, n_homotypic = n_homotypic,
      n_mitotic = n_mitotic,
      n_polyploid_like = n_polyploid_like, n_ambient = n_ambient,
      n_stress = n_stress,
      species = species, seed = seed
    )
    truth_col <- "truth_class"
  }

  sim_cluster_col <- if ("sim_cluster" %in% colnames(colData(sim))) "sim_cluster" else cluster_col
  adj <- md_adjudicate(
    sim, cluster_col = sim_cluster_col, sample_col = sample_col, species = species,
    run_scdblfinder = run_scdblfinder, return_sce = TRUE, verbose = FALSE, ...
  )

  cd <- .md_coldata_df(adj)
  truth <- as.character(cd[[truth_col]])
  pred <- as.character(cd$md_class)

  # Map simulator truth onto the adjudicator's class vocabulary.
  truth_mapped <- truth
  truth_mapped[truth_mapped == "original"] <- "singlet"
  singlet_family <- c("singlet", "mitotic_singlet", "polyploid_like_singlet")

  confusion <- table(truth = truth_mapped, predicted = pred)

  frac <- function(mask_num, mask_den) {
    den <- sum(mask_den)
    if (den == 0) return(NA_real_)
    sum(mask_num & mask_den) / den
  }

  is_truth_doublet <- truth_mapped == "doublet"
  truth_subclass <- if ("truth_subclass" %in% colnames(cd)) as.character(cd$truth_subclass) else truth_mapped
  is_truth_mito <- truth_mapped == "mitotic_singlet"
  is_truth_poly <- truth_mapped == "polyploid_like_singlet"
  is_truth_singletfam <- truth_mapped %in% singlet_family
  is_truth_biology <- is_truth_mito | is_truth_poly

  metrics <- c(
    doublet_sensitivity = frac(pred == "doublet", is_truth_doublet),
    heterotypic_doublet_sensitivity = frac(pred == "doublet", truth_subclass == "heterotypic_doublet"),
    homotypic_doublet_sensitivity = frac(pred == "doublet", truth_subclass == "homotypic_doublet"),
    doublet_false_positive_rate = frac(pred == "doublet", is_truth_singletfam),
    mitotic_recovered_as_singletfamily = frac(pred %in% singlet_family, is_truth_mito),
    mitotic_called_mitotic = frac(pred == "mitotic_singlet", is_truth_mito),
    polyploid_recovered_as_singletfamily = frac(pred %in% singlet_family, is_truth_poly),
    false_removal_rate_of_biology = frac(pred == "doublet", is_truth_biology),
    false_rescue_rate_of_doublets = frac(pred %in% singlet_family, is_truth_doublet),
    ambiguous_fraction = mean(pred == "ambiguous"),
    high_rna_marked_unresolved = frac(cd$md_identifiability == "high_rna_singlet_or_homotypic_doublet_unresolved", is_truth_poly)
  )

  if ("md_support_doublet" %in% colnames(cd)) {
    y <- is_truth_doublet
    s <- as.numeric(cd$md_support_doublet)
    n_pos <- sum(y); n_neg <- sum(!y)
    metrics["doublet_auroc"] <- if (n_pos && n_neg) {
      (sum(rank(s, ties.method = "average")[y]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    } else NA_real_
    ord_s <- order(s, decreasing = TRUE)
    y_ord <- y[ord_s]
    precision_at <- cumsum(y_ord) / seq_along(y_ord)
    metrics["doublet_average_precision"] <- if (n_pos) sum(precision_at[y_ord]) / n_pos else NA_real_
  }

  classes <- sort(unique(c(truth_mapped, pred)))
  per_class <- data.frame(
    class = classes,
    sensitivity = NA_real_,
    specificity = NA_real_,
    n_truth = NA_integer_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(classes)) {
    cl <- classes[i]
    tp <- sum(pred == cl & truth_mapped == cl)
    fn <- sum(pred != cl & truth_mapped == cl)
    fp <- sum(pred == cl & truth_mapped != cl)
    tn <- sum(pred != cl & truth_mapped != cl)
    per_class$sensitivity[i] <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    per_class$specificity[i] <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    per_class$n_truth[i] <- tp + fn
  }

  subtype_confusion <- table(truth_subclass = truth_subclass, predicted = pred)
  out <- list(confusion = confusion, subtype_confusion = subtype_confusion, metrics = metrics, per_class = per_class)
  if (isTRUE(return_object)) out$object <- adj
  out
}
