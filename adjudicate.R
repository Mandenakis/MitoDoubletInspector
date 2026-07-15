# Simulation engine ------------------------------------------------------------

#' Simulate adjudication stress-test data from an existing SCE
#'
#' Generates heterotypic and homotypic doublets with capture loss, heterogeneous
#' mitotic-like programmes, high-RNA/polyploid-like cells resampled from donor
#' proportions, stress states, and ambient contamination. This is a stress test,
#' not independent biological validation.
#'
#' @param sce Input SingleCellExperiment.
#' @param cluster_col Cluster column used to pick heterotypic doublets.
#' @param sample_col Optional sample/batch column. Both parents of a synthetic
#'   doublet are drawn from the same sample.
#' @param n_doublets Number of artificial doublets.
#' @param n_homotypic Number of artificial homotypic doublets.
#' @param n_mitotic Number of mitotic-like cells.
#' @param n_polyploid_like Number of polyploid-like cells.
#' @param n_ambient Number of ambient-contaminated cells.
#' @param n_stress Number of dissociation/stress-like cells.
#' @param doublet_capture_range Range of molecule-retention probabilities after
#'   summing two parent libraries.
#' @param species "human" or "mouse".
#' @param seed Random seed.
#' @return SingleCellExperiment with `truth_class` in colData.
#' @export
md_simulate_adjudication_dataset <- function(
  sce,
  cluster_col = NULL,
  sample_col = NULL,
  n_doublets = 200L,
  n_homotypic = 100L,
  n_mitotic = 200L,
  n_polyploid_like = 200L,
  n_ambient = 200L,
  n_stress = 200L,
  doublet_capture_range = c(0.55, 0.90),
  species = c("human", "mouse"),
  seed = 1L
) {
  species <- match.arg(species)
  if (length(doublet_capture_range) != 2L || any(!is.finite(doublet_capture_range)) ||
      doublet_capture_range[1] <= 0 || doublet_capture_range[2] > 1 ||
      doublet_capture_range[1] > doublet_capture_range[2]) {
    .md_stop("doublet_capture_range must be two increasing values in (0, 1].")
  }
  requested_n <- c(n_doublets, n_homotypic, n_mitotic, n_polyploid_like, n_ambient, n_stress)
  if (length(requested_n) != 6L || any(!is.finite(requested_n)) ||
      any(requested_n < 0) || any(abs(requested_n - round(requested_n)) > 0)) {
    .md_stop("All requested simulation counts must be non-negative finite integers.")
  }
  .md_with_seed(seed, {
  sce <- md_as_sce(sce)
  counts <- .md_sparse(.md_get_counts(sce))
  n <- ncol(sce)
  if (n < 4L) .md_stop("Need at least four cells to simulate from.")

  cluster_col <- cluster_col %||% .md_infer_col(sce, c("seurat_clusters", "cluster", "clusters", "celltype", "cell_type"))
  if (!is.null(sample_col) && !.md_col_exists(sce, sample_col)) {
    .md_stop("sample_col '", sample_col, "' was not found in colData.")
  }
  clusters <- if (!is.null(cluster_col) && .md_col_exists(sce, cluster_col)) {
    values <- as.character(colData(sce)[[cluster_col]])
    values[is.na(values) | values == ""] <- "unknown"
    factor(values)
  } else {
    factor(rep("global", n))
  }
  samples <- if (!is.null(sample_col) && .md_col_exists(sce, sample_col)) {
    values <- as.character(colData(sce)[[sample_col]])
    values[is.na(values) | values == ""] <- "unknown"
    factor(values)
  } else {
    factor(rep("all", n))
  }
  sample_values <- function(x, size = 1L, replace = FALSE) {
    if (size <= 0L) return(x[integer()])
    x[sample.int(length(x), size = size, replace = replace)]
  }

  sim_counts <- list(original = counts)
  truth <- rep("original", n)
  truth_subclass <- rep("original_singlet", n)
  sim_cluster <- as.character(clusters)
  sim_sample <- as.character(samples)
  parent_a <- parent_b <- rep(NA_character_, n)
  simulation_scale <- rep(1, n)

  thin_column <- function(v, probability) {
    v <- .md_sparse(v)
    if (length(v@x)) v@x <- stats::rbinom(length(v@x), size = round(v@x), prob = probability)
    v
  }

  # Heterotypic doublets where possible.
  if (n_doublets > 0L) {
    eligible <- which(vapply(seq_len(n), function(a) {
      any(samples == samples[a] & clusters != clusters[a])
    }, logical(1)))
    if (!length(eligible)) {
      .md_stop("No within-sample heterotypic parent pairs are available. Supply a valid cluster_col or set n_doublets = 0.")
    }
    dmat <- Matrix::Matrix(0, nrow = nrow(counts), ncol = n_doublets, sparse = TRUE)
    dclust <- character(n_doublets)
    dpa <- dpb <- dsample <- character(n_doublets)
    deff <- numeric(n_doublets)
    for (i in seq_len(n_doublets)) {
      a <- sample_values(eligible)
      possible <- which(samples == samples[a] & clusters != clusters[a])
      b <- sample_values(possible)
      eff <- stats::runif(1, doublet_capture_range[1], doublet_capture_range[2])
      dmat[, i] <- thin_column(counts[, a] + counts[, b], eff)
      dclust[i] <- as.character(clusters[a])
      dsample[i] <- as.character(samples[a])
      dpa[i] <- colnames(sce)[a]; dpb[i] <- colnames(sce)[b]; deff[i] <- eff
    }
    colnames(dmat) <- paste0("sim_doublet_", seq_len(n_doublets))
    sim_counts$doublet <- dmat
    truth <- c(truth, rep("doublet", n_doublets))
    truth_subclass <- c(truth_subclass, rep("heterotypic_doublet", n_doublets))
    sim_cluster <- c(sim_cluster, dclust)
    sim_sample <- c(sim_sample, dsample)
    parent_a <- c(parent_a, dpa); parent_b <- c(parent_b, dpb)
    simulation_scale <- c(simulation_scale, deff)
  }

  if (n_homotypic > 0L) {
    eligible <- which(vapply(seq_len(n), function(a) {
      any(seq_len(n) != a & samples == samples[a] & clusters == clusters[a])
    }, logical(1)))
    if (!length(eligible)) {
      .md_stop("No within-sample homotypic parent pairs are available. Supply a valid cluster_col or set n_homotypic = 0.")
    }
    hmat <- Matrix::Matrix(0, nrow = nrow(counts), ncol = n_homotypic, sparse = TRUE)
    hclust <- hpa <- hpb <- hsample <- character(n_homotypic)
    heff <- numeric(n_homotypic)
    for (i in seq_len(n_homotypic)) {
      a <- sample_values(eligible)
      possible <- which(seq_len(n) != a & samples == samples[a] & clusters == clusters[a])
      b <- sample_values(possible)
      eff <- stats::runif(1, doublet_capture_range[1], doublet_capture_range[2])
      hmat[, i] <- thin_column(counts[, a] + counts[, b], eff)
      hclust[i] <- as.character(clusters[a])
      hsample[i] <- as.character(samples[a])
      hpa[i] <- colnames(sce)[a]; hpb[i] <- colnames(sce)[b]; heff[i] <- eff
    }
    colnames(hmat) <- paste0("sim_homotypic_doublet_", seq_len(n_homotypic))
    sim_counts$homotypic_doublet <- hmat
    truth <- c(truth, rep("doublet", n_homotypic))
    truth_subclass <- c(truth_subclass, rep("homotypic_doublet", n_homotypic))
    sim_cluster <- c(sim_cluster, hclust)
    sim_sample <- c(sim_sample, hsample)
    parent_a <- c(parent_a, hpa); parent_b <- c(parent_b, hpb)
    simulation_scale <- c(simulation_scale, heff)
  }

  sets <- md_default_gene_sets(species)
  cycle_features <- unique(c(.md_match_genes(sce, sets$mitosis), .md_match_genes(sce, sets$cytokinesis)))

  if (n_mitotic > 0L) {
    idx <- sample_values(seq_len(n), n_mitotic, replace = TRUE)
    mmat <- counts[, idx, drop = FALSE]
    if (length(cycle_features)) {
      active <- sample_values(cycle_features, max(1L, ceiling(0.65 * length(cycle_features))))
      effects <- stats::rlnorm(length(active), meanlog = log(2.2), sdlog = 0.35)
      for (g in seq_along(active)) {
        mmat[active[g], ] <- round(mmat[active[g], , drop = FALSE] * effects[g] + stats::rpois(n_mitotic, 0.5))
      }
    }
    colnames(mmat) <- paste0("sim_mitotic_", seq_len(n_mitotic))
    sim_counts$mitotic <- mmat
    truth <- c(truth, rep("mitotic_singlet", n_mitotic))
    truth_subclass <- c(truth_subclass, rep("synthetic_mitotic_programme", n_mitotic))
    sim_cluster <- c(sim_cluster, as.character(clusters[idx]))
    sim_sample <- c(sim_sample, as.character(samples[idx]))
    parent_a <- c(parent_a, colnames(sce)[idx]); parent_b <- c(parent_b, rep(NA_character_, n_mitotic))
    simulation_scale <- c(simulation_scale, rep(NA_real_, n_mitotic))
  }

  if (n_polyploid_like > 0L) {
    idx <- sample_values(seq_len(n), n_polyploid_like, replace = TRUE)
    scales <- stats::runif(n_polyploid_like, 1.5, 2.3)
    pmat <- Matrix::Matrix(0, nrow = nrow(counts), ncol = n_polyploid_like, sparse = TRUE)
    for (i in seq_len(n_polyploid_like)) {
      v <- as.numeric(counts[, idx[i]])
      prob <- (v + 0.05) / sum(v + 0.05)
      target <- max(1L, round(sum(v) * scales[i]))
      pmat[, i] <- Matrix::Matrix(stats::rmultinom(1, target, prob), sparse = TRUE)
    }
    colnames(pmat) <- paste0("sim_polyploid_like_", seq_len(n_polyploid_like))
    sim_counts$polyploid <- pmat
    truth <- c(truth, rep("polyploid_like_singlet", n_polyploid_like))
    truth_subclass <- c(truth_subclass, rep("high_rna_same_identity_unresolved", n_polyploid_like))
    sim_cluster <- c(sim_cluster, as.character(clusters[idx]))
    sim_sample <- c(sim_sample, as.character(samples[idx]))
    parent_a <- c(parent_a, colnames(sce)[idx]); parent_b <- c(parent_b, rep(NA_character_, n_polyploid_like))
    simulation_scale <- c(simulation_scale, scales)
  }

  if (n_ambient > 0L) {
    idx <- sample_values(seq_len(n), n_ambient, replace = TRUE)
    ambient_profile <- Matrix::rowMeans(counts)
    amat <- counts[, idx, drop = FALSE]
    probs <- ambient_profile / sum(ambient_profile)
    for (i in seq_len(n_ambient)) {
      add_total <- max(1, round(sum(amat[, i]) * stats::runif(1, 0.03, 0.12)))
      add <- stats::rmultinom(1, size = add_total, prob = probs)
      amat[, i] <- amat[, i] + Matrix::Matrix(add, sparse = TRUE)
    }
    colnames(amat) <- paste0("sim_ambient_", seq_len(n_ambient))
    sim_counts$ambient <- amat
    truth <- c(truth, rep("stress_or_ambient", n_ambient))
    truth_subclass <- c(truth_subclass, rep("ambient_contamination", n_ambient))
    sim_cluster <- c(sim_cluster, as.character(clusters[idx]))
    sim_sample <- c(sim_sample, as.character(samples[idx]))
    parent_a <- c(parent_a, colnames(sce)[idx]); parent_b <- c(parent_b, rep(NA_character_, n_ambient))
    simulation_scale <- c(simulation_scale, rep(NA_real_, n_ambient))
  }

  if (n_stress > 0L) {
    idx <- sample_values(seq_len(n), n_stress, replace = TRUE)
    smat <- counts[, idx, drop = FALSE]
    stress_features <- .md_match_genes(sce, sets$stress)
    if (length(stress_features)) {
      active <- sample_values(stress_features, max(1L, ceiling(0.60 * length(stress_features))))
      smat[active, ] <- round(smat[active, , drop = FALSE] * stats::runif(1, 1.8, 3.5) + 1)
    }
    colnames(smat) <- paste0("sim_stress_", seq_len(n_stress))
    sim_counts$stress <- smat
    truth <- c(truth, rep("stress_or_ambient", n_stress))
    truth_subclass <- c(truth_subclass, rep("stress_programme", n_stress))
    sim_cluster <- c(sim_cluster, as.character(clusters[idx]))
    sim_sample <- c(sim_sample, as.character(samples[idx]))
    parent_a <- c(parent_a, colnames(sce)[idx]); parent_b <- c(parent_b, rep(NA_character_, n_stress))
    simulation_scale <- c(simulation_scale, rep(NA_real_, n_stress))
  }

  out_counts <- do.call(cbind, sim_counts)
  out <- SingleCellExperiment(assays = list(counts = out_counts))
  colData(out)$truth_class <- truth
  colData(out)$truth_subclass <- truth_subclass
  colData(out)$sim_cluster <- sim_cluster
  colData(out)$sim_sample <- sim_sample
  if (!is.null(sample_col)) colData(out)[[sample_col]] <- sim_sample
  colData(out)$sim_parent_a <- parent_a
  colData(out)$sim_parent_b <- parent_b
  colData(out)$sim_scale <- simulation_scale
  out <- .md_add_logcounts(out)
  out
  })
}
