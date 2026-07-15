# Independent synthetic reference dataset -------------------------------------

#' Simulate a typical multi-lineage scRNA-seq count matrix
#'
#' Generates a compact Gamma-Poisson UMI dataset with cardiomyocyte,
#' fibroblast, endothelial, myeloid, and rare neural-crest-like identities;
#' mitotic programmes; coherent high-RNA states; stress and ambient-contaminated
#' singlets; heterotypic and homotypic doublets; two batches; and complete
#' ground truth.
#' The generator uses latent mean programmes and count sampling rather than
#' editing normalized expression. It is suitable for smoke tests and failure
#' analysis, but it is not a substitute for experimentally labelled data.
#'
#' @param n_per_type Named integer vector of singlets per identity.
#' @param n_doublets,n_homotypic Numbers of heterotypic and homotypic doublets.
#' @param mitotic_fraction,high_rna_fraction,stress_fraction,ambient_fraction
#'   Fractions of base singlets assigned to these states. States are mutually
#'   exclusive.
#' @param n_background Number of non-marker background genes.
#' @param species Human or mouse symbol casing.
#' @param seed Random seed, restored on exit.
#' @return SingleCellExperiment with counts and truth metadata.
#' @export
md_simulate_typical_dataset <- function(
  n_per_type = c(
    cardiomyocyte = 300L, fibroblast = 250L, endothelial = 200L,
    myeloid = 150L, neural_crest_like = 50L
  ),
  n_doublets = 120L,
  n_homotypic = 60L,
  mitotic_fraction = 0.08,
  high_rna_fraction = 0.05,
  stress_fraction = 0.05,
  ambient_fraction = 0.05,
  n_background = 1200L,
  species = c("human", "mouse"),
  seed = 1L
) {
  species <- match.arg(species)
  if (!length(n_per_type) || is.null(names(n_per_type)) ||
      any(!is.finite(n_per_type)) || any(n_per_type < 5L)) {
    .md_stop("n_per_type must be a named vector with at least five cells per identity.")
  }
  doublet_n <- c(n_doublets, n_homotypic)
  if (length(doublet_n) != 2L || any(!is.finite(doublet_n)) ||
      any(abs(n_per_type - round(n_per_type)) > 0) ||
      any(doublet_n < 0) ||
      any(abs(c(n_doublets, n_homotypic) - round(c(n_doublets, n_homotypic))) > 0)) {
    .md_stop("Cell and doublet counts must be non-negative finite integers.")
  }
  if (n_doublets > 0L && length(n_per_type) < 2L) .md_stop("At least two identities are required for heterotypic doublets.")
  state_fractions <- c(mitotic_fraction, high_rna_fraction, stress_fraction, ambient_fraction)
  if (length(state_fractions) != 4L || any(!is.finite(state_fractions)) ||
      any(state_fractions < 0) || sum(state_fractions) >= 0.8) {
    .md_stop("State fractions must be finite, non-negative, and sum to less than 0.8.")
  }

  .md_with_seed(seed, {
    sample_values <- function(x, size) {
      if (size <= 0L) return(x[integer()])
      x[sample.int(length(x), size = size)]
    }
    markers <- list(
      cardiomyocyte = c("TNNT2", "ACTN2", "MYH6", "MYH7", "NKX2-5", "PLN", "RYR2", "TTN", "TNNI3", "DES"),
      fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "VIM", "COL6A1"),
      endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2", "ESAM", "CDH5"),
      myeloid = c("PTPRC", "LST1", "TYROBP", "FCER1G", "CTSS", "LYZ", "CSF1R", "AIF1"),
      neural_crest_like = c("SOX10", "S100B", "PLP1", "MPZ", "FOXD3", "TFAP2B", "ERBB3", "NGFR")
    )
    sets <- md_default_gene_sets("human")
    technical <- c("MALAT1", "RPLP0", "RPS18", "MT-CO1", "MT-ND1", "HBB", "HBA1", "ALB")
    genes <- unique(c(
      unlist(markers, use.names = FALSE), sets$S, sets$G2M, sets$mitosis,
      sets$cytokinesis, sets$stress, sets$ambient, technical,
      paste0("BG", sprintf("%04d", seq_len(as.integer(n_background))))
    ))
    if (species == "mouse") {
      genes <- vapply(genes, function(g) {
        if (grepl("^MT-", g)) sub("^MT-", "mt-", g) else paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))))
      }, character(1))
      markers <- lapply(markers, function(gs) {
        vapply(gs, function(g) paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g)))), character(1))
      })
      sets <- md_default_gene_sets("mouse")
    }

    cell_type <- rep(names(n_per_type), times = as.integer(n_per_type))
    n <- length(cell_type)
    cell_names <- paste0("cell_", seq_len(n))
    batch <- rep(c("batch_1", "batch_2"), length.out = n)
    batch <- sample(batch)

    state <- rep("singlet", n)
    state_n <- function(fraction) if (fraction > 0) max(1L, round(n * fraction)) else 0L
    n_mit <- state_n(mitotic_fraction)
    mit_idx <- sample_values(seq_len(n), n_mit)
    remaining <- setdiff(seq_len(n), mit_idx)
    n_high <- state_n(high_rna_fraction)
    high_idx <- sample_values(remaining, min(n_high, length(remaining)))
    remaining <- setdiff(remaining, high_idx)
    n_stress <- state_n(stress_fraction)
    stress_idx <- sample_values(remaining, min(n_stress, length(remaining)))
    remaining <- setdiff(remaining, stress_idx)
    n_ambient <- state_n(ambient_fraction)
    ambient_idx <- sample_values(remaining, min(n_ambient, length(remaining)))
    state[mit_idx] <- "mitotic_singlet"
    state[high_idx] <- "polyploid_like_singlet"
    state[stress_idx] <- "stress_or_ambient"
    state[ambient_idx] <- "stress_or_ambient"

    base_mean <- stats::rgamma(length(genes), shape = 0.7, rate = 2.5) + 0.01
    names(base_mean) <- genes
    means <- matrix(base_mean, nrow = length(genes), ncol = n,
                    dimnames = list(genes, cell_names))
    for (ct in names(markers)) {
      idx <- which(cell_type == ct)
      mg <- intersect(markers[[ct]], genes)
      means[mg, idx] <- means[mg, idx, drop = FALSE] * 9 + 0.6
    }

    mit_genes <- intersect(unique(c(sets$mitosis, sets$cytokinesis)), genes)
    if (length(mit_genes)) {
      active <- sample(mit_genes, max(5L, ceiling(0.70 * length(mit_genes))))
      fc <- stats::rlnorm(length(active), log(3), 0.25)
      means[active, mit_idx] <- means[active, mit_idx, drop = FALSE] * fc + 0.5
    }
    stress_genes <- intersect(sets$stress, genes)
    if (length(stress_genes) && length(stress_idx)) {
      active <- sample(stress_genes, max(3L, ceiling(0.70 * length(stress_genes))))
      means[active, stress_idx] <- means[active, stress_idx, drop = FALSE] * 4 + 0.4
    }
    lib_factor <- stats::rlnorm(n, meanlog = 0, sdlog = 0.35)
    lib_factor[high_idx] <- lib_factor[high_idx] * stats::runif(length(high_idx), 1.6, 2.2)
    batch_genes <- sample(genes, min(100L, length(genes)))
    means[batch_genes, batch == "batch_2"] <- means[batch_genes, batch == "batch_2", drop = FALSE] * 1.25
    means <- sweep(means, 2, lib_factor, "*")

    counts <- matrix(
      stats::rnbinom(length(means), mu = as.numeric(means), size = 1.5),
      nrow = nrow(means), ncol = ncol(means), dimnames = dimnames(means)
    )
    counts <- Matrix::Matrix(counts, sparse = TRUE)

    # Add a realistic low-fraction soup profile after count sampling. These
    # remain biological singlets; the truth label records the contamination.
    if (length(ambient_idx)) {
      soup <- Matrix::rowSums(counts) + 1
      soup <- as.numeric(soup / sum(soup))
      for (j in ambient_idx) {
        extra_n <- max(1L, round(sum(counts[, j]) * stats::runif(1, 0.05, 0.15)))
        counts[, j] <- counts[, j] + Matrix::Matrix(stats::rmultinom(1, extra_n, soup), sparse = TRUE)
      }
    }

    make_doublets <- function(number, homotypic) {
      out <- Matrix::Matrix(0, nrow = nrow(counts), ncol = number, sparse = TRUE)
      pa <- pb <- character(number)
      pcl <- pbatch <- character(number)
      eligible <- which(vapply(seq_len(n), function(a) {
        same_batch <- batch == batch[a]
        pool <- if (homotypic) {
          which(cell_type == cell_type[a] & same_batch)
        } else {
          which(cell_type != cell_type[a] & same_batch)
        }
        length(setdiff(pool, a)) > 0L
      }, logical(1)))
      if (number > 0L && !length(eligible)) {
        .md_stop("No valid within-batch parent pairs are available for the requested doublets.")
      }
      for (i in seq_len(number)) {
        a <- sample_values(eligible, 1)
        same_batch <- batch == batch[a]
        pool <- if (homotypic) {
          which(cell_type == cell_type[a] & same_batch)
        } else {
          which(cell_type != cell_type[a] & same_batch)
        }
        pool <- setdiff(pool, a)
        b <- sample_values(pool, 1)
        v <- .md_sparse(counts[, a] + counts[, b])
        if (length(v@x)) v@x <- stats::rbinom(length(v@x), round(v@x), stats::runif(1, 0.55, 0.9))
        out[, i] <- v
        pa[i] <- cell_names[a]; pb[i] <- cell_names[b]; pcl[i] <- cell_type[a]; pbatch[i] <- batch[a]
      }
      list(counts = out, parent_a = pa, parent_b = pb, cluster = pcl, batch = pbatch)
    }

    hetero <- make_doublets(as.integer(n_doublets), FALSE)
    homo <- make_doublets(as.integer(n_homotypic), TRUE)
    colnames(hetero$counts) <- paste0("heterotypic_doublet_", seq_len(n_doublets))
    colnames(homo$counts) <- paste0("homotypic_doublet_", seq_len(n_homotypic))
    all_counts <- cbind(counts, hetero$counts, homo$counts)

    truth_class <- c(state, rep("doublet", n_doublets + n_homotypic))
    singlet_subclass <- ifelse(state == "singlet", "ordinary_singlet", state)
    singlet_subclass[seq_len(n) %in% stress_idx] <- "stress_programme"
    singlet_subclass[seq_len(n) %in% ambient_idx] <- "ambient_contamination"
    truth_subclass <- c(
      singlet_subclass,
      rep("heterotypic_doublet", n_doublets), rep("homotypic_doublet", n_homotypic)
    )
    sim_cluster <- c(cell_type, hetero$cluster, homo$cluster)
    parent_a <- c(rep(NA_character_, n), hetero$parent_a, homo$parent_a)
    parent_b <- c(rep(NA_character_, n), hetero$parent_b, homo$parent_b)
    out_batch <- c(batch, hetero$batch, homo$batch)

    out <- SingleCellExperiment(assays = list(counts = all_counts))
    colData(out)$truth_class <- truth_class
    colData(out)$truth_subclass <- truth_subclass
    colData(out)$sim_cluster <- sim_cluster
    colData(out)$sample <- out_batch
    colData(out)$sim_parent_a <- parent_a
    colData(out)$sim_parent_b <- parent_b
    rowData(out)$symbol <- rownames(out)
    .md_add_logcounts(out)
  })
}
