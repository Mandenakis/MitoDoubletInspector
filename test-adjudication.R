# Internal utility functions ---------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

.md_msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) message(...)
}

.md_stop <- function(...) stop(..., call. = FALSE)

.md_require <- function(pkg, purpose = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    extra <- if (!is.null(purpose)) paste0(" for ", purpose) else ""
    .md_stop("Package '", pkg, "' is required", extra, ". Install it first.")
  }
  invisible(TRUE)
}

.md_warn_once <- local({
  seen <- new.env(parent = emptyenv())
  function(key, ...) {
    if (!exists(key, envir = seen, inherits = FALSE)) {
      assign(key, TRUE, envir = seen)
      warning(..., call. = FALSE)
    }
  }
})

.md_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.md_key <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  substr(x, 1, 60)
}

.md_is_seurat <- function(x) inherits(x, "Seurat")
.md_is_sce <- function(x) inherits(x, "SingleCellExperiment")
.md_is_matrix <- function(x) inherits(x, c("matrix", "dgCMatrix", "Matrix"))

.md_assay_names <- function(sce) assayNames(sce)

.md_get_assay <- function(sce, assay_name = NULL) {
  an <- .md_assay_names(sce)
  if (length(an) == 0L) .md_stop("Input object contains no assays.")
  if (!is.null(assay_name) && assay_name %in% an) return(assay(sce, assay_name))
  if ("counts" %in% an) return(assay(sce, "counts"))
  assay(sce, an[1])
}

.md_get_counts <- function(sce) {
  an <- .md_assay_names(sce)
  if ("counts" %in% an) return(assay(sce, "counts"))
  .md_warn_once("no_counts_assay", "No assay named 'counts' found; using the first assay as count-like input.")
  assay(sce, an[1])
}

.md_sparse <- function(x) {
  if (inherits(x, "sparseMatrix")) return(x)
  Matrix::Matrix(x, sparse = TRUE)
}

.md_add_logcounts <- function(sce, scale_factor = 1e4) {
  if ("logcounts" %in% .md_assay_names(sce)) return(sce)
  counts <- .md_sparse(.md_get_counts(sce))
  lib <- Matrix::colSums(counts)
  lib[is.na(lib) | lib <= 0] <- 1
  sf <- scale_factor / lib
  norm <- counts %*% Matrix::Diagonal(x = sf)
  norm@x <- log1p(norm@x)
  assay(sce, "logcounts") <- norm
  sce
}

.md_gene_symbols <- function(sce) {
  rd <- rowData(sce)
  candidate_cols <- c("symbol", "gene_symbol", "gene", "gene_name", "external_gene_name")
  for (cc in candidate_cols) {
    if (cc %in% colnames(rd)) {
      z <- as.character(rd[[cc]])
      z[is.na(z) | z == ""] <- rownames(sce)[is.na(z) | z == ""]
      return(z)
    }
  }
  rownames(sce)
}

.md_match_genes <- function(sce, genes) {
  if (length(genes) == 0L) return(character())
  sym <- .md_gene_symbols(sce)
  idx <- match(toupper(unique(genes)), toupper(sym))
  idx <- idx[!is.na(idx)]
  unique(rownames(sce)[idx])
}

.md_gene_pattern <- function(sce, pattern) {
  sym <- .md_gene_symbols(sce)
  rownames(sce)[grepl(pattern, sym, ignore.case = FALSE)]
}

.md_z <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(0, length(x)))
  med <- stats::median(x, na.rm = TRUE)
  x[is.na(x)] <- med
  s <- stats::sd(x)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

.md_robust_z <- function(x) {
  x <- as.numeric(x)
  finite <- is.finite(x)
  if (!any(finite)) return(rep(0, length(x)))
  med <- stats::median(x[finite])
  x[!finite] <- med
  s <- stats::mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s <= .Machine$double.eps) return(rep(0, length(x)))
  (x - med) / s
}

.md_percentile <- function(x, group = NULL, min_group_n = 20L) {
  x <- as.numeric(x)
  global <- function(z) {
    ok <- is.finite(z)
    out <- rep(0.5, length(z))
    if (sum(ok) > 1L && length(unique(z[ok])) > 1L) {
      out[ok] <- (rank(z[ok], ties.method = "average") - 0.5) / sum(ok)
    }
    out
  }
  if (is.null(group)) return(global(x))
  g <- as.character(group)
  out <- rep(NA_real_, length(x))
  for (lv in unique(g[!is.na(g)])) {
    idx <- which(g == lv)
    if (length(idx) >= min_group_n) out[idx] <- global(x[idx])
  }
  miss <- !is.finite(out)
  if (any(miss)) out[miss] <- global(x)[miss]
  out
}

.md_group_robust_z <- function(x, group = NULL, min_group_n = 20L) {
  x <- as.numeric(x)
  if (is.null(group)) return(.md_robust_z(x))
  g <- as.character(group)
  out <- rep(NA_real_, length(x))
  for (lv in unique(g[!is.na(g)])) {
    idx <- which(g == lv)
    if (length(idx) >= min_group_n) out[idx] <- .md_robust_z(x[idx])
  }
  miss <- !is.finite(out)
  if (any(miss)) out[miss] <- .md_robust_z(x)[miss]
  out
}

.md_minmax <- function(x) {
  x <- as.numeric(x)
  r <- range(x, finite = TRUE, na.rm = TRUE)
  if (any(!is.finite(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

.md_clip <- function(x, lower, upper) pmin(pmax(x, lower), upper)

.md_coldata_df <- function(sce) as.data.frame(colData(sce))

.md_col_exists <- function(sce, col) !is.null(col) && col %in% colnames(colData(sce))

.md_infer_col <- function(sce, candidates) {
  cd <- colnames(colData(sce))
  hit <- candidates[candidates %in% cd]
  if (length(hit)) hit[1] else NULL
}

.md_safe_quantile <- function(x, p, default = NA_real_) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 7))
}

.md_has_variation <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  length(x) > 1L && length(unique(x)) > 1L
}

.md_as_proportions <- function(x) {
  x <- .md_sparse(x)
  lib <- Matrix::colSums(x)
  lib[!is.finite(lib) | lib <= 0] <- 1
  x %*% Matrix::Diagonal(x = 1 / lib)
}

.md_softmax <- function(mat) {
  mat <- as.matrix(mat)
  out <- matrix(NA_real_, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  for (i in seq_len(nrow(mat))) {
    z <- mat[i, ]
    z[!is.finite(z)] <- -Inf
    m <- max(z)
    if (!is.finite(m)) {
      out[i, ] <- rep(1 / ncol(mat), ncol(mat))
    } else {
      ez <- exp(z - m)
      out[i, ] <- ez / sum(ez)
    }
  }
  out
}

.md_safe_colmeans <- function(x) {
  if (inherits(x, "sparseMatrix")) Matrix::colMeans(x) else colMeans(x)
}

.md_safe_rowmeans <- function(x) {
  if (inherits(x, "sparseMatrix")) Matrix::rowMeans(x) else rowMeans(x)
}

.md_row_variance_sparse <- function(x) {
  x <- .md_sparse(x)
  n <- ncol(x)
  if (n <= 1L) return(rep(0, nrow(x)))
  mu <- Matrix::rowMeans(x)
  ss <- Matrix::rowSums(x * x)
  v <- (ss - n * mu^2) / (n - 1)
  v[!is.finite(v) | v < 0] <- 0
  as.numeric(v)
}

.md_col_norms <- function(x) {
  x <- .md_sparse(x)
  sqrt(Matrix::colSums(x * x))
}

.md_profile_norm <- function(profile) {
  profile <- as.matrix(profile)
  nrm <- sqrt(colSums(profile^2))
  nrm[!is.finite(nrm) | nrm == 0] <- 1
  sweep(profile, 2, nrm, "/")
}

.md_subset_features <- function(sce, features, min_features = 10L) {
  features <- intersect(features, rownames(sce))
  if (length(features) < min_features) {
    .md_stop("Too few usable features after filtering: ", length(features), ".")
  }
  features
}

.md_return_object <- function(original, sce) {
  if (.md_is_seurat(original)) {
    .md_require("SeuratObject", "returning metadata to a Seurat object")
    md_cols <- grep("^md_", colnames(colData(sce)), value = TRUE)
    meta <- as.data.frame(colData(sce))[, md_cols, drop = FALSE]
    meta <- meta[colnames(original), , drop = FALSE]
    return(SeuratObject::AddMetaData(original, metadata = meta))
  }
  sce
}
