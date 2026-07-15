test_that("normalized matrices are rejected as counts", {
  x <- matrix(stats::runif(1000), nrow = 100)
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  expect_error(md_as_sce(x), "appears normalized")
})

test_that("cell-cycle scoring records gene coverage", {
  sim <- md_simulate_typical_dataset(
    n_per_type = c(cardiomyocyte = 12L, fibroblast = 12L),
    n_doublets = 0L, n_homotypic = 0L, n_background = 250L, seed = 3L
  )
  out <- md_score_cell_cycle(sim, species = "human")
  expect_true(all(c("md_mitosis_score", "md_cytokinesis_score") %in% colnames(SummarizedExperiment::colData(out))))
  coverage <- S4Vectors::metadata(out)$md_gene_set_coverage
  expect_true(coverage$md_mitosis_score[["matched"]] >= 3L)
})
