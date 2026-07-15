# mitoDoubletR

Mitosis-aware doublet adjudication for scRNA-seq (`miDAS`).

`mitoDoubletR` is a second-stage adjudicator. It does not replace `scDblFinder` or another primary caller. It asks whether a suspicious barcode is better explained by two independent transcriptomic identities, by one coherent identity carrying a mitotic/cytokinetic programme, by a coherent high-RNA state, or by stress/ambient contamination.

The final classes are:

- `singlet`
- `doublet`
- `mitotic_singlet`
- `polyploid_like_singlet`
- `stress_or_ambient`
- `ambiguous`

`polyploid_like_singlet` is deliberately not called “polyploid.” A same-identity homotypic doublet and a high-RNA/polyploid singlet can be observationally equivalent in RNA counts. Such cases carry `md_identifiability = "high_rna_singlet_or_homotypic_doublet_unresolved"`.

## What changed in 0.3.0

Version 0.2.0 contained valuable design improvements but the GitHub upload had flattened the package and permuted file contents. Version 0.3.0 reconstructs an installable package and changes five statistical decisions:

1. Full and no-cycle `scDblFinder` outputs are compared by within-sample percentile, because independently trained scores are not guaranteed to share an absolute scale.
2. Identity mixtures use mean per-cell count proportions, not exponentiated mean log-expression.
3. Mixture improvement is calibrated against an empirical singlet null with p/q values and a 95th-percentile search correction.
4. Class values are evidence supports, not probabilities. A doublet requires at least two of external-caller, lineage-conflict, and parent-mixture evidence.
5. Rare clusters are routed to uncertainty when evidence conflicts; rarity alone is never evidence of singlet status.

## Install

```r
install.packages("path/to/mitoDoubletR", repos = NULL, type = "source")
```

## Typical use

```r
library(mitoDoubletR)

obj <- md_adjudicate(
  obj,
  sample_col = "sample",
  cluster_col = "seurat_clusters",
  species = "mouse",
  run_scdblfinder = TRUE,
  protect_rare = TRUE,
  seed = 1234
)

md_summary(obj)
rescued <- md_rescue_candidates(obj)
audit <- md_reason_table(obj)
```

To use an existing caller instead:

```r
obj <- md_adjudicate(
  obj,
  doublet_score_col = "scDblFinder.score",
  doublet_class_col = "scDblFinder.class",
  run_scdblfinder = FALSE
)
```

An imported score alone cannot create the no-cycle counterfactual; `md_counterfactual_available` will therefore be false. Identity-mixture and coherence evidence still run.

## Reproducible stress test

```r
sim <- md_simulate_typical_dataset(seed = 11)
bench <- md_benchmark_adjudication(
  sim,
  cluster_col = "sim_cluster",
  truth_col = "truth_class",
  sample_col = "sample",
  run_scdblfinder = TRUE,
  seed = 11
)

bench$metrics
bench$subtype_confusion
```

The synthetic benchmark is a development test, not proof of biological validity. The external-validation vignette specifies orthogonal doublet, FUCCI cell-cycle, and DNA-ploidy datasets.

## Safe interpretation

- A `mitotic_singlet` call means a coherent identity plus an RNA mitotic programme and weak two-parent evidence. It does not prove that cytokinesis completed.
- A `polyploid_like_singlet` call means a coherent high-RNA profile. It does not establish DNA content and cannot exclude a homotypic doublet.
- Rescue candidates should be retained for downstream sensitivity analysis, not silently relabelled as ground-truth singlets.
- Benchmark sensitivity separately for heterotypic and homotypic doublets; the latter are intrinsically difficult from RNA alone.
