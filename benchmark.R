# Default gene sets ------------------------------------------------------------

#' Default gene sets used by mitoDoubletR
#'
#' Returns conservative RNA-only signatures for S phase, G2/M, mitosis,
#' cytokinesis, stress/dissociation, and a deliberately narrow set of
#' ambient-prone transcripts. Mitosis is not defined as the whole G2/M list:
#' this avoids calling every G2 cell a mitotic cell.
#'
#' @param species "human" or "mouse". Currently affects only returned casing;
#'   matching is case-insensitive internally.
#' @return Named list of character vectors.
#' @export
md_default_gene_sets <- function(species = c("human", "mouse")) {
  species <- match.arg(species)

  s_genes <- c(
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2",
    "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
    "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "BRIP1", "E2F8"
  )

  g2m_genes <- c(
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
    "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
    "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
  )

  mitosis <- c(
    "CCNB1", "CDK1", "PLK1", "AURKA", "AURKB", "BUB1", "BUB1B", "MAD2L1",
    "TTK", "CDC20", "PTTG1", "ESPL1", "NDC80", "NUF2", "SPC24", "SPC25",
    "CASC5", "KNL1", "ZWINT", "CENPE", "CENPF", "INCENP", "BIRC5", "CDCA8",
    "KIF2C", "KIF4A", "KIF11", "KIF18A", "TPX2", "ASPM", "NCAPG", "NCAPH"
  )

  cytokinesis <- c(
    "AURKB", "ANLN", "ECT2", "KIF23", "KIF20A", "KIF20B", "PRC1", "RACGAP1",
    "CIT", "PLK1", "CDK1", "CKAP2", "CKAP2L", "DLGAP5", "TACC3", "CEP55",
    "NUSAP1", "KIF4A", "KIF11", "KIF2C", "KIFC1", "INCENP", "BIRC5", "CDCA8"
  )

  stress <- c(
    "FOS", "JUN", "JUNB", "JUND", "ATF3", "DUSP1", "DUSP2", "EGR1", "EGR2",
    "BTG1", "BTG2", "IER2", "IER3", "NR4A1", "NR4A2", "HSPA1A", "HSPA1B",
    "HSP90AA1", "HSP90AB1", "DNAJA1", "DNAJB1", "HSPB1", "HSPH1", "HSPA8",
    "SOCS3", "KLF2", "KLF4", "ZFP36", "PPP1R15A", "GADD45B"
  )

  ambient <- c(
    "HBB", "HBA1", "HBA2", "HBD", "HBE1", "HBM", "ALB", "APOA1", "APOA2",
    "PRSS1", "PRSS2", "REG1A", "CPA1"
  )

  sets <- list(
    S = unique(s_genes),
    G2M = unique(g2m_genes),
    mitosis = unique(mitosis),
    cytokinesis = unique(cytokinesis),
    stress = unique(stress),
    ambient = unique(ambient)
  )

  if (identical(species, "mouse")) {
    sets <- lapply(sets, function(x) {
      vapply(x, function(g) {
        if (grepl("^MT-", g)) sub("^MT-", "mt-", g) else paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))))
      }, character(1))
    })
    sets$ambient <- unique(c(sets$ambient, "Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alb", "Apoa1", "Apoa2"))
  }
  sets
}
