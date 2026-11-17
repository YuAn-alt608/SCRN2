###############################################################################
# Step 21: colocalization prior sensitivity -- p12 = 1e-4 / 1e-5 / 1e-6
#
# Aim: PP.H4 is a Bayesian posterior and depends strongly on p12, the prior
#      probability that the two traits share a causal variant. The main analysis
#      uses the coloc default 1e-5; this step tests whether the 71 pairs with
#      PP.H4>0.8 are merely a product of that prior.
#
# Design: data preparation runs once (merge, allele alignment, MAF filter, beta
#         recovery); the same d1/d2 objects are fed to all three p12 values, so
#         any difference comes from the prior alone.
#         The p12=1e-5 output is also reconciled against results/coloc_results.csv
#         (the two must agree exactly).
#
# Input : data/eqtl_regions_all.txt.gz, data/region_*.tsv.gz, data/coloc_genes_all.tsv
# Output: results/coloc_p12_long.csv (long format: gene, outcome, p12, PP.H3, PP.H4, top_snp, top_snp_PP)
###############################################################################

suppressMessages({
  library(data.table)
  library(coloc)
})

BASE <- "D:/NewBeginning/MR1/DATA10"
setwd(BASE)

P12 <- c(1e-4, 1e-5, 1e-6)

OUTCOMES <- list(
  TSH  = list(file = "data/region_TSH.tsv.gz",  type = "quant", N = 482873,
              label = "TSH"),
  T4   = list(file = "data/region_T4.tsv.gz",   type = "quant", N = 191449,
              label = "Free T4"),
  HYPO = list(file = "data/region_HYPO.tsv.gz", type = "cc",    N = 1046991,
              s = 98805 / 1046991, label = "Hypothyroidism"),
  HASH = list(file = "data/region_HASH.tsv.gz", type = "cc",    N = 395640,
              s = 15654 / 395640, label = "Hashimoto thyroiditis")
)

eq <- fread("data/eqtl_regions_all.txt.gz", showProgress = FALSE)
setnames(eq,
         c("pval_eqtl", "snp", "snp_chr", "snp_pos", "z", "ea", "oa",
           "gene_id", "gene", "gene_chr", "gene_pos", "n_cohort", "n_eqtl", "fdr"))
eq[, snp_pos := as.integer(snp_pos)]
genes <- sort(unique(fread("data/coloc_genes_all.tsv")$gene))
eq <- eq[gene %in% genes]
cat(sprintf("eQTL records: %s rows, %d genes\n", format(nrow(eq), big.mark = ","),
            uniqueN(eq$gene)))

res <- list()

for (oc in names(OUTCOMES)) {
  cfg <- OUTCOMES[[oc]]
  o <- fread(cfg$file, showProgress = FALSE)
  setnames(o, c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                "beta", "se", "eaf", "pval", "rsid", "n"), skip_absent = TRUE)
  o[, base_pair_location := as.integer(base_pair_location)]
  setkey(o, chromosome, base_pair_location)
  cat(sprintf("\n=== %s: %s regional variants ===\n", cfg$label, format(nrow(o), big.mark = ",")))

  for (g in genes) {
    e <- eq[gene == g]
    if (nrow(e) < 50) next

    m <- merge(e, o, by.x = c("snp_chr", "snp_pos"),
               by.y = c("chromosome", "base_pair_location"), suffixes = c("", ".o"))
    if (nrow(m) < 50) next

    m[, flip := (ea == other_allele & oa == effect_allele)]
    m <- m[(ea == effect_allele & oa == other_allele) |
           (ea == other_allele & oa == effect_allele)]
    if (nrow(m) < 50) next
    m[, z_aligned := fifelse(flip, -z, z)]
    setorder(m, pval_eqtl)
    m <- m[!duplicated(snp)]
    if (nrow(m) < 50) next

    maf <- pmin(m$eaf, 1 - m$eaf)
    keep2 <- !is.na(maf) & maf > 0.01 & maf < 0.5
    m <- m[keep2]
    maf <- maf[keep2]
    if (nrow(m) < 50) next

    info <- 2 * maf * (1 - maf) * (m$n_eqtl + m$z^2)
    d1 <- list(snp = m$snp, beta = m$z_aligned / sqrt(info), varbeta = 1 / info,
               type = "quant", N = as.integer(median(m$n_eqtl)), MAF = maf, sdY = 1)
    d2 <- if (cfg$type == "cc") {
      list(snp = m$snp, beta = m$beta, varbeta = m$se^2, type = "cc",
           N = cfg$N, s = cfg$s, MAF = maf)
    } else {
      list(snp = m$snp, beta = m$beta, varbeta = m$se^2, type = "quant",
           N = cfg$N, MAF = maf, sdY = 1)
    }

    for (p in P12) {
      cr <- tryCatch(
        suppressMessages(suppressWarnings(coloc.abf(d1, d2, p12 = p))),
        error = function(err) NULL)
      if (is.null(cr)) next
      s <- cr$summary
      top <- cr$results[which.max(cr$results$SNP.PP.H4), ]
      res[[length(res) + 1]] <- data.table(
        gene = g, outcome = oc, p12 = p, nsnps = as.integer(s["nsnps"]),
        PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]],
        top_snp = top$snp, top_snp_PP = top$SNP.PP.H4)
    }
  }
  cat(sprintf("  %s done, %d rows accumulated\n", oc, length(res)))
  flush.console()
}

out <- rbindlist(res)
fwrite(out, "results/coloc_p12_long.csv")
cat(sprintf("\n-> results/coloc_p12_long.csv  %d rows\n", nrow(out)))

cat("\n===== pairs with PP.H4>0.8 under each prior =====\n")
print(out[, .(n_pairs = .N, H4_gt_0.8 = sum(PP.H4 > 0.8),
              H4_gt_0.5 = sum(PP.H4 > 0.5)), by = .(p12)][order(-p12)])

cat("\n===== by outcome =====\n")
print(dcast(out[PP.H4 > 0.8, .N, by = .(outcome, p12)], outcome ~ p12, value.var = "N"))

# reconcile against the main analysis
main <- fread("results/coloc_results.csv")
chk <- merge(out[p12 == 1e-5], main[, .(gene, outcome, PP.H4_main = PP.H4)],
             by = c("gene", "outcome"))
chk[, d := abs(PP.H4 - PP.H4_main)]
cat(sprintf("\nreconciliation: %d pairs, max absolute PP.H4 difference %.3g (should be 0)\n", nrow(chk), max(chk$d)))
