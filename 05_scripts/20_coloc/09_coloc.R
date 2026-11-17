###############################################################################
# Step 09: colocalization (coloc.abf) -- all significant genes
#
# Aim: test whether the eQTL and the thyroid outcome share the same causal variant
#      -- a question plain MR cannot answer but every reviewer asks
#
# Input : data/eqtl_regions_all.txt.gz           eQTLGen cis-eQTL (616 genes, +/-500 kb)
#         data/region_{TSH,T4,HYPO,HASH}.tsv.gz outcome statistics within each region
# Output: results/coloc_results.csv
#
# Note: eQTLGen releases Z scores only, not beta/MAF. Standard conversion:
#       beta = z / sqrt(2*maf*(1-maf)*(N+z^2)), varbeta = 1/(...)
#       MAF is taken from the EAF of the same-population outcome GWAS (after allele
#       alignment); the residual difference is sampling error only.
###############################################################################

suppressMessages({
  library(data.table)
  library(coloc)
})

BASE <- "D:/NewBeginning/MR1/DATA10"
setwd(BASE)
dir.create("results", showWarnings = FALSE)

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
gene_list <- fread("data/coloc_genes_all.tsv")$gene
eq <- eq[gene %in% gene_list]
cat(sprintf("eQTL records: %s rows, %d genes\n", format(nrow(eq), big.mark = ","),
            uniqueN(eq$gene)))

genes <- sort(unique(eq$gene))
res <- list()

ARGS <- commandArgs(trailingOnly = TRUE)
if (length(ARGS)) {
  OUTCOMES <- OUTCOMES[intersect(names(OUTCOMES), ARGS)]
  cat("Running outcomes only:", paste(names(OUTCOMES), collapse = ", "), "\n")
}

for (oc in names(OUTCOMES)) {
  cfg <- OUTCOMES[[oc]]
  if (!file.exists(cfg$file)) { cat("skip (file missing):", oc, "\n"); next }
  o <- fread(cfg$file, showProgress = FALSE)
  setnames(o, c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                "beta", "se", "eaf", "pval", "rsid", "n"),
           skip_absent = TRUE)
  o[, base_pair_location := as.integer(base_pair_location)]
  setkey(o, chromosome, base_pair_location)
  cat(sprintf("\n=== %s: %s regional variants ===\n", cfg$label, format(nrow(o), big.mark = ",")))

  for (g in genes) {
    e <- eq[gene == g]
    if (nrow(e) < 50) next

    m <- merge(e, o, by.x = c("snp_chr", "snp_pos"),
               by.y = c("chromosome", "base_pair_location"), suffixes = c("", ".o"))
    if (nrow(m) < 50) next

    # ---- allele alignment ----
    # The outcome file may hold several variants at the same physical position
    # (multi-allelic sites). Align alleles first, then keep a single variant per
    # position by eQTL significance, so coloc does not report duplicated snps.
    m[, flip := (ea == other_allele & oa == effect_allele)]
    m <- m[(ea == effect_allele & oa == other_allele) |
           (ea == other_allele & oa == effect_allele)]
    if (nrow(m) < 50) next
    m[, z_aligned := fifelse(flip, -z, z)]
    setorder(m, pval_eqtl)
    m <- m[!duplicated(snp)]
    if (nrow(m) < 50) next

    # ---- recover MAF and beta ----
    maf <- pmin(m$eaf, 1 - m$eaf)
    keep2 <- !is.na(maf) & maf > 0.01 & maf < 0.5
    m <- m[keep2]
    maf <- maf[keep2]
    if (nrow(m) < 50) next

    info <- 2 * maf * (1 - maf) * (m$n_eqtl + m$z^2)
    b_eqtl <- m$z_aligned / sqrt(info)
    vb_eqtl <- 1 / info

    d1 <- list(snp = m$snp, beta = b_eqtl, varbeta = vb_eqtl,
               type = "quant", N = as.integer(median(m$n_eqtl)),
               MAF = maf, sdY = 1)

    vb_out <- m$se^2
    d2 <- if (cfg$type == "cc") {
      list(snp = m$snp, beta = m$beta, varbeta = vb_out, type = "cc",
           N = cfg$N, s = cfg$s, MAF = maf)
    } else {
      list(snp = m$snp, beta = m$beta, varbeta = vb_out, type = "quant",
           N = cfg$N, MAF = maf, sdY = 1)
    }

    cr <- tryCatch(
      suppressMessages(suppressWarnings(coloc.abf(d1, d2, p12 = 1e-5))),
      error = function(err) {
        cat(sprintf("  [skip] %s / %s: %s\n", g, oc, conditionMessage(err)))
        NULL
      })
    if (is.null(cr)) next

    s <- cr$summary
    top <- cr$results[which.max(cr$results$SNP.PP.H4), ]
    res[[length(res) + 1]] <- data.table(
      gene = g, outcome = oc, nsnps = as.integer(s["nsnps"]),
      PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]], PP.H2 = s[["PP.H2.abf"]],
      PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]],
      top_snp = top$snp, top_snp_PP = top$SNP.PP.H4,
      pval_eqtl = min(m$pval_eqtl), pval_out = min(m$pval),
      maf_median = median(maf),
      n_snp_eqtl_sig = sum(m$pval_eqtl < 5e-8))
  }
  cat(sprintf("  done: %d genes\n", sum(sapply(res, function(x) x$outcome == oc))))
}

out <- rbindlist(res)
ARGS <- commandArgs(trailingOnly = TRUE)
suffix <- if (length(ARGS)) paste(ARGS, collapse = "_") else "all"
part_path <- sprintf("results/coloc_part_%s.csv", suffix)
fwrite(out, part_path)
cat(sprintf("\n-> %s\n", part_path))

cat("\n===== summary =====\n")
print(out[, .(genes = .N,
              H4_gt_0.8 = sum(PP.H4 > 0.8),
              H4_gt_0.5 = sum(PP.H4 > 0.5),
              H3_gt_0.8 = sum(PP.H3 > 0.8),
              median_H4 = round(median(PP.H4), 3)), by = outcome])

cat("\nstrong colocalization (PP.H4 > 0.8):\n")
print(out[PP.H4 > 0.8][order(-PP.H4),
      .(gene, outcome, nsnps, PP.H4 = round(PP.H4, 3),
        H3 = round(PP.H3, 3), top_snp, top_snp_PP = round(top_snp_PP, 3))])
