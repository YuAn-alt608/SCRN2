###############################################################################
# Step 21: coloc 先验敏感性 —— p12 = 1e-4 / 1e-5 / 1e-6
#
# 目的：PP.H4 是贝叶斯后验，强依赖"共享因果变异的先验概率" p12。
#       主分析用 coloc 默认 1e-5，本步检验 71 个 PP.H4>0.8 是否只是先验的产物。
#
# 做法：数据准备只做一次（合并、等位基因对齐、MAF 过滤、beta 还原），
#       同一份 d1/d2 分别喂给三个 p12，保证差异只来自先验。
#       同时用 p12=1e-5 的输出与 results/coloc_results.csv 对账（应完全一致）。
#
# 输入：data/eqtl_regions_all.txt.gz, data/region_*.tsv.gz, data/coloc_genes_all.tsv
# 输出：results/coloc_p12_long.csv（长表：gene, outcome, p12, PP.H3, PP.H4, top_snp, top_snp_PP）
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
cat(sprintf("eQTL 记录 %s 行，%d 个基因\n", format(nrow(eq), big.mark = ","),
            uniqueN(eq$gene)))

res <- list()

for (oc in names(OUTCOMES)) {
  cfg <- OUTCOMES[[oc]]
  o <- fread(cfg$file, showProgress = FALSE)
  setnames(o, c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                "beta", "se", "eaf", "pval", "rsid", "n"), skip_absent = TRUE)
  o[, base_pair_location := as.integer(base_pair_location)]
  setkey(o, chromosome, base_pair_location)
  cat(sprintf("\n=== %s: %s 条区域变异 ===\n", cfg$label, format(nrow(o), big.mark = ",")))

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
  cat(sprintf("  %s 完成，累计 %d 行\n", oc, length(res)))
  flush.console()
}

out <- rbindlist(res)
fwrite(out, "results/coloc_p12_long.csv")
cat(sprintf("\n-> results/coloc_p12_long.csv  %d 行\n", nrow(out)))

cat("\n===== 各先验下 PP.H4>0.8 的对数 =====\n")
print(out[, .(n_pairs = .N, H4_gt_0.8 = sum(PP.H4 > 0.8),
              H4_gt_0.5 = sum(PP.H4 > 0.5)), by = .(p12)][order(-p12)])

cat("\n===== 按结局 =====\n")
print(dcast(out[PP.H4 > 0.8, .N, by = .(outcome, p12)], outcome ~ p12, value.var = "N"))

# 与主分析对账
main <- fread("results/coloc_results.csv")
chk <- merge(out[p12 == 1e-5], main[, .(gene, outcome, PP.H4_main = PP.H4)],
             by = c("gene", "outcome"))
chk[, d := abs(PP.H4 - PP.H4_main)]
cat(sprintf("\n对账：%d 对，PP.H4 最大绝对差 %.3g（应为 0）\n", nrow(chk), max(chk$d)))
