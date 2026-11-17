# =============================================================================
#  Step 03: cis-eQTL 两样本孟德尔随机化
#  - 单 SNP: Wald ratio
#  - 多 SNP: IVW / MR-Egger / 加权中位数 / 加权众数 (TwoSampleMR)
#  - 逐基因输出 beta, se, p, F, R2；结局内 BH-FDR
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

BASE  <- "D:/NewBeginning/MR1/DATA10"
DATA  <- file.path(BASE, "data")
RES   <- file.path(BASE, "results")
dir.create(RES, showWarnings = FALSE)

N_EQTL <- 31684

ivs <- fread(file.path(DATA, "ivs.tsv"))
setnames(ivs, c("gene","snp","chr","pos","ea","nea","z","p","f_stat"))
cat(sprintf("工具变量: %d SNP / %d 基因\n", nrow(ivs), uniqueN(ivs$gene)))

OUTCOMES <- list(
  TSH  = list(file = "out_TSH.tsv",  type = "continuous", n = 482873,  label = "TSH"),
  T4   = list(file = "out_T4.tsv",   type = "continuous", n = 191449,  label = "Total T4"),
  HYPO = list(file = "out_HYPO.tsv", type = "binary",     ncase = 98805, ncont = 948186,
              label = "Hypothyroidism"),
  HASH = list(file = "out_HASH.tsv", type = "binary",     ncase = 15654, ncont = 379986,
              label = "Hashimoto thyroiditis")
)

run_one <- function(key, meta) {
  path <- file.path(DATA, meta$file)
  if (!file.exists(path)) { cat(sprintf("[skip] %s\n", key)); return(NULL) }

  o <- fread(path)
  o <- o[, .(rsid, ea_o = ea, nea_o = nea, beta_o = beta, se_o = se,
             pval_o = pval, eaf_o = eaf, n_o = n)]

  d <- merge(ivs, o, by.x = "snp", by.y = "rsid")
  cat(sprintf("[%s] 可匹配 SNP: %d\n", key, nrow(d)))
  if (nrow(d) == 0) return(NULL)

  # ---- 等位基因对齐 -----------------------------------------------------------
  same <- d$ea == d$ea_o & d$nea == d$nea_o
  flip <- d$ea == d$nea_o & d$nea == d$ea_o
  keep <- same | flip
  d <- d[keep]
  flip2 <- flip[keep]
  d[flip2, `:=`(beta_o = -beta_o, eaf_o = 1 - eaf_o)]
  d[, ea_o := NULL][, nea_o := NULL]

  # 回文位点（A/T、C/G）且频率接近 0.5 者剔除
  pal <- (toupper(d$ea) == "A" & toupper(d$nea) == "T") | (toupper(d$ea) == "T" & toupper(d$nea) == "A") |
         (toupper(d$ea) == "C" & toupper(d$nea) == "G") | (toupper(d$ea) == "G" & toupper(d$nea) == "C")
  d <- d[!(pal & eaf_o > 0.42 & eaf_o < 0.58)]
  cat(sprintf("[%s] 对齐后 SNP: %d\n", key, nrow(d)))

  # ---- 暴露效应量：由 Z 还原 beta/se -----------------------------------------
  # beta = Z / sqrt(2p(1-p)N)，p 取结局 EAF（与暴露同人群的近似）
  d[, p_use := pmin(pmax(eaf_o, 0.01), 0.99)]
  d[, scale := sqrt(2 * p_use * (1 - p_use) * N_EQTL)]
  d[, beta_x := z / scale]
  d[, se_x := 1 / scale]
  d[, r2 := f_stat / (f_stat + N_EQTL)]

  # ---- MR --------------------------------------------------------------------
  single <- d[, .N, by = gene][N == 1]$gene
  multi  <- d[, .N, by = gene][N >= 2]$gene

  res_single <- d[gene %in% single, {
    b <- beta_o / beta_x
    s <- se_o / abs(beta_x)
    list(nsnp = 1L, method = "Wald ratio", b = b, se = s,
         pval = 2 * pnorm(-abs(b / s)), f_stat = f_stat, r2 = r2)
  }, by = gene]

  res_multi <- NULL
  if (length(multi) > 0) {
    res_multi <- rbindlist(lapply(multi, function(g) {
      dd <- d[gene == g]
      ex <- data.table(SNP = dd$snp, beta.exposure = dd$beta_x, se.exposure = dd$se_x,
                       effect_allele.exposure = dd$ea, other_allele.exposure = dd$nea,
                       eaf.exposure = dd$p_use, exposure = g, id.exposure = g)
      ou <- data.table(SNP = dd$snp, beta.outcome = dd$beta_o, se.outcome = dd$se_o,
                       effect_allele.outcome = dd$ea, other_allele.outcome = dd$nea,
                       eaf.outcome = dd$p_use, outcome = meta$label, id.outcome = meta$label)
      h <- tryCatch(harmonise_data(ex, ou, action = 2), error = function(e) NULL)
      if (is.null(h) || nrow(h) < 2) return(NULL)
      m <- tryCatch(mr(h, method_list = c("mr_ivw", "mr_egger_regression",
                                          "mr_weighted_median", "mr_weighted_mode")),
                    error = function(e) NULL)
      if (is.null(m)) return(NULL)
      m$gene <- g
      m$f_stat <- max(dd$f_stat)
      m$r2 <- max(dd$r2)
      data.table(gene = m$gene, nsnp = m$nsnp, method = m$method,
                 b = m$b, se = m$se, pval = m$pval, f_stat = m$f_stat, r2 = m$r2)
    }), fill = TRUE)
  }

  out <- rbind(res_single, res_multi, fill = TRUE)
  out[, outcome := meta$label]
  out[, fdr := p.adjust(pval, method = "BH"), by = method]
  out
}

all_res <- list()
for (k in names(OUTCOMES)) {
  r <- run_one(k, OUTCOMES[[k]])
  if (!is.null(r)) all_res[[k]] <- r
}
res <- rbindlist(all_res, fill = TRUE)
fwrite(res, file.path(RES, "MR_all_results.csv"))

# 主分析（IVW）汇总
main <- res[method %in% c("Inverse variance weighted", "Wald ratio")]
fwrite(main[order(outcome, pval)], file.path(RES, "MR_main_IVW.csv"))

cat("\n===== 各结局 IVW 结果 =====\n")
print(main[, .(n_genes = .N,
               n_p05 = sum(pval < 0.05, na.rm = TRUE),
               n_fdr05 = sum(fdr < 0.05, na.rm = TRUE),
               min_p = suppressWarnings(min(pval, na.rm = TRUE)),
               n_na = sum(is.na(pval))), by = outcome])

# 显著基因（若目标文件被 Excel 等程序占用，自动降级为 _v2 文件名）
sig <- main[!is.na(fdr) & fdr < 0.05][order(outcome, pval)]
sig_path <- file.path(RES, "MR_significant_FDR05.csv")
ok <- tryCatch({ fwrite(sig, sig_path); TRUE },
               error = function(e) FALSE)
if (!ok) {
  sig_path <- file.path(RES, "MR_significant_FDR05_v2.csv")
  fwrite(sig, sig_path)
  cat("[提示] 原文件被占用，已写入:", basename(sig_path), "\n")
}
cat(sprintf("\nFDR<0.05 命中: %d 条 -> %s\n", nrow(sig), basename(sig_path)))
if (nrow(sig) > 0) print(head(sig[, .(outcome, gene, nsnp, b, se, pval, fdr)], 25))
