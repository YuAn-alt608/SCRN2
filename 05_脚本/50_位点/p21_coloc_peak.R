#!/usr/bin/env Rscript
# p21_coloc_peak.R —— 45,873,049 这个未归属甲减峰周围 ±200 kb 的基因，全部做窗口 coloc
#
# 背景：甲减在信号 2 内部另有一个峰在 45,873,049（P=3.09e-19，全区域次强），
# 之前只查了 SCRN2/TBKBP1/MRPL10 三个基因，都没 eQTL，于是标成"未归属"。
# 但那三个不是峰周围唯一的基因 —— ±200 kb 内有 17 个基因有 eQTL 数据。
# 本脚本把 17 个基因 × 2 组织 × 2 结局 × 2 个窗口全部跑完，看峰能不能被归属。
suppressMessages({library(coloc); library(data.table)})

BASE <- "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
IN <- file.path(BASE, "inputs"); OUT <- file.path(BASE, "tables")

WIN <- list(
  "峰±50kb(45.823-45.923)"  = c(45.823e6, 45.923e6),
  "峰±100kb(45.773-45.973)" = c(45.773e6, 45.973e6)
)

GENES <- c("NPEPPS", "RP11-580I16.2", "KPNB1", "TBKBP1", "TBX21", "OSBPL7",
           "MRPL10", "LRRC46", "SCRN2", "SP6", "RP11-6N17.4", "SP2",
           "AC003665.1", "RP11-6N17.3", "PNPO", "PRR15L", "CDK5RAP3")
OCS <- c("TSH", "HYPO")

rows <- list()
for (gene in GENES) for (oc in OCS) for (tissue in c("Thyroid", "WholeBlood")) {
  f <- file.path(IN, sprintf("%s__%s__%s.tsv", tissue, gene, oc))
  if (!file.exists(f)) next
  d <- fread(f, sep = "\t")
  d <- d[!is.na(beta_eqtl) & !is.na(se_eqtl) & se_eqtl > 0 &
         !is.na(se_out) & se_out > 0 & !is.na(maf)]
  if (nrow(d) < 50) next
  lead <- d[which.min(p_eqtl)]
  for (wn in names(WIN)) {
    r0 <- WIN[[wn]]
    s <- d[pos37 >= r0[1] & pos37 <= r0[2]]
    if (nrow(s) < 40) next
    d1 <- list(snp = as.character(s$pos37), beta = s$beta_eqtl, varbeta = s$se_eqtl^2,
               type = "quant", sdY = 1, MAF = s$maf, N = as.integer(median(s$n_eqtl)))
    d2 <- list(snp = as.character(s$pos37), beta = s$beta_out, varbeta = s$se_out^2,
               type = "quant", sdY = 1, MAF = s$maf,
               N = as.integer(median(s$n_out, na.rm = TRUE)))
    r <- tryCatch(suppressWarnings(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-4)),
                  error = function(e) NULL)
    if (is.null(r)) next
    x <- r$summary
    rows[[length(rows) + 1]] <- data.table(
      gene = gene, tissue = tissue, outcome = oc, window = wn, nsnp = as.integer(x["nsnps"]),
      H0 = x[["PP.H0.abf"]], H1 = x[["PP.H1.abf"]], H2 = x[["PP.H2.abf"]],
      H3 = x[["PP.H3.abf"]], H4 = x[["PP.H4.abf"]],
      eqtl_lead_pos = as.integer(lead$pos37), eqtl_lead_p = lead$p_eqtl,
      eqtl_lead_beta = lead$beta_eqtl)
  }
}

res <- rbindlist(rows)
fwrite(res, file.path(OUT, "26_峰周围基因_窗口coloc.csv"))
cat("完成，行数:", nrow(res), " 基因:", length(unique(res$gene)), "\n\n")
cat("=== 甲减，峰±50kb 窗口 ===\n")
print(dcast(res[outcome == "HYPO" & window == names(WIN)[1]],
            gene ~ tissue, value.var = "H4"), digits = 3)
cat("\n=== TSH，峰±50kb 窗口 ===\n")
print(dcast(res[outcome == "TSH" & window == names(WIN)[1]],
            gene ~ tissue, value.var = "H4"), digits = 3)
