#!/usr/bin/env Rscript
# p20_coloc_window.R —— 按信号窗口重跑 coloc.abf
#
# 动机：coloc.abf 假设每个性状只有一个因果变异。本区域甲减有三个信号，
# 全窗口跑时它会锁定最强的那一个（信号 1，45.80），于是 SCRN2 的 H4 被压到 0.00。
# 把窗口限制在信号 2 内（45.86–45.98 Mb）再看，才是"信号 2 与甲减共不共享"的正确问法。
suppressMessages({library(coloc); library(data.table)})

BASE <- "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
IN <- file.path(BASE, "inputs"); OUT <- file.path(BASE, "tables")

WIN <- list(
  "全窗口(45.34-46.54)" = c(45.34e6, 46.54e6),
  "信号1(45.72-45.86)"  = c(45.72e6, 45.86e6),
  "信号2(45.86-45.98)"  = c(45.86e6, 45.98e6),
  "信号3(46.00-46.12)"  = c(46.00e6, 46.12e6)
)
GENES <- c("SCRN2", "TBKBP1", "MRPL10", "COPZ2")
OCS <- c("TSH", "HYPO", "HASH", "T4")

rows <- list()
for (gene in GENES) for (oc in OCS) for (tissue in c("Thyroid", "WholeBlood")) {
  f <- file.path(IN, sprintf("%s__%s__%s.tsv", tissue, gene, oc))
  if (!file.exists(f)) next
  d <- fread(f, sep = "\t")
  d <- d[!is.na(beta_eqtl) & !is.na(se_eqtl) & se_eqtl > 0 & !is.na(se_out) & se_out > 0 & !is.na(maf)]
  if (nrow(d) < 50) next
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
      H3 = x[["PP.H3.abf"]], H4 = x[["PP.H4.abf"]])
  }
}

res <- rbindlist(rows)
fwrite(res, file.path(OUT, "24_按信号窗口的coloc.csv"))
cat("完成，行数:", nrow(res), "\n\n")
key <- res[gene == "SCRN2" & outcome %in% c("TSH", "HYPO")]
print(dcast(key, gene + tissue + outcome ~ window, value.var = "H4"), digits = 3)
cat("\n--- TBKBP1 ---\n")
key2 <- res[gene == "TBKBP1" & outcome %in% c("TSH", "HYPO")]
print(dcast(key2, gene + tissue + outcome ~ window, value.var = "H4"), digits = 3)
