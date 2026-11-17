#!/usr/bin/env Rscript
# p20_coloc_window.R -- re-run coloc.abf within each signal window
#
# Motivation: coloc.abf assumes a single causal variant per trait. Hypothyroidism
# has three signals in this region, so over the full window it locks onto the
# strongest one (signal 1, 45.80) and the SCRN2 H4 collapses to 0.00. Restricting
# the window to signal 2 (45.86-45.98 Mb) poses the right question: does signal 2
# share a causal variant with hypothyroidism?
suppressMessages({library(coloc); library(data.table)})

BASE <- "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
IN <- file.path(BASE, "inputs"); OUT <- file.path(BASE, "tables")

WIN <- list(
  "Full window (45.34-46.54)" = c(45.34e6, 46.54e6),
  "Signal 1 (45.72-45.86)"    = c(45.72e6, 45.86e6),
  "Signal 2 (45.86-45.98)"    = c(45.86e6, 45.98e6),
  "Signal 3 (46.00-46.12)"    = c(46.00e6, 46.12e6)
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
fwrite(res, file.path(OUT, "24_coloc_by_signal_window.csv"))
cat("done, rows:", nrow(res), "\n\n")
key <- res[gene == "SCRN2" & outcome %in% c("TSH", "HYPO")]
print(dcast(key, gene + tissue + outcome ~ window, value.var = "H4"), digits = 3)
cat("\n--- TBKBP1 ---\n")
key2 <- res[gene == "TBKBP1" & outcome %in% c("TSH", "HYPO")]
print(dcast(key2, gene + tissue + outcome ~ window, value.var = "H4"), digits = 3)
