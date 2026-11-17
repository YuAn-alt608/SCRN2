#!/usr/bin/env Rscript
# p21_coloc_peak.R -- window coloc for every gene within +/-200 kb of the
# unassigned hypothyroidism peak at 45,873,049
#
# Background: inside signal 2, hypothyroidism has a further peak at 45,873,049
# (P=3.09e-19, the second strongest in the region). Only SCRN2, TBKBP1 and MRPL10
# were examined before, and none of them has an eQTL there, so the peak was
# labelled "unassigned". Those three are not the only genes around the peak:
# 17 genes within +/-200 kb have eQTL data. This script runs all 17 genes x 2
# tissues x 2 outcomes x 2 windows to see whether the peak can be attributed.
suppressMessages({library(coloc); library(data.table)})

BASE <- "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
IN <- file.path(BASE, "inputs"); OUT <- file.path(BASE, "tables")

WIN <- list(
  "Peak +/-50 kb (45.823-45.923)"  = c(45.823e6, 45.923e6),
  "Peak +/-100 kb (45.773-45.973)" = c(45.773e6, 45.973e6)
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
fwrite(res, file.path(OUT, "26_peak_neighbour_genes_window_coloc.csv"))
cat("done, rows:", nrow(res), " genes:", length(unique(res$gene)), "\n\n")
cat("=== hypothyroidism, peak +/-50 kb window ===\n")
print(dcast(res[outcome == "HYPO" & window == names(WIN)[1]],
            gene ~ tissue, value.var = "H4"), digits = 3)
cat("\n=== TSH, peak +/-50 kb window ===\n")
print(dcast(res[outcome == "TSH" & window == names(WIN)[1]],
            gene ~ tissue, value.var = "H4"), digits = 3)
