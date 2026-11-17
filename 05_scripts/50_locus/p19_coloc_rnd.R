#!/usr/bin/env Rscript
# p19_coloc_rnd.R -- re-run coloc.abf with the Rand 2025 hypothyroidism outcome (same priors and LD conventions)
suppressMessages({library(coloc); library(data.table)})

BASE <- "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
IN   <- file.path(BASE, "inputs")
OUT  <- file.path(BASE, "tables")
P12  <- c(1e-6, 1e-5, 1e-4, 1e-3)

files <- list.files(IN, pattern = "__RHYPO\\.tsv$", full.names = TRUE)
cat("RHYPO input combinations:", length(files), "\n")

rows <- list()
for (f in files) {
  p <- strsplit(sub("\\.tsv$", "", basename(f)), "__")[[1]]
  tissue <- p[1]; gene <- p[2]
  d <- fread(f, sep = "\t")
  d <- d[!is.na(beta_eqtl) & !is.na(se_eqtl) & se_eqtl > 0 &
         !is.na(se_out) & se_out > 0 & !is.na(maf)]
  if (nrow(d) < 50) next
  d1 <- list(snp = as.character(d$pos37), beta = d$beta_eqtl, varbeta = d$se_eqtl^2,
             type = "quant", sdY = 1, MAF = d$maf, N = as.integer(median(d$n_eqtl)))
  d2 <- list(snp = as.character(d$pos37), beta = d$beta_out, varbeta = d$se_out^2,
             type = "quant", sdY = 1, MAF = d$maf, N = 1046991L)
  for (p12 in P12) {
    r <- tryCatch(suppressWarnings(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12)),
                  error = function(e) NULL)
    if (is.null(r)) next
    s <- r$summary
    rows[[length(rows) + 1]] <- data.table(
      tissue = tissue, gene = gene, outcome = "RHYPO", p12 = p12,
      nsnp = as.integer(s["nsnps"]),
      PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]], PP.H2 = s[["PP.H2.abf"]],
      PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]])
  }
}
res <- rbindlist(rows)
fwrite(res, file.path(OUT, "24_coloc_rand2025_hypothyroidism.csv"))
cat("done, rows:", nrow(res), "\n")
print(res[p12 == 1e-4 & gene %in% c("SCRN2", "TBKBP1", "COPZ2")])
