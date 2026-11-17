#!/usr/bin/env Rscript
# p02_coloc.R -- panel colocalization: coloc.abf (with a p12 prior grid) + coloc.susie (TSH)
suppressMessages({library(coloc); library(data.table)})

BASE <- "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
IN   <- file.path(BASE, "inputs")
OUT  <- file.path(BASE, "tables")
files <- list.files(IN, pattern = "\\.tsv$", full.names = TRUE)
cat("input combinations:", length(files), "\n")

P12 <- c(1e-6, 1e-5, 1e-4, 1e-3)
ld_cache <- list()
get_ld <- function(oc) {
  if (!is.null(ld_cache[[oc]])) return(ld_cache[[oc]])
  f <- file.path(BASE, "inputs_ld", paste0("ld_", oc, ".tsv.gz"))
  if (!file.exists(f)) return(NULL)
  m <- as.matrix(fread(f, sep = "\t", header = TRUE))
  rn <- as.integer(sub("^pos37", "", rownames(m)))
  if (any(is.na(rn))) rn <- as.integer(colnames(m))
  rownames(m) <- as.character(rn); colnames(m) <- as.character(rn)
  ld_cache[[oc]] <<- m
  m
}

rows <- list(); srows <- list()
for (f in files) {
  tag <- sub("\\.tsv$", "", basename(f))
  p <- strsplit(tag, "__")[[1]]
  tissue <- p[1]; gene <- p[2]; oc <- p[3]
  d <- fread(f, sep = "\t")
  d <- d[!is.na(beta_eqtl) & !is.na(se_eqtl) & se_eqtl > 0 & se_out > 0]
  if (nrow(d) < 50) next
  d$snp <- as.character(d$pos37)
  d1 <- list(snp = d$snp, beta = d$beta_eqtl, varbeta = d$se_eqtl^2,
             type = "quant", sdY = 1, MAF = d$maf, N = as.integer(median(d$n_eqtl)))
  d2 <- list(snp = d$snp, beta = d$beta_out, varbeta = d$se_out^2,
             type = "quant", sdY = 1, MAF = d$maf, N = as.integer(median(d$n_out, na.rm = TRUE)))
  if (is.na(d2$N)) d2$N <- 482873L
  for (p12 in P12) {
    r <- tryCatch(suppressWarnings(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12)),
                  error = function(e) NULL)
    if (is.null(r)) next
    s <- r$summary
    rows[[length(rows) + 1]] <- data.table(
      tissue = tissue, gene = gene, outcome = oc, p12 = p12, nsnp = as.integer(s["nsnps"]),
      PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]], PP.H2 = s[["PP.H2.abf"]],
      PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]])
  }
}

res <- rbindlist(rows)
fwrite(res, file.path(OUT, "03_panel_coloc_abf_prior_grid.csv"))

# ---- coloc.susie (TSH only: the primary outcome) ----
tsh <- files[grepl("__TSH\\.tsv$", files)]
for (f in tsh) {
  tag <- sub("\\.tsv$", "", basename(f)); p <- strsplit(tag, "__")[[1]]
  tissue <- p[1]; gene <- p[2]
  d <- fread(f, sep = "\t")
  d <- d[!is.na(beta_eqtl) & !is.na(se_eqtl) & se_eqtl > 0 & se_out > 0]
  if (nrow(d) < 200) next
  LD <- get_ld("TSH"); if (is.null(LD)) next
  d$snp <- as.character(d$pos37)
  com <- intersect(d$snp, rownames(LD))
  if (length(com) < 200) next
  d <- d[d$snp %in% com, ]; setkey(d, snp); d <- d[com, on = "snp"]
  L <- LD[com, com]
  d1 <- list(snp = d$snp, beta = d$beta_eqtl, varbeta = d$se_eqtl^2, type = "quant",
             sdY = 1, MAF = d$maf, N = as.integer(median(d$n_eqtl)), LD = L)
  d2 <- list(snp = d$snp, beta = d$beta_out, varbeta = d$se_out^2, type = "quant",
             sdY = 1, MAF = d$maf,
             N = ifelse(is.na(median(d$n_out, na.rm = TRUE)), 482873L,
                        as.integer(median(d$n_out, na.rm = TRUE))), LD = L)
  out <- tryCatch({
    s1 <- suppressWarnings(runsusie(d1)); s2 <- suppressWarnings(runsusie(d2))
    cr <- suppressWarnings(coloc.susie(s1, s2))
    if (is.null(cr$summary)) stop("no summary")
    sm <- as.data.table(cr$summary)
    sm[, `:=`(tissue = tissue, gene = gene, outcome = "TSH",
              n_cs1 = length(s1$sets$cs), n_cs2 = length(s2$sets$cs))]
    sm
  }, error = function(e) data.table(tissue = tissue, gene = gene, outcome = "TSH",
                                     n_cs1 = NA_integer_, n_cs2 = NA_integer_,
                                     err = conditionMessage(e)))
  srows[[length(srows) + 1]] <- out
  cat("susie:", tag, "ok\n")
}
if (length(srows)) {
  ss <- rbindlist(srows, fill = TRUE)
  fwrite(ss, file.path(OUT, "03b_panel_coloc_susie_TSH.csv"))
}
cat("done\n")
