#!/usr/bin/env Rscript
# p13_susie_aligned.R —— 用"等位基因对齐"的 LD 重跑 SuSiE + coloc.susie（甲状腺）
suppressMessages({library(coloc); library(data.table)})
BASE <- "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
OUT <- file.path(BASE, "tables")

rd_ld <- function(f) {
  x <- fread(file.path(BASE, "inputs_ld", f), sep = "\t", header = TRUE)
  rn <- as.character(x[[1]])
  M <- as.matrix(x[, -1, with = FALSE]); rownames(M) <- rn; colnames(M) <- rn
  M
}
L_out <- rd_ld("ld_TSH_aligned.tsv.gz")
L_eqtl <- rd_ld("ld_eqtl_aligned.tsv.gz")
cat("LD:", dim(L_out), dim(L_eqtl), "\n")

srows <- list(); csrows <- list()
for (g in c("SCRN2", "TBKBP1", "MRPL10")) {
  f <- file.path(BASE, sprintf("inputs/Thyroid__%s__TSH.tsv", g))
  d <- fread(f, sep = "\t")
  d <- d[!is.na(beta_eqtl) & !is.na(se_eqtl) & se_eqtl > 0 & se_out > 0]
  d$snp <- as.character(d$pos37)
  d <- d[as.numeric(d$pos37) >= 45800299 - 150000 & as.numeric(d$pos37) <= 45800299 + 150000]
  com <- Reduce(intersect, list(d$snp, rownames(L_out), rownames(L_eqtl)))
  cat(g, "共用 SNP:", length(com), "\n")
  d <- as.data.frame(d); d <- d[match(com, d$snp), ]
  Lo <- L_out[com, com]; Le <- L_eqtl[com, com]
  d1 <- list(snp = d$snp, beta = d$beta_eqtl, varbeta = d$se_eqtl^2, type = "quant",
             sdY = 1, MAF = d$maf, N = 574L, LD = Le)
  d2 <- list(snp = d$snp, beta = d$beta_out, varbeta = d$se_out^2, type = "quant",
             sdY = 1, MAF = d$maf, N = as.integer(median(d$n_out, na.rm = TRUE)), LD = Lo)
  s1 <- tryCatch(suppressWarnings(runsusie(d1, max_iter = 200)), error = function(e) NULL)
  s2 <- tryCatch(suppressWarnings(runsusie(d2, max_iter = 200)), error = function(e) NULL)
  if (is.null(s1) || is.null(s2)) { cat("susie fail", g, "\n"); next }
  for (src in c("eQTL", "TSH")) {
    s <- if (src == "eQTL") s1 else s2
    if (is.null(s$sets$cs)) next
    for (k in seq_along(s$sets$cs)) {
      idx <- s$sets$cs[[k]]; top <- idx[which.max(s$pip[idx])]
      csrows[[length(csrows) + 1]] <- data.table(
        gene = g, source = src, cs = names(s$sets$cs)[k], cs_size = length(idx),
        top_pos37 = as.numeric(d$snp[top]), top_pip = max(s$pip[idx]),
        top_p = if (src == "eQTL") d$p_eqtl[top] else d$p_out[top],
        top_beta = if (src == "eQTL") d$beta_eqtl[top] else d$beta_out[top])
    }
  }
  cr <- tryCatch(suppressWarnings(coloc.susie(s1, s2)), error = function(e) NULL)
  if (!is.null(cr$summary) && nrow(cr$summary) > 0) {
    sm <- as.data.table(cr$summary); sm[, gene := g]
    sm[, `:=`(n_cs_eqtl = length(s1$sets$cs), n_cs_tsh = length(s2$sets$cs))]
    srows[[length(srows) + 1]] <- sm
  } else {
    srows[[length(srows) + 1]] <- data.table(gene = g, n_cs_eqtl = length(s1$sets$cs),
      n_cs_tsh = length(s2$sets$cs), note = "无可配对信号")
  }
  cat("done", g, "\n")
}
ss <- rbindlist(srows, fill = TRUE); fwrite(ss, file.path(OUT, "03b_面板_coloc_susie_TSH.csv"))
cs <- rbindlist(csrows, fill = TRUE); fwrite(cs, file.path(OUT, "03c_可信集_TSH.csv"))
print(ss)
cat("\n--- 可信集（甲状腺，等位基因对齐的 LD）---\n")
print(cs[, .(gene, source, cs, cs_size, top_pos37, top_pip, top_p)][order(gene, source)])
