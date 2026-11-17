###############################################################################
# Step 12c: multi-signal colocalization with coloc.susie
#
# Aim: coloc.abf assumes one causal variant per trait, so regions carrying
#      several independent signals (HLA, NRG1 and the like) are systematically
#      misread as H3 ("distinct causal variants"). This step first fine-maps the
#      cis region with susie, then colocalizes signal pairs one by one.
#
# Input : tmp/susie/<gene>__<outcome>.tsv      aligned SNP statistics
#         tmp/susie/<gene>__<outcome>.ld.csv.gz 1000G EUR LD (r)
# Output: results/coloc_susie_results.csv
#
# LD source: OpenGWAS LD API (1000 Genomes phase 3, EUR)
# Stability: susie is very sensitive to inconsistency between z and LD, so two
#            passes are attempted:
#            (1) the full SNP set
#            (2) drop redundant SNPs with |r|>0.98 and rebuild both inputs
#            Only if both fail does the pair fall back to abf, tagged abf_fallback
###############################################################################

suppressMessages({
  library(data.table)
  library(coloc)
  library(susieR)
})

BASE <- "D:/NewBeginning/MR1/DATA10"
setwd(BASE)
OUTDIR <- "tmp/susie"
dir.create("results", showWarnings = FALSE)

CFG <- list(
  TSH  = list(type = "quant", N = 482873, label = "TSH"),
  T4   = list(type = "quant", N = 191449, label = "Free T4"),
  HYPO = list(type = "cc", N = 1046991, s = 98805 / 1046991, label = "Hypothyroidism"),
  HASH = list(type = "cc", N = 395640, s = 15654 / 395640, label = "Hashimoto thyroiditis")
)

make_pd <- function(R) {
  R <- (R + t(R)) / 2
  diag(R) <- 1
  ev <- eigen(R, symmetric = TRUE)
  if (min(ev$values) < 1e-6) {
    R <- cov2cor(ev$vectors %*% diag(pmax(ev$values, 1e-6)) %*% t(ev$vectors))
  }
  R
}

fit_susie <- function(d, lam = 0, L = 5, max_iter = 200) {
  # susie_rss is called directly rather than runsusie: the latter defaults to
  # L=10 and max_iter=100 and retries repeatedly on failure.
  # With lam > 0 the LD matrix is diagonally shrunk, which mitigates IBSS
  # non-convergence caused by inconsistency between the external LD reference
  # panel and the summary statistics (common for eQTLGen against 1000G EUR).
  R <- d$LD
  if (lam > 0) {
    R <- (1 - lam) * R + lam * diag(nrow(R))
    dimnames(R) <- dimnames(d$LD)
  }
  z <- d$beta / sqrt(d$varbeta)
  s <- tryCatch(suppressWarnings(susie_rss(z = z, R = R, n = d$N, L = L,
                                           estimate_residual_variance = FALSE,
                                           check_prior = FALSE, max_iter = max_iter)),
                error = function(e) NULL)
  if (!is.null(s) && isTRUE(s$converged)) return(s)
  NULL
}

shrink <- function(R, lam) {
  if (lam <= 0) return(R)
  Rs <- (1 - lam) * R + lam * diag(nrow(R))
  dimnames(Rs) <- dimnames(R)
  Rs
}

## greedily keep a low-redundancy SNP set by significance
prune_ld <- function(R, p, cut = 0.98) {
  o <- order(p)
  keep <- integer(0)
  for (i in o) {
    if (length(keep) == 0 || all(abs(R[i, keep]) <= cut)) keep <- c(keep, i)
  }
  sort(keep)
}

build_d1 <- function(d, R) {
  list(snp = d$snp, beta = d$beta_eqtl, varbeta = d$varbeta_eqtl, type = "quant",
       N = as.integer(median(d$n)), MAF = d$maf, sdY = 1, LD = R)
}
build_d2 <- function(d, R, cfg) {
  if (cfg$type == "cc") {
    list(snp = d$snp, beta = d$beta, varbeta = d$se^2, type = "cc",
         N = cfg$N, s = cfg$s, MAF = d$maf, LD = R)
  } else {
    list(snp = d$snp, beta = d$beta, varbeta = d$se^2, type = "quant",
         N = cfg$N, MAF = d$maf, sdY = 1, LD = R)
  }
}

## returns list(method, n1, n2, nsnps, sm) or NULL
## tries LD shrinkage 0 / 0.05 in turn; both traits must use the same shrunk
## matrix, otherwise coloc is not self-consistent
coloc_pair <- function(d, R, cfg) {
  lams <- c(0, 0.05)
  for (pass in 1:2) {
    if (pass == 2) {
      keep <- prune_ld(R, d$pmin, 0.98)
      if (length(keep) < 20) return(NULL)
      d <- d[keep]
      R <- R[keep, keep]
    }
    for (lam in lams) {
      Rs <- shrink(R, lam)
      d1 <- build_d1(d, Rs)
      d2 <- build_d2(d, Rs, cfg)
      s1 <- fit_susie(d1, lam)
      s2 <- fit_susie(d2, lam)
      if (is.null(s1) || is.null(s2)) next
      n1 <- if (is.null(s1$sets$cs)) 0L else length(s1$sets$cs)
      n2 <- if (is.null(s2$sets$cs)) 0L else length(s2$sets$cs)
      d1$susie.bf <- s1
      d2$susie.bf <- s2
      cs <- NULL
      tryCatch({
        invisible(capture.output(
          cs <- suppressWarnings(suppressMessages(coloc.susie(d1, d2, p12 = 1e-5)))))
      }, error = function(e) NULL)
      if (is.null(cs) || is.null(cs$summary)) next
      meth <- if (lam == 0) "susie" else sprintf("susie_shrink%.2f", lam)
      if (pass == 2) meth <- paste0(meth, "_pruned")
      return(list(method = meth, n1 = n1, n2 = n2, nsnps = nrow(d),
                  sm = as.data.frame(cs$summary)))
    }
  }
  NULL
}

files <- list.files(OUTDIR, pattern = "\\.tsv$", full.names = TRUE)
ARGS <- commandArgs(trailingOnly = TRUE)
CHUNK <- if (length(ARGS) >= 2) as.integer(ARGS[1]) else 0L
NCHUNK <- if (length(ARGS) >= 2) as.integer(ARGS[2]) else 1L
if (NCHUNK > 1) files <- files[seq_along(files) %% NCHUNK == (CHUNK %% NCHUNK)]
# ONLY="CD69__HASH" recomputes a single pair; OUT sets the output path so the
# official chunk files are not overwritten
ONLY <- Sys.getenv("ONLY", "")
if (nzchar(ONLY)) {
  files <- files[grepl(paste0("^", ONLY, "\\.tsv$"), basename(files))]
  cat("[ONLY]", ONLY, "->", length(files), "file(s)\n")
}
cat("to process:", length(files), "gene-outcome pairs (chunk", CHUNK, "/", NCHUNK, ")\n")
flush(stdout())

res <- list()
t_start <- Sys.time()
for (i in seq_along(files)) {
  f <- files[i]
  bn <- basename(f)
  g <- sub("__.*$", "", bn)
  oc <- sub("\\.tsv$", "", sub("^.*__", "", bn))
  cfg <- CFG[[oc]]
  if (is.null(cfg)) next

  ld_f <- sub("\\.tsv$", ".ld.csv.gz", f)
  if (!file.exists(ld_f)) next

  d <- tryCatch(fread(f), error = function(e) NULL)
  R <- tryCatch(as.matrix(fread(ld_f, header = FALSE)), error = function(e) NULL)
  if (is.null(d) || is.null(R)) next
  if (nrow(d) != nrow(R)) { cat("[dimension mismatch]", bn, "\n"); next }
  if (nrow(d) < 20) next

  d$pmin <- pmin(d$p, d$p_o)
  if (nrow(d) > 200) {                    # cap the susie problem by keeping the 200 most significant SNPs from both sides
    keep <- sort(order(d$pmin)[1:200])
    d <- d[keep]
    R <- R[keep, keep]
  }
  R <- make_pd(R)
  dimnames(R) <- list(d$snp, d$snp)

  setTimeLimit(cpu = Inf, elapsed = 25, transient = TRUE)
  cp <- tryCatch(coloc_pair(d, R, cfg), error = function(e) NULL)
  setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE)

  if (is.null(cp)) {
    ab <- suppressWarnings(suppressMessages(coloc.abf(
      build_d1(d, R)[c("snp", "beta", "varbeta", "type", "N", "MAF", "sdY")],
      build_d2(d, R, cfg)[intersect(c("snp", "beta", "varbeta", "type", "N", "MAF", "sdY", "s"),
                                    names(build_d2(d, R, cfg)))],
      p12 = 1e-5)))
    sm <- ab$summary
    res[[length(res) + 1]] <- data.table(
      gene = g, outcome = oc, method = "abf_fallback", nsnps = nrow(d),
      n_sig_eqtl = NA_integer_, n_sig_out = NA_integer_, n_pairs = 0L,
      PP.H0 = sm[["PP.H0.abf"]], PP.H1 = sm[["PP.H1.abf"]], PP.H2 = sm[["PP.H2.abf"]],
      PP.H3 = sm[["PP.H3.abf"]], PP.H4 = sm[["PP.H4.abf"]],
      top_snp = ab$results$snp[which.max(ab$results$SNP.PP.H4)])
    next
  }

  sm <- cp$sm
  h4col <- if ("PP.H4.abf" %in% names(sm)) "PP.H4.abf" else grep("^PP\\.H4", names(sm), value = TRUE)[1]
  h3col <- if ("PP.H3.abf" %in% names(sm)) "PP.H3.abf" else grep("^PP\\.H3", names(sm), value = TRUE)[1]
  k <- which.max(sm[[h4col]])
  pick <- function(kk) {
    cand <- if (paste0(kk, ".abf") %in% names(sm)) paste0(kk, ".abf") else kk
    if (cand %in% names(sm)) as.numeric(sm[[cand]][k]) else NA_real_
  }

  res[[length(res) + 1]] <- data.table(
    gene = g, outcome = oc, method = cp$method, nsnps = cp$nsnps,
    n_sig_eqtl = cp$n1, n_sig_out = cp$n2, n_pairs = nrow(sm),
    PP.H0 = pick("PP.H0"), PP.H1 = pick("PP.H1"), PP.H2 = pick("PP.H2"),
    PP.H3 = as.numeric(sm[[h3col]][k]), PP.H4 = as.numeric(sm[[h4col]][k]),
    top_snp = paste(unique(c(sm$hit1, sm$hit2)), collapse = ";"))

  if (i %% 10 == 0) {
    cat(sprintf("  %d/%d  %.0fs (%s/%s: %d x %d signals)\n", i, length(files),
                as.numeric(difftime(Sys.time(), t_start, units = "secs")),
                g, oc, cp$n1, cp$n2))
    flush(stdout())
  }
}

out <- rbindlist(res, fill = TRUE)
outfile <- Sys.getenv("OUT", "")
if (!nzchar(outfile)) {
  outfile <- if (NCHUNK > 1) sprintf("results/coloc_susie_part_%d-of-%d.csv", CHUNK, NCHUNK) else "results/coloc_susie_results.csv"
}
fwrite(out, outfile)
cat("written:", outfile, "(", nrow(out), "rows )\n")
cat(sprintf("\n-> %s  (%d rows)\n", outfile, nrow(out)))

cat("\n===== susie results =====\n")
print(out[, .(n = .N,
              H4_gt_0.8 = sum(PP.H4 > 0.8, na.rm = TRUE),
              H3_gt_0.8 = sum(PP.H3 > 0.8, na.rm = TRUE),
              median_H4 = round(median(PP.H4, na.rm = TRUE), 3)),
          by = method])
