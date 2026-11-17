###############################################################################
# Step 12c: 多信号共定位 coloc.susie
#
# 目的：coloc.abf 假设每性状仅 1 个因果变异，HLA、NRG1 这类存在多个独立
#       信号的区间会被系统性误判为 H3（"不同因果变异"）。本步在同一个 cis
#       区间内先用 susie 做细定位（fine-mapping），再逐对信号做共定位。
#
# 输入：tmp/susie/<gene>__<outcome>.tsv      对齐后的 SNP 统计量
#       tmp/susie/<gene>__<outcome>.ld.csv.gz 1000G EUR LD（r）
# 输出：results/coloc_susie_results.csv
#
# LD 来源：OpenGWAS LD API（1000 Genomes phase 3, EUR）
# 稳定性：susie 对 z 与 LD 的自洽性极敏感，故按两级尝试
#         (1) 全 SNP 集
#         (2) 剔除 |r|>0.98 的冗余 SNP 后两侧同时重建输入
#         都失败才退回 abf，并记为 abf_fallback
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
  T4   = list(type = "quant", N = 191449, label = "Total T4"),
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
  # 直接用 susie_rss 而非 runsusie：后者默认 L=10、max_iter=100 且失败时反复重试
  # lam > 0 时对 LD 做对角收缩，缓解外部 LD 参考面板与汇总统计量不自洽导致的
  # IBSS 不收敛（eQTLGen 与 1000G EUR 的常见问题）
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

## 按显著性贪心保留低冗余 SNP
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

## 返回 list(method, n1, n2, sm) 或 NULL
## 依次尝试 LD 收缩强度 0 / 0.05 / 0.15；两侧必须用同一收缩矩阵，否则 coloc 不自洽
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
# ONLY="CD69__HASH" 只补算单对；OUT 指定输出路径，避免污染正式分块
ONLY <- Sys.getenv("ONLY", "")
if (nzchar(ONLY)) {
  files <- files[grepl(paste0("^", ONLY, "\\.tsv$"), basename(files))]
  cat("[ONLY]", ONLY, "->", length(files), "个文件\n")
}
cat("待处理:", length(files), "个基因-结局对 (chunk", CHUNK, "/", NCHUNK, ")\n")
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
  if (nrow(d) != nrow(R)) { cat("[维度不符]", bn, "\n"); next }
  if (nrow(d) < 20) next

  d$pmin <- pmin(d$p, d$p_o)
  if (nrow(d) > 200) {                    # 控制 susie 规模，保留双侧最显著的 200 个
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
    cat(sprintf("  %d/%d  %.0fs（%s/%s: %d×%d 信号）\n", i, length(files),
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
cat("已写出:", outfile, "(", nrow(out), "行 )\n")
cat(sprintf("\n-> %s  (%d 行)\n", outfile, nrow(out)))

cat("\n===== susie 结论 =====\n")
print(out[, .(n = .N,
              H4_gt_0.8 = sum(PP.H4 > 0.8, na.rm = TRUE),
              H3_gt_0.8 = sum(PP.H3 > 0.8, na.rm = TRUE),
              median_H4 = round(median(PP.H4, na.rm = TRUE), 3)),
          by = method])
