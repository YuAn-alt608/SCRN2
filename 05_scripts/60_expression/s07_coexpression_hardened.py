#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s07_coexpression_hardened.py
============================
Purpose: three hardening tests of the co-expression result in s02b, each answering one objection.

Objection 1 (null distribution): s02b used ALB/INS/GFAP/MYH2 as negative controls. Those genes
        are barely expressed in thyroid (mean log2CPM 0.7-2.0), and low variance forces low
        correlation: that is a statistical necessity, not specificity. The correct null is a set
        of random genes matched on **expression level**.
Objection 2 (global axis): s02b controlled PC1-3 only. If a residual global axis remains, module
        correlations will dissolve. This script controls PC1-10 instead and reports the variance
        each PC explains.
Objection 3 (selection bias): s02b fixed the panel first and looked at correlations afterwards.
        This script does **unbiased discovery** instead: rank the whole transcriptome first,
        then see where the panel genes land.

Method:
  - partial Spearman: residualise gene ranks and covariate ranks, then correlate (as in s02b)
  - covariates: log10 library size, median expression over all genes, PC1..K of the top 2000
    variable genes
  - expression-matched null: bin genes by mean log2CPM (width 0.5), take every gene in the bin
    that contains SCRN2, giving the |rho| distribution attainable by equally expressed genes
  - module permutation: replace each module gene with a random gene from the same expression
    bin, rebuild the module score, repeat 2000 times for an empirical one-sided null

Output: results/s07_panel_pc10.csv, s07_modules_pc10.csv,
      s07_top_correlates_pc10.csv, s07_expression_matched_null.csv,
      s07_coexpression_hardened.txt
"""
import os, gzip
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1/DATA17_repro_checks"
RAW = os.path.join(BASE, "raw")
OUT = os.path.join(BASE, "results")
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(20261112)

PANEL = {
    "NKX2-1": "TF", "PAX8": "TF", "FOXE1": "TF",
    "TG": "EFF", "TPO": "EFF", "SLC5A5": "EFF", "DIO1": "EFF", "DIO2": "EFF", "TSHR": "EFF",
    "VAMP2": "EXO", "STXBP1": "EXO", "EXOC4": "EXO", "RAB27A": "EXO", "SYT1": "EXO",
    "ALB": "NEG", "INS": "NEG", "MYH2": "NEG", "GFAP": "NEG", "PTPRC": "NEG",
}
MODULES = {
    "TF": ["NKX2-1", "PAX8", "FOXE1"],
    "EFF": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR"],
    "EXO": ["VAMP2", "STXBP1", "EXOC4", "RAB27A", "SYT1"],
    "NEG": ["ALB", "INS", "MYH2", "GFAP", "PTPRC"],
}
TARGET = "SCRN2"


def load_gct(path):
    with gzip.open(path, "rt") as f:
        f.readline()
        ng, ns = map(int, f.readline().rstrip("\n").split("\t"))
        samples = f.readline().rstrip("\n").split("\t")[2:]
        rows, syms = [], []
        for line in f:
            p = line.rstrip("\n").split("\t")
            syms.append(p[1])
            rows.append(np.array(p[2:], dtype=float))
    return pd.DataFrame(np.vstack(rows), index=syms, columns=samples), ng, ns


def logcpm(mat):
    return np.log2(mat.div(mat.sum(axis=0), axis=1) * 1e6 + 1)


def build_cov(expr, npc):
    lib = expr.sum(axis=0)
    med = expr.median(axis=0)
    v = expr.var(axis=1).sort_values(ascending=False).head(2000).index
    X = expr.loc[v].values.T
    X = (X - X.mean(0)) / (X.std(0) + 1e-12)
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    var_expl = (S ** 2) / (S ** 2).sum()
    C = pd.DataFrame({"lib": np.log10(lib), "med": med}, index=expr.columns)
    for i in range(npc):
        C[f"PC{i+1}"] = U[:, i] * S[i]
    return C, var_expl


def residualise(values_by_gene, cov):
    """Residualise an (n_gene x n_sample) matrix on the covariates, in rank space."""
    n = cov.shape[0]
    C = cov.rank()
    Cm = np.column_stack([np.ones(n)] + [C[c].values for c in C.columns])
    H = Cm @ np.linalg.pinv(Cm)
    R = values_by_gene.rank(axis=1).values
    return R - R @ H.T, H


def partial_spearman_matrix(expr, target, cov):
    Rres, H = residualise(expr, cov)
    y = expr.loc[target].values
    y_rank = pd.Series(y).rank().values
    y_res = y_rank - H @ y_rank
    num = Rres @ y_res
    den = np.sqrt((Rres ** 2).sum(1) * (y_res ** 2).sum())
    r = num / np.where(den == 0, np.nan, den)
    dof = expr.shape[1] - 2 - (H.shape[1] - 1)
    t = r * np.sqrt(dof / np.clip(1 - r ** 2, 1e-12, None))
    p = 2 * stats.t.sf(np.abs(t), dof)
    return pd.DataFrame({"gene": expr.index, "rho": r, "p": p}), dof


def resid_vector(vec, cov):
    n = len(vec)
    C = cov.rank()
    Cm = np.column_stack([np.ones(n)] + [C[c].values for c in C.columns])
    H = Cm @ np.linalg.pinv(Cm)
    r = pd.Series(vec).rank().values
    return r - H @ r, H


def corr_resid(x_res, y_res):
    return float(x_res @ y_res / np.sqrt((x_res ** 2).sum() * (y_res ** 2).sum()))


log = ["=" * 80, "S07  co-expression hardening (PC1-10 + expression-matched null + unbiased discovery)", "=" * 80]

mat, ng, ns = load_gct(os.path.join(RAW, "gene_reads_v10_thyroid.gct.gz"))
expr = logcpm(mat)
expr = expr[~expr.index.duplicated()]
n = expr.shape[1]
scrn2_expr = expr.loc[TARGET]
mean_expr = expr.mean(axis=1)

cov3, ve3 = build_cov(expr, 3)
cov10, ve10 = build_cov(expr, 10)
log.append(f"\nthyroid: {expr.shape[0]} genes x {n} samples")
log.append(f"SCRN2 expression: mean log2(CPM+1)={scrn2_expr.mean():.2f}, "
           f"sd={scrn2_expr.std():.2f}, detection rate (log2CPM>1)={100*(scrn2_expr>1).mean():.1f}%")
log.append(f"variance explained by PC1-3 {100*ve3[:3].sum():.1f}%; by PC1-10 {100*ve10[:10].sum():.1f}%")

part3, dof = partial_spearman_matrix(expr, TARGET, cov3)
part10, _ = partial_spearman_matrix(expr, TARGET, cov10)
part3 = part3.set_index("gene")
part10 = part10.set_index("gene")

# ---------- 1. panel genes: PC1-3 vs PC1-10 ----------
rows = []
for g, mod in PANEL.items():
    if g not in expr.index:
        continue
    r3 = float(part3.loc[g, "rho"]); r10 = float(part10.loc[g, "rho"])
    rows.append(dict(gene=g, module=mod, mean_log2cpm=float(mean_expr[g]),
                     rho_pc3=r3, p_pc3=float(part3.loc[g, "p"]),
                     rho_pc10=r10, p_pc10=float(part10.loc[g, "p"])))
panel = pd.DataFrame(rows).sort_values("rho_pc10", ascending=False)
panel.to_csv(os.path.join(OUT, "s07_panel_pc10.csv"), index=False, encoding="utf-8-sig")

log.append("\n[1] panel genes: controlling PC1-3 (the original report) vs controlling PC1-10")
log.append(f"{'gene':<9}{'module':<6}{'meanCPM':>9}{'rho_PC3':>10}{'rho_PC10':>10}{'delta':>8}")
for _, r in panel.sort_values("rho_pc3", ascending=False).iterrows():
    log.append(f"{r.gene:<9}{r.module:<6}{r.mean_log2cpm:>9.2f}"
               f"{r.rho_pc3:>10.3f}{r.rho_pc10:>10.3f}{r.rho_pc10-r.rho_pc3:>+8.3f}")
log.append(f"  panel genes with |rho| >= 0.2: {int((panel.rho_pc3.abs()>=0.2).sum())} under PC1-3, "
           f"{int((panel.rho_pc10.abs()>=0.2).sum())} under PC1-10")

# ---------- 2. module scores: PC1-10 + expression-matched permutation ----------
zexpr = (expr.sub(expr.mean(axis=1), axis=0)).div(expr.std(axis=1) + 1e-12, axis=0)
y_res10, _ = resid_vector(scrn2_expr.values, cov10)
y_res3, _ = resid_vector(scrn2_expr.values, cov3)

# expression bins, used for the permutation
bins = (mean_expr // 0.25).astype(int)
bin_members = {b: np.where(bins.values == b)[0] for b in bins.unique()}

mod_rows = []
for name, genes in MODULES.items():
    genes = [g for g in genes if g in expr.index]
    score = zexpr.loc[genes].mean(axis=0).values
    s3, _ = resid_vector(score, cov3)
    s10, _ = resid_vector(score, cov10)
    r3 = corr_resid(s3, y_res3); r10 = corr_resid(s10, y_res10)
    # permutation: each gene replaced by a random gene from the same expression bin
    null = np.empty(2000)
    gene_idx = [expr.index.get_loc(g) for g in genes]
    for b in range(2000):
        picks = [RNG.choice(bin_members[bins.iloc[i]]) for i in gene_idx]
        sc = zexpr.values[picks].mean(axis=0)
        sc_res, _ = resid_vector(sc, cov10)
        null[b] = corr_resid(sc_res, y_res10)
    p_perm = (np.sum(null >= r10) + 1) / (2000 + 1)
    null_sd = null.std()
    mod_rows.append(dict(module=name, n_gene=len(genes), rho_pc3=r3, rho_pc10=r10,
                         perm_null_mean=null.mean(), perm_null_sd=null_sd,
                         perm_null_p95=np.quantile(null, 0.95),
                         z_vs_matched_null=(r10 - null.mean()) / null_sd,
                         p_perm_onesided=p_perm))
mods = pd.DataFrame(mod_rows)
mods.to_csv(os.path.join(OUT, "s07_modules_pc10.csv"), index=False, encoding="utf-8-sig")

log.append("\n[2] module scores: PC1-10 with an expression-matched random-module null (2000 permutations)")
log.append(f"{'module':<6}{'n':>4}{'rho_PC3':>10}{'rho_PC10':>10}{'null mean':>10}{'null SD':>9}"
           f"{'null 95%':>10}{'z':>8}{'p_perm':>10}")
for _, r in mods.iterrows():
    log.append(f"{r.module:<6}{r.n_gene:>4}{r.rho_pc3:>10.3f}{r.rho_pc10:>10.3f}"
               f"{r.perm_null_mean:>9.3f}{r.perm_null_sd:>8.3f}{r.perm_null_p95:>8.3f}"
               f"{r.z_vs_matched_null:>8.2f}{r.p_perm_onesided:>10.3g}")

# ---------- 3. unbiased discovery: whole-transcriptome ranking (PC1-10) ----------
disc = part10.reset_index().rename(columns={"rho": "rho_pc10", "p": "p_pc10"})
disc["rho_pc3"] = disc.gene.map(part3["rho"])
disc["mean_log2cpm"] = disc.gene.map(mean_expr)
disc = disc[disc.gene != TARGET]
disc["rank_pc10"] = disc.rho_pc10.rank(ascending=False).astype(int)
disc["rank_pc3"] = disc.rho_pc3.rank(ascending=False).astype(int)
disc = disc.sort_values("rho_pc10", ascending=False)
disc.to_csv(os.path.join(OUT, "s07_top_correlates_pc10.csv"), index=False, encoding="utf-8-sig")

log.append("\n[3] unbiased discovery: the top 15 genome-wide co-expression partners of SCRN2 (PC1-10)")
log.append(f"{'rank':>5}  {'gene':<14}{'rho_PC10':>10}{'rho_PC3':>10}{'meanCPM':>9}")
for i, r in enumerate(disc.head(15).itertuples(), 1):
    log.append(f"{i:>5}  {r.gene:<14}{r.rho_pc10:>10.3f}{r.rho_pc3:>10.3f}{r.mean_log2cpm:>9.2f}")

log.append("\n[4] where the panel genes sit in the unbiased ranking (out of %d genes)" % len(disc))
for g in ["NKX2-1", "PAX8", "FOXE1", "TG", "TPO", "SLC5A5", "RAB27A", "ALB"]:
    s = disc[disc.gene == g]
    if len(s):
        r = s.iloc[0]
        log.append(f"  {g:<9} PC1-10 rank {int(r.rank_pc10):>6}/{len(disc)}  "
                   f"(PC1-3 rank {int(r.rank_pc3):>6})  rho_PC10={r.rho_pc10:+.3f}")
    else:
        log.append(f"  {g:<9} not detected")

# ---------- 4. expression-matched null (single gene) ----------
sq = mean_expr[TARGET]
matched = mean_expr[(mean_expr - sq).abs() <= 0.5].index
matched = [g for g in matched if g != TARGET and g in part10.index]
mnull = part10.loc[matched, "rho"].dropna().values
pct = (np.sum(np.abs(mnull) >= 0) if False else None)

null_rows = []
for g in panel.gene:
    r10 = float(part10.loc[g, "rho"])
    same_side = np.sum(mnull >= r10) if r10 >= 0 else np.sum(mnull <= r10)
    p_emp = (same_side + 1) / (len(mnull) + 1)
    null_rows.append(dict(gene=g, rho_pc10=r10, n_matched=len(mnull),
                          null_abs_median=float(np.median(np.abs(mnull))),
                          null_abs_p95=float(np.quantile(np.abs(mnull), 0.95)),
                          p_emp_matched=p_emp))
nulldf = pd.DataFrame(null_rows)
nulldf.to_csv(os.path.join(OUT, "s07_expression_matched_null.csv"),
              index=False, encoding="utf-8-sig")

log.append(f"\n[5] expression-matched null: SCRN2 mean log2CPM={sq:.2f}, "
           f"{len(mnull)} genes within +/-0.5")
log.append(f"  null |rho|: median {np.median(np.abs(mnull)):.3f}, "
           f"95th percentile {np.quantile(np.abs(mnull),0.95):.3f}, "
           f"maximum {np.max(np.abs(mnull)):.3f}")
log.append(f"  (the negative control ALB in the original report reached only rho {float(part10.loc['ALB','rho']):+.3f}: "
           f"a low correlation is inevitable for a lowly expressed gene and is no noise floor)")
log.append(f"\n{'gene':<9}{'rho_PC10':>10}{'n matched':>10}{'null |rho| med':>15}"
           f"{'null |rho| 95%':>16}{'emp p':>9}")
for _, r in nulldf.sort_values("rho_pc10", ascending=False).iterrows():
    log.append(f"{r.gene:<9}{r.rho_pc10:>10.3f}{int(r.n_matched):>8}"
               f"{r.null_abs_median:>12.3f}{r.null_abs_p95:>12.3f}{r.p_emp_matched:>10.3g}")

log.append("\n[6] reading rules, fixed in advance so they cannot be rationalised afterwards")
log.append("  A. a module counts as carrying signal only if its rho exceeds the 95th percentile of the")
log.append("     expression-matched random-module null;")
log.append("  B. a panel gene deserves separate discussion only if it enters the top 1% of the unbiased")
log.append("     ranking (<579/57852);")
log.append("  C. if |rho| under PC1-10 falls below 0.2 with an empirical p>0.05, that gene's correlation")
log.append("     is not usable.")

txt = "\n".join(log)
open(os.path.join(OUT, "s07_coexpression_hardened.txt"), "w", encoding="utf-8").write(txt)
print(txt)
