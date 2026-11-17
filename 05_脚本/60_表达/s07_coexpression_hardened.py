#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s07_coexpression_hardened.py
============================
目的：对 s02b 的共表达结论做三项加固检验，回答三个具体质疑。

质疑 1（零分布）：s02b 用 ALB/INS/GFAP/MYH2 作阴性对照。这些基因在甲状腺
        近乎不表达（mean log2CPM 0.7–2.0），低方差必然低相关——这是统计必然而非
        特异性。正确的零分布应是**表达量匹配**的随机基因。
质疑 2（全局轴）：s02b 只控 PC1–3。若残余全局轴仍在，模块相关会随之消散。
        本脚本改控 PC1–10，并报告各 PC 解释的方差。
质疑 3（选择偏倚）：s02b 先定面板再看相关。本脚本改为**无偏发现**：
        先在全转录组排序，再看面板基因落在哪。

方法：
  - 偏 Spearman：对基因秩与协变量秩同时取残差后求相关（与 s02b 同法）
  - 协变量：log10 文库大小、全基因中位表达、top2000 变异基因的 PC1..K
  - 表达匹配零分布：把基因按 mean log2CPM 分箱（宽 0.5），在 SCRN2 所在箱内
    取所有基因，得到"同等表达量基因能达到的 |rho| 分布"
  - 模块置换：每个模块基因替换为同表达箱内的随机基因，重建模块分数，
    重复 2000 次，得到模块相关的经验零分布（单侧）

输出：results/s07_panel_pc10.csv, s07_modules_pc10.csv,
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
    """对 (n_gene × n_sample) 矩阵按协变量求残差（秩空间）。"""
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


log = ["=" * 80, "S07  共表达加固检验（PC1–10 + 表达匹配零分布 + 无偏发现）", "=" * 80]

mat, ng, ns = load_gct(os.path.join(RAW, "gene_reads_v10_thyroid.gct.gz"))
expr = logcpm(mat)
expr = expr[~expr.index.duplicated()]
n = expr.shape[1]
scrn2_expr = expr.loc[TARGET]
mean_expr = expr.mean(axis=1)

cov3, ve3 = build_cov(expr, 3)
cov10, ve10 = build_cov(expr, 10)
log.append(f"\n甲状腺：{expr.shape[0]} 基因 × {n} 样本")
log.append(f"SCRN2 表达：mean log2(CPM+1)={scrn2_expr.mean():.2f}，"
           f"sd={scrn2_expr.std():.2f}，检出率(log2CPM>1)={100*(scrn2_expr>1).mean():.1f}%")
log.append(f"PC1–3 解释方差 {100*ve3[:3].sum():.1f}%；PC1–10 解释方差 {100*ve10[:10].sum():.1f}%")

part3, dof = partial_spearman_matrix(expr, TARGET, cov3)
part10, _ = partial_spearman_matrix(expr, TARGET, cov10)
part3 = part3.set_index("gene")
part10 = part10.set_index("gene")

# ---------- 1. 面板基因：PC1-3 vs PC1-10 ----------
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

log.append("\n【1】面板基因：控 PC1–3（原报告）vs 控 PC1–10")
log.append(f"{'基因':<9}{'模块':<6}{'meanCPM':>9}{'rho_PC3':>10}{'rho_PC10':>10}{'Δ':>8}")
for _, r in panel.sort_values("rho_pc3", ascending=False).iterrows():
    log.append(f"{r.gene:<9}{r.module:<6}{r.mean_log2cpm:>9.2f}"
               f"{r.rho_pc3:>10.3f}{r.rho_pc10:>10.3f}{r.rho_pc10-r.rho_pc3:>+8.3f}")
log.append(f"  |rho| ≥ 0.2 的面板基因：PC1–3 下 {int((panel.rho_pc3.abs()>=0.2).sum())} 个，"
           f"PC1–10 下 {int((panel.rho_pc10.abs()>=0.2).sum())} 个")

# ---------- 2. 模块分数：PC1–10 + 表达匹配置换 ----------
zexpr = (expr.sub(expr.mean(axis=1), axis=0)).div(expr.std(axis=1) + 1e-12, axis=0)
y_res10, _ = resid_vector(scrn2_expr.values, cov10)
y_res3, _ = resid_vector(scrn2_expr.values, cov3)

# 表达分箱，供置换用
bins = (mean_expr // 0.25).astype(int)
bin_members = {b: np.where(bins.values == b)[0] for b in bins.unique()}

mod_rows = []
for name, genes in MODULES.items():
    genes = [g for g in genes if g in expr.index]
    score = zexpr.loc[genes].mean(axis=0).values
    s3, _ = resid_vector(score, cov3)
    s10, _ = resid_vector(score, cov10)
    r3 = corr_resid(s3, y_res3); r10 = corr_resid(s10, y_res10)
    # 置换：每个基因换成同表达箱内随机基因
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

log.append("\n【2】模块分数：PC1–10 + 表达匹配随机模块零分布（2000 次置换）")
log.append(f"{'模块':<6}{'n':>4}{'rho_PC3':>10}{'rho_PC10':>10}{'零均值':>9}{'零SD':>8}"
           f"{'零95%':>8}{'z':>8}{'p_perm':>10}")
for _, r in mods.iterrows():
    log.append(f"{r.module:<6}{r.n_gene:>4}{r.rho_pc3:>10.3f}{r.rho_pc10:>10.3f}"
               f"{r.perm_null_mean:>9.3f}{r.perm_null_sd:>8.3f}{r.perm_null_p95:>8.3f}"
               f"{r.z_vs_matched_null:>8.2f}{r.p_perm_onesided:>10.3g}")

# ---------- 3. 无偏发现：全转录组排序（PC1–10） ----------
disc = part10.reset_index().rename(columns={"rho": "rho_pc10", "p": "p_pc10"})
disc["rho_pc3"] = disc.gene.map(part3["rho"])
disc["mean_log2cpm"] = disc.gene.map(mean_expr)
disc = disc[disc.gene != TARGET]
disc["rank_pc10"] = disc.rho_pc10.rank(ascending=False).astype(int)
disc["rank_pc3"] = disc.rho_pc3.rank(ascending=False).astype(int)
disc = disc.sort_values("rho_pc10", ascending=False)
disc.to_csv(os.path.join(OUT, "s07_top_correlates_pc10.csv"), index=False, encoding="utf-8-sig")

log.append("\n【3】无偏发现：SCRN2 全转录组共表达 top 15（PC1–10 口径）")
log.append(f"{'rank':>5}  {'基因':<14}{'rho_PC10':>10}{'rho_PC3':>10}{'meanCPM':>9}")
for i, r in enumerate(disc.head(15).itertuples(), 1):
    log.append(f"{i:>5}  {r.gene:<14}{r.rho_pc10:>10.3f}{r.rho_pc3:>10.3f}{r.mean_log2cpm:>9.2f}")

log.append("\n【4】面板基因在无偏排序中的位置（总基因数 %d）" % len(disc))
for g in ["NKX2-1", "PAX8", "FOXE1", "TG", "TPO", "SLC5A5", "RAB27A", "ALB"]:
    s = disc[disc.gene == g]
    if len(s):
        r = s.iloc[0]
        log.append(f"  {g:<9} PC1–10 排名 {int(r.rank_pc10):>6}/{len(disc)}  "
                   f"(PC1–3 排名 {int(r.rank_pc3):>6})  rho_PC10={r.rho_pc10:+.3f}")
    else:
        log.append(f"  {g:<9} 未检出")

# ---------- 4. 表达匹配零分布（单基因） ----------
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

log.append(f"\n【5】表达匹配零分布：SCRN2 mean log2CPM={sq:.2f}，"
           f"±0.5 内共 {len(mnull)} 个基因")
log.append(f"  该零分布 |rho| 中位数 {np.median(np.abs(mnull)):.3f}，"
           f"95% 分位 {np.quantile(np.abs(mnull),0.95):.3f}，"
           f"最大 {np.max(np.abs(mnull)):.3f}")
log.append(f"  （原报告阴性对照 ALB 的 rho 仅 {float(part10.loc['ALB','rho']):+.3f}——"
           f"低表达基因的低相关是必然，不能当噪声地板）")
log.append(f"\n{'基因':<9}{'rho_PC10':>10}{'匹配数':>8}{'零|rho|中位':>12}"
           f"{'零|rho|95%':>12}{'经验p':>10}")
for _, r in nulldf.sort_values("rho_pc10", ascending=False).iterrows():
    log.append(f"{r.gene:<9}{r.rho_pc10:>10.3f}{int(r.n_matched):>8}"
               f"{r.null_abs_median:>12.3f}{r.null_abs_p95:>12.3f}{r.p_emp_matched:>10.3g}")

log.append("\n【6】判读规则（写死，避免事后合理化）")
log.append("  A. 模块 rho 须超过表达匹配随机模块零分布 95% 分位，才算该模块有信号；")
log.append("  B. 面板基因须在无偏排序中进入前 1%（<579/57852）才值得单独讨论；")
log.append("  C. 若 PC1–10 下 |rho| 掉到 0.2 以下且经验 p>0.05，则该基因的相关不可用。")

txt = "\n".join(log)
open(os.path.join(OUT, "s07_coexpression_hardened.txt"), "w", encoding="utf-8").write(txt)
print(txt)
