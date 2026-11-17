#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s08_tsh_t4_plane.py
===================
目的：给 SCRN2 在 (TSH, T4) 二维平面上找一个参照系，回答"这个签名是不是站得住"。

背景：
  每个基因的 MR 估计 b = beta_out / beta_eqtl。eQTL 的效应等位基因方向逐基因不同，
  因此每个点只定到一个整体正负号。**在 (b_TSH, b_T4) 同时翻号下不变的量是
  比率 r = b_T4 / b_TSH**，它才是可跨基因比较的量。
  - r ≈ -1.07（TPO 的斜率）  → 与合成通路同比例（激素产出改变 → 反馈改变 TSH）
  - r ≈ 0                     → 只动 TSH，不动 T4（"TSH-only"）
  - |r| → ∞                   → 只动 T4

三个问题：
  Q1  TSH 显著的基因里，r 的分布是什么样？是聚在 -1.07，还是散在 0？
  Q2  "TSH 显著而 T4 为零"是常态还是个例？→ 决定 SCRN2 的签名是否特殊
  Q3  SCRN2 的 r=+0.45 落在分布的哪个位置？

输出：results/s08_plane_genes.csv, s08_ratio_summary.csv, s08_plane.txt
      figures/fig_tsh_t4_plane.png
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1/DATA17_repro_checks"
D15 = "D:/NewBeginning/MR1/DATA15_thyroid_eQTL/results"
OUT = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
os.makedirs(OUT, exist_ok=True)
ZA = stats.norm.ppf(0.975)

d = pd.read_csv(os.path.join(D15, "gw_thyroid_mr.csv"))
w = d.pivot_table(index="gene", columns="outcome",
                  values=["b", "se", "p", "fdr", "F", "snp"], aggfunc="first")
w.columns = [f"{a}_{b}" for a, b in w.columns]
w = w.reset_index()
w = w.dropna(subset=["b_TSH", "b_T4"]).copy()
N = len(w)

# 比率及其 delta 法 se
w["r"] = w.b_T4 / w.b_TSH
w["se_r"] = np.abs(w.r) * np.sqrt((w.se_T4 / w.b_T4) ** 2 + (w.se_TSH / w.b_TSH) ** 2)
w["r_lo"] = w.r - ZA * w.se_r
w["r_hi"] = w.r + ZA * w.se_r
# T4 效应是否为 0（95% CI 含 0）
w["t4_null"] = (w.b_T4 - ZA * w.se_T4 <= 0) & (w.b_T4 + ZA * w.se_T4 >= 0)
w["tsh_sig"] = w.fdr_TSH < 0.05
w["t4_sig"] = w.fdr_T4 < 0.05

K_TPO = -1.0747
L = ["=" * 84, "S08  (TSH, T4) 二维平面：SCRN2 的参照系", "=" * 84]
L.append(f"\n可分析基因：{N}（同时有 TSH 与 T4 估计）")
L.append(f"FDR<0.05：TSH {int(w.tsh_sig.sum())} 个，T4 {int(w.t4_sig.sum())} 个，"
         f"两者兼备 {int((w.tsh_sig & w.t4_sig).sum())} 个")

# ---------- Q1/Q2：TSH 显著基因的 r 分布 ----------
A = w[w.tsh_sig & (w.b_TSH.abs() > 0.005)].copy()   # 排除 b_TSH≈0 导致比率不稳
L.append(f"\n【Q1】TSH 显著（FDR<0.05）且 |b_TSH|>0.005 的基因：{len(A)} 个")

med = A.r.median()
iqr = (A.r.quantile(0.25), A.r.quantile(0.75))
L.append(f"  r 中位数 {med:+.2f}，IQR [{iqr[0]:+.2f}, {iqr[1]:+.2f}]，"
         f"范围 [{A.r.min():+.2f}, {A.r.max():+.2f}]")
L.append(f"  参照：TPO 的 r = {K_TPO:+.2f}")

# 分箱
def cat(row):
    if row.r_lo <= K_TPO <= row.r_hi:
        return "合成同比例"
    if row.t4_null:
        return "TSH-only（T4 CI 含 0）"
    if row.r > 0:
        return "同向（T4 随 TSH 升）"
    return "部分合成（负但不同比例）"

A["pattern"] = A.apply(cat, axis=1)
vc = A.pattern.value_counts()
L.append("\n  模式构成：")
for k, v in vc.items():
    L.append(f"    {k:<24}{v:>4}  ({100*v/len(A):.1f}%)")

# 关键：TSH-only 是不是常态
n_tsh_only = int(A.t4_null.sum())
L.append(f"\n【Q2】TSH 显著基因中，T4 的 95% CI 含 0 的有 {n_tsh_only}/{len(A)} "
         f"（{100*n_tsh_only/len(A):.1f}%）")
L.append(f"  若 'TSH-only' 占比很高 → SCRN2 的签名是这一层的常态，不具特异性")
L.append(f"  若占比很低 → SCRN2 是少数派，签名有信息量")

# ---------- Q3：SCRN2 的位置 ----------
s = w[w.gene == "SCRN2"].iloc[0]
pct = 100 * (A.r < s.r).mean()
# 双侧尾部概率
tail = 2 * min((A.r >= s.r).mean(), (A.r <= s.r).mean())
L.append(f"\n【Q3】SCRN2：b_TSH={s.b_TSH:+.4f}, b_T4={s.b_T4:+.4f}, "
         f"r={s.r:+.3f}（95% CI [{s.r_lo:+.3f}, {s.r_hi:+.3f}]）")
L.append(f"  在 TSHSig 集合中的百分位：{pct:.1f}%（共 {len(A)} 个基因）")
L.append(f"  双侧经验 p（相对该分布）：{tail:.3f}")
L.append(f"  与 TPO 比例的一致性：{'一致' if s.r_lo <= K_TPO <= s.r_hi else '不一致'}")

# ---------- 检验：单一比例轴能否解释整个分布 ----------
L.append("\n【4】单一比例轴假设的检验")
L.append("  若所有基因都经同一通路起作用，则 r 应集中于 TPO 的斜率。")
inside = int(((A.r_lo <= K_TPO) & (K_TPO <= A.r_hi)).sum())
p_binom = stats.binomtest(inside, len(A), 0.95, alternative="less").pvalue
L.append(f"  r 的 95% CI 覆盖 {K_TPO:+.2f} 的基因：{inside}/{len(A)} = {100*inside/len(A):.1f}%")
L.append(f"  若真实比例为 95%（单一通路），观测到这么少的概率 p = {p_binom:.2e}")
L.append(f"  → 单一比例轴{'不成立' if p_binom < 1e-3 else '不能被排除'}")

# 加权相关
wx = 1 / (w.se_TSH ** 2); wy = 1 / (w.se_T4 ** 2)
sub = w[(w.b_TSH.abs() > 0.005)]
rho_s = stats.spearmanr(sub.b_TSH, sub.b_T4)
L.append(f"\n  全部 {len(sub)} 个基因的 b_TSH 与 b_T4 Spearman rho = {rho_s[0]:+.3f}（p={rho_s[1]:.2e}）")
L.append(f"  负相关 → 两轴共享反馈结构；但强度弱 → 大量基因只有一端变化")

w.to_csv(os.path.join(OUT, "s08_plane_genes.csv"), index=False, encoding="utf-8-sig")
A.to_csv(os.path.join(OUT, "s08_ratio_summary.csv"), index=False, encoding="utf-8-sig")

L.append("\n【5】TSH 显著基因中 T4 效应的分布（按 |z_T4| 排序，前 12）")
top = A.reindex(A.b_T4.abs().sort_values(ascending=False).index).head(12)
L.append(f"{'基因':<10}{'b_TSH':>9}{'b_T4':>9}{'r':>8}{'r_CI':>20}{'模式'}")
for _, x in top.iterrows():
    L.append(f"{x.gene:<10}{x.b_TSH:>9.4f}{x.b_T4:>9.4f}{x.r:>8.2f}"
             f"   [{x.r_lo:>6.2f},{x.r_hi:>6.2f}]  {x.pattern}")

# 已知甲状腺基因的定位
L.append("\n【6】已知甲状腺功能基因在平面上的位置")
KNOWN = ["TPO", "TG", "SLC5A5", "TSHR", "DUOX2", "DUOXA2", "IYD", "SLC26A4",
         "PAX8", "NKX2-1", "FOXE1", "DIO1", "DIO2", "SLC16A10", "THRA", "THRB",
         "TSHB", "SCRN2"]
L.append(f"{'基因':<10}{'b_TSH':>9}{'p_TSH':>10}{'b_T4':>9}{'p_T4':>10}{'r':>8}  判定")
for g in KNOWN:
    x = w[w.gene == g]
    if not len(x):
        continue
    x = x.iloc[0]
    verdict = ("TSH sig" if x.fdr_TSH < 0.05 else "") + (" / T4 sig" if x.fdr_T4 < 0.05 else "")
    verdict = verdict.strip(" /") or "都不是"
    L.append(f"{g:<10}{x.b_TSH:>9.4f}{x.p_TSH:>10.2e}{x.b_T4:>9.4f}{x.p_T4:>10.2e}"
             f"{x.r:>8.2f}  {verdict}")

txt = "\n".join(L)
open(os.path.join(OUT, "s08_plane.txt"), "w", encoding="utf-8").write(txt)
print(txt)
