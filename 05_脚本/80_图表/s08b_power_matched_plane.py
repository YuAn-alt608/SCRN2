#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s08b_power_matched_plane.py
===========================
"S08 说 68% 的 TSH 显著基因 T4 为零"——但 T4 GWAS（n≈191k）只有 TSH（n≈483k）的四成，
T4 的 se 大约是 TSH 的 1.6 倍。所以 68% 这个数可能纯粹是功效不足的产物。

本脚本把功效算进去，只在**对比例性效应有足够功效**的基因里比较。

步骤：
  1) 对每个 TSH 显著基因，按 TPO 斜率算出"若按比例起效"的 T4 预测值 pred=k*b_TSH
  2) 用该基因实际的 se_T4 计算检出 pred 的功效（α=0.05 双尾）
  3) 分层：有功效（≥80%）vs 无功效（<80%）
  4) 在有功效层里，看实际偏离比例性的比例，以及 SCRN2 的位置
  5) 期望检出数 vs 实际检出数 → 判断"T4 为零"是生物学还是功效

输出：results/s08b_power.csv, s08b_power.txt
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1/DATA17_repro_checks"
OUT = os.path.join(BASE, "results")
ZA = stats.norm.ppf(0.975)
K = -1.0747          # TPO 斜率
SE_K = 0.1884

w = pd.read_csv(os.path.join(OUT, "s08_plane_genes.csv"))
A = w[(w.fdr_TSH < 0.05) & (w.b_TSH.abs() > 0.005)].copy()

# 预测的 T4 效应（含斜率不确定度，用在该基因上）
A["pred_T4"] = K * A.b_TSH
A["se_pred"] = np.sqrt((A.b_TSH * SE_K) ** 2 + (K * A.se_TSH) ** 2)

# 检出 pred 的功效（相对零假设 b_T4=0）
ncp = A.pred_T4.abs() / A.se_T4
A["power"] = stats.norm.sf(ZA - ncp) + stats.norm.cdf(-ZA - ncp)

# 偏离比例性：pred 是否落在 b_T4 的 95% CI 外
A["ci_lo"] = A.b_T4 - ZA * A.se_T4
A["ci_hi"] = A.b_T4 + ZA * A.se_T4
A["departs"] = (A.pred_T4 < A.ci_lo) | (A.pred_T4 > A.ci_hi)

# T4 显著
A["t4_sig"] = A.fdr_T4 < 0.05

L = ["=" * 84, "S08b  功效匹配后的平面分析", "=" * 84]
L.append(f"\nTSH 显著基因 {len(A)} 个；T4 GWAS 与 TSH GWAS 样本量比约 0.40")
L.append(f"按比例模型预测的 T4 效应：|pred| 中位数 {A.pred_T4.abs().median():.4f}")
L.append(f"实际 se_T4 中位数 {A.se_T4.median():.4f}")
L.append(f"→ 中位基因检出比例性效应的功效仅 {100*A.power.median():.0f}%")

for thr in (0.80, 0.50):
    P = A[A.power >= thr]
    L.append(f"\n【功效 ≥{int(thr*100)}% 的亚组】n={len(P)}")
    if len(P) == 0:
        continue
    L.append(f"  T4 达 FDR<0.05：{int(P.t4_sig.sum())} 个（{100*P.t4_sig.mean():.1f}%）")
    L.append(f"  预期检出（若全部按比例起效）：{P.power.sum():.1f} 个")
    L.append(f"  → 实际/预期 = {P.t4_sig.sum()/max(P.power.sum(),1e-9):.2f}")
    L.append(f"  偏离比例性（pred 落在 T4 的 95% CI 外）：{int(P.departs.sum())}/{len(P)} "
             f"= {100*P.departs.mean():.1f}%")
    if (P.gene == "SCRN2").any():
        s = P[P.gene == "SCRN2"].iloc[0]
        pct = 100 * (P.departs & (P.b_T4 > P.pred_T4)).mean()
        L.append(f"  SCRN2 在此亚组：power={s.power:.2f}，pred={s.pred_T4:+.4f}，"
                 f"实测={s.b_T4:+.4f}，偏离={bool(s.departs)}")

# ---------- SCRN2 的"正向偏离"在同层里有多常见 ----------
P = A[A.power >= 0.80].copy()
if (P.gene == "SCRN2").any():
    s = P[P.gene == "SCRN2"].iloc[0]
    higher = P[(P.departs) & (P.b_T4 > P.pred_T4)]
    L.append(f"\n【SCRN2 的特殊性检验】")
    L.append(f"  SCRN2 属于'有功效但 T4 高于比例预测'的方向：pred {s.pred_T4:+.4f} → 实测 {s.b_T4:+.4f}")
    L.append(f"  同一亚组中同方向的基因 {len(higher)}/{len(P)}")
    if len(higher):
        L.append(f"  其中 b_T4 高于 SCRN2 的：{int((higher.b_T4 > s.b_T4).sum())} 个")
        L.append("  例：" + ", ".join(higher.sort_values('b_T4', ascending=False).gene.head(10)))
    p_two = 2 * min((P.b_T4 >= s.b_T4).mean(), (P.b_T4 <= s.b_T4).mean())
    L.append(f"  SCRN2 的 b_T4 在该亚组中的双侧经验 p = {p_two:.3f}")

A.to_csv(os.path.join(OUT, "s08b_power.csv"), index=False, encoding="utf-8-sig")

L.append("\n【结论判据】")
L.append("  若 实际/预期 ≈ 1  → 'T4 为零'主要是功效不足，SCRN2 的签名不具特异性")
L.append("  若 实际/预期 << 1 → 多数基因确实不动 T4，'TSH-only' 是常态")
P80 = A[A.power >= 0.80]
if len(P80):
    ratio = P80.t4_sig.sum() / max(P80.power.sum(), 1e-9)
    L.append(f"\n  本数据：功效≥80% 亚组实际/预期 = {ratio:.2f}")

txt = "\n".join(L)
open(os.path.join(OUT, "s08b_power.txt"), "w", encoding="utf-8").write(txt)
print(txt)
