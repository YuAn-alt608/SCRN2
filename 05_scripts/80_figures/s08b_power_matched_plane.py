#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s08b_power_matched_plane.py
===========================
"S08 says 68% of TSH-significant genes have a null T4 effect" -- but the T4 GWAS
(n~191k) has only about four tenths the sample size of TSH (n~483k), so the T4 se
is roughly 1.6 times larger. That 68% could therefore be an artefact of
insufficient power.

This script accounts for power and only compares genes that are **adequately
powered to detect a proportional effect**.

Steps:
  1) for each TSH-significant gene, compute the T4 value a proportional effect
     would predict from the TPO slope: pred = k * b_TSH
  2) use that gene's actual se_T4 to compute the power to detect pred (two-sided
     alpha = 0.05)
  3) stratify into adequately powered (>=80%) and underpowered (<80%)
  4) within the powered stratum, measure how often the data depart from
     proportionality, and where SCRN2 sits
  5) compare expected with observed detections -> decide whether "T4 is null" is
     biology or power

Output: 04_data/08_plane/s08b_power.csv and s08b_power.txt
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "04_data", "08_plane")
ZA = stats.norm.ppf(0.975)
K = -1.0747          # the TPO slope
SE_K = 0.1884

w = pd.read_csv(os.path.join(OUT, "s08_plane_genes.csv"))
A = w[(w.fdr_TSH < 0.05) & (w.b_TSH.abs() > 0.005)].copy()

# the predicted T4 effect (propagating the slope uncertainty for this gene)
A["pred_T4"] = K * A.b_TSH
A["se_pred"] = np.sqrt((A.b_TSH * SE_K) ** 2 + (K * A.se_TSH) ** 2)

# power to detect pred against the null b_T4 = 0
ncp = A.pred_T4.abs() / A.se_T4
A["power"] = stats.norm.sf(ZA - ncp) + stats.norm.cdf(-ZA - ncp)

# departure from proportionality: does pred fall outside the 95% CI of b_T4
A["ci_lo"] = A.b_T4 - ZA * A.se_T4
A["ci_hi"] = A.b_T4 + ZA * A.se_T4
A["departs"] = (A.pred_T4 < A.ci_lo) | (A.pred_T4 > A.ci_hi)

# T4 significant
A["t4_sig"] = A.fdr_T4 < 0.05

L = ["=" * 84, "S08b  the plane after power matching", "=" * 84]
L.append(f"\nTSH-significant genes: {len(A)}; T4 GWAS to TSH GWAS sample size ratio about 0.40")
L.append(f"T4 effect predicted by the proportional model: median |pred| {A.pred_T4.abs().median():.4f}")
L.append(f"Observed median se_T4: {A.se_T4.median():.4f}")
L.append(f"-> power to detect a proportional effect in the median gene is only {100*A.power.median():.0f}%")

for thr in (0.80, 0.50):
    P = A[A.power >= thr]
    L.append(f"\n[Subgroup with power >={int(thr*100)}%] n={len(P)}")
    if len(P) == 0:
        continue
    L.append(f"  T4 reaching FDR<0.05: {int(P.t4_sig.sum())} ({100*P.t4_sig.mean():.1f}%)")
    L.append(f"  Expected detections if every gene acted proportionally: {P.power.sum():.1f}")
    L.append(f"  -> observed/expected = {P.t4_sig.sum()/max(P.power.sum(),1e-9):.2f}")
    L.append(f"  Departure from proportionality (pred outside the T4 95% CI): {int(P.departs.sum())}/{len(P)} "
             f"= {100*P.departs.mean():.1f}%")
    if (P.gene == "SCRN2").any():
        s = P[P.gene == "SCRN2"].iloc[0]
        pct = 100 * (P.departs & (P.b_T4 > P.pred_T4)).mean()
        L.append(f"  SCRN2 in this subgroup: power={s.power:.2f}, pred={s.pred_T4:+.4f}, "
                 f"observed={s.b_T4:+.4f}, departs={bool(s.departs)}")

# ---------- how common an upward departure like the SCRN2 one is in the same stratum ----------
P = A[A.power >= 0.80].copy()
if (P.gene == "SCRN2").any():
    s = P[P.gene == "SCRN2"].iloc[0]
    higher = P[(P.departs) & (P.b_T4 > P.pred_T4)]
    L.append(f"\n[Test of what makes SCRN2 special]")
    L.append(f"  SCRN2 falls in the direction 'adequately powered but T4 above the proportional prediction': pred {s.pred_T4:+.4f} -> observed {s.b_T4:+.4f}")
    L.append(f"  Genes in the same subgroup pointing the same way: {len(higher)}/{len(P)}")
    if len(higher):
        L.append(f"  Of those, genes with b_T4 higher than SCRN2: {int((higher.b_T4 > s.b_T4).sum())}")
        L.append("  Examples: " + ", ".join(higher.sort_values('b_T4', ascending=False).gene.head(10)))
    p_two = 2 * min((P.b_T4 >= s.b_T4).mean(), (P.b_T4 <= s.b_T4).mean())
    L.append(f"  Two-sided empirical p for the SCRN2 b_T4 within that subgroup = {p_two:.3f}")

A.to_csv(os.path.join(OUT, "s08b_power.csv"), index=False, encoding="utf-8-sig")

L.append("\n[Decision rule]")
L.append("  If observed/expected ~ 1  -> 'T4 is zero' mostly reflects insufficient power and the SCRN2 signature is not specific")
L.append("  If observed/expected << 1 -> most genes genuinely do not move T4, and 'TSH-only' is the norm")
P80 = A[A.power >= 0.80]
if len(P80):
    ratio = P80.t4_sig.sum() / max(P80.power.sum(), 1e-9)
    L.append(f"\n  In these data: power>=80% subgroup observed/expected = {ratio:.2f}")

txt = "\n".join(L)
open(os.path.join(OUT, "s08b_power.txt"), "w", encoding="utf-8").write(txt + "\n")
print(txt)
