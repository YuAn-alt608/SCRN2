#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s08_tsh_t4_plane.py
===================
Aim: give SCRN2 a reference frame on the (TSH, T4) plane and answer "does this
signature hold up".

Background:
  The MR estimate for each gene is b = beta_out / beta_eqtl. The effect-allele
  direction of the eQTL differs from gene to gene, so each point is only defined up
  to an overall sign. **The quantity that is invariant under flipping both b_TSH
  and b_T4 is the ratio r = b_T4 / b_TSH**, and that is what can be compared
  across genes.
  - r ~ -1.07 (the TPO slope)  -> same ratio as the synthetic pathway (hormone output changes -> feedback changes TSH)
  - r ~ 0                      -> TSH moves, T4 does not ("TSH-only")
  - |r| -> infinity            -> only T4 moves

Three questions:
  Q1  among TSH-significant genes, what does the r distribution look like? Does it cluster at -1.07, or scatter around 0?
  Q2  is "TSH significant, T4 null" the norm or the exception? -> decides whether the SCRN2 signature is special
  Q3  where does the SCRN2 r=+0.45 sit in that distribution?

Output: results/s08_plane_genes.csv, s08_ratio_summary.csv, s08_plane.txt
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

# the ratio and its delta-method se
w["r"] = w.b_T4 / w.b_TSH
w["se_r"] = np.abs(w.r) * np.sqrt((w.se_T4 / w.b_T4) ** 2 + (w.se_TSH / w.b_TSH) ** 2)
w["r_lo"] = w.r - ZA * w.se_r
w["r_hi"] = w.r + ZA * w.se_r
# is the T4 effect null (95% CI includes 0)
w["t4_null"] = (w.b_T4 - ZA * w.se_T4 <= 0) & (w.b_T4 + ZA * w.se_T4 >= 0)
w["tsh_sig"] = w.fdr_TSH < 0.05
w["t4_sig"] = w.fdr_T4 < 0.05

K_TPO = -1.0747
L = ["=" * 84, "S08  the (TSH, T4) plane: a reference frame for SCRN2", "=" * 84]
L.append(f"\nAnalysable genes: {N} (TSH and T4 estimates both available)")
L.append(f"FDR<0.05: TSH {int(w.tsh_sig.sum())} genes, T4 {int(w.t4_sig.sum())} genes, "
         f"both {int((w.tsh_sig & w.t4_sig).sum())} genes")

# ---------- Q1/Q2: the r distribution among TSH-significant genes ----------
A = w[w.tsh_sig & (w.b_TSH.abs() > 0.005)].copy()   # exclude b_TSH~0, where the ratio is unstable
L.append(f"\n[Q1] Genes with a significant TSH effect (FDR<0.05) and |b_TSH|>0.005: {len(A)}")

med = A.r.median()
iqr = (A.r.quantile(0.25), A.r.quantile(0.75))
L.append(f"  median r {med:+.2f}, IQR [{iqr[0]:+.2f}, {iqr[1]:+.2f}], "
         f"range [{A.r.min():+.2f}, {A.r.max():+.2f}]")
L.append(f"  reference: TPO has r = {K_TPO:+.2f}")

# binning
def cat(row):
    if row.r_lo <= K_TPO <= row.r_hi:
        return "Synthetic, same ratio"
    if row.t4_null:
        return "TSH-only (T4 CI includes 0)"
    if row.r > 0:
        return "Concordant (T4 rises with TSH)"
    return "Partial, opposite at a different ratio"

A["pattern"] = A.apply(cat, axis=1)
vc = A.pattern.value_counts()
L.append("\n  Pattern composition:")
for k, v in vc.items():
    L.append(f"    {k:<40}{v:>4}  ({100*v/len(A):.1f}%)")

# the key question: is TSH-only the norm
n_tsh_only = int(A.t4_null.sum())
L.append(f"\n[Q2] Among TSH-significant genes, {n_tsh_only}/{len(A)} "
         f"({100*n_tsh_only/len(A):.1f}%) have a T4 95% CI that includes 0")
L.append(f"  If the 'TSH-only' share is high -> the SCRN2 signature is the norm at this layer, not specific")
L.append(f"  If the share is low -> SCRN2 is in the minority and the signature carries information")

# ---------- Q3: where SCRN2 sits ----------
s = w[w.gene == "SCRN2"].iloc[0]
pct = 100 * (A.r < s.r).mean()
# two-sided tail probability
tail = 2 * min((A.r >= s.r).mean(), (A.r <= s.r).mean())
L.append(f"\n[Q3] SCRN2: b_TSH={s.b_TSH:+.4f}, b_T4={s.b_T4:+.4f}, "
         f"r={s.r:+.3f} (95% CI [{s.r_lo:+.3f}, {s.r_hi:+.3f}])")
L.append(f"  percentile within the TSHSig set: {pct:.1f}% (of {len(A)} genes)")
L.append(f"  two-sided empirical p against that distribution: {tail:.3f}")
L.append(f"  agreement with the TPO ratio: {'yes' if s.r_lo <= K_TPO <= s.r_hi else 'no'}")

# ---------- test: can a single ratio axis explain the whole distribution ----------
L.append("\n[4] Test of the single-ratio-axis hypothesis")
L.append("  If every gene acted through one pathway, r would cluster around the TPO slope.")
inside = int(((A.r_lo <= K_TPO) & (K_TPO <= A.r_hi)).sum())
p_binom = stats.binomtest(inside, len(A), 0.95, alternative="less").pvalue
L.append(f"  Genes whose 95% CI for r covers {K_TPO:+.2f}: {inside}/{len(A)} = {100*inside/len(A):.1f}%")
L.append(f"  If the true proportion were 95% (a single pathway), the probability of observing this few is p = {p_binom:.2e}")
L.append(f"  -> the single-ratio axis {'does not hold' if p_binom < 1e-3 else 'cannot be excluded'}")

# weighted correlation
wx = 1 / (w.se_TSH ** 2); wy = 1 / (w.se_T4 ** 2)
sub = w[(w.b_TSH.abs() > 0.005)]
rho_s = stats.spearmanr(sub.b_TSH, sub.b_T4)
L.append(f"\n  Across all {len(sub)} genes, Spearman rho between b_TSH and b_T4 = {rho_s[0]:+.3f} (p={rho_s[1]:.2e})")
L.append(f"  Negative -> the two axes share feedback structure; weak -> many genes move at one end only")

w.to_csv(os.path.join(OUT, "s08_plane_genes.csv"), index=False, encoding="utf-8-sig")
A.to_csv(os.path.join(OUT, "s08_ratio_summary.csv"), index=False, encoding="utf-8-sig")

L.append("\n[5] Distribution of the T4 effect among TSH-significant genes (ordered by |z_T4|, top 12)")
top = A.reindex(A.b_T4.abs().sort_values(ascending=False).index).head(12)
L.append(f"{'gene':<10}{'b_TSH':>9}{'b_T4':>9}{'r':>8}{'r_CI':>20}  pattern")
for _, x in top.iterrows():
    L.append(f"{x.gene:<10}{x.b_TSH:>9.4f}{x.b_T4:>9.4f}{x.r:>8.2f}"
             f"   [{x.r_lo:>6.2f},{x.r_hi:>6.2f}]  {x.pattern}")

# where the known thyroid genes sit
L.append("\n[6] Where known thyroid genes sit on the plane")
KNOWN = ["TPO", "TG", "SLC5A5", "TSHR", "DUOX2", "DUOXA2", "IYD", "SLC26A4",
         "PAX8", "NKX2-1", "FOXE1", "DIO1", "DIO2", "SLC16A10", "THRA", "THRB",
         "TSHB", "SCRN2"]
L.append(f"{'gene':<10}{'b_TSH':>9}{'p_TSH':>10}{'b_T4':>9}{'p_T4':>10}{'r':>8}  call")
for g in KNOWN:
    x = w[w.gene == g]
    if not len(x):
        continue
    x = x.iloc[0]
    verdict = ("TSH sig" if x.fdr_TSH < 0.05 else "") + (" / T4 sig" if x.fdr_T4 < 0.05 else "")
    verdict = verdict.strip(" /") or "neither"
    L.append(f"{g:<10}{x.b_TSH:>9.4f}{x.p_TSH:>10.2e}{x.b_T4:>9.4f}{x.p_T4:>10.2e}"
             f"{x.r:>8.2f}  {verdict}")

txt = "\n".join(L)
open(os.path.join(OUT, "s08_plane.txt"), "w", encoding="utf-8").write(txt)
print(txt)
