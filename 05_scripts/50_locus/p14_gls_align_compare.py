#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
[Superseded by p15_ld_true_dosage.py; kept only as a methodological comparison]
p14_gls_align_compare.py
========================
Reconstructs how the three-instrument GLS estimate (0.0250 / 0.0080) in
s21_reference_frame.py was actually obtained, and answers one methodological
question:

  the s21 formulation  S = outer(se,se) * R * outer(sign(gamma), sign(gamma))
  flips the sign of the LD correlation matrix according to sign(gamma). If
  sign(gamma) happened to equal "is the effect allele ALT", this would be standard
  allele alignment; but the sign of gamma comes from the biological direction of
  the effect, not from the allele coding, and the two agree only when the effect
  allele is ALT and the effect is positive. That is therefore an implicit
  assumption, and it has to be checked against a genuinely REF/ALT-aligned LD
  matrix.

This script reports four algorithms:
  A  R = tool_R.npy (ALT dosage, includes negative correlations) + flip by sign(gamma)   <- the original s21 approach
  B  R = tool_ld_matrix.csv (LDlink, positive correlations) + no flipping                <- the other common formulation
  C  R = 1000G allele-aligned (ld_TSH_aligned), alternative triplet when rs4794224 is missing
  D  allele-aligned LD, per-instrument Wald + IVW (no Q weighting; reference only)

Output: tables/19_gls_sign_handling_comparison.csv and tables/19b_three_instruments_per_snp_wald.csv
"""
import math
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
T = f"{BASE}/tables"
D18 = "D:/NewBeginning/MR1/DATA18_manuscript/analysis"
REG = "D:/NewBeginning/MR1/DATA15_thyroid_eQTL/results/regional"

TOOLS3 = ["rs6503875", "rs62076102", "rs4794224"]        # the three used by s21
TOOLS3_ALT = ["rs6503875", "rs62076102", "rs11079797"]   # the three available in the 1000G panel


def wald(bg, seg, bo, seo):
    if bg == 0 or bo == 0:
        return np.nan, np.nan
    b = bo / bg
    se = abs(b) * math.sqrt((seo / bo) ** 2 + (seg / bg) ** 2)
    return b, se


def gls(b, se, R):
    S = np.outer(se, se) * R
    Si = np.linalg.inv(S)
    one = np.ones(len(b))
    est = float(one @ Si @ b / (one @ Si @ one))
    se_est = float(math.sqrt(1 / (one @ Si @ one)))
    Q = float((b - est) @ Si @ (b - est))
    df = len(b) - 1
    return est, se_est, Q, df, (stats.chi2.sf(Q, df) if df > 0 else np.nan)


# ---------- 1. per-instrument Wald ratio ----------
rows = []
store = {}
for oc in ["TSH", "T4"]:
    d = pd.read_csv(f"{REG}/SCRN2__{oc}.tsv", sep="\t")
    d = d.set_index("snp")
    for s in TOOLS3 + ["rs11079797"]:
        if s not in d.index:
            continue
        r = d.loc[s]
        b, se = wald(r.beta_eqtl, r.se_eqtl, r.beta, r.se)
        rows.append(dict(pheno=oc, snp=s, ea=r.ea, oa=r.oa,
                         beta_eqtl=r.beta_eqtl, se_eqtl=r.se_eqtl,
                         beta_out=r.beta, se_out=r.se,
                         F=(r.beta_eqtl / r.se_eqtl) ** 2,
                         b_wald=b, se_wald=se,
                         p_wald=2 * stats.norm.sf(abs(b / se))))
        store[(oc, s)] = (b, se, r.beta_eqtl)
wd = pd.DataFrame(rows)
wd.to_csv(f"{T}/19b_three_instruments_per_snp_wald.csv", index=False, encoding="utf-8-sig")
print("=== per-instrument Wald ratio (thyroid eQTL) ===")
print(wd[wd.pheno == "TSH"].round(5).to_string(index=False))

# ---------- 2. four LD handlings ----------
R_npy = np.load(f"{D18}/tool_R.npy")
nm = pd.read_csv(f"{D18}/tool_ld_matrix.csv", index_col=0)
R_csv = nm.values
names_csv = nm.index.tolist()

ld_al = pd.read_csv(f"{BASE}/inputs_ld/_deprecated_gwas_ea/ld_TSH_aligned.tsv.gz", sep="\t", index_col=0)
ld_al.index = ld_al.index.astype(int)
ld_al.columns = ld_al.columns.astype(int)
pos_of = {"rs6503875": 45915395, "rs62076102": 45887132, "rs11079797": 45912302}


def get_R(tools, mode):
    if mode == "npy":
        sub = pd.DataFrame(R_npy, index=names_csv, columns=names_csv).loc[tools, tools]
        return sub.values
    if mode == "csv":
        sub = nm.loc[tools, tools]
        return sub.values
    if mode == "aligned":
        p = [pos_of[s] for s in tools]
        return ld_al.loc[p, p].values
    raise ValueError(mode)


out = []
for pheno in ["TSH", "T4"]:
    for tools, tag in [(TOOLS3, "three instruments (s21: rs6503875+rs62076102+rs4794224)"),
                       (TOOLS3_ALT, "three instruments (1000G panel: rs6503875+rs62076102+rs11079797)"),
                       (["rs6503875", "rs62076102"], "two instruments (rs6503875+rs62076102)"),
                       (["rs6503875"], "single instrument (rs6503875)")]:
        bs = np.array([store[(pheno, s)][0] for s in tools])
        ses = np.array([store[(pheno, s)][1] for s in tools])
        gam = np.array([store[(pheno, s)][2] for s in tools])
        modes = ["npy", "csv", "aligned"] if not any(s == "rs4794224" for s in tools) else ["npy", "csv"]
        for mode in modes:
            R = get_R(tools, mode)
            for flip in ([False, True] if mode != "aligned" else [False]):
                Rf = R * (np.outer(np.sign(gam), np.sign(gam)) if flip else 1.0)
                try:
                    est, se_est, Q, df, pq = gls(bs, ses, Rf)
                except np.linalg.LinAlgError:
                    continue
                out.append(dict(pheno=pheno, tools=tag, k=len(tools),
                                ld_source=mode, sign_gamma_flip=flip,
                                b=est, se=se_est,
                                p=2 * stats.norm.sf(abs(est / se_est)),
                                Q=Q, df=df, p_Q=pq,
                                gamma_signs=";".join("+" if g > 0 else "-" for g in gam)))
res = pd.DataFrame(out)
res.to_csv(f"{T}/19_gls_sign_handling_comparison.csv", index=False, encoding="utf-8-sig")
print("\n=== GLS across the four LD handlings ===")
print(res[res.pheno == "TSH"].round(6).to_string(index=False))

print("\n=== comparison against s21 reference_frame.csv ===")
ref = pd.read_csv(f"{D18}/reference_frame.csv")
print(ref[ref.comparison.str.contains("GLS")][["comparison", "b_tsh", "se_tsh", "b_t4", "se_t4"]].to_string(index=False))
hit = res[(res.pheno == "TSH") & (np.abs(res.b - 0.025) < 0.0015)]
print("\nvariants that land near 0.025:")
print(hit.round(6).to_string(index=False) if len(hit) else "  none")

# ---------- 3. sensitivity: effect on the conclusions of the 0.0336 single-instrument vs the 0.025 GLS estimate ----------
print("\n=== effect on the section 12 conclusions ===")
for lab, b, s in [("single instrument 0.0336 (0.0088)", 0.0335877, 0.0088515)] + \
                 [(f"{r.tools}/{r.ld_source}/flip={r.sign_gamma_flip}", r.b, r.se)
                  for r in res[(res.pheno == "TSH") & (res.k == 3)].itertuples()]:
    req = -0.206 * b
    se_req = math.sqrt((-0.206 * s) ** 2 + (b * 0.013) ** 2)
    t4, se4 = 0.004976, 0.011956
    z = (t4 - req) / math.sqrt(se4 ** 2 + se_req ** 2)
    print(f"  {lab:52s} b={b:+.4f} req={req:+.4f} z={z:+.2f} p={2*stats.norm.sf(abs(z)):.3f}")
