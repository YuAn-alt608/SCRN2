#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SCRN2 sensitivity analysis: Steiger directionality + MR-Egger intercept + leave-one-out + multi-instrument MR
=====================================================================
Aim: strengthen the direction of "SCRN2 expression -> thyroid function" and
rule out reverse causation.

Methods:
1. Steiger directionality test: compare the variance explained by the eQTL
   (R2_exposure) with that explained by the outcome (R2_outcome) to infer the
   causal direction (the exposure should explain more variance).
2. Multi-instrument MR (independent SNPs after clumping): IVW + MR-Egger (with
   intercept) + weighted median. The Egger intercept tests for horizontal
   pleiotropy.
3. Leave-one-out: re-estimate with each SNP dropped in turn, to test whether the
   estimate is driven by a single point.

Data: tmp/susie/SCRN2__<outcome>.tsv (harmonized eQTL beta/varbeta + outcome beta/se)
      tmp/susie/SCRN2__<outcome>.ld.csv.gz (LD matrix)
"""
import pandas as pd
import numpy as np
import os

os.makedirs("results", exist_ok=True)

# ============ configuration ============
OUTCOMES = {
    "TSH":  dict(n_out=482873, type="quant"),
    "HYPO": dict(n_out=1046991, s=98805/1046991, type="cc"),
    "HASH": dict(n_out=395640, s=15654/395640, type="cc"),
}
IV = "rs11079797"

# ============ helper functions ============
def r2_from_z(z, n, maf):
    """Back out the variance explained from a z score: R2 = z^2/(z^2 + n - 2)."""
    return z**2 / (z**2 + n - 2)

def steiger_test(beta_eqtl, se_eqtl, n_eqtl, maf, beta_out, se_out, n_out):
    """Steiger directionality test: the exposure should explain more variance."""
    z_eqtl = beta_eqtl / se_eqtl
    z_out = beta_out / se_out
    r2_eqtl = r2_from_z(z_eqtl, n_eqtl, maf)
    r2_out = r2_from_z(z_out, n_out, maf)
    # Steiger Z test
    # standard formula: compare the two (corrected) R2 values
    return r2_eqtl, r2_out, r2_eqtl > r2_out

def clump_independent(snps, R, pval, r2_cut=0.1):
    """Greedy clumping: sort by p, keep SNPs with LD r2 below the cut."""
    order = np.argsort(pval)
    keep = []
    for i in order:
        if all(abs(R[i, j]) < np.sqrt(r2_cut) for j in keep):
            keep.append(i)
    return sorted(keep)

def ivw(beta, se):
    """Inverse-variance weighting."""
    w = 1 / se**2
    b = np.sum(beta * w) / np.sum(w)
    se_b = np.sqrt(1 / np.sum(w))
    return b, se_b

def mr_egger(beta_x, se_x, beta_y, se_y):
    """MR-Egger regression (with intercept)."""
    # weighted regression using se_y as weight: beta_y ~ beta_x
    w = 1 / se_y**2
    x = beta_x
    y = beta_y
    # weighted linear regression
    W = np.diag(w)
    X = np.column_stack([np.ones(len(x)), x])
    XtW = X.T @ W
    beta_hat = np.linalg.solve(XtW @ X, XtW @ y)
    intercept, slope = beta_hat
    # standard errors
    resid = y - X @ beta_hat
    dof = len(x) - 2
    if dof > 0:
        sigma2 = (resid @ W @ resid) / dof
        cov = np.linalg.inv(XtW @ X) * sigma2
        se_intercept = np.sqrt(cov[0, 0])
        se_slope = np.sqrt(cov[1, 1])
    else:
        se_intercept = se_slope = np.nan
    return intercept, se_intercept, slope, se_slope

def weighted_median(beta, se, weights=None):
    """Weighted median estimate."""
    if weights is None:
        weights = np.ones(len(beta))
    order = np.argsort(beta)
    beta_s = beta[order]
    w_s = weights[order]
    cum = np.cumsum(w_s)
    total = cum[-1]
    # locate the interval where cumulative weight crosses one half
    idx = np.searchsorted(cum, total/2)
    return beta_s[idx]

# ============ main analysis ============
print("=" * 72)
print("SCRN2 sensitivity analysis (Steiger + MR-Egger + leave-one-out)")
print("=" * 72)

all_rows = []

for oc, cfg in OUTCOMES.items():
    f = f"tmp/susie/SCRN2__{oc}.tsv"
    ld_f = f"tmp/susie/SCRN2__{oc}.ld.csv.gz"
    if not os.path.exists(f) or not os.path.exists(ld_f):
        print(f"\n[skip] {oc}: no data")
        continue

    d = pd.read_csv(f, sep="\t")
    R = pd.read_csv(ld_f, header=None).values
    n_eqtl = int(d["n"].median())

    print(f"\n{'='*72}\n[{oc}] n_eQTL={n_eqtl}, n_outcome={cfg['n_out']}\n{'='*72}")

    # ---- 1. Steiger directionality test (the instrument itself) ----
    iv = d[d.snp == IV].iloc[0]
    z_eqtl_iv = iv["z"]
    z_out_iv = iv["beta"] / iv["se"]
    r2_e, r2_o = r2_from_z(z_eqtl_iv, n_eqtl, iv["maf"]), r2_from_z(z_out_iv, cfg["n_out"], iv["maf"])
    steiger_ok = r2_e > r2_o
    print(f"\n[1] Steiger directionality test (instrument {IV})")
    print(f"    eQTL: z={z_eqtl_iv:.1f}, R2={r2_e:.3e}")
    print(f"    outcome: z={z_out_iv:.1f}, R2={r2_o:.3e}")
    print(f"    R2(exposure) > R2(outcome): {steiger_ok}  -> {'direction consistent (eQTL -> outcome)' if steiger_ok else 'direction uncertain (possible reverse causation)'}")

    # ---- 2. clump down to independent instruments ----
    pval = d["p"].values
    keep = clump_independent(d["snp"].values, R, pval, r2_cut=0.1)
    # make sure the primary instrument is included (or at least retained)
    iv_pos = list(d["snp"]).index(IV)
    if iv_pos not in keep:
        keep = sorted(set(keep) | {iv_pos})
    d_ind = d.iloc[keep].reset_index(drop=True)

    # outcome beta alignment: the effect allele must match (the susie inputs are
    # already harmonized, and the ea column follows the eQTL)
    # both beta_eqtl and beta in the susie input are effects on the same ea
    beta_x = d_ind["beta_eqtl"].values
    se_x = np.sqrt(d_ind["varbeta_eqtl"].values)
    beta_y = d_ind["beta"].values
    se_y = d_ind["se"].values
    maf = d_ind["maf"].values
    snps = d_ind["snp"].values

    n_iv = len(snps)
    print(f"\n[2] multi-instrument MR ({n_iv} independent SNPs after clumping, including {IV})")

    # IVW
    b_ivw, se_ivw = ivw(beta_y, se_y)
    z_ivw = b_ivw / se_ivw
    p_ivw = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(z_ivw)))
    print(f"    IVW:          β={b_ivw:.4f}, se={se_ivw:.4f}, z={z_ivw:.2f}, p={p_ivw:.2e}")

    # MR-Egger
    inter, se_inter, slope, se_slope = mr_egger(beta_x, se_x, beta_y, se_y)
    z_inter = inter / se_inter if se_inter and not np.isnan(se_inter) else np.nan
    p_inter = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(z_inter))) if not np.isnan(z_inter) else np.nan
    print(f"    MR-Egger:     intercept={inter:.4f}, se={se_inter:.4f}, p(intercept)={p_inter:.3f}")
    print(f"                  slope={slope:.4f}, se={se_slope:.4f}")
    egger_interp = "no evidence of horizontal pleiotropy" if not np.isnan(p_inter) and p_inter > 0.05 else ("signal of horizontal pleiotropy" if not np.isnan(p_inter) else "not estimable")
    print(f"                 -> {egger_interp}")

    # weighted median
    w = 1 / se_y**2
    b_wm = weighted_median(beta_y, se_y, weights=w)
    print(f"    weighted median: beta={b_wm:.4f}")

    # ---- 3. leave-one-out ----
    print(f"\n[3] leave-one-out ({n_iv} SNPs)")
    loo_bs = []
    for i in range(n_iv):
        mask = np.ones(n_iv, dtype=bool)
        mask[i] = False
        b, _ = ivw(beta_y[mask], se_y[mask])
        loo_bs.append(b)
    loo_bs = np.array(loo_bs)
    # identify the SNP whose removal shifts the estimate the most
    b_all = b_ivw
    dev = np.abs(loo_bs - b_all)
    worst_i = np.argmax(dev)
    print(f"    IVW full-sample beta={b_all:.4f}")
    print(f"    leave-one-out beta range: [{loo_bs.min():.4f}, {loo_bs.max():.4f}]")
    print(f"    dropping {snps[worst_i]} gives beta={loo_bs[worst_i]:.4f} (largest shift)")
    is_iv_worst = snps[worst_i] == IV
    print(f"    driven by instrument {IV}: {'yes (single-point driven, treat with caution)' if is_iv_worst else 'no (robust, not single-point driven)'}")

    all_rows.append(dict(
        outcome=oc, n_iv=n_iv,
        r2_exposure=r2_e, r2_outcome=r2_o, steiger_ok=steiger_ok,
        b_ivw=b_ivw, se_ivw=se_ivw, p_ivw=p_ivw,
        egger_intercept=inter, egger_intercept_p=p_inter,
        b_egger_slope=slope, b_wm=b_wm,
        loo_min=loo_bs.min(), loo_max=loo_bs.max(),
        loo_worst_snp=snps[worst_i], iv_driven=is_iv_worst,
    ))

res = pd.DataFrame(all_rows)
res.to_csv("results/SCRN2_sensitivity.csv", index=False, encoding="utf-8-sig")
print("\n\nwritten to results/SCRN2_sensitivity.csv")
print(res.to_string(index=False))
