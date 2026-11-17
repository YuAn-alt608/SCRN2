# -*- coding: utf-8 -*-
"""
Step 16: check the 733 FDR-significant hits in an independent sample, FinnGen R12

The outcomes in the original analysis come from the EBI GWAS Catalog and large thyroid-function
meta-analyses of UK Biobank lineage; FinnGen R12 is an independent cohort, so it can test whether
the instruments of a given gene give an effect in the same direction in another cohort.

The key control, which is more than repetition, is cross-trait directional consistency:
    higher TSH -> higher hypothyroidism risk (same direction)
    higher free T4 -> lower hypothyroidism risk (opposite direction)
If both expectations hold, the check result is not merely within-cohort correlation but is
consistent with thyroid-axis physiology.

Method:
  eQTLGen supplies only Z, so beta is recovered as beta = z / sqrt(2*maf*(1-maf)*(N+z^2)), with
  MAF approximated by FinnGen af_alt (both are European-ancestry samples). Alleles are harmonised
  to FinnGen. Multiple instruments are combined by IVW; a single instrument is the Wald ratio.

Output:
  results/finngen_replication.csv   per-hit replication results
  results/finngen_summary.txt       summary by outcome
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
TMP = os.path.join(BASE, "tmp")
DATA = os.path.join(BASE, "data")

N_EQTL = 31684  # eQTLGen sample size

# outcome code -> (FinnGen file, expected direction relative to the dichotomous outcome)
#   +1 : a higher value of this trait should raise hypothyroidism risk
#   -1 : a higher value of this trait should lower hypothyroidism risk
#    0 : this is itself a hypothyroidism endpoint, so the direction is the same by construction
MAP = {
    "HYPO": ("finngen_E4_HYTHY_AI_STRICT.tsv", 0, "hypothyroidism (autoimmune, independent cohort, direct check)"),
    "TSH":  ("finngen_E4_HYTHY_AI_STRICT.tsv", +1, "a TSH hit should move hypothyroidism risk the same way"),
    "T4":   ("finngen_E4_HYTHY_AI_STRICT.tsv", -1, "a free T4 hit should move hypothyroidism risk the other way"),
    "HASH": ("finngen_E4_THYROID.tsv", 0, "Hashimoto -> thyroid disease (same cohort, other endpoint, not independent)"),
}


def load_fg(fn):
    d = pd.read_csv(os.path.join(TMP, fn), sep="\t", dtype={"#chrom": str})
    d.columns = [c.lstrip("#").strip() for c in d.columns]
    d = d.rename(columns={"rsids": "snp", "alt": "ea", "ref": "nea",
                          "sebeta": "se", "af_alt": "af"})
    d = d[["snp", "ea", "nea", "beta", "se", "pval", "af"]].dropna(subset=["snp"])
    # a locus may appear on several rows (multiple rsids); keep the first causal-position row
    d = d.drop_duplicates("snp")
    return d.set_index("snp")


COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
fg_cache = {}


def replicate(genes_df, fg):
    """genes_df: one row per gene with snp/ea/nea/z. Returns the IVW estimate."""
    num = den = 0.0
    used = 0
    for r in genes_df.itertuples():
        if r.snp not in fg.index:
            continue
        f = fg.loc[r.snp]
        # harmonise to the eQTL effect allele
        if r.ea == f["ea"] or COMP.get(r.ea) == f["ea"]:
            sgn = 1.0
        elif r.ea == f["nea"] or COMP.get(r.ea) == f["nea"]:
            sgn = -1.0
        else:
            continue
        if not np.isfinite(f["beta"]) or not np.isfinite(f["se"]) or f["se"] <= 0:
            continue
        af = f["af"]
        if not np.isfinite(af):
            continue
        maf = min(af, 1 - af)
        if maf <= 0.001:
            continue
        info = 2 * maf * (1 - maf) * (N_EQTL + r.z ** 2)
        if info <= 0:
            continue
        b_exp = r.z / np.sqrt(info)
        if abs(b_exp) < 1e-12:
            continue
        b_out = sgn * f["beta"]
        se_out = f["se"]
        w = b_exp / se_out            # Wald ratio numerator
        den += (b_exp / se_out) ** 2  # IVW weight
        num += w * (b_out / se_out)
        used += 1
    if used == 0 or den <= 0:
        return np.nan, np.nan, 0
    b = num / den
    se = 1 / np.sqrt(den)
    return b, se, used


def main():
    mr = pd.read_csv(os.path.join(RES, "MR_annotated.csv"), encoding="utf-8-sig")
    mr.columns = [c.strip() for c in mr.columns]
    sig = mr[mr["fdr"].notna() & (mr["fdr"] < 0.05)].copy()
    iv = pd.read_csv(os.path.join(DATA, "ivs.tsv"), sep="\t")
    map_name = {"TSH": "TSH", "Free T4": "T4",
                "Hypothyroidism": "HYPO", "Hashimoto thyroiditis": "HASH"}
    sig["oc"] = sig["outcome"].map(map_name)
    print(f"hits to check: {len(sig)}")

    rows = []
    for oc, (fn, expect, note) in MAP.items():
        if fn not in fg_cache:
            fg_cache[fn] = load_fg(fn)
            print(f"  loading {fn}: {len(fg_cache[fn])} SNPs")
        fg = fg_cache[fn]

        sub = sig[sig["oc"] == oc]
        print(f"\n== {oc} ({len(sub)} hits) vs {fn} | {note}")
        for r in sub.itertuples():
            ins = iv[iv["gene"] == r.gene]
            b, se, n_used = replicate(ins, fg)
            ok = np.isfinite(b) and np.isfinite(se) and se > 0
            p = float(2 * (1 - norm.cdf(abs(b / se)))) if ok else np.nan
            rows.append({
                "gene": r.gene, "outcome": r.outcome, "oc": oc,
                "b_orig": r.b, "p_orig": r.pval,
                "snp": r.snp, "nsnp_orig": r.nsnp, "nsnp_used": n_used,
                "b_fg": b if ok else np.nan, "se_fg": se if ok else np.nan,
                "p_fg": p, "expect_sign": expect,
            })

    out = pd.DataFrame(rows)
    out["sign_ok"] = np.sign(out["b_fg"]) == np.sign(out["b_orig"])
    out["p_fg"] = pd.to_numeric(out["p_fg"], errors="coerce")
    # two different criteria:
    #   sign_ok      = the two betas share a sign (meaningful for HYPO/TSH/HASH only)
    #   physio_ok    = consistent with the thyroid axis (a T4 hit should point the other way)
    #   replicated   = same sign plus nominal significance (systematically understated for T4, so
    #                  the physio version is reported alongside)
    out["physio_ok"] = np.where(out["expect_sign"] < 0, ~out["sign_ok"], out["sign_ok"])
    out["replicated"] = out["sign_ok"] & (out["p_fg"] < 0.05)
    out["replicated_physio"] = out["physio_ok"] & (out["p_fg"] < 0.05)
    out.to_csv(os.path.join(RES, "finngen_replication.csv"), index=False, encoding="utf-8-sig")

    print("\n=== summary ===")
    lines = []
    for oc in MAP:
        s = out[(out["oc"] == oc) & out["b_fg"].notna()]
        if len(s) == 0:
            continue
        line = (f"{oc:5s} n={len(s):3d} | same dir. {int(s['sign_ok'].sum()):3d} "
                f"({100*s['sign_ok'].mean():.1f}%) | same direction and p<0.05 "
                f"{int(s['replicated'].sum()):3d} ({100*s['replicated'].mean():.1f}%)"
                f" | matches the physiological expectation and p<0.05 "
                f"{int(s['replicated_physio'].sum()):3d} ({100*s['replicated_physio'].mean():.1f}%)")
        print(line)
        lines.append(line)

    # cross-trait directional consistency: a TSH hit should agree, a T4 hit should oppose
    print("\n=== cross-trait directional consistency (against hypothyroidism risk) ===")
    for oc, expect in [("TSH", +1), ("T4", -1)]:
        s = out[(out["oc"] == oc) & out["b_fg"].notna()]
        if len(s) == 0:
            continue
        if expect > 0:
            agree = s["sign_ok"].mean()
        else:
            agree = (~s["sign_ok"]).mean()
        line = f"{oc:5s} fraction matching the physiological expectation {100*agree:.1f}% (n={len(s)})"
        print(line)
        lines.append(line)

    with open(os.path.join(RES, "finngen_summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\nwritten: results/finngen_replication.csv / finngen_summary.txt")

    # list the hits that replicated
    top = out[out["replicated"]].nsmallest(10, "p_fg")
    if len(top):
        print("\n=== the 10 strongest replicated hits ===")
        print(top[["gene", "outcome", "b_orig", "b_fg", "p_fg", "nsnp_used"]]
              .to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
