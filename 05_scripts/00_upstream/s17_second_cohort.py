#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s17_second_cohort.py
====================
A second, independent European TSH GWAS, used as a check.

Cohort
------
GCST90319321  Thyroid stimulating hormone levels, European, N=164,818, GWAMA
(local copy DATA17_repro_checks/raw/GCST90319321.h.tsv.gz, md5 13e0417c... verified)

What it does
------------
1) stream the chr17:45.5-46.0 Mb region out of the 1.1 GB genome-wide file
2) take the outcome effect of the three SCRN2 eQTL instrument SNPs, aligned by allele
3) Wald ratio = b_out / b_eqtl, standard error by the delta method
4) fixed-effect pooling with Rand 2025 and the South Asian GCST90448634

Output
------
second_cohort_region.tsv        full statistics for the region (for plotting)
second_cohort_mr.csv            Wald ratios of the three instruments in this cohort
three_cohort_meta.csv           the three-cohort pooled table
"""
import gzip, os, sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "..", "DATA17_repro_checks", "raw")
REG = os.path.join(HERE, "..", "..", "DATA15_thyroid_eQTL", "results", "regional")

FILE = os.path.join(RAW, "GCST90319321.h.tsv.gz")
SAS = os.path.join(RAW, "GCST90448634.h.tsv.gz")
LO, HI = 47_250_000, 48_450_000   # this harmonised file is actually GRCh38 (its meta says GRCh37, which the positions contradict)
TOOLS = ["rs6503875", "rs62076102", "rs4794224"]


def scan(fp, lo, hi, want=None):
    """Stream the requested chr17 interval; when `want` is given keep only those rsIDs."""
    rows = []
    with gzip.open(fp, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        need = ["base_pair_location", "effect_allele", "other_allele",
                "beta", "standard_error", "p_value", "rsid"]
        if "chromosome" in ix:
            need.append("chromosome")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue
            if "chromosome" in ix and p[ix["chromosome"]] not in ("17", "chr17"):
                continue
            try:
                pos = int(p[ix["base_pair_location"]])
            except ValueError:
                continue
            if pos < lo or pos > hi:
                continue
            rec = {k: p[ix[k]] for k in need if k in ix}
            if want and rec.get("rsid") not in want:
                continue
            rows.append(rec)
    d = pd.DataFrame(rows)
    for c in ["base_pair_location", "beta", "standard_error", "p_value"]:
        if c in d:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def eqtl_lookup():
    d = pd.read_csv(os.path.join(REG, "SCRN2__TSH.tsv"), sep="\t")
    d = d[d.iloc[:, 0].isin(TOOLS)]
    return d.set_index(d.columns[0])[["beta_eqtl", "se_eqtl", "ea", "oa", "pos37"]]


def align(b_out, ea_out, oa_out, ea_l, oa_l):
    """Align the outcome beta to the effect allele of the local eQTL record."""
    if ea_out == ea_l and oa_out == oa_l:
        return b_out, True
    if ea_out == oa_l and oa_out == ea_l:
        return -b_out, True
    return np.nan, False


def main():
    print("scanning the GCST90319321 region ...")
    reg = scan(FILE, LO, HI)
    print(f"  variants in the region: {len(reg)}")
    reg.to_csv(os.path.join(HERE, "second_cohort_region.tsv"), sep="\t", index=False)

    if len(reg):
        top = reg.nsmallest(5, "p_value")
        print("  five strongest in the region:")
        for _, r in top.iterrows():
            print(f"    {r.rsid:12s} pos={int(r.base_pair_location)} "
                  f"β={r.beta:+.4f} p={r.p_value:.2e} ({r.effect_allele}/{r.other_allele})")
        gws = reg[reg.p_value < 5e-8]
        print(f"  p<5e-8 in the region: {len(gws)}")

    eq = eqtl_lookup()
    print("\nthe three instrument SNPs in this cohort:")
    rows = []
    for s in TOOLS:
        if s not in eq.index:
            print(f"  {s}: absent from the eQTL table")
            continue
        e = eq.loc[s]
        hit = reg[reg.rsid == s]
        if not len(hit):
            print(f"  {s}: not reported in this cohort (probably QC-filtered)")
            rows.append(dict(snp=s, status="absent"))
            continue
        h = hit.iloc[0]
        b, ok = align(h.beta, h.effect_allele, h.other_allele, e.ea, e.oa)
        if not ok:
            print(f"  {s}: allele mismatch ({h.effect_allele}/{h.other_allele} vs {e.ea}/{e.oa})")
            rows.append(dict(snp=s, status="allele_mismatch"))
            continue
        wald = b / e.beta_eqtl
        se_w = h.standard_error / abs(e.beta_eqtl)
        rows.append(dict(snp=s, status="ok", beta_out=b, se_out=h.standard_error,
                         p_out=h.p_value, beta_eqtl=e.beta_eqtl, se_eqtl=e.se_eqtl,
                         wald=wald, se_wald=se_w, z=wald / se_w,
                         p_wald=2 * stats.norm.sf(abs(wald / se_w))))
        print(f"  {s}: β_out={b:+.4f} (se {h.standard_error:.4f}) p_out={h.p_value:.3g} "
              f"→ Wald={wald:+.4f} (se {se_w:.4f}) p={2*stats.norm.sf(abs(wald/se_w)):.3g}")

    mr = pd.DataFrame(rows)
    mr.to_csv(os.path.join(HERE, "second_cohort_mr.csv"), index=False)

    ok = mr[mr.status == "ok"]
    if len(ok):
        # within-cohort IVW (inverse-variance weighting of the Wald ratios)
        w = 1 / ok.se_wald ** 2
        est = (w * ok.wald).sum() / w.sum()
        se = np.sqrt(1 / w.sum())
        q = (w * (ok.wald - est) ** 2).sum()
        pq = stats.chi2.sf(q, len(ok) - 1)
        print(f"\nwithin-cohort IVW ({len(ok)} reported variants): {est:+.4f} "
              f"(95% CI {est-1.96*se:+.4f} ~ {est+1.96*se:+.4f}), "
              f"p={2*stats.norm.sf(abs(est/se)):.3g}, Q p={pq:.3f}")

        # ---- three-cohort check: always the lead SNP rs6503875, so the same instrument is compared ----
        LEAD = "rs6503875"
        lead = ok[ok.snp == LEAD]
        if not len(lead):
            print(f"[warn] this cohort does not report {LEAD}; the same-instrument pooling cannot be done")
            return
        lr = lead.iloc[0]
        cohorts = [dict(cohort="Rand2025_EUR", instrument=LEAD,
                        est=0.03359, se=0.00885, p=1.48e-4, n=482873)]
        cohorts.append(dict(cohort="GCST90319321_EUR", instrument=LEAD,
                            est=float(lr.wald), se=float(lr.se_wald),
                            p=float(lr.p_wald), n=164818))
        sas_fp = os.path.join(HERE, "..", "..", "DATA17_repro_checks", "results",
                              "s03_replication_lead.csv")
        if os.path.exists(sas_fp):
            s3 = pd.read_csv(sas_fp)
            print("\nthe South Asian file:")
            print(s3.to_string(index=False))
            cmap = {c.lower(): c for c in s3.columns}
            for _, r in s3.iterrows():
                wkey = cmap.get("wald") or cmap.get("b_mr_rep")
                skey = cmap.get("se") or cmap.get("se_mr_rep")
                wv = pd.to_numeric(r.get(wkey, np.nan), errors="coerce")
                sv = pd.to_numeric(r.get(skey, np.nan), errors="coerce")
                if np.isfinite(wv) and np.isfinite(sv):
                    cohorts.append(dict(cohort="GCST90448634_SAS", instrument=LEAD,
                                        est=float(wv), se=float(sv),
                                        p=float(2 * stats.norm.sf(abs(wv / sv))),
                                        n=38000))
        cm = pd.DataFrame(cohorts)
        ww = 1 / cm.se ** 2
        ce = (ww * cm.est).sum() / ww.sum()
        cse = np.sqrt(1 / ww.sum())
        cq = (ww * (cm.est - ce) ** 2).sum()
        cm.to_csv(os.path.join(HERE, "three_cohort_meta.csv"), index=False)
        print("\nthree-cohort fixed-effect pooling (same instrument rs6503875):")
        print(cm.to_string(index=False))
        i2 = max(0.0, (cq - (len(cm) - 1)) / cq) * 100 if cq > 0 else 0.0
        print(f"  pooled b={ce:+.4f} (95% CI {ce-1.96*cse:+.4f} to {ce+1.96*cse:+.4f}), "
              f"p={2*stats.norm.sf(abs(ce/cse)):.3g}, "
              f"Q p={stats.chi2.sf(cq, len(cm)-1):.3f}, I²={i2:.0f}%")


if __name__ == "__main__":
    main()
