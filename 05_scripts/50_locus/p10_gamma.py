#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p10_gamma.py
============
The same variant in the same tissue under two gamma estimates: eQTL Catalogue
(QTD000341 re-run) versus the official GTEx signif pairs. This quantifies the
effect on the cis-MR point estimate -- the point a reviewer is most likely to
attack. The whole-blood eQTLGen Z-based result is reported alongside as a third,
independent reference.
"""
import gzip
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1"
OUT = f"{BASE}/DATA21_complete_2026-11-14"
GEN = {"ENSG00000141295": "SCRN2", "ENSG00000198933": "TBKBP1", "ENSG00000159111": "MRPL10"}

# variant_id -> rsid (from the eQTL Catalogue extract)
cat = pd.read_csv(f"{OUT}/raw/thyroid_region_g38.tsv", sep="\t",
                  usecols=["variant", "rsid", "pvalue", "beta", "se", "gene_id"])
vid2rs = dict(zip(cat.variant, cat.rsid))
cat["vkey"] = cat.variant.str.replace("_b38", "", regex=False)

# rsid -> pos37
rs2pos = {}
for oc in ["TSH", "T4", "HYPO"]:
    d = pd.read_csv(f"{BASE}/DATA15_thyroid_eQTL/results/gwas_win/{oc}.tsv.gz", sep="\t",
                    usecols=["rs_id", "base_pair_location"])
    for rs, pos in zip(d.rs_id, d.base_pair_location):
        rs2pos.setdefault(rs, int(pos))

tsh = pd.read_csv(f"{BASE}/DATA15_thyroid_eQTL/results/gwas_win/TSH.tsv.gz", sep="\t")
tsh = tsh.rename(columns={"rs_id": "rsid", "base_pair_location": "pos37"})

# official GTEx
rows = []
with gzip.open(f"{BASE}/DATA15_thyroid_eQTL/tmp/GTEx_Analysis_v8_eQTL/"
               "Thyroid.v8.signif_variant_gene_pairs.txt.gz", "rt") as f:
    head = f.readline().rstrip("\n").split("\t")
    ix = {k: i for i, k in enumerate(head)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        gid = p[ix["gene_id"]].split(".")[0]
        if gid not in GEN:
            continue
        vid = p[ix["variant_id"]]
        rows.append(dict(gene=GEN[gid], gid=gid, variant=vid, alt=vid.split("_")[-2],
                         ref=vid.split("_")[-3],
                         p_gtex=float(p[ix["pval_nominal"]]),
                         b_gtex=float(p[ix["slope"]]), se_gtex=float(p[ix["slope_se"]]),
                         maf=float(p[ix["maf"]])))
gt = pd.DataFrame(rows)
gt["vkey"] = gt.variant.str.replace("_b38", "", regex=False)
gt["rsid"] = gt.vkey.map(vid2rs)
gt = gt[gt.rsid.notna()]

out = []
for g, sub in gt.groupby("gene"):
    i = sub.p_gtex.idxmin()
    r = sub.loc[i]
    t = tsh[tsh.rsid == r.rsid]
    if len(t) == 0:
        continue
    t = t.iloc[0]
    flip = t.effect_allele != "C" and t.other_allele == r.alt
    b_out = t.beta if t.effect_allele == r.alt else -t.beta if t.effect_allele == r.ref else np.nan
    if not np.isfinite(b_out):
        continue
    w_g = b_out / r.b_gtex
    se_g = abs(w_g) * np.sqrt((t.standard_error / b_out) ** 2 + (r.se_gtex / r.b_gtex) ** 2)
    c = cat[(cat.variant == r.vkey) & (cat.gene_id == r.gid)]
    if len(c) == 0:
        continue
    c = c.iloc[0]
    w_c = b_out / c.beta
    se_c = abs(w_c) * np.sqrt((t.standard_error / b_out) ** 2 + (c.se / c.beta) ** 2)
    out.append(dict(gene=g, lead_snp=r.rsid, pos37=rs2pos.get(r.rsid, np.nan),
                    alt=r.alt,
                    gamma_gtex=r.b_gtex, se_gtex=r.se_gtex, p_gtex=r.p_gtex,
                    gamma_eqtlcat=c.beta, se_eqtlcat=c.se, p_eqtlcat=c.pvalue,
                    gamma_ratio=abs(c.beta / r.b_gtex),
                    b_tsh=b_out, se_tsh=t.standard_error,
                    wald_gtex=w_g, se_wald_gtex=se_g,
                    p_wald_gtex=2 * stats.norm.sf(abs(w_g / se_g)),
                    wald_eqtlcat=w_c, se_wald_eqtlcat=se_c,
                    p_wald_eqtlcat=2 * stats.norm.sf(abs(w_c / se_c)),
                    wald_diff_pct=100 * (w_c - w_g) / w_g))
t = pd.DataFrame(out)
t.to_csv(f"{OUT}/tables/17_two_gamma_sensitivity.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 240)
print(t.round(5).to_string(index=False))
