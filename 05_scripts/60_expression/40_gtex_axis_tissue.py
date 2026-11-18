#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
40_gtex_axis_tissue.py - SCRN2 expression on both sides of the thyroid axis
===========================================================================
Reproduces 04_data/06_expression/SCRN2_axis_expression.csv, the table behind the
second and third sentences of section 2.9 (pituitary 59.1 TPM, fifth of 54; the
rs6503875 eQTL stronger in pituitary than in thyroid).

Why this table exists
---------------------
Section 2.9 previously reported thyroid and whole blood only. A claim about the
set-point needs the pituitary value as well, because a gene acting centrally has
to be expressed where the set-point is set. The GTEx v8 release carries all three
tissues, so the comparison is a lookup rather than an experiment.

What it does
------------
1. Asks the GTEx Portal API for the median TPM of one gene in every GTEx v8
   tissue (/expression/medianGeneExpression).
2. Ranks the 54 tissues by median TPM and keeps thyroid, pituitary and whole blood.
3. Joins the cis-eQTL record for rs6503875 from the deposited per-tissue table,
   so the effect size and the expression value come from the same release.

Usage
-----
    python 05_scripts/60_expression/40_gtex_axis_tissue.py
    python 05_scripts/60_expression/40_gtex_axis_tissue.py --gene SCRN2 --snp rs6503875

Needs network access for the API call (requests). The eQTL columns are read from
04_data/06_expression/SCRN2_rs6503875_GTEx_tissues.csv, which
39_gtex_cross_tissue_eqtl.py regenerates.
"""
import argparse
import os
import sys

import pandas as pd

BASE = "https://gtexportal.org/api/v2"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUTDIR = os.path.join(ROOT, "04_data", "06_expression")

KEEP = {"Thyroid": "Thyroid", "Pituitary": "Pituitary", "Whole_Blood": "Whole blood"}
ENSG = {"SCRN2": "ENSG00000141295.13"}


def get(session, path, **params):
    r = session.get(BASE + path, params=params, timeout=120)
    r.raise_for_status()
    return r.json()


def medians(session, gencode, dataset="gtex_v8"):
    js = get(session, "/expression/medianGeneExpression", gencodeId=gencode, datasetId=dataset)
    rows = js if isinstance(js, list) else js.get("data", [])
    return pd.DataFrame([(x["tissueSiteDetailId"], float(x["median"])) for x in rows],
                        columns=["tissue_id", "median_tpm"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default="SCRN2")
    ap.add_argument("--snp", default="rs6503875")
    ap.add_argument("--out", default=os.path.join(OUTDIR, "SCRN2_axis_expression.csv"))
    a = ap.parse_args()

    import requests

    gene = a.gene
    if gene not in ENSG:
        sys.exit("no Ensembl id on file for %s" % gene)

    s = requests.Session()
    m = medians(s, ENSG[gene])
    if m.empty:
        sys.exit("the API returned no median expression for %s" % gene)
    m = m.sort_values("median_tpm", ascending=False).reset_index(drop=True)
    m["rank"] = m.index + 1
    n = len(m)

    eqtl = pd.read_csv(os.path.join(OUTDIR, "SCRN2_rs6503875_GTEx_tissues.csv"))
    eqtl = eqtl[eqtl.snp == a.snp].set_index("tissue")

    rows = []
    for tid, label in KEEP.items():
        hit = m[m.tissue_id == tid]
        if hit.empty:
            continue
        r = hit.iloc[0]
        nes = p = ""
        if tid in eqtl.index:
            nes = round(float(eqtl.loc[tid, "nes"]), 4)
            p = float(eqtl.loc[tid, "p"])
        rows.append(dict(gene=gene, snp=a.snp, tissue=label, median_tpm=round(float(r.median_tpm), 2),
                         rank=int(r["rank"]), n_tissues=n, nes=nes, p=p))

    out = pd.DataFrame(rows)
    out.to_csv(a.out, index=False, encoding="utf-8-sig")
    pd.set_option("display.width", 160)
    print(out.to_string(index=False))
    print("\ntissues ranked: %d; pituitary percentile %.0f" %
          (n, 100 * (1 - (int(out.loc[out.tissue == "Pituitary", "rank"].iloc[0]) - 1) / n)))
    print("written", a.out)


if __name__ == "__main__":
    main()
