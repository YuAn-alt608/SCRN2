#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
39_gtex_cross_tissue_eqtl.py - cross-tissue cis-eQTL lookup for one variant
===========================================================================
Reproduces 04_data/06_expression/SCRN2_rs6503875_GTEx_tissues.csv, the table
behind section 2.9 of the manuscript ("significant in all 47 GTEx v8 tissues in
which it was testable") and Figure 8B.

What it does
------------
1. Asks the GTEx Portal API which v8 tissues carry cis-eQTL data
   (/dataset/tissueSiteDetail?datasetId=gtex_v8, field hasEGenes).
2. Queries /association/singleTissueEqtl for the variant and gene, paging
   1000 records at a time.
3. Writes one row per tissue: normalised effect size (NES) and P.

Why 47 and not 49
-----------------
GTEx v8 releases cis-eQTL data for 49 tissues. The API returns a single-tissue
record for rs6503875/SCRN2 in 47 of them; Brain_Hippocampus and
Brain_Hypothalamus carry no record for this variant. The manuscript states both
numbers (49 with data, 47 testable) and uses 47 as the denominator.

Usage
-----
    python 05_scripts/60_expression/39_gtex_cross_tissue_eqtl.py --snp rs6503875 --gene SCRN2
    python 05_scripts/60_expression/39_gtex_cross_tissue_eqtl.py --from-json cached.json

Needs network access for the API mode (requests). Output is written next to the
script unless --out is given.
"""
import argparse
import json
import os
import sys

BASE = "https://gtexportal.org/api/v2"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUTDIR = os.path.join(ROOT, "04_data", "06_expression")


def get(session, path, **params):
    r = session.get(BASE + path, params=params, timeout=120)
    r.raise_for_status()
    return r.json()


def eqtl_tissues(session, dataset="gtex_v8"):
    js = get(session, "/dataset/tissueSiteDetail", datasetId=dataset)
    rows = js if isinstance(js, list) else js.get("data", [])
    return sum(1 for x in rows if str(x.get("hasEGenes", "")).lower() in ("true", "1"))


def fetch(session, snp, gene, dataset="gtex_v8", page=1000):
    out, start = [], 0
    while True:
        js = get(session, "/association/singleTissueEqtl", geneSymbol=gene,
                 snpId=snp, datasetId=dataset, itemsPerPage=page, page=start)
        data = js.get("data", [])
        out.extend(data)
        if len(data) < page:
            break
        start += 1
    return out


def summarise(records, gene, snp):
    import pandas as pd
    df = pd.DataFrame([{"gene": gene, "snp": r["snpId"], "tissue": r["tissueSiteDetailId"],
                        "nes": round(r["nes"], 6), "p": r["pValue"]} for r in records])
    df = df.drop_duplicates("tissue").sort_values("nes").reset_index(drop=True)
    th = df.loc[df.tissue == "Thyroid", "nes"]
    others = df.loc[df.tissue != "Thyroid", "nes"]
    print("tissues returned          : %d" % len(df))
    print("significant (P < 0.05)    : %d/%d" % (int((df.p < 0.05).sum()), len(df)))
    print("same sign as thyroid      : %d/%d" % (int((df.nes < 0).sum()) - 1, len(df) - 1))
    if len(th):
        print("thyroid NES               : %.4f" % th.iloc[0])
        print("median NES, other tissues : %.4f" % others.median())
        print("percentile of thyroid NES : %.1f" % (100.0 * (others < th.iloc[0]).mean()))
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snp", default="rs6503875")
    ap.add_argument("--gene", default="SCRN2")
    ap.add_argument("--from-json", help="reuse a saved API dump instead of querying")
    ap.add_argument("--out", default=os.path.join(OUTDIR, "SCRN2_rs6503875_GTEx_tissues.csv"))
    a = ap.parse_args()

    if a.from_json:
        records = json.load(open(a.from_json))
    else:
        try:
            import requests
        except ImportError:
            sys.exit("requests is required for the API mode; use --from-json instead")
        s = requests.Session()
        print("tissues with cis-eQTL data: %d" % eqtl_tissues(s))
        records = fetch(s, a.snp, a.gene)

    df = summarise(records, a.gene, a.snp)
    df.to_csv(a.out, index=False)
    print("[written] %s" % a.out)


if __name__ == "__main__":
    main()
