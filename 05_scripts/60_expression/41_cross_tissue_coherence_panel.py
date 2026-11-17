#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
41_cross_tissue_coherence_panel.py - cross-tissue direction concordance for a gene panel
=======================================================================================
Reproduces 04_data/06_expression/cross_tissue_coherence_panel.csv, the table behind
layer 4 of the five-layer scheme (Figure 5) and its threshold sensitivity (Figure 5D).

What it does
------------
1. Rebuilds the candidate panel: genes that pass layer 1 (whole-blood cis-eQTL MR for
   TSH, FDR < 0.05) and layer 2 (thyroid cis-eQTL MR for TSH, FDR < 0.05) — 71 genes.
2. For each gene, uses the variant that instruments it in layer 2, i.e. the GTEx v8
   thyroid cis-eQTL that carries the thyroid MR estimate (rs6503875 for SCRN2). If that
   variant has no GTEx record, the whole-blood instruments are tried in order of
   decreasing |z|.
3. Pulls every tissue in which that variant is a cis-eQTL for the gene
   (/association/singleTissueEqtl) and counts, among tissues other than thyroid, the
   fraction whose normalised effect size has the same sign as the thyroid effect.

Layer 4 passes a gene when that fraction is >= 99% and at least MIN_TISSUES tissues
carry a significant record; genes below MIN_TISSUES are reported as not evaluable
rather than as passing on a two-tissue denominator.

What layer 4 does not do
------------------------
Direction concordance is close to universal for a cis-eQTL that is significant in many
tissues: most of the panel passes. Layer 4 is therefore a gate (it excludes a
tissue-idiosyncratic instrument) rather than a discriminator, and the manuscript says
so. The narrowing in Figure 5 comes from layers 1-3 and 5.

Usage
-----
    python 05_scripts/60_expression/41_cross_tissue_coherence_panel.py            # query the API
    python 05_scripts/60_expression/41_cross_tissue_coherence_panel.py --from-cache

Responses are cached under _rp/l4_cache/ so a rerun is offline after the first pass.
"""
import argparse
import json
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

BASE = "https://gtexportal.org/api/v2"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CACHE = os.path.join(ROOT, "_rp", "l4_cache")
OUT = os.path.join(ROOT, "04_data", "06_expression", "cross_tissue_coherence_panel.csv")

MIN_TISSUES = 10      # tissues with a significant record needed to call a gene
CONCORD = 0.99        # layer-4 threshold
MAX_INSTRUMENTS = 4   # fallback instruments tried per gene
WORKERS = 6

_local = threading.local()


def sess():
    import requests
    if not hasattr(_local, "s"):
        s = requests.Session()
        s.headers.update({"User-Agent": "SCRN2-analysis/1.0"})
        _local.s = s
    return _local.s


def get(path, **params):
    s = sess()
    for attempt in range(4):
        r = s.get(BASE + path, params=params, timeout=180)
        if r.status_code == 200:
            return r.json()
        if r.status_code in (429, 500, 502, 503):
            time.sleep(3 * (attempt + 1))
            continue
        r.raise_for_status()
    r.raise_for_status()


def panel():
    """candidate genes, with the layer-2 instrument first and whole-blood fallbacks after."""
    mr = pd.read_csv(os.path.join(ROOT, "04_data/01_screen/MR_significant_FDR05.csv"))
    thy = pd.read_csv(os.path.join(ROOT, "04_data/08_plane/gw_thyroid_mr.csv"))
    l1 = set(mr.loc[mr.outcome == "TSH", "gene"])
    sig2 = thy[(thy.outcome == "TSH") & (thy.fdr < 0.05)]
    l2 = set(sig2.gene)
    thyroid_iv = dict(zip(sig2.gene, sig2.snp))

    iv = pd.read_csv(os.path.join(ROOT, "04_data/01_screen/ivs.tsv"), sep="\t")
    iv["absz"] = iv.z.abs()
    iv = iv.sort_values(["gene", "absz"], ascending=[True, False])
    blood = {g: list(d.snp)[:MAX_INSTRUMENTS] for g, d in iv.groupby("gene")}

    plan = {}
    for g in sorted(l1 & l2):
        order = [thyroid_iv[g]] + [s for s in blood.get(g, []) if s != thyroid_iv.get(g)]
        plan[g] = order[: 1 + MAX_INSTRUMENTS]
    return plan


def records_for(gene, rsid):
    cache = os.path.join(CACHE, "%s__%s.json" % (gene, rsid))
    if os.path.exists(cache):
        return json.load(open(cache))
    js = get("/dataset/variant", snpId=rsid, datasetId="gtex_v8")
    rows = js.get("data", [])
    recs = []
    if rows:
        vid = rows[0]["variantId"]
        i = 0
        while True:
            js2 = get("/association/singleTissueEqtl", variantId=vid, datasetId="gtex_v8",
                      itemsPerPage=1000, page=i)
            data = js2.get("data", [])
            recs.extend({"geneSymbol": r["geneSymbol"],
                         "tissueSiteDetailId": r["tissueSiteDetailId"],
                         "nes": r["nes"], "pValue": r["pValue"]} for r in data)
            if len(data) < 1000:
                break
            i += 1
    json.dump(recs, open(cache, "w"))
    time.sleep(0.2)
    return recs


def one(gene, order):
    for rsid in order:
        recs = records_for(gene, rsid)
        thy = [x for x in recs if x["geneSymbol"] == gene and x["tissueSiteDetailId"] == "Thyroid"]
        if thy:
            others = [x for x in recs if x["geneSymbol"] == gene
                      and x["tissueSiteDetailId"] != "Thyroid"]
            sign = 1 if thy[0]["nes"] > 0 else -1
            sig = [x for x in others if x["pValue"] < 0.05]
            same = [x for x in sig if (1 if x["nes"] > 0 else -1) == sign]
            pct = 100.0 * len(same) / len(sig) if sig else float("nan")
            ok = len(sig) >= MIN_TISSUES
            return dict(gene=gene, instrument=rsid, n_tissue=len(recs),
                        n_sig_p05=len(sig), n_same_sign=len(same),
                        pct_same_sign=round(pct, 1) if sig else float("nan"),
                        nes_thyroid=round(thy[0]["nes"], 5), evaluable=ok,
                        layer4=bool(ok and pct >= 100 * CONCORD))
    return dict(gene=gene, instrument=";".join(order), n_tissue=0, n_sig_p05=0,
                n_same_sign=0, pct_same_sign=float("nan"), nes_thyroid=float("nan"),
                evaluable=False, layer4=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-cache", action="store_true")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    os.makedirs(CACHE, exist_ok=True)

    plan = panel()
    with ThreadPoolExecutor(WORKERS) as ex:
        rows = list(ex.map(lambda kv: one(*kv), plan.items()))

    df = pd.DataFrame(rows).sort_values("gene")
    df.to_csv(a.out, index=False)
    for _, r in df.iterrows():
        print("%-12s %-12s sig %2d same %2d %6.1f%% %s"
              % (r.gene, r.instrument, r.n_sig_p05, r.n_same_sign,
                 r.pct_same_sign if pd.notna(r.pct_same_sign) else float("nan"),
                 "PASS" if r.layer4 else ("n/e" if not r.evaluable else "")))
    print("\n[evaluable] %d/%d" % (int(df.evaluable.sum()), len(df)))
    print("[layer 4]   %d genes" % int(df.layer4.sum()))
    print("[written]   %s" % a.out)


if __name__ == "__main__":
    main()
