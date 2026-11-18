#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
42_five_layer_sensitivity.py - rebuild the five evidence layers and relax their thresholds
=========================================================================================
Reproduces the layer sizes, the region counts of Figure 5A and every configuration
plotted in Figure 5D, from the tables deposited in 04_data/. Nothing here is entered
by hand: each layer is a filter on a deposited table, so the per-configuration gene
lists can be read off directly.

Layers (all for the TSH branch)
-------------------------------
  L1  whole-blood cis-eQTL MR     04_data/01_screen/MR_significant_FDR05.csv   outcome TSH
  L2  thyroid cis-eQTL MR         04_data/08_plane/gw_thyroid_mr.csv  outcome TSH
  L3  colocalisation              04_data/02_coloc/coloc_results.csv        outcome TSH
  L4  cross-tissue coherence      04_data/06_expression/cross_tissue_coherence_panel.csv
  L5  eQTL target assignment      04_data/05_locus/25_eqtl_target_ranking.csv

Outputs
-------
  04_data/09_layers/five_layer_membership.csv            per gene, per configuration
  04_data/09_layers/five_layer_threshold_sensitivity.csv per configuration, layer sizes
                                                       and the five-layer intersection

Usage
-----
    python 05_scripts/90_figures/42_five_layer_sensitivity.py
"""
import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
P = lambda *a: os.path.join(ROOT, *a)

SENTINELS = ["rs6503875", "rs62076102", "rs4794224"]   # the three instrument SNPs named in the text

CONFIGS = [
    ("baseline",     dict(h4=0.80, top=3, concord=1.00)),
    ("PP.H4>=0.70",  dict(h4=0.70, top=3, concord=1.00)),
    ("top-5 targets", dict(h4=0.80, top=5, concord=1.00)),
    ("concordance>=90%", dict(h4=0.80, top=3, concord=0.90)),
    ("all three relaxed", dict(h4=0.70, top=5, concord=0.90)),
]


def load():
    mr = pd.read_csv(P("04_data/01_screen/MR_significant_FDR05.csv"))
    thy = pd.read_csv(P("04_data/08_plane/gw_thyroid_mr.csv"))
    co = pd.read_csv(P("04_data/02_coloc/coloc_results.csv"))
    panel = pd.read_csv(P("04_data/06_expression/cross_tissue_coherence_panel.csv"))
    tgt = pd.read_csv(P("04_data/05_locus/25_eqtl_target_ranking.csv"))
    return mr, thy, co, panel, tgt


def layers(cfg, mr, thy, co, panel, tgt):
    l1 = set(mr.loc[mr.outcome == "TSH", "gene"])
    l2 = set(thy.loc[(thy.outcome == "TSH") & (thy.fdr < 0.05), "gene"])
    l3 = set(co.loc[(co.outcome == "TSH") & (co["PP.H4"] >= cfg["h4"]), "gene"])
    p = panel[panel.evaluable]
    l4 = set(p.loc[p.pct_same_sign >= 100 * cfg["concord"], "gene"])
    t = tgt[(tgt.snp.isin(SENTINELS)) & (tgt["rank"] <= cfg["top"])]
    l5 = set(t.gene)
    return [l1, l2, l3, l4, l5]


def main():
    mr, thy, co, panel, tgt = load()
    all_genes = sorted(set(mr.gene) | set(thy.gene) | set(co.gene) | set(panel.gene) | set(tgt.gene))

    print("layer sizes by configuration")
    print("%-22s %5s %5s %5s %5s %5s   %s" % ("configuration", "L1", "L2", "L3", "L4", "L5", "five-layer intersection"))
    rows, membership = [], {}
    for name, cfg in CONFIGS:
        L = layers(cfg, mr, thy, co, panel, tgt)
        five = set.intersection(*[set(x) for x in L])
        print("%-22s %5d %5d %5d %5d %5d   %s"
              % (name, len(L[0]), len(L[1]), len(L[2]), len(L[3]), len(L[4]),
                 ", ".join(sorted(five)) if five else "(empty)"))
        rows.append(dict(configuration=name, PP_H4_min=cfg["h4"], target_rank_max=cfg["top"],
                         concordance_min=cfg["concord"],
                         layer1_whole_blood_MR=len(L[0]), layer2_thyroid_MR=len(L[1]),
                         layer3_colocalisation=len(L[2]), layer4_cross_tissue=len(L[3]),
                         layer5_target=len(L[4]), n_five_layer=len(five),
                         five_layer_genes=";".join(sorted(five))))
        for g in all_genes:
            bit = sum(1 << i for i, x in enumerate(L) if g in x)
            membership.setdefault(g, {})[name] = bit

    # bitmask table for the baseline, comparable with table_region_counts.csv
    base = pd.DataFrame([(g, v["baseline"]) for g, v in membership.items()],
                        columns=["gene", "bitmask"])
    cnt = base[base.bitmask > 0].bitmask.value_counts().sort_index()
    print("\nbaseline regions (bitmask = L1 + 2*L2 + 4*L3 + 8*L4 + 16*L5)")
    for bit, n in cnt.items():
        genes = base.loc[base.bitmask == bit, "gene"].tolist()
        print("  %2d  n=%2d  %s" % (bit, n, ", ".join(sorted(genes)[:8]) + (" ..." if n > 8 else "")))

    # the cumulative chain quoted in Table 3A and Figure 5B
    L = layers(CONFIGS[0][1], mr, thy, co, panel, tgt)
    chain = [len(L[0]), len(L[0] & L[1]), len(L[0] & L[1] & L[2]),
             len(L[0] & L[1] & L[2] & L[3]), len(set.intersection(*[set(x) for x in L]))]
    print("\ncumulative chain (layers 1..5): %s" % " -> ".join(map(str, chain)))
    print("occupied regions: %d of 31" % len(cnt))

    # region counts in the deposited form
    out = pd.DataFrame({"bitmask": cnt.index.astype(int), "n_genes": cnt.values})
    out["genes"] = [",".join(sorted(base.loc[base.bitmask == b, "gene"])) for b in out.bitmask]
    out.to_csv(P("04_data/09_layers/table_region_counts.csv"), index=False)
    json.dump({str(b): {"n": int(n), "genes": sorted(base.loc[base.bitmask == b, "gene"])}
               for b, n in cnt.items()},
              open(P("04_data/09_layers/region_counts.json"), "w"), indent=1)

    mem = pd.DataFrame([dict(gene=g, **d) for g, d in membership.items()],
                       columns=["gene"] + [c[0] for c in CONFIGS])
    mem = mem.sort_values(["baseline", "gene"], ascending=[False, True])
    mem.to_csv(P("04_data/09_layers/five_layer_membership.csv"), index=False)
    pd.DataFrame(rows).to_csv(P("04_data/09_layers/five_layer_threshold_sensitivity.csv"), index=False)
    print("\n[written] 04_data/09_layers/five_layer_membership.csv")
    print("[written] 04_data/09_layers/five_layer_threshold_sensitivity.csv")

    # the numbers the manuscript quotes
    b = rows[0]
    assert b["layer1_whole_blood_MR"] == 337, b
    assert b["layer2_thyroid_MR"] == 180, b
    assert b["layer3_colocalisation"] == 26, b
    assert b["layer4_cross_tissue"] == 36, b
    assert b["layer5_target"] == 4, b
    assert b["five_layer_genes"] == "SCRN2", b
    assert all(r["five_layer_genes"] == "SCRN2" for r in rows), rows
    print("\n[check] baseline 337/180/26/36/4, SCRN2 sole member in every configuration")


if __name__ == "__main__":
    main()
