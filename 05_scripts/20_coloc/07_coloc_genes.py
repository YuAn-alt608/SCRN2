# -*- coding: utf-8 -*-
"""
Step 07: choose target genes for colocalization and define the analysis windows

Three-tier gene selection logic:
  tier1 = genes consistent across outcomes (replicated across outcomes, strongest evidence)
  tier2 = top-N significant genes per outcome (ranked by p)
  tier3 = canonical thyroid positive controls (run regardless of significance, for methodological validation)

Output: data/coloc_genes.tsv
"""
import csv, os, collections

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")

TOP_N = 25
PATTERN = ["TPO", "TG", "TSHR", "PDE8B", "FOXE1", "SLC26A4", "DIO1", "DIO2",
           "HLA-DRB1", "CD40", "CTLA4", "PTPN22", "TGFB1", "NRXN3", "VAV3",
           "ITM2A", "CAPZB", "AADAT", "GLIS3", "BABAM2"]

rows = list(csv.DictReader(open(os.path.join(DATA, "results_placeholder")
                              if False else os.path.join(BASE, "results", "MR_annotated.csv"),
                              encoding="utf-8-sig")))
sig = [r for r in rows if float(r["fdr"]) < 0.05]

by_gene = collections.defaultdict(dict)
for r in sig:
    by_gene[r["gene"]][r["outcome"]] = r

info = {}
for g, d in by_gene.items():
    best = min(d.values(), key=lambda x: float(x["pval"]))
    info[g] = {
        "gene": g,
        "n_outcome": len(d),
        "outcomes": "+".join(sorted(d)),
        "best_pval": float(best["pval"]),
        "min_pval": min(float(x["pval"]) for x in d.values()),
        "chr": best["chr"], "pos": best["pos"], "snp": best["snp"],
        "tier": "",
    }

# tier1: multiple outcomes
for g in info:
    if info[g]["n_outcome"] >= 2:
        info[g]["tier"] = "T1_multi"

# tier2: top genes per outcome
for oc in ["TSH", "Free T4", "Hypothyroidism"]:
    sub = sorted([r for r in sig if r["outcome"] == oc], key=lambda x: float(x["pval"]))[:TOP_N]
    for r in sub:
        g = r["gene"]
        if g in info and not info[g]["tier"]:
            info[g]["tier"] = f"T2_top_{oc}"

# tier3: positive controls (keep an existing tier if present, and tag as control)
for g in PATTERN:
    if g in info:
        if not info[g]["tier"]:
            info[g]["tier"] = "Tier3_control"
        info[g]["tier"] = info[g]["tier"] + "|control"

# keep only the selected genes
sel = [v for v in info.values() if v["tier"]]
sel.sort(key=lambda x: (x["tier"], x["min_pval"]))

out = os.path.join(DATA, "coloc_genes.tsv")
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene", "tier", "n_outcome", "outcomes", "min_pval", "chr", "pos", "lead_snp", "window_kb"])
    for v in sel:
        w.writerow([v["gene"], v["tier"], v["n_outcome"], v["outcomes"],
                    f"{v['min_pval']:.3e}", v["chr"], v["pos"], v["snp"], 500])

print(f"{len(sel)} target genes")
c = collections.Counter(v["tier"] for v in sel)
for k, v in sorted(c.items()):
    print(f"  {k:<26} {v}")
print(f"-> {out}")
