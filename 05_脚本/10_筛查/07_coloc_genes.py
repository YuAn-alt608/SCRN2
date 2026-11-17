# -*- coding: utf-8 -*-
"""
Step 07: 为共定位分析选定目标基因并生成区域定义

选基因的三层逻辑：
  tier1 = 多结局一致基因（跨结局重复，证据最强）
  tier2 = 各结局 top-N 显著基因（按 p 排序）
  tier3 = 甲状腺经典阳性对照（无论是否显著都做，用于方法学验证）

输出: data/coloc_genes.tsv
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

# tier1: 多结局
for g in info:
    if info[g]["n_outcome"] >= 2:
        info[g]["tier"] = "T1_multi"

# tier2: 各结局 top
for oc in ["TSH", "Total T4", "Hypothyroidism"]:
    sub = sorted([r for r in sig if r["outcome"] == oc], key=lambda x: float(x["pval"]))[:TOP_N]
    for r in sub:
        g = r["gene"]
        if g in info and not info[g]["tier"]:
            info[g]["tier"] = f"T2_top_{oc}"

# tier3: 阳性对照（若曾显著则保留已有 tier，另标 control）
for g in PATTERN:
    if g in info:
        if not info[g]["tier"]:
            info[g]["tier"] = "T3_control"
        info[g]["tier"] = info[g]["tier"] + "|control"

# 只保留被选中的
sel = [v for v in info.values() if v["tier"]]
sel.sort(key=lambda x: (x["tier"], x["min_pval"]))

out = os.path.join(DATA, "coloc_genes.tsv")
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene", "tier", "n_outcome", "outcomes", "min_pval", "chr", "pos", "lead_snp", "window_kb"])
    for v in sel:
        w.writerow([v["gene"], v["tier"], v["n_outcome"], v["outcomes"],
                    f"{v['min_pval']:.3e}", v["chr"], v["pos"], v["snp"], 500])

print(f"目标基因 {len(sel)} 个")
c = collections.Counter(v["tier"] for v in sel)
for k, v in sorted(c.items()):
    print(f"  {k:<26} {v}")
print(f"-> {out}")
