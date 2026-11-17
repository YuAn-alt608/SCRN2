# -*- coding: utf-8 -*-
"""
Step 10: MR 与共定位的联合证据分级

单一 MR 只能说明"遗传预测的表达与结局相关"，无法排除 LD 混杂；
共定位 PP.H4 高才说明两者共享同一个因果变异。两者结合才能定级。

分级规则（对每个 基因×结局 组合）：
  A   MR FDR<0.05 + coloc PP.H4>0.8 + 该基因在≥2 个结局显著   → 最强
  B   MR FDR<0.05 + coloc PP.H4>0.8                            → 强
  C   MR FDR<0.05 + coloc PP.H3>0.8（不同因果变异）             → MR 结果需谨慎
  D   MR FDR<0.05 + coloc H4<0.5 且 H3<0.5（无信息）            → 未获支持
  E   MR FDR<0.05 + 区域 SNP 不足，未做 coloc                   → 不可评估
输出: results/evidence_tiers.csv
"""
import csv, os
from collections import defaultdict

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")

mr = list(csv.DictReader(open(os.path.join(RES, "MR_annotated.csv"), encoding="utf-8-sig")))
sig = [r for r in mr if r["fdr"] not in ("", "nan") and float(r["fdr"]) < 0.05]

coloc = {}
if os.path.exists(os.path.join(RES, "coloc_results.csv")):
    for r in csv.DictReader(open(os.path.join(RES, "coloc_results.csv"), encoding="utf-8-sig")):
        coloc[(r["gene"], r["outcome"])] = r

# 每个基因覆盖的结局数
n_out = defaultdict(set)
for r in sig:
    n_out[r["gene"]].add(r["outcome"])

LABEL = {
    "TSH": "TSH", "Total T4": "Total T4",
    "Hypothyroidism": "Hypothyroidism", "Hashimoto thyroiditis": "Hashimoto",
}
# coloc 里结局用短码
CODE = {"TSH": "TSH", "Total T4": "T4", "Hypothyroidism": "HYPO", "Hashimoto thyroiditis": "HASH"}

rows = []
for r in sig:
    g, oc = r["gene"], r["outcome"]
    c = coloc.get((g, CODE[oc]))
    nsig = len(n_out[g])
    if c is None:
        tier, h4, h3, nsnp, topsnp = "E_unassessable", "", "", "", ""
    else:
        h4, h3 = float(c["PP.H4"]), float(c["PP.H3"])
        nsnp, topsnp = c["nsnps"], c["top_snp"]
        if h4 > 0.8:
            tier = "A_multi" if nsig >= 2 else "B_strong"
        elif h3 > 0.8:
            tier = "C_distinct"
        else:
            tier = "D_unsupported"
    rows.append({
        "gene": g, "outcome": oc, "b": r["b"], "pval": r["pval"], "fdr": r["fdr"],
        "nsnp_iv": r["nsnp"], "f_stat": r["f_stat"], "chr": r["chr"], "pos": r["pos"],
        "snp": r["snp"], "locus": r["locus"], "n_edc_chem": r["n_edc_chem"],
        "edc_classes": r["edc_classes"],
        "n_outcome_sig": nsig,
        "PP.H4": h4, "PP.H3": h3, "coloc_nsnps": nsnp, "coloc_top_snp": topsnp,
        "tier": tier,
    })

order = {"A_multi": 0, "B_strong": 1, "C_distinct": 2, "D_unsupported": 3, "E_unassessable": 4}
rows.sort(key=lambda x: (order[x["tier"]], -x["n_outcome_sig"], float(x["pval"])))

out = os.path.join(RES, "evidence_tiers.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

from collections import Counter
print(f"MR 显著记录 {len(rows)} 条")
c = Counter(r["tier"] for r in rows)
for k in ["A_multi", "B_strong", "C_distinct", "D_unsupported", "E_unassessable"]:
    print(f"  {k:<16} {c.get(k, 0):>4}")

print("\n== A 级（多结局 + 强共定位）==")
for r in [x for x in rows if x["tier"] == "A_multi"][:30]:
    print(f"  {r['gene']:<10} {LABEL[r['outcome']]:<16} b={float(r['b']):>7.3f} "
          f"p={float(r['pval']):.1e} H4={float(r['PP.H4']):.2f} "
          f"edc={r['n_edc_chem']:<3} {r['edc_classes'][:40]}")

print("\n== B 级 top10 ==")
for r in [x for x in rows if x["tier"] == "B_strong"][:10]:
    print(f"  {r['gene']:<10} {LABEL[r['outcome']]:<16} b={float(r['b']):>7.3f} "
          f"p={float(r['pval']):.1e} H4={float(r['PP.H4']):.2f}")

# 分结局统计
print("\n== 各结局的共定位支持率 ==")
agg = defaultdict(lambda: [0, 0])
for r in rows:
    a = agg[r["outcome"]]
    a[0] += 1
    if r["tier"] in ("A_multi", "B_strong"):
        a[1] += 1
for k, v in agg.items():
    print(f"  {k:<24} {v[1]}/{v[0]} ({100*v[1]/v[0]:.1f}%)")
