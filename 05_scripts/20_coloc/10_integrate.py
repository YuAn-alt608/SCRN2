# -*- coding: utf-8 -*-
"""
Step 10: joint evidence grading from MR and colocalization

MR alone only shows that "genetically predicted expression correlates with the
outcome" and cannot exclude LD confounding; only a high colocalization PP.H4
shows the two traits share the same causal variant. Grading needs both.

Grading rules (per gene x outcome pair):
  A   MR FDR<0.05 + coloc PP.H4>0.8 + gene significant in >=2 outcomes  -> strongest
  B   MR FDR<0.05 + coloc PP.H4>0.8                                     -> strong
  C   MR FDR<0.05 + coloc PP.H3>0.8 (distinct causal variants)          -> MR result needs care
  D   MR FDR<0.05 + coloc H4<0.5 and H3<0.5 (uninformative)             -> not supported
  E   MR FDR<0.05 + too few regional SNPs, coloc not run                -> not assessable
Output: results/evidence_tiers.csv
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

# number of outcomes covered by each gene
n_out = defaultdict(set)
for r in sig:
    n_out[r["gene"]].add(r["outcome"])

LABEL = {
    "TSH": "TSH", "Free T4": "Free T4",
    "Hypothyroidism": "Hypothyroidism", "Hashimoto thyroiditis": "Hashimoto",
}
# coloc uses short outcome codes
CODE = {"TSH": "TSH", "Free T4": "T4", "Hypothyroidism": "HYPO", "Hashimoto thyroiditis": "HASH"}

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
print(f"{len(rows)} significant MR records")
c = Counter(r["tier"] for r in rows)
for k in ["A_multi", "B_strong", "C_distinct", "D_unsupported", "E_unassessable"]:
    print(f"  {k:<16} {c.get(k, 0):>4}")

print("\n== tier A (multiple outcomes + strong colocalization) ==")
for r in [x for x in rows if x["tier"] == "A_multi"][:30]:
    print(f"  {r['gene']:<10} {LABEL[r['outcome']]:<16} b={float(r['b']):>7.3f} "
          f"p={float(r['pval']):.1e} H4={float(r['PP.H4']):.2f} "
          f"edc={r['n_edc_chem']:<3} {r['edc_classes'][:40]}")

print("\n== tier B, top 10 ==")
for r in [x for x in rows if x["tier"] == "B_strong"][:10]:
    print(f"  {r['gene']:<10} {LABEL[r['outcome']]:<16} b={float(r['b']):>7.3f} "
          f"p={float(r['pval']):.1e} H4={float(r['PP.H4']):.2f}")

# per-outcome statistics
print("\n== colocalization support rate by outcome ==")
agg = defaultdict(lambda: [0, 0])
for r in rows:
    a = agg[r["outcome"]]
    a[0] += 1
    if r["tier"] in ("A_multi", "B_strong"):
        a[1] += 1
for k, v in agg.items():
    print(f"  {k:<24} {v[1]}/{v[0]} ({100*v[1]/v[0]:.1f}%)")
