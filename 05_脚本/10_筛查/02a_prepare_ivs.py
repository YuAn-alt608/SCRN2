# -*- coding: utf-8 -*-
"""
Step 02a: 从 eQTLGen 显著 cis-eQTL 中挑选工具变量
  - 限制在 Module A 的候选基因
  - P < 5e-8
  - F = Z^2 > 10
  - 基因内按物理距离 1 Mb 贪心去重（近似 LD 独立）
输出: data/ivs.tsv, data/iv_summary.csv
"""
import gzip, os, csv
from collections import defaultdict

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, DATA = os.path.join(BASE, "raw"), os.path.join(BASE, "data")

P_THR, F_MIN, WINDOW = 5e-8, 10.0, 1_000_000
N_EQTL = 31684

with open(os.path.join(DATA, "edc_thyroid_genes.txt"), encoding="utf-8") as f:
    cand = {l.strip() for l in f if l.strip()}
print(f"候选基因: {len(cand)}")

rows = []
with gzip.open(os.path.join(RAW, "eQTLGen_cis_FDR05.txt.gz"), "rt",
               encoding="utf-8", errors="ignore") as f:
    header = f.readline().rstrip("\n").split("\t")
    H = {k: i for i, k in enumerate(header)}
    for line in f:
        c = line.rstrip("\n").split("\t")
        g = c[H["GeneSymbol"]]
        if g not in cand:
            continue
        p = float(c[H["Pvalue"]])
        if p >= P_THR:
            continue
        z = float(c[H["Zscore"]])
        if z * z < F_MIN:
            continue
        rows.append((g, c[H["SNP"]], c[H["SNPChr"]], int(c[H["SNPPos"]]),
                     c[H["AssessedAllele"]], c[H["OtherAllele"]], z, p))

print(f"通过 P<5e-8 的 SNP-基因对: {len(rows):,}")

# 基因内距离去重：按 P 升序，保留与已选 SNP 距离 > 1Mb 者
by_gene = defaultdict(list)
for r in rows:
    by_gene[r[0]].append(r)

kept = []
for g, lst in by_gene.items():
    lst.sort(key=lambda x: x[7])
    chosen = []
    for r in lst:
        if all(abs(r[3] - s[3]) > WINDOW for s in chosen):
            chosen.append(r)
    kept.extend(chosen)

print(f"去重后 SNP-基因对: {len(kept):,}；基因数: {len(set(k[0] for k in kept)):,}")

with open(os.path.join(DATA, "ivs.tsv"), "w", encoding="utf-8") as f:
    f.write("gene\tsnp\tchr\tpos\tea\tnea\tz\tp\tf_stat\n")
    for r in kept:
        f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]:.6f}\t{r[7]:.6e}\t{r[6]**2:.2f}\n")

# 汇总
import statistics
n_by_gene = defaultdict(int)
f_by_gene = defaultdict(float)
for r in kept:
    n_by_gene[r[0]] += 1
    f_by_gene[r[0]] += r[6] ** 2

with open(os.path.join(DATA, "iv_summary.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["gene", "n_snp", "max_f", "sum_f"])
    for g in sorted(n_by_gene):
        w.writerow([g, n_by_gene[g], round(max(r[6] ** 2 for r in kept if r[0] == g), 1),
                    round(f_by_gene[g], 1)])

print(f"单 SNP 基因: {sum(1 for v in n_by_gene.values() if v==1)}")
print(f"≥2 SNP 基因: {sum(1 for v in n_by_gene.values() if v>=2)}")
print(f"F 中位数: {statistics.median([r[6]**2 for r in kept]):.1f}")
