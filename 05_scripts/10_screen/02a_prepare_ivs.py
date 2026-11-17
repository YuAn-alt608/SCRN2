# -*- coding: utf-8 -*-
"""
Step 02a: select instruments from the significant cis-eQTLs in eQTLGen
  - restricted to the Module A candidate genes
  - P < 5e-8
  - F = Z^2 > 10
  - greedy de-duplication within a gene at 1 Mb physical distance (approximate LD independence)
Output: data/ivs.tsv, data/iv_summary.csv
"""
import gzip, os, csv
from collections import defaultdict

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, DATA = os.path.join(BASE, "raw"), os.path.join(BASE, "data")

P_THR, F_MIN, WINDOW = 5e-8, 10.0, 1_000_000
N_EQTL = 31684

with open(os.path.join(DATA, "edc_thyroid_genes.txt"), encoding="utf-8") as f:
    cand = {l.strip() for l in f if l.strip()}
print(f"candidate genes: {len(cand)}")

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

print(f"SNP-gene pairs passing P<5e-8: {len(rows):,}")

# within-gene distance pruning: sort by P and keep SNPs more than 1 Mb from those already taken
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

print(f"SNP-gene pairs after pruning: {len(kept):,}; genes: {len(set(k[0] for k in kept)):,}")

with open(os.path.join(DATA, "ivs.tsv"), "w", encoding="utf-8") as f:
    f.write("gene\tsnp\tchr\tpos\tea\tnea\tz\tp\tf_stat\n")
    for r in kept:
        f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]:.6f}\t{r[7]:.6e}\t{r[6]**2:.2f}\n")

# summary
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

print(f"genes with one SNP: {sum(1 for v in n_by_gene.values() if v==1)}")
print(f"genes with >=2 SNPs: {sum(1 for v in n_by_gene.values() if v>=2)}")
print(f"median F: {statistics.median([r[6]**2 for r in kept]):.1f}")
