# -*- coding: utf-8 -*-
"""
Step 01b: 构建 EDC 响应基因集（对应方案 Module A）
逻辑复刻源论文：
  1. EDC 化学物库（01a 产出）
  2. CTD 实验验证的人源化学物–基因互作
  3. 以「靶向该基因的 EDC 化学物数量」为度中心性
  4. 度值 >= 2 x 中位数 者为高响应基因
  5. 与 GTEx 甲状腺表达基因取交集 → 甲状腺可检验基因集
输出:
  data/gene_degree_all.csv        全部基因的度值
  data/edc_thyroid_genes.txt      最终待检验基因
  data/moduleA_stats.txt
"""
import gzip, csv, os, statistics
from collections import defaultdict

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, DATA = os.path.join(BASE, "raw"), os.path.join(BASE, "data")

# 1. EDC 化学物
edc = set()
with open(os.path.join(DATA, "edc_chemicals.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        edc.add(r["chemical_id"].replace("MESH:","").replace("MESH",""))
print(f"[1] EDC 化学物: {len(edc)}")

# 2. 甲状腺表达基因
with open(os.path.join(DATA, "thyroid_expressed.txt"), encoding="utf-8") as f:
    thyroid = {l.strip() for l in f if l.strip()}
print(f"[2] GTEx 甲状腺 TPM>=1 基因: {len(thyroid)}")

# 3. 扫描 CTD 互作，统计每个基因被多少种 EDC 靶向
chems_of_gene = defaultdict(set)
ixn_count = defaultdict(int)
total_pairs = 0
with gzip.open(os.path.join(RAW, "CTD_chem_gene_ixns.tsv.gz"), "rt",
               encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 8:
            continue
        cid, gene, organism = c[1].strip(), c[3].strip(), c[6].strip()
        if organism != "Homo sapiens":
            continue
        if cid not in edc:
            continue
        total_pairs += 1
        chems_of_gene[gene].add(cid)
        ixn_count[gene] += 1

print(f"[3] EDC–基因互作记录: {total_pairs:,}；涉及基因: {len(chems_of_gene):,}")

# 4. 度中心性阈值
degrees = {g: len(s) for g, s in chems_of_gene.items()}
med = statistics.median(degrees.values())
thr = 2 * med
high = {g: d for g, d in degrees.items() if d >= thr}
print(f"[4] 度值中位数 {med:.1f}，阈值 {thr:.1f}，高响应基因 {len(high):,}")

# 5. 与甲状腺表达取交集
final = sorted(g for g in high if g in thyroid)
print(f"[5] ∩ 甲状腺表达基因 = {len(final):,}")

# 输出
with open(os.path.join(DATA, "gene_degree_all.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["gene", "n_edc_chemicals", "n_interactions", "thyroid_expressed", "high_responder"])
    for g in sorted(degrees, key=lambda x: -degrees[x]):
        w.writerow([g, degrees[g], ixn_count[g], int(g in thyroid), int(g in high)])

with open(os.path.join(DATA, "edc_thyroid_genes.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(final))

with open(os.path.join(DATA, "moduleA_stats.txt"), "w", encoding="utf-8") as f:
    f.write(f"EDC_chemicals\t{len(edc)}\n")
    f.write(f"thyroid_expressed_genes\t{len(thyroid)}\n")
    f.write(f"edc_gene_pairs\t{total_pairs}\n")
    f.write(f"genes_with_any_edc\t{len(degrees)}\n")
    f.write(f"degree_median\t{med}\n")
    f.write(f"degree_threshold\t{thr}\n")
    f.write(f"high_responder_genes\t{len(high)}\n")
    f.write(f"thyroid_candidate_genes\t{len(final)}\n")

print("\n前 20 位高响应基因:")
for g in sorted(degrees, key=lambda x: -degrees[x])[:20]:
    print(f"  {g:<12} EDC={degrees[g]:>5}  interactions={ixn_count[g]:>6}  thyroid={g in thyroid}")
