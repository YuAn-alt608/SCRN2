# -*- coding: utf-8 -*-
"""
Step 01b: build the EDC-responsive gene set (Module A of the design)
Logic follows the source paper:
  1. the EDC chemical library from 01a
  2. experimentally validated human chemical-gene interactions in CTD
  3. degree centrality = number of EDC chemicals targeting the gene
  4. genes with degree >= 2 x median count as high-response genes
  5. intersect with genes expressed in GTEx thyroid -> the testable thyroid gene set
Output:
  data/gene_degree_all.csv        degree of every gene
  data/edc_thyroid_genes.txt      the genes finally taken forward
  data/moduleA_stats.txt
"""
import gzip, csv, os, statistics
from collections import defaultdict

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, DATA = os.path.join(BASE, "raw"), os.path.join(BASE, "data")

# 1. EDC chemicals
edc = set()
with open(os.path.join(DATA, "edc_chemicals.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        edc.add(r["chemical_id"].replace("MESH:","").replace("MESH",""))
print(f"[1] EDC chemicals: {len(edc)}")

# 2. thyroid-expressed genes
with open(os.path.join(DATA, "thyroid_expressed.txt"), encoding="utf-8") as f:
    thyroid = {l.strip() for l in f if l.strip()}
print(f"[2] GTEx thyroid genes with TPM>=1: {len(thyroid)}")

# 3. scan the CTD interactions and count how many distinct EDCs target each gene
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

print(f"[3] EDC-gene interaction records: {total_pairs:,}; genes involved: {len(chems_of_gene):,}")

# 4. degree-centrality threshold
degrees = {g: len(s) for g, s in chems_of_gene.items()}
med = statistics.median(degrees.values())
thr = 2 * med
high = {g: d for g, d in degrees.items() if d >= thr}
print(f"[4] median degree {med:.1f}, threshold {thr:.1f}, high-response genes {len(high):,}")

# 5. intersect with thyroid expression
final = sorted(g for g in high if g in thyroid)
print(f"[5] intersect thyroid-expressed genes = {len(final):,}")

# output
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

print("\ntop 20 high-response genes:")
for g in sorted(degrees, key=lambda x: -degrees[x])[:20]:
    print(f"  {g:<12} EDC={degrees[g]:>5}  interactions={ixn_count[g]:>6}  thyroid={g in thyroid}")
