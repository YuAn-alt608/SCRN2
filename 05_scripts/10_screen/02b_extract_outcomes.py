# -*- coding: utf-8 -*-
"""
Step 02b: extract the instrument loci from the full outcome GWAS files
Read in a stream so that 10M+ rows are never held in memory.
Output: data/out_<outcome>.tsv
"""
import gzip, os, csv, sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, DATA = os.path.join(BASE, "raw"), os.path.join(BASE, "data")

ivs = set()
with open(os.path.join(DATA, "ivs.tsv"), encoding="utf-8") as f:
    next(f)
    for line in f:
        ivs.add(line.split("\t")[1])
print(f"target SNPs: {len(ivs)}")

FILES = {
    "TSH":   ("TSH.tsv.gz",   True),
    "T4":    ("T4.tsv.gz",    True),
    "HYPO":  ("HYPO.tsv.gz",  True),
    "LT":    ("LT.tsv",       False),
}

ALIAS = {
    "rs_id": "rsid", "rsid": "rsid", "variant_id": "rsid", "SNP": "rsid", "snpid": "rsid",
    "effect_allele": "ea", "EA": "ea", "allele1": "ea", "alt": "ea",
    "other_allele": "nea", "non_effect_allele": "nea", "NEA": "nea", "allele2": "nea", "ref": "nea",
    "beta": "beta", "standard_error": "se", "se": "se",
    "p_value": "pval", "pval": "pval", "p": "pval",
    "effect_allele_frequency": "eaf", "eaf": "eaf", "AF": "eaf",
    "chromosome": "chr", "chr": "chr",
    "base_pair_location": "pos", "pos": "pos", "BP": "pos",
    "n": "n", "N": "n", "total_sample_size": "n",
}

ONLY = set(sys.argv[1:])
for name, (fn, gz) in FILES.items():
    if ONLY and name not in ONLY:
        continue
    path = os.path.join(RAW, fn)
    if not os.path.exists(path):
        print(f"[skip] {name} not downloaded")
        continue
    op = gzip.open(path, "rt", encoding="utf-8", errors="ignore") if gz else open(path, encoding="utf-8", errors="ignore")
    with op as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {}
        for i, h in enumerate(header):
            key = ALIAS.get(h.strip())
            if key and key not in col:
                col[key] = i
        need = ["rsid", "ea", "nea", "beta", "se", "pval"]
        if any(k not in col for k in need):
            print(f"[skip] {name} missing columns: {header}")
            continue

        out_path = os.path.join(DATA, f"out_{name}.tsv")
        n_hit = 0
        with open(out_path, "w", encoding="utf-8") as o:
            o.write("rsid\tea\tnea\tbeta\tse\tpval\teaf\tn\tchr\tpos\n")
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) <= max(col[k] for k in need):
                    continue
                rs = c[col["rsid"]]
                if rs not in ivs:
                    continue
                def g(k, d=""):
                    return c[col[k]] if k in col and col[k] < len(c) else d
                o.write("\t".join([rs, g("ea"), g("nea"), g("beta"), g("se"),
                                   g("pval"), g("eaf"), g("n"), g("chr"), g("pos")]) + "\n")
                n_hit += 1
    print(f"[ok] {name}: matched {n_hit} / {len(ivs)} SNPs  -> {out_path}")
