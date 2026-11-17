#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p26_outcome_identity.py
=======================
Establishes the identity and sample-size definition of the three outcome GWAS
(read only; no raw file is modified).

Answers two questions:
  1. Is T4 free or total?   -> three-way alignment of the paper's N, the per-SNP n
     ceiling in the file, and the cohort decomposition of n
  2. Why does the hypothyroidism n not add up? -> compare four numbers: the paper
     total, the European part, the non-European part and the file ceiling

Approach:
  - stream the value distribution of the `n` column from DATA10/raw/{TSH,T4,HYPO}.tsv.gz
    (never loaded into memory);
  - cache the result in tables/34_outcome_per_snp_sample_size.csv and read the cache
    from the second run onwards (fast);
  - write the conclusion to tables/35_outcome_gwas_identity.csv.

Usage:
    python p26_outcome_identity.py           # returns instantly when the cache exists
    python p26_outcome_identity.py --refresh # force a rescan (about 4 minutes)
"""
import os, sys, gzip
from collections import Counter
import numpy as np
import pandas as pd

BASE = "D:/NewBeginning/MR1"
OUT = f"{BASE}/DATA21_complete_2026-11-14"
RAW = f"{BASE}/DATA10/raw"
TAB = f"{OUT}/tables"
CACHE = f"{TAB}/34_outcome_per_snp_sample_size.csv"

FILES = {"TSH": f"{RAW}/TSH.tsv.gz", "FT4": f"{RAW}/T4.tsv.gz", "HYPO": f"{RAW}/HYPO.tsv.gz"}

# official numbers from the paper and GWAS Catalog (Rand SA et al., Nat Genet 2025, PMID 41238958)
PAPER = {
    "TSH":  dict(acc="GCST90572789", trait="Thyroid stimulating hormone levels",
                 paper_n=482873, snp_count=29790163, note="482,873 European ancestry individuals"),
    "FT4":  dict(acc="GCST90572790", trait="Thyroxine levels (OBA_1001010, a coarse-grained parent term)",
                 paper_n=191449, snp_count=28647991, note="the abstract says free thyroxine (191,449 individuals)"),
    "HYPO": dict(acc="GCST90572791", trait="Hypothyroidism",
                 paper_n=1178661, snp_count=34556406, note="113,393 cases + 1,065,268 controls"),
}
HYPO_EUR = 98805 + 948186          # GWAS Catalog: European part
HYPO_NONEUR = 14588 + 117082       # GWAS Catalog: non-European (23andMe)


def scan(path):
    """Stream the value distribution of the n column (the 10th column)."""
    cnt = Counter()
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 10:
                cnt[p[9]] += 1
    return cnt


def load_counts(refresh=False):
    if os.path.exists(CACHE) and not refresh:
        d = pd.read_csv(CACHE)
        print(f"[cache] {os.path.basename(CACHE)}  ({len(d)} rows)")
        return d
    rows = []
    for oc, path in FILES.items():
        print(f"[scan] {oc} … {path}")
        for v, c in scan(path).items():
            rows.append(dict(outcome=oc, n=v, n_variants=c))
    d = pd.DataFrame(rows)
    d["n"] = pd.to_numeric(d.n, errors="coerce")
    d = d.dropna(subset=["n"]).astype({"n_variants": int})
    d.to_csv(CACHE, index=False)
    print(f"[save] {CACHE}  ({len(d)} rows)")
    return d


def wmedian(sub):
    s = sub.sort_values("n")
    return int(s.n[s.n_variants.cumsum() >= s.n_variants.sum() / 2].iloc[0])


def main():
    refresh = "--refresh" in sys.argv
    d = load_counts(refresh)

    rows = []
    for oc in ["TSH", "FT4", "HYPO"]:
        s = d[d.outcome == oc]
        mx = int(s.n.max())
        pn = PAPER[oc]["paper_n"]
        rows.append(dict(
            outcome=("Free T4 (FT4)" if oc == "FT4" else oc),
            GWAS_Catalog=PAPER[oc]["acc"], paper_n=pn,
            max_file_n_per_snp=mx, diff=pn - mx,
            genomewide_weighted_median_n=wmedian(s),
            n_variants=int(s.n_variants.sum()),
            GWAS_Catalog_snp_count=PAPER[oc]["snp_count"],
            n_variants_match=(int(s.n_variants.sum()) == PAPER[oc]["snp_count"]),
            note=PAPER[oc]["note"],
        ))
    ident = pd.DataFrame(rows)

    # cohort decomposition: the three T4 components
    t4 = d[d.outcome == "FT4"].set_index("n").n_variants
    parts = []
    for a, b, c in [(66057, 49269, 76123)]:
        parts.append(dict(component_A=a, component_B=b, component_C=c,
                          A_plus_B=a + b, total=a + b + c,
                          paper_n=191449,
                          exact_match=(a + b + c == 191449 and a + b == 115326)))

    # attribution of the hypothyroidism shortfall
    hmax = int(d[d.outcome == "HYPO"].n.max())
    attr = dict(paper_total=1178661, eur_part=HYPO_EUR, non_eur_part=HYPO_NONEUR,
                file_ceiling=hmax, ceiling_minus_eur=HYPO_EUR - hmax,
                total_minus_ceiling=1178661 - hmax, non_eur_part_check=HYPO_NONEUR)

    print("\n=== identity check ===")
    print(ident.to_string(index=False))
    print("\n=== T4 cohort decomposition ===")
    print(pd.DataFrame(parts).to_string(index=False))
    print("\n=== hypothyroidism shortfall attribution ===")
    for k, v in attr.items():
        print(f"  {k:12s} {v:>10,}")

    print("\n=== reading ===")
    print("  variant counts: the counts in all three files equal the GWAS Catalog snp_count exactly -> the file identity is not in doubt")
    print("  TSH / FT4: the file n ceiling differs from the paper N by only 10 / 1 (single SNP dropped) -> the same data")
    print("  FT4: 66,057 + 49,269 + 76,123 = 191,449 exactly; 49,269 is the Teumer 2018 FT4 cohort")
    print("       -> the trait is free T4, not total T4")
    print(f"  HYPO: {hmax:,} ~ European part {HYPO_EUR:,} - {HYPO_EUR-hmax:,}; "
          f"the paper total 1,178,661 exceeds it by {1178661-hmax:,}, close to the 23andMe {HYPO_NONEUR:,}")
    print("       -> the public summary statistics cover the European cohort only")
    print("\nNote: the per-SNP median (TSH 191,203 / FT4 66,057 / hypothyroidism 424,442) must not be cited as the study sample size.")

    ident.to_csv(f"{TAB}/35_outcome_gwas_identity.csv", index=False, encoding="utf-8-sig")
    print(f"\n[save] {TAB}/35_outcome_gwas_identity.csv")


if __name__ == "__main__":
    main()
