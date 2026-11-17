#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Path A: from the GTEx v8 median TPM matrix across tissues, assess the thyroid expression level
and tissue specificity of SCRN2/DDX5/PSRC1.
Output: results/genes_expression_profile.csv plus a text summary
"""
import gzip
import pandas as pd
import numpy as np

GTEX = "raw/GTEx_median_tpm.gct.gz"
OUT = "results/genes_expression_profile.csv"

# target genes: core candidates + paralogues + positive controls + pleiotropy controls
GENES = ["SCRN2", "DDX5", "PSRC1",  # core candidates
         "SCRN1", "SCRN3",          # secernin paralogues
         "TG", "TPO", "TSHR",       # classical thyroid positive controls
         "SORT1", "CELSR2"]         # 1p13 pleiotropy controls

# read the GCT (two metadata rows, then the column names, then the data)
rows = []
with gzip.open(GTEX, "rt") as f:
    lines = f
    header = None
    for i, line in enumerate(lines):
        if i == 2:
            header = line.rstrip("\n").split("\t")
            break
    for line in lines:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 3:
            continue
        sym = cols[1]
        if sym in GENES:
            rows.append(cols)

# columns: Name(0), Description(1), tissues from index 2 onwards
tissues = header[2:]
df = pd.DataFrame(rows, columns=["ensg", "symbol"] + tissues)
for t in tissues:
    df[t] = pd.to_numeric(df[t], errors="coerce")

print(f"genes matched {len(df)} / {len(GENES)}")
print("missing genes:", [g for g in GENES if g not in set(df.symbol)])

# thyroid and whole-blood columns
thy = "Thyroid"
blood = "Whole Blood"

# tissue specificity tau: tau = sum(1 - x_i/max) / (n-1)
def tau_of(row_vals):
    v = row_vals.values.astype(float)
    m = v.max()
    if m == 0 or np.isnan(m):
        return np.nan
    n = len(v)
    return np.sum(1 - v / m) / (n - 1)

recs = []
for _, r in df.iterrows():
    vals = r[tissues].astype(float)
    thy_val = r[thy]
    blood_val = r[blood]
    # rank of thyroid expression among the 54 tissues (descending, 1 = highest)
    rank = int((vals > thy_val).sum()) + 1
    # tissue with the highest expression
    top_tissue = tissues[int(np.argmax(vals))]
    recs.append(dict(
        gene=r["symbol"],
        ensg=r["ensg"],
        thyroid_tpm=round(thy_val, 3),
        blood_tpm=round(blood_val, 3),
        thyroid_blood_ratio=round(thy_val / blood_val, 3) if blood_val > 0 else np.nan,
        thyroid_rank=rank,          # 1 = thyroid is the highest-expressing tissue
        n_tissues=len(tissues),
        top_tissue=top_tissue,
        top_tpm=round(vals.max(), 3),
        tau=round(tau_of(vals), 3),
    ))

out = pd.DataFrame(recs).sort_values("thyroid_tpm", ascending=False)
out.to_csv(OUT, index=False, encoding="utf-8-sig")
print(f"\nwritten to {OUT}\n")

# print
print("=" * 78)
print("thyroid expression profile (GTEx v8 median TPM, 54 tissues)")
print("=" * 78)
for _, r in out.iterrows():
    pct = r["thyroid_rank"] / r["n_tissues"] * 100
    tag = ""
    if r["gene"] in ("SCRN2", "DDX5", "PSRC1"):
        tag = "  <<< core candidate"
    print(f"{r['gene']:<7} thyroid TPM={r['thyroid_tpm']:>8.3f}  blood TPM={r['blood_tpm']:>8.3f}  "
          f"thyroid/blood={r['thyroid_blood_ratio']:>6.2f}  "
          f"rank={r['thyroid_rank']}/{r['n_tissues']}({pct:.0f}th pct)  top tissue={r['top_tissue']}  tau={r['tau']}{tag}")
