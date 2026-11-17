#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Path B, extended: from the CellxGene "Human Cell Landscape" (Han 2020, microwell-seq) h5ad,
extract the thyroid subset and ask where SCRN2 is expressed across thyroid cell types,
whether in follicular epithelium (thyrocytes) or in immune and stromal cells.

Output: results/SCRN2_thyroid_celltype.csv plus a text summary
"""
import anndata as ad
import pandas as pd
import numpy as np
import sys

H5AD = "tmp/scrn2_sc/hcl_thyroid.h5ad"
OUT = "results/SCRN2_thyroid_celltype.csv"

print("reading h5ad ...", flush=True)
adata = ad.read_h5ad(H5AD)
print(f"cells={adata.n_obs}, genes={adata.n_vars}", flush=True)

# 1) list the obs columns to identify the tissue and cell-type fields
print("\nobs fields:", list(adata.obs.columns), flush=True)

# 2) locate the tissue field (it may be tissue / tissue_original_ontology_term_id / organ / compartment)
tissue_col = None
for c in adata.obs.columns:
    lc = c.lower()
    if "tissue" in lc and "ontology" not in lc:
        tissue_col = c
        break
if tissue_col is None:
    for c in adata.obs.columns:
        if "tissue" in c.lower():
            tissue_col = c
            break

ct_col = None
for c in adata.obs.columns:
    if c in ("cell_type", "cell_type_original", "celltype", "CellType", "cell_type_name"):
        ct_col = c
        break
if ct_col is None:
    for c in adata.obs.columns:
        if "cell" in c.lower() and "type" in c.lower():
            ct_col = c
            break

print(f"tissue field={tissue_col}, cell-type field={ct_col}", flush=True)

if tissue_col is None:
    print("no tissue field found; printing sample unique values per field for manual inspection")
    for c in adata.obs.columns[:20]:
        uniq = adata.obs[c].astype(str).unique()
        print(f"  {c}: {list(uniq[:6])}")
    sys.exit(0)

# 3) print the unique tissue values that mention thyroid
tvals = adata.obs[tissue_col].astype(str).unique()
thy_vals = [v for v in tvals if "thyroid" in v.lower()]
print(f"\ntissue values containing thyroid: {thy_vals}", flush=True)

if not thy_vals:
    print("the tissue field of this dataset has no thyroid entry; check whether an ontology mapping is needed")
    sys.exit(0)

# 4) extract the thyroid subset
mask = adata.obs[tissue_col].astype(str).str.lower().str.contains("thyroid", na=False)
sub = adata[mask].copy()
print(f"thyroid cells: {sub.n_obs}", flush=True)

# 5) locate SCRN2 in var (it may be named by symbol or by ENSG)
var_names = list(sub.var_names)
scr_idx = None
for i, name in enumerate(var_names):
    if str(name).upper() == "SCRN2":
        scr_idx = i
        break
# also try any gene-symbol column in var
gene_sym_col = None
if scr_idx is None:
    for c in sub.var.columns:
        if c in ("feature_name", "gene_symbol", "symbol", "gene_symbols"):
            gene_sym_col = c
            break
    if gene_sym_col:
        for i, v in enumerate(sub.var[gene_sym_col]):
            if str(v).upper() == "SCRN2":
                scr_idx = i
                break

print(f"SCRN2 index in var: {scr_idx}, gene symbol column used={gene_sym_col}", flush=True)

if scr_idx is None:
    print("!!! SCRN2 not found in this dataset (it may be undetected, or named by ENSG)")
    print("first 10 names in var:", [str(v) for v in var_names[:10]])
    sys.exit(0)

# 6) pull the SCRN2 expression vector (X may be sparse or dense)
X = sub.X
if hasattr(X, "toarray"):
    X = X.toarray()
scr_expr = X[:, scr_idx]

# 7) summarise by cell type
if ct_col is None:
    print("no cell-type field found; writing SCRN2 expression per cell instead")
    sub.obs["SCRN2_expr"] = scr_expr
    print(sub.obs[["SCRN2_expr"]].describe())
    sys.exit(0)

sub.obs["SCRN2_expr"] = scr_expr
g = sub.obs.groupby(ct_col).agg(
    n_cells=("SCRN2_expr", "size"),
    mean_expr=("SCRN2_expr", "mean"),
    pct_expressing=("SCRN2_expr", lambda x: (x > 0).mean() * 100),
    max_expr=("SCRN2_expr", "max"),
).sort_values("mean_expr", ascending=False)

print("\n===== SCRN2 expression by thyroid cell type =====", flush=True)
print(g.to_string(), flush=True)

# 8) also score a few follicular-epithelial marker genes as a reference frame for annotation quality
markers = ["TG", "TPO", "TSHR", "PAX8", "NKX2-1", "PTPRC"]
marker_map = {}
for m in markers:
    idx = None
    for i, name in enumerate(var_names):
        if str(name).upper() == m:
            idx = i
            break
    if idx is not None:
        marker_map[m] = X[:, idx]

if marker_map:
    print("\n===== thyroid marker gene expression (annotation quality check) =====", flush=True)
    for m, vec in marker_map.items():
        # cell type with the highest expression
        sub.obs["_m"] = vec
        top = sub.obs.groupby(ct_col)["_m"].mean().sort_values(ascending=False)
        print(f"{m}: highest-expressing cell type = {top.index[0]} (mean={top.iloc[0]:.3f})", flush=True)
        sub.obs.drop(columns="_m", inplace=True)

# 9) save
g.to_csv(OUT, encoding="utf-8-sig")
print(f"\nwritten to {OUT}", flush=True)
