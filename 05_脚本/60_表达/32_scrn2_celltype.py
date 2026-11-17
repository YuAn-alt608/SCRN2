#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""B 路径深化：从 CellxGene "Human Cell Landscape"（Han 2020, microwell-seq）h5ad
中提取甲状腺组织子集，查 SCRN2 在甲状腺各细胞类型的表达分布，
判断它到底落在滤泡上皮（thyrocyte）还是免疫/间质细胞。

输出：results/SCRN2_thyroid_celltype.csv + 文本汇总
"""
import anndata as ad
import pandas as pd
import numpy as np
import sys

H5AD = "tmp/scrn2_sc/hcl_thyroid.h5ad"
OUT = "results/SCRN2_thyroid_celltype.csv"

print("读取 h5ad ...", flush=True)
adata = ad.read_h5ad(H5AD)
print(f"细胞数={adata.n_obs}, 基因数={adata.n_vars}", flush=True)

# 1) 打印 obs 的字段名，确认组织字段与细胞类型字段
print("\nobs 字段:", list(adata.obs.columns), flush=True)

# 2) 找出组织字段（可能叫 tissue / tissue_original_ontology_term_id / organ / compartment 等）
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

print(f"组织字段={tissue_col}, 细胞类型字段={ct_col}", flush=True)

if tissue_col is None:
    print("未找到组织字段，打印各字段唯一值样本以便人工判断")
    for c in adata.obs.columns[:20]:
        uniq = adata.obs[c].astype(str).unique()
        print(f"  {c}: {list(uniq[:6])}")
    sys.exit(0)

# 3) 打印组织字段里含 thyroid 的唯一值
tvals = adata.obs[tissue_col].astype(str).unique()
thy_vals = [v for v in tvals if "thyroid" in v.lower()]
print(f"\n组织字段含 thyroid 的值: {thy_vals}", flush=True)

if not thy_vals:
    print("该数据集组织字段里没有 thyroid，检查是否需要 ontology 映射")
    sys.exit(0)

# 4) 提取甲状腺子集
mask = adata.obs[tissue_col].astype(str).str.lower().str.contains("thyroid", na=False)
sub = adata[mask].copy()
print(f"甲状腺细胞数: {sub.n_obs}", flush=True)

# 5) 找 SCRN2 基因在 var 中的位置（可能用 symbol 或 ENSG）
var_names = list(sub.var_names)
scr_idx = None
for i, name in enumerate(var_names):
    if str(name).upper() == "SCRN2":
        scr_idx = i
        break
# 也尝试 var 里的 gene symbol 列
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

print(f"SCRN2 在 var 的索引: {scr_idx}, 用到的 gene symbol 列={gene_sym_col}", flush=True)

if scr_idx is None:
    print("!!! 该数据集中未找到 SCRN2 基因（可能未被检测到，或用 ENSG 命名）")
    print("var 前 10 个名字:", [str(v) for v in var_names[:10]])
    sys.exit(0)

# 6) 提取 SCRN2 表达向量（X 可能是稀疏矩阵或 dense）
X = sub.X
if hasattr(X, "toarray"):
    X = X.toarray()
scr_expr = X[:, scr_idx]

# 7) 按细胞类型分组统计
if ct_col is None:
    print("未找到细胞类型字段，直接输出 SCRN2 在每个细胞的表达")
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

print("\n===== SCRN2 在甲状腺各细胞类型的表达 =====", flush=True)
print(g.to_string(), flush=True)

# 8) 也统计几个甲状腺滤泡上皮标记基因，作为参照系（判断细胞注释质量）
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
    print("\n===== 甲状腺标记基因表达（验证注释质量）=====", flush=True)
    for m, vec in marker_map.items():
        # 找表达最高的细胞类型
        sub.obs["_m"] = vec
        top = sub.obs.groupby(ct_col)["_m"].mean().sort_values(ascending=False)
        print(f"{m}: 最高表达细胞类型 = {top.index[0]} (mean={top.iloc[0]:.3f})", flush=True)
        sub.obs.drop(columns="_m", inplace=True)

# 9) 保存
g.to_csv(OUT, encoding="utf-8-sig")
print(f"\n已写入 {OUT}", flush=True)
