#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""A 路径：用 GTEx v8 全组织中位 TPM 矩阵，评估 SCRN2/DDX5/PSRC1 的甲状腺表达水平与组织特异性。
输出：results/genes_expression_profile.csv  +  文本汇总
"""
import gzip
import pandas as pd
import numpy as np

GTEX = "raw/GTEx_median_tpm.gct.gz"
OUT = "results/genes_expression_profile.csv"

# 目标基因：核心候选 + 旁系同源 + 阳性对照 + 多效性对照
GENES = ["SCRN2", "DDX5", "PSRC1",  # 核心候选
         "SCRN1", "SCRN3",          # secernin 旁系
         "TG", "TPO", "TSHR",       # 甲状腺经典阳性对照
         "SORT1", "CELSR2"]         # 1p13 多效性对照

# 读 GCT（前2行元数据，第3行列名，之后数据）
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

# 列：Name(0) Description(1) 组织从 index 2 开始
tissues = header[2:]
df = pd.DataFrame(rows, columns=["ensg", "symbol"] + tissues)
for t in tissues:
    df[t] = pd.to_numeric(df[t], errors="coerce")

print(f"匹配到基因 {len(df)} / {len(GENES)}")
print("缺失基因:", [g for g in GENES if g not in set(df.symbol)])

# 甲状腺 vs 全血 列
thy = "Thyroid"
blood = "Whole Blood"

# 组织特异性 tau：tau = sum(1 - x_i/max) / (n-1)
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
    # 甲状腺表达在 54 组织中的排名（降序，1=最高）
    rank = int((vals > thy_val).sum()) + 1
    # 最高表达组织
    top_tissue = tissues[int(np.argmax(vals))]
    recs.append(dict(
        gene=r["symbol"],
        ensg=r["ensg"],
        thyroid_tpm=round(thy_val, 3),
        blood_tpm=round(blood_val, 3),
        thyroid_blood_ratio=round(thy_val / blood_val, 3) if blood_val > 0 else np.nan,
        thyroid_rank=rank,          # 1 = 甲状腺是最高表达组织
        n_tissues=len(tissues),
        top_tissue=top_tissue,
        top_tpm=round(vals.max(), 3),
        tau=round(tau_of(vals), 3),
    ))

out = pd.DataFrame(recs).sort_values("thyroid_tpm", ascending=False)
out.to_csv(OUT, index=False, encoding="utf-8-sig")
print(f"\n已写入 {OUT}\n")

# 打印
print("=" * 78)
print("基因甲状腺表达档案（GTEx v8 中位 TPM，54 组织）")
print("=" * 78)
for _, r in out.iterrows():
    pct = r["thyroid_rank"] / r["n_tissues"] * 100
    tag = ""
    if r["gene"] in ("SCRN2", "DDX5", "PSRC1"):
        tag = "  <<< 核心候选"
    print(f"{r['gene']:<7} 甲状腺TPM={r['thyroid_tpm']:>8.3f}  全血TPM={r['blood_tpm']:>8.3f}  "
          f"甲状腺/全血={r['thyroid_blood_ratio']:>6.2f}  "
          f"组织排名={r['thyroid_rank']}/{r['n_tissues']}({pct:.0f}%位)  最高组织={r['top_tissue']}  tau={r['tau']}{tag}")
