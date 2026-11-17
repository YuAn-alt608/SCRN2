#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p07_locus.py
============
必做 4（位点身份）+ 必做 6（样本重叠）+ 必做 5 补全（GRCh38 坐标）。

位点身份三件事：
  (a) 本地位点内四个结局的最小 P，以及是否过 5e-8；
  (b) GWAS Catalog 已报告关联（v2 REST，本地缓存于 tables/12）；
  (c) 明确写清本地缺 Rand 2025 / Williams 2023 的哨兵表（无法核对其独立位点计数）。

样本重叠：给出三队列的估计与异质性 Q；分队列与去 UKB 需要 Rand 补充表，本地无。
"""
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1"
OUT = f"{BASE}/DATA21_complete_2026-11-14"
CENTER = 45915395

# ---------------------------------------------------------------- (a) 区域最小 P
rows = []
for oc in ["TSH", "T4", "HYPO", "HASH"]:
    d = pd.read_csv(f"{BASE}/DATA15_thyroid_eQTL/results/gwas_win/{oc}.tsv.gz", sep="\t")
    d = d[(d.chromosome.astype(str) == "17") &
          (d.base_pair_location.between(CENTER - 500000, CENTER + 500000))]
    d = d.assign(p=pd.to_numeric(d.p_value, errors="coerce"))
    i = d.p.idxmin()
    rs = d.loc[i, "rs_id"] if "rs_id" in d and pd.notna(d.loc[i, "rs_id"]) else "-"
    rows.append(dict(outcome=oc, n_snp=len(d), min_p=float(d.loc[i, "p"]),
                     min_p_snp=rs, pos_grch37=int(d.loc[i, "base_pair_location"]),
                     passes_5e8=bool(d.loc[i, "p"] < 5e-8),
                     n_snps_below_5e8=int((d.p < 5e-8).sum())))
lt = pd.DataFrame(rows)
lt.to_csv(f"{OUT}/tables/13_位点身份_区域最小P.csv", index=False, encoding="utf-8-sig")
print(lt.to_string(index=False))

# ---------------------------------------------------------------- (b) 已报告关联
try:
    gc = pd.read_csv(f"{OUT}/tables/12_GWAS目录_位点身份.csv")
    print("\nGWAS Catalog 命中（本窗口 4 个变异）：")
    print(gc[["rsid", "p", "efo", "study", "genes", "author"]].to_string(index=False))
except Exception as e:
    print("GWAS Catalog 表缺失", e)

# ---------------------------------------------------------------- (c) 三队列 + 异质性
tc = pd.DataFrame([
    dict(cohort="Rand2025_EUR (discovery)", n=482873, est=0.03359, se=0.00885, p=1.48e-4),
    dict(cohort="GCST90319321_EUR (= Figueredo)", n=164818, est=0.04038, se=0.01444, p=5.16e-3),
    dict(cohort="GCST90448634_SAS (Jacobs)", n=38000, est=0.02882, se=0.02413, p=0.232),
])
w = 1 / tc.se ** 2
b_pool = float((w * tc.est).sum() / w.sum())
se_pool = float(np.sqrt(1 / w.sum()))
Q = float((w * (tc.est - b_pool) ** 2).sum())
pq = float(stats.chi2.sf(Q, len(tc) - 1))
tc["weight"] = w / w.sum()
meta = pd.DataFrame([dict(cohort="fixed-effect pooled", n=int(tc.n.sum()), est=b_pool,
                          se=se_pool, p=2 * stats.norm.sf(abs(b_pool / se_pool)),
                          Q=Q, df=len(tc) - 1, p_Q=pq, I2=max(0, (Q - (len(tc) - 1)) / Q) * 100)])
tc.to_csv(f"{OUT}/tables/14_三队列_样本重叠.csv", index=False, encoding="utf-8-sig")
meta.to_csv(f"{OUT}/tables/14b_三队列_合并.csv", index=False, encoding="utf-8-sig")
print("\n", tc.to_string(index=False))
print(meta.to_string(index=False))

# ---------------------------------------------------------------- 补 GRCh38 坐标
vt = pd.read_csv(f"{OUT}/tables/11_变异信息表.csv")
th = pd.read_csv(f"{OUT}/raw/thyroid_region_g38.tsv", sep="\t", usecols=["rsid", "position"])
th = th.dropna().drop_duplicates("rsid").rename(columns={"position": "pos_grch38_new"})
vt = vt.drop(columns=["pos_grch38"]).merge(th, left_on="rsid", right_on="rsid", how="left")
vt = vt.rename(columns={"pos_grch38_new": "pos_grch38"})
vt.to_csv(f"{OUT}/tables/11_变异信息表.csv", index=False, encoding="utf-8-sig")
print("\n", vt.to_string(index=False))
