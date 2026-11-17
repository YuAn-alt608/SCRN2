#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p08_summary.py
==============
汇总成一张决策表 + 工作簿 + 热图：
  基因 × 组织 × 结局：lead SNP、eQTL P/F、cis-MR（Wald/IVW）、coloc（H3/H4，p12=1e-4）、
  SMR、HEIDI、有效工具数。并给出 TSH 下的基因排序。
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
T = f"{BASE}/tables"

mr = pd.read_csv(f"{T}/04_cisMR_SMR_HEIDI.csv")
co = pd.read_csv(f"{T}/03_面板_coloc_abf_先验网格.csv")
genes_in = pd.read_csv(f"{T}/01_基因清单.csv").gene_name.tolist()

co4 = co[co.p12 == 1e-4][["tissue", "gene", "outcome", "PP.H3", "PP.H4", "nsnp"]]
co4 = co4.rename(columns={"PP.H3": "coloc_H3", "PP.H4": "coloc_H4",
                          "nsnp": "coloc_nsnp"})
d = mr.merge(co4, on=["tissue", "gene", "outcome"], how="left")
d = d[d.gene.isin(genes_in)]
d["F_lead"] = d.F_lead.round(1)
d["heidi_pass"] = d.heidi_p > 0.01

# 排序分数：TSH，要求 F>10
def score(r):
    s = 0.0
    if r.lead_p_eqtl < 5e-8: s += 1
    if r.p_wald < 0.05: s += 1
    if r.coloc_H4 > 0.8: s += 1
    if r.p_smr < 0.05 and (r.heidi_p > 0.01 if r.heidi_p == r.heidi_p else False): s += 1
    return s

d["证据分_0to4"] = d.apply(score, axis=1)
d.to_csv(f"{T}/15_决策表_全基因.csv", index=False, encoding="utf-8-sig")

tsh = d[d.outcome == "TSH"].sort_values(["证据分_0to4", "p_wald"], ascending=[False, True])
tsh.to_csv(f"{T}/15b_决策表_TSH.csv", index=False, encoding="utf-8-sig")
cols = ["tissue", "gene", "lead_snp", "lead_p_eqtl", "F_lead", "b_wald", "se_wald", "p_wald",
        "n_iv", "k_eff", "coloc_H3", "coloc_H4", "p_smr", "heidi_p", "证据分_0to4"]
pd.set_option("display.width", 260, "display.max_columns", 40)
print(tsh[cols].round(4).to_string(index=False))

# 热图：coloc H4
piv = (d[d.outcome == "TSH"].pivot_table(index="gene", columns="tissue",
                                         values="coloc_H4"))
piv = piv.reindex(piv.max(axis=1).sort_values(ascending=False).index)
fig, ax = plt.subplots(figsize=(4.6, 7))
im = ax.imshow(piv.values, cmap="RdYlBu_r", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=30, ha="right", fontsize=9)
ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index, fontsize=8)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v = piv.values[i, j]
        if v == v:
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if v > 0.6 or v < 0.15 else "black")
ax.set_title("coloc PP.H4 with TSH (p12=1e-4)", fontsize=10)
fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
fig.tight_layout()
fig.savefig(f"{BASE}/fig_面板coloc热图.png", dpi=200)
print("热图已存")

# 工作簿
with pd.ExcelWriter(f"{T}/../MR1_补做结果_2026-11-14.xlsx", engine="openpyxl") as xw:
    tsh[cols].to_excel(xw, sheet_name="TSH决策表", index=False)
    d.to_excel(xw, sheet_name="全结局决策表", index=False)
    co.to_excel(xw, sheet_name="coloc先验网格", index=False)
    mr.to_excel(xw, sheet_name="cisMR_SMR_HEIDI", index=False)
    pd.read_csv(f"{T}/05_rs7219420_边际.csv").to_excel(xw, sheet_name="rs7219420边际", index=False)
    pd.read_csv(f"{T}/06_两SNP联合模型.csv").to_excel(xw, sheet_name="两SNP联合", index=False)
    pd.read_csv(f"{T}/07_两SNP的cis_eQTL.csv").to_excel(xw, sheet_name="两SNP_eQTL", index=False)
    pd.read_csv(f"{T}/11_变异信息表.csv").to_excel(xw, sheet_name="变异信息表", index=False)
    pd.read_csv(f"{T}/13_位点身份_区域最小P.csv").to_excel(xw, sheet_name="位点身份", index=False)
    pd.read_csv(f"{T}/14_三队列_样本重叠.csv").to_excel(xw, sheet_name="样本重叠", index=False)
    for f in ["12_GWAS目录_位点身份.csv"]:
        if os.path.exists(f"{T}/{f}"):
            pd.read_csv(f"{T}/{f}").to_excel(xw, sheet_name="GWAS目录", index=False)
    if os.path.exists(f"{T}/16_变异_cCRE注释.csv"):
        pd.read_csv(f"{T}/16_变异_cCRE注释.csv").to_excel(xw, sheet_name="cCRE注释", index=False)
    if os.path.exists(f"{T}/17_两套gamma敏感性.csv"):
        pd.read_csv(f"{T}/17_两套gamma敏感性.csv").to_excel(xw, sheet_name="两套gamma", index=False)
    if os.path.exists(f"{T}/18_工具数与GLS_IVW.csv"):
        pd.read_csv(f"{T}/18_工具数与GLS_IVW.csv").to_excel(xw, sheet_name="工具数扫描", index=False)
    # 2026-11-14 修正：TSH 侧 SuSiE 结果不可靠（等位基因对齐后仍报 prior 方差过大，
    # 且给出的可信集首位 SNP 与边际 P 矛盾），已删除 03b/03c 的 TSH 版本，
    # 只保留 eQTL 侧（QTD000341 官方 credible_sets）。
    if os.path.exists(f"{T}/03c_可信集_eQTL侧.csv"):
        pd.read_csv(f"{T}/03c_可信集_eQTL侧.csv").to_excel(xw, sheet_name="可信集_eQTL侧", index=False)
    for f, sn in [("19_GLS方向处理对比.csv", "GLS工具集对比"),
                  ("19b_三工具_逐工具Wald.csv", "逐工具Wald"),
                  ("19c_LD目标等位基因.csv", "LD等位基因"),
                  ("20_模拟验证_signγ.csv", "GLS符号模拟验证"),
                  ("21_已发表哨兵_17q21.csv", "已发表哨兵"),
                  ("22_哨兵基因映射唯一性.csv", "哨兵映射唯一性"),
                  ("23_信号解析归属总表.csv", "信号解析归属"),
                  ("24_按信号窗口的coloc.csv", "按窗口coloc"),
                  ("25_HYPO数据来源核验.csv", "甲减来源核验"),
                  ("26_峰周围基因_窗口coloc.csv", "邻近基因coloc"),
                  ("27_肩部基因归属表.csv", "肩部归属"),
                  ("28_峰位点eQTL全扫描.csv", "峰位点eQTL扫描"),
                  ("29_肩部位点逐基因Wald.csv", "肩部逐基因Wald"),
                  ("30_肩部GTEx跨组织.csv", "肩部GTEx跨组织"),
                  ("31_联合条件模型.csv", "联合条件模型")]:
        p = f"{T}/{f}"
        if os.path.exists(p):
            pd.read_csv(p).to_excel(xw, sheet_name=sn, index=False)
    if os.path.exists(f"{T}/08_区域图数据.csv"):
        pd.read_csv(f"{T}/08_区域图数据.csv").to_excel(xw, sheet_name="区域图数据", index=False)
    for f, sn in [("09b_逐组织_最强eQTL_p.csv", "逐组织eQTL"), ("10_逐组织_关键SNP的eQTL.csv", "关键SNP逐组织")]:
        p = f"{T}/{f}"
        if os.path.exists(p):
            pd.read_csv(p).to_excel(xw, sheet_name=sn, index=False)
    for ws in xw.book.worksheets:
        ws.freeze_panes = "A2"
print("工作簿已存")
