#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p22_shoulder_attribution.py —— 45,873,049 肩部（原记作"未归属甲减峰"）的基因归属

背景
----
甲减在 17q21.32 的 P 值剖面（10 kb 分箱）显示：主峰 45.80–45.81（P≈1e-33），
其右侧 45.83–45.87 有一段平台（P 恒在 3e-19），此前被当作"第二个未归属的甲减峰"。
本脚本把该肩部 ±200 kb 内的全部基因（Ensembl GRCh37 端点取的 23 个）纳入，
能测的 17 个跑完：eQTL lead → 与肩部/主峰的 LD → cis-MR → 窗口 coloc → 裁决。

坐标口径
--------
全部 GRCh37（pos37）。偏移 1,922,634 由 Ensembl GRCh37 归档端点逐位点直接核对：
  rs6503875 45,915,395；rs7219420 45,800,299；rs1989291 45,807,547；rs4794063 45,804,494
"""
import os
import numpy as np
import pandas as pd

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
T = f"{BASE}/tables"
PEAK = 45_873_049          # 肩部最优点（甲减 P=3.09e-19）
LEADS = {                  # 三个信号的代表位点
    "肩部": PEAK,
    "主峰_rs1989291": 45_807_547,
    "TSH峰_rs62076102": 45_887_132,
    "SCRN2_rs6503875": 45_915_395,
}

ld = pd.read_csv(f"{BASE}/inputs_ld/ld_TSH.tsv.gz", sep="\t", index_col=0)
ld.index = ld.index.astype(int)
ld.columns = ld.columns.astype(int)

genes = pd.read_csv(f"{BASE}/raw/ens_genes_17_peak_b37.csv")
genes["name"] = genes.gene_name.fillna(genes.gene_id)
genes["tss37"] = np.where(genes.strand == 1, genes.start, genes.end)
genes["dist_peak_kb"] = (genes.tss37 - PEAK) / 1000

# 别名：AC003665.1 是 SP2-AS1 的旧名（同一 ENSG00000234494），面板里两个名字各存了一份
ALIAS = {"SP2-AS1": "AC003665.1"}

mr = pd.read_csv(f"{T}/04_cisMR_SMR_HEIDI.csv")
co = pd.read_csv(f"{T}/26_峰周围基因_窗口coloc.csv")

rows = []
for _, g in genes.iterrows():
    name = g["name"]
    for tissue, lab in [("Thyroid", "甲状腺"), ("WholeBlood", "全血")]:
        f = f"{BASE}/inputs/{tissue}__{name}__TSH.tsv"
        if not os.path.exists(f):
            alt = [k for k, v in ALIAS.items() if v == name]
            if alt:
                f = f"{BASE}/inputs/{tissue}__{alt[0]}__TSH.tsv"
            if not os.path.exists(f):
                continue
        d = pd.read_csv(f, sep="\t", usecols=["pos37", "p_eqtl", "beta_eqtl"])
        row = dict(基因=name, gene_id=g.gene_id, 组织=lab, tss37=int(g.tss37),
                   距肩部_kb=round(g.dist_peak_kb, 1))
        lead = d.loc[d.p_eqtl.idxmin()]
        row.update(eQTL_lead=int(lead.pos37), eQTL_lead_p=lead.p_eqtl,
                   距肩部_eQTLlead_kb=round((int(lead.pos37) - PEAK) / 1000, 1),
                   eQTL_at_peak=float(d.loc[d.pos37 == PEAK, "p_eqtl"].iloc[0])
                   if (d.pos37 == PEAK).any() else np.nan)
        # 与各信号 lead 的 LD
        for k, v in LEADS.items():
            if v in ld.index and int(lead.pos37) in ld.index:
                row[f"r_{k}"] = round(float(ld.loc[int(lead.pos37), v]), 3)
            else:
                row[f"r_{k}"] = np.nan
        # cis-MR
        for oc in ["HYPO", "TSH"]:
            q = mr[(mr.gene == name) & (mr.tissue == tissue) & (mr.outcome == oc)]
            if not len(q) and name in ALIAS:
                q = mr[(mr.gene == ALIAS[name]) & (mr.tissue == tissue) & (mr.outcome == oc)]
            if len(q):
                row[f"MR_{oc}_b"] = float(q.b_wald.iloc[0])
                row[f"MR_{oc}_p"] = float(q.p_wald.iloc[0])
                row[f"MR_{oc}_F"] = float(q.F_lead.iloc[0])
                row[f"MR_{oc}_niv"] = int(q.n_iv.iloc[0])
        # 窗口 coloc（峰±50kb）
        q2 = co[(co.gene == name) & (co.tissue == tissue) &
                (co.outcome == "HYPO") & (co.window.str.contains("50kb"))]
        row["coloc_H3_HYPO"] = float(q2.H3.iloc[0]) if len(q2) else np.nan
        row["coloc_H4_HYPO"] = float(q2.H4.iloc[0]) if len(q2) else np.nan
        rows.append(row)

s = pd.DataFrame(rows)


s["该位点组织内排名"] = (s.groupby("组织")["eQTL_at_peak"]
                            .rank(method="min").astype("Int64"))


def verdict(r):
    if pd.isna(r.get("eQTL_lead_p")):
        return "无 eQTL 数据"
    at_peak = r.get("eQTL_at_peak", np.nan)
    lead_dist = abs(r.get("距肩部_eQTLlead_kb", np.nan))
    h4 = r.get("coloc_H4_HYPO", np.nan)
    h3 = r.get("coloc_H3_HYPO", np.nan)

    # 只有"该位点上该组织的首位 eQTL"才算位点层面重合；否则任何强 eQTL 基因
    # 都会因单倍型连锁在这个位点"显著"，无法区分（实测全血 11 个基因 P<1e-5）
    if r.get("该位点组织内排名", 99) == 1 and at_peak < 1e-3 and lead_dist > 20:
        return ("位点层面重合（该变异本身是该基因 eQTL，P=%.1e），"
                "但主 eQTL 在 %+.0f kb 外，窗口 coloc 判 H3=%.2f（单因果假设被违反）→ 候选"
                % (at_peak, r["距肩部_eQTLlead_kb"], h3 if not pd.isna(h3) else np.nan))
    if not pd.isna(h4) and h4 > 0.75:
        return "共享（coloc H4=%.2f）" % h4
    base = ("不共享（H3=%.2f）" % h3) if (not pd.isna(h3) and h3 > 0.75) else "不可判定"
    if not pd.isna(r.get("MR_HYPO_p", np.nan)) and r.MR_HYPO_p < 0.01:
        base += "；但甲减 cis-MR 显著（P=%.1e）" % r.MR_HYPO_p
    return base


s["裁决"] = [verdict(r) for _, r in s.iterrows()]
s = s.sort_values(["组织", "coloc_H4_HYPO"], ascending=[True, False])
s.to_csv(f"{T}/27_肩部基因归属表.csv", index=False, encoding="utf-8-sig")

print(f"基因-组织行数 {len(s)}｜基因 {s.基因.nunique()}")
print("\n=== 甲减：局部窗口 coloc（峰±50kb）H4 排序 ===")
c = [x for x in ["基因", "组织", "距肩部_kb", "eQTL_lead", "eQTL_lead_p", "r_肩部",
                 "MR_HYPO_b", "MR_HYPO_p", "coloc_H3_HYPO", "coloc_H4_HYPO"] if x in s.columns]
print(s[c].sort_values("coloc_H4_HYPO", ascending=False).head(14).round(3).to_string(index=False))

print("\n=== 有甲减 cis-MR 信号的行（P<0.01）===")
q = s[(s.MR_HYPO_p < 0.01)].sort_values("MR_HYPO_p")
print(q[["基因", "组织", "eQTL_lead", "r_肩部", "r_主峰_rs1989291", "r_SCRN2_rs6503875",
         "MR_HYPO_b", "MR_HYPO_p", "MR_HYPO_F", "MR_HYPO_niv"]].round(4).to_string(index=False))
