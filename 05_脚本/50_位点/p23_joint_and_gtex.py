#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p23_joint_and_gtex.py —— 两件事：
  (1) 肩部位点（chr17:45,873,049）在 49 个 GTEx v8 组织里是不是某个基因的 eQTL
  (2) 用 GLS 联合模型（b_j = (R'V⁻¹R)⁻¹R'V⁻¹b，V = diag(se)·R·diag(se)）判定
      45.800 / 45.807 / 45.873 / 45.887 / 45.915 这五个位点里，哪些在条件化后仍然独立

为什么用 GLS 而不是 z 尺度的 R⁻¹z：
  后者隐含"各 SNP 的 SE 相同"。本区间 SE 从 0.0023 到 0.0093 差了 4 倍，
  且 rs62076102 与 rs6503875 的 r=0.767 近共线，R⁻¹z 会把 r≈0.03 的两位点也"吸收"掉
  （实测 rs6503875 z 从 −4.22 塌到 −1.15），这与 LD 结构矛盾。GLS 才是正确形式。
"""
import glob
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
T = f"{BASE}/tables"
SHOULDER = 45_873_049
SHOULDER_B38 = 47_795_683          # GRCh37 + 1,922,634

# ---------------- (1) GTEx 跨组织 ----------------
GTEX = "D:/NewBeginning/MR1/DATA10/tmp/gtex/GTEx_Analysis_v8_eQTL"
rows = []
if os.path.exists(f"{T}/30_肩部GTEx跨组织.csv"):
    print("（复用已生成的 30_肩部GTEx跨组织.csv）")
for f in [] if os.path.exists(f"{T}/30_肩部GTEx跨组织.csv") else sorted(glob.glob(f"{GTEX}/*.v8.signif_variant_gene_pairs.txt.gz")):
    tis = os.path.basename(f).split(".v8.signif")[0]
    d = pd.read_csv(f, sep="\t", usecols=["variant_id", "gene_id", "pval_nominal", "slope"])
    v = d.variant_id.str.split("_", expand=True)
    d = d[(v[0] == "chr17") & (v[1].astype(int) == SHOULDER_B38)]
    if not len(d):
        continue
    for _, r in d.iterrows():
        rows.append(dict(组织=tis, 基因=r.gene_id, P=r.pval_nominal, slope=r.slope))
if not rows:
    g = pd.read_csv(f"{T}/30_肩部GTEx跨组织.csv")
    print(f"=== 肩部位点 GTEx eQTL：{len(g)} 条，{g.组织.nunique()} 组织，{g.基因.nunique()} 基因 ===")
    print(g.groupby("基因名").agg(组织数=("组织", "nunique"), 最强P=("P", "min")).to_string())
else:
    g = pd.DataFrame(rows)
if rows:
    # 补基因名（窗口内的 23 个 + 窗口外但该位点仍命中的 2 个）
    ens = pd.read_csv(f"{BASE}/raw/ens_genes_17_peak_b37.csv")[["gene_id", "gene_name"]]
    name_map = dict(zip(ens.gene_id, ens.gene_name))
    name_map.update({"ENSG00000178852": "EFCAB13", "ENSG00000228782": "MRPL45P2"})
    g["基因名"] = g.基因.str.split(".").str[0].map(name_map)
    # 面板的甲状腺 eQTL 就是 GTEx v8 甲状腺（QTD000341, n=574）——同源，不能算独立验证
    g["与面板同源"] = g.组织 == "Thyroid"
    g = g.sort_values("P")
    g.to_csv(f"{T}/30_肩部GTEx跨组织.csv", index=False, encoding="utf-8-sig")
print(f"=== 肩部位点在 GTEx 中的显著 eQTL：{len(g)} 条，{g.组织.nunique()} 个组织，{g.基因.nunique()} 个基因 ===")
if len(g):
    print(g.head(15).round(6).to_string(index=False))
print("\n按基因汇总（组织数 / 最强 P）:")
print(g.groupby("基因名").agg(组织数=("组织", "nunique"), 最强P=("P", "min"),
                            最强组织=("P", lambda s: g.loc[s.idxmin(), "组织"])).to_string())

# ---------------- (2) 联合条件模型 ----------------
ld = pd.read_csv(f"{BASE}/inputs_ld/ld_TSH.tsv.gz", sep="\t", index_col=0)
ld.index = ld.index.astype(int)
ld.columns = ld.columns.astype(int)

BASE_SNP = {"rs7219420": 45_800_299, "rs1989291": 45_807_547, "肩部": SHOULDER,
            "rs62076102": 45_887_132, "rs6503875": 45_915_395}
MODELS = {
    "五SNP": dict(BASE_SNP),
    "四SNP(去rs62076102)": {k: v for k, v in BASE_SNP.items() if k != "rs62076102"},
    "对:rs7219420+rs1989291": {"rs7219420": 45_800_299, "rs1989291": 45_807_547},
    "对:rs1989291+肩部": {"rs1989291": 45_807_547, "肩部": SHOULDER},
    "对:肩部+rs62076102": {"肩部": SHOULDER, "rs62076102": 45_887_132},
    "对:rs62076102+rs6503875": {"rs62076102": 45_887_132, "rs6503875": 45_915_395},
    "对:rs7219420+rs6503875": {"rs7219420": 45_800_299, "rs6503875": 45_915_395},
    "对:肩部+rs6503875": {"肩部": SHOULDER, "rs6503875": 45_915_395},
}


def gls_joint(d, snps):
    pos = list(snps.values())
    R = ld.loc[pos, pos].values
    b = np.array([d.loc[p, "beta_out"] for p in pos])
    se = np.array([d.loc[p, "se_out"] for p in pos])
    V = np.outer(se, se) * R
    Vi = np.linalg.inv(V)
    M = np.linalg.inv(R.T @ Vi @ R)
    bj = M @ (R.T @ Vi @ b)
    return bj, np.sqrt(np.diag(M)), R


out = []
for oc in ["TSH", "T4", "HYPO", "HASH"]:
    f = f"{BASE}/inputs/Thyroid__SCRN2__{oc}.tsv"
    d = pd.read_csv(f, sep="\t", usecols=["pos37", "beta_out", "se_out"])
    d = d.drop_duplicates("pos37").set_index("pos37")
    for mname, snps in MODELS.items():
        bj, sej, R = gls_joint(d, snps)
        for i, (n, p) in enumerate(snps.items()):
            z_m = d.loc[p, "beta_out"] / d.loc[p, "se_out"]
            out.append(dict(结局=oc, 模型=mname, SNP=n, pos37=p,
                            边际b=d.loc[p, "beta_out"], 边际se=d.loc[p, "se_out"],
                            边际z=z_m, 边际P=2 * stats.norm.sf(abs(z_m)),
                            联合b=bj[i], 联合se=sej[i], 联合z=bj[i] / sej[i],
                            联合P=2 * stats.norm.sf(abs(bj[i] / sej[i])),
                            与模型首SNP_r=R[0, i] if i else 1.0))
j = pd.DataFrame(out)
j.to_csv(f"{T}/31_联合条件模型.csv", index=False, encoding="utf-8-sig")

print("\n注：GTEx 的 Thyroid 与面板甲状腺数据同源（QTD000341），不构成独立验证；其余 22 个组织独立。")
print("\n=== 联合模型：TSH ===")
q = j[(j.结局 == "TSH")]
print(q[["模型", "SNP", "边际z", "联合z", "联合P"]].round(4).to_string(index=False))
print("\n=== 联合模型：甲减 ===")
q = j[(j.结局 == "HYPO")]
print(q[["模型", "SNP", "边际z", "联合z", "联合P"]].round(4).to_string(index=False))
