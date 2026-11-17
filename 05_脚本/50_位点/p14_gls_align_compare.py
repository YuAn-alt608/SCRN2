#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
[已被 p15_ld_true_dosage.py 取代，仅留作方法学对照]
p14_gls_align_compare.py
========================
核对 s21_reference_frame.py 的 GLS 三工具合并估计（0.0250 / 0.0080）到底是怎么算出来的，
并回答一个方法学问题：

  s21 的做法  S = outer(se,se) * R * outer(sign(gamma), sign(gamma))
  这是把 LD 相关矩阵按 sign(γ) 翻符号。如果 sign(γ) 恰好等于"效应等位基因是否 ALT"，
  那它就是标准的等位基因对齐；但 γ 的符号来自生物学效应方向，不是等位基因编码方向，
  两者只有在"效应等位基因 = ALT 且效应为正"时才一致。所以这是一个隐式假设，
  需要用真正的 REF/ALT 对齐 LD 来验证。

本脚本输出 4 种算法：
  A  R=tool_R.npy（ALT 剂量，含负相关） + sign(γ) 翻符号   ← s21 原做法
  B  R=tool_ld_matrix.csv（LDlink 正相关） + 不翻符号       ← 另一种常见写法
  C  R=1000G 等位基因对齐（ld_TSH_aligned）  三工具 rs4794224 缺失时的替代组合
  D  等位基因对齐 LD，逐工具 Wald + IVW（无 Q 加权，仅作参考）

输出 tables/19_GLS方向处理对比.csv 与 tables/19b_三工具_逐工具Wald.csv
"""
import math
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
T = f"{BASE}/tables"
D18 = "D:/NewBeginning/MR1/DATA18_manuscript/analysis"
REG = "D:/NewBeginning/MR1/DATA15_thyroid_eQTL/results/regional"

TOOLS3 = ["rs6503875", "rs62076102", "rs4794224"]        # s21 用的三个
TOOLS3_ALT = ["rs6503875", "rs62076102", "rs11079797"]   # 1000G 面板内可用的三个


def wald(bg, seg, bo, seo):
    if bg == 0 or bo == 0:
        return np.nan, np.nan
    b = bo / bg
    se = abs(b) * math.sqrt((seo / bo) ** 2 + (seg / bg) ** 2)
    return b, se


def gls(b, se, R):
    S = np.outer(se, se) * R
    Si = np.linalg.inv(S)
    one = np.ones(len(b))
    est = float(one @ Si @ b / (one @ Si @ one))
    se_est = float(math.sqrt(1 / (one @ Si @ one)))
    Q = float((b - est) @ Si @ (b - est))
    df = len(b) - 1
    return est, se_est, Q, df, (stats.chi2.sf(Q, df) if df > 0 else np.nan)


# ---------- 1. 逐工具 Wald ----------
rows = []
store = {}
for oc in ["TSH", "T4"]:
    d = pd.read_csv(f"{REG}/SCRN2__{oc}.tsv", sep="\t")
    d = d.set_index("snp")
    for s in TOOLS3 + ["rs11079797"]:
        if s not in d.index:
            continue
        r = d.loc[s]
        b, se = wald(r.beta_eqtl, r.se_eqtl, r.beta, r.se)
        rows.append(dict(pheno=oc, snp=s, ea=r.ea, oa=r.oa,
                         beta_eqtl=r.beta_eqtl, se_eqtl=r.se_eqtl,
                         beta_out=r.beta, se_out=r.se,
                         F=(r.beta_eqtl / r.se_eqtl) ** 2,
                         b_wald=b, se_wald=se,
                         p_wald=2 * stats.norm.sf(abs(b / se))))
        store[(oc, s)] = (b, se, r.beta_eqtl)
wd = pd.DataFrame(rows)
wd.to_csv(f"{T}/19b_三工具_逐工具Wald.csv", index=False, encoding="utf-8-sig")
print("=== 逐工具 Wald（甲状腺 eQTL）===")
print(wd[wd.pheno == "TSH"].round(5).to_string(index=False))

# ---------- 2. 四种 LD 处理 ----------
R_npy = np.load(f"{D18}/tool_R.npy")
nm = pd.read_csv(f"{D18}/tool_ld_matrix.csv", index_col=0)
R_csv = nm.values
names_csv = nm.index.tolist()

ld_al = pd.read_csv(f"{BASE}/inputs_ld/_deprecated_gwas_ea/ld_TSH_aligned.tsv.gz", sep="\t", index_col=0)
ld_al.index = ld_al.index.astype(int)
ld_al.columns = ld_al.columns.astype(int)
pos_of = {"rs6503875": 45915395, "rs62076102": 45887132, "rs11079797": 45912302}


def get_R(tools, mode):
    if mode == "npy":
        sub = pd.DataFrame(R_npy, index=names_csv, columns=names_csv).loc[tools, tools]
        return sub.values
    if mode == "csv":
        sub = nm.loc[tools, tools]
        return sub.values
    if mode == "aligned":
        p = [pos_of[s] for s in tools]
        return ld_al.loc[p, p].values
    raise ValueError(mode)


out = []
for pheno in ["TSH", "T4"]:
    for tools, tag in [(TOOLS3, "三工具(s21: rs6503875+rs62076102+rs4794224)"),
                       (TOOLS3_ALT, "三工具(1000G面板: rs6503875+rs62076102+rs11079797)"),
                       (["rs6503875", "rs62076102"], "两工具(rs6503875+rs62076102)"),
                       (["rs6503875"], "单工具(rs6503875)")]:
        bs = np.array([store[(pheno, s)][0] for s in tools])
        ses = np.array([store[(pheno, s)][1] for s in tools])
        gam = np.array([store[(pheno, s)][2] for s in tools])
        modes = ["npy", "csv", "aligned"] if not any(s == "rs4794224" for s in tools) else ["npy", "csv"]
        for mode in modes:
            R = get_R(tools, mode)
            for flip in ([False, True] if mode != "aligned" else [False]):
                Rf = R * (np.outer(np.sign(gam), np.sign(gam)) if flip else 1.0)
                try:
                    est, se_est, Q, df, pq = gls(bs, ses, Rf)
                except np.linalg.LinAlgError:
                    continue
                out.append(dict(pheno=pheno, tools=tag, k=len(tools),
                                ld_source=mode, sign_gamma_flip=flip,
                                b=est, se=se_est,
                                p=2 * stats.norm.sf(abs(est / se_est)),
                                Q=Q, df=df, p_Q=pq,
                                gamma_signs=";".join("+" if g > 0 else "-" for g in gam)))
res = pd.DataFrame(out)
res.to_csv(f"{T}/19_GLS方向处理对比.csv", index=False, encoding="utf-8-sig")
print("\n=== GLS 四种 LD 处理对比 ===")
print(res[res.pheno == "TSH"].round(6).to_string(index=False))

print("\n=== 与 s21 reference_frame.csv 对照 ===")
ref = pd.read_csv(f"{D18}/reference_frame.csv")
print(ref[ref.comparison.str.contains("GLS")][["comparison", "b_tsh", "se_tsh", "b_t4", "se_t4"]].to_string(index=False))
hit = res[(res.pheno == "TSH") & (np.abs(res.b - 0.025) < 0.0015)]
print("\n接近 0.025 的算法：")
print(hit.round(6).to_string(index=False) if len(hit) else "  无")

# ---------- 3. 敏感性：用对齐 LD 的 0.0336 单工具 vs 0.025 GLS 对结论的影响 ----------
print("\n=== 对 §12 结论的影响 ===")
for lab, b, s in [("单工具 0.0336(0.0088)", 0.0335877, 0.0088515)] + \
                 [(f"{r.tools}/{r.ld_source}/flip={r.sign_gamma_flip}", r.b, r.se)
                  for r in res[(res.pheno == "TSH") & (res.k == 3)].itertuples()]:
    req = -0.206 * b
    se_req = math.sqrt((-0.206 * s) ** 2 + (b * 0.013) ** 2)
    t4, se4 = 0.004976, 0.011956
    z = (t4 - req) / math.sqrt(se4 ** 2 + se_req ** 2)
    print(f"  {lab:52s} b={b:+.4f} req={req:+.4f} z={z:+.2f} p={2*stats.norm.sf(abs(z)):.3f}")
