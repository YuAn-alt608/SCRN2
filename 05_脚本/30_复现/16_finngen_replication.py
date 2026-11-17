# -*- coding: utf-8 -*-
"""
Step 16: 用 FinnGen R12 对 733 条 FDR 显著命中做独立样本复现

原分析结局来自 EBI GWAS Catalog / 大型甲功 meta（UKB 系），FinnGen R12 为
独立队列，可用于检验"同一基因的工具变量在另一个队列里是否给出同向效应"。

关键对照（不是简单重复）= 跨表型方向一致性：
    TSH 升高 -> 甲减风险升高（同向）
    总 T4 升高 -> 甲减风险降低（反向）
若这两个预期都成立，说明复现结果不只是"同队列相关"，而是符合甲状腺轴生理。

方法：
  eQTLGen 只给 Z 值，按 beta = z / sqrt(2*maf*(1-maf)*(N+z^2)) 还原，
  MAF 用 FinnGen 的 af_alt 近似（同为欧洲人群）。等位基因按 FinnGen 对齐。
  多工具变量用 IVW 合并，单工具变量即 Wald ratio。

输出:
  results/finngen_replication.csv   逐条复现结果
  results/finngen_summary.txt       分结局汇总
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
TMP = os.path.join(BASE, "tmp")
DATA = os.path.join(BASE, "data")

N_EQTL = 31684  # eQTLGen 样本量

# 结局代码 -> (FinnGen 文件, 该结局与原二分类结局的预期方向)
#   +1 : 该结局升高应与甲减风险同向
#   -1 : 该结局升高应与甲减风险反向
#    0 : 本身就是甲减类结局，直接同向
MAP = {
    "HYPO": ("finngen_E4_HYTHY_AI_STRICT.tsv", 0, "甲减（自身免疫型，独立队列，直接复现）"),
    "TSH":  ("finngen_E4_HYTHY_AI_STRICT.tsv", +1, "TSH 命中 -> 甲减风险应同向"),
    "T4":   ("finngen_E4_HYTHY_AI_STRICT.tsv", -1, "总 T4 命中 -> 甲减风险应反向"),
    "HASH": ("finngen_E4_THYROID.tsv", 0, "桥本 -> 甲状腺病（同队列不同终点，非独立）"),
}


def load_fg(fn):
    d = pd.read_csv(os.path.join(TMP, fn), sep="\t", dtype={"#chrom": str})
    d.columns = [c.lstrip("#").strip() for c in d.columns]
    d = d.rename(columns={"rsids": "snp", "alt": "ea", "ref": "nea",
                          "sebeta": "se", "af_alt": "af"})
    d = d[["snp", "ea", "nea", "beta", "se", "pval", "af"]].dropna(subset=["snp"])
    # 一个位点可能有多行（多 rsid），保留因果位点首个
    d = d.drop_duplicates("snp")
    return d.set_index("snp")


COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
fg_cache = {}


def replicate(genes_df, fg):
    """genes_df: 每个基因一行，含 snp/ea/nea/z。返回 IVW 估计。"""
    num = den = 0.0
    used = 0
    for r in genes_df.itertuples():
        if r.snp not in fg.index:
            continue
        f = fg.loc[r.snp]
        # 对齐到 eQTL 的效应等位基因
        if r.ea == f["ea"] or COMP.get(r.ea) == f["ea"]:
            sgn = 1.0
        elif r.ea == f["nea"] or COMP.get(r.ea) == f["nea"]:
            sgn = -1.0
        else:
            continue
        if not np.isfinite(f["beta"]) or not np.isfinite(f["se"]) or f["se"] <= 0:
            continue
        af = f["af"]
        if not np.isfinite(af):
            continue
        maf = min(af, 1 - af)
        if maf <= 0.001:
            continue
        info = 2 * maf * (1 - maf) * (N_EQTL + r.z ** 2)
        if info <= 0:
            continue
        b_exp = r.z / np.sqrt(info)
        if abs(b_exp) < 1e-12:
            continue
        b_out = sgn * f["beta"]
        se_out = f["se"]
        w = b_exp / se_out            # Wald ratio 分子
        den += (b_exp / se_out) ** 2  # IVW 权重
        num += w * (b_out / se_out)
        used += 1
    if used == 0 or den <= 0:
        return np.nan, np.nan, 0
    b = num / den
    se = 1 / np.sqrt(den)
    return b, se, used


def main():
    mr = pd.read_csv(os.path.join(RES, "MR_annotated.csv"), encoding="utf-8-sig")
    mr.columns = [c.strip() for c in mr.columns]
    sig = mr[mr["fdr"].notna() & (mr["fdr"] < 0.05)].copy()
    iv = pd.read_csv(os.path.join(DATA, "ivs.tsv"), sep="\t")
    map_name = {"TSH": "TSH", "Free T4": "T4",
                "Hypothyroidism": "HYPO", "Hashimoto thyroiditis": "HASH"}
    sig["oc"] = sig["outcome"].map(map_name)
    print(f"待复现命中: {len(sig)} 条")

    rows = []
    for oc, (fn, expect, note) in MAP.items():
        if fn not in fg_cache:
            fg_cache[fn] = load_fg(fn)
            print(f"  载入 {fn}: {len(fg_cache[fn])} 个 SNP")
        fg = fg_cache[fn]

        sub = sig[sig["oc"] == oc]
        print(f"\n== {oc} ({len(sub)} 条) vs {fn} | {note}")
        for r in sub.itertuples():
            ins = iv[iv["gene"] == r.gene]
            b, se, n_used = replicate(ins, fg)
            ok = np.isfinite(b) and np.isfinite(se) and se > 0
            p = float(2 * (1 - norm.cdf(abs(b / se)))) if ok else np.nan
            rows.append({
                "gene": r.gene, "outcome": r.outcome, "oc": oc,
                "b_orig": r.b, "p_orig": r.pval,
                "snp": r.snp, "nsnp_orig": r.nsnp, "nsnp_used": n_used,
                "b_fg": b if ok else np.nan, "se_fg": se if ok else np.nan,
                "p_fg": p, "expect_sign": expect,
            })

    out = pd.DataFrame(rows)
    out["sign_ok"] = np.sign(out["b_fg"]) == np.sign(out["b_orig"])
    out["p_fg"] = pd.to_numeric(out["p_fg"], errors="coerce")
    # 关键区分：
    #   sign_ok      = 两个 beta 同号（仅对 HYPO/TSH/HASH 有意义）
    #   physio_ok    = 符合甲状腺轴预期（T4 命中应为反向）
    #   replicated   = 同号 + 名义显著（对 T4 会被系统性低估，故并列给出 physio 版）
    out["physio_ok"] = np.where(out["expect_sign"] < 0, ~out["sign_ok"], out["sign_ok"])
    out["replicated"] = out["sign_ok"] & (out["p_fg"] < 0.05)
    out["replicated_physio"] = out["physio_ok"] & (out["p_fg"] < 0.05)
    out.to_csv(os.path.join(RES, "finngen_replication.csv"), index=False, encoding="utf-8-sig")

    print("\n=== 汇总 ===")
    lines = []
    for oc in MAP:
        s = out[(out["oc"] == oc) & out["b_fg"].notna()]
        if len(s) == 0:
            continue
        line = (f"{oc:5s} n={len(s):3d} | 同向 {int(s['sign_ok'].sum()):3d} "
                f"({100*s['sign_ok'].mean():.1f}%) | 同向且 p<0.05 "
                f"{int(s['replicated'].sum()):3d} ({100*s['replicated'].mean():.1f}%)"
                f" | 符合生理预期且 p<0.05 "
                f"{int(s['replicated_physio'].sum()):3d} ({100*s['replicated_physio'].mean():.1f}%)")
        print(line)
        lines.append(line)

    # 跨表型方向一致性：TSH 命中应同向、T4 命中应反向
    print("\n=== 跨表型方向一致性（对甲减风险）===")
    for oc, expect in [("TSH", +1), ("T4", -1)]:
        s = out[(out["oc"] == oc) & out["b_fg"].notna()]
        if len(s) == 0:
            continue
        if expect > 0:
            agree = s["sign_ok"].mean()
        else:
            agree = (~s["sign_ok"]).mean()
        line = f"{oc:5s} 符合生理预期的比例 {100*agree:.1f}% (n={len(s)})"
        print(line)
        lines.append(line)

    with open(os.path.join(RES, "finngen_summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n已写 results/finngen_replication.csv / finngen_summary.txt")

    # 逐条列出复现成功的例子
    top = out[out["replicated"]].nsmallest(10, "p_fg")
    if len(top):
        print("\n=== 复现最强的前 10 条 ===")
        print(top[["gene", "outcome", "b_orig", "b_fg", "p_fg", "nsnp_used"]]
              .to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
