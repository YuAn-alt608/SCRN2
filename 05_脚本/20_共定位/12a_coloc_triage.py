#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
12a_coloc_triage.py
-------------------
在跑 coloc.susie 之前，先做一次"是否需要 susie"的分诊（不依赖 LD）。

逻辑：coloc.abf 假设每个性状只有 1 个因果变异。若某个基因-结局对的
eQTL 侧或结局侧存在 >=2 个距离上独立的信号峰（>250 kb），abf 就可能把
"两个性状各自有信号但峰不重合"误判成 H3，此时必须用 susie 复核。
只有 1 个峰的区间，susie 与 abf 结论一致，无需重复计算。

输出 results/susie_triage.csv
"""
import sys, os
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
DATA = os.path.join(BASE, "data")

CLUMP_KB = 250
PV_SIG = 1e-5


def clump_peaks(df, pcol="p"):
    """距离聚类：按 P 升序贪心取峰，剔除 250 kb 内其他 SNP。返回 (n_peaks, lead_positions, lead_ps)"""
    if df.empty:
        return 0, [], []
    d = df.dropna(subset=[pcol, "pos"]).sort_values(pcol)
    d = d[d[pcol] < PV_SIG]
    if d.empty:
        return 0, [], []
    leads, lpos, lp = [], [], []
    for pos, p in zip(d["pos"].values, d[pcol].values):
        if all(abs(pos - q) > CLUMP_KB * 1000 for q in lpos):
            lpos.append(pos)
            lp.append(p)
    return len(lpos), lpos, lp


def main():
    coloc = pd.read_csv(os.path.join(RES, "coloc_results.csv"), encoding="utf-8-sig")
    genes = pd.read_csv(os.path.join(DATA, "coloc_genes_all.tsv"), sep="\t")
    print("coloc pairs:", len(coloc), "| genes:", len(genes), "| gene cols:", genes.columns.tolist())

    gcol = "gene" if "gene" in genes.columns else genes.columns[0]
    ccol = "chr" if "chr" in genes.columns else ("chrom" if "chrom" in genes.columns else None)
    pcol = "pos" if "pos" in genes.columns else ("tss" if "tss" in genes.columns else None)

    # MR 显著性，用于排序优先级
    mr = pd.read_csv(os.path.join(RES, "MR_significant_FDR05_v2.csv"), encoding="utf-8-sig")
    OMAP = {"TSH": "TSH", "Free T4": "T4", "Hypothyroidism": "HYPO", "Hashimoto thyroiditis": "HASH"}
    mr["outcome"] = mr["outcome"].map(OMAP)
    mrkey = mr.groupby(["gene", "outcome"])["pval"].min().rename("p_mr").reset_index()

    # eQTL 区域
    eq = pd.read_csv(os.path.join(DATA, "eqtl_regions_all.txt.gz"), sep="\t",
                     usecols=["Pvalue", "SNP", "SNPChr", "SNPPos", "GeneSymbol"],
                     dtype={"Pvalue": float, "SNPPos": int})
    eq = eq.rename(columns={"Pvalue": "p", "SNPPos": "pos", "GeneSymbol": "gene",
                            "SNPChr": "chr"})
    eq["chr"] = eq["chr"].astype(str).str.replace("chr", "", regex=False)
    # cis 窗口以 eQTL 实际跨度为准（eQTLGen 取基因 ±1 Mb，跨度达 1 Mb 时不能用 TSS 当中心）
    span = eq.groupby("gene").agg(chr=("chr", "first"), lo=("pos", "min"), hi=("pos", "max")).reset_index()
    span = {r.gene: (r.chr, int(r.lo), int(r.hi)) for r in span.itertuples()}
    eq = eq[eq["p"] < PV_SIG]
    eq_g = {k: v for k, v in eq.groupby("gene")}

    triage = []
    for out in ["TSH", "T4", "HYPO", "HASH"]:
        f = os.path.join(DATA, f"region_{out}.tsv.gz")
        if not os.path.exists(f):
            print("missing", f)
            continue
        oo = pd.read_csv(f, sep="\t",
                         usecols=["chromosome", "base_pair_location", "p_value", "rs_id"],
                         dtype={"chromosome": str, "base_pair_location": int, "p_value": float})
        oo = oo.rename(columns={"chromosome": "chr", "base_pair_location": "pos", "p_value": "p", "rs_id": "snp"})
        oo["chr"] = oo["chr"].str.replace("chr", "", regex=False)
        oo = oo[oo["p"] < PV_SIG]
        print(out, "outcome SNPs P<1e-5:", len(oo))
        for _, g in genes.iterrows():
            span_g = span.get(g[gcol])
            if span_g is None:
                continue
            c_g, lo_g, hi_g = span_g
            sub = oo[(oo["chr"] == c_g) & (oo["pos"] >= lo_g) & (oo["pos"] <= hi_g)]
            n_out, pos_out, p_out = clump_peaks(sub)
            e = eq_g.get(g[gcol])
            n_eqtl, pos_eqtl, p_eqtl = clump_peaks(e) if e is not None else (0, [], [])
            triage.append({"gene": g[gcol], "outcome": out,
                           "n_sig_eqtl": n_eqtl, "n_sig_out": n_out,
                           "lead_eqtl_pos": pos_eqtl[0] if pos_eqtl else np.nan,
                           "lead_out_pos": pos_out[0] if pos_out else np.nan,
                           "dist_lead_kb": (abs(pos_eqtl[0] - pos_out[0]) / 1000.0
                                            if (pos_eqtl and pos_out) else np.nan)})
        del oo

    t = pd.DataFrame(triage)
    t = t.merge(coloc[["gene", "outcome", "PP.H3", "PP.H4", "PP.H1"]], on=["gene", "outcome"], how="left")
    t = t.merge(mrkey, on=["gene", "outcome"], how="left")
    t["multi_signal"] = (t["n_sig_eqtl"] >= 2) | (t["n_sig_out"] >= 2)
    # 需要 susie 复核的三种情形
    t["need_susie"] = (
        (t["multi_signal"] & (t["PP.H3"] > 0.5)) |          # 疑似假 H3
        (t["multi_signal"] & (t["PP.H4"] > 0.5)) |          # 需确认共享信号
        (~t["multi_signal"] & (t["PP.H3"] > 0.5) & (t["PP.H4"] > 0.2))  # 单峰但结论不稳
    )
    t = t.sort_values(["need_susie", "p_mr"], ascending=[False, True])
    t.to_csv(os.path.join(RES, "susie_triage.csv"), index=False, encoding="utf-8-sig")

    print("\n=== 分诊结果 ===")
    print(t.groupby("outcome").agg(n=("gene", "size"), multi=("multi_signal", "sum"),
                                   need=("need_susie", "sum")))
    print("need_susie 合计:", int(t["need_susie"].sum()))
    print(t[t["need_susie"]].head(15)[["gene", "outcome", "n_sig_eqtl", "n_sig_out", "PP.H3", "PP.H4", "p_mr"]].to_string())


if __name__ == "__main__":
    main()
