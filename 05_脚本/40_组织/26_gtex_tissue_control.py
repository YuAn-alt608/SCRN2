# -*- coding: utf-8 -*-
"""
Step 26: GTEx 多组织对照（本地 v8 全组织，与第 06 节同一数据源、同一口径）

问题：第 06 节只做了甲状腺，34% 检出率无法判断高低。
      如果同样的工具变量在肝、肌肉、肺里也有三成检出，"甲状腺验证"就不成立。

做法：对 A/B 级全部 (基因, 工具变量) 对，在 49 个组织的
      *.v8.signif_variant_gene_pairs.txt.gz 中按 GRCh38 坐标 + ENSG 匹配，
      统计各组织的检出率（p<0.05）与方向一致率，看甲状腺排第几。

      与 06 节的区别：06 节只查甲状腺、用 64 个可映射 SNP 作分母；
      本步覆盖全部 49 个组织、70 个 (基因,SNP) 对作分母，口径更严。

输入：tmp/gtex/GTEx_Analysis_v8_eQTL/*.signif_variant_gene_pairs.txt.gz
输出：
  results/gtex_tissue_control.csv   长表（gene, snp, tissue, p, slope）
  results/gtex_tissue_rates.csv     逐组织检出率与一致率
  results/gtex_tissue_summary.txt
"""
import os, json, glob, gzip
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
DATA = os.path.join(BASE, "data")
TMP = os.path.join(BASE, "tmp")
GDIR = os.path.join(TMP, "gtex", "GTEx_Analysis_v8_eQTL")
PAL = {"AT", "TA", "CG", "GC"}


def main():
    tier = pd.read_csv(os.path.join(RES, "evidence_tiers_v2.csv"), encoding="utf-8-sig")
    tier.columns = [c.strip() for c in tier.columns]
    genes = sorted(tier.loc[tier["tier"].isin(["A_multi", "B_strong"]), "gene"].unique())
    iv = pd.read_csv(os.path.join(DATA, "ivs.tsv"), sep="\t")
    pairs = iv[iv["gene"].isin(genes)][["gene", "snp", "ea", "z"]].drop_duplicates(["gene", "snp"])

    vcache = json.load(open(os.path.join(TMP, "ensembl_var_cache.json"), encoding="utf-8"))
    gcache = json.load(open(os.path.join(TMP, "mygene_cache.json"), encoding="utf-8"))
    ensg2sym = {v: k for k, v in gcache.items() if v}

    key2pair = {}
    for r in pairs.itertuples():
        m = vcache.get(r.snp)
        if m:
            key2pair[(str(m["chr"]), str(m["pos"]), r.gene)] = r
    print(f"目标 (基因,SNP) 对 {len(pairs)} 个，坐标可用 {len(key2pair)} 个", flush=True)

    files = sorted(glob.glob(os.path.join(GDIR, "*.v8.signif_variant_gene_pairs.txt.gz")))
    rows = []
    for f in files:
        tissue = os.path.basename(f).split(".v8")[0]
        with gzip.open(f, "rt") as fh:
            next(fh)
            for line in fh:
                x = line.rstrip("\n").split("\t")
                q = x[0].split("_")
                ensg = x[1].split(".")[0]
                sym = ensg2sym.get(ensg)
                if sym is None:
                    continue
                k = (q[0].replace("chr", ""), q[1], sym)
                if k not in key2pair:
                    continue
                # 列序：0 variant_id 1 gene_id 2 tss_distance 3 ma_samples
                #       4 ma_count 5 maf 6 pval_nominal 7 slope 8 slope_se
                rows.append(dict(gene=sym, snp=key2pair[k].snp, tissue=tissue,
                                 p=float(x[6]), slope=float(x[7]),
                                 ref=q[2], alt=q[3]))
        print(f"  {tissue:<40s} 累计 {len(rows)}", flush=True)

    d = pd.DataFrame(rows).drop_duplicates(["gene", "snp", "tissue"])

    meta = pairs.set_index(["gene", "snp"])
    d["ea"] = [meta.loc[(g, s), "ea"] for g, s in zip(d["gene"], d["snp"])]
    d["z"] = [meta.loc[(g, s), "z"] for g, s in zip(d["gene"], d["snp"])]
    # 方向：GTEx slope 相对 ALT，与 eQTLGen 的 ea 直接比对（回文位点单独标注）
    sgn = []
    for r in d.itertuples():
        sgn.append(1 if r.ea == r.alt else (-1 if r.ea == r.ref else 0))
    d["same"] = ((d["slope"] * pd.Series(sgn, index=d.index) * d["z"]) > 0)
    d.to_csv(os.path.join(RES, "gtex_tissue_control.csv"), index=False,
             encoding="utf-8-sig")

    N = len(pairs)
    sig = d[d["p"] < 0.05]
    per = sig.groupby("tissue").agg(n_sig=("snp", "size"), n_same=("same", "sum"))
    per["detect_pct"] = 100 * per["n_sig"] / N
    per["same_pct"] = 100 * per["n_same"] / per["n_sig"]
    per = per.sort_values("n_sig", ascending=False)
    per.to_csv(os.path.join(RES, "gtex_tissue_rates.csv"), encoding="utf-8-sig")

    thy = per.loc["Thyroid"]
    rank = int((per["n_sig"] > thy["n_sig"]).sum() + 1)
    above = list(per.index[:rank - 1])
    lines = [
        f"A/B 级 (基因,SNP) 对 {N} 个；GTEx v8 组织 {len(files)} 个",
        f"甲状腺检出 {int(thy['n_sig'])}/{N}（{thy['detect_pct']:.1f}%），"
        f"方向一致 {int(thy['n_same'])}/{int(thy['n_sig'])}（{thy['same_pct']:.1f}%）",
        f"甲状腺排名 第 {rank}/{len(per)}；检出高于它的组织: "
        + (", ".join(above) if above else "无"),
        f"全体组织中位检出 {per['n_sig'].median():.0f}/{N}"
        f"（{100*per['n_sig'].median()/N:.1f}%）",
        "",
        "检出率最高的 10 个组织:",
    ]
    for t, r in per.head(10).iterrows():
        lines.append(f"  {t:<40s} {int(r['n_sig']):>3d}/{N} ({r['detect_pct']:>4.1f}%)"
                     f"  同向 {r['same_pct']:>4.1f}%")
    lines += ["", "检出率最低的 5 个组织:"]
    for t, r in per.tail(5).iterrows():
        lines.append(f"  {t:<40s} {int(r['n_sig']):>3d}/{N} ({r['detect_pct']:>4.1f}%)"
                     f"  同向 {r['same_pct']:>4.1f}%")
    txt = "\n".join(lines)
    print("\n" + txt)
    open(os.path.join(RES, "gtex_tissue_summary.txt"), "w", encoding="utf-8").write(txt + "\n")


if __name__ == "__main__":
    main()
