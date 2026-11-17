# -*- coding: utf-8 -*-
"""
Step 24: GTEx 多组织对照 —— 甲状腺的"组织验证"到底有多特异

问题：第 06 节只做了甲状腺一个组织，34% 的检出率无法判断高低。
      如果同样的工具变量在肝、肌肉、胰腺里也有 30% 检出，"甲状腺验证"就不成立。

做法：对 A/B 级全部 70 个（基因, 工具变量）对，向 GTEx v8 查询该变异在所有
      组织的 eQTL 记录（singleTissueEqtl，按 geneSymbol 过滤到目标基因），
      统计每个组织的检出率（p<0.05）与方向一致率，再看甲状腺排在什么位置。

数据源：GTEx Portal API v2，datasetId=gtex_v8。
      已核对：RERE/rs301802 的甲状腺记录为 p=7.04189e-09、nes=0.148，
      与本地 Thyroid.v8.signif_variant_gene_pairs 完全一致，说明 API 与 v8 同源。
      注意：v2 默认返回 gtex_v10，必须显式传 datasetId=gtex_v8。

变异 ID：API 只认 variantId（chr_pos_ref_alt_b38）。用 NCBI 缓存构造，
        失败时尝试 REF/ALT 互换（多等位基因位点两种写法都可能存在）。

输出：
  results/gtex_tissue_control.csv    长表（gene, snp, tissue, p, nes）
  results/gtex_tissue_summary.txt    各组织检出率与甲状腺排名
"""
import os, json, time, subprocess
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
DATA = os.path.join(BASE, "data")
TMP = os.path.join(BASE, "tmp")
CACHE = os.path.join(TMP, "gtex_api_cache.json")
VAR_CACHE = os.path.join(TMP, "ensembl_var_cache.json")
API = "https://gtexportal.org/api/v2/association/singleTissueEqtl"


def query(vid, tries=3):
    url = f"{API}?variantId={vid}&datasetId=gtex_v8&itemsPerPage=1000"
    for _ in range(tries):
        try:
            p = subprocess.run(["curl", "-sS", "--max-time", "30", url],
                               capture_output=True, text=True, timeout=40)
            if p.returncode == 0 and p.stdout.strip():
                return json.loads(p.stdout).get("data", [])
        except Exception:
            pass
        time.sleep(0.5)
    return None


def main():
    tier = pd.read_csv(os.path.join(RES, "evidence_tiers_v2.csv"), encoding="utf-8-sig")
    tier.columns = [c.strip() for c in tier.columns]
    genes = sorted(tier.loc[tier["tier"].isin(["A_multi", "B_strong"]), "gene"].unique())

    iv = pd.read_csv(os.path.join(DATA, "ivs.tsv"), sep="\t")
    pairs = iv[iv["gene"].isin(genes)][["gene", "snp"]].drop_duplicates()
    print(f"A/B 级 {len(genes)} 个基因，(基因,SNP) 对 {len(pairs)} 个", flush=True)

    vcache = json.load(open(VAR_CACHE, encoding="utf-8"))
    cache = json.load(open(CACHE, encoding="utf-8")) if os.path.exists(CACHE) else {}

    rows, fail = [], []
    for i, r in enumerate(pairs.itertuples(), 1):
        key = f"{r.gene}|{r.snp}"
        if key in cache:
            recs = cache[key]
        else:
            m = vcache.get(r.snp)
            recs = None
            if m:
                for ref, alt in [(m["ref"], m["alt"]), (m["alt"], m["ref"])]:
                    vid = f"chr{m['chr']}_{m['pos']}_{ref}_{alt}_b38"
                    recs = query(vid)
                    if recs:
                        break
            cache[key] = recs
            if (i % 10) == 0:
                json.dump(cache, open(CACHE, "w", encoding="utf-8"))
                print(f"  {i}/{len(pairs)}", flush=True)
        if recs is None:
            fail.append(key)
            continue
        for x in recs:
            if x.get("geneSymbol") != r.gene:
                continue
            rows.append(dict(gene=r.gene, snp=r.snp,
                             tissue=x["tissueSiteDetailId"],
                             p=x["pValue"], nes=x["nes"]))
    json.dump(cache, open(CACHE, "w", encoding="utf-8"))

    d = pd.DataFrame(rows).drop_duplicates(["gene", "snp", "tissue"])
    d.to_csv(os.path.join(RES, "gtex_tissue_control.csv"), index=False,
             encoding="utf-8-sig")

    iv2 = iv.drop_duplicates("snp").set_index("snp")
    d["ea"] = d["snp"].map(iv2["ea"])
    d["z"] = d["snp"].map(iv2["z"])
    d["same"] = (d["nes"] * d["z"]) > 0

    n_pairs = len(pairs)
    sig = d[d["p"] < 0.05]
    per = sig.groupby("tissue").agg(n_sig=("snp", "size"),
                                    n_same=("same", "sum"))
    per["detect_pct"] = 100 * per["n_sig"] / n_pairs
    per["same_pct"] = 100 * per["n_same"] / per["n_sig"]
    per = per.sort_values("n_sig", ascending=False)
    per.to_csv(os.path.join(RES, "gtex_tissue_rates.csv"), encoding="utf-8-sig")

    thy = per.loc["Thyroid"] if "Thyroid" in per.index else None
    rank = (per["n_sig"] > thy["n_sig"]).sum() + 1 if thy is not None else None
    n_tis = d["tissue"].nunique()
    med = per["n_sig"].median()

    lines = [
        f"A/B 级 (基因,SNP) 对 {n_pairs} 个；查询失败 {len(fail)} 个" + (f"（{fail[:5]}）" if fail else ""),
        f"GTEx v8 覆盖组织 {n_tis} 个；表中出现显著 eQTL 的组织 {len(per)} 个",
        "",
    ]
    if thy is not None:
        lines += [
            f"甲状腺：检出 {int(thy['n_sig'])}/{n_pairs}（{thy['detect_pct']:.1f}%），"
            f"方向一致 {int(thy['n_same'])}/{int(thy['n_sig'])}（{thy['same_pct']:.1f}%）",
            f"甲状腺检出率排名：第 {rank}/{len(per)}（并列按上界计）",
            f"所有组织的中位检出数：{med:.0f}/{n_pairs}（{100*med/n_pairs:.1f}%）",
            "",
            "检出数最高的 8 个组织:",
        ]
        for t, r in per.head(8).iterrows():
            lines.append(f"  {t:<38s} {int(r['n_sig']):>3d}/{n_pairs} "
                         f"({r['detect_pct']:.0f}%)  同向 {r['same_pct']:.0f}%")
        lines += ["", "检出数最低的 5 个组织:"]
        for t, r in per.tail(5).iterrows():
            lines.append(f"  {t:<38s} {int(r['n_sig']):>3d}/{n_pairs} "
                         f"({r['detect_pct']:.0f}%)  同向 {r['same_pct']:.0f}%")
        above = per[per["n_sig"] > thy["n_sig"]]
        if len(above):
            lines += ["", f"检出高于甲状腺的组织 {len(above)} 个: " + ", ".join(above.index[:12])]
    txt = "\n".join(lines)
    print("\n" + txt)
    open(os.path.join(RES, "gtex_tissue_summary.txt"), "w", encoding="utf-8").write(txt + "\n")


if __name__ == "__main__":
    main()
