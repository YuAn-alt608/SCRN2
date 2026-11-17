# -*- coding: utf-8 -*-
"""
Step 17: 用 GTEx v8 甲状腺组织 eQTL 验证 A/B 级命中

背景：全部证据都建立在"全血 eQTL + 甲状腺疾病"上，存在组织错配。本步检查
A/B 级命中的工具变量 SNP 在**甲状腺组织**里是否为同一基因的 eQTL，且方向是否
与全血一致。GTEx 甲状腺 n=574，功效有限，定位为"方向验证"而非独立复现。

注意坐标系：GTEx v8 为 GRCh38，本项目 eQTLGen 为 GRCh37，需用 Ensembl 换算。

输出:
  results/gtex_thyroid_validation.csv
  results/gtex_summary.txt
"""
import os, json, time, gzip, subprocess, urllib.request, urllib.parse
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
DATA = os.path.join(BASE, "data")
TMP = os.path.join(BASE, "tmp")
GTEX = os.path.join(TMP, "gtex", "GTEx_Analysis_v8_eQTL",
                    "Thyroid.v8.signif_variant_gene_pairs.txt.gz")
CACHE = os.path.join(TMP, "ensembl_var_cache.json")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


# NCBI RefSeq 染色体登录号 -> 染色体名
ACC2CHR = {f"NC_{i:06d}": str(i) for i in range(1, 23)}
ACC2CHR.update({"NC_000023": "X", "NC_000024": "Y", "NC_000017": "17"})


def ensembl_pos(rsid, cache):
    """rsID -> GRCh38 坐标。

    原先用 Ensembl REST，但该服务对本机 IP 稳定限流（实测连续 20 s 超时），
    因此改用 NCBI Variation Services（限速 3 req/s，无需 key）。
    SPDI 的 position 是 0-based，转 1-based 需 +1。
    """
    if rsid in cache:
        return cache[rsid]
    out = None
    if not rsid.startswith("rs"):
        cache[rsid] = None
        return None
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid[2:]}"
    for _ in range(2):
        try:
            p = subprocess.run(["curl", "-sS", "--max-time", "25", url],
                               capture_output=True, text=True, timeout=35)
            if p.returncode != 0 or not p.stdout.strip():
                time.sleep(0.4)
                continue
            d = json.loads(p.stdout)
            pl = d.get("primary_snapshot_data", {}).get("placements_with_allele", [])
            for x in pl:
                # GRCh38 的 RefSeq 版本号各染色体不同（chr1=.11, chr7=.14, chr12=.12…），
                # 不能用后缀判断，必须看 assembly 字段
                traits = (x.get("placement_annot", {})
                           .get("seq_id_traits_by_assembly", []))
                if not any(str(t.get("assembly_name", "")).startswith("GRCh38")
                           for t in traits):
                    continue
                acc = x.get("seq_id", "")
                base = acc.split(".")[0]
                chrom = ACC2CHR.get(base)
                if chrom is None:
                    continue
                for a in x.get("alleles", []):
                    spdi = a.get("allele", {}).get("spdi", {})
                    ref = spdi.get("deleted_sequence", "")
                    alt = spdi.get("inserted_sequence", "")
                    if len(ref) == 1 and len(alt) == 1 and ref != alt:
                        out = {"chr": chrom, "pos": int(spdi["position"]) + 1,
                               "ref": ref, "alt": alt}
                        break
                if out:
                    break
            break
        except Exception:
            time.sleep(0.4)
    cache[rsid] = out
    return out


def ensg_of(symbol):
    q = f"https://mygene.info/v3/query?q=symbol:{urllib.parse.quote(symbol)}&species=human&fields=ensembl.gene&size=1"
    try:
        p = subprocess.run(["curl", "-sS", "--max-time", "20", q],
                           capture_output=True, text=True, timeout=30)
        d = json.loads(p.stdout)
        e = d["hits"][0].get("ensembl", {}).get("gene")
        if isinstance(e, list):
            e = e[0]
        return e
    except Exception:
        return None


def main():
    t = pd.read_csv(os.path.join(RES, "evidence_tiers_v2.csv"), encoding="utf-8-sig")
    t.columns = [c.strip() for c in t.columns]
    ab = t[t["tier"].isin(["A_multi", "B_strong"])].copy()
    print(f"A/B 级命中 {len(ab)} 条，涉及 {ab['gene'].nunique()} 个基因", flush=True)

    iv = pd.read_csv(os.path.join(DATA, "ivs.tsv"), sep="\t")
    genes = sorted(ab["gene"].unique())

    ensg = {}
    for g in genes:
        ensg[g] = ensg_of(g)
        time.sleep(0.12)
    miss = [g for g, v in ensg.items() if not v]
    print(f"ENSG 映射成功 {len(genes)-len(miss)}/{len(genes)}" + (f"，失败: {miss}" if miss else ""), flush=True)
    ensg2sym = {v.split(".")[0]: k for k, v in ensg.items() if v}

    # 收集 A/B 级基因的全部工具变量 SNP
    sub_iv = iv[iv["gene"].isin(genes)]
    snps = sorted(sub_iv["snp"].unique())
    print(f"工具变量 SNP {len(snps)} 个，开始换算 GRCh38 坐标…", flush=True)

    cache = {}
    if os.path.exists(CACHE):
        cache = json.load(open(CACHE, encoding="utf-8"))
    for i, s in enumerate(snps):
        ensembl_pos(s, cache)
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(snps)}")
    json.dump(cache, open(CACHE, "w", encoding="utf-8"))
    ok = sum(1 for s in snps if cache.get(s))
    print(f"坐标换算成功 {ok}/{len(snps)}", flush=True)

    # 目标键：(chr, pos)
    key2snp = {}
    for s in snps:
        m = cache.get(s)
        if m:
            key2snp[(m["chr"], m["pos"])] = s
    print(f"目标位点 {len(key2snp)} 个，扫描 GTEx 甲状腺…", flush=True)

    frames = []
    for ch in pd.read_csv(GTEX, sep="\t", chunksize=500_000,
                          usecols=["variant_id", "gene_id", "maf",
                                   "pval_nominal", "slope", "slope_se"]):
        ch["chr"] = (ch["variant_id"].str.split("_", expand=True)[0]
                     .str.replace("chr", "", regex=False))
        ch["pos"] = pd.to_numeric(ch["variant_id"].str.split("_", expand=True)[1],
                                  errors="coerce")
        ch["ensg"] = ch["gene_id"].str.split(".").str[0]
        ch = ch[ch["ensg"].isin(ensg2sym)]
        ch = ch[[(c, p) in key2snp for c, p in zip(ch["chr"], ch["pos"])]]
        if len(ch):
            frames.append(ch)
    print(f"命中 GTEx 记录 {sum(len(f) for f in frames) if frames else 0} 行")

    if not frames:
        print("未在甲状腺中检出任何 A/B 级工具变量的 eQTL 记录")
        return

    g = pd.concat(frames, ignore_index=True)
    g["snp"] = [key2snp[(c, p)] for c, p in zip(g["chr"], g["pos"])]
    g["gene"] = g["ensg"].map(ensg2sym)
    g["alt_g"] = g["variant_id"].str.split("_").str[3]

    # 与全血方向比较：eQTLGen 的 z 是相对 ivs.tsv 的 ea
    iv2 = iv.set_index("snp")
    rows = []
    for r in g.itertuples():
        if r.snp not in iv2.index:
            continue
        e = iv2.loc[r.snp]
        if isinstance(e, pd.DataFrame):
            e = e.iloc[0]
        ea = e["ea"]
        # GTEx slope 相对 alt
        if ea == r.alt_g or COMP.get(ea) == r.alt_g:
            sgn = 1
        else:
            sgn = -1
        rows.append({"gene": r.gene, "snp": r.snp, "tissue": "Thyroid",
                     "maf": r.maf, "p_gtex": r.pval_nominal,
                     "slope_gtex": r.slope * sgn,
                     "z_eqtlgen": e["z"],
                     "same_dir": (r.slope * sgn * e["z"]) > 0})
    v = pd.DataFrame(rows)
    v.to_csv(os.path.join(RES, "gtex_thyroid_validation.csv"),
             index=False, encoding="utf-8-sig")

    # 每个基因：是否有甲状腺 eQTL、是否同向
    per = v.groupby("gene").agg(n_snp=("snp", "size"),
                                n_sig=("p_gtex", lambda x: int((x < 0.05).sum())),
                                n_same=("same_dir", "sum"))
    per["validated"] = per["n_sig"] > 0
    per["dir_ok"] = per["n_same"] >= per["n_sig"] / 2

    lines = []
    lines.append(f"A/B 级基因 {len(genes)} 个；坐标换算成功 {ok}/{len(snps)} 个 SNP")
    lines.append(f"GTEx 甲状腺中检出 eQTL 记录的 SNP {v['snp'].nunique()} 个，涉及基因 {v['gene'].nunique()} 个")
    sig = v[v["p_gtex"] < 0.05]
    lines.append(f"其中 p<0.05 的 {len(sig)} 条记录，涉及 {sig['gene'].nunique()} 个基因")
    if len(sig):
        lines.append(f"方向与全血一致 {int(sig['same_dir'].sum())}/{len(sig)} "
                     f"({100*sig['same_dir'].mean():.1f}%)")
    lines.append(f"至少 1 个显著 eQTL 的基因 {int(per['validated'].sum())}/{len(per)}")
    print("\n".join(lines))
    open(os.path.join(RES, "gtex_summary.txt"), "w", encoding="utf-8").write(
        "\n".join(lines) + "\n")
    print("\n已写 results/gtex_thyroid_validation.csv / gtex_summary.txt")


if __name__ == "__main__":
    main()
