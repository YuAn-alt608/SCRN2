# -*- coding: utf-8 -*-
"""
29_thyroid_3source_summary.py
汇总三个甲状腺 eQTL 数据源对 A/B 级 70 个工具变量的组织特异性重验证，
输出对比表 + 文本汇总。

三个来源：
  来源1 Zenodo GTEx v8 full nominal (Head et al. 2026, salmon/GENCODE v38)
  来源2 eQTL Catalogue GTEx 甲状腺 (QTD000341, 标准化)
  来源3 GTEx 官方 API (gtex_api_cache.json)
"""
import os, json
import pandas as pd
import numpy as np

BASE = "D:/NewBeginning/MR1/DATA10"
RES = os.path.join(BASE, "results")

def main():
    iv = pd.read_csv(os.path.join(BASE, "data/ivs.tsv"), sep="\t")
    ab = pd.read_csv(os.path.join(BASE, "tmp/eqtl_sources/ab_genes.csv"))
    iv = iv[iv.gene.isin(set(ab.gene))]

    # 读三来源结果
    d1 = source1_zenodo()
    d2 = source2_eqtlcat()
    d3 = source3_api()

    base = iv[["gene", "snp", "chr", "pos", "ea", "nea", "z"]].copy()
    base = base.merge(d1, on=["gene", "snp"], how="left")
    base = base.merge(d2, on=["gene", "snp"], how="left")
    base = base.merge(d3, on=["gene", "snp"], how="left")

    out = os.path.join(RES, "thyroid_3source_validation.csv")
    base.to_csv(out, index=False, encoding="utf-8-sig")
    print(f"已写 {out}", flush=True)

    # 汇总统计
    print("\n" + "=" * 60)
    print("三来源对比汇总（70 个工具变量，统一口径：该 SNP 对该基因，显著阈值 p<0.05）")
    print("=" * 60)
    for name, sigc, samec, hasc in [
        ("来源1 Zenodo(salmon)", "s1_sig", "s1_same", "s1_p"),
        ("来源2 eQTL Catalogue", "s2_sig", "s2_same", "s2_p"),
        ("来源3 GTEx官方API", "s3_sig", "s3_same", "s3_p"),
    ]:
        has = base[hasc].notna()
        sig = base[sigc].fillna(False)
        sig_has = has & sig
        same_all = base[samec][has].sum()
        same_all_denom = base[samec][has].notna().sum()
        same_sig = base[samec][sig_has].sum()
        same_sig_denom = base[samec][sig_has].notna().sum()
        print(f"\n{name}:")
        print(f"  可匹配(有记录): {int(has.sum())}/70")
        print(f"  显著(p<0.05): {int(sig.sum())}/70  (占可匹配 {int(sig.sum())}/{int(has.sum())}={100*int(sig.sum())/max(int(has.sum()),1):.0f}%)")
        print(f"  方向一致(全部可匹配): {int(same_all)}/{int(same_all_denom)}")
        print(f"  方向一致(仅显著子集): {int(same_sig)}/{int(same_sig_denom)}")

    # 三来源一致显著的核心 SNP
    print("\n" + "=" * 60)
    print("三来源都显著且方向一致的 SNP（最强组织验证证据）")
    print("=" * 60)
    all_sig = base[base[["s1_sig", "s2_sig", "s3_sig"]].fillna(False).all(axis=1)]
    print(f"三来源都显著: {len(all_sig)} 个")
    if len(all_sig):
        print(all_sig[["gene", "snp", "s1_p", "s2_p", "s3_p"]].to_string(index=False))

    # 只有全血显著、甲状腺不显著的（组织错配证据）
    print("\n" + "=" * 60)
    print("甲状腺里不显著（组织错配）的 SNP 数")
    print("=" * 60)
    for name, sigc, hasc in [
        ("来源1 Zenodo", "s1_sig", "s1_p"),
        ("来源2 eQTL Catalogue", "s2_sig", "s2_p"),
        ("来源3 GTEx官方API", "s3_sig", "s3_p"),
    ]:
        has = base[hasc].notna()
        notsig = (has & ~base[sigc].fillna(False)).sum()
        nomatch = (~has).sum()
        print(f"{name}: 有记录但不显著 {int(notsig)}，无记录 {int(nomatch)}")


def source1_zenodo():
    iv = _load_iv()
    path = os.path.join(BASE, "tmp/eqtl_sources/thyroid_zenodo/Thyroid",
                        "gene_qtls_nominal_Thyroid.GENCODE_v38.salmon.txt.gz")
    if not os.path.exists(path):
        return pd.DataFrame(dict(gene=iv.gene, snp=iv.snp,
                                 s1_sig=False, s1_p=np.nan, s1_same=np.nan))
    mygene = _mygene()
    ensg_of = {r.gene: mygene.get(r.gene) for r in iv.itertuples()}
    iv_rs = {r.snp: r for r in iv.itertuples()}
    cache = _refalt_cache()
    hits = {}
    import gzip
    with gzip.open(path, "rt") as f:
        for line in f:
            x = line.split()
            if len(x) < 14:
                continue
            phe = x[0]; rsid = x[7]
            ivr = iv_rs.get(rsid)
            if ivr is None:
                continue
            if ensg_of.get(ivr.gene) != phe:
                continue
            p = float(x[11]); slope = float(x[13])
            if rsid not in hits or p < hits[rsid]["p"]:
                hits[rsid] = dict(p=p, slope=slope)
    rows = []
    for r in iv.itertuples():
        if r.snp in hits:
            h = hits[r.snp]
            vm = cache.get(r.snp) or {}
            ref, alt = vm.get("ref"), vm.get("alt")
            sgn = _sign_of(r.ea, ref, alt) if (ref and alt) else 1
            aligned = h["slope"] * sgn
            rows.append(dict(gene=r.gene, snp=r.snp, s1_sig=h["p"] < 0.05,
                             s1_p=h["p"], s1_same=(aligned * r.z > 0)))
        else:
            rows.append(dict(gene=r.gene, snp=r.snp, s1_sig=False,
                             s1_p=np.nan, s1_same=np.nan))
    return pd.DataFrame(rows)


def source2_eqtlcat():
    iv = _load_iv()
    path = os.path.join(BASE, "tmp/eqtl_sources/QTD000341_thyroid.all.tsv.gz")
    if not os.path.exists(path) or os.path.getsize(path) < 3_000_000_000:
        return pd.DataFrame(dict(gene=iv.gene, snp=iv.snp,
                                 s2_sig=False, s2_p=np.nan, s2_same=np.nan))
    mygene = _mygene()
    ensg_of = {r.gene: mygene.get(r.gene) for r in iv.itertuples()}
    iv_rs = {r.snp: r for r in iv.itertuples()}
    hits = {}
    for chunk in pd.read_csv(path, sep="\t", chunksize=2_000_000,
                             usecols=["ref", "alt", "rsid", "pvalue", "beta",
                                      "gene_id", "molecular_trait_id"],
                             dtype={"rsid": str}, low_memory=False):
        m = chunk[chunk.rsid.isin(set(iv.snp))]
        for r in m.itertuples():
            ivr = iv_rs.get(r.rsid)
            if ivr is None:
                continue
            ensg = ensg_of.get(ivr.gene)
            if ensg and r.gene_id != ensg and r.molecular_trait_id != ensg:
                continue
            key = r.rsid
            if key not in hits or r.pvalue < hits[key]["p"]:
                sgn = _sign_of(ivr.ea, r.ref, r.alt)
                hits[key] = dict(p=r.pvalue, beta=r.beta, sgn=sgn)
    rows = []
    for r in iv.itertuples():
        if r.snp in hits:
            h = hits[r.snp]
            aligned = h["beta"] * h["sgn"]
            rows.append(dict(gene=r.gene, snp=r.snp, s2_sig=h["p"] < 0.05,
                             s2_p=h["p"], s2_same=(aligned * r.z > 0)))
        else:
            rows.append(dict(gene=r.gene, snp=r.snp, s2_sig=False,
                             s2_p=np.nan, s2_same=np.nan))
    return pd.DataFrame(rows)


def source3_api():
    iv = _load_iv()
    cache_path = os.path.join(BASE, "tmp/gtex_api_cache.json")
    if not os.path.exists(cache_path):
        return pd.DataFrame(dict(gene=iv.gene, snp=iv.snp,
                                 s3_sig=False, s3_p=np.nan, s3_same=np.nan))
    cache = json.load(open(cache_path))
    rows = []
    for r in iv.itertuples():
        recs = cache.get(f"{r.gene}|{r.snp}") or []
        # 只取"该基因 + 甲状腺组织"的记录
        same = [x for x in recs if x.get("geneSymbol") == r.gene
                and x.get("tissueSiteDetailId") == "Thyroid"]
        if same:
            best = min(same, key=lambda x: x.get("pValue", 1))
            p = best.get("pValue"); nes = best.get("nes")
            vid = best.get("variantId", "")
            q = vid.split("_")
            ref = q[2] if len(q) > 2 else None
            alt = q[3] if len(q) > 3 else None
            sgn = _sign_of(r.ea, ref, alt) if (ref and alt) else 1
            aligned = nes * sgn if nes is not None else np.nan
            same_flag = (aligned * r.z > 0) if not np.isnan(aligned) else np.nan
            rows.append(dict(gene=r.gene, snp=r.snp, s3_sig=p < 0.05,
                             s3_p=p, s3_same=same_flag))
        else:
            rows.append(dict(gene=r.gene, snp=r.snp, s3_sig=False,
                             s3_p=np.nan, s3_same=np.nan))
    return pd.DataFrame(rows)


# ---- helpers ----
def _load_iv():
    iv = pd.read_csv(os.path.join(BASE, "data/ivs.tsv"), sep="\t")
    ab = pd.read_csv(os.path.join(BASE, "tmp/eqtl_sources/ab_genes.csv"))
    return iv[iv.gene.isin(set(ab.gene))]


def _mygene():
    try:
        return json.load(open(os.path.join(BASE, "tmp/mygene_cache.json")))
    except Exception:
        return {}


def _refalt_cache():
    try:
        return json.load(open(os.path.join(BASE, "tmp/ensembl_var_cache.json")))
    except Exception:
        return {}


def _sign_of(ea, ref, alt):
    COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if ea == alt:
        return 1
    if ea == ref:
        return -1
    if ea == COMP.get(alt) and COMP.get(alt) == ref:
        return 1
    if ea == COMP.get(ref) and COMP.get(ref) == alt:
        return -1
    return 0


if __name__ == "__main__":
    main()
