# -*- coding: utf-8 -*-
"""
28_thyroid_3source_validate.py
对 A/B 级 65 个基因的 70 个工具变量 SNP，在三个甲状腺 eQTL 数据源中做组织特异性重验证：
  来源1  Zenodo GTEx v8 甲状腺 full nominal (Head et al. 2026)
  来源2  eQTL Catalogue GTEx 甲状腺 (QTD000341)
  来源3  本地 GTEx v8 官方 signif_variant_gene_pairs（显著位点）
统一输出：每个 SNP 在三个来源里是否为显著甲状腺 eQTL、beta 方向是否与 eQTLGen 全血一致。
"""
import gzip, os, io, sys, json
import pandas as pd
import numpy as np

BASE = "D:/NewBeginning/MR1/DATA10"
RES = os.path.join(BASE, "results")
TMP = os.path.join(BASE, "tmp/eqtl_sources")

# 互补碱基
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def load_iv():
    iv = pd.read_csv(os.path.join(BASE, "data/ivs.tsv"), sep="\t")
    ab = pd.read_csv(os.path.join(TMP, "ab_genes.csv"))
    iv = iv[iv.gene.isin(set(ab.gene))]
    return iv

def sign_of(ea, ref, alt):
    """判断 eQTL 记录的 effect allele 与 eQTLGen 的 ea 是否同向。
    返回 +1（同向）/ -1（反向）/ 0（无法判断，非回文且不匹配）"""
    if ea == alt:
        return 1
    if ea == ref:
        return -1
    # 回文位点：ea 是 alt 的互补（即 alt 与 ref 互补），方向取决于正链假设
    if ea == COMP.get(alt) and COMP.get(alt) == ref:
        return 1  # 回文，按正链假设同向
    if ea == COMP.get(ref) and COMP.get(ref) == alt:
        return -1
    return 0

# ---------- 来源3：GTEx 官方 API 缓存 ----------
def source3_gtx_signif():
    """用 tmp/gtex_api_cache.json。每个键 gene|snp 存该 SNP 在甲状腺关联的所有基因。
    提取 geneSymbol == 该基因 的记录，取最小 pValue + 对应 nes。"""
    iv = load_iv()
    cache_path = os.path.join(BASE, "tmp/gtex_api_cache.json")
    if not os.path.exists(cache_path):
        return None
    cache = json.load(open(cache_path))
    out = []
    for r in iv.itertuples():
        recs = cache.get(f"{r.gene}|{r.snp}") or []
        # 筛选 geneSymbol == 该基因 的记录
        same = [x for x in recs if x.get("geneSymbol") == r.gene]
        if same:
            best = min(same, key=lambda x: x.get("pValue", 1))
            p = best.get("pValue")
            nes = best.get("nes")
            vid = best.get("variantId", "")
            q = vid.split("_")
            ref = q[2] if len(q) > 2 else None
            alt = q[3] if len(q) > 3 else None
            sgn = sign_of(r.ea, ref, alt) if (ref and alt) else 1
            aligned = nes * sgn if nes is not None else np.nan
            out.append(dict(gene=r.gene, snp=r.snp, ea=r.ea, z=r.z,
                            s3_sig=(p < 5e-8), s3_p=p, s3_nes=nes,
                            s3_same=(aligned * r.z > 0) if not np.isnan(aligned) else np.nan))
        else:
            # 该基因在甲状腺无记录（该 SNP 不是该基因的 eQTL，或未达阈值）
            out.append(dict(gene=r.gene, snp=r.snp, ea=r.ea, z=r.z,
                            s3_sig=False, s3_p=np.nan, s3_nes=np.nan,
                            s3_same=np.nan))
    return pd.DataFrame(out)

# ---------- 来源2：eQTL Catalogue ----------
def source2_eqtlcat():
    """QTD000341.all.tsv.gz 全 nominal；effect allele = ALT
    列: molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se ... gene_id ... rsid"""
    iv = load_iv()
    path = os.path.join(TMP, "QTD000341_thyroid.all.tsv.gz")
    if not os.path.exists(path) or os.path.getsize(path) < 1e8:
        return None
    # 建索引：按 (gene_id 或 gene symbol, pos)
    # gene 列是 gene symbol；eQTL Catalogue 用 gene_id=ENSG，但 molecular_trait_id 也是 ENSG。
    # 我们用 rsid 或 chr:pos 匹配。
    iv_rs = dict(zip(iv.snp, iv.itertuples(index=False)))
    iv_pos = {}
    for r in iv.itertuples():
        iv_pos[(int(r.chr), int(r.pos))] = r
    # 需要把 A/B 基因的 symbol -> ENSG（精确匹配该基因）
    try:
        mygene = json.load(open(os.path.join(BASE, "tmp/mygene_cache.json")))
    except Exception:
        mygene = {}
    ensg_of = {r.gene: mygene.get(r.gene) for r in iv.itertuples()}
    iv_rs = dict(zip(iv.snp, iv.itertuples(index=False)))
    hits = {}
    for chunk in pd.read_csv(path, sep="\t", chunksize=2_000_000,
                             usecols=["chromosome", "position", "ref", "alt", "rsid",
                                      "pvalue", "beta", "gene_id", "molecular_trait_id"],
                             dtype={"rsid": str}, low_memory=False):
        m = chunk[chunk.rsid.isin(set(iv.snp))]
        for r in m.itertuples():
            ivr = iv_rs.get(r.rsid)
            if ivr is None:
                continue
            # 精确匹配该基因（gene_id 或 molecular_trait_id 都是 ENSG）
            ensg = ensg_of.get(ivr.gene)
            if ensg and r.gene_id != ensg and r.molecular_trait_id != ensg:
                continue
            key = r.rsid
            if key not in hits or r.pvalue < hits[key]["p"]:
                sgn = sign_of(ivr.ea, r.ref, r.alt)
                hits[key] = dict(gene=ivr.gene, snp=ivr.snp, ea=ivr.ea, z=ivr.z,
                                 ref=r.ref, alt=r.alt, p=r.pvalue, beta=r.beta,
                                 sgn=sgn)
    out = []
    for r in iv.itertuples():
        if r.snp in hits:
            h = hits[r.snp]
            aligned = h["beta"] * h["sgn"]
            out.append(dict(gene=r.gene, snp=r.snp, ea=r.ea, z=r.z,
                            s2_sig=(h["p"] < 5e-8), s2_p=h["p"],
                            s2_aligned_slope=aligned,
                            s2_same=(aligned * r.z > 0)))
        else:
            out.append(dict(gene=r.gene, snp=r.snp, ea=r.ea, z=r.z,
                            s2_sig=False, s2_p=np.nan,
                            s2_aligned_slope=np.nan, s2_same=np.nan))
    return pd.DataFrame(out)

# ---------- 来源1：Zenodo full nominal ----------
def source1_zenodo():
    """解压 Thyroid.tar 后，取 ge/Salmon 或 ge/STAR 的 nominal 文件。
    Zenodo 列: phe_id phe_chr phe_from phe_to phe_strd n_var_in_cis dist_phe_var
                var_id var_chr var_from var_to nom_pval r_squared slope best_hit
    var_id 形如 chrX_POS_REF_ALT 或 chrX_POS_REF_ALT_b38
    phe_id 是 ENSG 编号，需映射 gene symbol。"""
    iv = load_iv()
    tar_dir = os.path.join(TMP, "thyroid_zenodo")
    # 找一个 nominal ge 文件
    cand = []
    if os.path.isdir(tar_dir):
        for root, _, files in os.walk(tar_dir):
            for fn in files:
                if fn.endswith(".txt.gz") and "ge" in fn.lower():
                    cand.append(os.path.join(root, fn))
    if not cand:
        return None
    # 用第一个 ge 文件（可指定 Salmon 或 STAR）
    path = cand[0]
    # Zenodo 无表头，**空格分隔**，15 列：
    # 0 phe_id(ENSG) 1 phe_chr 2 phe_from 3 phe_to 4 phe_strd 5 n_var_in_cis
    # 6 dist_phe_var 7 rsid 8 var_chr 9 var_from 10 var_to
    # 11 nom_pval 12 r_squared 13 slope 14 best_hit
    iv_rs = {}
    for r in iv.itertuples():
        iv_rs[r.snp] = r
    # 基因 symbol -> ENSG（精确匹配该基因的 eQTL，而非该 SNP 的所有基因）
    try:
        mygene = json.load(open(os.path.join(BASE, "tmp/mygene_cache.json")))
    except Exception:
        mygene = {}
    # 反向：ENSG -> symbol（因为 Zenodo phe_id 是 ENSG）
    ensg_of = {r.gene: mygene.get(r.gene) for r in iv.itertuples()}
    hits = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            x = line.split()  # 空格分隔
            if len(x) < 14:
                continue
            phe_id = x[0]      # ENSG
            rsid = x[7]
            ivr = iv_rs.get(rsid)
            if ivr is None:
                continue
            # 只保留该 SNP 对应"该基因"的记录
            if ensg_of.get(ivr.gene) != phe_id:
                continue
            pval = float(x[11]); slope = float(x[13])
            if rsid not in hits or pval < hits[rsid]["p"]:
                hits[rsid] = dict(gene=ivr.gene, snp=ivr.snp, ea=ivr.ea, z=ivr.z,
                                  p=pval, slope=slope)
    out = []
    # Zenodo QTLtools slope 对应 ALT allele（hg38 reference）；用缓存 ref/alt 对齐到 eQTLGen ea
    try:
        cache = json.load(open(os.path.join(BASE, "tmp/ensembl_var_cache.json")))
    except Exception:
        cache = {}
    for r in iv.itertuples():
        if r.snp in hits:
            h = hits[r.snp]
            varmeta = cache.get(r.snp) or {}
            ref = varmeta.get("ref"); alt = varmeta.get("alt")
            sgn = 1
            if ref and alt:
                sgn = sign_of(r.ea, ref, alt)  # 假设 Zenodo effect allele = ALT
            aligned = h["slope"] * sgn
            out.append(dict(gene=r.gene, snp=r.snp, ea=r.ea, z=r.z,
                            s1_sig=(h["p"] < 5e-8), s1_p=h["p"],
                            s1_aligned_slope=aligned, s1_same=(aligned * r.z > 0)))
        else:
            out.append(dict(gene=r.gene, snp=r.snp, ea=r.ea, z=r.z,
                            s1_sig=False, s1_p=np.nan,
                            s1_aligned_slope=np.nan, s1_same=np.nan))
    return pd.DataFrame(out)

def main():
    iv = load_iv()
    print(f"A/B 级基因 {iv.gene.nunique()} 个 / 工具变量 SNP {len(iv)} 个", flush=True)

    base = iv[["gene", "snp", "chr", "pos", "ea", "nea", "z"]].copy()

    # 来源3（本地，必定有）
    s3 = source3_gtx_signif()
    base = base.merge(s3, on=["gene", "snp", "ea", "z"], how="left")
    print(f"来源3 本地GTEx显著位点：完成", flush=True)

    # 来源2（eQTL Catalogue）
    s2 = source2_eqtlcat()
    if s2 is not None:
        base = base.merge(s2, on=["gene", "snp", "ea", "z"], how="left")
        print(f"来源2 eQTL Catalogue：完成", flush=True)
    else:
        print("来源2 eQTL Catalogue：文件未就绪，跳过", flush=True)

    # 来源1（Zenodo）
    s1 = source1_zenodo()
    if s1 is not None:
        base = base.merge(s1, on=["gene", "snp", "ea", "z"], how="left")
        print(f"来源1 Zenodo：完成", flush=True)
    else:
        print("来源1 Zenodo：未解压，跳过", flush=True)

    out_path = os.path.join(RES, "thyroid_3source_validation.csv")
    base.to_csv(out_path, index=False, encoding="utf-8-sig")
    print(f"已写 {out_path}", flush=True)

    # 汇总
    for src, sigc, samec in [("s1", "s1_sig", "s1_same"),
                              ("s2", "s2_sig", "s2_same"),
                              ("s3", "s3_sig", "s3_same")]:
        if sigc in base:
            sig = base[sigc].sum()
            same = base[samec].sum() if samec in base else np.nan
            print(f"  {src}: 显著 {int(sig)}/{len(base)} 同向 {int(same) if not np.isnan(same) else 'NA'}", flush=True)

if __name__ == "__main__":
    main()
