# -*- coding: utf-8 -*-
"""
29_thyroid_3source_summary.py
Summarise the tissue-specificity retest of the 70 A/B-tier instruments against three thyroid
eQTL sources, writing a comparison table plus a text summary.

Three sources:
  source 1 Zenodo GTEx v8 full nominal (Head et al. 2026, salmon/GENCODE v38)
  source 2 eQTL Catalogue GTEx thyroid (QTD000341, harmonised)
  source 3 GTEx official API (gtex_api_cache.json)
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

    # read the three sources
    d1 = source1_zenodo()
    d2 = source2_eqtlcat()
    d3 = source3_api()

    base = iv[["gene", "snp", "chr", "pos", "ea", "nea", "z"]].copy()
    base = base.merge(d1, on=["gene", "snp"], how="left")
    base = base.merge(d2, on=["gene", "snp"], how="left")
    base = base.merge(d3, on=["gene", "snp"], how="left")

    out = os.path.join(RES, "thyroid_3source_validation.csv")
    base.to_csv(out, index=False, encoding="utf-8-sig")
    print(f"written {out}", flush=True)

    # summary statistics
    print("\n" + "=" * 60)
    print("three-source comparison (70 instruments; uniform definition: each SNP against its own gene, p<0.05)")
    print("=" * 60)
    for name, sigc, samec, hasc in [
        ("source 1 Zenodo (salmon)", "s1_sig", "s1_same", "s1_p"),
        ("source 2 eQTL Catalogue", "s2_sig", "s2_same", "s2_p"),
        ("source 3 GTEx API", "s3_sig", "s3_same", "s3_p"),
    ]:
        has = base[hasc].notna()
        sig = base[sigc].fillna(False)
        sig_has = has & sig
        same_all = base[samec][has].sum()
        same_all_denom = base[samec][has].notna().sum()
        same_sig = base[samec][sig_has].sum()
        same_sig_denom = base[samec][sig_has].notna().sum()
        print(f"\n{name}:")
        print(f"  matched (has a record): {int(has.sum())}/70")
        print(f"  significant (p<0.05): {int(sig.sum())}/70  (of matched {int(sig.sum())}/{int(has.sum())}={100*int(sig.sum())/max(int(has.sum()),1):.0f}%)")
        print(f"  same direction (all matched): {int(same_all)}/{int(same_all_denom)}")
        print(f"  same direction (significant only): {int(same_sig)}/{int(same_sig_denom)}")

    # core SNPs significant in all three sources
    print("\n" + "=" * 60)
    print("SNPs significant and directionally concordant in all three sources (the strongest evidence)")
    print("=" * 60)
    all_sig = base[base[["s1_sig", "s2_sig", "s3_sig"]].fillna(False).all(axis=1)]
    print(f"significant in all three: {len(all_sig)}")
    if len(all_sig):
        print(all_sig[["gene", "snp", "s1_p", "s2_p", "s3_p"]].to_string(index=False))

    # significant in whole blood but not in thyroid (the tissue-mismatch evidence)
    print("\n" + "=" * 60)
    print("SNPs not significant in thyroid (tissue mismatch)")
    print("=" * 60)
    for name, sigc, hasc in [
        ("source 1 Zenodo", "s1_sig", "s1_p"),
        ("source 2 eQTL Catalogue", "s2_sig", "s2_p"),
        ("source 3 GTEx API", "s3_sig", "s3_p"),
    ]:
        has = base[hasc].notna()
        notsig = (has & ~base[sigc].fillna(False)).sum()
        nomatch = (~has).sum()
        print(f"{name}: has a record but not significant {int(notsig)}, no record {int(nomatch)}")


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
        # keep only records for this gene in thyroid
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
