# -*- coding: utf-8 -*-
"""
28_thyroid_3source_validate.py
Re-test the 70 instrument SNPs of the 65 A/B-tier genes for tissue specificity against three
thyroid eQTL sources:
  source 1  Zenodo GTEx v8 thyroid full nominal (Head et al. 2026)
  source 2  eQTL Catalogue GTEx thyroid (QTD000341)
  source 3  local GTEx v8 official signif_variant_gene_pairs (significant pairs)
Common output: for every SNP, whether it is a significant thyroid eQTL in each source and
whether the beta direction agrees with whole blood in eQTLGen.
"""
import gzip, os, io, sys, json
import pandas as pd
import numpy as np

BASE = "D:/NewBeginning/MR1/DATA10"
RES = os.path.join(BASE, "results")
TMP = os.path.join(BASE, "tmp/eqtl_sources")

# complementary bases
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def load_iv():
    iv = pd.read_csv(os.path.join(BASE, "data/ivs.tsv"), sep="\t")
    ab = pd.read_csv(os.path.join(TMP, "ab_genes.csv"))
    iv = iv[iv.gene.isin(set(ab.gene))]
    return iv

def sign_of(ea, ref, alt):
    """Decide whether the effect allele of an eQTL record points the same way as the eQTLGen ea.
    Returns +1 (same), -1 (opposite) or 0 (undecidable: not palindromic and no match)."""
    if ea == alt:
        return 1
    if ea == ref:
        return -1
    # palindromic site: ea is the complement of alt (so alt and ref are complementary) and the
    # direction depends on the strand assumption
    if ea == COMP.get(alt) and COMP.get(alt) == ref:
        return 1  # palindromic, same direction under the plus-strand assumption
    if ea == COMP.get(ref) and COMP.get(ref) == alt:
        return -1
    return 0

# ---------- source 3: cached GTEx official API ----------
def source3_gtx_signif():
    """Uses tmp/gtex_api_cache.json. Each key gene|snp holds every gene the SNP is associated with
    in thyroid. Records with geneSymbol == the target gene are kept, taking the smallest pValue
    and its nes."""
    iv = load_iv()
    cache_path = os.path.join(BASE, "tmp/gtex_api_cache.json")
    if not os.path.exists(cache_path):
        return None
    cache = json.load(open(cache_path))
    out = []
    for r in iv.itertuples():
        recs = cache.get(f"{r.gene}|{r.snp}") or []
        # keep records whose geneSymbol is the target gene
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
            # no thyroid record for this gene (the SNP is not its eQTL, or falls below threshold)
            out.append(dict(gene=r.gene, snp=r.snp, ea=r.ea, z=r.z,
                            s3_sig=False, s3_p=np.nan, s3_nes=np.nan,
                            s3_same=np.nan))
    return pd.DataFrame(out)

# ---------- source 2: eQTL Catalogue ----------
def source2_eqtlcat():
    """QTD000341.all.tsv.gz, all nominal pairs; effect allele = ALT.
    Columns: molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se
    ... gene_id ... rsid"""
    iv = load_iv()
    path = os.path.join(TMP, "QTD000341_thyroid.all.tsv.gz")
    if not os.path.exists(path) or os.path.getsize(path) < 1e8:
        return None
    # build an index on (gene_id or gene symbol, pos). The gene column holds a symbol; eQTL
    # Catalogue uses gene_id=ENSG, and molecular_trait_id is an ENSG as well, so matching is done
    # on rsid or chr:pos.
    iv_rs = dict(zip(iv.snp, iv.itertuples(index=False)))
    iv_pos = {}
    for r in iv.itertuples():
        iv_pos[(int(r.chr), int(r.pos))] = r
    # map A/B-tier symbols to ENSG, so the match is to that gene and not to any gene at the locus
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
            # exact match to the gene (both gene_id and molecular_trait_id are ENSG)
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

# ---------- source 1: Zenodo full nominal ----------
def source1_zenodo():
    """After unpacking Thyroid.tar, take the nominal file from ge/Salmon or ge/STAR.
    Zenodo columns: phe_id phe_chr phe_from phe_to phe_strd n_var_in_cis dist_phe_var
                var_id var_chr var_from var_to nom_pval r_squared slope best_hit
    var_id looks like chrX_POS_REF_ALT or chrX_POS_REF_ALT_b38
    phe_id is an ENSG, which has to be mapped to a gene symbol."""
    iv = load_iv()
    tar_dir = os.path.join(TMP, "thyroid_zenodo")
    # pick a nominal ge file
    cand = []
    if os.path.isdir(tar_dir):
        for root, _, files in os.walk(tar_dir):
            for fn in files:
                if fn.endswith(".txt.gz") and "ge" in fn.lower():
                    cand.append(os.path.join(root, fn))
    if not cand:
        return None
    # use the first ge file (Salmon or STAR can be requested)
    path = cand[0]
    # Zenodo has no header, is **space separated**, and has 15 columns:
    # 0 phe_id(ENSG) 1 phe_chr 2 phe_from 3 phe_to 4 phe_strd 5 n_var_in_cis
    # 6 dist_phe_var 7 rsid 8 var_chr 9 var_from 10 var_to
    # 11 nom_pval 12 r_squared 13 slope 14 best_hit
    iv_rs = {}
    for r in iv.itertuples():
        iv_rs[r.snp] = r
    # gene symbol -> ENSG, so the match is to that gene's eQTL and not to every gene at the SNP
    try:
        mygene = json.load(open(os.path.join(BASE, "tmp/mygene_cache.json")))
    except Exception:
        mygene = {}
    # and the reverse, ENSG -> symbol, because the Zenodo phe_id is an ENSG
    ensg_of = {r.gene: mygene.get(r.gene) for r in iv.itertuples()}
    hits = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            x = line.split()  # space separated
            if len(x) < 14:
                continue
            phe_id = x[0]      # ENSG
            rsid = x[7]
            ivr = iv_rs.get(rsid)
            if ivr is None:
                continue
            # keep only the records in which this SNP belongs to this gene
            if ensg_of.get(ivr.gene) != phe_id:
                continue
            pval = float(x[11]); slope = float(x[13])
            if rsid not in hits or pval < hits[rsid]["p"]:
                hits[rsid] = dict(gene=ivr.gene, snp=ivr.snp, ea=ivr.ea, z=ivr.z,
                                  p=pval, slope=slope)
    out = []
    # the Zenodo QTLtools slope refers to the ALT allele against the hg38 reference; the cached
    # ref/alt pair aligns it to the eQTLGen ea
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
                sgn = sign_of(r.ea, ref, alt)  # assuming the Zenodo effect allele is ALT
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
    print(f"{iv.gene.nunique()} A/B-tier genes / {len(iv)} instrument SNPs", flush=True)

    base = iv[["gene", "snp", "chr", "pos", "ea", "nea", "z"]].copy()

    # source 3 (local, always available)
    s3 = source3_gtx_signif()
    base = base.merge(s3, on=["gene", "snp", "ea", "z"], how="left")
    print(f"source 3 local GTEx significant pairs: done", flush=True)

    # source 2 (eQTL Catalogue)
    s2 = source2_eqtlcat()
    if s2 is not None:
        base = base.merge(s2, on=["gene", "snp", "ea", "z"], how="left")
        print(f"source 2 eQTL Catalogue: done", flush=True)
    else:
        print("source 2 eQTL Catalogue: file not present, skipped", flush=True)

    # source 1 (Zenodo)
    s1 = source1_zenodo()
    if s1 is not None:
        base = base.merge(s1, on=["gene", "snp", "ea", "z"], how="left")
        print(f"source 1 Zenodo: done", flush=True)
    else:
        print("source 1 Zenodo: not unpacked, skipped", flush=True)

    out_path = os.path.join(RES, "thyroid_3source_validation.csv")
    base.to_csv(out_path, index=False, encoding="utf-8-sig")
    print(f"written {out_path}", flush=True)

    # summary
    for src, sigc, samec in [("s1", "s1_sig", "s1_same"),
                              ("s2", "s2_sig", "s2_same"),
                              ("s3", "s3_sig", "s3_same")]:
        if sigc in base:
            sig = base[sigc].sum()
            same = base[samec].sum() if samec in base else np.nan
            print(f"  {src}: significant {int(sig)}/{len(base)} same direction {int(same) if not np.isnan(same) else 'NA'}", flush=True)

if __name__ == "__main__":
    main()
