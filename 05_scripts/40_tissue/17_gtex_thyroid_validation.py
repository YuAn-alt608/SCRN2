# -*- coding: utf-8 -*-
"""
Step 17: validate the A/B-tier hits against GTEx v8 thyroid eQTLs

Background: every piece of evidence rests on whole-blood eQTLs against thyroid disease, so a
tissue mismatch is possible. This step asks whether the instrument SNPs behind the A/B-tier
hits are eQTLs for the same gene in **thyroid tissue**, and in the same direction as in whole
blood. GTEx thyroid has n=574, so power is limited: this is a direction check, not a replication.

Note the coordinate systems: GTEx v8 is GRCh38 and the eQTLGen data used here are GRCh37, so
positions have to be lifted over.

Output:
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


# NCBI RefSeq chromosome accession -> chromosome name
ACC2CHR = {f"NC_{i:06d}": str(i) for i in range(1, 23)}
ACC2CHR.update({"NC_000023": "X", "NC_000024": "Y", "NC_000017": "17"})


def ensembl_pos(rsid, cache):
    """rsID -> GRCh38 position.

    Ensembl REST was tried first but throttles this IP reliably (20 s of consecutive timeouts),
    so NCBI Variation Services is used instead (3 requests/s, no key). SPDI positions are
    0-based, so converting to 1-based adds one.
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
                # the GRCh38 RefSeq version differs by chromosome (chr1=.11, chr7=.14, chr12=.12...),
                # so the suffix cannot decide it; the assembly field has to be read
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
    print(f"{len(ab)} A/B-tier hits covering {ab['gene'].nunique()} genes", flush=True)

    iv = pd.read_csv(os.path.join(DATA, "ivs.tsv"), sep="\t")
    genes = sorted(ab["gene"].unique())

    ensg = {}
    for g in genes:
        ensg[g] = ensg_of(g)
        time.sleep(0.12)
    miss = [g for g, v in ensg.items() if not v]
    print(f"ENSG mapped {len(genes)-len(miss)}/{len(genes)}" + (f", failed: {miss}" if miss else ""), flush=True)
    ensg2sym = {v.split(".")[0]: k for k, v in ensg.items() if v}

    # collect every instrument SNP of the A/B-tier genes
    sub_iv = iv[iv["gene"].isin(genes)]
    snps = sorted(sub_iv["snp"].unique())
    print(f"{len(snps)} instrument SNPs, lifting over to GRCh38 ...", flush=True)

    cache = {}
    if os.path.exists(CACHE):
        cache = json.load(open(CACHE, encoding="utf-8"))
    for i, s in enumerate(snps):
        ensembl_pos(s, cache)
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(snps)}")
    json.dump(cache, open(CACHE, "w", encoding="utf-8"))
    ok = sum(1 for s in snps if cache.get(s))
    print(f"positions lifted {ok}/{len(snps)}", flush=True)

    # lookup key: (chr, pos)
    key2snp = {}
    for s in snps:
        m = cache.get(s)
        if m:
            key2snp[(m["chr"], m["pos"])] = s
    print(f"{len(key2snp)} target loci, scanning GTEx thyroid ...", flush=True)

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
    print(f"GTEx records matched {sum(len(f) for f in frames) if frames else 0}")

    if not frames:
        print("no eQTL record for any A/B-tier instrument in thyroid")
        return

    g = pd.concat(frames, ignore_index=True)
    g["snp"] = [key2snp[(c, p)] for c, p in zip(g["chr"], g["pos"])]
    g["gene"] = g["ensg"].map(ensg2sym)
    g["alt_g"] = g["variant_id"].str.split("_").str[3]

    # compare with the whole-blood direction: the eQTLGen z is relative to ea in ivs.tsv
    iv2 = iv.set_index("snp")
    rows = []
    for r in g.itertuples():
        if r.snp not in iv2.index:
            continue
        e = iv2.loc[r.snp]
        if isinstance(e, pd.DataFrame):
            e = e.iloc[0]
        ea = e["ea"]
        # the GTEx slope is relative to alt
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

    # per gene: is there a thyroid eQTL, and is it in the same direction
    per = v.groupby("gene").agg(n_snp=("snp", "size"),
                                n_sig=("p_gtex", lambda x: int((x < 0.05).sum())),
                                n_same=("same_dir", "sum"))
    per["validated"] = per["n_sig"] > 0
    per["dir_ok"] = per["n_same"] >= per["n_sig"] / 2

    lines = []
    lines.append(f"{len(genes)} A/B-tier genes; positions lifted for {ok}/{len(snps)} SNPs")
    lines.append(f"{v['snp'].nunique()} SNPs with an eQTL record in GTEx thyroid, covering {v['gene'].nunique()} genes")
    sig = v[v["p_gtex"] < 0.05]
    lines.append(f"of which {len(sig)} records reach p<0.05, covering {sig['gene'].nunique()} genes")
    if len(sig):
        lines.append(f"same direction as whole blood {int(sig['same_dir'].sum())}/{len(sig)} "
                     f"({100*sig['same_dir'].mean():.1f}%)")
    lines.append(f"genes with at least one significant eQTL {int(per['validated'].sum())}/{len(per)}")
    print("\n".join(lines))
    open(os.path.join(RES, "gtex_summary.txt"), "w", encoding="utf-8").write(
        "\n".join(lines) + "\n")
    print("\nwritten: results/gtex_thyroid_validation.csv / gtex_summary.txt")


if __name__ == "__main__":
    main()
