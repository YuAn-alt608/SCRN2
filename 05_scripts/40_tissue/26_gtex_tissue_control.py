# -*- coding: utf-8 -*-
"""
Step 26: GTEx multi-tissue control (local v8, all tissues; same source and definition as section 06)

Question: section 06 looked at thyroid alone, so a 34% detection rate cannot be judged high or
      low. If the same instruments were also detected in a third of liver, muscle and lung
      samples, the "thyroid validation" would not hold.

Approach: for every A/B-tier (gene, instrument) pair, match on GRCh38 position plus ENSG in the
      *.v8.signif_variant_gene_pairs.txt.gz files of the 49 tissues, then compute the detection
      rate (p<0.05) and directional agreement per tissue and see where thyroid ranks.

      Difference from section 06: that section queried thyroid only and used the 64 mappable SNPs
      as the denominator; this step covers all 49 tissues with all 70 (gene, SNP) pairs as the
      denominator, which is the stricter definition.

Input: tmp/gtex/GTEx_Analysis_v8_eQTL/*.signif_variant_gene_pairs.txt.gz
Output:
  results/gtex_tissue_control.csv   long table (gene, snp, tissue, p, slope)
  results/gtex_tissue_rates.csv     detection and agreement rate per tissue
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
    print(f"{len(pairs)} target (gene, SNP) pairs, {len(key2pair)} with usable positions", flush=True)

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
                # column order: 0 variant_id, 1 gene_id, 2 tss_distance, 3 ma_samples
                #       4 ma_count 5 maf 6 pval_nominal 7 slope 8 slope_se
                rows.append(dict(gene=sym, snp=key2pair[k].snp, tissue=tissue,
                                 p=float(x[6]), slope=float(x[7]),
                                 ref=q[2], alt=q[3]))
        print(f"  {tissue:<40s} cumulative {len(rows)}", flush=True)

    d = pd.DataFrame(rows).drop_duplicates(["gene", "snp", "tissue"])

    meta = pairs.set_index(["gene", "snp"])
    d["ea"] = [meta.loc[(g, s), "ea"] for g, s in zip(d["gene"], d["snp"])]
    d["z"] = [meta.loc[(g, s), "z"] for g, s in zip(d["gene"], d["snp"])]
    # direction: the GTEx slope is relative to ALT, compared directly with the eQTLGen ea
    # (palindromic sites are flagged separately)
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
        f"A/B-tier (gene, SNP) pairs: {N}; GTEx v8 tissues: {len(files)}",
        f"detected in thyroid {int(thy['n_sig'])}/{N} ({thy['detect_pct']:.1f}%),"
        f" same direction {int(thy['n_same'])}/{int(thy['n_sig'])} ({thy['same_pct']:.1f}%)",
        f"thyroid rank {rank} of {len(per)}; tissues detecting more: "
        + (", ".join(above) if above else "none"),
        f"median detection across tissues {per['n_sig'].median():.0f}/{N}"
        f"({100*per['n_sig'].median()/N:.1f}%)",
        "",
        "10 tissues with the highest detection rate:",
    ]
    for t, r in per.head(10).iterrows():
        lines.append(f"  {t:<40s} {int(r['n_sig']):>3d}/{N} ({r['detect_pct']:>4.1f}%)"
                     f"  same dir. {r['same_pct']:>4.1f}%")
    lines += ["", "5 tissues with the lowest detection rate:"]
    for t, r in per.tail(5).iterrows():
        lines.append(f"  {t:<40s} {int(r['n_sig']):>3d}/{N} ({r['detect_pct']:>4.1f}%)"
                     f"  same dir. {r['same_pct']:>4.1f}%")
    txt = "\n".join(lines)
    print("\n" + txt)
    open(os.path.join(RES, "gtex_tissue_summary.txt"), "w", encoding="utf-8").write(txt + "\n")


if __name__ == "__main__":
    main()
