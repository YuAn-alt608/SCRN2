#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
12a_coloc_triage.py
-------------------
Before running coloc.susie, triage which pairs actually need susie (no LD required).

Logic: coloc.abf assumes a single causal variant per trait. If the eQTL side or
the outcome side of a gene-outcome pair carries >=2 distance-separated signal
peaks (>250 kb apart), abf can misread "both traits have a signal but the peaks
do not coincide" as H3; such pairs must be re-checked with susie.
Regions with a single peak give the same answer under susie and abf, so they are
not recomputed.

Output: results/susie_triage.csv
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
    """Distance clustering: greedy peak picking in ascending P, dropping SNPs within 250 kb. Returns (n_peaks, lead_positions, lead_ps)."""
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

    # MR significance, used to prioritise the queue
    mr = pd.read_csv(os.path.join(RES, "MR_significant_FDR05_v2.csv"), encoding="utf-8-sig")
    OMAP = {"TSH": "TSH", "Free T4": "T4", "Hypothyroidism": "HYPO", "Hashimoto thyroiditis": "HASH"}
    mr["outcome"] = mr["outcome"].map(OMAP)
    mrkey = mr.groupby(["gene", "outcome"])["pval"].min().rename("p_mr").reset_index()

    # eQTL regions
    eq = pd.read_csv(os.path.join(DATA, "eqtl_regions_all.txt.gz"), sep="\t",
                     usecols=["Pvalue", "SNP", "SNPChr", "SNPPos", "GeneSymbol"],
                     dtype={"Pvalue": float, "SNPPos": int})
    eq = eq.rename(columns={"Pvalue": "p", "SNPPos": "pos", "GeneSymbol": "gene",
                            "SNPChr": "chr"})
    eq["chr"] = eq["chr"].astype(str).str.replace("chr", "", regex=False)
    # the cis window follows the actual eQTL span (eQTLGen uses gene +/-1 Mb, so the TSS cannot serve as centre once the span reaches 1 Mb)
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
    # three situations that need a susie re-check
    t["need_susie"] = (
        (t["multi_signal"] & (t["PP.H3"] > 0.5)) |          # suspected spurious H3
        (t["multi_signal"] & (t["PP.H4"] > 0.5)) |          # shared signal to confirm
        (~t["multi_signal"] & (t["PP.H3"] > 0.5) & (t["PP.H4"] > 0.2))  # single peak but unstable call
    )
    t = t.sort_values(["need_susie", "p_mr"], ascending=[False, True])
    t.to_csv(os.path.join(RES, "susie_triage.csv"), index=False, encoding="utf-8-sig")

    print("\n=== triage result ===")
    print(t.groupby("outcome").agg(n=("gene", "size"), multi=("multi_signal", "sum"),
                                   need=("need_susie", "sum")))
    print("need_susie total:", int(t["need_susie"].sum()))
    print(t[t["need_susie"]].head(15)[["gene", "outcome", "n_sig_eqtl", "n_sig_out", "PP.H3", "PP.H4", "p_mr"]].to_string())


if __name__ == "__main__":
    main()
