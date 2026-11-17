#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p22_shoulder_attribution.py -- gene attribution for the 45,873,049 shoulder
(previously recorded as an "unassigned hypothyroidism peak")

Background
----------
The hypothyroidism P-value profile across 17q21.32 (10 kb bins) shows a main peak
at 45.80-45.81 (P~1e-33) and, to its right, a plateau at 45.83-45.87 (P pinned at
3e-19) that was previously treated as a second, unassigned hypothyroidism peak.
This script takes every gene within +/-200 kb of that shoulder (23 genes by
Ensembl GRCh37 endpoints) and runs the 17 that are testable through the whole
chain: eQTL lead -> LD with the shoulder and main peaks -> cis-MR -> window
coloc -> verdict.

Coordinate conventions
----------------------
All GRCh37 (pos37). The offset 1,922,634 was checked position by position against
the Ensembl GRCh37 archive endpoints:
  rs6503875 45,915,395; rs7219420 45,800,299; rs1989291 45,807,547; rs4794063 45,804,494
"""
import os
import numpy as np
import pandas as pd

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
T = f"{BASE}/tables"
PEAK = 45_873_049          # shoulder optimum (hypothyroidism P=3.09e-19)
LEADS = {                  # representative variant of each signal
    "shoulder": PEAK,
    "main_peak_rs1989291": 45_807_547,
    "tsh_peak_rs62076102": 45_887_132,
    "scrn2_rs6503875": 45_915_395,
}

ld = pd.read_csv(f"{BASE}/inputs_ld/ld_TSH.tsv.gz", sep="\t", index_col=0)
ld.index = ld.index.astype(int)
ld.columns = ld.columns.astype(int)

genes = pd.read_csv(f"{BASE}/raw/ens_genes_17_peak_b37.csv")
genes["name"] = genes.gene_name.fillna(genes.gene_id)
genes["tss37"] = np.where(genes.strand == 1, genes.start, genes.end)
genes["dist_peak_kb"] = (genes.tss37 - PEAK) / 1000

# alias: AC003665.1 is the former name of SP2-AS1 (same ENSG00000234494); the panel
# carries both names
ALIAS = {"SP2-AS1": "AC003665.1"}

mr = pd.read_csv(f"{T}/04_cisMR_SMR_HEIDI.csv")
co = pd.read_csv(f"{T}/26_peak_neighbour_genes_window_coloc.csv")

rows = []
for _, g in genes.iterrows():
    name = g["name"]
    for tissue in ["Thyroid", "WholeBlood"]:
        f = f"{BASE}/inputs/{tissue}__{name}__TSH.tsv"
        if not os.path.exists(f):
            alt = [k for k, v in ALIAS.items() if v == name]
            if alt:
                f = f"{BASE}/inputs/{tissue}__{alt[0]}__TSH.tsv"
            if not os.path.exists(f):
                continue
        d = pd.read_csv(f, sep="\t", usecols=["pos37", "p_eqtl", "beta_eqtl"])
        row = dict(gene=name, gene_id=g.gene_id, tissue=tissue, tss37=int(g.tss37),
                   dist_shoulder_kb=round(g.dist_peak_kb, 1))
        lead = d.loc[d.p_eqtl.idxmin()]
        row.update(eQTL_lead=int(lead.pos37), eQTL_lead_p=lead.p_eqtl,
                   dist_shoulder_eqtl_lead_kb=round((int(lead.pos37) - PEAK) / 1000, 1),
                   eQTL_at_peak=float(d.loc[d.pos37 == PEAK, "p_eqtl"].iloc[0])
                   if (d.pos37 == PEAK).any() else np.nan)
        # LD with the lead variant of each signal
        for k, v in LEADS.items():
            if v in ld.index and int(lead.pos37) in ld.index:
                row[f"r_{k}"] = round(float(ld.loc[int(lead.pos37), v]), 3)
            else:
                row[f"r_{k}"] = np.nan
        # cis-MR
        for oc in ["HYPO", "TSH"]:
            q = mr[(mr.gene == name) & (mr.tissue == tissue) & (mr.outcome == oc)]
            if not len(q) and name in ALIAS:
                q = mr[(mr.gene == ALIAS[name]) & (mr.tissue == tissue) & (mr.outcome == oc)]
            if len(q):
                row[f"MR_{oc}_b"] = float(q.b_wald.iloc[0])
                row[f"MR_{oc}_p"] = float(q.p_wald.iloc[0])
                row[f"MR_{oc}_F"] = float(q.F_lead.iloc[0])
                row[f"MR_{oc}_niv"] = int(q.n_iv.iloc[0])
        # window coloc (peak +/-50 kb)
        q2 = co[(co.gene == name) & (co.tissue == tissue) &
                (co.outcome == "HYPO") & (co.window.str.contains("50 kb"))]
        row["coloc_H3_HYPO"] = float(q2.H3.iloc[0]) if len(q2) else np.nan
        row["coloc_H4_HYPO"] = float(q2.H4.iloc[0]) if len(q2) else np.nan
        rows.append(row)

s = pd.DataFrame(rows)


s["rank_within_tissue_at_locus"] = (s.groupby("tissue")["eQTL_at_peak"]
                                       .rank(method="min").astype("Int64"))


def verdict(r):
    if pd.isna(r.get("eQTL_lead_p")):
        return "no eQTL data"
    at_peak = r.get("eQTL_at_peak", np.nan)
    lead_dist = abs(r.get("dist_shoulder_eqtl_lead_kb", np.nan))
    h4 = r.get("coloc_H4_HYPO", np.nan)
    h3 = r.get("coloc_H3_HYPO", np.nan)

    # only the top eQTL of that gene in that tissue counts as locus-level overlap;
    # otherwise any strongly associated gene looks "significant" at this position
    # through haplotype linkage and cannot be told apart (11 whole-blood genes reach
    # P<1e-5 in practice)
    if r.get("rank_within_tissue_at_locus", 99) == 1 and at_peak < 1e-3 and lead_dist > 20:
        return ("locus-level overlap (this variant is itself an eQTL for the gene, P=%.1e), "
                "but the main eQTL lies %+.0f kb away and the window coloc gives H3=%.2f "
                "(the single-causal-variant assumption is violated) -> candidate"
                % (at_peak, r["dist_shoulder_eqtl_lead_kb"], h3 if not pd.isna(h3) else np.nan))
    if not pd.isna(h4) and h4 > 0.75:
        return "shared (coloc H4=%.2f)" % h4
    base = ("not shared (H3=%.2f)" % h3) if (not pd.isna(h3) and h3 > 0.75) else "undetermined"
    if not pd.isna(r.get("MR_HYPO_p", np.nan)) and r.MR_HYPO_p < 0.01:
        base += "; but the hypothyroidism cis-MR is significant (P=%.1e)" % r.MR_HYPO_p
    return base


s["verdict"] = [verdict(r) for _, r in s.iterrows()]
s = s.sort_values(["tissue", "coloc_H4_HYPO"], ascending=[True, False])
s.to_csv(f"{T}/27_shoulder_gene_attribution.csv", index=False, encoding="utf-8-sig")

print(f"gene-tissue rows: {len(s)} | genes: {s.gene.nunique()}")
print("\n=== hypothyroidism: local window coloc (peak +/-50 kb), ranked by H4 ===")
c = [x for x in ["gene", "tissue", "dist_shoulder_kb", "eQTL_lead", "eQTL_lead_p", "r_shoulder",
                 "MR_HYPO_b", "MR_HYPO_p", "coloc_H3_HYPO", "coloc_H4_HYPO"] if x in s.columns]
print(s[c].sort_values("coloc_H4_HYPO", ascending=False).head(14).round(3).to_string(index=False))

print("\n=== rows with a hypothyroidism cis-MR signal (P<0.01) ===")
q = s[(s.MR_HYPO_p < 0.01)].sort_values("MR_HYPO_p")
print(q[["gene", "tissue", "eQTL_lead", "r_shoulder", "r_main_peak_rs1989291", "r_scrn2_rs6503875",
         "MR_HYPO_b", "MR_HYPO_p", "MR_HYPO_F", "MR_HYPO_niv"]].round(4).to_string(index=False))
