#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p01_build_panel.py
==================
Locus gene panel: pair every gene with an eQTL within +/-500 kb, in two tissues
-- thyroid (GTEx v8 / eQTL Catalogue QTD000341) and whole blood (eQTLGen) -- with
the four outcomes (TSH / T4 / HYPO / HASH), producing:

  tables/01_gene_list.csv            all genes in the window (including those without an eQTL) with coordinates
  tables/02_panel_lead.csv           lead eQTL SNP and both effect sets per gene x tissue x outcome
  inputs/{tissue}__{gene}__{oc}.tsv  coloc input (snp, pos37, both beta/se sets, maf)
  inputs_ld/ld_{oc}.tsv.gz           1000G EUR r matrix (indexed by pos37)

Coordinate conventions: pos37 = GRCh37 used by the GWAS window; pos38 = GRCh38 of the eQTL Catalogue variant_id.
Allele conventions: each row carries eqtl_ea / gwas_ea, and all betas are aligned to eqtl_ea.
"""
import os, sys, json
import numpy as np
import pandas as pd

BASE = "D:/NewBeginning/MR1"
OUT = f"{BASE}/DATA21_complete_2026-11-14"
RAW = f"{OUT}/raw"
os.makedirs(f"{OUT}/inputs", exist_ok=True)
os.makedirs(f"{OUT}/inputs_ld", exist_ok=True)
os.makedirs(f"{OUT}/tables", exist_ok=True)

sys.path.insert(0, f"{BASE}/DATA15_thyroid_eQTL/scripts")
from ld_utils import get_region

CENTER37 = 45915395          # rs6503875, GRCh37
HALF = 500_000
N_EQTLGEN = 31684
OCS = ["TSH", "T4", "HYPO", "HASH"]

# ---------------------------------------------------------------- 1. outcome windows
REN = {"base_pair_location": "pos37", "rs_id": "snp", "rsid": "snp",
       "effect_allele": "gwas_ea", "other_allele": "gwas_oa",
       "standard_error": "se_out", "p_value": "p_out",
       "effect_allele_frequency": "eaf_out", "n": "n_out"}
gwas = {}
for oc in OCS:
    d = pd.read_csv(f"{BASE}/DATA15_thyroid_eQTL/results/gwas_win/{oc}.tsv.gz", sep="\t")
    d = d.rename(columns=REN)
    d = d[(d.chromosome.astype(str).str.replace("chr", "") == "17") &
          (d.pos37.between(CENTER37 - HALF, CENTER37 + HALF))].copy()
    if "snp" not in d:                       # the HASH window has no rs_id column
        d["snp"] = "chr17:" + d.pos37.astype(str)
    if "n_out" not in d:
        d["n_out"] = np.nan
    d["snp"] = d.snp.fillna("chr17:" + d.pos37.astype(str)).astype(str)
    d = d.drop_duplicates("pos37")
    gwas[oc] = d[["snp", "pos37", "gwas_ea", "gwas_oa", "beta", "se_out", "p_out", "eaf_out", "n_out"]]
    print(f"[GWAS] {oc}: {len(d)} SNP, pos37 {d.pos37.min()}-{d.pos37.max()}")

# region-wide rsID <-> GRCh37 position map (union of the three windows that carry
# rsIDs), used to place the HASH variants
rs2pos = pd.concat([gwas[o][["snp", "pos37"]] for o in ["TSH", "T4", "HYPO"]])
rs2pos = rs2pos[~rs2pos.snp.str.startswith("chr17:")].drop_duplicates("snp")
print(f"[MAP] rsID->pos37: {len(rs2pos)} entries")

# ---------------------------------------------------------------- 2. LD panel
ld = get_region("17", CENTER37 - HALF, CENTER37 + HALF)
pos_ld = np.asarray(ld["pos"], dtype=int)
Graw = np.asarray(ld["G"], dtype=float)          # (n_snp, n_sample)
G = Graw - Graw.mean(1, keepdims=True)
sd = G.std(1, ddof=0)
keep = sd > 0
pos_ld, G, sd, Graw = pos_ld[keep], G[keep], sd[keep], Graw[keep]
R = (G / sd[:, None]) @ (G / sd[:, None]).T / G.shape[1]
maf = np.minimum(Graw.mean(1), 1 - Graw.mean(1))
maf = np.clip(maf, 1e-4, 0.5)
# keep only the first variant at each position (indels / multi-allelic sites) so
# that pos37 is unique
dup = pd.Series(pos_ld).duplicated().values
if dup.any():
    print(f"[LD] dropped {int(dup.sum())} duplicate-position variants")
    pos_ld, maf = pos_ld[~dup], maf[~dup]
    R = R[np.ix_(~dup, ~dup)]
print(f"[LD] {len(pos_ld)} SNP, {G.shape[1]} samples")

# SNP subset actually used by each outcome (one LD file per outcome, for coloc.susie)
for oc in OCS:
    snps = list(dict.fromkeys([p for p in gwas[oc].pos37 if p in set(pos_ld)]))
    if not snps:
        continue
    idx = pd.Series(np.arange(len(pos_ld)), index=pos_ld).loc[snps].values
    sub = pd.DataFrame(R[np.ix_(idx, idx)], index=snps, columns=snps)
    sub.index.name = "pos37"
    sub.round(4).to_csv(f"{OUT}/inputs_ld/ld_{oc}.tsv.gz", sep="\t", compression="gzip")
    print(f"[LD] {oc}: {len(snps)}-SNP submatrix written")

ld_map = pd.DataFrame({"pos37": pos_ld, "maf1000g": maf})
ld_map.to_csv(f"{OUT}/tables/00_ld_panel_loci.csv", index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------- 3. gene coordinates
eg = pd.read_csv(f"{BASE}/DATA15_thyroid_eQTL/tmp/GTEx_Analysis_v8_eQTL/Thyroid.v8.egenes.txt.gz",
                 sep="\t", usecols=["gene_id", "gene_name", "gene_chr", "gene_start", "gene_end",
                                    "strand", "pval_nominal"])
eg["gene_id"] = eg.gene_id.str.split(".").str[0]
eg = eg.drop_duplicates("gene_id")
# GRCh38 centre (rs6503875: 45,915,395 GRCh37 -> 47,838,029 GRCh38)
CENTER38 = 47838029
genes = eg[(eg.gene_chr.astype(str).str.replace("chr", "") == "17") &
           (eg.gene_end >= CENTER38 - HALF) & (eg.gene_start <= CENTER38 + HALF)].copy()
genes["dist_to_center_kb"] = ((genes.gene_start + genes.gene_end) / 2 - CENTER38) / 1000
genes.to_csv(f"{OUT}/tables/01_gene_list.csv", index=False, encoding="utf-8-sig")
print(f"[GENES] {len(genes)} genes in the window: " + ", ".join(genes.gene_name.tolist()))

# ---------------------------------------------------------------- 4. eQTLs in both tissues
th = pd.read_csv(f"{RAW}/thyroid_region_g38.tsv", sep="\t")
th = th.rename(columns={"gene_id": "gene_id", "pvalue": "p_eqtl", "beta": "beta_eqtl",
                        "se": "se_eqtl", "maf": "maf_eqtl", "rsid": "snp",
                        "position": "pos38", "ref": "ref", "alt": "alt"})
th = th[th.snp.notna() & (th.snp != ".")]
th["eqtl_ea"] = th.alt           # eQTL Catalogue: beta is relative to ALT
th["eqtl_oa"] = th.ref
th["tissue"] = "Thyroid"
th["n_eqtl"] = 574

wb = pd.read_csv(f"{RAW}/eqtlgen_region_all.tsv", sep="\t")
wb = wb.rename(columns={"SNP": "snp", "SNPPos": "pos37", "AssessedAllele": "eqtl_ea",
                        "OtherAllele": "eqtl_oa", "GeneSymbol": "gene_name"})
wb["gene_id"] = wb["Gene"]
wb["p_eqtl"] = wb.Pvalue.astype(float)
wb["z_eqtl"] = wb.Zscore.astype(float)
wb["n_eqtl"] = wb.NrSamples
wb["tissue"] = "WholeBlood"
# eQTLGen releases Z only: se = 1/sqrt(2p(1-p)(N+z^2)), with p from 1000G EUR
wb = wb.merge(ld_map, on="pos37", how="inner")
wb["se_eqtl"] = 1.0 / np.sqrt(2 * wb.maf1000g * (1 - wb.maf1000g) * (wb.n_eqtl + wb.z_eqtl ** 2))
wb["beta_eqtl"] = wb.z_eqtl * wb.se_eqtl
wb["maf_eqtl"] = wb.maf1000g

# add pos37 to the thyroid table (rsID match against the GWAS tables)
th = th.merge(rs2pos, on="snp", how="inner")
print(f"[EQTL] thyroid {len(th)} rows / {th.gene_id.nunique()} genes; "
      f"whole blood {len(wb)} rows / {wb.gene_id.nunique()} genes")

# ---------------------------------------------------------------- 5. per gene x tissue x outcome
def harmonise(e, g):
    m = e.merge(g.drop(columns=["snp"]), on="pos37", how="inner")
    ok = ((m.eqtl_ea == m.gwas_ea) & (m.eqtl_oa == m.gwas_oa)) | \
         ((m.eqtl_ea == m.gwas_oa) & (m.eqtl_oa == m.gwas_ea))
    m = m[ok].copy()
    flip = m.eqtl_ea == m.gwas_oa
    m["beta_out"] = np.where(flip, -m.beta, m.beta)
    m["z_out"] = m.beta_out / m.se_out
    m = m.drop_duplicates("pos37")
    m["snp"] = "chr17:" + m.pos37.astype(str)
    return m

lead_rows, panel_rows = [], []
for tissue, e in [("Thyroid", th), ("WholeBlood", wb)]:
    for gid, sub in e.groupby("gene_id"):
        sname = genes.loc[genes.gene_id == gid, "gene_name"]
        sname = sname.iloc[0] if len(sname) else sub.gene_name.iloc[0] if "gene_name" in sub else gid
        if "gene_name" not in sub:
            sub = sub.assign(gene_name=sname)
        for oc in OCS:
            m = harmonise(sub, gwas[oc])
            if len(m) < 30:
                continue
            m["maf"] = m.maf1000g if tissue == "WholeBlood" else np.minimum(m.maf_eqtl, 1 - m.maf_eqtl)
            m = m[(m.maf > 0.01) & m.beta_eqtl.notna() & (m.se_eqtl > 0)]
            if len(m) < 30:
                continue
            tag = f"{tissue}__{sname}__{oc}"
            m[["snp", "pos37", "eqtl_ea", "eqtl_oa", "gwas_ea", "beta_eqtl", "se_eqtl", "p_eqtl",
               "beta_out", "se_out", "p_out", "maf", "n_eqtl", "n_out"]].to_csv(
                f"{OUT}/inputs/{tag}.tsv", sep="\t", index=False)
            i = m.p_eqtl.idxmin()
            panel_rows.append(dict(tissue=tissue, gene=sname, gene_id=gid, outcome=oc,
                                   nsnp=len(m), lead_snp=m.loc[i, "snp"], lead_pos37=int(m.loc[i, "pos37"]),
                                   lead_p_eqtl=m.loc[i, "p_eqtl"], lead_beta_eqtl=m.loc[i, "beta_eqtl"],
                                   lead_se_eqtl=m.loc[i, "se_eqtl"], lead_maf=m.loc[i, "maf"],
                                   lead_beta_out=m.loc[i, "beta_out"], lead_se_out=m.loc[i, "se_out"],
                                   lead_p_out=m.loc[i, "p_out"]))
print(f"[INPUT] wrote {len(panel_rows)} gene x tissue x outcome combinations")

pn = pd.DataFrame(panel_rows)
pn.to_csv(f"{OUT}/tables/02_panel_lead.csv", index=False, encoding="utf-8-sig")

# overview: the best eQTL per tissue and gene
summary = (pn.sort_values("lead_p_eqtl").groupby(["tissue", "gene"]).first().reset_index()
             [["tissue", "gene", "lead_snp", "lead_p_eqtl", "lead_beta_eqtl", "lead_se_eqtl"]])
summary.to_csv(f"{OUT}/tables/02b_top_eqtl_per_gene.csv", index=False, encoding="utf-8-sig")
print(summary.to_string(index=False))
