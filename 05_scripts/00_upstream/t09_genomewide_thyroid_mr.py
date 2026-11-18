#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
t09_genomewide_thyroid_mr.py
============================
Genome-wide MR of the EDC candidate genes with thyroid-tissue eQTLs
(replaces the earlier B06 step).

  1. take the lead thyroid cis-eQTL of every EDC candidate gene (p < 5e-8,
     instrument grade)
  2. pull those SNPs out of the four full thyroid-outcome GWAS files
     (TSH/T4/HYPO by rsID; HASH carries no rsID, so by GRCh37 position)
  3. Wald ratio, then FDR within each outcome

Output
  results/gw_thyroid_ivs.csv
  results/gw_thyroid_mr.csv
  results/gw_thyroid_mr_summary.txt
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1/DATA15_thyroid_eQTL"
D10 = "D:/NewBeginning/MR1/DATA10"
R = os.path.join(BASE, "results")
TMP = os.path.join(BASE, "tmp")

OC = {
    "TSH":  ("TSH.tsv.gz", 482873),
    "T4":   ("T4.tsv.gz", 191449),
    "HYPO": ("HYPO.tsv.gz", 1046991),
    "HASH": ("HASH_full.tsv.gz", 395640),
}
KEEP = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
        "beta", "standard_error", "p_value", "rs_id"]


def build_ivs():
    d = pd.read_csv(os.path.join(TMP, "thyroid_gw_sig.tsv"), sep="\t")
    for c in ["position", "maf", "pvalue", "beta", "se", "an"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["rsid", "pvalue", "beta", "se"])
    d = d[d.se > 0]
    d = d.sort_values("pvalue").drop_duplicates(["molecular_trait_id", "position"])
    lead = d.sort_values("pvalue").groupby("molecular_trait_id", as_index=False).first()
    lead = lead[lead.pvalue < 5e-8].copy()
    lead["F"] = (lead.beta / lead.se) ** 2
    m = pd.read_csv(os.path.join(TMP, "cand_ensg.tsv"), sep="\t")
    lead = lead.merge(m, left_on="molecular_trait_id", right_on="ensg", how="left")
    out = pd.DataFrame({
        "gene": lead.gene_name, "ensg": lead.molecular_trait_id,
        "rsid": lead.rsid, "chr": lead.chromosome,
        "pos38": lead.position.astype("Int64"),
        "ref": lead.ref, "alt": lead.alt,
        "b_e": lead.beta, "se_e": lead.se, "p_e": lead.pvalue,
        "maf": lead.maf, "F": lead.F,
    })
    out.to_csv(os.path.join(R, "gw_thyroid_ivs.csv"), index=False, encoding="utf-8-sig")
    return out


def scan(oc, rs, want_pos=None):
    fn = OC[oc][0]
    path = os.path.join(D10, "raw", fn)
    parts, tot = [], 0
    for ch in pd.read_csv(path, sep="\t", chunksize=2_000_000,
                          dtype={"chromosome": str}, low_memory=False,
                          na_values=["NA", ""]):
        tot += len(ch)
        cols = set(ch.columns)
        if "rs_id" in cols:
            sel = ch[ch["rs_id"].isin(rs)]
        elif want_pos is not None:
            p = pd.to_numeric(ch["base_pair_location"], errors="coerce")
            sel = ch[p.isin(want_pos)]
        else:
            sel = ch.iloc[0:0]
        if len(sel):
            parts.append(sel[[c for c in KEEP if c in sel.columns]].copy())
    d = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=KEEP)
    print(f"  {oc}: scanned {tot:,} rows -> {len(d):,} matched", flush=True)
    return d


def main():
    print("building thyroid instruments ...", flush=True)
    lead = build_ivs()
    print(f"  candidate genes with a p<5e-8 thyroid cis-eQTL: {len(lead):,}", flush=True)

    rs = set(lead.rsid)
    gws, posmap = {}, {}
    for oc in ["TSH", "T4", "HYPO"]:
        d = scan(oc, rs)
        if "rs_id" not in d.columns:
            continue
        d = d.drop_duplicates("rs_id")
        posmap.update(dict(zip(d["rs_id"],
                               pd.to_numeric(d["base_pair_location"], errors="coerce"))))
        gws[oc] = d
    lead["pos37"] = lead.rsid.map(posmap)
    print(f"  rsID -> GRCh37 mappings: {len(posmap):,}", flush=True)

    want = set(lead.pos37.dropna().astype(int))
    gws["HASH"] = scan("HASH", set(), want)

    rows = []
    for oc in OC:
        d = gws.get(oc)
        if d is None or d.empty:
            continue
        d = d.copy()
        d["bp"] = pd.to_numeric(d["base_pair_location"], errors="coerce")
        for c in ["beta", "standard_error"]:
            d[c] = pd.to_numeric(d[c], errors="coerce")
        if oc == "HASH":
            mm = lead.dropna(subset=["pos37"]).merge(
                d, left_on="pos37", right_on="bp", how="inner")
        else:
            mm = lead.merge(d, left_on="rsid", right_on="rs_id", how="inner")
        same = (mm["alt"] == mm["effect_allele"]) & (mm["ref"] == mm["other_allele"])
        flip = (mm["alt"] == mm["other_allele"]) & (mm["ref"] == mm["effect_allele"])
        mm = mm[same | flip].copy()
        sgn = np.where(mm["alt"] == mm["effect_allele"], 1.0, -1.0)
        by = mm["beta"].values * sgn
        sy = mm["standard_error"].values
        bx, sx = mm["b_e"].values, mm["se_e"].values
        b = by / bx
        se = np.sqrt(sy ** 2 / bx ** 2 + by ** 2 * sx ** 2 / bx ** 4)
        z = b / se
        rows.append(pd.DataFrame(dict(
            gene=mm["gene"].values, ensg=mm["ensg"].values, outcome=oc,
            snp=mm["rsid"].values, F=mm["F"].values,
            b=b, se=se, z=z, p=2 * stats.norm.sf(np.abs(z)))))

    out = pd.concat(rows, ignore_index=True)
    out["fdr"] = np.nan
    for oc in OC:
        m = out.outcome == oc
        if not m.any():
            continue
        pv = out.loc[m, "p"].values
        ok = ~np.isnan(pv)
        q = np.full(len(pv), np.nan)
        if ok.any():
            q[ok] = stats.false_discovery_control(pv[ok])
        out.loc[m, "fdr"] = q
    out = out.sort_values(["outcome", "p"]).reset_index(drop=True)
    out.to_csv(os.path.join(R, "gw_thyroid_mr.csv"), index=False, encoding="utf-8-sig")

    L = ["Genome-wide MR with thyroid eQTLs (EDC candidate genes)", "=" * 72,
         f"Candidate genes with an instrument-grade (p<5e-8) thyroid cis-eQTL: {len(lead):,}",
         f"Gene-outcome pairs with an MR estimate: {len(out):,}", ""]
    for oc in OC:
        s = out[out.outcome == oc]
        if s.empty:
            continue
        L.append(f"{oc:<5} n={len(s):>5,}  nominal p<0.05: {int((s.p<0.05).sum()):>4}"
                 f"  FDR<0.05: {int((s.fdr<0.05).sum()):>4}"
                 f"  FDR<0.10: {int((s.fdr<0.10).sum()):>4}")
    L.append("")
    for oc in OC:
        s = out[(out.outcome == oc) & (out.fdr < 0.05)].sort_values("p")
        if s.empty:
            continue
        L.append(f"{oc} FDR<0.05 ({len(s)}): " + ", ".join(s.gene.head(30).astype(str)))
    txt = "\n".join(L)
    open(os.path.join(R, "gw_thyroid_mr_summary.txt"), "w", encoding="utf-8").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
