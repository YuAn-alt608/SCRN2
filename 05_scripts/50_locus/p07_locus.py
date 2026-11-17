#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p07_locus.py
============
Required analysis 4 (locus identity) + required analysis 6 (sample overlap) +
completion of required analysis 5 (GRCh38 coordinates).

Locus identity covers three things:
  (a) the minimum P of the four outcomes inside the local locus, and whether it
      passes 5e-8;
  (b) associations already reported in GWAS Catalog (v2 REST, cached locally in
      tables/12);
  (c) an explicit note that the sentinel tables of Rand 2025 / Williams 2023 are
      not available locally, so their independent-locus counts cannot be checked.

Sample overlap: reported as the three-cohort estimate and heterogeneity Q;
splitting by cohort and dropping UK Biobank would require the Rand supplementary
tables, which are not available locally.
"""
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1"
OUT = f"{BASE}/DATA21_complete_2026-11-14"
CENTER = 45915395

# ---------------------------------------------------------------- (a) regional minimum P
rows = []
for oc in ["TSH", "T4", "HYPO", "HASH"]:
    d = pd.read_csv(f"{BASE}/DATA15_thyroid_eQTL/results/gwas_win/{oc}.tsv.gz", sep="\t")
    d = d[(d.chromosome.astype(str) == "17") &
          (d.base_pair_location.between(CENTER - 500000, CENTER + 500000))]
    d = d.assign(p=pd.to_numeric(d.p_value, errors="coerce"))
    i = d.p.idxmin()
    rs = d.loc[i, "rs_id"] if "rs_id" in d and pd.notna(d.loc[i, "rs_id"]) else "-"
    rows.append(dict(outcome=oc, n_snp=len(d), min_p=float(d.loc[i, "p"]),
                     min_p_snp=rs, pos_grch37=int(d.loc[i, "base_pair_location"]),
                     passes_5e8=bool(d.loc[i, "p"] < 5e-8),
                     n_snps_below_5e8=int((d.p < 5e-8).sum())))
lt = pd.DataFrame(rows)
lt.to_csv(f"{OUT}/tables/13_locus_identity_regional_min_p.csv", index=False, encoding="utf-8-sig")
print(lt.to_string(index=False))

# ---------------------------------------------------------------- (b) previously reported associations
try:
    gc = pd.read_csv(f"{OUT}/tables/12_gwas_catalog_locus_identity.csv")
    print("\nGWAS Catalog hits (4 variants in this window):")
    print(gc[["rsid", "p", "efo", "study", "genes", "author"]].to_string(index=False))
except Exception as e:
    print("GWAS Catalog table missing", e)

# ---------------------------------------------------------------- (c) three cohorts + heterogeneity
tc = pd.DataFrame([
    dict(cohort="Rand2025_EUR (discovery)", n=482873, est=0.03359, se=0.00885, p=1.48e-4),
    dict(cohort="GCST90319321_EUR (= Figueredo)", n=164818, est=0.04038, se=0.01444, p=5.16e-3),
    dict(cohort="GCST90448634_SAS (Jacobs)", n=38000, est=0.02882, se=0.02413, p=0.232),
])
w = 1 / tc.se ** 2
b_pool = float((w * tc.est).sum() / w.sum())
se_pool = float(np.sqrt(1 / w.sum()))
Q = float((w * (tc.est - b_pool) ** 2).sum())
pq = float(stats.chi2.sf(Q, len(tc) - 1))
tc["weight"] = w / w.sum()
meta = pd.DataFrame([dict(cohort="fixed-effect pooled", n=int(tc.n.sum()), est=b_pool,
                          se=se_pool, p=2 * stats.norm.sf(abs(b_pool / se_pool)),
                          Q=Q, df=len(tc) - 1, p_Q=pq, I2=max(0, (Q - (len(tc) - 1)) / Q) * 100)])
tc.to_csv(f"{OUT}/tables/14_three_cohort_sample_overlap.csv", index=False, encoding="utf-8-sig")
meta.to_csv(f"{OUT}/tables/14b_three_cohort_pooled.csv", index=False, encoding="utf-8-sig")
print("\n", tc.to_string(index=False))
print(meta.to_string(index=False))

# ---------------------------------------------------------------- fill in GRCh38 coordinates
vt = pd.read_csv(f"{OUT}/tables/11_variant_info.csv")
th = pd.read_csv(f"{OUT}/raw/thyroid_region_g38.tsv", sep="\t", usecols=["rsid", "position"])
th = th.dropna().drop_duplicates("rsid").rename(columns={"position": "pos_grch38_new"})
vt = vt.drop(columns=["pos_grch38"]).merge(th, left_on="rsid", right_on="rsid", how="left")
vt = vt.rename(columns={"pos_grch38_new": "pos_grch38"})
vt.to_csv(f"{OUT}/tables/11_variant_info.csv", index=False, encoding="utf-8-sig")
print("\n", vt.to_string(index=False))
