#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p04_joint.py
============
Required analysis 2: marginal beta/SE/P of rs7219420 (per cohort) + cis-MR for
TBKBP1 + joint modelling of the two SNPs.

Joint model: with marginal estimates b=(b1,b2), diagonal SE matrix S and the
two-SNP LD correlation matrix R, the marginal covariance is V_m = S R S; the
joint effects are b_j = R^-1 b and V_j = R^-1 V_m R^-T.
When r2 is close to 0, b_j is close to b; the point is to quantify exactly that:
how far the estimates move after conditioning and how much the SE inflates.
"""
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1"
OUT = f"{BASE}/DATA21_complete_2026-11-14"
POS_A = 45915395     # rs6503875 GRCh37
POS_B = 45800299     # rs7219420 GRCh37

ld = pd.read_csv(f"{OUT}/inputs_ld/ld_TSH.tsv.gz", sep="\t", index_col=0)
ld.index = ld.index.astype(int); ld.columns = ld.columns.astype(int)
r = float(ld.loc[POS_A, POS_B])
R = np.array([[1.0, r], [r, 1.0]])
print(f"rs6503875-rs7219420: r={r:.4f}, r2={r**2:.5f}")

# ---------------------------------------------------------------- marginal estimates
disc = pd.read_csv(f"{BASE}/DATA15_thyroid_eQTL/results/gwas_win/TSH.tsv.gz", sep="\t")
row = lambda d, rs: d[d.rs_id == rs].iloc[0]
a_d = row(disc, "rs6503875"); b_d = row(disc, "rs7219420")

rep = pd.read_csv(f"{BASE}/DATA18_manuscript/analysis/second_cohort_region.tsv", sep="\t")
a_r = rep[rep.rsid == "rs6503875"].iloc[0]; b_r = rep[rep.rsid == "rs7219420"].iloc[0]


def align(beta, ea, ref_ea):
    return beta if ea == ref_ea else -beta


# harmonise to the C allele of rs6503875 and the C allele of rs7219420
marg = pd.DataFrame([
    dict(cohort="Discovery EUR (n=482,873)", snp="rs6503875",
         ea="C", beta=align(a_d.beta, a_d.effect_allele, "C"), se=a_d.standard_error,
         p=a_d.p_value, eaf=(a_d.effect_allele_frequency if a_d.effect_allele == "C"
                             else 1 - a_d.effect_allele_frequency), n=int(a_d.n)),
    dict(cohort="Discovery EUR (n=482,873)", snp="rs7219420",
         ea="C", beta=align(b_d.beta, b_d.effect_allele, "C"), se=b_d.standard_error,
         p=b_d.p_value, eaf=(b_d.effect_allele_frequency if b_d.effect_allele == "C"
                             else 1 - b_d.effect_allele_frequency), n=int(b_d.n)),
    dict(cohort="Replication EUR GCST90319321 (n=164,818)", snp="rs6503875",
         ea="C", beta=align(a_r.beta, a_r.effect_allele, "C"), se=a_r.standard_error,
         p=a_r.p_value, eaf=np.nan, n=164818),
    dict(cohort="Replication EUR GCST90319321 (n=164,818)", snp="rs7219420",
         ea="C", beta=align(b_r.beta, b_r.effect_allele, "C"), se=b_r.standard_error,
         p=b_r.p_value, eaf=np.nan, n=164818),
    dict(cohort="Replication SAS GCST90448634 (n=38,000)", snp="rs6503875",
         ea="C", beta=np.nan, se=np.nan, p=np.nan, eaf=np.nan, n=38000),
    dict(cohort="Replication SAS GCST90448634 (n=38,000)", snp="rs7219420",
         ea="C", beta=np.nan, se=np.nan, p=np.nan, eaf=np.nan, n=38000),
])
marg["z"] = marg.beta / marg.se
marg.to_csv(f"{OUT}/tables/05_rs7219420_marginal.csv", index=False, encoding="utf-8-sig")
print(marg.to_string(index=False))

# ---------------------------------------------------------------- joint / conditional
rows = []
for coh, b1, s1, b2, s2 in [
    ("Discovery EUR", marg.query("snp=='rs6503875' and cohort.str.startswith('Discovery')").beta.iloc[0],
     marg.query("snp=='rs6503875' and cohort.str.startswith('Discovery')").se.iloc[0],
     marg.query("snp=='rs7219420' and cohort.str.startswith('Discovery')").beta.iloc[0],
     marg.query("snp=='rs7219420' and cohort.str.startswith('Discovery')").se.iloc[0]),
    ("Replication EUR", marg.query("snp=='rs6503875' and cohort.str.startswith('Replication EUR')").beta.iloc[0],
     marg.query("snp=='rs6503875' and cohort.str.startswith('Replication EUR')").se.iloc[0],
     marg.query("snp=='rs7219420' and cohort.str.startswith('Replication EUR')").beta.iloc[0],
     marg.query("snp=='rs7219420' and cohort.str.startswith('Replication EUR')").se.iloc[0]),
]:
    b = np.array([b1, b2]); S = np.diag([s1, s2])
    Vm = S @ R @ S
    A = np.linalg.inv(R)
    bj = A @ b
    Vj = A @ Vm @ A.T
    sej = np.sqrt(np.diag(Vj))
    chi2 = float(bj @ np.linalg.inv(Vj) @ bj)
    pj = stats.chi2.sf(chi2, 2)
    z1 = b1 / s1; z2 = b2 / s2
    rows.append(dict(cohort=coh, r=r, r2=r ** 2,
                     b_marg_rs6503875=b1, se_marg_rs6503875=s1, p_marg_rs6503875=2 * stats.norm.sf(abs(z1)),
                     b_marg_rs7219420=b2, se_marg_rs7219420=s2, p_marg_rs7219420=2 * stats.norm.sf(abs(z2)),
                     b_joint_rs6503875=bj[0], se_joint_rs6503875=sej[0],
                     b_joint_rs7219420=bj[1], se_joint_rs7219420=sej[1],
                     shift_pct_rs6503875=100 * (bj[0] - b1) / b1,
                     shift_pct_rs7219420=100 * (bj[1] - b2) / b2,
                     se_inflation=1 / np.sqrt(1 - r ** 2),
                     chi2_joint=chi2, df=2, p_joint=pj))
jt = pd.DataFrame(rows)
jt.to_csv(f"{OUT}/tables/06_two_snp_joint_model.csv", index=False, encoding="utf-8-sig")
print()
print(jt.T.to_string())

# ---------------------------------------------------------------- cis-eQTL values at these SNPs in the locus
eq = []
for tissue in ["Thyroid", "WholeBlood"]:
    f = f"{OUT}/inputs/{tissue}__TBKBP1__TSH.tsv"
    d = pd.read_csv(f, sep="\t")
    for gene in ["TBKBP1", "SCRN2", "MRPL10"]:
        g = pd.read_csv(f"{OUT}/inputs/{tissue}__{gene}__TSH.tsv", sep="\t")
        for pos, rs in [(POS_B, "rs7219420"), (POS_A, "rs6503875")]:
            x = g[g.pos37 == pos]
            if len(x) == 0:
                continue
            x = x.iloc[0]
            eq.append(dict(tissue=tissue, gene=gene, snp=rs, pos37=pos,
                           eqtl_ea=x.eqtl_ea, beta_eqtl=x.beta_eqtl, se_eqtl=x.se_eqtl,
                           p_eqtl=x.p_eqtl, F=(x.beta_eqtl / x.se_eqtl) ** 2,
                           beta_tsh=x.beta_out, se_tsh=x.se_out,
                           wald=(x.beta_out / x.beta_eqtl if x.beta_eqtl else np.nan),
                           se_wald=(abs(x.beta_out / x.beta_eqtl) *
                                    np.sqrt((x.se_out / x.beta_out) ** 2 +
                                            (x.se_eqtl / x.beta_eqtl) ** 2)
                                    if x.beta_eqtl and x.beta_out else np.nan)))
eqt = pd.DataFrame(eq)
eqt["z_wald"] = eqt.wald / eqt.se_wald
eqt["p_wald"] = 2 * stats.norm.sf(abs(eqt.z_wald))
eqt.to_csv(f"{OUT}/tables/07_cis_eqtl_at_two_snps.csv", index=False, encoding="utf-8-sig")
print()
print(eqt.round(5).to_string(index=False))
