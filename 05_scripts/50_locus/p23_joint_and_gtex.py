#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p23_joint_and_gtex.py -- two tasks:
  (1) is the shoulder variant (chr17:45,873,049) an eQTL for some gene in the 49
      GTEx v8 tissues?
  (2) use a GLS joint model (b_j = (R'V^-1 R)^-1 R'V^-1 b, V = diag(se).R.diag(se))
      to decide which of the five positions 45.800 / 45.807 / 45.873 / 45.887 /
      45.915 stay independent after conditioning

Why GLS rather than R^-1 z on the z scale:
  the latter implicitly assumes equal SE across SNPs. In this interval the SE
  ranges from 0.0023 to 0.0093, a four-fold spread, and rs62076102 with
  rs6503875 is near-collinear at r=0.767; R^-1 z would also absorb the pair at
  r~0.03 (rs6503875 z collapses from -4.22 to -1.15 in practice), which
  contradicts the LD structure. GLS is the correct form.
"""
import glob
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
T = f"{BASE}/tables"
SHOULDER = 45_873_049
SHOULDER_B38 = 47_795_683          # GRCh37 + 1,922,634

# ---------------- (1) GTEx across tissues ----------------
GTEX = "D:/NewBeginning/MR1/DATA10/tmp/gtex/GTEx_Analysis_v8_eQTL"
rows = []
if os.path.exists(f"{T}/30_shoulder_gtex_cross_tissue.csv"):
    print("(reusing the existing 30_shoulder_gtex_cross_tissue.csv)")
for f in [] if os.path.exists(f"{T}/30_shoulder_gtex_cross_tissue.csv") else sorted(glob.glob(f"{GTEX}/*.v8.signif_variant_gene_pairs.txt.gz")):
    tis = os.path.basename(f).split(".v8.signif")[0]
    d = pd.read_csv(f, sep="\t", usecols=["variant_id", "gene_id", "pval_nominal", "slope"])
    v = d.variant_id.str.split("_", expand=True)
    d = d[(v[0] == "chr17") & (v[1].astype(int) == SHOULDER_B38)]
    if not len(d):
        continue
    for _, r in d.iterrows():
        rows.append(dict(tissue=tis, gene=r.gene_id, P=r.pval_nominal, slope=r.slope))
if not rows:
    g = pd.read_csv(f"{T}/30_shoulder_gtex_cross_tissue.csv")
    print(f"=== shoulder variant as a GTEx eQTL: {len(g)} records, {g.tissue.nunique()} tissues, {g.gene.nunique()} genes ===")
    print(g.groupby("gene_name").agg(n_tissues=("tissue", "nunique"), min_P=("P", "min")).to_string())
else:
    g = pd.DataFrame(rows)
if rows:
    # add gene names (the 23 inside the window plus the 2 outside that still hit here)
    ens = pd.read_csv(f"{BASE}/raw/ens_genes_17_peak_b37.csv")[["gene_id", "gene_name"]]
    name_map = dict(zip(ens.gene_id, ens.gene_name))
    name_map.update({"ENSG00000178852": "EFCAB13", "ENSG00000228782": "MRPL45P2"})
    g["gene_name"] = g.gene.str.split(".").str[0].map(name_map)
    # the panel's thyroid eQTL is GTEx v8 thyroid (QTD000341, n=574) -- the same
    # source, so it cannot count as independent validation
    g["same_source_as_panel"] = g.tissue == "Thyroid"
    g = g.sort_values("P")
    g.to_csv(f"{T}/30_shoulder_gtex_cross_tissue.csv", index=False, encoding="utf-8-sig")
print(f"=== significant GTEx eQTLs at the shoulder variant: {len(g)} records, {g.tissue.nunique()} tissues, {g.gene.nunique()} genes ===")
if len(g):
    print(g.head(15).round(6).to_string(index=False))
print("\nsummary per gene (number of tissues / strongest P):")
print(g.groupby("gene_name").agg(n_tissues=("tissue", "nunique"), min_P=("P", "min"),
                                 best_tissue=("P", lambda s: g.loc[s.idxmin(), "tissue"])).to_string())

# ---------------- (2) joint conditional model ----------------
ld = pd.read_csv(f"{BASE}/inputs_ld/ld_TSH.tsv.gz", sep="\t", index_col=0)
ld.index = ld.index.astype(int)
ld.columns = ld.columns.astype(int)

BASE_SNP = {"rs7219420": 45_800_299, "rs1989291": 45_807_547, "shoulder": SHOULDER,
            "rs62076102": 45_887_132, "rs6503875": 45_915_395}
MODELS = {
    "five-SNP": dict(BASE_SNP),
    "four-SNP (minus rs62076102)": {k: v for k, v in BASE_SNP.items() if k != "rs62076102"},
    "pair:rs7219420+rs1989291": {"rs7219420": 45_800_299, "rs1989291": 45_807_547},
    "pair:rs1989291+shoulder": {"rs1989291": 45_807_547, "shoulder": SHOULDER},
    "pair:shoulder+rs62076102": {"shoulder": SHOULDER, "rs62076102": 45_887_132},
    "pair:rs62076102+rs6503875": {"rs62076102": 45_887_132, "rs6503875": 45_915_395},
    "pair:rs7219420+rs6503875": {"rs7219420": 45_800_299, "rs6503875": 45_915_395},
    "pair:shoulder+rs6503875": {"shoulder": SHOULDER, "rs6503875": 45_915_395},
}


def gls_joint(d, snps):
    pos = list(snps.values())
    R = ld.loc[pos, pos].values
    b = np.array([d.loc[p, "beta_out"] for p in pos])
    se = np.array([d.loc[p, "se_out"] for p in pos])
    V = np.outer(se, se) * R
    Vi = np.linalg.inv(V)
    M = np.linalg.inv(R.T @ Vi @ R)
    bj = M @ (R.T @ Vi @ b)
    return bj, np.sqrt(np.diag(M)), R


out = []
for oc in ["TSH", "T4", "HYPO", "HASH"]:
    f = f"{BASE}/inputs/Thyroid__SCRN2__{oc}.tsv"
    d = pd.read_csv(f, sep="\t", usecols=["pos37", "beta_out", "se_out"])
    d = d.drop_duplicates("pos37").set_index("pos37")
    for mname, snps in MODELS.items():
        bj, sej, R = gls_joint(d, snps)
        for i, (n, p) in enumerate(snps.items()):
            z_m = d.loc[p, "beta_out"] / d.loc[p, "se_out"]
            out.append(dict(outcome=oc, model=mname, SNP=n, pos37=p,
                            marginal_b=d.loc[p, "beta_out"], marginal_se=d.loc[p, "se_out"],
                            marginal_z=z_m, marginal_P=2 * stats.norm.sf(abs(z_m)),
                            joint_b=bj[i], joint_se=sej[i], joint_z=bj[i] / sej[i],
                            joint_P=2 * stats.norm.sf(abs(bj[i] / sej[i])),
                            r_with_first_snp=R[0, i] if i else 1.0))
j = pd.DataFrame(out)
j.to_csv(f"{T}/31_joint_conditional_model.csv", index=False, encoding="utf-8-sig")

print("\nNote: GTEx Thyroid and the panel thyroid data come from the same source (QTD000341), so it is not independent validation; the other 22 tissues are independent.")
print("\n=== joint model: TSH ===")
q = j[(j.outcome == "TSH")]
print(q[["model", "SNP", "marginal_z", "joint_z", "joint_P"]].round(4).to_string(index=False))
print("\n=== joint model: hypothyroidism ===")
q = j[(j.outcome == "HYPO")]
print(q[["model", "SNP", "marginal_z", "joint_z", "joint_P"]].round(4).to_string(index=False))
