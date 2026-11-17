#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p05_regional.py
===============
Required analyses 3 and 5:
  1) regional plot data: TSH beta/SE/EAF/N plus the eQTLs of SCRN2/TBKBP1/MRPL10
     (thyroid, whole blood) over the same window, on a common LD reference
     (1000G EUR), together with r2 against the lead SNP.
  2) variant information table: rsID, GRCh37/38, effect/reference alleles, both
     EAF and N sets.
  3) plotting: a three-panel regional plot (TSH -log10P coloured by LD; the
     eQTLs of two genes).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = "D:/NewBeginning/MR1"
OUT = f"{BASE}/DATA21_complete_2026-11-14"
LEAD37 = 45915395
LEAD38 = 47838029

tsh = pd.read_csv(f"{BASE}/DATA15_thyroid_eQTL/results/gwas_win/TSH.tsv.gz", sep="\t")
tsh = tsh[(tsh.chromosome.astype(str) == "17") &
          (tsh.base_pair_location.between(LEAD37 - 500000, LEAD37 + 500000))].copy()
tsh = tsh.rename(columns={"base_pair_location": "pos37", "rs_id": "rsid"})
tsh = tsh.drop_duplicates("pos37")

ld = pd.read_csv(f"{OUT}/inputs_ld/ld_TSH.tsv.gz", sep="\t", index_col=0)
ld.index = ld.index.astype(int); ld.columns = ld.columns.astype(int)
r2 = (ld[LEAD37] ** 2).rename("r2_lead")

reg = tsh[["pos37", "rsid", "effect_allele", "other_allele", "beta", "standard_error",
           "effect_allele_frequency", "p_value", "n"]].copy()
reg = reg.rename(columns={"beta": "beta_tsh", "standard_error": "se_tsh",
                          "effect_allele_frequency": "eaf_tsh", "p_value": "p_tsh",
                          "effect_allele": "tsh_ea", "other_allele": "tsh_oa", "n": "n_tsh"})
reg = reg.merge(r2, left_on="pos37", right_index=True, how="left")

# eQTLs for three genes in two tissues
for tissue, tag in [("Thyroid", "thyroid"), ("WholeBlood", "whole blood")]:
    for gene in ["SCRN2", "TBKBP1", "MRPL10"]:
        f = f"{OUT}/inputs/{tissue}__{gene}__TSH.tsv"
        g = pd.read_csv(f, sep="\t")[["pos37", "eqtl_ea", "beta_eqtl", "se_eqtl", "p_eqtl"]]
        g = g.rename(columns={"eqtl_ea": f"{gene}_ea_{tissue}",
                              "beta_eqtl": f"{gene}_beta_{tissue}",
                              "se_eqtl": f"{gene}_se_{tissue}",
                              "p_eqtl": f"{gene}_p_{tissue}"})
        reg = reg.merge(g, on="pos37", how="left")

reg.to_csv(f"{OUT}/tables/08_regional_plot_data.csv", index=False, encoding="utf-8-sig")
print("regional data:", reg.shape)
print(reg.dropna(subset=["p_tsh"]).p_tsh.min())

# ---------------------------------------------------------------- variant information table
KEY = {"rs6503875": 45915395, "rs7219420": 45800299, "rs62076102": 45887132,
       "rs4794224": 45957463, "rs11079797": None, "rs1989291": None,
       "rs62076103": None, "rs3047616": None}
POS38 = {"rs6503875": 47838029, "rs7219420": 47722933, "rs62076102": 47809766,
         "rs4794224": 47880097}
rows = []
for rs, p37 in KEY.items():
    if p37 is None:
        x = reg[reg.rsid == rs]
        if len(x) == 0:
            continue
        p37 = int(x.pos37.iloc[0])
    x = reg[reg.pos37 == p37]
    if len(x) == 0:
        continue
    x = x.iloc[0]
    rows.append(dict(rsid=rs, pos_grch37=int(p37), pos_grch38=POS38.get(rs, np.nan),
                     tsh_ea=x.tsh_ea, tsh_oa=x.tsh_oa,
                     eaf_tsh_of_tsh_ea=round(float(x.eaf_tsh), 4) if pd.notna(x.eaf_tsh) else np.nan,
                     n_tsh=int(x.n_tsh) if pd.notna(x.n_tsh) else np.nan,
                     beta_tsh=x.beta_tsh, se_tsh=x.se_tsh, p_tsh=x.p_tsh,
                     r2_with_rs6503875=x.r2_lead))
vt = pd.DataFrame(rows)
vt.to_csv(f"{OUT}/tables/11_variant_info.csv", index=False, encoding="utf-8-sig")
print(vt.to_string(index=False))

# ---------------------------------------------------------------- figure
def manh(ax, sub, ycol, title, color_by_ld=True, ylab="-log10 P"):
    sub = sub.dropna(subset=[ycol])
    y = -np.log10(np.clip(sub[ycol].values, 1e-320, None))
    c = sub.r2_lead.fillna(0).values
    order = np.argsort(c)
    if color_by_ld:
        sc = ax.scatter(sub.pos37.values[order] / 1e6, y[order], c=c[order],
                        cmap="viridis_r", s=9, linewidths=0, vmin=0, vmax=1)
    else:
        ax.scatter(sub.pos37.values / 1e6, y, s=9, color="#c0392b", linewidths=0)
    ax.set_title(title, fontsize=10, loc="left")
    ax.set_ylabel(ylab, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    return sc if color_by_ld else None


fig, axes = plt.subplots(3, 1, figsize=(9, 8.6), sharex=True,
                         gridspec_kw=dict(height_ratios=[1.5, 1, 1], hspace=0.12))
sc = manh(axes[0], reg, "p_tsh", "A  TSH GWAS (n = 482,873)")
for ax in axes:
    for pos, lab, cc in [(45.800299, "rs7219420", "#2c3e50"), (45.915395, "rs6503875", "#8e44ad")]:
        ax.axvline(pos, ls=":", lw=1, color=cc, alpha=0.9)
axes[0].annotate("rs7219420", (45.800299, 6.9), xytext=(3, 4), textcoords="offset points",
                 fontsize=8, color="#2c3e50")
axes[0].annotate("rs6503875", (45.915395, 4.7), xytext=(3, 4), textcoords="offset points",
                 fontsize=8, color="#8e44ad")
cb = fig.colorbar(sc, ax=axes[0], pad=0.01, fraction=0.03)
cb.set_label("r² with rs6503875", fontsize=8); cb.ax.tick_params(labelsize=7)

for ax, gene, tissue, ttl in [(axes[1], "TBKBP1", "Thyroid", "B  TBKBP1 eQTL, thyroid"),
                              (axes[2], "SCRN2", "Thyroid", "C  SCRN2 eQTL, thyroid")]:
    col = f"{gene}_p_{tissue}"
    manh(ax, reg, col, ttl, color_by_ld=False)
axes[2].set_xlabel("Position on chr17 (Mb, GRCh37)", fontsize=9)
fig.suptitle("17q21.32: TSH and cis-eQTL for TBKBP1 / SCRN2", fontsize=11, y=0.94)
fig.savefig(f"{OUT}/fig_regional.png", dpi=200, bbox_inches="tight")
print("figure saved")
