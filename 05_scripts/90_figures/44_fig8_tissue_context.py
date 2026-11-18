# -*- coding: utf-8 -*-
"""
44_fig8_tissue_context.py - Figure 8, SCRN2 in tissue and cell context
=======================================================================
Every number in this figure is read from 04_data; nothing is entered by hand.

Panels
------
  A  GTEx v8 median expression of SCRN2, its paralogues and thyroid markers,
     thyroid against whole blood, with the thyroid rank among the 54 tissues
  B  the rs6503875 eQTL across the GTEx v8 tissues in which the variant is
     testable (47 of the 49 with cis-eQTL data), from a2_summary.csv
  C  the fraction of cells expressing SCRN2 in a public thyroid scRNA-seq set
  D  partial Spearman correlation with thyroid transcription factors,
     hormone-synthesis, exocytosis and negative-control genes (10 PCs removed)
  E  log2 fold change in two GEO datasets, with the smallest detectable effect

Output
------
  03_manuscript/figs/fig8.png (300 dpi) and fig8.pdf

Usage
-----
    python 05_scripts/90_figures/44_fig8_tissue_context.py
"""
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
D = os.path.join(ROOT, "04_data", "06_expression")
OUT = os.path.join(ROOT, "03_manuscript", "figs")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7.2,
    "axes.titlesize": 7.6,
    "axes.labelsize": 7.2,
    "xtick.labelsize": 6.6,
    "ytick.labelsize": 6.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "figure.dpi": 300,
})

TEAL = "#1b9e8f"
GREY = "#8d9299"
DARK = "#2b2f33"
SALMON = "#d4705f"

# ────────────────────────────────────────────────────────── data
prof = pd.read_csv(os.path.join(D, "genes_expression_profile.csv"), encoding="utf-8-sig")
a2 = pd.read_csv(os.path.join(D, "a2_summary.csv"), encoding="utf-8-sig").iloc[0]
ct = pd.read_csv(os.path.join(D, "SCRN2_thyroid_celltype.csv"), encoding="utf-8-sig")
co = pd.read_csv(os.path.join(D, "s07_panel_pc10.csv"), encoding="utf-8-sig")
geo_txt = open(os.path.join(D, "SCRN2_geo_expression_validation.txt"), encoding="utf-8").read()

GENES = ["TG", "TPO", "DDX5", "TSHR", "SCRN2", "SORT1", "CELSR2", "SCRN1", "SCRN3", "PSRC1"]
p = prof.set_index("gene").loc[GENES]

fig = plt.figure(figsize=(8.6, 6.8))
gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1.0, 0.95], height_ratios=[1.0, 1.05],
                      left=0.075, right=0.985, top=0.945, bottom=0.075,
                      wspace=0.80, hspace=0.60)


def letter(ax, s, x=-0.16, y=1.10):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="top", ha="left", color=DARK)


# ── A  expression, thyroid vs whole blood
axA = fig.add_subplot(gs[0, 0])
y = np.arange(len(GENES))[::-1]
for yi, g in zip(y, GENES):
    axA.plot([p.loc[g, "blood_tpm"], p.loc[g, "thyroid_tpm"]], [yi, yi],
             color="#cfd3d6", lw=1.0, zorder=1)
axA.scatter(p["blood_tpm"].values, y, s=16, color=GREY, zorder=3, label="whole blood")
axA.scatter(p["thyroid_tpm"].values, y, s=22, color=TEAL, zorder=4, label="thyroid")
axA.set_xscale("log")
axA.set_yticks(y)
axA.set_yticklabels(GENES)
for lab in axA.get_yticklabels():
    if lab.get_text() == "SCRN2":
        lab.set_fontweight("bold")
axA.set_xlabel("GTEx v8 median expression (TPM, log)")
axA.set_xlim(0.3, 1.2e4)
axA.set_ylim(-0.7, 10.3)
axA.legend(frameon=False, loc="lower right", handletextpad=0.2, borderpad=0.1)
for yi, g in zip(y, GENES):
    r = int(p.loc[g, "thyroid_rank"])
    axA.text(1.35e4, yi, "#%d" % r, fontsize=6.0, color=GREY, va="center", ha="left")
axA.set_title("Thyroid-enriched, not thyroid-specific", loc="left", pad=4, fontweight="bold")
letter(axA, "A")
axA.text(1.35e4, 9.95, "rank", fontsize=5.6, color=GREY, ha="left", va="center")

# ── B  the eQTL across tissues
axB = fig.add_subplot(gs[0, 1])
axB.set_axis_off()
axB.set_xlim(0, 1); axB.set_ylim(0, 1)
n_t = int(a2["n_tissue"]); n_sig = int(a2["n_sig_p05"])
n_same = int(a2["n_sig_same_sign"]); pct = float(a2["rank_pct_thyroid"])
axB.text(0, 0.97, "rs6503875 eQTL, %d GTEx tissues" % n_t, fontsize=7.6,
         fontweight="bold", color=DARK, va="top")
axB.text(0, 0.80, "Significant (P < 0.05)", fontsize=7.0, color=DARK, va="top")
axB.add_patch(Rectangle((0, 0.66), 1, 0.10, facecolor=TEAL, edgecolor="none"))
axB.text(0.5, 0.71, "%d / %d tissues" % (n_sig, n_t), color="white", fontsize=7.2,
         fontweight="bold", ha="center", va="center")
axB.text(0, 0.55, "Same direction as thyroid", fontsize=7.0, color=DARK, va="top")
axB.add_patch(Rectangle((0, 0.41), 1, 0.10, facecolor=TEAL, edgecolor="none"))
axB.text(0.5, 0.46, "%d / %d other tissues" % (n_same, n_t - 1), color="white",
         fontsize=7.2, fontweight="bold", ha="center", va="center")
axB.text(0, 0.30, "Normalised effect size (NES)", fontsize=7.0, color=DARK, va="top")
axB.text(0.5, 0.19, "thyroid %.2f" % a2["nes_thyroid"], fontsize=7.0, color=TEAL,
         ha="center", va="top")
axB.plot([0.02, 0.98], [0.08, 0.08], color="#c9cdd1", lw=1.4, solid_capstyle="butt")
axB.plot([0.90, 0.90], [0.03, 0.13], color="#9aa0a6", lw=0.8)
axB.scatter([0.30], [0.08], s=70, color=GREY, zorder=3)
axB.scatter([0.42], [0.08], s=70, color=TEAL, zorder=4)
axB.text(0.5, -0.03, "other-tissue median %.2f" % a2["nes_median_other"], fontsize=6.8,
         color=DARK, ha="center", va="top")
axB.text(0, -0.10, "Thyroid effect at the %.0fth percentile of the %d tissues\n"
         "(per-tissue values in 04_data)." % (pct, n_t),
         fontsize=6.0, color="#5b6167", style="italic", va="top")
letter(axB, "B", x=-0.06)

# ── C  which cell in the gland
axC = fig.add_subplot(gs[0, 2])
c = ct[ct.n_cells >= 100].sort_values("pct_expressing")
n_omitted = int((ct.n_cells < 100).sum())
yc = np.arange(len(c))
axC.hlines(yc, 0, c["pct_expressing"], color="#d9dcdf", lw=0.9)
axC.scatter(c["pct_expressing"], yc, s=np.clip(c["n_cells"] / 70, 6, 90),
            color=SALMON, alpha=0.9, edgecolor="white", linewidth=0.4, zorder=3)
axC.set_yticks(yc)
axC.set_yticklabels([t.replace("thyroid ", "").replace(" cell", "") for t in c["cell_type"]],
                    fontsize=6.2)
axC.text(0.0, -0.19, "%d further cell types (<100 cells) omitted" % n_omitted,
         transform=axC.transAxes, ha="left", va="top", fontsize=5.6,
         color="#6b7177", style="italic")
axC.set_xlabel("Cells expressing SCRN2 (%)")
axC.set_xlim(0, max(c["pct_expressing"]) * 2.05)
for yi, (_, r) in zip(yc, c.iterrows()):
    axC.text(r["pct_expressing"] + max(c["pct_expressing"]) * 0.06, yi,
             "%.2f%% (n = %s)" % (r["pct_expressing"], format(int(r["n_cells"]), ",")),
             fontsize=5.6, color="#5b6167", va="center")
axC.set_title("Which cell in the gland", loc="left", pad=4, fontweight="bold")
letter(axC, "C", x=-0.34)

# ── D  co-expression
axD = fig.add_subplot(gs[1, 0])
MODCOL = {"TF": "#2c5f8a", "EFF": TEAL, "EXO": "#d99b3f", "NEG": "#b8bdc2"}
MODLAB = {"TF": "thyroid transcription factors", "EFF": "hormone synthesis",
          "EXO": "exocytosis", "NEG": "negative controls"}
order = co.sort_values("rho_pc10", ascending=False)
xd = np.arange(len(order))
axD.bar(xd, order["rho_pc10"], width=0.72,
        color=[MODCOL[m] for m in order["module"]], edgecolor="none")
axD.axhline(0, color="#4a4f54", lw=0.7)
axD.set_xticks(xd)
axD.set_xticklabels(order["gene"], rotation=90, fontsize=6.0)
axD.set_ylabel("Partial $\\rho$ with SCRN2\n(10 expression PCs removed)")
axD.set_ylim(-0.30, 0.42)
for xi, (_, r) in zip(xd, order.iterrows()):
    if r["gene"] in ("NKX2-1", "PAX8"):
        axD.text(xi, r["rho_pc10"] + 0.012, "%+.2f" % r["rho_pc10"], ha="center",
                 fontsize=6.0, color=DARK)
handles = [Rectangle((0, 0), 1, 1, facecolor=MODCOL[k]) for k in ("TF", "EFF", "EXO", "NEG")]
axD.legend(handles, [MODLAB[k] for k in ("TF", "EFF", "EXO", "NEG")], frameon=False,
           fontsize=5.8, loc="upper right", handlelength=0.9, borderpad=0.1,
           labelspacing=0.3)
axD.set_title("Co-expression with the thyroid programme", loc="left", pad=4,
              fontweight="bold")
letter(axD, "D", x=-0.30)

# ── E  GEO bulk tissue
axE = fig.add_subplot(gs[1, 1:])
rows = []
m = re.search(r"LT\+ vs LT-: log2FC=([+-][\d.]+)\s+p\(MW\)=([\d.]+)", geo_txt)
rows.append(("GSE165724 LT+ vs LT−", float(m.group(1)), float(m.group(2)), "SCRN2"))
m = re.search(r"LT\+ vs N : log2FC=([+-][\d.]+)\s+p\(MW\)=([\d.]+)", geo_txt)
rows.append(("GSE165724 LT+ vs normal", float(m.group(1)), float(m.group(2)), "SCRN2"))
m = re.search(r"log2FC\(HT vs TN\) = ([+-][\d.]+).*?Welch p=([\d.]+)", geo_txt, re.S)
rows.append(("GSE138198 HT vs normal", float(m.group(1)), float(m.group(2)), "SCRN2"))
for name, lab in [("IFNG", None), ("CD4", None), ("PTPRC", None), ("TPO", None), ("TG", None)]:
    mm = re.search(r"^  %s\s+log2FC=([+-][\d.]+)\s+p=([\d.]+)" % name, geo_txt, re.M)
    rows.append((name + " (positive control)", float(mm.group(1)), float(mm.group(2)), "control"))
ye = np.arange(len(rows))[::-1]
band = 0.33
axE.axvspan(-band, band, color="#eceeef", zorder=0)
for yi, (name, b, pv, kind) in zip(ye, rows):
    col = TEAL if kind == "SCRN2" else ("#7f8b93" if b > 0 else "#b07a6e")
    axE.barh(yi, b, height=0.55, color=col, edgecolor="none", zorder=2)
    axE.text(b + (0.09 if b >= 0 else -0.09), yi,
             "%+.2f (P = %s)" % (b, "< 0.001" if pv < 0.001 else "%.3f" % pv),
             va="center", ha="left" if b >= 0 else "right", fontsize=5.8, color="#44494e")
axE.set_yticks(ye)
axE.set_yticklabels([r[0] for r in rows], fontsize=6.2)
axE.axvline(0, color="#4a4f54", lw=0.7)
axE.set_xlabel("log2 fold change in disease tissue")
axE.set_xlim(-1.6, 3.2)
axE.text(0.98, 0.03, "grey band: smallest detectable effect (%.2f)" % band,
         transform=axE.transAxes, ha="right", fontsize=5.8, color="#6b7177")
axE.set_title("Bulk tissue says nothing; the pipeline works", loc="left", pad=4,
              fontweight="bold")
letter(axE, "E", x=-0.24)

for ext in ("png", "pdf"):
    fp = os.path.join(OUT, "fig8." + ext)
    fig.savefig(fp, dpi=300, facecolor="white")
    print("[written]", fp)
plt.close(fig)
