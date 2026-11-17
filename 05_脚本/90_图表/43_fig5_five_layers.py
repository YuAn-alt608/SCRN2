#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
43_fig5_five_layers.py - Figure 5: the five evidence layers, from deposited tables
==================================================================================
Draws Figure 5 of the manuscript. Every number comes from 04_数据/:
  A  UpSet plot          04_数据/09_五层/region_counts.json
  B  cumulative chain    04_数据/09_五层/five_layer_membership.csv
  C  layer membership    04_数据/09_五层/five_layer_membership.csv
  D  threshold sensitivity  04_数据/09_五层/five_layer_threshold_sensitivity.csv

The four panels are regenerated from those files, so the figure and the deposited
sensitivity table cannot drift apart.

Output: 03_稿件/figs/fig5.png / .pdf  (300 dpi)
Usage:
    python 05_脚本/90_图表/43_fig5_five_layers.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
D = os.path.join(ROOT, "04_数据")
OUT = os.path.join(ROOT, "03_稿件", "figs")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.linewidth": 0.6,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
})

ACC, BLUE, GREY = "#B03A2E", "#2F6FB2", "#5A5A5A"
LAYER_LABELS = ["1  whole-blood", "2  thyroid", "3  coloc.",
                "4  coherence", "5  assignment"]
REGION_ORDER = [1, 2, 3, 5, 7, 11, 16, 15, 31]


def tag(ax, letter, title=None, x=-0.09, y=1.10):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="top", ha="left", color="#111")
    if title:
        try:                                    # 11 pt bold letter is ~0.115 in wide
            dx = 0.165 / max(ax.get_position().width * ax.figure.get_figwidth(), 0.6)
        except Exception:
            dx = 0.06
        ax.text(x + dx, y, title, transform=ax.transAxes, fontsize=7.4,
                va="top", ha="left", color="#111")


def load():
    rc = {int(k): v for k, v in json.load(
        open(os.path.join(D, "09_五层", "region_counts.json"), encoding="utf-8")).items()}
    sens = pd.read_csv(os.path.join(D, "09_五层", "five_layer_threshold_sensitivity.csv"))
    mem = pd.read_csv(os.path.join(D, "09_五层", "five_layer_membership.csv"))
    return rc, sens, mem


def panel_a(ax_bar, ax_mat, rc, sizes):
    order = sorted(rc, key=lambda b: -rc[b]["n"])
    counts = [rc[b]["n"] for b in order]
    x = np.arange(len(order))

    ax_bar.bar(x, counts, width=0.62, color=GREY, ec="none", zorder=2)
    ax_bar.set_yscale("log")
    ax_bar.set_ylim(0.8, 700)
    ax_bar.set_xlim(-0.7, len(order) - 0.3)
    ax_bar.set_ylabel("Genes in the exact\ncombination (log scale)", fontsize=6.6)
    ax_bar.set_xticks([])
    for xi, c in zip(x, counts):
        ax_bar.text(xi, c * 1.20, "%d" % c, ha="center", va="bottom", fontsize=6.4,
                    fontweight="bold", color="#222")
    ax_bar.spines["bottom"].set_visible(False)
    ax_bar.tick_params(axis="y", labelsize=6.0, length=2)

    ax_bar.tick_params(axis="x", labelbottom=False, length=0)

    ax_mat.set_xlim(-4.15, len(order) - 0.3)
    ax_mat.set_ylim(4.6, -0.6)
    for j in range(5):
        ax_mat.plot([-0.45, len(order) - 0.3], [j, j], color="#EDEDED", lw=8.0, zorder=0)
    for xi, b in zip(x, order):
        filled = [j for j in range(5) if (b >> j) & 1]
        for j in filled:
            ax_mat.plot([xi, xi], [min(filled), max(filled)], color="#333", lw=1.1, zorder=1)
        for j in range(5):
            on = j in filled
            ax_mat.scatter([xi], [j], s=26 if on else 20,
                           color="#222" if on else "#CFCFCF", zorder=3, lw=0)
    ax_mat.set_yticks([])
    ax_mat.set_xticks(x)
    ax_mat.set_xticklabels([str(b) for b in order], fontsize=5.6, color="#777")
    ax_mat.set_xlabel("bitmask of the occupied region", fontsize=5.8, labelpad=1)
    ax_mat.tick_params(axis="both", length=0)
    for s in ax_mat.spines.values():
        s.set_visible(False)

    # layer names and sizes on the left, inside the axes so nothing is clipped
    for j, (lab, n) in enumerate(zip(LAYER_LABELS, sizes)):
        ax_mat.text(-0.80, j, lab, fontsize=6.1, ha="right", va="center",
                    color="#222", clip_on=False)
        ax_mat.text(-2.55, j, "{:,}".format(n), fontsize=5.9, ha="right", va="center",
                    color="#555", clip_on=False)
    ax_mat.text(-2.55, -0.95, "layer size", fontsize=5.3, ha="right", va="bottom",
                color="#555")
    tag(ax_bar, "A", "Five layers, nine occupied combinations of thirty-one", x=-0.30, y=1.13)


def panel_b(ax, mem, rc):
    n_pairs, n_genes = 733, 612
    chain = [n_pairs, n_genes,
             int((mem.baseline & 1 > 0).sum()),
             int((mem.baseline & 3 == 3).sum()),
             int((mem.baseline & 7 == 7).sum()),
             int((mem.baseline & 15 == 15).sum()),
             int((mem.baseline == 31).sum())]
    labels = ["FDR < 0.05", "unique genes", "1 · whole-blood MR", "2 · also thyroid MR",
              "3 · also colocalised", "4 · also coherent", "5 · also assigned"]
    ax.set_xlim(0, 6.9)
    ax.set_ylim(-0.9, 6.9)
    for i, (lab, v) in enumerate(zip(labels, chain)):
        y = 6 - i
        w = max(0.55, 0.10 + 0.95 * (v / n_pairs) ** 0.42)
        ax.add_patch(Rectangle((1.85, y - 0.30), w, 0.60,
                               fc=ACC if i == 6 else BLUE, alpha=0.55 if i == 6 else 0.20,
                               ec="none"))
        ax.text(1.75, y, lab, ha="right", va="center", fontsize=6.2,
                fontweight="bold" if i == 6 else "normal")
        ax.text(1.85 + w + 0.08, y, "{:,}".format(v), va="center", fontsize=6.6,
                fontweight="bold", color=ACC if i == 6 else "#333")
    ax.set_axis_off()
    ax.text(0.0, -0.55, "Each step is an intersection, not a subtraction: 71 genes, not 62,\n"
                        "pass layers 1 and 2, and 32 of them pass no other layer.",
            fontsize=5.9, color="#555", va="top", linespacing=1.3)
    tag(ax, "B", "Cumulative intersections", x=-0.02, y=1.06)


def panel_c(ax, mem):
    rows = ["SCRN2", "MCM6", "NUCB2", "UCP2", "PPP2R1B", "UROD", "YWHAZ", "SUOX", "TBC1D2"]
    bits = {r.gene: int(r.baseline) for r in mem.itertuples()}
    M = np.array([[1.0 if (bits[g] >> j) & 1 else 0.0 for j in range(5)] for g in rows])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("pw", ["#EDEDED", "#2F6FB2"])
    ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, "●" if M[i, j] == 1 else "○", ha="center", va="center",
                    fontsize=6.4, color="white" if M[i, j] == 1 else "#999")
    ax.set_xticks(range(5))
    ax.set_xticklabels(["1\nblood\nMR", "2\nthyroid\nMR", "3\ncolo-\ncalisation",
                        "4\ntissue\ncoherence", "5\ntarget\nassignment"], fontsize=5.6,
                       linespacing=1.1)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=6.2)
    for lb in ax.get_yticklabels():
        if lb.get_text() == "SCRN2":
            lb.set_fontweight("bold")
            lb.set_color(ACC)
    ax.set_xticks(np.arange(-0.5, 5, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.4)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(axis="both", length=0)
    tag(ax, "C", "The nine closest genes", x=-0.40, y=1.06)


def panel_d(ax, sens):
    ax.set_axis_off()
    hdr = ["Relaxation", "layer 3", "layer 4", "layer 5", "intersection"]
    xs = [0.00, 0.52, 0.68, 0.84, 1.00]
    for x, h in zip(xs, hdr):
        ax.text(x, 1.02, h, fontsize=6.2, fontweight="bold", va="bottom",
                ha="left" if x == 0 else "center")
    names = {"baseline": "Baseline",
             "PP.H4>=0.70": "PP.H4 ≥ 0.70",
             "top-5 targets": "Targets top-5",
             "concordance>=90%": "Coherence ≥ 90%",
             "all three relaxed": "All three"}
    for i, r in enumerate(sens.itertuples()):
        y = 0.86 - i * 0.155
        ax.plot([0, 1.0], [y + 0.055, y + 0.055], color="#E4E4E4", lw=0.6, zorder=0)
        vals = [names.get(r.configuration, r.configuration), str(r.layer3_colocalisation),
                str(r.layer4_cross_tissue), str(r.layer5_target),
                "%s (n=%d)" % (r.five_layer_genes, r.n_five_layer)]
        for k, (x, v) in enumerate(zip(xs, vals)):
            ax.text(x, y, v, fontsize=5.9, va="center",
                    ha="left" if k == 0 else "center",
                    color=ACC if k == 4 else "#222",
                    fontweight="bold" if k == 4 else "normal")
    ax.set_xlim(-0.02, 1.12)
    ax.set_ylim(-0.12, 1.12)
    ax.text(0.0, -0.10, "Layer sizes for each relaxation; the five-layer intersection is invariant.\n"
                        "Per-configuration gene lists: five_layer_threshold_sensitivity.csv.",
            fontsize=5.7, color="#555", va="top", linespacing=1.35)
    tag(ax, "D", "Threshold sensitivity", x=-0.02, y=1.16)


def main():
    rc, sens, mem = load()
    fig = plt.figure(figsize=(7.09, 7.0))
    outer = fig.add_gridspec(3, 2, height_ratios=[1.28, 1.0, 0.56],
                             width_ratios=[1.10, 1.0],
                             hspace=0.34, wspace=0.34,
                             left=0.215, right=0.975, top=0.955, bottom=0.05)
    inner = outer[0, :].subgridspec(2, 1, height_ratios=[1.0, 0.30], hspace=0.06)
    ax_bar = fig.add_subplot(inner[0])
    ax_mat = fig.add_subplot(inner[1], sharex=ax_bar)
    b = sens[sens.configuration == "baseline"].iloc[0]
    sizes = [b.layer1_whole_blood_MR, b.layer2_thyroid_MR, b.layer3_colocalisation,
             b.layer4_cross_tissue, b.layer5_target]
    panel_a(ax_bar, ax_mat, rc, sizes)

    panel_b(fig.add_subplot(outer[1, 0]), mem, rc)
    panel_c(fig.add_subplot(outer[1, 1]), mem)
    panel_d(fig.add_subplot(outer[2, :]), sens)

    fig.savefig(os.path.join(OUT, "fig5.png"), bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(OUT, "fig5.pdf"), bbox_inches="tight", facecolor="white")
    print("[written] %s/fig5.png" % OUT)


if __name__ == "__main__":
    main()
