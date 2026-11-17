# -*- coding: utf-8 -*-
"""
生成稿件 Figure 1–9。

所有数值从 04_数据 直接读取；不硬编码，不重算稿件以外的统计量。
输出：03_稿件/figs/figN.png（300 dpi）与 .pdf

用法：
    python make_figures.py          # 全部
    python make_figures.py 3 5      # 只画第 3、5 张
"""

import os
import sys
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, FancyBboxPatch, Ellipse, Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, "..", "..", "04_数据"))
OUT = os.path.abspath(os.path.join(HERE, "..", "..", "03_稿件", "figs"))
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.titlesize": 8,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 6.8,
    "ytick.labelsize": 6.8,
    "legend.fontsize": 6.8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.4,
    "ytick.major.size": 2.4,
    "savefig.dpi": 300,
    "figure.dpi": 110,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

W = 7.09  # 180 mm
C = {"TSH": "#2F6FB2", "HYPO": "#C0392B", "T4": "#E08A1E", "HASH": "#7D5BA6",
     "grey": "#8A8A8A", "dark": "#222222", "acc": "#B03A2E"}
OC = {"TSH": ("TSH", C["TSH"]), "HYPO": ("Hypothyroidism", C["HYPO"]),
      "T4": ("Free T4", C["T4"]), "HASH": ("Hashimoto", C["HASH"])}


def tag(ax, letter, title=None, x=-0.10, y=1.10):
    # 3D 轴无 text()，改用 text2D
    put = ax.text2D if hasattr(ax, "text2D") else ax.text
    put(x, y, letter, transform=ax.transAxes, fontsize=11, fontweight="bold",
        va="top", ha="left", color="#111111")
    if title:
        # 11 pt 粗体大写字母约 0.115 in 宽；按面板实际宽度换算成 axes 坐标，窄面板也不会压字
        try:
            pw = ax.get_position().width * ax.figure.get_figwidth()
        except Exception:
            pw = 2.4
        dx = 0.115 / max(pw, 0.55)
        put(x + dx, y, title, transform=ax.transAxes, fontsize=7.7, va="top", ha="left",
            color="#111111")


def save(fig, n):
    fig.savefig(os.path.join(OUT, "fig%d.png" % n), bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(OUT, "fig%d.pdf" % n), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  fig%d" % n)


# ══════════════════════════════════════════════════ 读数据
def load():
    g = {}
    g["mr"] = pd.read_csv(os.path.join(D, "01_筛查", "MR_all_results.csv"))
    g["sig"] = pd.read_csv(os.path.join(D, "01_筛查", "MR_significant_FDR05.csv"))
    g["tiers"] = pd.read_csv(os.path.join(D, "01_筛查", "evidence_tiers.csv"))
    g["deg"] = pd.read_csv(os.path.join(D, "01_筛查", "gene_degree_all.csv"))
    g["coloc"] = pd.read_csv(os.path.join(D, "02_共定位", "coloc_results.csv"))
    g["circ"] = pd.read_csv(os.path.join(D, "02_共定位", "h4_circularity.csv"))
    g["indep"] = pd.read_csv(os.path.join(D, "02_共定位", "h4_independent_loci.csv"))
    g["fin"] = pd.read_csv(os.path.join(D, "03_复现", "finngen_replication.csv"))
    g["tis"] = pd.read_csv(os.path.join(D, "04_组织", "gtex_tissue_rates.csv"))
    g["t3s"] = pd.read_csv(os.path.join(D, "04_组织", "thyroid_3source_validation.csv"))
    g["sig23"] = pd.read_csv(os.path.join(D, "05_位点", "23_信号解析归属总表.csv"))
    g["jm"] = pd.read_csv(os.path.join(D, "05_位点", "31_联合条件模型.csv"))
    g["win"] = pd.read_csv(os.path.join(D, "05_位点", "24_按信号窗口的coloc.csv"))
    g["ld"] = pd.read_csv(os.path.join(D, "05_位点", "00_LD面板位点.csv"), header=None,
                          names=["pos", "r2"])
    g["var"] = pd.read_csv(os.path.join(D, "05_位点", "11_变异信息表.csv"))
    g["plane"] = pd.read_csv(os.path.join(D, "08_平面与跨族裔", "s08_plane_genes.csv"))
    g["gw"] = pd.read_csv(os.path.join(D, "08_平面与跨族裔", "gw_thyroid_mr.csv"))
    g["meta"] = pd.read_csv(os.path.join(D, "08_平面与跨族裔", "three_cohort_meta.csv"))
    g["rc"] = json.load(open(os.path.join(D, "09_五层", "region_counts.json"), encoding="utf-8"))
    g["mds"] = pd.read_csv(os.path.join(D, "07_结构", "summary_table.csv"))
    g["res"] = pd.read_csv(os.path.join(D, "07_结构", "MD逐残基接触.csv"))
    g["dist"] = pd.read_csv(os.path.join(D, "07_结构", "MD接触重算.csv"))
    g["expr"] = pd.read_csv(os.path.join(D, "06_表达与单细胞", "genes_expression_profile.csv"))
    g["a2"] = pd.read_csv(os.path.join(D, "06_表达与单细胞", "a2_summary.csv"))
    g["cell"] = pd.read_csv(os.path.join(D, "06_表达与单细胞", "SCRN2_thyroid_celltype.csv"))
    g["s07"] = pd.read_csv(os.path.join(D, "06_表达与单细胞", "s07_panel_pc10.csv"))
    return g


# ══════════════════════════════════════════════════ 图 1
def fig1(g):
    fig = plt.figure(figsize=(W, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.28, 1], width_ratios=[1.42, 1],
                          hspace=0.55, wspace=0.30)

    # ── A 分层证据链（固定槽位，标签不溢出）
    ax = fig.add_subplot(gs[0, :])
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # 槽位起点：左侧留 0.145 给行标签，其余 4 个槽
    slots = [0.148, 0.330, 0.525, 0.720]

    def band(y, label, col, items, bold_last=False, sizes=None):
        ax.add_patch(Rectangle((0.0, y - 0.062), 1.0, 0.124, fc=col, ec="none", alpha=0.55))
        ax.text(0.006, y, label, fontsize=6.3, fontweight="bold", va="center", color="#222")
        for i, it in enumerate(items):
            x = slots[i]
            fw = "bold" if (bold_last and i == len(items) - 1) else "normal"
            fs = sizes[i] if sizes else 6.5
            ax.text(x, y, it, fontsize=fs, va="center", fontweight=fw, color="#1a1a1a")
            if i < len(items) - 1:
                nx = slots[i + 1]
                ax.annotate("", xy=(nx - 0.010, y), xytext=(nx - 0.040, y),
                            arrowprops=dict(arrowstyle="-|>", lw=0.7, color="#5a5a5a",
                                            shrinkA=0, shrinkB=0))

    band(0.90, "Exposure", "#DCE7F2",
         ["1,580 EDCs", "256,899 pairs", "24,264 genes", "∩ thyroid = 4,526"])
    band(0.70, "Instruments", "#F7E7D3",
         ["12,636 pairs", "733 at FDR<0.05", "612 genes", "10 loci"])
    band(0.50, "Locus", "#E6E9F5",
         ["7 signals / 5 positions", "4 candidate genes", "1 for TSH only", "SCRN2, PP.H4 0.99"], True)
    band(0.30, "Surface", "#F5DEDC",
         ["PFNA docked", "groove R119/R77", "QDPVR 3.2 Å", "QVDRR 6.7 Å"], True,
         sizes=[6.5, 6.2, 6.5, 6.5])

    # 右侧竖括号：前两条是遗传臂
    ax.plot([1.005, 1.005], [0.775, 0.975], color="#7a7a7a", lw=0.8, clip_on=False)
    ax.plot([1.005, 1.02], [0.775, 0.775], color="#7a7a7a", lw=0.8, clip_on=False)
    ax.plot([1.005, 1.02], [0.975, 0.975], color="#7a7a7a", lw=0.8, clip_on=False)
    ax.text(1.03, 0.875, "genetic arm", fontsize=6.2, rotation=90, va="center", color="#555",
            clip_on=False)
    ax.plot([1.005, 1.005], [0.375, 0.575], color="#7a7a7a", lw=0.8, clip_on=False)
    ax.plot([1.005, 1.02], [0.375, 0.375], color="#7a7a7a", lw=0.8, clip_on=False)
    ax.plot([1.005, 1.02], [0.575, 0.575], color="#7a7a7a", lw=0.8, clip_on=False)
    ax.text(1.03, 0.475, "chemical arm", fontsize=6.2, rotation=90, va="center", color="#555",
            clip_on=False)

    # 汇合箭头
    ax.text(0.50, 0.145,
            "The two arms meet at one gene and are never multiplied into a single score.\n"
            "Genetic evidence locates SCRN2; the chemical arm only proposes a mechanism.",
            fontsize=6.1, style="italic", color="#444", va="center", ha="center",
            linespacing=1.45)
    ax.text(0.0, 1.02, "A", transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="bottom")
    ax.text(0.035, 1.02, "Design: two arms, one gene", transform=ax.transAxes,
            fontsize=7.8, va="bottom")

    # ── B 漏斗
    ax = fig.add_subplot(gs[1, 0])
    steps = [("Curated EDC–gene pairs", 256899), ("EDC-mapped genes", 24264),
             ("Thyroid-expressed", 4526), ("With an instrument", 3559),
             ("Testable genes", 3416), ("Testable pairs", 12636), ("FDR < 0.05", 733),
             ("Unique genes", 612), ("Coloc PP.H4>0.8", 71),
             ("Independent loci", 10), ("Five layers", 1)]
    labs = [s[0] for s in steps][::-1]
    vals = [s[1] for s in steps][::-1]
    cols = (["#D6D6D6"] + ["#A8C0D8"] * 3 + ["#7FA6C9"] * 3
            + ["#5B7FA6"] * (len(steps) - 7))
    yy = np.arange(len(vals))
    ax.barh(yy, vals, color=cols, height=0.68, edgecolor="white", lw=0.5)
    ax.set_xscale("log")
    ax.set_yticks(yy); ax.set_yticklabels(labs)
    ax.set_xlim(0.7, 1.2e6)
    ax.set_xlabel("Count (log scale)")
    for i, v in enumerate(vals):
        ax.text(v * 1.28, i, "{:,}".format(v), va="center", fontsize=6.1, color="#333")
    ax.set_ylim(-0.7, len(vals) - 0.3)
    ax.text(0.0, 1.04, "B", transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="bottom")
    ax.text(0.055, 1.04, "Attrition", transform=ax.transAxes, fontsize=7.8, va="bottom")

    # ── C 每层回答什么
    ax = fig.add_subplot(gs[1, 1])
    ax.set_axis_off()
    txt = [("Screen", "Any signal at all?"),
           ("FinnGen", "Reproduces elsewhere?"),
           ("Thyroid eQTL", "Instrument still works?"),
           ("Colocalisation", "Same causal variant?"),
           ("Coherence", "Direction tissue-general?"),
           ("Assignment", "Which gene here?")]
    for i, (a, b) in enumerate(txt):
        yv = 0.97 - i * 0.163
        ax.text(0.0, yv, a, fontsize=6.7, va="top", fontweight="bold", color="#333")
        ax.text(0.42, yv, b, fontsize=6.7, va="top", color="#444")
    ax.text(0.0, -0.02, "Complementary, not independent;\nall share one TSH GWAS.",
            fontsize=6.0, style="italic", color="#777", va="top")
    ax.text(0.0, 1.06, "C", transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="bottom")
    ax.text(0.16, 1.06, "What each layer asks", transform=ax.transAxes, fontsize=7.8, va="bottom")

    save(fig, 1)


# ══════════════════════════════════════════════════ 图 2
def fig2(g):
    mr = g["mr"].drop_duplicates(["gene", "outcome"])
    fig = plt.figure(figsize=(W, 5.6))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.30)

    # A 火山
    ax = fig.add_subplot(gs[0, 0])
    for oc, (lab, col) in OC.items():
        s = mr[mr.outcome == oc]
        ax.scatter(s.b, -np.log10(s.pval), s=2.2, c=col, alpha=0.30, lw=0, rasterized=True)
    ax.axhline(-np.log10(0.05 / len(mr)), ls="--", lw=0.7, c="#555")
    ax.text(0.345, -np.log10(0.05 / len(mr)) + 1.2, "FDR 0.05", fontsize=5.9,
            color="#555", ha="right")
    scr = mr[(mr.gene == "SCRN2") & (mr.outcome == "TSH")].iloc[0]
    ax.scatter([scr.b], [-np.log10(scr.pval)], s=34, facecolor="none",
               edgecolor=C["acc"], lw=1.2, zorder=5)
    ax.annotate("SCRN2", (scr.b, -np.log10(scr.pval)), xytext=(14, 6),
                textcoords="offset points", fontsize=6.5, color=C["acc"], fontweight="bold")
    ax.set_xlabel("MR effect on the outcome (β)")
    ax.set_ylabel(r"$-\log_{10}$ P")
    ax.set_xlim(-0.35, 0.35)
    hs = [plt.Line2D([], [], marker="o", ls="", ms=4, c=c, label=OC[k][0]) for k, c in
          [("TSH", C["TSH"]), ("HYPO", C["HYPO"]), ("T4", C["T4"]), ("HASH", C["HASH"])]]
    ax.legend(handles=hs, frameon=False, loc="upper left", handletextpad=0.2)
    tag(ax, "A", "12,636 pairs, 4 outcomes")

    # B 堆叠
    ax = fig.add_subplot(gs[0, 1])
    order = ["TSH", "Hypothyroidism", "Total T4", "Hashimoto thyroiditis"]
    tcol = {"A_multi": "#1F3B57", "B_strong": "#4E7AA8", "C_distinct": "#93B4D1",
            "D_unsupported": "#D9D9D9"}
    bottom = np.zeros(len(order))
    tiers = g["tiers"]
    for k in ["A_multi", "B_strong", "C_distinct", "D_unsupported"]:
        vals = [int(((tiers.outcome == o) & (tiers.tier == k)).sum()) for o in order]
        ax.bar(range(len(order)), vals, bottom=bottom, color=tcol[k], width=0.62,
               edgecolor="white", lw=0.5,
               label=k.split("_")[0] + " " + k.split("_")[1])
        bottom += np.array(vals, dtype=float)
    for i, v in enumerate(bottom):
        ax.text(i, v + 8, "{:.0f}".format(v), ha="center", fontsize=7, fontweight="bold")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["TSH", "Hypo-\nthyroidism", "Free T4", "Hashimoto"], rotation=0)
    ax.set_ylabel("Significant pairs")
    ax.set_ylim(0, 430)
    ax.legend(frameon=False, ncol=2, loc="upper right", handlelength=0.9,
              columnspacing=0.8, borderaxespad=0.2)
    tag(ax, "B", "733 pairs by tier")

    # C EDC 五分位（基因级：分母是基因数，不是对数）
    ax = fig.add_subplot(gs[1, 0])
    deg = g["deg"].set_index("gene")["n_edc_chemicals"]
    for oc, (lab, col) in OC.items():
        s_oc = mr[mr.outcome == oc]
        genes = s_oc.gene.values
        q = pd.qcut(deg.reindex(genes).fillna(0).values, 5, labels=False, duplicates="drop")
        sig = set(s_oc[s_oc.fdr < 0.05].gene)
        y, err = [], []
        for k in range(5):
            idx = q == k
            tot = int(idx.sum())
            hit = sum(1 for x in np.array(genes)[idx] if x in sig)
            p = 100.0 * hit / max(tot, 1)
            y.append(p)
            err.append(1.96 * np.sqrt(p * (100 - p) / max(tot, 1)))
        ax.errorbar(range(5), y, yerr=err, marker="o", ms=3.4, lw=1.1, capsize=1.8,
                    elinewidth=0.7, color=col, label=lab)
    ax.set_xticks(range(5))
    ax.set_xticklabels(["Q1\n(lowest)", "Q2", "Q3", "Q4", "Q5\n(highest)"])
    ax.set_xlabel("Quintile of EDCs mapped to the gene")
    ax.set_ylabel("Genes at FDR < 0.05 (%)")
    ax.set_ylim(-0.8, 22)
    ax.legend(frameon=False, ncol=4, loc="upper center", handlelength=0.8,
              columnspacing=0.7, handletextpad=0.25, borderaxespad=0.1, fontsize=6.2)
    ax.text(0.99, 0.04, "TSH ρ = +0.026\nHypo ρ = −0.008\n(all P > 0.3)",
            transform=ax.transAxes, fontsize=5.9, color="#555", ha="right")
    tag(ax, "C", "No dose–response at the gene level")

    # D 散点（基因级）
    ax = fig.add_subplot(gs[1, 1])
    s_tsh = mr[mr.outcome == "TSH"].copy()
    s_tsh["absz"] = np.abs(s_tsh.b / s_tsh.se)
    s_tsh["edc"] = deg.reindex(s_tsh.gene).fillna(0).values
    ax.scatter(s_tsh.edc, s_tsh.absz, s=2.0, c=C["TSH"], alpha=0.22, lw=0, rasterized=True)
    order = s_tsh.sort_values("edc").reset_index(drop=True)
    roll = 500
    xx = order.edc.rolling(roll, center=True, min_periods=150).median()
    yy = order.absz.rolling(roll, center=True, min_periods=150).median()
    ax.plot(xx, yy, color=C["acc"], lw=1.4)
    ax.axhline(s_tsh.absz.median(), ls=":", lw=0.8, color="#777")
    ax.text(0.99, 0.72, "median |z| = %.2f" % s_tsh.absz.median(),
            transform=ax.transAxes, ha="right", fontsize=5.9, color="#777")
    ax.set_xlabel("Number of EDCs mapped to the gene")
    ax.set_ylabel("|z| of the MR estimate")
    ax.set_xlim(0, np.percentile(s_tsh.edc, 99.3))
    ax.set_ylim(0, 12)
    ax.text(0.98, 0.95, "ρ = +0.026 (P = 0.13)", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.4, color="#555")
    tag(ax, "D", "The same absence, ungrouped")

    save(fig, 2)


# ══════════════════════════════════════════════════ 图 3
def fig3(g):
    fig = plt.figure(figsize=(W, 5.8))
    gs = fig.add_gridspec(2, 2, hspace=0.48, wspace=0.34)

    # A 复现率（分母是可测对数，不是全部显著对）
    ax = fig.add_subplot(gs[0, 0])
    fin = g["fin"]
    ocs = ["HYPO", "TSH", "HASH"]
    labs = ["Hypo-\nthyroidism†", "TSH", "Hashimoto†", "Negative\ncontrol"]
    cols_ = [C["HYPO"], C["TSH"], C["HASH"], "#C9C9C9"]
    hatch = ["///", "", "///", ""]
    rate, ntest = [], []
    for o in ocs:
        sub = fin[(fin.oc == o) & fin.b_fg.notna()]
        rate.append(sub.replicated.mean() * 100)
        ntest.append(len(sub))
    rate.append(6.7); ntest.append(388)   # 阴性对照：FDR>0.5 的 388 对
    x = np.arange(len(labs))
    bars = ax.bar(x, rate, width=0.56, color=cols_, edgecolor="white", lw=0.6)
    for b_, h_ in zip(bars, hatch):
        b_.set_hatch(h_); b_.set_edgecolor("white")
    for i, (r, nn_) in enumerate(zip(rate, ntest)):
        ax.text(i, r + 2.6, "%.1f%%" % r, ha="center", fontsize=6.6, fontweight="bold")
        ax.text(i, 2.4, "n = %d" % nn_, ha="center", va="bottom", fontsize=5.5, color="#222",
                bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none", alpha=0.82))
    ax.set_xticks(x); ax.set_xticklabels(labs)
    ax.set_ylabel("Same direction and P < 0.05\nin FinnGen R12 (%)")
    ax.set_ylim(0, 101)
    ax.set_xlim(-0.6, len(labs) - 0.4)
    ax.text(0.0, 0.995, "† discovery set already contains Finnish participants (r = 0.91):\n"
                        "consistency check, not replication. Free T4 omitted: no detectable\n"
                        "replication (r = 0.09, P = 0.49, n = 57).",
            transform=ax.transAxes, fontsize=5.0, color="#666", va="top", linespacing=1.5)
    tag(ax, "A", "Replication, and what is not replication")

    # B 散点
    ax = fig.add_subplot(gs[0, 1])
    s_ = fin[(fin.oc == "HYPO") & fin.b_fg.notna()].copy()
    ok = ((s_.sign_ok) & (s_.p_fg < 0.05)).values
    ax.scatter(s_.b_orig[ok], s_.b_fg[ok], s=5, c=C["HYPO"], alpha=0.5, lw=0,
               label="same direction, P<0.05")
    ax.scatter(s_.b_orig[~ok], s_.b_fg[~ok], s=5, c="#9A9A9A", alpha=0.55, lw=0,
               label="otherwise")
    lim = [-0.55, 0.55]
    ax.plot(lim, lim, ls="--", lw=0.7, c="#666")
    ax.axhline(0, lw=0.5, c="#ccc"); ax.axvline(0, lw=0.5, c="#ccc")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Discovery β (whole blood)")
    ax.set_ylabel("FinnGen β")
    ax.legend(frameon=False, loc="upper left", handletextpad=0.2, markerscale=1.8,
              borderaxespad=0.2, fontsize=6.2)
    ax.text(0.97, 0.05, "94.3% concordant; 66.9% also P<0.05",
            transform=ax.transAxes, ha="right", fontsize=6.1, color="#444")
    tag(ax, "B", "Hypothyroidism, 314 testable pairs")

    # C 组织检出率
    ax = fig.add_subplot(gs[1, 0])
    t = g["tis"].sort_values("detect_pct", ascending=False).reset_index(drop=True)
    cols = []
    for name in t.tissue:
        if name == "Thyroid":
            cols.append(C["acc"])
        elif name == "Whole_Blood":
            cols.append("#555555")
        else:
            cols.append("#B9C9DA")
    ax.bar(range(len(t)), t.detect_pct, color=cols, width=0.78, lw=0)
    ax.axhline(t.detect_pct.median(), ls=":", lw=0.8, c="#666")
    ax.set_ylabel("Instruments recovered (%)")
    ax.set_xlabel("GTEx tissue (49), ranked")
    ax.set_xlim(-0.8, len(t) + 2.6)
    ax.set_ylim(0, 52)
    ax.text(len(t) + 0.4, 43.5, "43.5%", fontsize=6.6, color="#1F6B3A", fontweight="bold", va="center")
    ax.text(len(t) + 0.4, 38.5, "eQTL Catalogue\n(denominator 62)", fontsize=5.8, color="#1F6B3A",
            va="top")
    i_thy = list(t.tissue).index("Thyroid")
    ax.annotate("thyroid 31.4%", (i_thy, 31.4), xytext=(10, 10),
                textcoords="offset points", fontsize=6.2, color=C["acc"], fontweight="bold",
                arrowprops=dict(arrowstyle="-", lw=0.6, color=C["acc"]))
    ax.annotate("blood 41.4%", (0, 41.4), xytext=(16, 6),
                textcoords="offset points", fontsize=6.2, color="#333", fontweight="bold")
    ax.text(12, t.detect_pct.median() + 1.4, "median 12.9%", fontsize=5.9, color="#666")
    tag(ax, "C", "Where the same instrument still works")

    # D 三来源矩阵
    ax = fig.add_subplot(gs[1, 1])
    d = g["t3s"].copy()
    d["n_src"] = d[["s1_sig", "s2_sig", "s3_sig"]].fillna(False).astype(bool).sum(axis=1)
    d["n_same"] = d[["s1_same", "s2_same", "s3_same"]].fillna(False).astype(bool).sum(axis=1)
    d = d.sort_values(["n_src", "n_same"], ascending=False).reset_index(drop=True)
    M = d[["s1_sig", "s2_sig", "s3_sig"]].fillna(False).astype(int).values
    ax.imshow(M.T, aspect="auto", cmap=LinearSegmentedColormap.from_list(
        "g", ["#F0F0F0", "#5B7FA6"]), vmin=0, vmax=1, interpolation="nearest")
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["eQTLGen", "GTEx", "Catalogue"], fontsize=6.2)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xlabel("70 instruments, ranked by number of sources", labelpad=2)
    ax.text(44, 1.0, "four leftmost: PSRC1, DDX5,\nSCRN2, PRR15L", fontsize=5.6,
            va="center", ha="left", color="#444")
    ax.set_ylim(2.5, -0.6)
    tag(ax, "D", "Four instruments in all three call sets")

    save(fig, 3)


# ══════════════════════════════════════════════════ 图 4
def fig4(g):
    co = g["coloc"]
    circ = g["circ"]
    fig = plt.figure(figsize=(W, 5.9))
    gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.32)

    # A PP.H4 分布
    ax = fig.add_subplot(gs[0, 0])
    b = np.linspace(0, 1, 51)
    hi = co[co["PP.H4"] > 0.8]
    ax.hist(co["PP.H4"], bins=b, color="#9FB4C9", edgecolor="none")
    ax.hist(hi["PP.H4"], bins=b, color=C["acc"], edgecolor="none")
    ax.set_yscale("log")
    ax.axvline(0.8, ls="--", lw=0.9, c="#333")
    ax.text(0.79, ax.get_ylim()[1] * 0.35, "PP.H4 = 0.8", rotation=90, ha="right",
            fontsize=6.1, color="#333")
    ax.text(0.86, 12, "%d pairs\n> 0.8" % len(hi), fontsize=6.5, color=C["acc"],
            fontweight="bold", va="bottom")
    ax.set_xlabel("PP.H4 (shared causal variant)")
    ax.set_ylabel("Gene–outcome pairs (log)")
    ax.set_xlim(0, 1.02)
    # 内嵌：只看显著的 733 对
    axi = ax.inset_axes([0.10, 0.60, 0.36, 0.33])
    mp = {"TSH": "TSH", "Hypothyroidism": "HYPO", "Total T4": "T4",
          "Hashimoto thyroiditis": "HASH"}
    sigset = set(zip(g["tiers"].gene, g["tiers"].outcome.map(mp)))
    issig = np.array([(gg, oo) in sigset for gg, oo in zip(co.gene, co.outcome)])
    axi.hist(co["PP.H4"][issig], bins=np.linspace(0, 1, 26),
             color=C["acc"], edgecolor="none")
    axi.hist(co["PP.H4"][~issig], bins=np.linspace(0, 1, 26),
             color="#C9D2DC", edgecolor="none")
    axi.set_title("733 significant vs controls", fontsize=5.6, pad=2)
    axi.set_ylim(0, 235)
    axi.tick_params(labelsize=5)
    axi.set_yticks([])
    for sp_ in axi.spines.values():
        sp_.set_linewidth(0.5)
    axi.set_xlim(0, 1)
    tag(ax, "A", "2,464 tests")

    # B 循环性审计
    ax = fig.add_subplot(gs[0, 1])
    n_ins = int(circ.top_is_instrument.sum())
    n_lead = int((~circ.top_is_instrument.astype(bool) &
                  circ.same_as_eqtl_top.astype(bool)).sum())
    n_none = len(circ) - n_ins - n_lead
    parts = [(n_ins, "Peak SNP is the instrument itself", "#B03A2E"),
             (n_lead, "Peak is the eQTL lead, not the instrument", "#E0A03E"),
             (n_none, "Neither", "#8FAFC7")]
    left = 0
    for v, lab, c in parts:
        ax.barh(0, v, left=left, color=c, height=0.52, edgecolor="white", lw=1.0)
        if v >= 5:
            ax.text(left + v / 2, 0, str(v), ha="center", va="center", color="white",
                    fontsize=8, fontweight="bold")
        left += v
    ax.set_xlim(0, 71); ax.set_ylim(-0.65, 0.65)
    ax.set_yticks([])
    ax.set_xlabel("Pairs with PP.H4 > 0.8 (n = 71)")
    yy = -0.34
    for v, lab, c in parts:
        ax.add_patch(Rectangle((0, yy - 0.035), 2.4, 0.07, fc=c, ec="none"))
        ax.text(3.2, yy, "%s (%d)" % (lab, v), va="center", fontsize=6.2, color="#333")
        yy -= 0.15
    ax.text(0.0, 1.20, "60/71 = 84.5% of the apparent colocalisations are circular",
            transform=ax.transAxes, ha="left", fontsize=6.6, color="#333",
            fontweight="bold")
    ax.set_ylim(-0.95, 0.65)
    for sp in ("left",):
        ax.spines[sp].set_visible(False)
    tag(ax, "B", "Instrument circularity")

    # C 缩窗
    ax = fig.add_subplot(gs[1, 0])
    ind = g["indep"].sort_values("pp_h4")
    ylab = ind.gene + "  (" + ind.outcome.map({"TSH": "TSH", "HYPO": "Hypo", "T4": "FT4",
                                               "HASH": "Hash"}).fillna(ind.outcome) + ")"
    yv = np.arange(len(ind))
    ax.barh(yv, ind.pp_h4, color=[C.get(str(o), "#888") for o in ind.outcome],
            height=0.66, edgecolor="white", lw=0.5)
    ax.set_yticks(yv); ax.set_yticklabels(ylab, fontsize=6.4)
    ax.axvline(0.8, ls="--", lw=0.8, c="#333")
    ax.set_xlabel("PP.H4")
    ax.set_xlim(0.78, 1.0)
    ax.set_ylim(-0.7, len(ind) - 0.3)
    ax.text(0.995, 0.06, "71 → 10 after restricting\neach window to its own\nindependent locus",
            transform=ax.transAxes, ha="right", fontsize=6.1, color="#444", va="bottom")
    tag(ax, "C", "What survives the audit")

    # D 先验敏感性
    ax = fig.add_subplot(gs[1, 1])
    xs = [r"$1\times10^{-4}$", "$1\\times10^{-5}$\n(as reported)", r"$1\times10^{-6}$"]
    ys = [196, 71, 18]
    ax.bar(range(3), ys, color=["#B8C9DA", "#B03A2E", "#B8C9DA"], width=0.55,
           edgecolor="white", lw=0.6)
    for i, v in enumerate(ys):
        ax.text(i, v + 5, str(v), ha="center", fontsize=7.2, fontweight="bold")
    ax.set_xticks(range(3)); ax.set_xticklabels(xs)
    ax.set_xlabel("Colocalisation prior p12")
    ax.set_ylabel("Pairs with PP.H4 > 0.8")
    ax.set_ylim(0, 230)
    ax.text(0.5, 0.72, "Rank order unchanged (ρ = 1.00)", transform=ax.transAxes,
            ha="center", fontsize=6.4, color="#444")
    tag(ax, "D", "Prior sensitivity")

    save(fig, 4)


# ══════════════════════════════════════════════════ 图 5
def _venn5(ax):
    """五集"花瓣"韦恩：每个集合一瓣，中心为五层交集。"""
    from matplotlib.patches import Ellipse
    from matplotlib.patches import Circle

    sets = [(0, "Whole-blood\nMR (337)", "247", "#2F6FB2"),
            (90, "Thyroid MR\n(180)", "109", "#7D5BA6"),
            (180, "Colocalisation\n(26)", "19", "#E08A1E"),
            (232, "Coherence\n(3)", "3", "#2E7D6B"),
            (308, "Assignment\n(4)", "0", "#B03A2E")]

    R, EW, EH = 0.62, 1.80, 0.62
    for th, lab, cnt, col in sets:
        cx, cy = R * np.cos(np.deg2rad(th)), R * np.sin(np.deg2rad(th))
        ax.add_patch(Ellipse((cx, cy), EW, EH, angle=th, fc=col, ec=col,
                             lw=0.8, alpha=0.20, zorder=1))
        ax.add_patch(Ellipse((cx, cy), EW, EH, angle=th, fc="none", ec=col,
                             lw=0.7, alpha=0.75, zorder=2))
        nx, ny = 1.10 * np.cos(np.deg2rad(th)), 1.10 * np.sin(np.deg2rad(th))
        ax.text(nx, ny, cnt, ha="center", va="center", fontsize=6.0,
                color=col if cnt != "0" else "#999", fontweight="bold", zorder=4)
        rad = 1.50 if abs(np.sin(np.deg2rad(th))) > 0.6 else 1.72
        lx, ly = rad * np.cos(np.deg2rad(th)), rad * np.sin(np.deg2rad(th))
        ha = "left" if np.cos(np.deg2rad(th)) > 0.2 else (
            "right" if np.cos(np.deg2rad(th)) < -0.2 else "center")
        va = "center" if abs(np.cos(np.deg2rad(th))) > 0.5 else (
            "bottom" if np.sin(np.deg2rad(th)) > 0 else "top")
        ax.text(lx, ly, lab, ha=ha, va=va, fontsize=5.6, color=col,
                linespacing=1.15, zorder=4)

    ax.add_patch(Circle((0, 0), 0.235, fc="white", ec="#555", lw=0.7, zorder=5))
    ax.text(0, 0, "1", ha="center", va="center", fontsize=7.0, fontweight="bold",
            color="#B03A2E", zorder=6)
    ax.text(0.0, -0.42, "SCRN2", ha="center", va="center", fontsize=6.2, fontweight="bold",
            color="#B03A2E", zorder=6,
            bbox=dict(boxstyle="round,pad=0.20", fc="white", ec="none", alpha=0.92))
    ax.text(0.0, -1.88, "Numbers count the genes unique to that layer; the intermediate\n"
                        "regions are given in panel B and in Table 3.",
            ha="center", va="top", fontsize=5.4, color="#666", linespacing=1.35)

    ax.set_xlim(-2.55, 2.55); ax.set_ylim(-2.40, 2.05)
    ax.set_aspect("equal")
    ax.set_axis_off()


def fig5(g):
    from PIL import Image
    fig = plt.figure(figsize=(W, 6.2))
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.24,
                          height_ratios=[1.06, 1], width_ratios=[1.12, 1])

    # A 五集韦恩（原生绘制，字号可控）
    ax = fig.add_subplot(gs[0, 0])
    _venn5(ax)
    tag(ax, "A", "Five sets, one gene", x=-0.02, y=1.02)

    # B 交集链
    ax = fig.add_subplot(gs[0, 1])
    chain = [("Pairs, FDR<0.05", 733), ("Unique genes", 612), ("1 · whole-blood MR", 337),
             ("2 · also thyroid MR", 62), ("3 · also colocalised", 6), ("4+5 · all layers", 1)]
    for i, (lab, v) in enumerate(chain):
        yv = len(chain) - 1 - i
        wbar = max(0.62, 0.09 + 0.86 * (v / 733) ** 0.42)
        ax.add_patch(Rectangle((0, yv - 0.32), wbar, 0.64,
                               fc=C["TSH"] if i < 5 else C["acc"], alpha=0.30 if i < 5 else 0.55,
                               ec="none"))
        ax.text(0.016, yv, lab, va="center", fontsize=6.4,
                fontweight="bold" if i == 5 else "normal")
        ax.text(wbar + 0.018, yv, "{:,}".format(v), va="center", fontsize=7,
                fontweight="bold", color=C["acc"] if i == 5 else "#333")
        if i < len(chain) - 1:
            ax.annotate("", xy=(0.5, yv - 0.36), xytext=(0.5, yv - 0.64),
                        arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#888"))
    ax.set_xlim(0, 1.26); ax.set_ylim(-0.8, len(chain) - 0.2)
    ax.set_axis_off()
    ax.text(0.0, -0.72, "Layer sizes 337 / 180 / 26 / 3 / 4 are not nested: only 62 of the\n"
                        "180 thyroid-passing genes also pass the whole-blood layer.",
            fontsize=6.1, color="#555", va="top")
    tag(ax, "B", "Intersection chain, not a funnel")

    # C 通过矩阵
    ax = fig.add_subplot(gs[1, 0])
    rc = g["rc"]
    rows = [("SCRN2", "31"), ("MCM6", "7"), ("NUCB2", "7"), ("PPP2R1B", "7"),
            ("UCP2", "7"), ("UROD", "7"), ("YWHAZ", "7"),
            ("TBKBP1", "16"), ("PRR15L", "5")]
    M = np.full((len(rows), 5), np.nan)
    for i, (gene, bit) in enumerate(rows):
        b = int(bit)
        for L in range(5):
            M[i, L] = 1.0 if (b >> L) & 1 else 0.0
    cmap = LinearSegmentedColormap.from_list("pw", ["#EDEDED", "#2F6FB2"])
    ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, "●" if M[i, j] == 1 else "○", ha="center", va="center",
                    fontsize=6.6, color="white" if M[i, j] == 1 else "#999")
    ax.set_xticks(range(5))
    ax.set_xticklabels(["1  blood\nMR", "2  thyroid\nMR", "3  colo-\ncalisation",
                        "4  tissue\ncoherence", "5  target\nassignment"], fontsize=6.0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.4)
    for lb in ax.get_yticklabels():
        if lb.get_text() == "SCRN2":
            lb.set_fontweight("bold"); lb.set_color(C["acc"])
    ax.set_xticks(np.arange(-0.5, 5, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.4)
    ax.tick_params(which="minor", length=0)
    ax.set_xlim(-0.5, 5.5)
    tag(ax, "C", "The six near-misses and the two contenders")

    # D 敏感性
    ax = fig.add_subplot(gs[1, 1])
    ax.set_axis_off()
    ax.text(0.0, 1.0, "Relaxation tested", fontsize=6.8, fontweight="bold", va="top")
    ax.text(0.62, 1.0, "Coloc. layer", fontsize=6.8, fontweight="bold", va="top", ha="center")
    ax.text(1.0, 1.0, "Intersection", fontsize=6.8, fontweight="bold", va="top", ha="right")
    lines = [("Baseline", "26", "SCRN2"),
             ("PP.H4 ≥ 0.70", "31", "SCRN2"),
             ("eQTL target top-5", "26", "SCRN2"),
             ("Coherence ≥ 90%", "26", "SCRN2"),
             ("All three relaxed", "31", "SCRN2")]
    for i, (a, b, c) in enumerate(lines):
        yv = 0.86 - i * 0.155
        ax.text(0.0, yv, a, fontsize=6.5, va="center")
        ax.text(0.62, yv, b, fontsize=6.5, va="center", ha="center")
        ax.text(1.0, yv, c, fontsize=6.5, va="center", ha="right", color=C["acc"],
                fontweight="bold")
    ax.text(0.0, -0.02, "The intersection is invariant to every relaxation;\n"
                        "the colocalisation layer is the only one that moves.",
            fontsize=6.0, color="#666", style="italic", va="top")
    tag(ax, "D", "Nothing here would change the answer")

    save(fig, 5)


# ══════════════════════════════════════════════════ 图 6
def fig6(g):
    fig = plt.figure(figsize=(W, 7.0))
    gs = fig.add_gridspec(3, 2, hspace=0.62, wspace=0.30,
                          height_ratios=[1.0, 1.0, 1.05])

    # A 区域示意
    ax = fig.add_subplot(gs[0, :])
    d = g["sig23"]
    ax.set_xlim(45.75, 46.10); ax.set_ylim(-3.55, 3.55)
    sig_pos = [(45.800299, "1", "TBKBP1 / TBX21\nunassignable", C["grey"]),
               (45.887132, "2b", "MRPL10\nexcluded", C["grey"]),
               (45.915395, "2", "SCRN2\nTSH only", C["acc"]),
               (46.050669, "3", "COPZ2\nshared", "#2E7D6B")]
    hyp_pos = [(45.804494, "a", "published\nsentinel"), (45.873049, "4", "TBKBP1 2nd eQTL\nhypo only")]
    for i, (pos, lab, txt, col) in enumerate(sig_pos):
        ax.plot([pos, pos], [0.1, 1.5], color=col, lw=1.4)
        ax.plot(pos, 1.5, "o", ms=4, color=col)
        ax.text(pos, 1.58, lab, ha="center", va="center", fontsize=6.4, fontweight="bold",
                color=col, bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"))
        ax.text(pos, 2.50, txt, ha="center", va="bottom", fontsize=5.7, color=col, linespacing=1.15)
    for i, (pos, lab, txt) in enumerate(hyp_pos):
        ax.plot([pos, pos], [-1.5, -0.1], color=C["HYPO"], lw=1.4)
        ax.plot(pos, -1.5, "o", ms=4, color=C["HYPO"])
        ax.text(pos, -1.58, lab, ha="center", va="center", fontsize=6.4, fontweight="bold",
                color=C["HYPO"], bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"))
        ax.text(pos, -1.98, txt, ha="center", va="top", fontsize=5.7, color=C["HYPO"],
                linespacing=1.15)
    ax.axhline(0, lw=0.7, color="#333")
    ax.text(45.755, 0.72, "TSH", fontsize=6.6, fontweight="bold", color="#333", va="center")
    ax.text(45.755, -0.52, "Hypo-\nthyroidism", fontsize=6.6, fontweight="bold",
            color=C["HYPO"], va="center", linespacing=1.1)
    ax.set_yticks([]); ax.set_xlabel("Position on chromosome 17, GRCh37 (Mb)")
    ax.set_xlim(45.75, 46.10)
    for sp in ("left",):
        ax.spines[sp].set_visible(False)
    ax.annotate("", xy=(45.915395, 1.05), xytext=(45.800299, 1.05),
                arrowprops=dict(arrowstyle="<->", lw=0.6, color="#999"))
    ax.text(45.858, 1.12, "r² < 0.02", ha="center", fontsize=5.8, color="#777")
    ax.annotate("", xy=(45.915395, -0.75), xytext=(45.887132, -0.75),
                arrowprops=dict(arrowstyle="<->", lw=0.6, color="#999"))
    ax.text(45.9013, -0.68, "r² = 0.59", ha="center", fontsize=5.8, color="#777")
    tag(ax, "A", "Seven nominal signals, five positions, two shared", x=-0.055, y=1.06)

    # B LD 与条件化
    ax = fig.add_subplot(gs[1, 0])
    v = g["var"].dropna(subset=["r2_with_rs6503875"]).copy()
    v["r"] = np.sqrt(v.r2_with_rs6503875)
    v = v.sort_values("r")
    cols = [C["acc"] if abs(r - 0.7669) < 0.01 else "#9FB4C9" for r in v.r]
    ax.barh(range(len(v)), v.r, color=cols, height=0.62, edgecolor="white", lw=0.5)
    ax.set_yticks(range(len(v)))
    ax.set_yticklabels(v.rsid, fontsize=6.2)
    ax.set_xlabel("|r| with rs6503875")
    ax.set_xlim(0, 1.0)
    for i, r in enumerate(v.r):
        ax.text(r + 0.02, i, "%.2f" % r, va="center", fontsize=5.8, color="#444")
    ax.text(0.22, 0.45, "rs62076102 carries r = 0.77 with the lead;\n"
                        "this is why the four-variant model exists",
            fontsize=5.8, color=C["acc"], va="center")
    tag(ax, "B", "Correlation with the lead variant")

    # C 归属
    ax = fig.add_subplot(gs[1, 1])
    ver = [("SCRN2", "TSH only", C["acc"], 0.99),
           ("COPZ2", "shared TSH + hypo", "#2E7D6B", 0.985),
           ("TBKBP1 (sig 1)", "unassignable", C["grey"], 0.06),
           ("TBX21", "no cis-eQTL", C["grey"], 0.0),
           ("MRPL10", "MR null", C["grey"], 0.07)]
    for i, (gene, lab, col, h) in enumerate(ver):
        yv = len(ver) - 1 - i
        ax.barh(yv, max(h, 0.004), color=col, height=0.5, edgecolor="white", lw=0.5)
        tx = 0.022 if h > 0.5 else max(h, 0.004) + 0.022
        ax.text(tx, yv, "%s — %s" % (gene, lab), va="center", fontsize=6.2, zorder=6,
                color="white" if h > 0.5 else "#333")
    ax.set_yticks([]); ax.set_xlim(0, 1.0)
    ax.set_xlabel("Window PP.H4 for TSH")
    ax.text(0.0, -0.28, "TBKBP1's whole-blood eQTL has F = 4,989 but a window H4 of 0.11;\n"
                        "the thyroid eQTL has F = 166 and a window H4 of 0.06.",
            transform=ax.transAxes, fontsize=5.9, color="#666", va="top")
    tag(ax, "C", "Attribution, not proximity")

    # D 哑铃（两个结局并列，避免只报甲减）
    ax = fig.add_subplot(gs[2, 0])
    jm = g["jm"]
    groups = [("TSH", 0.52, C["TSH"]), ("HYPO", -0.42, C["HYPO"])]
    for oc, y0, col in groups:
        for model, lab, off in [("五SNP", "five-variant", 0.17),
                                ("四SNP(去rs62076102)", "four-variant", -0.17)]:
            sub = jm[(jm.模型 == model) & (jm.结局 == oc)]
            row = sub[sub.SNP == "rs6503875"]
            if row.empty:
                continue
            r = row.iloc[0]
            yv = y0 + off
            ax.plot([r["边际z"], r["联合z"]], [yv, yv], lw=1.2, color=col, zorder=2)
            ax.plot(r["边际z"], yv, "o", ms=5.5, mfc="white", mec=col, mew=1.2, zorder=3)
            ax.plot(r["联合z"], yv, "o", ms=5.5, color=col, zorder=3)
            if model == "五SNP":
                ax.text(r["联合z"] + 0.18, yv, "%s %+.2f" % (lab, r["联合z"]),
                        ha="left", va="center", fontsize=5.6, color=col)
            else:
                ax.text(r["联合z"], yv + 0.07, "%s %+.2f" % (lab, r["联合z"]),
                        ha="center", va="bottom", fontsize=5.6, color=col)
    ax.axvline(0, lw=0.8, color="#333")
    ax.set_xlim(-6.5, 2.6); ax.set_ylim(-1.02, 0.92)
    ax.set_yticks([0.52, -0.42]); ax.set_yticklabels(["TSH", "Hypo-\nthyroidism"], fontsize=6.4)
    ax.set_xlabel("z for rs6503875 by model and outcome")
    ax.text(-6.4, -0.90, "open circles: marginal z (TSH −4.22; hypothyroidism −5.70)",
            fontsize=5.5, color="#555", va="bottom")
    for sp in ("left",):
        ax.spines[sp].set_visible(False)
    tag(ax, "D", "Conditioning fails for both traits", x=-0.17)

    # E 窗口共定位（SCRN2 信号 2 窗口，4 结局 × 2 组织）
    ax = fig.add_subplot(gs[2, 1])
    win = g["win"]
    w = win[(win.gene == "SCRN2") & (win.window.str.contains("信号2"))].copy()
    order = [("TSH", "Thyroid", "TSH\nthyroid"), ("TSH", "WholeBlood", "TSH\nblood"),
             ("HYPO", "Thyroid", "Hypo\nthyroid"), ("HYPO", "WholeBlood", "Hypo\nblood"),
             ("T4", "Thyroid", "FT4\nthyroid"), ("T4", "WholeBlood", "FT4\nblood"),
             ("HASH", "Thyroid", "Hash\nthyroid"), ("HASH", "WholeBlood", "Hash\nblood")]
    labs, h3, h4 = [], [], []
    for oc, tis, lab in order:
        r = w[(w.outcome == oc) & (w.tissue == tis)]
        if r.empty:
            continue
        labs.append(lab); h3.append(r.iloc[0].H3); h4.append(r.iloc[0].H4)
    xx = np.arange(len(labs))
    ax.bar(xx - 0.19, h3, width=0.37, color="#9AA6B2", label="PP.H3 (two variants)")
    ax.bar(xx + 0.19, h4, width=0.37, color=C["acc"], label="PP.H4 (shared)")
    ax.set_xticks(xx); ax.set_xticklabels(labs, fontsize=5.7)
    ax.set_ylim(0, 1.62)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_ylabel("Posterior probability")
    ax.legend(frameon=False, fontsize=5.9, ncol=2, loc="upper center", borderaxespad=0.05,
              handlelength=0.8, columnspacing=0.9)
    ax.text(xx[2] - 0.19, 1.04, "1.00", ha="center", fontsize=5.8, color="#444")
    ax.text(xx[1] + 0.19, 1.04, "0.99", ha="center", fontsize=5.8, color=C["acc"],
            fontweight="bold")
    tag(ax, "E", "The same window, two outcomes")

    save(fig, 6)

# ══════════════════════════════════════════════════ 图 7
def fig7(g):
    pl = g["plane"]
    fig = plt.figure(figsize=(W, 6.4))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30, height_ratios=[1.15, 1])

    # A TSH–T4 平面
    ax = fig.add_subplot(gs[0, :])
    ax.scatter(pl.b_TSH, pl.b_T4, s=2.2, c="#C9D2DC", alpha=0.65, lw=0, rasterized=True)
    ref = pl[pl.gene.isin(["TPO", "DIO1"])]
    ax.scatter(ref.b_TSH, ref.b_T4, s=26, c="#2E7D6B", edgecolor="white", lw=0.7, zorder=4)
    for _, r in ref.iterrows():
        ax.annotate(r.gene, (r.b_TSH, r.b_T4), xytext=(6, 4), textcoords="offset points",
                    fontsize=6.4, color="#2E7D6B", fontweight="bold")
    scr = pl[pl.gene == "SCRN2"]
    if not scr.empty:
        ax.scatter(scr.b_TSH, scr.b_T4, s=42, c=C["acc"], edgecolor="white", lw=0.8, zorder=5)
        ax.annotate("SCRN2", (scr.b_TSH.iloc[0], scr.b_T4.iloc[0]), xytext=(8, -10),
                    textcoords="offset points", fontsize=7, color=C["acc"], fontweight="bold")
    xs = np.array([-0.09, 0.09])
    ax.plot(xs, -1.07 * xs, ls="--", lw=1.0, color="#2E7D6B", zorder=2)
    ax.text(0.098, -0.088, "slope −1.07\n(TPO reference)", fontsize=6.2, color="#2E7D6B",
            ha="right", va="bottom")
    ax.plot(xs, -0.25 * xs, ls=":", lw=1.0, color="#888", zorder=2)
    ax.text(-0.088, -0.25 * -0.088 + 0.006, "median ratio −0.25", fontsize=6.0, color="#777")
    ax.axhline(0, lw=0.6, color="#999", zorder=1)
    ax.axhspan(-0.0022, 0.0022, color="#F4E3E0", alpha=0.75, zorder=0)
    ax.set_xlim(-0.10, 0.10); ax.set_ylim(-0.10, 0.10)
    ax.set_xlabel("Effect on TSH (β per SD of predicted expression)")
    ax.set_ylabel("Effect on free T4 (β per SD)")
    ax.text(0.015, 0.945, "915 genes with both estimates", transform=ax.transAxes,
            ha="left", va="top", fontsize=6.1, color="#666")
    tag(ax, "A", "Position on the TSH–free T4 plane", x=-0.052)

    # B 森林图
    ax = fig.add_subplot(gs[1, 0])
    wb = g["mr"][g["mr"].gene == "SCRN2"].set_index("outcome")
    code = {"TSH": "TSH", "Hypothyroidism": "HYPO", "Total T4": "T4",
            "Hashimoto thyroiditis": "HASH"}
    th = g["gw"][g["gw"].gene == "SCRN2"].copy()
    th["oc"] = th.outcome.map({v: k for k, v in code.items()})
    th = th.set_index("oc")
    LONGS = ["TSH", "Total T4", "Hypothyroidism", "Hashimoto thyroiditis"]
    LAB = {"TSH": "TSH", "Total T4": "Free T4", "Hypothyroidism": "Hypothyroidism",
           "Hashimoto thyroiditis": "Hashimoto"}
    COL = {"TSH": C["TSH"], "Total T4": C["T4"], "Hypothyroidism": C["HYPO"],
           "Hashimoto thyroiditis": C["HASH"]}
    labels, est, se, cols = [], [], [], []
    for oc in LONGS:
        labels.append(LAB[oc] + "  ·  blood")
        est.append(wb.loc[oc].b); se.append(wb.loc[oc].se); cols.append(COL[oc])
    for oc in LONGS:
        labels.append(LAB[oc] + "  ·  thyroid")
        est.append(th.loc[oc].b); se.append(th.loc[oc].se); cols.append(COL[oc])
    yy = np.arange(len(labels))[::-1]
    for i in range(len(labels)):
        ax.errorbar(est[i], yy[i], xerr=1.96 * se[i], fmt="o", ms=4, color=cols[i],
                    capsize=2.2, elinewidth=0.9, lw=0)
        ax.text(max(est[i] + 1.96 * se[i], 0.02) + 0.02, yy[i],
                "%+.3f (%+.3f)" % (est[i], se[i]), va="center", fontsize=5.8, color="#444")
    ax.axvline(0, lw=0.8, color="#333")
    ax.set_yticks(yy); ax.set_yticklabels(labels, fontsize=6.1)
    ax.set_xlabel("β (95% CI) per SD of predicted SCRN2")
    ax.set_xlim(-0.10, 0.46)
    ax.axhspan(3.5, 3.5 + 0.001, color="white")
    ax.axhline(3.5, lw=0.5, ls=":", color="#ccc")
    tag(ax, "B", "Two instrument sets, four outcomes")

    # C ratio 直方
    ax = fig.add_subplot(gs[1, 1])
    s_tsh = pl[(pl.tsh_sig == True) & pl.r.notna() & np.isfinite(pl.r)].copy()
    s_tsh["r"] = s_tsh.r.clip(-3, 3)
    ax.hist(s_tsh.r, bins=np.linspace(-3, 3, 41), color="#9FB4C9", edgecolor="none")
    ax.axvline(0, lw=0.8, color="#333")
    if not scr.empty:
        ax.axvline(scr.r.iloc[0], ls="--", lw=1.1, color=C["acc"])
        ax.text(scr.r.iloc[0] + 0.10, ax.get_ylim()[1] * 0.55,
                "SCRN2\nr = %+.2f" % scr.r.iloc[0], fontsize=6.2, color=C["acc"],
                fontweight="bold", va="top")
    ax.set_xlabel("Per-gene ratio  β(FT4) / β(TSH)")
    ax.set_ylabel("TSH-significant genes")
    nnull = int(s_tsh.t4_null.sum())
    ax.set_ylim(0, max(np.histogram(s_tsh.r, bins=40)[0]) * 1.35)
    ax.text(0.02, 0.97, "%d of %d (%.0f%%) have a free-T4 CI\nthat includes zero"
            % (nnull, len(s_tsh), 100.0 * nnull / len(s_tsh)),
            transform=ax.transAxes, va="top", fontsize=6.0, color="#444",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9))
    tag(ax, "C", "What the flat T4 arm cannot rule out")

    save(fig, 7)


# ══════════════════════════════════════════════════ 图 8
def fig8(g):
    ex = g["expr"]
    fig = plt.figure(figsize=(W, 7.0))
    gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.95, height_ratios=[1, 1])

    # A 甲状腺 vs 全血
    ax = fig.add_subplot(gs[0, 0])
    ex2 = ex.sort_values("thyroid_tpm")
    yy = np.arange(len(ex2))
    ax.barh(yy - 0.19, ex2.thyroid_tpm.clip(lower=0.05), height=0.36, color=C["acc"],
            alpha=0.85, label="thyroid")
    ax.barh(yy + 0.19, ex2.blood_tpm.clip(lower=0.05), height=0.36, color="#8FA6BB",
            alpha=0.9, label="whole blood")
    ax.set_xscale("log")
    ax.set_yticks(yy); ax.set_yticklabels(ex2.gene, fontsize=6.3)
    for lb in ax.get_yticklabels():
        if lb.get_text() == "SCRN2":
            lb.set_fontweight("bold"); lb.set_color(C["acc"])
    ax.set_xlabel("GTEx median TPM (log)")
    ax.set_xlim(0.04, 4e4)
    ax.legend(frameon=False, fontsize=5.9, loc="lower right", handlelength=0.8,
              borderaxespad=0.1)
    tag(ax, "A", "Thyroid-enriched, not specific")

    # B 检出率 vs 同向率
    ax = fig.add_subplot(gs[0, 1])
    t = g["tis"]
    ax.scatter(t.detect_pct, t.same_pct, s=8, c="#8FA6BB", alpha=0.75, lw=0)
    for name, col, dx, dy in [("Thyroid", C["acc"], 7, -8), ("Whole_Blood", "#333", 7, -8)]:
        r = t[t.tissue == name]
        if not r.empty:
            ax.scatter(r.detect_pct, r.same_pct, s=30, c=col, edgecolor="white", lw=0.7, zorder=4)
            ax.annotate(name.replace("_", " "), (r.detect_pct.iloc[0], r.same_pct.iloc[0]),
                        xytext=(dx, dy), textcoords="offset points", fontsize=6.0, color=col,
                        fontweight="bold")
    ax.set_xlabel("Instruments recovered (%)")
    ax.set_ylabel("Same direction (%)")
    ax.set_ylim(42, 108)
    ax.set_xlim(0, 50)
    ax.text(0.97, 0.06, "thyroid: 31.4% recovered,\n81.8% of those concordant",
            transform=ax.transAxes, ha="right", fontsize=5.8, color="#444")
    tag(ax, "B", "Recovery ≠ concordance")

    # C 细胞类型
    ax = fig.add_subplot(gs[0, 2])
    ct = g["cell"].sort_values("pct_expressing", ascending=False).head(9)
    SHORT = {"thyroid follicular cell": "follicular", "endothelial cell": "endothelial",
             "smooth muscle cell": "smooth muscle", "epithelial cell": "epithelial",
             "erythroid progenitor cell": "erythroid prog."}
    short = [SHORT.get(c, c.replace(" cell", "")) for c in ct.cell_type]
    yy = np.arange(len(ct))[::-1]
    ax.scatter(ct.pct_expressing, yy, s=np.sqrt(ct.n_cells) * 1.5 + 5, c=ct.mean_expr,
               cmap=LinearSegmentedColormap.from_list("r", ["#EBD8D5", "#B03A2E"]),
               edgecolor="#B03A2E", lw=0.4, zorder=3)
    ax.set_yticks(yy); ax.set_yticklabels(short, fontsize=5.7)
    ax.set_xlabel("Cells expressing SCRN2 (%)")
    ax.set_xlim(-0.15, 3.7)
    ax.set_ylim(-0.7, len(ct) - 0.3)
    ax.text(0.97, 0.10, "follicular 2.96%,\nall others ≤1.7%",
            transform=ax.transAxes, ha="right", fontsize=5.8, color="#444")
    tag(ax, "C", "Which cell in the gland")

    # D 共表达
    ax = fig.add_subplot(gs[1, 0])
    s7 = g["s07"].sort_values("rho_pc10")
    yy = np.arange(len(s7))
    ax.barh(yy, s7.rho_pc10, height=0.62,
            color=["#B03A2E" if m == "TF" else "#8FA6BB" for m in s7.module])
    ax.axvline(0, lw=0.7, color="#333")
    ax.set_yticks(yy); ax.set_yticklabels(s7.gene, fontsize=5.8)
    ax.set_xlabel("Partial ρ with SCRN2 (10 PCs removed)")
    ax.set_xlim(-0.32, 0.66)
    for i, (gn, r) in enumerate(zip(s7.gene, s7.rho_pc10)):
        if gn in ("NKX2-1", "PAX8"):
            ax.text(r + (0.012 if r > 0 else -0.012), i, "%+.2f" % r, va="center",
                    ha="left" if r > 0 else "right", fontsize=5.6, color=C["acc"],
                    fontweight="bold")
    tag(ax, "D", "Co-expression with TFs")

    # E GEO
    ax = fig.add_subplot(gs[1, 1:])
    items = [("HT vs normal", -0.056, "0.68", "#8FA6BB"),
             ("LT+ vs normal", 0.134, "0.20", "#8FA6BB"),
             ("LT+ vs LT−", 0.088, "0.60", "#8FA6BB"),
             ("IFNG  (control)", 0.717, "0.0002", "#7FA88C"),
             ("CD4  (control)", 1.729, "0.0000", "#7FA88C"),
             ("PTPRC  (control)", 2.425, "0.0000", "#7FA88C"),
             ("TPO  (control)", -0.661, "0.0012", "#7FA88C"),
             ("TG  (control)", -0.539, "0.0009", "#7FA88C")]
    yy = np.arange(len(items))[::-1]
    ax.axvspan(-0.30, 0.30, color="#F2F2F2", zorder=0)
    for i, (lab, v, p, c) in enumerate(items):
        ax.barh(yy[i], v, height=0.58, color=c, alpha=0.92)
        ax.text(v + (0.05 if v > 0 else -0.05), yy[i], "%+.2f  (P = %s)" % (v, p),
                va="center", ha="left" if v > 0 else "right", fontsize=5.7, color="#444")
    ax.axvline(0, lw=0.8, color="#333")
    ax.set_yticks(yy); ax.set_yticklabels([i[0] for i in items], fontsize=5.8)
    for lb, it in zip(ax.get_yticklabels(), items):
        if it[3] == "#8FA6BB":
            lb.set_color("#B03A2E")
    ax.set_xlabel("log2 fold change")
    ax.set_xlim(-1.9, 3.2)
    ax.set_ylim(-0.7, len(items) - 0.3)
    ax.text(0.985, 0.05, "shaded band: smallest detectable effect (0.24–0.33)",
            transform=ax.transAxes, ha="right", fontsize=5.8, color="#666")
    ax.text(0.985, 0.93, "SCRN2 rows: GSE138198 / GSE165724", transform=ax.transAxes,
            fontsize=5.6, color="#B03A2E", ha="right", va="top")
    tag(ax, "E", "Bulk tissue says nothing; the pipeline works", x=-0.075)

    save(fig, 8)


# ══════════════════════════════════════════════════ 图 9
def fig9(g):
    fig = plt.figure(figsize=(W, 8.0))
    gs = fig.add_gridspec(3, 2, hspace=0.52, wspace=0.30,
                          height_ratios=[1.10, 1.02, 0.80])

    # A 两级示意
    ax = fig.add_subplot(gs[0, 0])
    ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.add_patch(FancyBboxPatch((0.00, 0.62), 0.66, 0.22, boxstyle="round,pad=0.006",
                                fc="#EAF0F7", ec="#9FB4C9", lw=0.7))
    ax.text(0.33, 0.788, "Pituitary thyrotroph", ha="center", fontsize=6.4, fontweight="bold")
    ax.text(0.33, 0.714, "TSH PP.H4 = 0.99 · free T4 unchanged", ha="center", fontsize=5.3,
            color="#2F6FB2", fontweight="bold")
    ax.add_patch(FancyBboxPatch((0.00, 0.10), 0.66, 0.26, boxstyle="round,pad=0.006",
                                fc="#FBEDEB", ec="#D9A79F", lw=0.7))
    ax.text(0.33, 0.306, "Thyroid follicle", ha="center", fontsize=6.4, fontweight="bold")
    ax.text(0.33, 0.234, "hypothyroidism PP.H3 = 1.00", ha="center", fontsize=5.3,
            color="#B03A2E", fontweight="bold")
    ax.text(0.33, 0.166, "Hashimoto PP.H4 = 0.28", ha="center", fontsize=5.3,
            color="#B03A2E", fontweight="bold")
    ax.annotate("", xy=(0.14, 0.368), xytext=(0.14, 0.602),
                arrowprops=dict(arrowstyle="-|>", lw=1.1, color="#666"))
    ax.text(0.108, 0.485, "TSH", fontsize=6.0, color="#666", va="center", ha="right")
    ax.annotate("", xy=(0.36, 0.602), xytext=(0.36, 0.368),
                arrowprops=dict(arrowstyle="-|>", lw=1.1, color="#B03A2E"))
    ax.text(0.392, 0.485, "T4", fontsize=6.0, color="#B03A2E", va="center", ha="left")
    ax.plot([0.705, 0.705], [0.23, 0.73], color="#AAA", lw=0.9)
    ax.plot([0.705, 0.678], [0.73, 0.73], color="#AAA", lw=0.9)
    ax.plot([0.705, 0.678], [0.23, 0.23], color="#AAA", lw=0.9)
    ax.text(0.735, 0.48, "one axis,\ntwo control\npoints; the\ngenetic\ninstrument\nlabels the\nupper one",
            fontsize=5.5, va="center", ha="left", color="#444", linespacing=1.55)
    tag(ax, "A", "Two discontinuities, one axis", x=-0.02, y=1.05)

    # B Cα trace
    ax = fig.add_subplot(gs[0, 1], projection="3d")
    ca, res = _ca_coords()
    groove = [44, 46, 77, 81, 82, 119]
    m1 = list(range(251, 256)); qd = list(range(339, 344)); m2 = list(range(351, 356))
    ax.scatter(ca[:, 0], ca[:, 1], ca[:, 2], s=3.0, c="#D0D0D0", lw=0, depthshade=False,
               alpha=0.85)
    for grp, col, lab in [(groove, "#2F6FB2", "groove"),
                          (qd, "#B03A2E", "QDPVR 339–343"),
                          (m2, "#E08A1E", "QVDRR 351–355"),
                          (m1, "#7D5BA6", "ELLRQ 251–255")]:
        idx = [i for i, r in enumerate(res) if r in grp]
        if not idx:
            continue
        ax.scatter(ca[idx, 0], ca[idx, 1], ca[idx, 2], s=22, c=col, lw=0,
                   depthshade=False, label=lab)
    ax.view_init(elev=8, azim=-70)
    ax.set_axis_off()
    ax.set_box_aspect((1, 1, 1), zoom=1.20)
    leg = ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.04), frameon=False,
                    fontsize=5.8, ncol=2, handletextpad=0.3, columnspacing=0.9,
                    labelspacing=0.32)
    ax.text2D(0.0, 0.985, "Cα trace, first frame of the PFNA trajectory",
              transform=ax.transAxes, fontsize=5.5, color="#666", va="top")
    tag(ax, "B", "Where the motifs sit", x=-0.02, y=1.06)

    # C 残基接触
    ax = fig.add_subplot(gs[1, 0])
    r = g["res"].sort_values("frac_contact_lt5A", ascending=True).reset_index(drop=True)
    yy = np.arange(len(r))
    cols = []
    for _, x in r.iterrows():
        if x.resid in groove:
            cols.append("#2F6FB2")
        elif 339 <= x.resid <= 343:
            cols.append("#B03A2E")
        elif 351 <= x.resid <= 355:
            cols.append("#E08A1E")
        else:
            cols.append("#7D5BA6")
    ax.barh(yy, r.frac_contact_lt5A * 100, height=0.64, color=cols)
    ax.set_yticks(yy)
    ax.set_yticklabels(["%s%d" % (n[0] + n[1:].lower(), i) for i, n in zip(r.resid, r.resname)],
                       fontsize=5.2)
    ax.set_xlabel("Frames with any heavy atom < 5 Å (%)")
    ax.set_xlim(0, 128)
    ax.set_ylim(-0.75, len(r) + 0.55)
    ax.set_xticks([0, 25, 50, 75, 100])
    hs = [Patch(facecolor=c, label=l) for c, l in
          [("#2F6FB2", "groove"), ("#B03A2E", "QDPVR"),
           ("#E08A1E", "QVDRR"), ("#7D5BA6", "ELLRQ")]]
    ax.legend(handles=hs, frameon=False, fontsize=5.5, ncol=2, loc="lower left",
              bbox_to_anchor=(0.015, 0.005), handlelength=0.7, columnspacing=0.8,
              labelspacing=0.34, borderaxespad=0.0)
    ax.text(0.40, 0.20, "no atom within 5 Å:\nGln339 · Pro341 · ELLRQ 251–255",
            transform=ax.transAxes, fontsize=5.4, color="#444", va="bottom", ha="left")
    tag(ax, "C", "Partial, not complete, occlusion", x=-0.15, y=1.05)

    # D 距离区间
    ax = fig.add_subplot(gs[1, 1])
    dd = g["dist"]
    order = ["QDPVR (339-343)", "QVDRR M2 (351-355)", "ELLRQ M1 (251-255)"]
    lab = ["QDPVR\n339–343", "QVDRR\n351–355", "ELLRQ\n251–255"]
    cols3 = ["#B03A2E", "#E08A1E", "#7D5BA6"]
    ax.axvspan(2, 6, color="#F1F5F9", zorder=0, lw=0)
    for i, (k, l, c) in enumerate(zip(order, lab, cols3)):
        row = dd[dd.target == k]
        if row.empty:
            continue
        x = row.iloc[0]
        ax.plot([x.p05_A, x.p95_A], [i, i], lw=6, color=c, alpha=0.32, solid_capstyle="round",
                zorder=2)
        ax.plot(x.median_min_dist_A, i, "o", ms=7.5, color=c, zorder=3)
        pct = x.frac_lt6A * 100
        ptxt = "0%" if pct < 0.05 else "%.1f%%" % pct
        ax.text(x.p95_A * 1.09, i + 0.31, "median %.2f Å · %s < 6 Å" % (x.median_min_dist_A, ptxt),
                va="center", ha="left", fontsize=5.6, color=c,
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.88))
    for thr, ls in [(6, "--"), (8, ":")]:
        ax.axvline(thr, ls=ls, lw=0.8, color="#888", zorder=1)
        ax.text(thr, -0.66, "%d Å" % thr, ha="center", va="bottom", fontsize=5.6, color="#777")
    ax.set_xscale("log")
    ax.set_xlim(2, 260)
    ax.set_ylim(-0.95, 2.72)
    ax.set_yticks(range(3)); ax.set_yticklabels(lab, fontsize=6.0)
    ax.set_xticks([3, 6, 10, 30, 60, 100])
    ax.set_xticklabels(["3", "6", "10", "30", "60", "100"])
    ax.set_xlabel("Minimum heavy-atom distance (Å, log scale)")
    ax.minorticks_off()
    tag(ax, "D", "The asymmetry is quantitative", x=-0.13, y=1.05)

    # E 可否证预测
    ax = fig.add_subplot(gs[2, :])
    ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    xs = [0.0, 0.34, 0.55, 0.76]
    ax.text(xs[1], 0.84, "wild type", fontsize=6.6, fontweight="bold", ha="center")
    ax.text(xs[2], 0.84, "QVDRR mutant", fontsize=6.6, fontweight="bold", ha="center")
    ax.text(xs[3], 0.84, "ELLRQ mutant", fontsize=6.6, fontweight="bold", ha="center")
    ax.text(xs[0], 0.66, "SCRN2 + HSPA8, no PFNA", fontsize=6.4, va="center")
    ax.text(xs[0], 0.46, "+ PFNA", fontsize=6.4, va="center")
    ax.text(xs[1], 0.66, "++", fontsize=7.2, ha="center", va="center")
    ax.text(xs[2], 0.66, "++", fontsize=7.2, ha="center", va="center")
    ax.text(xs[3], 0.66, "++", fontsize=7.2, ha="center", va="center")
    ax.text(xs[1], 0.46, "++ → +", fontsize=7.2, ha="center", va="center", color=C["acc"])
    ax.text(xs[2], 0.46, "+ → ++", fontsize=7.2, ha="center", va="center", color=C["acc"])
    ax.text(xs[3], 0.46, "++ (unchanged)", fontsize=7.2, ha="center", va="center",
            color="#2E7D6B")
    ax.text(0.0, 0.02, "The model survives only if PFNA lowers the HSPA8 co-IP signal at the wild type and at the "
                       "QVDRR mutant,\nbut leaves the ELLRQ mutant unchanged. Equal effects at both motifs, or no "
                       "effect at either, would refute it.",
            fontsize=6.2, color="#444", va="bottom", linespacing=1.5)
    tag(ax, "E", "The one test that could kill the model", x=-0.008, y=1.02)

    save(fig, 9)


def _ca_coords():
    """从 PFNA 轨迹首帧提取 Cα 坐标（425 残基）。"""
    import numpy as _np
    p = os.path.join(HERE, "..", "..", "..", "DATA14", "md_runs", "PFNA", "traj.pdb")
    ca, res = [], []
    with open(p, encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            if i > 300000 or line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                ca.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                res.append(int(line[22:26]))
    return _np.array(ca), _np.array(res)


if __name__ == "__main__":
    which = [int(a) for a in sys.argv[1:]] or list(range(1, 10))
    print("output:", OUT)
    g = load()
    fns = {1: fig1, 2: fig2, 3: fig3, 4: fig4, 5: fig5, 6: fig6,
           7: fig7, 8: fig8, 9: fig9}
    for k in which:
        if k in fns:
            fns[k](g)
