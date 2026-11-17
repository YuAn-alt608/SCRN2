#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
38_scrn2_ctd_edc.py
回到 CTD 数据库，挖掘 SCRN2 的 EDC 关联证据

背景：本研究"EDC 响应基因 × 甲状腺 eQTL-MR"框架中，SCRN2 是主推候选。
      但收尾阶段 SCRN2 的叙事越来越偏"甲状腺遗传学"，EDC 这一层被稀释。
      本脚本回到 CTD chem-gene interactions，把 SCRN2 的 EDC 关联证据系统挖出，
      回答"这个研究还和 EDC 有什么关系"。

关键逻辑：
- SCRN2 作为 EDC 响应基因：CTD 中 98 条交互记录，其中 35 条属 EDC 类化学物
- EDC 暴露 -> SCRN2 表达下调（人类 mRNA 31/42 下调）-> 与 MR 负向效应自洽
  （SCRN2 低表达 -> 甲减/桥本风险增高）
- 因此完整因果链：EDC 暴露 -> SCRN2 表达下调 -> 甲状腺功能减退/桥本风险增高

输出：
- results/SCRN2_CTD_EDC_associations.csv  (EDC 类全部交互)
- results/SCRN2_CTD_EDC_human.csv          (聚焦人类证据)
- results/SCRN2_CTD_EDC_summary.txt        (文本报告)
- figures/Fig09_SCRN2_EDC_chain.png/pdf    (EDC->SCRN2->甲状腺 因果链图)
"""

import gzip
import os
import re
import sys
import numpy as np
import pandas as pd
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as fs

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(BASE, "raw")
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
TMP = os.path.join(BASE, "tmp", "geo_validation")

CTD_FILE = os.path.join(RAW, "CTD_chem_gene_ixns.tsv.gz")


# EDC 分类规则（与 01a_build_chemical_set.py 完全一致）
EDC_RULES = [
    ("PFAS", r"perfluor|polyfluor|\bPFOA\b|\bPFOS\b|\bPFNA\b|\bPFHxS\b|\bPFDA\b|\bPFBA\b|\bPFHpA\b|\bPFDoA\b"),
    ("Phthalate", r"phthalat"),
    ("Bisphenol", r"bisphenol"),
    ("PCB", r"polychlorinated biphenyl|\bPCB\s?\d|\bPCB\b"),
    ("PBDE", r"polybrominated diphenyl ether|\bPBDE\b|\btetrabromobisphenol\b|diphenyl ether"),
    ("Organochlorine", r"\bDDT\b|dichlorodiphenyltrichloroethane|\bDDE\b|\bDDD\b|hexachlorocyclohexane|lindane|chlordane|heptachlor|\baldrin\b|dieldrin|endosulfan|toxaphene|mirex|hexachlorobenzene"),
    ("Paraben", r"paraben"),
    ("Antimicrobial", r"triclosan|triclocarban"),
    ("Dioxin", r"\bdioxin\b|dibenzodioxin|dibenzofuran|tetrachlorodibenzodioxin|\bTCDD\b"),
    ("Metal", r"\bcadmium\b|\blead\b|\bmercury\b|\bmethylmercury\b|\barsenic\b|arsenite|arsenate|manganese|chromium|nickel|\bcopper\b"),
    ("Perchlorate", r"perchlorate|thiocyanate"),
    ("Alkylphenol", r"octylphenol|nonylphenol|alkylphenol|benzophenone"),
    ("Pesticide", r"glyphosate|atrazine|chlorpyrifos|malathion|permethrin|cypermethrin|imidacloprid|fipronil|diazinon|parathion|acetamiprid|fenthion|decamethrin|deltamethrin"),
]


def classify(name):
    n = name.lower()
    for cls, pat in EDC_RULES:
        if re.search(pat, n):
            return cls
    return None


def direction(actions):
    if "decreases^expression" in actions:
        return "down"
    if "increases^expression" in actions:
        return "up"
    if "affects" in actions:
        return "affects"
    return "other"


def load_scrN2_ctd():
    """提取 SCRN2 全部 CTD 交互记录"""
    rows = []
    with gzip.open(CTD_FILE, "rt", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 8 and p[3] == "SCRN2":
                rows.append({
                    "chemical": p[0], "mesh": p[1], "cas": p[2],
                    "gene": p[3], "gene_id": p[4],
                    "molecule": p[5], "organism": p[6], "taxid": p[7],
                    "interaction": p[8],
                    "actions": p[9] if len(p) > 9 else "",
                    "pmid": p[10] if len(p) > 10 else "",
                })
    return rows


def main():
    rows = load_scrN2_ctd()
    print(f"SCRN2 CTD 总交互记录: {len(rows)}")

    # 标注 EDC 类别与方向
    for r in rows:
        r["edc_class"] = classify(r["chemical"])
        r["direction"] = direction(r["actions"])
        r["is_human"] = r["organism"] == "Homo sapiens"

    edc_rows = [r for r in rows if r["edc_class"] is not None]
    print(f"EDC 类交互记录: {len(edc_rows)} / {len(rows)}")

    human_edc = [r for r in edc_rows if r["is_human"]]
    print(f"人类 EDC 交互记录: {len(human_edc)}")

    # ---------- 保存完整 EDC 关联表 ----------
    cols = ["chemical", "mesh", "cas", "edc_class", "molecule", "organism",
            "direction", "interaction", "pmid"]
    df_edc = pd.DataFrame([{k: r[k] for k in cols} for r in edc_rows])
    df_edc.to_csv(os.path.join(RES, "SCRN2_CTD_EDC_associations.csv"), index=False)

    # 人类证据子表
    df_human = pd.DataFrame([{k: r[k] for k in cols} for r in human_edc])
    df_human.to_csv(os.path.join(RES, "SCRN2_CTD_EDC_human.csv"), index=False)

    # ---------- 汇总统计 ----------
    # 人类 mRNA 方向分布（全部化学物，含非EDC，用于说明"EDC下调是普遍趋势"）
    human_mrna = [r for r in rows if r["is_human"] and r["molecule"] == "mRNA"]
    human_mrna_down = sum(1 for r in human_mrna if r["direction"] == "down")
    human_mrna_up = sum(1 for r in human_mrna if r["direction"] == "up")

    # 人类 EDC mRNA 方向
    human_edc_mrna = [r for r in human_edc if r["molecule"] == "mRNA"]
    human_edc_down = sum(1 for r in human_edc_mrna if r["direction"] == "down")
    human_edc_up = sum(1 for r in human_edc_mrna if r["direction"] == "up")

    # ---------- 文本报告 ----------
    L = []
    L.append("=" * 78)
    L.append("SCRN2 的 EDC 关联证据（回到 CTD 数据库）")
    L.append("=" * 78)
    L.append("")
    L.append("核心问题：本研究还和 EDC 有什么关系？")
    L.append("答案：SCRN2 正是以「EDC 响应基因」身份进入候选集的。CTD 证据确凿存在，")
    L.append("      但必须诚实报告一个方向性问题（见第五节）。")
    L.append("")
    L.append("-" * 78)
    L.append(f"一、SCRN2 在 CTD 中的交互证据总量")
    L.append(f"  总交互记录: {len(rows)} 条")
    L.append(f"  其中 EDC 类化学物: {len(edc_rows)} 条")
    L.append(f"  其中人类 EDC: {len(human_edc)} 条")
    L.append("")
    L.append("-" * 78)
    L.append("二、EDC 暴露对 SCRN2 表达的方向（人类证据）")
    L.append(f"  人类 mRNA 全部化学物: 下调 {human_mrna_down}/{len(human_mrna)}, "
             f"上调 {human_mrna_up}/{len(human_mrna)}")
    L.append(f"  人类 EDC mRNA: 下调 {human_edc_down}, 上调 {human_edc_up}")
    L.append("  → 绝大多数 EDC 暴露导致 SCRN2 表达【下调】")
    L.append("")
    L.append("-" * 78)
    L.append("三、最核心的 EDC 人类证据（直接、单一化学物、mRNA、下调）")
    key = [r for r in human_edc if r["molecule"] == "mRNA" and r["direction"] == "down"]
    for r in key:
        L.append(f"  [{r['edc_class']}] {r['chemical']} → SCRN2 mRNA ↓  "
                 f"(PMID {r['pmid']})")
    L.append("")
    L.append("-" * 78)
    L.append("四、PFAS 证据（本研究首要 EDC 类别，PMID 25812627 同一研究 3 种）")
    pf = [r for r in human_edc if r["edc_class"] == "PFAS"]
    for r in pf:
        L.append(f"  {r['chemical']} → {r['direction']} (PMID {r['pmid']})")
    L.append("")
    L.append("-" * 78)
    L.append("五、【必须诚实报告】方向性矛盾")
    L.append("  MR 主结果（Wald ratio）b 值全为【正】：")
    L.append("    SCRN2 表达↑ → TSH↑ (+0.027) / 甲减↑ (+0.070) / 桥本↑ (+0.176)")
    L.append("    即：SCRN2 高表达 = 疾病风险增高，低表达 = 保护")
    L.append("")
    L.append("  CTD 方向：EDC 暴露 → SCRN2 表达【下调】")
    L.append("")
    L.append("  两者结合 → EDC 下调 SCRN2 → 甲减风险【降低】")
    L.append("  这与「EDC 危害甲状腺」的公共卫生常识相矛盾。")
    L.append("")
    L.append("  三种诚实解读（论文中应如实讨论，而非强行闭环）：")
    L.append("   (1) CTD 多为高剂量/体外/异种的急性应激反应，")
    L.append("       不代表人类体内慢性低剂量暴露的净效应；")
    L.append("   (2) EDC 的甲状腺毒性【可能不通过】SCRN2 介导——")
    L.append("       SCRN2 是 EDC 响应基因，但其「响应」与「致病」方向解耦；")
    L.append("   (3) 甲状腺反馈调节使单向表达变化难以线性外推到疾病。")
    L.append("")
    L.append("  → 结论：SCRN2 的「EDC 响应基因」身份成立（CTD 证据充分），")
    L.append("      但其「介导 EDC 甲状腺毒性」的因果角色【尚不能】由现有证据支持，")
    L.append("      需在论文中降级为「候选机制线索」，而非确证因果链。")
    L.append("=" * 78)

    txt = "\n".join(L)
    with open(os.path.join(RES, "SCRN2_CTD_EDC_summary.txt"), "w", encoding="utf-8") as f:
        f.write(txt)
    print("\n" + txt)

    # ---------- 因果链图（如实反映方向矛盾）----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5.6))
    ax.axis("off")

    # 三节点：EDC → SCRN2 下调 → 但 MR 是"高表达致病"
    boxes = [
        (0.03, 0.55, "EDC exposure\nPFAS · Bisphenol A\nArsenic · PBDE",
         fs.C_ORANGE),
        (0.36, 0.55, "SCRN2 expression\nDOWN-regulated\n(CTD: 31/42 human mRNA ↓)",
         fs.C_RED),
        (0.69, 0.55, "Thyroid dysfunction\nMR: high SCRN2 →\nTSH↑ · Hypo↑ · Hashimoto↑",
         fs.C_PURPLE),
    ]
    for x, y, text, color in boxes:
        ax.add_patch(plt.Rectangle((x, y), 0.27, 0.30, facecolor=color,
                                   alpha=0.16, edgecolor=color, linewidth=1.8,
                                   transform=ax.transAxes))
        ax.text(x + 0.135, y + 0.15, text, ha="center", va="center",
                fontsize=9.5, transform=ax.transAxes, weight="bold")

    # 箭头1：EDC → SCRN2 下调（实线，方向确证）
    ax.annotate("", xy=(0.36, 0.70), xytext=(0.30, 0.70),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color=fs.C_DARK, lw=2.0))
    ax.text(0.33, 0.745, "down-regulates\n(consistent)", ha="center", fontsize=8.5,
            transform=ax.transAxes, style="italic", color=fs.C_DARK)

    # 箭头2：SCRN2 下调 → 疾病（虚线 + 问号，表示方向矛盾）
    ax.annotate("", xy=(0.69, 0.70), xytext=(0.63, 0.70),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color=fs.C_GRAY, lw=2.0,
                                linestyle="--"))
    ax.text(0.66, 0.745, "DIRECTION\nCONFLICT?", ha="center", fontsize=8.5,
            transform=ax.transAxes, style="italic", color=fs.C_GRAY,
            weight="bold")

    # 底部说明框
    ax.add_patch(plt.Rectangle((0.03, 0.10), 0.94, 0.30, facecolor=fs.C_GRAY,
                               alpha=0.08, edgecolor=fs.C_GRAY, linewidth=1.0,
                               transform=ax.transAxes))
    ax.text(0.05, 0.34,
            "Honest reading: CTD confirms SCRN2 is an EDC-responsive gene, but the "
            "down-regulation direction CONFLICTS with\nMR (where HIGH expression raises "
            "disease risk). EDC thyroid toxicity may NOT be mediated through SCRN2 —\n"
            "SCRN2's 'EDC-response' and 'disease-causation' roles are directionally "
            "uncoupled. Treat as candidate clue, not proven chain.",
            ha="left", va="center", fontsize=8.8, transform=ax.transAxes,
            color=fs.C_DARK, linespacing=1.6)

    # EDC 类别小图例（左下角外）
    legend_items = [
        ("PFAS", fs.C_BLUE), ("Bisphenol", fs.C_CYAN), ("Arsenic", fs.C_GREEN),
        ("PBDE", fs.C_YELLOW), ("Phthalate", fs.C_PURPLE),
    ]
    lx = 0.03
    for name, c in legend_items:
        ax.add_patch(plt.Rectangle((lx, 0.035), 0.028, 0.05, facecolor=c,
                                   edgecolor=c, transform=ax.transAxes))
        ax.text(lx + 0.034, 0.06, name, fontsize=7.5, va="center",
                transform=ax.transAxes)
        lx += 0.085

    fig.suptitle("SCRN2 as an EDC-responsive gene: evidence and the directionality caveat",
                 fontsize=11.5, y=0.99)
    fig.tight_layout()
    out_png = os.path.join(FIG, "Fig09_SCRN2_EDC_chain.png")
    out_pdf = os.path.join(FIG, "Fig09_SCRN2_EDC_chain.pdf")
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    print(f"\n图已保存: {out_png}")
    print(f"图已保存: {out_pdf}")

    print("\n=== EDC 类别分布 ===")
    for k, v in Counter(r["edc_class"] for r in edc_rows).most_common():
        print(f"  {k}: {v}")
    print("\n=== 人类 EDC 方向分布 ===")
    for k, v in Counter(r["direction"] for r in human_edc).most_common():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
