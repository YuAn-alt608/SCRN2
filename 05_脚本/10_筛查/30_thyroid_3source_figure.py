# -*- coding: utf-8 -*-
"""
30_thyroid_3source_figure.py
生成三来源甲状腺 eQTL 验证的对比图。
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "D:/NewBeginning/MR1/DATA10"
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False

def main():
    d = pd.read_csv(os.path.join(RES, "thyroid_3source_validation.csv"),
                    encoding="utf-8-sig")
    n = len(d)

    sources = [
        ("来源1\nZenodo(salmon)", "s1_p", "s1_same"),
        ("来源2\neQTL Catalogue", "s2_p", "s2_same"),
        ("来源3\nGTEx官方API", "s3_p", "s3_same"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # 面板A：三个来源的匹配/显著/方向
    labels = ["可匹配", "显著", "方向一致"]
    for i, (name, pcol, samecol) in enumerate(sources):
        ax = axes[0]
        has = d[pcol].notna().sum()
        sig = (d[pcol] < 0.05).sum()
        same = d[samecol].sum()
        vals = [has, sig, same]
        x = np.arange(3) + i * 0.25
        ax.bar(x, vals, width=0.22, label=name.replace("\n", " "))
    axes[0].set_xticks(np.arange(3) + 0.25)
    axes[0].set_xticklabels(labels)
    axes[0].set_ylabel("SNP 数")
    axes[0].set_title(f"A. 三来源检出对比（共 {n} 个工具变量）")
    axes[0].legend(fontsize=8)
    axes[0].axhline(n, ls="--", color="gray", alpha=0.5)

    # 面板B：方向一致率（分全部可匹配 vs 仅显著子集）
    ax = axes[1]
    rates_all, rates_sig = [], []
    for name, pcol, samecol in sources:
        has = d[pcol].notna()
        sig = has & (d[pcol] < 0.05)
        denom_all = d[samecol][has].notna().sum()
        denom_sig = d[samecol][sig].notna().sum()
        same_all = d[samecol][has].sum()
        same_sig = d[samecol][sig].sum()
        rates_all.append((same_all / denom_all * 100) if denom_all else 0)
        rates_sig.append((same_sig / denom_sig * 100) if denom_sig else 0)
    x = np.arange(3)
    w = 0.36
    b1 = ax.bar(x - w/2, rates_all, w, label="全部可匹配", color="#4C72B0")
    b2 = ax.bar(x + w/2, rates_sig, w, label="仅显著子集", color="#C44E52")
    ax.set_xticks(x)
    ax.set_xticklabels([s[0].split("\n")[1] for s in sources])
    ax.set_ylabel("方向一致率 (%)")
    ax.set_ylim(0, 100)
    ax.set_title("B. 与全血方向一致性")
    ax.legend(fontsize=8)
    for b, r in zip(b1, rates_all):
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 2,
                f"{r:.0f}%", ha="center", fontsize=8)
    for b, r in zip(b2, rates_sig):
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 2,
                f"{r:.0f}%", ha="center", fontsize=8)

    # 面板C：三来源显著重叠（Venn 精确计数）
    ax = axes[2]
    s1 = set(d[d.s1_p < 0.05].snp)
    s2 = set(d[d.s2_p < 0.05].snp)
    s3 = set(d[d.s3_p < 0.05].snp)
    three = len(s1 & s2 & s3)
    only1 = (len(s1 - s2 - s3) + len(s2 - s1 - s3) + len(s3 - s1 - s2))
    two = (len(s1 & s2) + len(s1 & s3) + len(s2 & s3)
           - 3 * three)
    cats = ["仅1个来源", "2个来源", "3个来源"]
    ax.bar(range(3), [only1, two, three], color=["#DD8452", "#8C8C8C", "#C44E52"])
    ax.set_xticks(range(3))
    ax.set_xticklabels(cats)
    ax.set_ylabel("显著 SNP 数")
    ax.set_title("C. 三来源显著重叠")
    for i, v in enumerate([only1, two, three]):
        ax.text(i, v + 0.3, str(v), ha="center", fontsize=11)

    fig.suptitle("甲状腺组织 eQTL 三来源验证：A/B 级工具变量的组织特异性",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(FIG, "Fig11_thyroid_3source.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"已写 {out}")

    # 文本摘要
    summary = [
        "三来源甲状腺 eQTL 验证（A/B 级 70 个工具变量，统一口径：该 SNP 对该基因）",
        "=" * 60,
    ]
    for name, pcol, samecol in sources:
        has = d[pcol].notna()
        sig = has & (d[pcol] < 0.05)
        same_all = d[samecol][has].sum()
        denom_all = d[samecol][has].notna().sum()
        same_sig = d[samecol][sig].sum()
        denom_sig = d[samecol][sig].notna().sum()
        summary.append(
            f"{name.replace(chr(10),' ')}: 可匹配 {int(has.sum())}/70, "
            f"显著 {int(sig.sum())}, 方向一致(全部) {int(same_all)}/{int(denom_all)}, "
            f"方向一致(显著) {int(same_sig)}/{int(denom_sig)}")
    summary.append("")
    summary.append("核心结论：全血强 eQTL 在甲状腺里大部分不再是原基因的显著 eQTL，")
    summary.append("组织错配得到量化证实；即使甲状腺里仍显著的位点，方向一致性中等。")
    with open(os.path.join(RES, "thyroid_3source_summary.txt"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(summary))
    print("已写", os.path.join(RES, "thyroid_3source_summary.txt"))


if __name__ == "__main__":
    main()
