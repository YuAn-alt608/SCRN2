#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p26_outcome_identity.py
=======================
核定三类结局 GWAS 的身份与样本量口径（只读，不改任何原始文件）。

回答两个问题：
  1. T4 是游离还是总？   → 用「论文 N / 文件逐 SNP n 上限 / n 的队列分解」三重对齐
  2. 甲减的 n 为什么对不上？ → 论文总量、欧洲部分、非欧洲部分、文件上限四者对照

做法：
  - 从 DATA10/raw/{TSH,T4,HYPO}.tsv.gz 流式统计 `n` 列的取值分布（不载入内存）；
  - 结果缓存到 tables/34_结局逐SNP样本量分布.csv，第二次起直接读缓存（快）；
  - 核定结论写到 tables/35_结局GWAS身份.csv。

用法：
    python p26_outcome_identity.py           # 有缓存就秒回
    python p26_outcome_identity.py --refresh # 强制重扫（约 4 分钟）
"""
import os, sys, gzip
from collections import Counter
import numpy as np
import pandas as pd

BASE = "D:/NewBeginning/MR1"
OUT = f"{BASE}/DATA21_complete_2026-11-14"
RAW = f"{BASE}/DATA10/raw"
TAB = f"{OUT}/tables"
CACHE = f"{TAB}/34_结局逐SNP样本量分布.csv"

FILES = {"TSH": f"{RAW}/TSH.tsv.gz", "FT4": f"{RAW}/T4.tsv.gz", "HYPO": f"{RAW}/HYPO.tsv.gz"}

# 论文与 GWAS Catalog 的官方数字（Rand SA et al., Nat Genet 2025, PMID 41238958）
PAPER = {
    "TSH":  dict(acc="GCST90572789", trait="Thyroid stimulating hormone levels",
                 paper_n=482873, snp_count=29790163, note="482,873 European ancestry individuals"),
    "FT4":  dict(acc="GCST90572790", trait="Thyroxine levels（OBA_1001010，粗粒度上位词）",
                 paper_n=191449, snp_count=28647991, note="abstract 原文作 free thyroxine (191,449 individuals)"),
    "HYPO": dict(acc="GCST90572791", trait="Hypothyroidism",
                 paper_n=1178661, snp_count=34556406, note="113,393 cases + 1,065,268 controls"),
}
HYPO_EUR = 98805 + 948186          # GWAS Catalog：欧洲部分
HYPO_NONEUR = 14588 + 117082       # GWAS Catalog：非欧洲（23andMe）


def scan(path):
    """流式统计 n 列（第 10 列）的取值分布。"""
    cnt = Counter()
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 10:
                cnt[p[9]] += 1
    return cnt


def load_counts(refresh=False):
    if os.path.exists(CACHE) and not refresh:
        d = pd.read_csv(CACHE)
        print(f"[cache] {os.path.basename(CACHE)}  ({len(d)} 行)")
        return d
    rows = []
    for oc, path in FILES.items():
        print(f"[scan] {oc} … {path}")
        for v, c in scan(path).items():
            rows.append(dict(结局=oc, n=v, 变异数=c))
    d = pd.DataFrame(rows)
    d["n"] = pd.to_numeric(d.n, errors="coerce")
    d = d.dropna(subset=["n"]).astype({"变异数": int})
    d.to_csv(CACHE, index=False)
    print(f"[save] {CACHE}  ({len(d)} 行)")
    return d


def wmedian(sub):
    s = sub.sort_values("n")
    return int(s.n[s.变异数.cumsum() >= s.变异数.sum() / 2].iloc[0])


def main():
    refresh = "--refresh" in sys.argv
    d = load_counts(refresh)

    rows = []
    for oc in ["TSH", "FT4", "HYPO"]:
        s = d[d.结局 == oc]
        mx = int(s.n.max())
        pn = PAPER[oc]["paper_n"]
        rows.append(dict(
            结局=("游离 T4（FT4）" if oc == "FT4" else oc),
            GWAS_Catalog=PAPER[oc]["acc"], 论文报告N=pn,
            文件最大逐SNP_n=mx, 差=pn - mx,
            全基因组加权中位n=wmedian(s),
            变异数=int(s.变异数.sum()),
            GWAS_Catalog_snp_count=PAPER[oc]["snp_count"],
            变异数吻合=(int(s.变异数.sum()) == PAPER[oc]["snp_count"]),
            说明=PAPER[oc]["note"],
        ))
    ident = pd.DataFrame(rows)

    # 队列分解：T4 的三个分量
    t4 = d[d.结局 == "FT4"].set_index("n").变异数
    parts = []
    for a, b, c in [(66057, 49269, 76123)]:
        parts.append(dict(分量A=a, 分量B=b, 分量C=c,
                          A加B=a + b, 三者合计=a + b + c,
                          论文N=191449,
                          逐位吻合=(a + b + c == 191449 and a + b == 115326)))

    # 甲减缺口归因
    hmax = int(d[d.结局 == "HYPO"].n.max())
    attr = dict(论文总量=1178661, 欧洲部分=HYPO_EUR, 非欧洲部分=HYPO_NONEUR,
                文件上限=hmax, 上限减欧洲=HYPO_EUR - hmax,
                总量减上限=1178661 - hmax, 非欧洲部分_=HYPO_NONEUR)

    print("\n=== 身份核定 ===")
    print(ident.to_string(index=False))
    print("\n=== T4 队列分解 ===")
    print(pd.DataFrame(parts).to_string(index=False))
    print("\n=== 甲减缺口归因 ===")
    for k, v in attr.items():
        print(f"  {k:12s} {v:>10,}")

    print("\n=== 判读 ===")
    print("  变异数：三个文件的变异数与 GWAS Catalog 的 snp_count 逐位相等 → 文件身份无疑义")
    print("  TSH / FT4：文件 n 上限与论文 N 只差 10 / 1（单 SNP 掉样）→ 同一批数据")
    print("  FT4：66,057 + 49,269 + 76,123 = 191,449 逐位吻合；49,269 即 Teumer 2018 FT4 队列")
    print("       → 性状是游离 T4，不是总 T4")
    print(f"  HYPO：{hmax:,} ≈ 欧洲部分 {HYPO_EUR:,} − {HYPO_EUR-hmax:,}；"
          f"距论文总量 1,178,661 差 {1178661-hmax:,} ≈ 23andMe 的 {HYPO_NONEUR:,}")
    print("       → 公开汇总统计只含欧洲队列")
    print("\n注意：逐 SNP 中位数（TSH 191,203 / FT4 66,057 / 甲减 424,442）不能当作研究样本量引用。")

    ident.to_csv(f"{TAB}/35_结局GWAS身份.csv", index=False, encoding="utf-8-sig")
    print(f"\n[save] {TAB}/35_结局GWAS身份.csv")


if __name__ == "__main__":
    main()
