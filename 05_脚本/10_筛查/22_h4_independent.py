# -*- coding: utf-8 -*-
"""
Step 22: 把 71 个 PP.H4>0.8 位点按"循环性"分类，并导出真正独立的清单

循环性的定义（与报告 07 节一致）：
  A 类 工具变量循环   —— coloc 的 top SNP 同时是 eQTL 最强位点和 MR 工具变量（60 对）
  B 类 eQTL 顶位点重叠 —— top SNP 不是工具变量，但仍是 eQTL 最强位点（1 对，LCN2）
  C 类 独立           —— top SNP 既不是工具变量、也不是 eQTL 最强位点（10 对）

只有 C 类提供"独立于工具变量选择"的共定位证据。本步把 C 类单独导出，
并附上 PP.H4、top SNP 的置信度、分级、FinnGen 复现、GTEx 验证状态。

输出：
  results/h4_loci_classes.csv     71 行，带 class 列
  results/h4_independent_loci.csv 10 行，C 类清单（核心增量证据）
  results/h4_independent_summary.txt
"""
import os
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")

OC_NAME = {"TSH": "TSH", "T4": "Total T4", "HYPO": "Hypothyroidism",
           "HASH": "Hashimoto thyroiditis"}


def main():
    c = pd.read_csv(os.path.join(RES, "h4_circularity.csv"))
    c = c.rename(columns={"PP.H4": "pp_h4", "PP.H3": "pp_h3"})
    tier = pd.read_csv(os.path.join(RES, "evidence_tiers_v2.csv"), encoding="utf-8-sig")
    tier.columns = [x.strip() for x in tier.columns]
    tier = tier.rename(columns={"PP.H4": "pp_h4", "PP.H3": "pp_h3"})
    fg = pd.read_csv(os.path.join(RES, "finngen_replication.csv"), encoding="utf-8-sig")
    gx = pd.read_csv(os.path.join(RES, "gtex_thyroid_validation.csv"), encoding="utf-8-sig")

    c["top_is_instrument"] = c["top_is_instrument"].astype(str).str.lower() == "true"
    c["same_as_eqtl_top"] = c["same_as_eqtl_top"].astype(str).str.lower() == "true"

    def cls(r):
        if r.top_is_instrument:
            return "A_工具变量循环"
        if r.same_as_eqtl_top:
            return "B_eQTL顶位点"
        return "C_独立"

    c["locus_class"] = c.apply(cls, axis=1)

    tmap = {(r.gene, r.outcome): r for r in tier.itertuples()}
    rows = []
    for r in c.itertuples():
        t = tmap.get((r.gene, OC_NAME.get(r.outcome, r.outcome)))
        f = fg[(fg["gene"] == r.gene) & (fg["oc"] == r.outcome)]
        g = gx[gx["gene"] == r.gene]
        if len(g) and r.top_snp in set(g["snp"]):
            gstat = "同一 SNP"
        elif len(g):
            gstat = "基因有验证"
        else:
            gstat = "未验证"
        rows.append(dict(
            gene=r.gene, outcome=r.outcome, locus_class=r.locus_class,
            # 分级表只收录 FDR<0.05 的命中；H4>0.8 但 MR 未达门槛的记为"未达FDR"
            tier=(getattr(t, "tier", "") if t is not None else "") or "未达FDR",
            pp_h4=r.pp_h4, pp_h3=r.pp_h3,
            top_snp=r.top_snp, top_snp_PP=r.top_snp_PP,
            nsnps=r.nsnps, n_snp_eqtl_sig=r.n_snp_eqtl_sig,
            pval_out=r.pval_out, maf_median=r.maf_median,
            mr_beta=getattr(t, "b", None) if t is not None else None,
            mr_p=getattr(t, "pval", None) if t is not None else None,
            finngen_replicated=(bool(f["replicated_physio"].iloc[0]) if len(f) else None),
            gtex_thyroid=gstat,
        ))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(RES, "h4_loci_classes.csv"), index=False, encoding="utf-8-sig")

    ind = d[d["locus_class"] == "C_独立"].sort_values("pp_h4", ascending=False)
    ind.to_csv(os.path.join(RES, "h4_independent_loci.csv"), index=False,
               encoding="utf-8-sig")

    n = len(d)
    n_a = int((d["locus_class"] == "A_工具变量循环").sum())
    n_b = int((d["locus_class"] == "B_eQTL顶位点").sum())
    n_c = len(ind)
    strong = ind[ind["top_snp_PP"] > 0.5]
    lines = [
        f"PP.H4>0.8 共 {n} 对",
        f"  A 类 工具变量循环（top SNP = eQTL 最强 = 工具变量）: {n_a} ({100*n_a/n:.1f}%)",
        f"  B 类 top SNP 是 eQTL 最强位点、但不是工具变量: {n_b} ({100*n_b/n:.1f}%)",
        f"  C 类 独立（top SNP 与两者都不同）: {n_c} ({100*n_c/n:.1f}%)",
        "",
        ("C 类中 top SNP 后验 >0.5（单个变异可定位）: "
         + str(len(strong)) + " 对 —— "
         + ", ".join(f"{r.gene}x{r.outcome}(H4={r.pp_h4:.3f})" for r in strong.itertuples())
         if len(strong) else "C 类中无高置信 top SNP"),
        "",
        "C 类清单（按 PP.H4 降序）:",
    ]
    for r in ind.itertuples():
        fgs = ("复现" if r.finngen_replicated else "未复现") \
            if r.finngen_replicated is not None else "未入复现表"
        lines.append(
            f"  {r.gene:<8s} {r.outcome:<5s} H4={r.pp_h4:.3f} top={r.top_snp}"
            f"(PP={r.top_snp_PP:.2f}) tier={r.tier} FinnGen={fgs}")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(RES, "h4_independent_summary.txt"), "w", encoding="utf-8").write(txt + "\n")


if __name__ == "__main__":
    main()
