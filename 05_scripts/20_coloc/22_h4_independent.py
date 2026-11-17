# -*- coding: utf-8 -*-
"""
Step 22: classify the 71 loci with PP.H4>0.8 by circularity and export the
genuinely independent list

Circularity classes (identical to section 07 of the report):
  class A  instrument circularity -- the coloc top SNP is both the strongest
           eQTL signal and the MR instrument (60 pairs)
  class B  eQTL top-site overlap  -- the top SNP is not the instrument but is
           still the strongest eQTL signal (1 pair, LCN2)
  class C  independent            -- the top SNP is neither the instrument nor
           the strongest eQTL signal (10 pairs)

Only class C offers colocalization evidence that is independent of instrument
selection. This step exports class C with its PP.H4, top-SNP confidence, tier,
FinnGen replication and GTEx validation status.

Output:
  results/h4_loci_classes.csv     71 rows with a class column
  results/h4_independent_loci.csv 10 rows, the class C list (core added evidence)
  results/h4_independent_summary.txt
"""
import os
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")

OC_NAME = {"TSH": "TSH", "T4": "Free T4", "HYPO": "Hypothyroidism",
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
            return "A_instrument_circularity"
        if r.same_as_eqtl_top:
            return "B_eqtl_top_site"
        return "C_independent"

    c["locus_class"] = c.apply(cls, axis=1)

    tmap = {(r.gene, r.outcome): r for r in tier.itertuples()}
    rows = []
    for r in c.itertuples():
        t = tmap.get((r.gene, OC_NAME.get(r.outcome, r.outcome)))
        f = fg[(fg["gene"] == r.gene) & (fg["oc"] == r.outcome)]
        g = gx[gx["gene"] == r.gene]
        if len(g) and r.top_snp in set(g["snp"]):
            gstat = "same_snp"
        elif len(g):
            gstat = "gene_validated"
        else:
            gstat = "not_validated"
        rows.append(dict(
            gene=r.gene, outcome=r.outcome, locus_class=r.locus_class,
            # the tier table only lists FDR<0.05 hits; H4>0.8 without MR
            # significance is recorded as below_FDR
            tier=(getattr(t, "tier", "") if t is not None else "") or "below_FDR",
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

    ind = d[d["locus_class"] == "C_independent"].sort_values("pp_h4", ascending=False)
    ind.to_csv(os.path.join(RES, "h4_independent_loci.csv"), index=False,
               encoding="utf-8-sig")

    n = len(d)
    n_a = int((d["locus_class"] == "A_instrument_circularity").sum())
    n_b = int((d["locus_class"] == "B_eqtl_top_site").sum())
    n_c = len(ind)
    strong = ind[ind["top_snp_PP"] > 0.5]
    lines = [
        f"PP.H4>0.8 in {n} pairs",
        f"  class A instrument circularity (top SNP = strongest eQTL = instrument): {n_a} ({100*n_a/n:.1f}%)",
        f"  class B top SNP is the strongest eQTL site but not the instrument: {n_b} ({100*n_b/n:.1f}%)",
        f"  class C independent (top SNP differs from both): {n_c} ({100*n_c/n:.1f}%)",
        "",
        ("class C with top-SNP posterior >0.5 (single variant localizable): "
         + str(len(strong)) + " pairs -- "
         + ", ".join(f"{r.gene}x{r.outcome}(H4={r.pp_h4:.3f})" for r in strong.itertuples())
         if len(strong) else "class C has no high-confidence top SNP"),
        "",
        "class C list (by descending PP.H4):",
    ]
    for r in ind.itertuples():
        fgs = ("replicated" if r.finngen_replicated else "not replicated") \
            if r.finngen_replicated is not None else "not in replication table"
        lines.append(
            f"  {r.gene:<8s} {r.outcome:<5s} H4={r.pp_h4:.3f} top={r.top_snp}"
            f"(PP={r.top_snp_PP:.2f}) tier={r.tier} FinnGen={fgs}")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(RES, "h4_independent_summary.txt"), "w", encoding="utf-8").write(txt + "\n")


if __name__ == "__main__":
    main()
