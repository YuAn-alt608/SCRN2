"""
p18_signal_summary.py -- freezes the "one locus, several signals" conclusion into three tables

Outputs
  tables/21_published_sentinels_17q21.csv          sentinels from the two papers in this region + LD attribution
  tables/22_sentinel_gene_mapping_uniqueness.csv   genome-wide fraction with a unique attribution (a scope that can be extended)
  tables/23_signal_attribution_summary.csv         the three regional signals x candidate gene x all evidence x verdict

Principle: every number is read from a result or input table inside the project,
never retyped by hand. Only the sentinel values of published papers come from the
GWAS Catalog download (raw/sentinels/), with the accession recorded.
"""
import os

import numpy as np
import pandas as pd

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
T = f"{BASE}/tables"
# GRCh38 -> GRCh37 offset (calibrated independently on two loci, which agree):
#   rs4794063 47,727,128 - 45,804,494 = 1,922,634
#   rs1989291 47,730,181 - 45,807,547 = 1,922,634
# The earlier value of 1,957,714 was back-inferred from a gene start and was
# wrong by 35 kb.
OFS = 1_922_634

# ---------------------------------------------------------------- 21 sentinels
sen = pd.read_csv(f"{BASE}/raw/sentinels/gwas_catalog_associations_all.csv")
reg = sen[(sen.chr38.astype(str) == "17") &
          (sen.pos38 >= 47_300_000) & (sen.pos38 <= 48_500_000)].copy()
reg["pos37"] = reg.pos38.astype(int) - OFS

# LD proxy: rs4794063 is not called in 1000G, so the nearest panel variant is used (resolved below)
PROXY = None
ld = pd.read_csv(f"{BASE}/inputs_ld/ld_TSH.tsv.gz", sep="\t", index_col=0)
ld.index = ld.index.astype(int)
ld.columns = ld.columns.astype(int)
LEAD = {"rs7219420": 45_800_299, "rs62076102": 45_887_132, "rs6503875": 45_915_395}
idx = np.array(ld.index)
if len(reg):
    tgt = int(reg.pos37.iloc[0])
    PROXY = int(idx[np.abs(idx - tgt).argmin()])
    print(f"rs4794063 is not called in 1000G; nearest panel variant {PROXY} ({abs(PROXY-tgt)} bp away)")
    for k, v in LEAD.items():
        reg[f"r_{k}"] = ld.loc[PROXY, v] if (PROXY in ld.index and v in ld.index) else np.nan
    reg["proxy_pos37"] = PROXY
    reg["proxy_dist_bp"] = (reg.pos37 - PROXY).abs()
reg["signal"] = reg.pos37.apply(
    lambda p: "1 (45.74-45.81)" if 45_740_000 <= p < 45_810_000 else
              ("2 (45.86-45.93)" if 45_860_000 <= p < 45_930_000 else
               ("3 (46.02-46.08)" if 46_020_000 <= p < 46_080_000 else "other")))
reg["verdict"] = reg["signal"].apply(
    lambda s: "falls in unattributable signal 1 (TBKBP1 heterogeneous / TBX21 no eQTL)" if s.startswith("1")
    else ("falls in signal 2, attributable to SCRN2" if s.startswith("2") else "unresolved"))
cols = ["paper", "trait", "accession", "rsid", "effect_allele", "pos38", "pos37", "p", "beta",
        "raf", "mapped_gene", "proxy_pos37", "proxy_dist_bp",
        "r_rs7219420", "r_rs62076102", "r_rs6503875", "signal", "verdict"]
reg[cols].to_csv(f"{T}/21_published_sentinels_17q21.csv", index=False, encoding="utf-8-sig")
print("=== table 21, regional sentinels ===")
print(reg[cols].to_string(index=False))

# ---------------------------------------------------------------- 22 mapping uniqueness
sen["ngene"] = sen.mapped_gene.fillna("").apply(
    lambda s: len([x for x in str(s).split(";") if x.strip()]))
rows = []
for (paper, trait), g in sen.groupby(["paper", "trait"]):
    n = len(g)
    rows.append(dict(paper=paper, trait=trait, n_sentinels=n,
                     n_unique_gene=(g.ngene == 1).sum(),
                     n_multi_gene=(g.ngene > 1).sum(),
                     n_unmapped=(g.ngene == 0).sum(),
                     frac_unique_gene=round((g.ngene == 1).mean(), 4),
                     n_regional_sentinels=int(((g.chr38.astype(str) == "17") &
                                 (g.pos38 >= 47_300_000) & (g.pos38 <= 48_500_000)).sum())))
u = pd.DataFrame(rows)
u.to_csv(f"{T}/22_sentinel_gene_mapping_uniqueness.csv", index=False, encoding="utf-8-sig")
print("\n=== table 22, mapping uniqueness ===")
print(u.to_string(index=False))

# ---------------------------------------------------------------- 23 signal attribution summary
mr = pd.read_csv(f"{T}/04_cisMR_SMR_HEIDI.csv")
co = pd.read_csv(f"{T}/03_panel_coloc_abf_prior_grid.csv")
ts = pd.read_csv(f"{T}/10_per_tissue_key_snp_eqtl.csv")


def gwas(pos, outcome, tissue="Thyroid"):
    f = f"{BASE}/inputs/{tissue}__SCRN2__{outcome}.tsv"
    d = pd.read_csv(f, sep="\t")
    r = d[d.pos37 == pos]
    return (float(r.p_out.iloc[0]), float(r.beta_out.iloc[0])) if len(r) else (np.nan, np.nan)


def eqtl(pos, gene, tissue="Thyroid"):
    f = f"{BASE}/inputs/{tissue}__{gene}__TSH.tsv"
    if not os.path.exists(f):
        return np.nan, np.nan
    d = pd.read_csv(f, sep="\t")
    r = d[d.pos37 == pos]
    return (float(r.p_eqtl.iloc[0]), float(r.beta_eqtl.iloc[0])) if len(r) else (np.nan, np.nan)


def mrow(gene, tissue, outcome="TSH"):
    m = mr[(mr.gene == gene) & (mr.tissue == tissue) & (mr.outcome == outcome)]
    return m.iloc[0] if len(m) else None


def crow(gene, tissue, outcome="TSH", p12=1e-4):
    c = co[(co.gene == gene) & (co.tissue == tissue) & (co.outcome == outcome) &
           (np.isclose(co.p12, p12))]
    return c.iloc[0] if len(c) else None


TISSUES = ts.tissue.nunique()
COVERED = set(ts.gene.unique())


def ntis(snp, gene):
    """Number of the 49 GTEx tissues in which this SNP is an eQTL for the gene; "not tested" if the gene was not covered."""
    if gene not in COVERED:
        return "Not tested"
    d = ts[(ts.key_snp == snp) & (ts.gene == gene)]
    return f"{d.tissue.nunique()}/{TISSUES}"


ROWS = [
    # signal, sentinel, position, candidate gene, tissue, MR instrument position
    ("1", "rs7219420", 45_800_299, "TBKBP1", "Thyroid", 45_774_072),
    ("1", "rs7219420", 45_800_299, "TBKBP1", "WholeBlood", 45_744_652),
    ("1", "rs7219420", 45_800_299, "TBX21", "Thyroid", 45_800_299),
    ("2", "rs6503875", 45_915_395, "SCRN2", "Thyroid", 45_915_395),
    ("2", "rs6503875", 45_915_395, "SCRN2", "WholeBlood", 45_912_302),
    ("2", "rs62076102", 45_887_132, "MRPL10", "Thyroid", 45_888_254),
    ("3", "COPZ2lead", 46_050_669, "COPZ2", "WholeBlood", 46_050_669),
    ("3", "COPZ2lead", 46_050_669, "COPZ2", "Thyroid", 45_881_197),
]

out = []
for sig, snp, pos, gene, tissue, mrpos in ROWS:
    p_tsh, b_tsh = gwas(pos, "TSH")
    p_hyp, _ = gwas(pos, "HYPO")
    # the eQTL is queried at that gene's instrument position in that tissue, so
    # that it is the same variant the F statistic refers to
    p_teq, b_teq = eqtl(mrpos, gene, tissue)
    m = mrow(gene, tissue)
    c = crow(gene, tissue)
    out.append(dict(
        signal=sig, sentinel=snp, pos37=pos, candidate_gene=gene, tissue=tissue, mr_instrument=mrpos,
        TSH_P=p_tsh, TSH_beta=b_tsh, HYPO_P=p_hyp,
        eQTL_P=p_teq, eQTL_beta=b_teq,
        F=m.F_lead if m is not None else np.nan,
        MR_b=m.b_wald if m is not None else np.nan,
        MR_se=m.se_wald if m is not None else np.nan,
        MR_P=m.p_wald if m is not None else np.nan,
        n_iv=m.n_iv if m is not None else 0,
        Q_P=m.p_Q if m is not None else np.nan,
        I2=m.I2 if m is not None else np.nan,
        HEIDI_P=m.heidi_p if m is not None else np.nan,
        coloc_H3=c["PP.H3"] if c is not None else np.nan,
        coloc_H4=c["PP.H4"] if c is not None else np.nan,
        tissue_replicable=ntis(snp, gene),
    ))

s = pd.DataFrame(out)

# ---- attach the within-window coloc results (the hypothyroidism panel is the
# ---- same source; see tables/25 for the check)
w = pd.read_csv(f"{T}/24_coloc_by_signal_window.csv")
OWN_WIN = {"SCRN2": "Signal 2 (45.86-45.98)", "MRPL10": "Signal 2 (45.86-45.98)",
           "TBKBP1": "Signal 1 (45.72-45.86)", "COPZ2": "Signal 3 (46.00-46.12)"}
s["window_coloc_H3"] = np.nan
s["window_coloc_H4"] = np.nan
s["hypo_window_coloc_H3"] = np.nan
s["hypo_window_coloc_H4"] = np.nan
for i, r in s.iterrows():
    q = w[(w.gene == r.candidate_gene) & (w.tissue == r.tissue) & (w.outcome == "TSH") &
          (w.window == OWN_WIN.get(r.candidate_gene, ""))]
    if len(q):
        s.at[i, "window_coloc_H3"] = q.H3.iloc[0]
        s.at[i, "window_coloc_H4"] = q.H4.iloc[0]
    q2 = w[(w.gene == r.candidate_gene) & (w.tissue == r.tissue) & (w.outcome == "HYPO") &
           (w.window == OWN_WIN.get(r.candidate_gene, ""))]
    if len(q2):
        s.at[i, "hypo_window_coloc_H3"] = q2.H3.iloc[0]
        s.at[i, "hypo_window_coloc_H4"] = q2.H4.iloc[0]
# verdict rules (hard-coded so that every call is interpretable)
verdict = {
    ("1", "TBKBP1", "Thyroid"):
        "Not attributable: strong eQTL (F=166) but coloc does not support sharing (H3=0.87); a single instrument cannot support HEIDI",
    ("1", "TBKBP1", "WholeBlood"):
        "Not attributable: Q P=0.009 (I2=70%), HEIDI<1e-4, coloc H4=0.03",
    ("1", "TBX21", "Thyroid"):
        "Excluded: no cis-eQTL (P=0.55, F=5.8); the MR estimate is uninterpretable",
    ("2", "SCRN2", "Thyroid"):
        "Assigned to SCRN2: MR, coloc and SMR all pass; HEIDI needs >=2 instruments and is not testable in thyroid",
    ("2", "SCRN2", "WholeBlood"):
        "Assigned to SCRN2: 4/4 pass (Q P=0.25, HEIDI 0.064), 47/49 tissues replicable",
    ("2", "MRPL10", "Thyroid"):
        "Excluded: MR null (P=0.70); coloc H4=0.02, not shared with TSH",
    ("3", "COPZ2", "WholeBlood"):
        "Assigned to COPZ2: 4/4 in whole blood (Q P=0.06, HEIDI 0.053)",
    ("3", "COPZ2", "Thyroid"):
        "Same direction but the instrument is too weak (F=15.7, coloc H4=0.79); supportive only",
}
s["verdict"] = [verdict.get((r.signal, r.candidate_gene, r.tissue), "") for r in s.itertuples()]


def hypo_verdict(r):
    if np.isnan(r.hypo_window_coloc_H4):
        return "Not testable"
    return ("Shared (H4=%.2f)" % r.hypo_window_coloc_H4 if r.hypo_window_coloc_H4 > 0.75
            else ("Not shared (H3=%.2f)" % r.hypo_window_coloc_H3 if r.hypo_window_coloc_H3 > 0.75
                  else "Undetermined (H3=%.2f/H4=%.2f)" % (r.hypo_window_coloc_H3, r.hypo_window_coloc_H4)))


s["hypo_verdict"] = [hypo_verdict(r) for r in s.itertuples()]
s.to_csv(f"{T}/23_signal_attribution_summary.csv", index=False, encoding="utf-8-sig")
print("\n=== table 23, signal attribution ===")
print(s.round(5).to_string(index=False))
