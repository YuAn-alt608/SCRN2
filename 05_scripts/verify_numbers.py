#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
verify_numbers.py - regenerate every statistic reported in the manuscript
==========================================================================
Usage:  python 05_scripts/verify_numbers.py
Output: 05_scripts/verification_report.txt  (also printed to stdout)

Design
------
1. Every number is recomputed from a table in 04_data/ (or 06_supplementary/).
   No intermediate document and no hand-typed value is used as a source.
2. Every percentage is printed with its numerator and denominator, so a
   reader can check the arithmetic without re-running anything.
3. Where a value is read from a deposited file rather than recomputed here,
   or comes from an external dataset (cohort GWAS), the entry names the file
   or accession it is read from.

Requires: python 3.9+, pandas, numpy, scipy.  No network access.
"""
import io
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # repository root
D = os.path.join(ROOT, "04_data")
SUP = os.path.join(ROOT, "06_supplementary")

LINES = []


def say(s=""):
    print(s)
    LINES.append(str(s))


def title(s):
    say("")
    say("=" * 78)
    say(s)
    say("=" * 78)


def f(x, n=2):
    try:
        return ("%." + str(n) + "f") % x
    except Exception:
        return str(x)


def pct(a, b, n=1):
    return "%s%% (%s/%s)" % (f(100.0 * a / b, n), a, b) if b else "n/a"


def P(name):
    return os.path.join(D, name)


title("SCRN2 / thyroid set-point - numerical verification report")
say("repository root : %s" % ROOT)
say("generated       : %s" % pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"))

# ------------------------------------------------------------------ 1 screen
title("1  Screen input: EDCs -> genes -> tested pairs")
chem = pd.read_csv(P("01_screen/edc_chemicals.csv"))
genes = [l.strip() for l in open(P("01_screen/edc_thyroid_genes.txt"), encoding="utf-8") if l.strip()]
mr = pd.read_csv(P("01_screen/MR_all_results.csv"))
uniq = mr.drop_duplicates(["gene", "outcome"])
tiers = pd.read_csv(P("01_screen/evidence_tiers.csv"))
deg = pd.read_csv(P("01_screen/gene_degree_all.csv")).set_index("gene")
ivs = pd.read_csv(P("01_screen/iv_summary.csv"))

say("EDC chemicals                            : %d" % len(chem))
say("EDC-responsive genes expressed in thyroid : %d" % len(genes))
say("MR table rows                            : %d  (49 pairs carry a second estimator)" % len(mr))
say("unique gene-outcome pairs                : %d  <- the denominator used in the text" % len(uniq))
say("duplicate rows                           : %d" % (len(mr) - len(uniq)))
say("genes / outcomes tested                  : %d / %d" % (uniq.gene.nunique(), uniq.outcome.nunique()))
say("outcomes per gene                        : %s"
    % uniq.groupby("gene").outcome.nunique().value_counts().sort_index().to_dict())
ivar = pd.read_csv(P("01_screen/ivs.tsv"), sep=None, engine="python")
say("instrumented genes                       : %d  (variants: %d)" % (len(ivs), int(ivs.n_snp.sum())))
say("genes with a single instrument           : %s" % pct(int((ivs.n_snp == 1).sum()), len(ivs)))
say("median F over the %d instrumenting variants : %s  <- text says 338.5"
    % (len(ivar), f(ivar.f_stat.median(), 1)))

sig = uniq[uniq.fdr < 0.05]
say("")
say("FDR < 0.05                               : %d pairs, %d genes" % (len(sig), sig.gene.nunique()))
say("share of all tested pairs                : %s  <- text says 5.8%%" % pct(len(sig), len(uniq), 2))
say("evidence_tiers rows                      : %d  (matches significant pairs: %s)"
    % (len(tiers), set(zip(tiers.gene, tiers.outcome)) == set(zip(sig.gene, sig.outcome))))
say("tier distribution                        : %s" % tiers.tier.value_counts().to_dict())
say("significant pairs by outcome             : %s" % tiers.groupby("outcome").size().to_dict())

# ------------------------------------------------------- 1b EDC burden
title("1b EDC burden vs association strength")
mr2 = uniq.copy()
mr2["z"] = mr2.b / mr2.se
degv = deg.loc["SCRN2", "n_edc_chemicals"] if "SCRN2" in deg.index else float("nan")
say("EDCs mapped to SCRN2                     : %s" % degv)
for th in [3, 5, 10, 15, 20]:
    s = deg[(deg.n_edc_chemicals >= th) & (deg.thyroid_expressed == 1)]
    say("  threshold >=%2d : %5d genes ; SCRN2 retained: %s" % (th, len(s), "SCRN2" in set(s.index)))
say("")
for oc in ["TSH", "Free T4", "Hypothyroidism", "Hashimoto thyroiditis"]:
    d = mr2[mr2.outcome == oc].copy()
    d["deg"] = deg.reindex(d.gene).n_edc_chemicals.values
    d = d.dropna(subset=["deg", "z"])
    rho, pp = stats.spearmanr(d.deg, d.z.abs())
    q = d.groupby(pd.qcut(d.deg, 5, duplicates="drop"), observed=True).apply(
        lambda x: 100.0 * (x.fdr < 0.05).mean())
    hit, non = d[d.fdr < 0.05].deg, d[d.fdr >= 0.05].deg
    say("  %-22s n=%4d  rho=%+.3f (P=%.3f)  median EDC %g vs %g (MWU P=%.3g)"
        % (oc, len(d), rho, pp, hit.median(), non.median(), stats.mannwhitneyu(hit, non).pvalue))
    say("  %-22s quintile %% significant: %s"
        % ("", ", ".join(f(v, 1) for v in q.values)))

# ------------------------------------------------------------ 2 five layers
title("2  Five evidence layers (TSH branch)")
rc = json.load(open(P("09_layers/region_counts.json")))
for k in sorted(rc, key=int):
    g = rc[k]["genes"] if isinstance(rc[k], dict) else rc[k]
    say("  bitmask %-3s n=%-4d %s" % (k, len(g), ", ".join(g[:8]) + ("..." if len(g) > 8 else "")))
say("")
say("layer sizes reported in the text : 337 / 180 / 26 / 36 / 4")
say("intersection of all five layers  : %s" % rc["31"]["genes"])
say("pass layers 1-3, fail layer 4    : %s" % rc["7"]["genes"])
say("pass layers 1-4, fail layer 5    : %s" % rc["15"]["genes"])
n_occ = sum(1 for k, v in rc.items() if (v["genes"] if isinstance(v, dict) else v))
say("non-empty combinations occupied  : %d of 31" % n_occ)

# --------------------------------------------------------------- 3 FinnGen
title("3  FinnGen R12 consistency check")
fin = pd.read_csv(P("03_replication/finngen_replication.csv"))
say("rows: %d ; by outcome: %s" % (len(fin), fin.groupby("outcome").size().to_dict()))
say("")
say(open(P("03_replication/finngen_null_control.txt"), encoding="utf-8").read().strip())
say("")
for oc, g in fin.groupby("outcome"):
    testable = g[g.p_fg.notna()] if "p_fg" in g else g.iloc[:0]
    if not len(testable):
        continue
    same = int(testable.physio_ok.sum())
    ok = int(testable.replicated_physio.sum())
    say("  %-24s testable %4d ; expected direction %s ; recovered %s"
        % (oc, len(testable), pct(same, len(testable)), pct(ok, len(testable))))
r = {}
for oc, g in fin.dropna(subset=["b_orig", "b_fg"]).groupby("outcome"):
    if len(g) > 2:
        r[oc] = stats.pearsonr(g.b_orig, g.b_fg).statistic
say("")
say("discovery-FinnGen correlation (r) by outcome: %s" % {k: f(v, 2) for k, v in r.items()})

say("")
say("pooled over outcomes, by evidence tier (section 2.2):")
mp = {"TSH": "TSH", "Hypothyroidism": "HYPO", "Free T4": "T4", "Hashimoto thyroiditis": "HASH"}
tt = tiers[["gene", "outcome", "tier"]].copy()
tt["oc2"] = tt.outcome.map(mp)
ff = fin.rename(columns={"oc": "oc2"})[["gene", "oc2", "p_fg", "physio_ok", "replicated_physio"]]
mm = ff.merge(tt, on=["gene", "oc2"], how="inner")
mm = mm[mm.p_fg.notna()]
for tier, g in mm.groupby("tier"):
    say("  %-14s testable %3d ; recovered %s ; expected direction %s"
        % (tier, len(g), pct(int(g.replicated_physio.sum()), len(g)),
           pct(int(g.physio_ok.sum()), len(g))))

# ------------------------------------------------------------- 4 coloc
title("4  Colocalisation and the circularity audit")
co = pd.read_csv(P("02_coloc/coloc_results.csv"))
say("pairs tested                : %d (genes %d, outcomes %d)"
    % (len(co), co.gene.nunique(), co.outcome.nunique()))
h4 = co[co["PP.H4"] > 0.8].copy()
mp = {"TSH": "TSH", "Hypothyroidism": "HYPO", "Free T4": "T4", "Hashimoto thyroiditis": "HASH"}
sigset = set(zip(tiers.gene, tiers.outcome.map(mp)))
co["issig"] = [(g, o) in sigset for g, o in zip(co.gene, co.outcome)]
h4 = co[co["PP.H4"] > 0.8]
say("pairs with PP.H4 > 0.8      : %d  (%d in the FDR set, %d not)"
    % (len(h4), int((h4.issig).sum()), int((~h4.issig).sum())))
say("control pairs without FDR   : %d" % int((~co.issig).sum()))

cy = pd.read_csv(P("02_coloc/h4_circularity.csv"))
for c, lab in [("top_is_instrument", "peak variant is the MR instrument"),
               ("same_as_eqtl_top", "peak variant equals the eQTL lead"),
               ("eqtl_top_is_instrument", "eQTL lead is the MR instrument")]:
    if c in cy.columns:
        say("  %-34s %s" % (lab, pct(int(cy[c].sum()), len(cy))))
ind = pd.read_csv(P("02_coloc/h4_independent_loci.csv"))
say("pairs still > 0.8 in own-locus windows: %d (%s)" % (len(ind), ", ".join(ind.gene)))
say("")
say(open(P("02_coloc/coloc_p12_summary.txt"), encoding="utf-8").read()
    .split("retention of the 71 primary pairs")[0].strip())

# -------------------------------------------------- 5 tissue transferability
title("5  Tissue transferability of the instruments")
ts = pd.read_csv(P("04_tissue/thyroid_3source_validation.csv"))
say("instruments (tier A/B genes) : %d" % len(ts))
for c, lab in [("s1", "Zenodo salmon"), ("s2", "eQTL Catalogue thyroid"), ("s3", "GTEx portal API")]:
    rec = int(ts[c + "_p"].notna().sum())
    sg = int(ts[c + "_sig"].sum())
    same = int(ts.loc[ts[c + "_p"].notna(), c + "_same"].sum())
    same_sig = int(ts.loc[ts[c + "_sig"] == 1, c + "_same"].sum())
    say("  %-24s matchable %2d/%d ; significant %2d (%s of 70; %s of matchable) ; same sign %s of matchable, %s of significant"
        % (lab, rec, len(ts), sg, f(100.0 * sg / len(ts), 1),
           f(100.0 * sg / rec, 1) if rec else "n/a", pct(same, rec), pct(same_sig, sg)))
allsrc = ts[ts.s1_sig & ts.s2_sig & ts.s3_sig]
say("significant in all three     : %d (%s)" % (len(allsrc), ", ".join(allsrc.gene)))
for _, row in allsrc.iterrows():
    say("  %-8s same sign s1/s2/s3 = %s/%s/%s"
        % (row.gene, row.s1_same, row.s2_same, row.s3_same))
gt = pd.read_csv(P("04_tissue/gtex_tissue_rates.csv"))
say("")
say("GTEx v8 tissues scanned      : %d (denominator 70 in every tissue)" % len(gt))
for t in ["Whole_Blood", "Thyroid"]:
    row = gt[gt.tissue == t].iloc[0]
    say("  %-12s significant %s ; same sign %d" % (t, pct(int(row.n_sig), 70), int(row.n_same)))
say("  median tissue detection    : %s" % pct(int(gt.n_sig.median()), 70))
thy = gt.loc[gt.tissue == "Thyroid", "detect_pct"].iloc[0]
say("  thyroid rank               : %d of %d" % (int((gt.detect_pct > thy).sum()) + 1, len(gt)))

# --------------------------------------------------------------- 6 SCRN2
title("6  SCRN2: two instrument sets, four outcomes")
say("[A] whole-blood eQTLGen instrument rs11079797 (Wald ratio)")
s = mr[mr.gene == "SCRN2"]
say(s[["outcome", "method", "nsnp", "b", "se", "pval", "fdr"]].to_string(index=False))
say("")
say("[B] thyroid GTEx v8 instrument rs6503875 (Wald ratio)")
gw = pd.read_csv(P("08_plane/gw_thyroid_mr.csv"), encoding="utf-8-sig")
say(gw[gw.gene == "SCRN2"][["outcome", "snp", "F", "b", "se", "z", "p", "fdr"]].to_string(index=False))
say("")
say("[C] region-wide colocalisation (all pairs, default priors)")
say(co[co.gene == "SCRN2"][["outcome", "PP.H3", "PP.H4"]].to_string(index=False))
say("")
say("[D] window-restricted colocalisation (04_data/05_locus/24_coloc_by_signal_window.csv)")
wc = pd.read_csv(P("05_locus/24_coloc_by_signal_window.csv"), encoding="utf-8-sig")
say(wc[wc.gene == "SCRN2"][["tissue", "outcome", "window", "nsnp", "H3", "H4"]].to_string(index=False))

# ------------------------------------------------------- 7 locus dissection
title("7  Locus dissection and conditional models")
at = pd.read_csv(P("05_locus/23_signal_attribution_summary.csv"), encoding="utf-8-sig")
cols = ["signal", "sentinel", "candidate_gene", "tissue", "F", "coloc_H3", "coloc_H4",
        "window_coloc_H3", "window_coloc_H4", "verdict"]
say("signal attribution table (%d rows x %d columns)" % at.shape)
say(at[[c for c in cols if c in at.columns]].to_string(index=False))
say("")
cp = pd.read_csv(P("05_locus/31_joint_conditional_model.csv"), encoding="utf-8-sig")
q = cp[(cp["outcome"].isin(["TSH", "HYPO"])) & (cp.SNP.isin(["rs6503875", "shoulder"]))]
say("conditional z at rs6503875 and at the 45.873 Mb shoulder:")
say(q[["outcome", "model", "SNP", "pos37", "marginal_z", "marginal_P", "joint_z", "joint_P", "r_with_first_snp"]].to_string(index=False))
say("")
say("reading: in the five-variant model the hypothyroid estimate at rs6503875 is absorbed")
say("         (z -5.70 -> +0.16); in the four-variant model that drops rs62076102 a residual")
say("         z of -2.98 remains. The claim is 'explained by the neighbouring signal',")
say("         not 'no association with hypothyroidism'.")

# ---------------------------------------------------- 8 axis and ancestry
title("8  TSH-free T4 plane and cross-ancestry comparison")
plane = pd.read_csv(P("08_plane/s08_plane_genes.csv"))
say("genes with both estimates (thyroid instrument): %d" % len(plane))
lr = stats.linregress(plane.b_TSH, plane.b_T4)
say("Pearson r = %s ; OLS slope = %s"
    % (f(stats.pearsonr(plane.b_TSH, plane.b_T4).statistic, 3), f(lr.slope, 3)))
say("median per-gene ratio r = beta_FT4/beta_TSH : %s" % f(plane.r.median(), 2))
say("median ratio among TSH-significant genes    : %s"
    % f(plane.loc[plane.tsh_sig == 1, "r"].median(), 2))
say("TSH-significant genes with a T4 interval covering 0: %s"
    % pct(int(plane.loc[plane.tsh_sig == 1, "t4_null"].sum()),
          int((plane.tsh_sig == 1).sum()), 0))
tpo = plane[plane.gene == "TPO"]
if len(tpo):
    say("TPO ratio (glandular reference)             : %s" % f(tpo.r.iloc[0], 2))
if "SCRN2" in set(plane.gene):
    say("SCRN2 ratio                                 : %s" % f(plane.loc[plane.gene == "SCRN2", "r"].iloc[0], 2))
say("")
say("power-matched subset (04_data/08_plane/s08b_power.txt):")
say(open(P("08_plane/s08b_power.txt"), encoding="utf-8").read().strip())
say("")
say("cross-ancestry (04_data/08_plane/x02_transancestry.txt):")
# first block of the file only: pooling and the wording of the conclusion.
# The Bayesian section below it is reported by its own numbers in the manuscript and
# is not repeated here, to keep this section to the point estimates being verified.
for line in open(P("08_plane/x02_transancestry.txt"), encoding="utf-8"):
    s = line.strip()
    if s.startswith("====="):
        continue
    if s.startswith("Bayesian weight"):
        break
    if s:
        say("  " + s)

# ------------------------------------------------------------ 9 expression
title("9  Expression, single cell and co-expression")
a2 = pd.read_csv(P("06_expression/a2_summary.csv"))
say("[E] cross-tissue eQTL for the thyroid instrument (47 of the 49 GTEx v8 tissues with cis-eQTL data)")
say(a2[a2.gene == "SCRN2"].to_string(index=False))
say("    per-tissue source: 04_data/06_expression/SCRN2_rs6503875_GTEx_tissues.csv")
tk = pd.read_csv(P("06_expression/SCRN2_rs6503875_GTEx_tissues.csv"))
say("    tissues returned %d ; significant %d ; same sign %d of the other %d"
    % (len(tk), int((tk.p < 0.05).sum()), int((tk.nes < 0).sum()) - 1, len(tk) - 1))
pr = pd.read_csv(P("06_expression/genes_expression_profile.csv"))
say(pr[pr.gene == "SCRN2"][["thyroid_tpm", "blood_tpm", "thyroid_blood_ratio",
                            "thyroid_rank", "n_tissues", "top_tissue", "tau"]].to_string(index=False))
ct = pd.read_csv(P("06_expression/SCRN2_thyroid_celltype.csv"))
ct = ct[ct.n_cells >= 100].sort_values("pct_expressing", ascending=False)
say("")
say("single cell (cell types with >=100 cells):")
say(ct.head(6).to_string(index=False))
co_ex = pd.read_csv(P("06_expression/s07_panel_pc10.csv"))
say("")
say("partial Spearman rho after 10 expression PCs:")
say(co_ex[["gene", "module", "rho_pc10"]].to_string(index=False))

# -------------------------------------------------------------- 10 exposure
title("10 Exposure arm (Supplementary Table S1)")
s1 = pd.read_csv(os.path.join(SUP, "Table_S3_EDC_exposure_arm_MR.csv"))
say("exposures instrumented : %d ; testable pairs : %d" % (s1["EDC exposure"].nunique(), len(s1)))
say("instrument entries summed over pairs : %d  (each variant is reused across outcomes)"
    % int(s1["n SNP"].sum()))
say("distinct variants     : %d  (one to three per exposure; see the Instruments sheet)"
    % int(s1.groupby("EDC exposure")["n SNP"].max().sum()))
say("nominally significant : %d of %d ; smallest FDR q = %s"
    % (int((s1.P < 0.05).sum()), len(s1), f(s1["FDR q"].min(), 3)))
say(s1[s1.P < 0.05][["EDC exposure", "Outcome", "n SNP", "beta", "SE", "P", "FDR q"]].to_string(index=False))
say("")
h = pd.read_csv(P("06_expression/SCRN2_CTD_EDC_human.csv"))
say("curated human EDC records for SCRN2 : %d" % len(h))
say("  by class     : %s" % h.edc_class.value_counts().to_dict())
say("  by molecule  : %s" % h.molecule.value_counts().to_dict())
say("  by direction : %s" % h.direction.value_counts().to_dict())
say("  mRNA         : %s" % h[h.molecule == "mRNA"].direction.value_counts().to_dict())
say("  protein      : %s" % h[h.molecule == "protein"].direction.value_counts().to_dict())

# ------------------------------------------------------------- 11 structure
title("11 Structural arm: docking and MD")
dp = pd.read_csv(P("07_structure/dock_PFNA.csv"))
say("best blind-docking score : PFNA %s kcal/mol; mean of top ten %s kcal/mol"
    % (f(dp.energy.iloc[0]), f(dp.energy.iloc[:10].mean())))
C = dp[["cx", "cy", "cz"]].values[:10]
_dbest = np.linalg.norm(C - C[0], axis=1)[1:]
iu = np.triu_indices(len(C), 1)
_dpair = np.linalg.norm(C[iu[0]] - C[iu[1]], axis=1)
say("pose convergence        : all ten centroids within %s A of the best pose;"
    " widest separation between any two %s A"
    % (f(_dbest.max(), 2), f(max(_dpair), 2)))
say("")
say(pd.read_csv(P("07_structure/summary_table.csv")).to_string(index=False))
say("")
say("per-target minimum distance over the 10-ns trajectory (1001 frames, heavy atoms, 5 A):")
say(pd.read_csv(P("07_structure/md_contacts_recomputed.csv")).to_string(index=False))
say("")
say("per-residue contact frequency (<5 A):")
say(pd.read_csv(P("07_structure/md_per_residue_contacts.csv")).to_string(index=False))
mo = json.load(open(P("07_structure/motif_occlusion.json")))
say("")
say("motif occlusion in the docked pose (total buried %.1f A^2):" % mo.get("total_buried", 0))
for k, v in sorted(mo.get("motifs", {}).items()):
    say("  %-3s %s (%d-%d) buried %6.1f A^2 of %6.1f = %5.1f%% ; min distance %.1f A"
        % (k, v.get("seq"), v["range"][0], v["range"][1], v.get("buried", 0), v.get("apo", 0),
           100 * v.get("frac", 0), v.get("dmin", float("nan"))))

# ------------------------------------------------------------ 12 not derived
title("12 Values not recomputed by this script")
say("Each entry names where the value comes from, so it can be checked at source.")
say("  power-matched counts (Figure 7D)        : 04_data/08_plane/s08b_power.txt")
say("  GEO log2 fold changes (Figure 8E)       : 04_data/06_expression/SCRN2_geo_expression_validation.txt")
say("  threshold sensitivity (Figure 5D)       : recomputed offline by 05_scripts/90_figures/42_five_layer_sensitivity.py")
say("  South Asian and second European cohort  : external GWAS, accession in Table 1")
say("  MD trajectory-derived metrics           : 04_data/07_structure/ (trajectories not redistributed)")

title("End of report")
say("Generated by verify_numbers.py. If a number here disagrees with the manuscript,")
say("this report is the one computed from the deposited tables - correct the manuscript.")

with io.open(os.path.join(HERE, "verification_report.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(LINES) + "\n")
print("\n[written] 05_scripts/verification_report.txt")
