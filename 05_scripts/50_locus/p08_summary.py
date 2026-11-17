#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p08_summary.py
==============
Collapses everything into one decision table, one workbook and one heatmap:
  gene x tissue x outcome with lead SNP, eQTL P/F, cis-MR (Wald/IVW),
  coloc (H3/H4 at p12=1e-4), SMR, HEIDI and the effective number of instruments,
  plus the gene ranking under TSH.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
T = f"{BASE}/tables"

mr = pd.read_csv(f"{T}/04_cisMR_SMR_HEIDI.csv")
co = pd.read_csv(f"{T}/03_panel_coloc_abf_prior_grid.csv")
genes_in = pd.read_csv(f"{T}/01_gene_list.csv").gene_name.tolist()

co4 = co[co.p12 == 1e-4][["tissue", "gene", "outcome", "PP.H3", "PP.H4", "nsnp"]]
co4 = co4.rename(columns={"PP.H3": "coloc_H3", "PP.H4": "coloc_H4",
                          "nsnp": "coloc_nsnp"})
d = mr.merge(co4, on=["tissue", "gene", "outcome"], how="left")
d = d[d.gene.isin(genes_in)]
d["F_lead"] = d.F_lead.round(1)
d["heidi_pass"] = d.heidi_p > 0.01

# ranking score: TSH, requiring F>10
def score(r):
    s = 0.0
    if r.lead_p_eqtl < 5e-8: s += 1
    if r.p_wald < 0.05: s += 1
    if r.coloc_H4 > 0.8: s += 1
    if r.p_smr < 0.05 and (r.heidi_p > 0.01 if r.heidi_p == r.heidi_p else False): s += 1
    return s

d["evidence_score_0to4"] = d.apply(score, axis=1)
d.to_csv(f"{T}/15_decision_table_all_genes.csv", index=False, encoding="utf-8-sig")

tsh = d[d.outcome == "TSH"].sort_values(["evidence_score_0to4", "p_wald"], ascending=[False, True])
tsh.to_csv(f"{T}/15b_decision_table_TSH.csv", index=False, encoding="utf-8-sig")
cols = ["tissue", "gene", "lead_snp", "lead_p_eqtl", "F_lead", "b_wald", "se_wald", "p_wald",
        "n_iv", "k_eff", "coloc_H3", "coloc_H4", "p_smr", "heidi_p", "evidence_score_0to4"]
pd.set_option("display.width", 260, "display.max_columns", 40)
print(tsh[cols].round(4).to_string(index=False))

# heatmap: coloc H4
piv = (d[d.outcome == "TSH"].pivot_table(index="gene", columns="tissue",
                                         values="coloc_H4"))
piv = piv.reindex(piv.max(axis=1).sort_values(ascending=False).index)
fig, ax = plt.subplots(figsize=(4.6, 7))
im = ax.imshow(piv.values, cmap="RdYlBu_r", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=30, ha="right", fontsize=9)
ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index, fontsize=8)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v = piv.values[i, j]
        if v == v:
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if v > 0.6 or v < 0.15 else "black")
ax.set_title("coloc PP.H4 with TSH (p12=1e-4)", fontsize=10)
fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
fig.tight_layout()
fig.savefig(f"{BASE}/fig_panel_coloc_heatmap.png", dpi=200)
print("heatmap saved")

# workbook
with pd.ExcelWriter(f"{T}/../MR1_supplementary_results_2026-11-14.xlsx", engine="openpyxl") as xw:
    tsh[cols].to_excel(xw, sheet_name="TSH_decision", index=False)
    d.to_excel(xw, sheet_name="all_outcomes_decision", index=False)
    co.to_excel(xw, sheet_name="coloc_prior_grid", index=False)
    mr.to_excel(xw, sheet_name="cisMR_SMR_HEIDI", index=False)
    pd.read_csv(f"{T}/05_rs7219420_marginal.csv").to_excel(xw, sheet_name="rs7219420_marginal", index=False)
    pd.read_csv(f"{T}/06_two_snp_joint_model.csv").to_excel(xw, sheet_name="two_snp_joint", index=False)
    pd.read_csv(f"{T}/07_cis_eqtl_at_two_snps.csv").to_excel(xw, sheet_name="two_snp_eqtl", index=False)
    pd.read_csv(f"{T}/11_variant_info.csv").to_excel(xw, sheet_name="variant_info", index=False)
    pd.read_csv(f"{T}/13_locus_identity_regional_min_p.csv").to_excel(xw, sheet_name="locus_identity", index=False)
    pd.read_csv(f"{T}/14_three_cohort_sample_overlap.csv").to_excel(xw, sheet_name="sample_overlap", index=False)
    for f in ["12_gwas_catalog_locus_identity.csv"]:
        if os.path.exists(f"{T}/{f}"):
            pd.read_csv(f"{T}/{f}").to_excel(xw, sheet_name="gwas_catalog", index=False)
    if os.path.exists(f"{T}/16_variant_ccre_annotation.csv"):
        pd.read_csv(f"{T}/16_variant_ccre_annotation.csv").to_excel(xw, sheet_name="ccre_annotation", index=False)
    if os.path.exists(f"{T}/17_two_gamma_sensitivity.csv"):
        pd.read_csv(f"{T}/17_two_gamma_sensitivity.csv").to_excel(xw, sheet_name="two_gamma", index=False)
    if os.path.exists(f"{T}/18_instrument_count_gls_ivw.csv"):
        pd.read_csv(f"{T}/18_instrument_count_gls_ivw.csv").to_excel(xw, sheet_name="instrument_count_scan", index=False)
    # 2026-11-14 fix: the TSH-side SuSiE result is unreliable (it still reports an
    # excessive prior variance after allele alignment, and its lead credible-set
    # SNP contradicts the marginal P). The TSH versions of 03b/03c were removed and
    # only the eQTL side (the official QTD000341 credible_sets) is kept.
    if os.path.exists(f"{T}/03c_credible_set_eqtl_side.csv"):
        pd.read_csv(f"{T}/03c_credible_set_eqtl_side.csv").to_excel(xw, sheet_name="credible_set_eqtl_side", index=False)
    for f, sn in [("19_gls_sign_handling_comparison.csv", "gls_sign_handling"),
                  ("19b_three_instruments_per_snp_wald.csv", "per_snp_wald"),
                  ("19c_ld_target_allele.csv", "ld_target_allele"),
                  ("20_simulation_sign_gamma.csv", "gls_sign_simulation"),
                  ("21_published_sentinels_17q21.csv", "published_sentinels"),
                  ("22_sentinel_gene_mapping_uniqueness.csv", "sentinel_mapping_uniqueness"),
                  ("23_signal_attribution_summary.csv", "signal_attribution"),
                  ("24_coloc_by_signal_window.csv", "coloc_by_window"),
                  ("25_hypo_data_source_check.csv", "hypo_source_check"),
                  ("26_peak_neighbour_genes_window_coloc.csv", "neighbour_gene_coloc"),
                  ("27_shoulder_gene_attribution.csv", "shoulder_attribution"),
                  ("28_peak_variant_eqtl_scan.csv", "peak_variant_eqtl_scan"),
                  ("29_shoulder_locus_per_gene_wald.csv", "shoulder_per_gene_wald"),
                  ("30_shoulder_gtex_cross_tissue.csv", "shoulder_gtex_cross_tissue"),
                  ("31_joint_conditional_model.csv", "joint_conditional_model")]:
        p = f"{T}/{f}"
        if os.path.exists(p):
            pd.read_csv(p).to_excel(xw, sheet_name=sn, index=False)
    if os.path.exists(f"{T}/08_regional_plot_data.csv"):
        pd.read_csv(f"{T}/08_regional_plot_data.csv").to_excel(xw, sheet_name="regional_plot_data", index=False)
    for f, sn in [("09b_per_tissue_top_eqtl_p.csv", "per_tissue_eqtl"), ("10_per_tissue_key_snp_eqtl.csv", "key_snp_per_tissue")]:
        p = f"{T}/{f}"
        if os.path.exists(p):
            pd.read_csv(p).to_excel(xw, sheet_name=sn, index=False)
    for ws in xw.book.worksheets:
        ws.freeze_panes = "A2"
print("workbook saved")
