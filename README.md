# SCRN2 / thyroid set-point — analysis code and derived result tables

Code and derived result tables for the manuscript *"Colocalisation-resolved expression-QTL
screening of endocrine-disrupting chemical-responsive genes nominates SCRN2 as a candidate
thyroid set-point gene"*.

This repository contains **no individual-level data** and **no raw GWAS or eQTL summary
statistics**. It contains the analysis scripts and every derived table from which the
numbers in the manuscript are taken.

## Author

Zizheng Yuan — University of Lisbon, Lisbon, Portugal

## Citation

If you use this code or these tables, please cite the manuscript above, or this repository:

```bibtex
@software{yuan_scrn2_2026,
  author = {Yuan, Zizheng},
  title  = {SCRN2 / thyroid set-point: analysis code and derived result tables},
  year   = {2026},
  url    = {https://github.com/YuAn-alt608/SCRN2}
}
```

## Layout

```
04_data/                 derived result tables, one folder per analysis block
  01_screen/               EDC–gene universe, instruments, MR screen, evidence tiers
  02_coloc/             coloc.abf posteriors, circularity audit, independent-locus reruns
  03_replication/               FinnGen R12 consistency check and negative controls
  04_tissue/               GTEx tissue panel and the three thyroid eQTL processings
  05_locus/               locus dissection: LD panel, joint GLS models, window colocalisation
  06_expression/       expression, cross-tissue eQTL, single-cell, co-expression, CTD records, GEO checks
                         SCRN2_rs6503875_GTEx_tissues.csv = the 47-tissue lookup behind section 2.9
                         cross_tissue_coherence_panel.csv = the 71-gene panel behind layer 4 of Figure 5
  07_structure/               docking scores, MD metrics, per-residue contacts, MM-GBSA, pockets
  08_plane/       TSH–FT4 plane, power-matched subset, cross-ancestry meta-analysis
  09_layers/               five-layer membership, region counts and threshold sensitivity

05_scripts/                 analysis code
  00_upstream/              three scripts that build the upstream MR tables from the raw GWAS
                            downloads, plus their own README explaining why they cannot be
                            re-run from this repository alone
  10_screen … 90_figures/     numbered by analysis block, mirroring 04_data
  verify_numbers.py      regenerates every reported percentage and prints numerator/denominator
  verification_report.txt  output of verify_numbers.py
  recompute_md_contacts.py  per-residue contact frequencies from the MD trajectories

06_supplementary/             Supplementary tables S1-S3: one typeset docx, plus S1/S2 as CSV + XLSX

07_structure_models/             structures and simulation settings for the structural arm
  pdb/                   AlphaFold2 SCRN2 model, the PFNA docked complex, ligand poses
  mdp/                   em / nvt / npt / production .mdp; raw/ holds the mdout.mdp dump
```

## Supplementary material

The manuscript cites three numbered supplementary items, all typeset in one
submission file, `06_supplementary/Supplementary_Material_SCRN2.docx`, built by
`06_supplementary/_build_supplementary.py`:

- **Supplementary Table S1** — layer membership and threshold sensitivity for the
  five evidence layers (Figure 5D, Table 3C, section 2.5 and 4.8):
  `06_supplementary/Table_S1_five_layer_membership_and_sensitivity.xlsx`, with the
  five sensitivity rows also in `Table_S1_sensitivity.csv`
- **Supplementary Table S2** — PFNA at the SCRN2 surface: docking, 10 ns all-atom
  MD and MM-GBSA (section 2.10). Sixteen numbers, assembled from
  `04_data/07_structure/summary_table.csv`, `dock_PFNA.csv` and `mmgbsa_PFNA.json`.
- **Supplementary Table S3** — per-pair MR estimates for the nine instrumentable
  EDC exposures (section 2.10):
  `06_supplementary/Table_S3_EDC_exposure_arm_MR.{csv,xlsx}`

The tables are numbered in the order the manuscript cites them (S1 in section
2.5, then S2 and S3 in section 2.10); v20 renumbered them, old S2 → S1, old S3 →
S2, old S1 → S3.

The docx is not in the repository (the package ships code and result tables, not
drafts); it is rebuilt in seconds from the two CSVs and the PFNA block. Column
definitions are in the README in that folder.

## Reproducing the reported numbers

`05_scripts/verify_numbers.py` reads the tables in `04_data/` and prints, for each reported
statistic, the numerator and the denominator. Running it is the fastest way to check any
number in the manuscript:

```bash
python 05_scripts/verify_numbers.py        # writes 05_scripts/verification_report.txt
```

It needs only pandas, numpy and scipy, reads nothing outside this repository, and does
not touch the network. If a number in the report disagrees with the manuscript, the
report is the one computed from the deposited tables.

Other entry points are self-contained in the same way:
`05_scripts/60_expression/39_gtex_cross_tissue_eqtl.py` regenerates the per-tissue table behind
section 2.9 and Figure 8B (47 tissues, 49 with cis-eQTL data),
`05_scripts/60_expression/40_gtex_axis_tissue.py` rebuilds `SCRN2_axis_expression.csv`, the
thyroid, pituitary and whole-blood comparison behind the set-point paragraph of section 2.9,
`05_scripts/80_figures/s08b_power_matched_plane.py` regenerates `04_data/08_plane/s08b_power.csv`
and the text next to it byte for byte from `s08_plane_genes.csv`, and
`05_scripts/recompute_md_contacts.py` recomputes per-residue contacts from a trajectory.

The three scripts in `05_scripts/00_upstream/` are the exception. Each one reads a full GWAS
summary-statistics file of several GB, which is not redistributed; they are deposited so the
derivation of `gw_thyroid_mr.csv`, `three_cohort_meta.csv` and `x02_transancestry.txt` is on
record. `x02_transancestry.py` needs no input at all and reproduces its deposited text exactly,
so it can be run here as a check.

The remaining scripts are numbered in the order they were run within each block and are
archived as the analysis record rather than as one-click pipelines: several of them read
the author's working directories (for example `D:/NewBeginning/MR1/DATA10`), because they
were written against the internal archive from which this release was assembled. Those
directories held the raw GWAS/eQTL downloads and are not redistributed; where a script
matters for a specific table, its `BASE` line is the only edit needed, and every value it
produced is already in `04_data/`.

## Environment

- R 4.4.2 — coloc 5.2.3, TwoSampleMR 0.6.x, susieR 0.12.x
- Python 3.11 — pandas 2.x, numpy 1.26, scipy 1.11, MDAnalysis 2.7, RDKit, PyMOL (for figures)
- GROMACS 2022.5 for the molecular dynamics runs in `07_structure/`

## External data sources

All input data are public and are identified by accession in Table 1 of the manuscript:
CTD (curated human chemical–gene interactions), eQTLGen, GTEx v8, eQTL Catalogue
(QTD000341), FinnGen R12 (E4_HYTHY_AI_STRICT), GWAS Catalog (GCST90572789, GCST90572790,
GCST90572791, GCST90018855, GCST90448634, GCST90319321) and 1000 Genomes phase 3.

Every number in Figure 5 is derived from `04_data/09_layers/`, with the per-configuration
gene lists in Supplementary Table S1 and the regeneration script in
`05_scripts/90_figures/42_five_layer_sensitivity.py`; nothing in that figure is entered by hand.

Every panel of Figure 7 is read from a deposited table: panel A from `s08_plane_genes.csv`,
panel B from the MR tables, panel C from the same plane table, and panel D from
`s08b_power.csv`, which `05_scripts/80_figures/s08b_power_matched_plane.py` regenerates from
the plane table. Figure 8 is drawn by `05_scripts/90_figures/44_fig8_tissue_context.py` from
`04_data/06_expression/`; panel E parses the summary text that the GEO re-analysis left behind
(`SCRN2_geo_expression_validation.txt`) rather than recomputing it from the GEO series matrix,
which is public but not redistributed here.

One analysis arm that appeared in earlier drafts is not part of this version: the T3 docking and
MD series was dropped, and Supplementary Table S3 reports PFNA only. S3 is typeset with S1 and
S2 in the one supplementary file, so the manuscript itself carries no appendix and five tables.
