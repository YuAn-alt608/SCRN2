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
  10_screen … 90_figures/     numbered by analysis block, mirroring 04_data
  verify_numbers.py      regenerates every reported percentage and prints numerator/denominator
  verification_report.txt  output of verify_numbers.py
  recompute_md_contacts.py  per-residue contact frequencies from the MD trajectories

06_supplementary/             Supplementary Tables S1 (exposure-arm MR) and S2 (five layers); CSV + XLSX

07_structure_models/             structures and simulation settings for the structural arm
  pdb/                   AlphaFold2 SCRN2 model, the PFNA docked complex, ligand poses
  mdp/                   em / nvt / npt / production .mdp; raw/ holds the mdout.mdp dump
```

## Supplementary material

The manuscript cites two numbered supplementary items, and both are here:

- **Supplementary Table S1** — per-pair MR estimates for the nine instrumentable
  EDC exposures (section 2.10):
  `06_supplementary/Table_S1_EDC_exposure_arm_MR.{csv,xlsx}`
- **Supplementary Table S2** — layer membership and threshold sensitivity for the
  five evidence layers (Figure 5D, Table 3C, section 4.8):
  `06_supplementary/Table_S2_five_layer_membership_and_sensitivity.xlsx`

Nothing in the text cites an S3 or a Figure S1, so the numbering matches
one-to-one. Column definitions are in the README in that folder.

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

Two other entry points are self-contained in the same way:
`05_scripts/60_expression/39_gtex_cross_tissue_eqtl.py` regenerates the per-tissue table behind
section 2.9 and Figure 8B (47 tissues, 49 with cis-eQTL data), and
`05_scripts/recompute_md_contacts.py` recomputes per-residue contacts from a trajectory.

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
gene lists in Supplementary Table S2 and the regeneration script in
`05_scripts/90_figures/42_five_layer_sensitivity.py`; nothing in that figure is entered by hand.

Two panels still carry values entered in the plotting scripts rather than derived from a
deposited table — the power-matched counts (Figure 7D) and the GEO fold changes (Figure 8E).
They are listed here as the remaining gap between the figures and the deposited tables.
