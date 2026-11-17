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
04_数据/                 derived result tables, one folder per analysis block
  01_筛查/               EDC–gene universe, instruments, MR screen, evidence tiers
  02_共定位/             coloc.abf posteriors, circularity audit, independent-locus reruns
  03_复现/               FinnGen R12 consistency check and negative controls
  04_组织/               GTEx tissue panel and the three thyroid eQTL processings
  05_位点/               locus dissection: LD panel, joint GLS models, window colocalisation
  06_表达与单细胞/       expression, single-cell, co-expression, CTD records, GEO checks
  07_结构/               docking scores, MD metrics, per-residue contacts, MM-GBSA, pockets
  08_平面与跨族裔/       TSH–FT4 plane, power-matched subset, cross-ancestry meta-analysis
  09_五层/               five-layer region counts and bitmask combinations

05_脚本/                 analysis code
  10_筛查 … 90_图表/     numbered by analysis block, mirroring 04_数据
  verify_numbers.py      regenerates every reported percentage and prints numerator/denominator
  核对报告.txt           output of verify_numbers.py
  recompute_md_contacts.py  per-residue contact frequencies from the MD trajectories

06_补充材料/             Supplementary Table S1 (exposure-arm MR, 32 pairs; CSV + XLSX)

07_结构模型/             structures and simulation settings for the structural arm
  pdb/                   AlphaFold2 SCRN2 model, PFNA and T3 docked complexes, ligand poses
  mdp/                   em / nvt / npt / production .mdp; raw/ holds the two mdout.mdp dumps
```

## Supplementary material

The manuscript cites one numbered supplementary item, **Supplementary Table S1**
(section 2.10, per-pair estimates for the nine instrumentable EDC exposures). It
is `06_补充材料/Table_S1_EDC_exposure_arm_MR.{csv,xlsx}`; see the README in that
folder for the column definitions. Nothing in the text cites an S2 or a Figure S1,
so the numbering matches one-to-one.

## Reproducing the reported numbers

`05_脚本/verify_numbers.py` reads the tables in `04_数据/` and prints, for each reported
statistic, the numerator and the denominator. Running it is the fastest way to check any
number in the manuscript:

```bash
python 05_脚本/verify_numbers.py > 05_脚本/核对报告.txt
```

Scripts are numbered in the order they were run within each block. Paths are resolved
relative to the repository root, so the repository can be cloned anywhere.

## Environment

- R 4.4.2 — coloc 5.2.3, TwoSampleMR 0.6.x, susieR 0.12.x
- Python 3.11 — pandas 2.x, numpy 1.26, scipy 1.11, MDAnalysis 2.7, RDKit, PyMOL (for figures)
- GROMACS 2022.5 for the molecular dynamics runs in `07_结构/`

## External data sources

All input data are public and are identified by accession in Table 1 of the manuscript:
CTD (curated human chemical–gene interactions), eQTLGen, GTEx v8, eQTL Catalogue
(QTD000341), FinnGen R12 (E4_HYTHY_AI_STRICT), GWAS Catalog (GCST90572789, GCST90572790,
GCST90572791, GCST90018855, GCST90448634, GCST90319321) and 1000 Genomes phase 3.

Two figures depend on values that are entered in the plotting scripts rather than derived
from a deposited table — the threshold-sensitivity result (Figure 5D), the power-matched
counts (Figure 7D) and the GEO fold changes (Figure 8E). These are labelled
"author-reported" in the manuscript and in the figure legends.
