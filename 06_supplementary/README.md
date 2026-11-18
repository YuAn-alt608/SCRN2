# 06_supplementary — Supplementary Tables S1–S3

The manuscript cites three numbered supplementary items, and all three are
typeset in one submission file:

- **`Supplementary_Material_SCRN2.docx`** (+ `.pdf`) — the three tables, numbered
  in the order the manuscript cites them, one table per page. A referee reads
  S1 (the layer-threshold table, cited in section 2.5), S2 (the PFNA simulation
  summary, cited in section 2.10) and S3 (32 rows of exposure-arm MR estimates,
  cited in section 2.10) without opening three attachments. Built by
  `_build_supplementary.py`.
- **Supplementary Table S1** — five-layer membership and threshold sensitivity
  (Figure 5D and Table 3C), cited in section 2.5 and Methods 4.8.
- **Supplementary Table S2** — PFNA at the SCRN2 surface: docking, 10 ns all-atom
  MD and MM-GBSA (section 2.10). Sixteen numbers; the underlying rows are
  `04_data/07_structure/summary_table.csv`, `dock_PFNA.csv` and
  `mmgbsa_PFNA.json`.
- **Supplementary Table S3** — the exposure-arm MR estimates (section 2.10).

The numbering changed in v20. Until v19 the three tables were S1 = exposure-arm
MR, S2 = layers, S3 = PFNA, which is not the order in which the manuscript cites
them; S2 is cited first (section 2.5), then S3 (section 2.10), then S1 (section
2.10). v20 renumbered them old S2 → S1, old S3 → S2, old S1 → S3, and the
machine-readable files were renamed with them. The PFNA table was typeset in a
manuscript appendix up to v15 and moved here in v16, so that the manuscript
carries five tables and one figure set and the supplement carries the three
supplementary tables. The docx is built on the manuscript file, which is why it
inherits the same stylesheet and page numbering.

## Files

| File | Use |
|---|---|
| `Supplementary_Material_SCRN2.docx` | Submission file: S1, S2 and S3 in one document (4 pages) |
| `Supplementary_Material_SCRN2.pdf` | Same, exported for review |
| `_build_supplementary.py` | Builds the docx from the two CSVs and the PFNA block of the manuscript |
| `Table_S1_five_layer_membership_and_sensitivity.xlsx` | Workbook version of S1: notes, sensitivity table, per-gene membership in every configuration, and the 71-gene coherence panel |
| `Table_S1_sensitivity.csv` | Machine-readable sensitivity table (the numbers in Figure 5D) |
| `Table_S3_EDC_exposure_arm_MR.csv` | Machine-readable version of S3; this is what `verify_numbers.py` reads |
| `Table_S3_EDC_exposure_arm_MR.xlsx` | Workbook version: title, header, notes, and a second sheet with the per-exposure instrument counts |

The two S3 files contain the same 32 rows. The docx prints S3 to three decimals
and S1 to the precision of the CSV; the spreadsheets keep the six-significant-digit
values, so round-off can differ in the last digit of `MDE80` only.

## Table S3 — what it contains

One row per testable exposure–outcome pair (9 EDC exposures × 4 outcomes, of
which 32 pairs had at least one instrument):

| Column | Meaning |
|---|---|
| EDC exposure | MECPP, MEHHP, MEHP, MEOHP, MMP, MethylParaben, PFOA, PFOS, PropylParaben |
| Outcome | TSH, Free T4, hypothyroidism, Hashimoto thyroiditis |
| n SNP | Number of instruments for that pair |
| Method | Inverse-variance weighted (11 pairs) or Wald ratio (21 pairs) |
| beta, SE, P | MR estimate, standard error, P value |
| FDR q | Benjamini–Hochberg across the 32 pairs |
| Cochran Q, P (Q), I² | Between-instrument heterogeneity, where ≥2 instruments |
| MDE80 | Smallest effect detectable at 80% power given the observed instrument strength |
| beta (WM), SE (WM), P (WM) | Weighted-median sensitivity analysis |

Three pairs reach nominal P < 0.05 — MECPP–hypothyroidism (β = 0.162,
P = 0.0079), MEHHP–hypothyroidism (β = −0.539, P = 0.027) and PFOS–FT4
(β = −0.077, P = 0.013) — and none survives FDR correction (minimum q = 0.21).
The PFOS–FT4 pair is the one whose direction matches the toxicological
literature. These are the numbers quoted in section 2.10 and in the Discussion.

### Counting the instruments

The 32 pairs come from 13 variants, because most exposures were instrumented by
one to three variants and each variant is reused across the four outcomes:

| Exposure | Testable pairs | Distinct variants |
|---|---|---|
| MECPP | 4 | 3 |
| MEHHP | 4 | 1 |
| MEHP | 4 | 1 |
| MEOHP | 4 | 1 |
| MMP | 1 | 1 |
| MethylParaben | 3 | 1 |
| PFOA | 4 | 1 |
| PFOS | 4 | 2 |
| PropylParaben | 4 | 2 |
| **Total** | **32** | **13** |

Summing `n SNP` over the 32 pairs gives 47, the figure the sheet note reports.
The two numbers count different things and both are correct; the Methods of the
manuscript uses the 13-variant count for the instrument set and this table for
the per-pair estimates.

## Table S1 — what it contains

Five layers on the TSH branch: whole-blood MR (337 genes), thyroid MR (180),
colocalisation (26), cross-tissue coherence (36) and target assignment (4).
Membership is encoded as a bitmask, `bitmask = L1 + 2*L2 + 4*L3 + 8*L4 + 16*L5`.

| Sheet | Content |
|---|---|
| `Notes` | What the table is and how the bitmask is defined |
| `Sensitivity` | Five configurations (baseline, PP.H4 ≥ 0.70, targets top-5, coherence ≥ 90%, all three) with the size of each layer and the five-layer intersection |
| `Membership` | Every gene that appears in at least one layer, with its bitmask in each configuration |
| `Coherence panel` | The 71 genes passing layers 1 and 2, with the variant used, the number of testable tissues, the concordance and the layer-4 call |

Layer 4 is stored as a concordance of 1.00 — every tissue in which the instrument is a
significant cis-eQTL must agree in sign. A gene has at most 48 such tissues (49 GTEx v8
tissues with cis-eQTL data, minus thyroid), so 0.99 and 1.00 select the same genes; the
earlier 0.99 value was a threshold that could not bind, and the relaxed configuration uses
0.90, which admits one discordant tissue in ten.

In every configuration the five-layer intersection contains one gene, SCRN2.
The layer sizes are 337 / 180 / 26 / 36 / 4 at baseline; relaxing colocalisation
moves layer 3 to 31, relaxing the target criterion moves layer 5 to 7, and
relaxing direction concordance to 90% moves layer 4 to 43.

## Regenerating

```bash
python 05_scripts/60_expression/41_cross_tissue_coherence_panel.py   # coherence panel (GTEx API)
python 05_scripts/90_figures/42_five_layer_sensitivity.py          # membership + sensitivity
python 06_supplementary/_build_supplementary.py                    # the docx
python 06_supplementary/_export.py                                 # the pdf
```

`41_…` needs network access on first run and caches its GTEx responses; `42_…`
is offline and reads only `04_data/`. The S3 numbers are re-printed, with
numerators and denominators, by section 10 of `05_scripts/verify_numbers.py`.
