# 06_补充材料 — Supplementary Table S1

The manuscript cites exactly one numbered supplementary item, **Supplementary
Table S1**, in the last paragraph of section 2.10:

> "Only nine EDC exposures could be instrumented at all, giving 13 distinct
> variants and 32 testable exposure–outcome pairs (47 instrument entries once
> each variant is reused across the four outcomes) … The per-pair estimates are
> given in Supplementary Table S1."

This folder holds that table, under exactly that number. No other supplementary
table or figure is cited anywhere in the text, so there is no S2 or Figure S1 to
match, and the numbering cannot collide.

## Files

| File | Use |
|---|---|
| `Table_S1_EDC_exposure_arm_MR.csv` | Machine-readable version; this is what `verify_numbers.py` reads |
| `Table_S1_EDC_exposure_arm_MR.xlsx` | Submission version: title, header, notes, and a second sheet with the per-exposure instrument counts |

Both contain the same 32 rows.

## What the table contains

One row per testable exposure–outcome pair (9 EDC exposures × 4 outcomes, of
which 32 pairs had at least one instrument):

| Column | Meaning |
|---|---|
| EDC exposure | MECPP, MEHHP, MEHP, MEOHP, MMP, MethylParaben, PFOA, PFOS, PropylParaben |
| Outcome | TSH, FT4, hypothyroidism, Hashimoto thyroiditis |
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

## Counting the instruments

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

## Regenerating

The table is a direct subset of the exposure-arm MR output (all rows with at
least one instrument). That intermediate file is not part of this release; the
published table is the authoritative record of the estimates. The block is
re-printed, with numerators and denominators, by section 10 of
`05_脚本/verify_numbers.py`.
