# 00_upstream — the three generators that were not in the earlier deposit

Three tables in `04_data/` were produced before this release was assembled, by
scripts that lived in the working archives rather than here. They are included
because without them those tables cannot be traced back to code:

| Script | Writes | Read by |
|---|---|---|
| `t09_genomewide_thyroid_mr.py` | `04_data/08_plane/gw_thyroid_mr.csv` | `41_cross_tissue_coherence_panel.py`, `42_five_layer_sensitivity.py`, `90_figures/make_figures.py`, Table 5A |
| `s17_second_cohort.py` | `04_data/08_plane/three_cohort_meta.csv` | `90_figures/make_figures.py`, Table 5B |
| `x02_transancestry.py` | `04_data/08_plane/x02_transancestry.txt`, `.csv` | `verify_numbers.py`, section 2.8 |

None of the three can be re-run from this repository alone: they stream
multi-gigabyte GWAS files (the four thyroid outcome files, GCST90319321 and
GCST90448634) that are not redistributed here. They are included as a record of
how the deposited tables were made, not as a pipeline you can execute. Their
absolute input paths point at the working archive.

`x02_transancestry.py` is the one exception to "cannot be checked": it reads
nothing, and running it reproduces `04_data/08_plane/x02_transancestry.txt`
byte for byte. That is the cheapest way to confirm that the deposited summary
and the code behind it still agree.

Everything else in the deposit — the screen, colocalisation, replication,
tissue, locus-resolved, expression and structure analyses — is in the numbered
directories of `05_scripts/` and runs from `04_data/` alone.
