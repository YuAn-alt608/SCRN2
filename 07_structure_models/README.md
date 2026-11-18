# 07_structure_models — docked complexes and GROMACS parameter files

Geometry and simulation settings for the structural arm of the paper (Supplementary
Table S3, Figure 9, and the last Methods paragraph). The numbers these files produce —
docking scores, contact counts, residence fractions, diffusion, MM-GBSA — are
in `04_data/07_structure/`; this folder holds the structures themselves and the
run settings needed to reproduce them.

## `pdb/`

| File | Content | Heavy atoms |
|---|---|---|
| `AF-SCRN2-Q96FV2.pdb` | AlphaFold2 model of human SCRN2 (UniProt Q96FV2), chain A, residues 1–425 | 3,277 |
| `PFNA_complex.pdb` | **The docked complex**: same receptor + PFNA as chain B (residue `PFN`, residue 1). Geometry used for the PLIP contact analysis and as the MD starting structure | 3,277 + 29 |
| `PFNA_docked_pose.pdb` | Ligand-only coordinates of the top-scoring blind-docking pose (no receptor) | 29 |
| `PFNA_ligand.pdb` | Input ligand geometry (RDKit ETKDGv3 → MMFF94, AM1-BCC charges, GAFF2 typing) | 28 |

Notes on the coordinates:

- No `CRYST1` record: these are in the docking frame, without a simulation box.
  The box and solvation of the MD system are not stored here; the MD system was
  53,696 atoms (6,489 solute: 6,461 protein, 28 ligand) in a triclinic box, as
  reported in the Methods.
- Hydrogens are absent from the PDB files: they were added by GROMACS
  (`pdb2gmx`, amber99sb-ildn) at system build time, which is why the protein
  contributes more atoms to the simulation than the 3,277 heavy atoms in
  `AF-SCRN2-Q96FV2.pdb`. For the same reason the ligand enters the simulation as
  the 28-atom carboxylate of `PFNA_ligand.pdb`, not as the 29-atom neutral acid
  of the docking pose.
- The two ligand counts are not a contradiction: `PFNA_docked_pose.pdb` carries
  the carboxylic proton (29 atoms) because the blind docking was run on the
  neutral molecule, and the AM1-BCC/GAFF2 parameterisation used for MD was run
  on the deprotonated carboxylate (28 atoms, net charge −1), which is what
  Supplementary Table S3 reports.
- Residue numbering is that of the AlphaFold2 model (1–425), not of a
  crystallographic reference.

## `mdp/`

| File | Stage | Length |
|---|---|---|
| `em.mdp` | Energy minimisation (steepest descent) | 20,000 steps |
| `nvt.mdp` | NVT equilibration, position-restrained | 100 ps (2 fs × 50,000) |
| `npt.mdp` | NPT equilibration, position-restrained | 200 ps (2 fs × 100,000) |
| `md_production.mdp` | Production, unrestrained | 10 ns (2 fs × 5,000,000) |

All stages: GROMACS 2022.5, amber99sb-ildn + TIP3P, PME, Verlet cutoff scheme,
h-bonds constrained with LINCS, V-rescale thermostat at **310 K**, C-rescale
barostat at 1 bar, coupling groups `Solute` / `Solvent`.

`mdp/raw/` holds the `mdout.mdp` file that GROMACS writes at start-up: the
complete, verbatim parameter set of the PFNA production run, which can be used to
check any line of the condensed `md_production.mdp`.

### Provenance, and one superseded template

The condensed `.mdp` files were rebuilt from the run logs (`em.log`, `nvt.log`,
`npt.log` and the `mdout.mdp` dump), because the original `.mdp` files were
not kept in the working archive. Each file says which log it came from; every
value is one that the log records, and unlisted parameters stayed at GROMACS
defaults.

The docking archive also contains `md_setup/md.mdp`, a **template that was not
used for the reported trajectories**: it specifies 50,000,000 steps (100 ns),
Nose–Hoover and Parrinello–Rahman coupling at 300 K. The runs behind Supplementary Table S3
used 5,000,000 steps (10 ns) with V-rescale and C-rescale at 310 K, as stated in
the Methods and as recorded in `mdp/raw/`. Do not use that template to reproduce
the paper.

## Reproducing

```bash
# build the system (protein + ligand + TIP3P + 0.15 M NaCl)
gmx pdb2gmx -f pdb/PFNA_complex.pdb -ff amber99sb-ildn -water tip3p -ignh
# then, in order
gmx grompp -f mdp/em.mdp  -c box.gro -p topol.top -o em.tpr   && gmx mdrun -deffnm em
gmx grompp -f mdp/nvt.mdp -c em.gro  -r em.gro -p topol.top -o nvt.tpr && gmx mdrun -deffnm nvt
gmx grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr && gmx mdrun -deffnm npt
gmx grompp -f mdp/md_production.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr && gmx mdrun -deffnm md
```

Contact, residence and diffusion statistics from the resulting trajectory are
computed by `05_scripts/recompute_md_contacts.py`.
