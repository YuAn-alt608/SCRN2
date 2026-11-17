#!/usr/bin/env python
"""01 - ligand preparation: 3D structure of PFNA + GAFF2/AM1-BCC parameters -> GROMACS topology

Charge state (a modelling decision that has to be stated):
  PFNA  perfluorononanoic acid, pKa ~0.5, fully deprotonated at physiological pH -> carboxylate, net charge -1

Only one ligand is simulated, so there is no comparator ligand and no charge-state mismatch to align.
"""
from pathlib import Path
import subprocess, os, sys, shutil

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

ROOT = Path(__file__).resolve().parents[1]
LIG = ROOT / "ligands"
AMBER = Path("D:/Tools/ChimeraX/bin/amber20")
BIN = AMBER / "bin"
# antechamber is a cygwin build that assembles an sh command line from AMBERHOME; backslashes are eaten by the shell,
# so paths must be passed with forward slashes or it dies with "Cannot properly run .../bondtype"
AMBERHOME_POSIX = AMBER.as_posix()

SPECIES = {
    # name: (SMILES, net charge, note)
    "PFNA": ("[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", -1,
             "perfluorononanoic acid, carboxylate"),
}


def build_pdb(name, smiles, seed=0xC0FFEE):
    m = Chem.MolFromSmiles(smiles)
    assert m is not None, name
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = seed
    ps.useRandomCoords = True
    assert AllChem.EmbedMolecule(m, ps) == 0, f"embed failed {name}"
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    p = LIG / f"{name}.pdb"
    Chem.MolToPDBFile(m, str(p))
    return p, Description_ok(m)


def Description_ok(m):
    from rdkit.Chem import rdMolDescriptors
    return dict(formula=rdMolDescriptors.CalcMolFormula(m),
                mw=round(Descriptors.MolWt(m), 2),
                rotb=rdMolDescriptors.CalcNumRotatableBonds(m),
                tpsa=round(Descriptors.TPSA(m), 1),
                logp=round(Descriptors.MolLogP(m), 2),
                natoms=m.GetNumAtoms())


def run(cmd, cwd):
    env = dict(os.environ)
    env["AMBERHOME"] = AMBERHOME_POSIX
    env["PATH"] = str(BIN) + os.pathsep + env.get("PATH", "")
    r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, shell=False)
    if r.returncode != 0:
        print(r.stdout[-3000:]); print(r.stderr[-3000:])
        raise SystemExit(f"FAILED: {' '.join(map(str, cmd))}")
    return r.stdout


def main():
    for name, (smi, nc, note) in SPECIES.items():
        pdb, info = build_pdb(name, smi)
        work = LIG / f"work_{name}"
        work.mkdir(exist_ok=True)
        shutil.copy(pdb, work / f"{name}.pdb")
        print(f"=== {name} ({note}) ===")
        print("   ", info, "net charge", nc)

        # antechamber: GAFF2 atom types + AM1-BCC charges (sqm semi-empirical; a 30-50 atom ligand takes minutes)
        run([BIN.joinpath("antechamber.exe").as_posix(), "-i", f"{name}.pdb", "-fi", "pdb",
             "-o", f"{name}.mol2", "-fo", "mol2", "-c", "bcc", "-nc", str(nc),
             "-s", "2", "-at", "gaff2", "-rn", name[:3].upper(), "-pf", "y"], work)
        # parmchk2: fills in the bond/angle/dihedral parameters GAFF2 does not carry
        run([BIN.joinpath("parmchk2.exe").as_posix(), "-i", f"{name}.mol2", "-f", "mol2",
             "-o", f"{name}.frcmod", "-s", "gaff2"], work)
        # acpype: mol2+frcmod -> GROMACS itp (charges are read from the mol2, not recomputed)
        run([sys.executable, "-m", "acpype", "-i", f"{name}.mol2", "-b", name,
             "-c", "user", "-n", str(nc), "-a", "gaff2", "-o", "gmx"], work)

        for f in (f"{name}.itp", f"{name}_GMX.itp", f"{name}_atomtypes.itp", f"{name}.top"):
            src = work / f
            if src.exists():
                shutil.copy(src, LIG / f)
        print("   ->", sorted(x.name for x in LIG.glob(f"{name}*")))


if __name__ == "__main__":
    main()
