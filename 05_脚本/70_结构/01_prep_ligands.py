#!/usr/bin/env python
"""01 · 配体准备：PFNA 的 3D 结构 + GAFF2/AM1-BCC 参数 -> GROMACS 拓扑

电荷态选择（这是必须交代的建模决定）：
  PFNA  全氟壬酸，pKa ~0.5，生理 pH 下完全去质子化 -> 羧酸根，净电荷 -1

只模拟一个配体，因此没有对照配体，也没有需要对齐的电荷态差异。
"""
from pathlib import Path
import subprocess, os, sys, shutil

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

ROOT = Path(__file__).resolve().parents[1]
LIG = ROOT / "ligands"
AMBER = Path("D:/Tools/ChimeraX/bin/amber20")
BIN = AMBER / "bin"
# antechamber 是 cygwin 编译的，它会用 AMBERHOME 拼出 sh 命令行；反斜杠会被 shell 吃掉，
# 必须传正斜杠，否则报 "Cannot properly run .../bondtype"
AMBERHOME_POSIX = AMBER.as_posix()

SPECIES = {
    # name: (SMILES, net charge, 说明)
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

        # antechamber: GAFF2 原子类型 + AM1-BCC 电荷（sqm 半经验，30-50 原子几分钟内可完成）
        run([BIN.joinpath("antechamber.exe").as_posix(), "-i", f"{name}.pdb", "-fi", "pdb",
             "-o", f"{name}.mol2", "-fo", "mol2", "-c", "bcc", "-nc", str(nc),
             "-s", "2", "-at", "gaff2", "-rn", name[:3].upper(), "-pf", "y"], work)
        # parmchk2: 补齐 GAFF2 缺失的键/角/二面角参数
        run([BIN.joinpath("parmchk2.exe").as_posix(), "-i", f"{name}.mol2", "-f", "mol2",
             "-o", f"{name}.frcmod", "-s", "gaff2"], work)
        # acpype: mol2+frcmod -> GROMACS itp（电荷取自 mol2，不再重算）
        run([sys.executable, "-m", "acpype", "-i", f"{name}.mol2", "-b", name,
             "-c", "user", "-n", str(nc), "-a", "gaff2", "-o", "gmx"], work)

        for f in (f"{name}.itp", f"{name}_GMX.itp", f"{name}_atomtypes.itp", f"{name}.top"):
            src = work / f
            if src.exists():
                shutil.copy(src, LIG / f)
        print("   ->", sorted(x.name for x in LIG.glob(f"{name}*")))


if __name__ == "__main__":
    main()
