"""05 - build the MD system: put the docked pose back on the protein, make a periodic box,
solvate and ionise, then energy-minimise

Key points:
  * the [atomtypes] block of the ligand itp is split out into LIG_atomtypes.itp and placed
    after the force-field include and before any [moleculetype], since GROMACS requires atom
    types to be defined first
  * atom type names carry a G prefix to avoid clashing with c3/ca and friends in amber99sb-ildn
  * a 400 kJ/mol/nm^2 positional restraint on the ligand during equilibration, released in production
  * custom index groups Solute (protein + ligand) and Solvent (water + ions) for thermostat coupling
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
LIG = ROOT / "ligands"
RES = ROOT / "results"
PREP = ROOT / "md_prep"
GMX = "gmx"


def sh(cmd, cwd, inp=None, check=True):
    r = subprocess.run(cmd, cwd=str(cwd), input=inp, capture_output=True, text=True, shell=False)
    if check and r.returncode != 0:
        print(r.stdout[-2500:]); print(r.stderr[-1500:])
        raise SystemExit(f"FAILED: {' '.join(map(str, cmd))}")
    return r


def read_pdb(path):
    xyz, name, resname, el = [], [], [], []
    for ln in Path(path).read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            name.append(ln[12:16].strip())
            resname.append(ln[17:20].strip())
            el.append((ln[76:78].strip() or ln[12:16].strip()[0]))
    return np.array(xyz), name, resname, el


def read_gro(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[1])
    return lines[0], n, lines[2:2 + n]


def fmt_gro(resnum, resname, aname, anum, xyz):
    return (f"{resnum % 100000:5d}{resname:<5s}{aname:>5s}{anum % 100000:5d}"
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}")


def apply_pose(name):
    """Apply the docking result (centre + quaternion) to the ligand coordinates."""
    import pandas as pd
    d = pd.read_csv(RES / f"dock_{name}.csv")
    row = d.sort_values("energy").iloc[0]          # lowest-energy pose
    c = np.array([row.cx, row.cy, row.cz], dtype=float)
    q = np.array([float(v) for v in str(row.quat).split(",")])
    xyz, anames, resname, el = read_pdb(LIG / f"{name}.pdb")
    heavy = np.array([e != "H" for e in el])
    com0 = xyz[heavy].mean(0)              # the docking centroid used heavy atoms only
    R = Rotation.from_quat(q).as_matrix()
    new = (xyz - com0) @ R.T + c
    return new, anames, el, resname[0]


def relieve_clash(xyz, heavy, tree, target=2.2, step=0.15, max_iter=60):
    """Push the ligand out of the protein as a rigid body.

    The direction is the mean unit vector from each clashing ligand atom to its nearest
    protein atom (the local outward normal), which is steadier than a centroid direction:
    a ligand wedged into a groove does not slide along the tangent. Translation only, so
    the pose is preserved.
    """
    xyz = xyz.copy()
    for _ in range(max_iter):
        d, idx = tree.query(xyz[heavy])
        if d.min() >= target:
            break
        near = tree.data[idx]
        vec = xyz[heavy] - near                       # each heavy atom relative to its nearest protein atom
        clash = d < target
        vec = vec[clash]
        n = np.linalg.norm(vec, axis=1, keepdims=True)
        vec = vec / np.where(n < 1e-6, 1.0, n)
        xyz += step * vec.mean(0)
    return xyz


def build_one(name):
    work = ROOT / "md_runs" / name
    work.mkdir(parents=True, exist_ok=True)
    print(f"\n===== {name} =====")

    lig_xyz, anames, el, ligres = apply_pose(name)
    nl = len(lig_xyz)
    heavy = np.array([e != "H" for e in el])

    # post-docking check: the ligand must not sit inside the protein (heavy atoms only; pdb2gmx added H)
    prot_xyz, _, _, _ = read_pdb(ROOT / "structures" / "SCRN2_AF.pdb")
    from scipy.spatial import cKDTree
    tree = cKDTree(prot_xyz)
    dmin = tree.query(lig_xyz[heavy])[0].min()
    print(f"  ligand {nl} atoms, closest protein heavy atom {dmin:.2f} A")
    if dmin < 2.2:                       # the score tolerates 2.2 A, but poses often interpenetrate within 1 A
        lig_xyz = relieve_clash(lig_xyz, heavy, tree, 2.2)
        dmin2 = tree.query(lig_xyz[heavy])[0].min()
        print(f"  clash relief: {dmin:.2f} -> {dmin2:.2f} A")

    # protein gro (already written by pdb2gmx, hydrogens included)
    title, np_, prot_lines = read_gro(PREP / "prot.gro")

    # complex gro: protein first, ligand second, matching the topology order
    out = [f"{name} complex", f"{np_ + nl:5d}"] + prot_lines
    lig_nm = lig_xyz / 10.0          # gro is in nm; the ligand coordinates come from a PDB in A
    for k in range(nl):
        out.append(fmt_gro(1, ligres[:5], anames[k][:5], np_ + k + 1, lig_nm[k]))
    # box vectors: whatever pdb2gmx wrote is meaningless (an AF model has no unit cell), so a box
    # is given from the coordinate range + 2 nm and rebuilt by editconf -d
    allx = np.vstack([np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])]
                                for l in prot_lines]), lig_nm])
    box = allx.max(0) - allx.min(0) + 2.0        # nm
    out.append("   " + "".join(f"{float(v):10.5f}" for v in box))
    (work / "complex.gro").write_text("\n".join(out) + "\n", encoding="ascii")

    # split the ligand topology: atomtypes get their own file
    itp = (LIG / f"{name}.itp").read_text(encoding="utf-8")
    m = re.search(r"\[ atomtypes \]\n(.*?)\n\n", itp, re.S)
    if not m:
        raise SystemExit("ligand itp lacks [ atomtypes ]")
    (work / "LIG_atomtypes.itp").write_text(
        "[ atomtypes ]\n" + m.group(1) + "\n", encoding="ascii")
    body = itp[:m.start()] + itp[m.end():]
    body = re.sub(r"^;.*\n", "", body, count=1)

    # ligand positional restraints (active during equilibration only)
    pr = ["[ position_restraints ]", "; atom  type  fx  fy  fz"]
    for i in range(1, nl + 1):
        pr.append(f"{i:6d} 1 400.0 400.0 400.0")
    (work / "LIG_posre.itp").write_text("\n".join(pr) + "\n", encoding="ascii")
    body += '\n#ifdef POSRES_LIG\n#include "LIG_posre.itp"\n#endif\n'
    (work / f"{name}.itp").write_text(body, encoding="utf-8")

    # top-level topology
    top = (PREP / "topol_prot.top").read_text()
    top = top.replace('#include "amber99sb-ildn.ff/forcefield.itp"',
                      '#include "amber99sb-ildn.ff/forcefield.itp"\n'
                      '#include "LIG_atomtypes.itp"', 1)
    top = top.replace('#include "amber99sb-ildn.ff/tip3p.itp"',
                      f'#include "{name}.itp"\n'
                      '#include "amber99sb-ildn.ff/tip3p.itp"', 1)
    top = re.sub(r"(Protein_chain_A\s+1\n)", r"\1" + f"{name}                1\n", top, count=1)
    (work / "topol.top").write_text(top)
    shutil.copy(PREP / "posre_prot.itp", work / "posre_prot.itp")

    # box + solvent + ions
    sh([GMX, "editconf", "-f", "complex.gro", "-o", "box.gro",
        "-bt", "triclinic", "-d", "1.0", "-c"], work)
    sh([GMX, "solvate", "-cp", "box.gro", "-cs", "spc216.gro",
        "-o", "solv.gro", "-p", "topol.top"], work)
    sh([GMX, "grompp", "-f", str(ROOT / "mdp" / "ions.mdp"), "-c", "solv.gro",
        "-p", "topol.top", "-o", "ions.tpr", "-maxwarn", "2"], work)
    r = sh([GMX, "genion", "-s", "ions.tpr", "-o", "solv_ions.gro", "-p", "topol.top",
            "-pname", "NA", "-nname", "CL", "-conc", "0.15", "-neutral"], work,
           inp="SOL\n")
    print("  " + [ln for ln in r.stdout.splitlines() if "charge" in ln.lower()][-1:].__str__())
    print("  " + [ln for ln in r.stdout.splitlines() if "replaced" in ln.lower()][-1:].__str__())

    # index: Solute = protein + ligand. GROMACS has a short line buffer, so atom numbers wrap (<=15 per line)
    _, ntot, _ = read_gro(work / "solv_ions.gro")
    nsolute = np_ + nl

    def grp(tag, lo, hi):
        out = [f"[ {tag} ]"]
        ids = list(range(lo, hi + 1))
        for i in range(0, len(ids), 15):
            out.append(" ".join(str(v) for v in ids[i:i + 15]))
        return out

    idx = grp("Solute", 1, nsolute) + grp("Solvent", nsolute + 1, ntot) + grp("System", 1, ntot)
    (work / "index.ndx").write_text("\n".join(idx) + "\n", encoding="ascii")
    print(f"  {ntot} atoms in total ({nsolute} solute)")

    sh([GMX, "grompp", "-f", str(ROOT / "mdp" / "em.mdp"), "-c", "solv_ions.gro",
        "-r", "solv_ions.gro", "-p", "topol.top", "-n", "index.ndx",
        "-o", "em.tpr", "-pp", "processed.top", "-maxwarn", "2"], work)
    sh([GMX, "mdrun", "-deffnm", "em", "-ntmpi", "1"], work)
    tail = (work / "em.log").read_text(errors="ignore").splitlines()
    print("  " + [ln for ln in tail if "Maximum force" in ln or "Potential Energy" in ln][-2:].__str__())


if __name__ == "__main__":
    for nm in (sys.argv[1:] or ["PFNA"]):
        build_one(nm)
