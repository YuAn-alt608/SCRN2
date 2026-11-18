"""mdlib - parse GROMACS topologies and multi-frame PDB trajectories for 07/08

Why the parsing is hand-rolled:
  * neither MDAnalysis nor mdtraj is available locally (no wheel for Python 3.13)
  * the gmx tools select groups through interactive stdin, which scripts get wrong
  * MM-GBSA needs per-frame charges, LJ parameters and connectivity, all of which
    exist only in the topology
"""
from pathlib import Path

import numpy as np

# Atom-name sets used to isolate backbone atoms; the element table is only needed
# to label atom types that the force field stores as an atomic number.
BACKBONE = {"N", "CA", "C", "O"}

ELEM_BY_Z = {1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 11: "Na",
             15: "P", 16: "S", 17: "Cl", 35: "Br", 53: "I", 19: "K"}


def read_top(path, n_moltypes=2):
    """Read [atomtypes] and [atoms]/[bonds] from the grompp-expanded processed.top.

    Only the first n_moltypes moleculetypes are kept (0 = protein, 1 = ligand);
    solvent and ions are skipped. The returned dict follows exactly the atom order
    of the gro/pdb files, so indices can be used interchangeably.
    """
    lines = Path(path).read_text(errors="ignore").splitlines()
    at = {}
    elem, charge, mass, aname, resname = [], [], [], [], []
    bonds = []
    sec = None
    mt = -1
    offset = 0
    local_n = 0
    for ln in lines:
        s = ln.split(";")[0].strip()
        if not s:
            continue
        if s.startswith("["):
            sec = s.strip("[] ").lower().split()[0]
            if sec == "moleculetype":
                mt += 1
                offset += local_n
                local_n = 0
            continue
        t = s.split()
        if sec == "atomtypes" and len(t) >= 7:
            # name  bond_type  mass  charge  ptype  V(sigma)  W(eps)
            # note: in the amber force fields column 2 is an atomic number, whereas
            # hand-written GAFF blocks put an element symbol there
            bt = t[1]
            e = ELEM_BY_Z.get(int(bt), "X") if bt.isdigit() else bt
            at[t[0]] = (e, float(t[2]), float(t[5]), float(t[6]))
        elif sec == "atoms" and mt < n_moltypes and len(t) >= 8 and t[0].isdigit():
            e, m, sig, eps = at[t[1]]
            elem.append(e); mass.append(m); charge.append(float(t[6]))
            aname.append(t[4]); resname.append(t[3])
            local_n += 1
        elif sec == "bonds" and mt < n_moltypes and len(t) >= 2:
            bonds.append((offset + int(t[0]) - 1, offset + int(t[1]) - 1))

    types = []
    # Second pass: collect each atom's atomtype name, then map it to sigma/eps
    mt = -1
    sec = None
    for ln in lines:
        s = ln.split(";")[0].strip()
        if not s:
            continue
        if s.startswith("["):
            sec = s.strip("[] ").lower().split()[0]
            if sec == "moleculetype":
                mt += 1
            continue
        t = s.split()
        if sec == "atoms" and mt < n_moltypes and len(t) >= 8 and t[0].isdigit():
            types.append(t[1])

    sigma = np.array([at[t][2] for t in types])
    eps = np.array([at[t][3] for t in types])
    return dict(
        elem=np.array(elem), charge=np.array(charge), mass=np.array(mass),
        sigma=sigma, eps=eps, aname=np.array(aname), resname=np.array(resname),
        types=np.array(types),
        bonds=np.array(bonds, dtype=int).reshape(-1, 2),
    )


def read_pdb_frames(path):
    """Stream a multi-frame PDB written by trjconv; returns (coords[nf,na,3] nm, boxes[nf,3,3] nm)."""
    frames, boxes = [], []
    xyz = []
    box = np.eye(3) * 8.0
    with open(path, errors="ignore") as fh:
        for ln in fh:
            rec = ln[:6]
            if rec == "CRYST1":
                a, b, c = float(ln[6:15]), float(ln[15:24]), float(ln[24:33])
                al, be, ga = float(ln[33:40]), float(ln[40:47]), float(ln[47:54])
                box = cell_to_vectors(a, b, c, al, be, ga) / 10.0    # Å -> nm
            elif rec in ("ATOM  ", "HETATM"):
                xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
            elif rec.startswith("ENDMDL"):
                frames.append(np.array(xyz) / 10.0)
                boxes.append(box.copy())
                xyz = []
    if xyz:
        frames.append(np.array(xyz) / 10.0)
        boxes.append(box.copy())
    return np.array(frames), np.array(boxes)


def load_frames(work, pdb_name="traj.pdb"):
    """Load trajectory frames: frames.npz if present (float32, about 1/6 the size of
    the PDB), otherwise fall back to the multi-frame PDB.

    frames.npz is written by 00_export_frames.py to cut bandwidth when the analysis
    is uploaded to a cloud (Colab) runtime. Both paths return nm coordinates of
    shape (nf, na, 3).
    """
    work = Path(work)
    npz = work / "frames.npz"
    if npz.exists():
        z = np.load(npz)
        return z["frames"].astype(float), z["boxes"].astype(float)
    return read_pdb_frames(work / pdb_name)


def cell_to_vectors(a, b, c, al, be, ga):
    """Rebuild the GROMACS box vectors (as rows) from the cell parameters."""
    al, be, ga = np.radians([al, be, ga])
    v1 = np.array([a, 0.0, 0.0])
    v2 = np.array([b * np.cos(ga), b * np.sin(ga), 0.0])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz = np.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return np.array([v1, v2, [cx, cy, cz]])


def fit_transform(mobile, ref):
    """Kabsch superposition (row-vector convention): return R, t minimising
    |mobile @ R + t - ref|.

    With H = mobile_c^T @ ref_c, the optimal rotation is R = U D V^T, not V D U^T:
    the latter is its transpose, i.e. a rotation the other way round, and that mistake
    inflates the RMSD several-fold while still producing numbers that look plausible.
    """
    mc, rc = mobile.mean(0), ref.mean(0)
    H = (mobile - mc).T @ (ref - rc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1.0, 1.0, d])
    R = U @ D @ Vt
    return R, rc - mc @ R


def unwrap_ligand(lig, prot_c, box):
    """Translate the ligand as a rigid body to its nearest periodic image of the
    protein (after -pbc mol it can still sit in a neighbouring cell)."""
    shift_best = np.zeros(3)
    best = np.linalg.norm(lig.mean(0) - prot_c)
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                s = i * box[0] + j * box[1] + k * box[2]
                d = np.linalg.norm(lig.mean(0) + s - prot_c)
                if d < best:
                    best, shift_best = d, s
    return lig + shift_best
