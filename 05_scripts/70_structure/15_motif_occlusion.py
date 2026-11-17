"""15_motif_occlusion.py - does PFNA occlude the KFERQ-like motifs of SCRN2?

Question: the report in Adv Sci 2025 (PMID 39836524) states that SCRN2 carries two KFERQ-like
motifs (M1 = 251-ELLRQ, M2 = 351-QVDRR) and that mutating either one weakens HSPA8/LAMP2A
binding. Any claim that EDCs act on SCRN2 through CMA therefore runs through one geometric
question: does the PFNA binding site cover, or occlude, these two motifs?

Method: Shrake-Rupley solvent-accessible surface area (probe 1.4 A, 92 points), computed
per atom for the apo protein and for the protein+PFNA complex, then differenced to give the
area buried by ligand binding (dSASA > 0 means hidden by the ligand).

The closest heavy-atom distance from the ligand to each motif residue is reported alongside,
together with the list of interface residues the ligand actually buries.

Usage:
    python 05_scripts/70_structure/15_motif_occlusion.py     # defaults built into the script
    python scripts/15_motif_occlusion.py --lig <mol2>
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "F": 1.47,
       "P": 1.85, "CL": 1.75, "BR": 1.85, "I": 1.98, "H": 1.20}
PROBE = 1.4
N_SPHERE = 92


def sphere_points(n=N_SPHERE):
    """Uniform points on a sphere by golden-spiral sampling (the usual Shrake-Rupley implementation)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


SPH = sphere_points()


def read_pdb(path):
    elem, resn, resi, aname, xyz = [], [], [], [], []
    for l in open(path, errors="ignore"):
        if l.startswith(("ATOM", "HETATM")):
            e = l[76:78].strip().upper()
            if not e:
                e = l[12:16].strip()[0].upper()
            if e == "H":
                continue                      # heavy-atom SASA is enough
            elem.append(e)
            aname.append(l[12:16].strip())
            resn.append(l[17:20].strip())
            resi.append(int(l[22:26]))
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return dict(elem=np.array(elem), resn=np.array(resn), resi=np.array(resi),
                aname=np.array(aname), xyz=np.array(xyz, dtype=float))


def read_mol2(path):
    elem, aname, xyz = [], [], []
    sec = None
    for l in open(path, errors="ignore"):
        if l.startswith("@<TRIPOS>"):
            sec = l.strip()
            continue
        if sec != "@<TRIPOS>ATOM":
            continue
        t = l.split()
        if len(t) < 6:
            continue
        typ = t[5].split(".")[0].upper()
        elem.append(typ)
        aname.append(t[1])
        xyz.append([float(t[2]), float(t[3]), float(t[4])])
    return dict(elem=np.array(elem), aname=np.array(aname),
                xyz=np.array(xyz, dtype=float))


def sasa(xyz, rad, blockers_xyz=None, blockers_rad=None):
    """Per-atom SASA (A^2). Atoms in blockers occlude but contribute no area."""
    if blockers_xyz is not None and len(blockers_xyz):
        allx = np.vstack([xyz, blockers_xyz])
        allr = np.concatenate([rad, blockers_rad])
        n_own = len(xyz)
    else:
        allx, allr = xyz, rad
        n_own = len(xyz)
    tree = cKDTree(allx)
    rmax = allr.max()
    out = np.zeros(len(xyz))
    for i in range(len(xyz)):
        ri = allr[i] + PROBE
        idx = tree.query_ball_point(allx[i], ri + rmax + 0.1)
        nb = [j for j in idx if j != i]
        if not nb:
            out[i] = 4.0 * np.pi * ri ** 2
            continue
        nb = np.array(nb)
        p = xyz[i] + ri * SPH                      # (92,3)
        d = np.linalg.norm(p[:, None, :] - allx[nb][None, :, :], axis=2)
        blocked = (d < (allr[nb] + PROBE)[None, :]).any(1)
        out[i] = (4.0 * np.pi * ri ** 2) * (~blocked).mean()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lig", default=str(
        Path("D:/NewBeginning/MR1/DATA10/docking/gromacs/PFNA_docked.mol2")))
    ap.add_argument("--pdb", default=str(ROOT / "structures" / "SCRN2_AF.pdb"))
    a = ap.parse_args()

    P = read_pdb(a.pdb)
    L = read_mol2(a.lig)
    rad_p = np.array([VDW.get(e, 1.70) for e in P["elem"]])
    rad_l = np.array([VDW.get(e, 1.70) for e in L["elem"]])

    print(f"protein heavy atoms {len(P['xyz'])}; ligand atoms {len(L['xyz'])} "
          f"({', '.join(sorted(set(L['elem'])))})")
    print(f"ligand file: {a.lig}")

    s_apo = sasa(P["xyz"], rad_p)
    s_holo = sasa(P["xyz"], rad_p, L["xyz"], rad_l)
    dsasa = s_apo - s_holo

    # per-residue totals
    res_ids = sorted(set(P["resi"].tolist()))
    per_res = {}
    for r in res_ids:
        m = P["resi"] == r
        per_res[r] = dict(resn=str(P["resn"][m][0]),
                          apo=float(s_apo[m].sum()),
                          holo=float(s_holo[m].sum()),
                          buried=float(dsasa[m].sum()))

    MOTIFS = {"M1": (251, 255, "ELLRQ"), "M2": (351, 355, "QVDRR")}

    print(f"\ntotal buried area {dsasa.sum():.1f} A^2  "
          f"(apo {s_apo.sum():.0f} -> holo {s_holo.sum():.0f} A^2)")

    report = {"ligand": a.lig, "total_buried": float(dsasa.sum()),
              "apo_total": float(s_apo.sum()), "holo_total": float(s_holo.sum()),
              "motifs": {}, "interface": []}

    print("\n=== how far the ligand buries each KFERQ-like motif ===")
    for key, (lo, hi, seq) in MOTIFS.items():
        rows = [(r, per_res[r]) for r in range(lo, hi + 1) if r in per_res]
        apo = sum(x["apo"] for _, x in rows)
        bur = sum(x["buried"] for _, x in rows)
        frac = bur / apo if apo > 0 else 0.0
        # closest heavy-atom distance
        dmin = min(np.linalg.norm(P["xyz"][P["resi"] == r][:, None, :]
                                  - L["xyz"][None, :, :], axis=2).min()
                   for r, _ in rows)
        print(f"\n{key} = {lo}-{hi} {seq}   apo SASA {apo:6.1f} A^2  "
              f"buried {bur:5.1f} A^2 ({frac*100:4.1f}%)  closest heavy atom {dmin:5.2f} A")
        for r, x in rows:
            flag = " <-- buried >5A^2" if x["buried"] > 5 else ""
            print(f"   {x['resn']}{r:<4d} apo {x['apo']:6.1f}  holo {x['holo']:6.1f}"
                  f"  Δ {x['buried']:5.1f} A^2{flag}")
        report["motifs"][key] = dict(range=[lo, hi], seq=seq, apo=apo,
                                     buried=bur, frac=frac, dmin=float(dmin),
                                     per_residue={str(r): x for r, x in rows})

    # interface residues (the 15 most buried)
    top = sorted(per_res.items(), key=lambda kv: -kv[1]["buried"])[:15]
    print("\n=== residues most buried by PFNA (the interface) ===")
    for r, x in top:
        if x["buried"] < 2:
            continue
        inmotif = " [M1]" if 251 <= r <= 255 else (" [M2]" if 351 <= r <= 355 else "")
        print(f"   {x['resn']}{r:<4d} Δ {x['buried']:6.1f} A^2  "
              f"({x['buried']/max(x['apo'],1e-6)*100:4.0f}% of its own){inmotif}")
        report["interface"].append(dict(resi=r, resn=x["resn"],
                                        buried=x["buried"],
                                        frac=x["buried"] / max(x["apo"], 1e-6)))

    OUT.mkdir(exist_ok=True)
    (OUT / "motif_occlusion.json").write_text(json.dumps(report, indent=1))
    print(f"\n-> {OUT / 'motif_occlusion.json'}")


if __name__ == "__main__":
    main()
