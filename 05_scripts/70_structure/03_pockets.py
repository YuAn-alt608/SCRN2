"""03 - SCRN2 pocket survey (LIGSITE style)

Premise: SCRN2 has no known physiological ligand and no co-crystal structure, so the binding site
can only be located structurally. The first question is therefore what cavities SCRN2 actually has:
if there is not even one enclosed cavity able to hold a polycyclic hydrophobic molecule, the MD that
follows samples a surface groove rather than a pocket.

Criteria (a simplification of LIGSITE/PASS):
  1. 0.8 A grid; a grid point within 1.8 A of any heavy atom counts as "protein-occupied"
  2. along the six +/-x/y/z directions, test whether protein is occupied within 12 A ahead
  3. a "void point" with >=4 of the six directions blocked is a pocket point (a flat surface blocks 1-2)
  4. 26-connected clustering; report volume, the three axis spans, lining residues, mean blocked directions
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
AF = ROOT / "structures" / "SCRN2_AF.pdb"
GRID = 0.8
PROBE = 2.6          # void: farther than 2.6 A from the nearest atom
OCC = 1.8            # occupied: closer than 1.8 A to the nearest atom
WALL = 15            # 12 Å / 0.8 Å


def load_pdb(path):
    xyz, res, b = [], [], []
    for ln in Path(path).read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            el = (ln[76:78].strip() or ln[12:16].strip()[0])
            if el == "H":
                continue
            xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            res.append((ln[17:20].strip(), int(ln[22:26])))
            b.append(float(ln[60:66]))
    return np.array(xyz), res, np.array(b)


def wall_along(o, axis, k, sign):
    """Along the +/- directions of an axis, is there a True within k steps ahead of o[i]"""
    if sign < 0:
        o = np.flip(o, axis=axis)
    n = o.shape[axis]
    cs = np.cumsum(o, axis=axis, dtype=np.int32)
    shape = list(cs.shape)
    shape[axis] = 1
    cs2 = np.concatenate([np.zeros(shape, np.int32), cs], axis=axis)
    lo = np.take(cs2, np.arange(n) + 1, axis=axis)
    hi = np.take(cs2, np.clip(np.arange(n) + 1 + k, 0, n), axis=axis)
    out = (hi - lo) > 0
    return np.flip(out, axis=axis) if sign < 0 else out


def main():
    xyz, res, bfac = load_pdb(AF)
    df = pd.DataFrame(dict(res=[r[1] for r in res], plddt=bfac))
    per = df.groupby("res").plddt.mean()
    print(f"protein: {len(xyz)} heavy atoms / {len(set(res))} residues")
    print(f"pLDDT: mean {bfac.mean():.1f}; {int((per < 50).sum())} residues below 50,"
          f" {int((per < 70).sum())} below 70")

    lo, hi = xyz.min(0) - 3, xyz.max(0) + 3
    axes = [np.arange(lo[i], hi[i] + GRID, GRID) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    flat = grid.reshape(-1, 3)
    dist, _ = cKDTree(xyz).query(flat, workers=-1)
    dist = dist.reshape(grid.shape[:3])
    void, occ = dist > PROBE, dist < OCC
    print(f"grid {grid.shape[:3]}; void points {int(void.sum())}; occupied {int(occ.sum())}")

    cnt = np.zeros(occ.shape, np.uint8)
    for ax in range(3):
        for sg in (1, -1):
            cnt += wall_along(occ, ax, WALL, sg)
    pocket = void & (cnt >= 4)
    print(f"pocket points (>=4 directions blocked): {int(pocket.sum())}; deep pockets (>=5): {int((void & (cnt >= 5)).sum())}")

    tree = cKDTree(xyz)
    lab, n = ndimage.label(pocket, structure=np.ones((3, 3, 3)))
    rows = []
    for cid in range(1, n + 1):
        idx = np.argwhere(lab == cid)
        if len(idx) < 60:
            continue
        pts = grid[idx[:, 0], idx[:, 1], idx[:, 2]]
        c = pts.mean(0)
        w, v = np.linalg.eigh(np.cov((pts - c).T))
        v = v[:, np.argsort(w)[::-1]]
        proj = (pts - c) @ v
        ext = proj.max(0) - proj.min(0)
        near = tree.query_ball_point(pts, 4.5)
        uid = sorted({res[j][1] for lst in near for j in lst})
        rows.append(dict(id=cid, vol=round(len(pts) * GRID ** 3, 1),
                         walls=round(float(cnt[idx[:, 0], idx[:, 1], idx[:, 2]].mean()), 2),
                         x=round(c[0], 1), y=round(c[1], 1), z=round(c[2], 1),
                         L1=round(ext[0], 1), L2=round(ext[1], 1), L3=round(ext[2], 1),
                         n_res=len(uid), res=",".join(map(str, uid[:16]))))

    t = pd.DataFrame(rows).sort_values("vol", ascending=False).reset_index(drop=True)
    t.to_csv(ROOT / "results" / "pockets_ligsite.csv", index=False)
    pd.set_option("display.width", 220, "display.max_colwidth", 70)
    print(f"\n{len(t)} pockets found (top 10):")
    print(t.head(10)[["id", "vol", "walls", "L1", "L2", "L3", "n_res", "res"]].to_string(index=False))
    fit = t[(t.L1 >= 11) & (t.L2 >= 5)]
    print(f"\npockets with a long axis >=11 A and a short axis >=5 A (able to hold a polycyclic hydrophobic molecule): {len(fit)}")
    if len(fit):
        print(fit[["id", "vol", "walls", "L1", "L2", "L3", "res"]].to_string(index=False))


if __name__ == "__main__":
    main()
