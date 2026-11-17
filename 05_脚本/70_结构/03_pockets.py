"""03 · SCRN2 口袋普查（LIGSITE 式）

前提问题：SCRN2 没有已知的生理配体，也没有共晶结构，结合位点只能从结构上找。
所以先量清楚 SCRN2 到底有哪些凹腔：若连一个能容纳多环疏水分子的封闭空腔都没有，
后面的 MD 就只能是对表面沟槽的采样，而不是对口袋的采样。

判据（LIGSITE/PASS 的简化）：
  1. 0.8 Å 网格；格点落在任一重原子 1.8 Å 内记为"蛋白占据"
  2. 对 ±x/±y/±z 六个方向，看前方 12 Å 内是否有蛋白占据
  3. 六个方向中被堵住 >=4 个的"空点"= 口袋点（平坦表面只有 1–2 个方向被堵）
  4. 26-连通聚类，输出体积、三轴跨度、衬里残基、平均被堵方向数
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
AF = ROOT / "structures" / "SCRN2_AF.pdb"
GRID = 0.8
PROBE = 2.6          # 空点：距最近原子 > 2.6 Å
OCC = 1.8            # 占据：距最近原子 < 1.8 Å
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
    """沿 axis 的 ±方向，o[i] 前方 k 步内是否有 True"""
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
    print(f"蛋白: {len(xyz)} 重原子 / {len(set(res))} 残基")
    print(f"pLDDT: 均值 {bfac.mean():.1f}；<50 残基 {int((per < 50).sum())} 个，"
          f"<70 残基 {int((per < 70).sum())} 个")

    lo, hi = xyz.min(0) - 3, xyz.max(0) + 3
    axes = [np.arange(lo[i], hi[i] + GRID, GRID) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    flat = grid.reshape(-1, 3)
    dist, _ = cKDTree(xyz).query(flat, workers=-1)
    dist = dist.reshape(grid.shape[:3])
    void, occ = dist > PROBE, dist < OCC
    print(f"网格 {grid.shape[:3]}；空点 {int(void.sum())}；占据 {int(occ.sum())}")

    cnt = np.zeros(occ.shape, np.uint8)
    for ax in range(3):
        for sg in (1, -1):
            cnt += wall_along(occ, ax, WALL, sg)
    pocket = void & (cnt >= 4)
    print(f"口袋点(>=4 方向被堵): {int(pocket.sum())}；深口袋(>=5): {int((void & (cnt >= 5)).sum())}")

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
    print(f"\n检出 {len(t)} 个口袋（前 10）：")
    print(t.head(10)[["id", "vol", "walls", "L1", "L2", "L3", "n_res", "res"]].to_string(index=False))
    fit = t[(t.L1 >= 11) & (t.L2 >= 5)]
    print(f"\n主轴 >=11 Å 且次轴 >=5 Å（可容纳多环疏水分子）的口袋: {len(fit)} 个")
    if len(fit):
        print(fit[["id", "vol", "walls", "L1", "L2", "L3", "res"]].to_string(index=False))


if __name__ == "__main__":
    main()
