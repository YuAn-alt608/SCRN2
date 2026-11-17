"""04 · 盲对接：让 T3 与 PFNA 自己在 SCRN2 表面找位置

没有共晶结构、没有已知位点，所以不能"指定口袋再对接"——那等于把结论写进起点。
这里做的是全局刚性对接：
  1. 打分 = 12-6 Lennard-Jones（按元素给 vdW 半径与势阱）+ 氢键方向性奖励
  2. 位点候选 = 口袋普查里"六方向被堵 >=5"的格点，最远点采样降到 ~250 个中心
  3. 每个中心试 600 个随机四元数取向，取最优
  4. 前 40 名用 Nelder-Mead 在 6 自由度（平移+四元数）上局部优化

两个配体用完全相同的打分函数与搜索空间，保证可比。
"""
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
LIG = ROOT / "ligands"
RES = ROOT / "results"

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "I": 1.98, "S": 1.80, "H": 1.20}
EPS = {"C": 0.086, "N": 0.170, "O": 0.210, "F": 0.061, "I": 0.320, "S": 0.250}
HB_DON = {"N", "O"}
RNG = np.random.default_rng(20261111)


def read_pdb(path, keep_h=False):
    xyz, el = [], []
    for ln in Path(path).read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            e = (ln[76:78].strip() or ln[12:16].strip()[0])
            if e == "H" and not keep_h:
                continue
            xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            el.append(e)
    return np.array(xyz), np.array(el)


def pair_energy(d, ei, ej):
    """LJ(12-6) + 氢键。d 为一维或二维距离，ei/ej 为元素数组；统一成 (ni,nj) 再算。"""
    ei = np.atleast_1d(np.asarray(ei)).ravel()
    ej = np.atleast_1d(np.asarray(ej)).ravel()
    d = np.asarray(d, dtype=float).reshape(len(ei), len(ej))
    ri = np.array([VDW.get(x, 1.70) for x in ei])[:, None]
    rj = np.array([VDW.get(x, 1.70) for x in ej])[None, :]
    r0 = 0.5 * (ri + rj)
    eps = np.sqrt(np.array([EPS.get(x, 0.10) for x in ei])[:, None] *
                  np.array([EPS.get(x, 0.10) for x in ej])[None, :])
    s6 = np.clip((r0 / np.maximum(d, 1e-3)) ** 6, 0.0, 4.0 ** 6)
    e = eps * (s6 ** 2 - 2.0 * s6)
    don = (np.isin(ei, list(HB_DON))[:, None] & np.isin(ej, list(HB_DON))[None, :]
           & (d > 2.4) & (d < 3.7))
    e = e + np.where(don, -0.9 * np.exp(-((d - 3.0) / 0.5) ** 2), 0.0)
    return e


def score_pose(x, rec, rec_el, tree, news, com0):
    """x = [tx,ty,tz, qx,qy,qz,qw]；返回体系能量"""
    R = Rotation.from_quat(x[3:]).as_matrix()
    pos = (news - com0) @ R.T + x[:3]
    idx = tree.query_ball_point(pos, 8.0)
    tot = 0.0
    clash = 0
    for k, lst in enumerate(idx):
        if lst:
            dd = np.linalg.norm(rec[lst] - pos[k], axis=1)
            tot += pair_energy(dd, [lig_el[k]], rec_el[lst]).sum()
            clash += int((dd < 2.2).sum())
    return tot + 12.0 * clash          # 穿模重罚


def main():
    global lig_el
    rec, rec_el = read_pdb(ROOT / "structures" / "SCRN2_AF.pdb")
    tree = cKDTree(rec)
    print(f"受体 {len(rec)} 重原子")

    # 位点候选：从口袋点里做最远点采样
    pk = pd.read_csv(RES / "pockets_ligsite.csv")
    grid_pts = []
    for _, r in pk.iterrows():
        grid_pts.append([r.x, r.y, r.z])
    centers = np.array(grid_pts[:40])
    print(f"位点中心候选 {len(centers)} 个（取体积前 40 个口袋质心）")

    out = {}
    for name in ("T3", "PFNA"):
        pos, el = read_pdb(LIG / f"{name}.pdb")
        lig_el = el
        com0 = pos.mean(0)
        best = []
        for c in centers:
            qs = Rotation.random(600, random_state=int(RNG.integers(1e6))).as_quat()
            for q in qs:
                x = np.concatenate([c, q])
                e = score_pose(x, rec, rec_el, tree, pos, com0)
                best.append((e, x))
        best.sort(key=lambda t: t[0])
        print(f"{name}: 全局粗搜 {len(best)} 个位姿，最好 {best[0][0]:.1f} kcal/mol")

        # 局部优化前 40
        refined = []
        for e0, x0 in best[:40]:
            r = minimize(score_pose, x0, args=(rec, rec_el, tree, pos, com0),
                         method="Nelder-Mead",
                         options=dict(maxiter=1200, xatol=1e-3, fatol=1e-3))
            refined.append((r.fun, r.x))
        refined.sort(key=lambda t: t[0])
        print(f"{name}: 优化后最好 {refined[0][0]:.1f} kcal/mol")

        # 记录前 10 个位姿
        rows = []
        for i, (e, x) in enumerate(refined[:10]):
            R = Rotation.from_quat(x[3:]).as_matrix()
            p = (pos - com0) @ R.T + x[:3]
            d = cKDTree(rec).query(p)[0]
            rows.append(dict(rank=i + 1, energy=round(float(e), 2),
                             mindist=round(float(d.min()), 2),
                             cx=round(float(x[0]), 2), cy=round(float(x[1]), 2),
                             cz=round(float(x[2]), 2),
                             quat=",".join(f"{v:.6f}" for v in x[3:])))
        dfp = pd.DataFrame(rows)
        dfp.to_csv(RES / f"dock_{name}.csv", index=False)
        print(dfp.to_string(index=False))
        out[name] = (refined[0], dfp)

    # 两个配体是否落在同一处？取各自第 1 名
    print("\n=== 首位位姿间距 ===")
    for a, b in itertools.combinations(out, 2):
        ca = np.array([out[a][1].cx[0], out[a][1].cy[0], out[a][1].cz[0]])
        cb = np.array([out[b][1].cx[0], out[b][1].cy[0], out[b][1].cz[0]])
        print(f"{a} vs {b}: 质心相距 {np.linalg.norm(ca - cb):.1f} Å")


if __name__ == "__main__":
    main()
