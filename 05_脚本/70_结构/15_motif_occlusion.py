"""15_motif_occlusion.py — PFNA 是否遮蔽 SCRN2 的 KFERQ 样基序？

问题：论文（Adv Sci 2025, PMID 39836524）报告 SCRN2 含两个 KFERQ 样基序
（M1 = 251-ELLRQ，M2 = 351-QVDRR），突变其一即削弱 HSPA8/LAMP2A 结合。
若要论证"EDC 通过 CMA 间接影响 SCRN2"，必须先回答一个几何问题：
PFNA 的结合位点是否覆盖（遮蔽）了这两个基序？

方法：Shrake-Rupley 溶剂可及表面积（probe 1.4 A，92 点），
分别计算 脱辅基蛋白（apo）与 蛋白+PFNA 复合物（holo）的逐原子 SASA，
作差得到配体结合引起的埋藏面积 ΔSASA（>0 表示被配体挡住）。

同时给出配体到每个基序残基的最近重原子距离，
以及配体实际埋藏的界面残基清单（用于判断"谁被挡住了"）。

用法：
    python scripts/15_motif_occlusion.py                 # 默认 DATA10 位点A
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
    """Golden-spiral 均匀球面点（Shrake-Rupley 常用实现）。"""
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
                continue                      # 重原子 SASA 即可
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
    """逐原子 SASA（A^2）。blockers 参与遮挡但不输出面积。"""
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

    print(f"蛋白重原子 {len(P['xyz'])}；配体原子 {len(L['xyz'])} "
          f"({', '.join(sorted(set(L['elem'])))})")
    print(f"配体文件：{a.lig}")

    s_apo = sasa(P["xyz"], rad_p)
    s_holo = sasa(P["xyz"], rad_p, L["xyz"], rad_l)
    dsasa = s_apo - s_holo

    # 逐残基汇总
    res_ids = sorted(set(P["resi"].tolist()))
    per_res = {}
    for r in res_ids:
        m = P["resi"] == r
        per_res[r] = dict(resn=str(P["resn"][m][0]),
                          apo=float(s_apo[m].sum()),
                          holo=float(s_holo[m].sum()),
                          buried=float(dsasa[m].sum()))

    MOTIFS = {"M1": (251, 255, "ELLRQ"), "M2": (351, 355, "QVDRR")}

    print(f"\n总埋藏面积 {dsasa.sum():.1f} A^2  "
          f"（apo {s_apo.sum():.0f} → holo {s_holo.sum():.0f} A^2）")

    report = {"ligand": a.lig, "total_buried": float(dsasa.sum()),
              "apo_total": float(s_apo.sum()), "holo_total": float(s_holo.sum()),
              "motifs": {}, "interface": []}

    print("\n=== KFERQ 样基序被配体埋藏的程度 ===")
    for key, (lo, hi, seq) in MOTIFS.items():
        rows = [(r, per_res[r]) for r in range(lo, hi + 1) if r in per_res]
        apo = sum(x["apo"] for _, x in rows)
        bur = sum(x["buried"] for _, x in rows)
        frac = bur / apo if apo > 0 else 0.0
        # 最近重原子距离
        dmin = min(np.linalg.norm(P["xyz"][P["resi"] == r][:, None, :]
                                  - L["xyz"][None, :, :], axis=2).min()
                   for r, _ in rows)
        print(f"\n{key} = {lo}-{hi} {seq}   apo SASA {apo:6.1f} A^2  "
              f"被埋藏 {bur:5.1f} A^2 ({frac*100:4.1f}%)  最近重原子距 {dmin:5.2f} A")
        for r, x in rows:
            flag = " <-- 埋藏>5A^2" if x["buried"] > 5 else ""
            print(f"   {x['resn']}{r:<4d} apo {x['apo']:6.1f}  holo {x['holo']:6.1f}"
                  f"  Δ {x['buried']:5.1f} A^2{flag}")
        report["motifs"][key] = dict(range=[lo, hi], seq=seq, apo=apo,
                                     buried=bur, frac=frac, dmin=float(dmin),
                                     per_residue={str(r): x for r, x in rows})

    # 界面残基（被埋藏最多的 15 个）
    top = sorted(per_res.items(), key=lambda kv: -kv[1]["buried"])[:15]
    print("\n=== 被 PFNA 埋藏最多的残基（界面）===")
    for r, x in top:
        if x["buried"] < 2:
            continue
        inmotif = " [M1]" if 251 <= r <= 255 else (" [M2]" if 351 <= r <= 355 else "")
        print(f"   {x['resn']}{r:<4d} Δ {x['buried']:6.1f} A^2  "
              f"({x['buried']/max(x['apo'],1e-6)*100:4.0f}% 自身){inmotif}")
        report["interface"].append(dict(resi=r, resn=x["resn"],
                                        buried=x["buried"],
                                        frac=x["buried"] / max(x["apo"], 1e-6)))

    OUT.mkdir(exist_ok=True)
    (OUT / "motif_occlusion.json").write_text(json.dumps(report, indent=1))
    print(f"\n-> {OUT / 'motif_occlusion.json'}")


if __name__ == "__main__":
    main()
