"""17_hspa8_insertion_test.py — 刚体插入 / 位阻冲突检验（纯几何，不做 MD）

问题：折叠态的 SCRN2，能不能把基序（QVDRR 351-355 / QDPVR 339-343）送进
HSPA8 的底物结合裂隙？PFNA 结合在 DATA14 位点上会不会挡住 SBD 的接近路径？

流程：
 1. 参照几何：DnaK SBD + 底物肽 NRLLLTG 的晶体结构（1DKX），提供"裂隙已含伸展肽"的真实坐标。
    （注：DATA13+ 的 complex_M2_WT.pdb 里 HSPA8 链是 510-616 = SBDα 螺旋盖，
      不是底物裂隙所在的 SBDβ 394-509，因此不能用作参照。）
 2. 把 AlphaFold 人 HSPA8 全长按 SBDβ（394-509）叠合到 1DKX 的 DnaK SBD 上。
 3. 把 SCRN2（DATA14 MD 末帧，含 PFNA）按基序主链叠合到 1DKX 肽段的 5 残基窗口（遍历 register）。
 4. 量冲突：HSPA8(SBDβ/SBDα/NBD) ↔ SCRN2、HSPA8(SBDβ) ↔ PFNA，以及对照（天然肽 vs DnaK）。

判据（重原子对距离）：
   severe  < 2.0 A   物理上不可能（范德华半径和 ~2.8-3.2 A）
   bad     < 2.5 A
   contact < 3.0 A
"""
import json
import os
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
# 参照结构随项目走（云端 / 他人复现时无需改路径）；仍可用环境变量覆盖
DKX = Path(os.environ.get("DKX_PDB", ROOT / "structures" / "1dkx.pdb"))
AFH = Path(os.environ.get("AFH_PDB", ROOT / "structures" / "AF-HSPA8-P11142.pdb"))

AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
       'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}


def read_pdb(path):
    el, ch, resn, resi, an, xyz = [], [], [], [], [], []
    for l in open(path, errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        e = l[76:78].strip().upper()
        if not e:
            e = l[12:16].strip().lstrip("0123456789")[0].upper()
        el.append(e); ch.append(l[21]); resn.append(l[17:20].strip())
        resi.append(int(l[22:26])); an.append(l[12:16].strip())
        xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return dict(elem=np.array(el), chain=np.array(ch), resn=np.array(resn),
                resi=np.array(resi), aname=np.array(an),
                xyz=np.array(xyz, dtype=float))


def read_gro(path):
    L = open(path, errors="ignore").read().splitlines()
    n = int(L[1])
    resid = np.array([int(l[:5]) for l in L[2:2 + n]])
    rn = np.array([l[5:10].strip() for l in L[2:2 + n]])
    an = np.array([l[10:15].strip() for l in L[2:2 + n]])
    xyz = np.array([[float(l[20 + 8 * i:28 + 8 * i]) for i in range(3)]
                    for l in L[2:2 + n]])
    el = []
    for a, r in zip(an, rn):
        a = a.lstrip("0123456789")
        if r == "PFNA":
            el.append("F" if a.startswith("F") else ("O" if a.startswith("O") else "C"))
        else:
            el.append(a[0] if a[0] in "CNOSPF" else "C")
    return dict(elem=np.array(el), resn=rn, resi=resid, aname=an, xyz=xyz)


def subset(d, m):
    return {k: v[m] for k, v in d.items()}


def kabsch(P, Q):
    """行向量：返回 R(3x3),t 使 P@R.T+t ≈ Q。"""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, Q.mean(0) - P.mean(0) @ R.T


def apply(x, R, t):
    return x @ R.T + t


def rmsd(P, Q, R, t):
    return float(np.sqrt(((apply(P, R, t) - Q) ** 2).sum(1).mean()))


def clash(A, B, cutoff=3.0):
    """原子对级别的冲突统计（与全矩阵广播等价，但用 KD 树，内存 O(邻居数)）。

    返回的 severe/bad/contact 是**原子对数**（不是残基数），
    与旧的全矩阵实现逐位一致，只是不必为 3300x1500 建临时数组。
    """
    if len(A) == 0 or len(B) == 0:
        return dict(min=99.0, severe=0, bad=0, contact=0)
    from scipy.spatial import cKDTree
    tree = cKDTree(B)
    nbrs = tree.query_ball_point(A, r=cutoff)
    ai, bi = [], []
    for i, idx in enumerate(nbrs):
        if idx:
            ai.extend([i] * len(idx))
            bi.extend(idx)
    if not bi:
        # 没有任何接触，最近距离单独查一次
        d, _ = tree.query(A, k=1)
        return dict(min=float(d.min()), severe=0, bad=0, contact=0)
    d = np.linalg.norm(A[np.asarray(ai)] - B[np.asarray(bi)], axis=1)
    return dict(min=float(d.min()), severe=int((d < 2.0).sum()),
                bad=int((d < 2.5).sum()), contact=int((d < cutoff).sum()))


def paired_bb(d1, d2, pairs):
    P, Q = [], []
    for a, b in pairs:
        for an in ("N", "CA", "C", "O"):
            i = np.where((d1["resi"] == a) & (d1["aname"] == an))[0]
            j = np.where((d2["resi"] == b) & (d2["aname"] == an))[0]
            if len(i) and len(j):
                P.append(d1["xyz"][i[0]]); Q.append(d2["xyz"][j[0]])
    return np.array(P), np.array(Q)


def ca_seq(d, lo, hi):
    r = [x for x in sorted(set(d["resi"])) if lo <= x <= hi]
    s = []
    for x in r:
        m = (d["resi"] == x) & (d["aname"] == "CA")
        s.append(AA3.get(d["resn"][m][0], "X") if m.any() else "X")
    return "".join(s), r


def main():
    af = read_pdb(AFH)
    dk = read_pdb(DKX)
    keep = dk["resn"] != "HOH"
    dk_a = subset(dk, keep & (dk["chain"] == "A"))
    dk_b = subset(dk, keep & (dk["chain"] == "B"))

    hs, hr = ca_seq(af, 394, 509)
    ds, dr = ca_seq(dk_a, 389, 520)
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -10, -0.5, "local"
    b = al.align(hs, ds)[0]
    ia, ib = [], []
    for (s1, e1), (s2, e2) in zip(b.aligned[0], b.aligned[1]):
        ia += list(range(s1, e1)); ib += list(range(s2, e2))
    pairs = [(hr[a], dr[c]) for a, c in zip(ia, ib)]
    idn = 100 * sum(hs[a] == ds[c] for a, c in zip(ia, ib)) / len(ia)

    P, Q = paired_bb(af, dk_a, pairs)
    R, t = kabsch(P, Q)
    rms = rmsd(P, Q, R, t)
    print(f"SBDβ 比对 HSPA8(394-509)↔DnaK：{len(pairs)} 残基，一致性 {idn:.1f}%，"
          f"主链叠合 RMSD {rms:.2f} A")

    afx = apply(af["xyz"], R, t)
    dk_prot = subset(dk_a, dk_a["elem"] != "H")
    pep = subset(dk_b, dk_b["elem"] != "H")

    # 检查叠合质量：变换后 AF-HSPA8 的 SBDβ 到 DnaK SBDβ 的距离
    # （已在上面打印 RMSD）

    print("\n[对照] 1DKX 天然肽 NRLLLTG ↔ DnaK SBD：", clash(pep["xyz"], dk_prot["xyz"]))

    g = read_gro(ROOT / "md_runs" / "PFNA" / "md.gro")
    scr = subset(g, (g["resn"] != "PFNA") & (g["elem"] != "H"))
    pfna = g["xyz"][(g["resn"] == "PFNA") & (g["elem"] != "H")]

    def dom(lo, hi):
        return afx[(af["resi"] >= lo) & (af["resi"] <= hi) & (af["elem"] != "H")]
    sbd_b, sbd_a, nbd = dom(394, 509), dom(510, 616), dom(1, 383)

    pep_res = sorted(set(pep["resi"]))
    pep_seq = "".join(AA3.get(pep["resn"][pep["resi"] == r][0], "X") for r in pep_res)
    print(f"\n参照肽段：{pep_seq}（残基编号 {pep_res}）")

    MOTIFS = {"QVDRR": (351, 355), "QDPVR": (339, 343)}
    rep = dict(sbd_fit_rmsd=rms, n_pairs=len(pairs), identity_pct=idn,
               control_pep_vs_dnak=clash(pep["xyz"], dk_prot["xyz"]),
               placements={})

    for tag, (lo, hi) in MOTIFS.items():
        print(f"\n=== {tag} {lo}-{hi} ===")
        rows = []
        for start in range(0, len(pep_res) - 4):
            win = pep_res[start:start + 5]
            P, Q = paired_bb(scr, pep, list(zip(range(lo, hi + 1), win)))
            if len(P) < 15:
                continue
            R, t = kabsch(P, Q)
            r = rmsd(P, Q, R, t)
            X = apply(scr["xyz"], R, t)
            c_b, c_a = clash(X, sbd_b), clash(X, sbd_a)
            c_n = clash(X, nbd)
            c_l = clash(apply(pfna, R, t), sbd_b)
            rows.append(dict(register=[int(win[0]), int(win[-1])],
                             motif_rmsd=round(r, 2), vs_SBDbeta=c_b,
                             vs_SBDalpha=c_a, vs_NBD=c_n,
                             PFNA_vs_SBDbeta=c_l))
            print(f"  register {win[0]}-{win[-1]}（基序 RMSD {r:.2f} A）: "
                  f"SBDβ severe {c_b['severe']:5d} bad {c_b['bad']:5d} <3A {c_b['contact']:5d}"
                  f" 最近 {c_b['min']:.2f} A | PFNA↔SBDβ 最近 {c_l['min']:.2f} A"
                  f"（severe {c_l['severe']}）")
        rows.sort(key=lambda x: (x["vs_SBDbeta"]["severe"], -x["vs_SBDbeta"]["min"]))
        rep["placements"][tag] = rows
        if rows:
            z = rows[0]
            print(f"  --> 最有利 register {z['register']}：severe {z['vs_SBDbeta']['severe']}，"
                  f"bad {z['vs_SBDbeta']['bad']}，最近 {z['vs_SBDbeta']['min']:.2f} A；"
                  f"PFNA↔SBDβ 最近 {z['PFNA_vs_SBDbeta']['min']:.2f} A")

    OUT.mkdir(exist_ok=True)
    (OUT / "hspa8_insertion.json").write_text(json.dumps(rep, indent=1))
    print(f"\n-> {OUT / 'hspa8_insertion.json'}")


if __name__ == "__main__":
    main()
