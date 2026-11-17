"""05 · 建 MD 体系：把对接位姿放回蛋白，做周期盒、加溶剂与离子、能量最小化

关键处理：
  * 配体 itp 里的 [atomtypes] 单独拆成 LIG_atomtypes.itp，放在力场 include 之后、
    任何 [moleculetype] 之前（GROMACS 要求原子类型先定义）
  * 原子类型名统一加 G 前缀，避免与 amber99sb-ildn 的 c3/ca 等同名冲突
  * 平衡段给配体 400 kJ/mol/nm^2 的位置限制，生产段完全放开
  * 自定义索引 Solute(蛋白+配体) / Solvent(水+离子)，供温控分组
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
    """把对接结果（中心 + 四元数）施加到配体坐标上"""
    import pandas as pd
    d = pd.read_csv(RES / f"dock_{name}.csv")
    row = d.sort_values("energy").iloc[0]          # 取能量最低的位姿
    c = np.array([row.cx, row.cy, row.cz], dtype=float)
    q = np.array([float(v) for v in str(row.quat).split(",")])
    xyz, anames, resname, el = read_pdb(LIG / f"{name}.pdb")
    heavy = np.array([e != "H" for e in el])
    com0 = xyz[heavy].mean(0)              # 对接时的质心只算重原子
    R = Rotation.from_quat(q).as_matrix()
    new = (xyz - com0) @ R.T + c
    return new, anames, el, resname[0]


def relieve_clash(xyz, heavy, tree, target=2.2, step=0.15, max_iter=60):
    """把穿插进蛋白的配体整体推开。

    方向取「每个冲突配体原子 - 其最近蛋白原子」的平均单位向量（局部外法线），
    比用质心方向更稳：即使配体卡在沟槽里也不会沿切线乱滑。刚体平移，不破坏位姿。
    """
    xyz = xyz.copy()
    for _ in range(max_iter):
        d, idx = tree.query(xyz[heavy])
        if d.min() >= target:
            break
        near = tree.data[idx]
        vec = xyz[heavy] - near                       # 每个重原子相对最近蛋白原子的方向
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

    # 对接后自检：配体不能穿进蛋白（只比较重原子，H 已经由 pdb2gmx 加到蛋白上）
    prot_xyz, _, _, _ = read_pdb(ROOT / "structures" / "SCRN2_AF.pdb")
    from scipy.spatial import cKDTree
    tree = cKDTree(prot_xyz)
    dmin = tree.query(lig_xyz[heavy])[0].min()
    print(f"  配体 {nl} 原子，与蛋白最近重原子 {dmin:.2f} Å")
    if dmin < 2.2:                       # 对接打分允许 2.2 Å，实际位姿常在 1 Å 内穿插
        lig_xyz = relieve_clash(lig_xyz, heavy, tree, 2.2)
        dmin2 = tree.query(lig_xyz[heavy])[0].min()
        print(f"  穿模修整：{dmin:.2f} -> {dmin2:.2f} Å")

    # 蛋白 gro（pdb2gmx 已生成，含氢）
    title, np_, prot_lines = read_gro(PREP / "prot.gro")

    # 复合物 gro：蛋白在前，配体在后（与拓扑顺序一致）
    out = [f"{name} complex", f"{np_ + nl:5d}"] + prot_lines
    lig_nm = lig_xyz / 10.0          # gro 用 nm，配体坐标来自 PDB(Å)
    for k in range(nl):
        out.append(fmt_gro(1, ligres[:5], anames[k][:5], np_ + k + 1, lig_nm[k]))
    # 盒向量：pdb2gmx 写出的盒子可能无意义（AF 模型没有晶胞），
    # 这里直接按坐标范围 + 2 nm 给一个，交给 editconf -d 重新构建
    allx = np.vstack([np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])]
                                for l in prot_lines]), lig_nm])
    box = allx.max(0) - allx.min(0) + 2.0        # 单位 nm
    out.append("   " + "".join(f"{float(v):10.5f}" for v in box))
    (work / "complex.gro").write_text("\n".join(out) + "\n", encoding="ascii")

    # 配体拓扑拆分：atomtypes 单独一份
    itp = (LIG / f"{name}.itp").read_text(encoding="utf-8")
    m = re.search(r"\[ atomtypes \]\n(.*?)\n\n", itp, re.S)
    if not m:
        raise SystemExit("ligand itp lacks [ atomtypes ]")
    (work / "LIG_atomtypes.itp").write_text(
        "[ atomtypes ]\n" + m.group(1) + "\n", encoding="ascii")
    body = itp[:m.start()] + itp[m.end():]
    body = re.sub(r"^;.*\n", "", body, count=1)

    # 配体位置限制（只在平衡段生效）
    pr = ["[ position_restraints ]", "; atom  type  fx  fy  fz"]
    for i in range(1, nl + 1):
        pr.append(f"{i:6d} 1 400.0 400.0 400.0")
    (work / "LIG_posre.itp").write_text("\n".join(pr) + "\n", encoding="ascii")
    body += '\n#ifdef POSRES_LIG\n#include "LIG_posre.itp"\n#endif\n'
    (work / f"{name}.itp").write_text(body, encoding="utf-8")

    # 顶层拓扑
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

    # 盒 + 溶剂 + 离子
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

    # 索引：Solute = 蛋白+配体。GROMACS 的行缓冲有限，原子号必须折行（每行 <=15 个）
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
    print(f"  总原子 {ntot}（溶质 {nsolute}）")

    sh([GMX, "grompp", "-f", str(ROOT / "mdp" / "em.mdp"), "-c", "solv_ions.gro",
        "-r", "solv_ions.gro", "-p", "topol.top", "-n", "index.ndx",
        "-o", "em.tpr", "-pp", "processed.top", "-maxwarn", "2"], work)
    sh([GMX, "mdrun", "-deffnm", "em", "-ntmpi", "1"], work)
    tail = (work / "em.log").read_text(errors="ignore").splitlines()
    print("  " + [ln for ln in tail if "Maximum force" in ln or "Potential Energy" in ln][-2:].__str__())


if __name__ == "__main__":
    for nm in (sys.argv[1:] or ["T3", "PFNA"]):
        build_one(nm)
