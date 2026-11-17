"""02 · 配体拓扑：RDKit 建结构 -> ChimeraX 赋 AM1-BCC 电荷与 GAFF 类型 -> parmchk2 取全参数 -> GROMACS itp

为什么绕这么远：
  * antechamber 是 cygwin 程序，内部用 `sh -c "rm ..."` 清理中间文件；本机对 rm 的
    拦截会让它在"键入"环节中断，输出的芳香碳全被误判成 sp3(c3)。不能直接用。
  * acpype 需要 tleap，本机 ChimeraX 自带的 amber20 是精简版，没有 tleap。
  * ChimeraX 自带同一套 GAFF 类型规则 + sqm(AM1-BCC)，且是原生 Windows 程序，
    不受上述 shell 问题影响；parmchk2 仍可用，用来取完整的 BOND/ANGLE/DIHE/NONBON。

产物：ligands/<NAME>.itp（含 [atomtypes]，类型名加 G 前缀避免与 amber99sb 冲突）
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
LIG = ROOT / "ligands"
AMBER = Path("D:/Tools/ChimeraX/bin/amber20")
CX = Path("D:/Tools/ChimeraX/bin/ChimeraX.exe")

# 净电荷：PFNA 生理 pH 下完全去质子化(-1)
SPECIES = {
    "PFNA": ("[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", -1),
}

MASS = {"C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
        "I": 126.904, "H": 1.008, "S": 32.06, "P": 30.974}


# ---------------------------------------------------------------- 1. 写干净 PDB
def write_pdb(mol, path, resname="UNL"):
    lines = []
    conf = mol.GetConformer()
    for i, a in enumerate(mol.GetAtoms(), 1):
        p = conf.GetAtomPosition(i - 1)
        el = a.GetSymbol()
        name = f"{el}{i}"[:4]
        chg = ""
        if a.GetFormalCharge() == 1:
            chg = "1+"
        elif a.GetFormalCharge() == -1:
            chg = "1-"
        lines.append(
            f"HETATM{i:5d} {name:<4s} {resname:<3s} A   1    "
            f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {el:>2s}{chg:>2s}")
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n", encoding="ascii")


# ------------------------------------------------- 2. ChimeraX: AM1-BCC + GAFF
CX_SCRIPT = r'''
import json, os
from chimerax.core.commands import run
inp, outp = os.environ["LIG_IN"], os.environ["LIG_OUT"]
resname, netq = os.environ["LIG_RES"], os.environ["LIG_Q"]
run(session, 'open "%s"' % inp)
run(session, "addcharge nonstd #1 %s %s method am1-bcc" % (resname, netq))
m = session.models.list()[0]
atoms = [dict(serial=a.serial_number, name=a.name, element=a.element.name,
              gaff=getattr(a, "gaff_type", None),
              charge=float(a.charge) if a.charge is not None else 0.0,
              xyz=[round(float(x), 4) for x in a.coord]) for a in m.atoms]
bonds = [[b.atoms[0].serial_number, b.atoms[1].serial_number,
          int(getattr(b, "order", 1) or 1)] for b in m.bonds]
data = dict(resname=m.residues[0].name, atoms=atoms, bonds=bonds,
            qtot=sum(a["charge"] for a in atoms))
json.dump(data, open(outp, "w"), indent=1)
print("WROTE", len(atoms), round(data["qtot"], 4))
'''


def run_chimerax(inp, outp, resname, netq):
    sc = ROOT / "scripts" / "_cx_charge.py"
    sc.write_text(CX_SCRIPT, encoding="utf-8")
    env = dict(os.environ, LIG_IN=str(inp), LIG_OUT=str(outp),
               LIG_RES=resname, LIG_Q=str(netq))
    # ChimeraX 打不开深层目录的绝对路径脚本（"Error opening python file"），
    # 但 cwd + 相对路径可以，所以统一在项目根目录下用相对路径调用
    r = subprocess.run([str(CX), "--nogui", "--exit", "--script", "scripts/_cx_charge.py"],
                       cwd=str(ROOT), env=env, capture_output=True, text=True)
    out = (r.stdout or "") + (r.stderr or "")
    if not Path(outp).exists():
        tail = "\n".join(out.strip().splitlines()[-12:])
        raise SystemExit(f"ChimeraX addcharge failed:\n{tail}")


# ---------------------------------------------------- 3. 类型兜底（备用规则）
def fallback_types(mol):
    """若 ChimeraX 未给出 gaff_type，用规则从 RDKit 化学环境推断（只覆盖本研究所用化学）。"""
    out = []
    for a in mol.GetAtoms():
        s, hyb, arom = a.GetSymbol(), str(a.GetHybridization()), a.GetIsAromatic()
        nb = a.GetNeighbors()
        if s == "H":
            n = nb[0]
            if n.GetSymbol() == "O":
                t = "ho"
            elif n.GetSymbol() == "N":
                t = "hn"
            elif n.GetIsAromatic() or str(n.GetHybridization()) == "SP2":
                t = "ha"
            else:
                elec = sum(1 for x in n.GetNeighbors()
                           if x.GetSymbol() in ("N", "O", "F", "Cl", "Br", "I", "S"))
                t = "hc" if elec == 0 else ("h1" if elec == 1 else "h2")
        elif s == "C":
            if arom:
                t = "ca"
            elif hyb == "SP2":
                nO = sum(1 for x in nb if x.GetSymbol() == "O")
                t = "c" if nO >= 2 else "c2"
            else:
                t = "c3"
        elif s == "O":
            if hyb == "SP2" and len(nb) == 1:
                t = "o"
            elif all(x.GetIsAromatic() for x in nb) and len(nb) == 2:
                t = "os"
            elif len(nb) == 2 and any(x.GetIsAromatic() for x in nb) and \
                    any(x.GetSymbol() == "C" and not x.GetIsAromatic() for x in nb):
                t = "os"
            elif any(x.GetSymbol() == "H" for x in nb):
                t = "oh"
            else:
                t = "o"
            if len(nb) == 2 and nb[0].GetIsAromatic() and nb[1].GetIsAromatic() and \
                    all(x.GetSymbol() == "C" for x in nb) and a.GetTotalNumHs() == 0 and \
                    a.GetFormalCharge() == 0 and hyb == "SP3":
                t = "os"
        elif s == "N":
            t = "n4" if len(nb) + a.GetTotalNumHs() >= 4 else ("n3" if hyb == "SP3" else "n")
        elif s == "F":
            t = "f"
        elif s == "I":
            t = "i"
        elif s == "S":
            t = "s"
        else:
            t = s.lower()
        out.append(t)
    return out


# ------------------------------------------------------- 4. 解析 frcmod & 输出
def parse_frcmod(path):
    """把 frcmod 解析成 {BOND:[(key,vals)], ...}。

    frcmod 的空格很随意：`o -c3-os  113.500  122.300` 这种键名里带空格，
    而 DIHE 是 4 个数(idivf,K,phase,path)、IMPROPER 只有 3 个。所以先把键名内的
    空格压掉，再取第一段连续数字。
    """
    rec = {"BOND": [], "ANGLE": [], "DIHE": [], "IMPROPER": [], "NONBON": [], "MASS": []}
    sec = None
    for ln in Path(path).read_text(errors="ignore").splitlines():
        s = ln.strip()
        if not s or s.startswith("Remark"):
            continue
        up = s.upper()
        if up in ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON", "END"):
            sec = up
            continue
        if sec is None:
            continue
        flat = re.sub(r"\s*-\s*", "-", s)
        tok = flat.split()
        if len(tok) < 2:
            continue
        k = tuple(tok[0].split("-"))
        nums = []
        for t in tok[1:]:
            try:
                nums.append(float(t))
            except ValueError:
                break
        if nums:
            rec[sec].append((k, nums))
    return rec


def assign_2(rec, ta, tb):
    for k, v in rec["BOND"]:
        if len(k) == 2 and ((k[0], k[1]) == (ta, tb) or (k[0], k[1]) == (tb, ta)):
            return v[0], v[1]
    return None


def assign_3(rec, ta, tb, tc):
    for k, v in rec["ANGLE"]:
        if len(k) == 3 and (k == (ta, tb, tc) or k == (tc, tb, ta)):
            return v[0], v[1]
    return None


def _match(k, t):
    """键名与四体元组匹配，GT 通配符 X 可匹配任意类型"""
    if len(k) != 4:
        return False
    return all(kk in ("X", "x") or kk == cc for kk, cc in zip(k, t))


def assign_4(rec, section, t):
    """返回该四体对应的所有条目（同一键可能有多重大小不同的项）"""
    outs = []
    for k, v in rec[section]:
        if _match(k, t) or _match(k, t[::-1]):
            if section == "DIHE":
                if len(v) >= 4:
                    outs.append((int(v[0]), v[1], v[2], int(round(v[3]))))
            else:
                if len(v) >= 3:
                    outs.append((1, v[0], v[1], int(round(v[2]))))
    return outs


def main():
    for name, (smi, netq) in SPECIES.items():
        print(f"=== {name} (net {netq}) ===")
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 20261111
        AllChem.EmbedMolecule(m, ps)
        AllChem.MMFFOptimizeMolecule(m, maxIters=5000)
        pdb = LIG / f"{name}.pdb"
        write_pdb(m, pdb)

        js = LIG / f"{name}_cx.json"
        if not js.exists():
            run_chimerax(pdb, js, "UNL", netq)
        d = json.loads(js.read_text())
        atoms = d["atoms"]
        print(f"   ChimeraX: {len(atoms)} 原子, 总电荷 {d['qtot']:+.4f}")

        # 类型：优先 ChimeraX 的 gaff_type；缺失则用规则兜底
        types = [a["gaff"] for a in atoms]
        n_missing = sum(1 for t in types if not t)
        if n_missing:
            fb = fallback_types(m)
            types = [t if t else fb[i] for i, t in enumerate(types)]
            print(f"   {n_missing} 个原子类型缺失，用规则补齐")

        # 写回带 GAFF 类型的 mol2，交给 parmchk2 取完整参数
        mol2 = LIG / f"{name}_typed.mol2"
        with open(mol2, "w", encoding="ascii") as f:
            f.write(f"@<TRIPOS>MOLECULE\n{name}\n{len(atoms)} {len(d['bonds'])} 1 0 0\n"
                    "SMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
            for a, t in zip(atoms, types):
                f.write(f"{a['serial']:7d} {a['name']:<6s}{a['xyz'][0]:10.4f}{a['xyz'][1]:10.4f}"
                        f"{a['xyz'][2]:10.4f} {t:<6s} 1 {name:<4s}{a['charge']:9.4f}\n")
            f.write("@<TRIPOS>BOND\n")
            for i, b in enumerate(d["bonds"], 1):
                f.write(f"{i:6d}{b[0]:6d}{b[1]:6d} {b[2] if b[2] in (1, 2, 3) else 1}\n")
            f.write("@<TRIPOS>SUBSTRUCTURE\n1 %s 1\n" % name)

        frc = LIG / f"{name}_full.frcmod"
        env = dict(os.environ, AMBERHOME=AMBER.as_posix())
        r = subprocess.run([(AMBER / "bin" / "parmchk2.exe").as_posix(),
                            "-i", mol2.name, "-f", "mol2", "-o", frc.name,
                            "-s", "2", "-a", "Y"],
                           cwd=LIG, env=env, capture_output=True, text=True)
        if r.returncode != 0 or not frc.exists():
            raise SystemExit("parmchk2 failed\n" + (r.stdout or "")[-1500:] + (r.stderr or "")[-800:])
        par = parse_frcmod(frc)
        print(f"   parmchk2: BOND {len(par['BOND'])} / ANGLE {len(par['ANGLE'])} / "
              f"DIHE {len(par['DIHE'])} / IMPROPER {len(par['IMPROPER'])} / NONBON {len(par['NONBON'])}")

        # ---------------- 生成 GROMACS itp ----------------
        nb = {}
        for k, v in par["NONBON"]:
            if len(k) == 1 and len(v) >= 2:
                nb[k[0]] = (v[0], v[1])                # Rmin/2 (Å), epsilon (kcal/mol)

        used = sorted(set(types))
        mass = {a["element"]: MASS.get(a["element"], 12.011) for a in atoms}
        el_of_type = {}
        for a, t in zip(atoms, types):
            el_of_type.setdefault(t, a["element"])

        lines = [f"; {name} - GAFF2 parameters (AM1-BCC charges), converted from parmchk2 full set",
                 "[ atomtypes ]"]
        R2S = 2.0 / 2.0 ** (1.0 / 6.0)   # sigma = 2*R*/2^(1/6)（与 GROMACS amber 端口一致）
        for t in used:
            if t not in nb:
                raise SystemExit(f"类型 {t} 缺 NONBON 参数")
            rmin2, eps = nb[t]
            lines.append(f"G{t:<6s} {el_of_type[t]:>2s} {mass[el_of_type[t]]:8.3f} 0.000  A "
                         f"{rmin2 * R2S / 10.0:9.5f} {eps * 4.184:9.5f}")
        lines += ["", f"[ moleculetype ]", f"{name}   3", "", "[ atoms ]"]
        for i, (a, t) in enumerate(zip(atoms, types), 1):
            lines.append(f"{i:5d} G{t:<7s} 1  {name:<4s} {a['name']:<5s} {i:5d} "
                         f"{a['charge']:9.4f} {MASS.get(a['element'], 12.011):8.3f}")

        # 键/角/二面角：按连通性枚举，从 frcmod 取参数
        adj = {i: set() for i in range(1, len(atoms) + 1)}
        for b in d["bonds"]:
            adj[b[0]].add(b[1]); adj[b[1]].add(b[0])
        tmap = {a["serial"]: t for a, t in zip(atoms, types)}

        # 注意 GROMACS 的列序与 frcmod / prmtop 不同，务必按 GROMACS 手册写：
        #   [ bonds ]      ai aj funct  b0(nm)   kb(kJ/mol/nm^2)      <- frcmod 是 (K, b0)
        #   [ angles ]     ai aj ak funct  th0(deg)  cth(kJ/mol/rad^2) <- frcmod 是 (K, th0)
        #   [ dihedrals ]  funct 9:  phi0(deg)  kphi  n              <- frcmod 是 (idivf, Vn, phase, n)
        #   [ dihedrals ]  funct 4:  phi0(deg)  kphi  n              <- frcmod 是 (Vn, phase, n=2)
        miss = [0, 0, 0]
        lines += ["", "[ bonds ]", ";  ai    aj  funct       b0(nm)          kb"]
        for i in range(1, len(atoms) + 1):
            for j in sorted(x for x in adj[i] if x > i):
                v = assign_2(par, tmap[i], tmap[j])
                if v is None:
                    miss[0] += 1; continue
                lines.append(f"{i:5d} {j:5d} 1 {v[1] / 10.0:12.5f} {v[0] * 4.184 * 100:13.3f}")

        lines += ["", "[ angles ]",
                  ";  ai    aj    ak  funct      th0(deg)          cth"]
        for j in range(1, len(atoms) + 1):
            nbrs = sorted(adj[j])
            for ai in range(len(nbrs)):
                for bi in range(ai + 1, len(nbrs)):
                    i, k = nbrs[ai], nbrs[bi]
                    v = assign_3(par, tmap[i], tmap[j], tmap[k])
                    if v is None:
                        miss[1] += 1; continue
                    lines.append(f"{i:5d} {j:5d} {k:5d} 1 "
                                 f"{v[1]:12.3f} {v[0] * 2 * 4.184:13.3f}")

        # 1-4 对：只写 3 列，让 grompp 按 [ atomtypes ] 组合规则 + fudgeLJ 自动填 c6/c12
        pairs = set()
        for j in range(1, len(atoms) + 1):
            for k in sorted(x for x in adj[j] if x > j):
                for i in sorted(adj[j] - {k}):
                    for l in sorted(adj[k] - {j}):
                        pairs.add((min(i, l), max(i, l)))
        pairs = sorted(pairs)
        lines += ["", "[ pairs ]", ";  ai    aj  funct"]
        for i, l in pairs:
            lines.append(f"{i:5d} {l:5d} 1")

        lines += ["", "[ dihedrals ]   ; type 9 = proper periodic",
                  ";  ai    aj    ak    al  funct      phi0(deg)        kphi   n"]
        for j in range(1, len(atoms) + 1):
            for k in sorted(x for x in adj[j] if x > j):
                for i in sorted(adj[j] - {k}):
                    for l in sorted(adj[k] - {j}):
                        v = assign_4(par, "DIHE", (tmap[i], tmap[j], tmap[k], tmap[l]))
                        if not v:
                            miss[2] += 1; continue
                        for idivf, K, ph, path in v:
                            lines.append(f"{i:5d} {j:5d} {k:5d} {l:5d} 9 "
                                         f"{ph:12.3f} {K / (2.0 * idivf) * 4.184:11.4f} "
                                         f"{abs(path):3d}")

        lines += ["", "[ dihedrals ]   ; type 4 = improper periodic, central atom is 3rd",
                  ";  ai    aj    ak    al  funct      phi0(deg)        kphi   n"]
        n_imp = 0
        for j in range(1, len(atoms) + 1):
            if len(adj[j]) != 3:
                continue
            a_ = m.GetAtomWithIdx(j - 1)
            if not (a_.GetIsAromatic() or str(a_.GetHybridization()) == "SP2"):
                continue
            i, k, l = sorted(adj[j])
            v = assign_4(par, "IMPROPER", (tmap[i], tmap[j], tmap[k], tmap[l]))
            if not v:
                continue
            for _idivf, K, ph, path in v:
                # Amber 约定中心原子在第 3 位 -> GROMACS 也按 i-j-k-l 顺序算二面角
                lines.append(f"{i:5d} {k:5d} {j:5d} {l:5d} 4 "
                             f"{ph:12.3f} {K / 2.0 * 4.184:11.4f} {abs(path):3d}")
                n_imp += 1
                break

        itp = LIG / f"{name}.itp"
        itp.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"   缺失参数: bond {miss[0]}, angle {miss[1]}, dihedral {miss[2]}；improper {n_imp} 条")
        print(f"   -> {itp}  ({len(lines)} 行, 电荷和 {sum(a['charge'] for a in atoms):+.4f})")
        sanity(itp)


def sanity(itp):
    """列序自检：把 itp 读回来，检查每个参数的物理量级是否合理。

    这一关专门用来抓「b0 与 kb 写反」这类静默错误——写反时 b0 会变成 1e5 量级、
    kb 变成 0.1 量级，GROMACS 不会报错，但 EM 会直接炸到 1e10 kJ/mol。
    """
    sec = None
    bad = []
    for ln in Path(itp).read_text(encoding="utf-8").splitlines():
        s = ln.split(";")[0].strip()
        if not s:
            continue
        if s.startswith("["):
            sec = s.strip("[] ").lower().split()[0]
            kind = s
            continue
        if sec not in ("bonds", "angles", "dihedrals"):
            continue
        t = s.split()
        if "type 4" in kind or "type 9" in kind:
            pass
        if sec == "bonds" and len(t) >= 5:
            b0, kb = float(t[3]), float(t[4])
            if not (0.08 <= b0 <= 0.30 and 1e4 <= kb <= 1e7):
                bad.append(("bond", b0, kb))
        elif sec == "angles" and len(t) >= 6:
            th0, cth = float(t[4]), float(t[5])
            if not (60.0 <= th0 <= 180.0 and 10.0 <= cth <= 5000.0):
                bad.append(("angle", th0, cth))
        elif sec == "dihedrals" and len(t) >= 8:
            phi0, kphi = float(t[5]), float(t[6])
            if not (abs(phi0) <= 360.0 and 0.0 <= kphi <= 100.0):
                bad.append(("dihedral", phi0, kphi))
    if bad:
        raise SystemExit(f"!! {itp.name} 参数列序异常 {len(bad)} 条，例如 {bad[:3]}")
    print(f"   自检通过：键/角/二面角参数量级正常")


if __name__ == "__main__":
    main()
