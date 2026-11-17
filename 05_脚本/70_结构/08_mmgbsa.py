"""08 · MM-GBSA（单轨迹法，OBC2 隐式溶剂，无熵项）

必须用 ChimeraX 自带的 Python 运行（本机唯一带 openmm / parmed 的解释器）：
    D:\\Tools\\ChimeraX\\bin\\python.exe scripts\\08_mmgbsa.py T3 PFNA

公式：
    dG_bind = <dE_vdW> + <dE_elec> + <dG_GB> + <dG_np>
  * dE_vdW / dE_elec：复合物与两个单体之间的**分子间**作用能，直接对配体-蛋白原子对求和，
    无截断（MM-GBSA 气相项惯例）。分子内项在同一帧、同一参数下完全抵消，故不需要真的
    拆三个体系去算。
  * dG_GB = G_pol(复合物) - G_pol(受体) - G_pol(配体)，用 OpenMM GBSAOBCForce（OBC-II，
    igb=5），ε_sol=78.3，ε_solute=1。半径取 **mbondi2**（igb=5 的配套集），筛选因子 0.8。
  * dG_np = -gamma * dSASA_buried，gamma=0.0072 kcal/mol/A^2，SASA 用 Shrake-Rupley
    (probe 1.4 A)；由于埋藏面积只涉及界面原子，只对 10 A 内的界面原子求和（等价且快）。
  * 不计算 -TdS（构象熵），因此结果是"有效结合自由能"，用于**相对比较**而非绝对预测。
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import openmm
from openmm import Platform, System, Context, VerletIntegrator, GBSAOBCForce, unit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import mdlib  # noqa: E402

RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.50, "F": 1.50,
         "P": 1.85, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98}
SCALE = 0.8
GAMMA = 0.0072         # kcal/mol/A^2
KCOUL = 332.0637       # kcal*A/mol/e^2
KJ2KCAL = 0.239006
EPS_SOLV = 78.3
IFACE_CUT = 10.0       # A，界面原子判定
N_SPHERE = 92


def mbondi2_radii(elem, bonds):
    """Amber mbondi2 半径集——igb=5 (OBC-II) 的配套半径。

    必须用它，不能用纯 Bondi：mbondi2 把与氧相连的氢缩到 0.80 A、与氮相连的 1.30 A。
    T3 有两个酚羟基、一个铵基，PFNA 一个都没有，用 Bondi 会让两者去溶剂化惩罚
    出现不对称偏差，直接污染 ΔG_GB 的可比性。
    """
    from collections import defaultdict
    rad = np.array([RADII.get(str(e), 1.70) for e in elem], dtype=float)
    nbr = defaultdict(list)
    for a, b in bonds:
        nbr[int(a)].append(int(b))
        nbr[int(b)].append(int(a))
    for i, e in enumerate(elem):
        if str(e) != "H":
            continue
        heavy = [str(elem[j]) for j in nbr[i] if str(elem[j]) != "H"]
        if not heavy:
            continue
        rad[i] = 0.80 if heavy[0] == "O" else (1.30 if heavy[0] == "N" else 1.20)
    return rad


def sphere_points(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi), np.cos(phi)], 1)


SPH = sphere_points(N_SPHERE)


def mm_interaction(xyz_a, xyz_b, sig_a, eps_a, q_a, sig_b, eps_b, q_b):
    """配体-蛋白分子间 LJ + 库仑（kcal/mol），无截断。坐标单位 A。"""
    r = np.linalg.norm(xyz_a[:, None, :] - xyz_b[None, :, :], axis=2)
    inv = 1.0 / np.maximum(r, 1e-6)
    sig = 0.5 * (sig_a[:, None] + sig_b[None, :])
    eps = np.sqrt(eps_a[:, None] * eps_b[None, :])
    s6 = (sig * inv) ** 6
    e_vdw = float((4.0 * eps * (s6 * s6 - s6)).sum())
    e_ele = float((KCOUL * np.outer(q_a, q_b) * inv).sum())
    return e_vdw, e_ele


def sasa_subset(xyz_all, rad_all, subset, probe=1.4):
    """只对 subset 里的原子算 SASA，但遮挡由全部原子提供。坐标 A。"""
    tree = cKDTree(xyz_all)
    rmax = rad_all.max()
    total = 0.0
    for i in subset:
        ri = rad_all[i] + probe
        nb = np.array(tree.query_ball_point(xyz_all[i], ri + rmax), dtype=int)
        nb = nb[(nb != i)]
        pts = xyz_all[i] + ri * SPH
        if len(nb):
            d = np.linalg.norm(pts[:, None, :] - xyz_all[nb][None, :, :], axis=2)
            exposed = (d >= (rad_all[nb] + probe)[None, :] - 1e-9).all(1)
        else:
            exposed = np.ones(len(pts), bool)
        total += exposed.sum() / len(pts) * 4.0 * np.pi * ri * ri
    return total


def gb_energy(xyz_nm, q, rad_nm, mass, platform):
    sysm = System()
    for m in mass:
        sysm.addParticle(float(m))
    gb = GBSAOBCForce()
    gb.setSolventDielectric(EPS_SOLV)
    gb.setSoluteDielectric(1.0)
    for i in range(len(q)):
        gb.addParticle(float(q[i]), float(rad_nm[i]), SCALE)
    sysm.addForce(gb)
    ctx = Context(sysm, VerletIntegrator(0.001), platform)
    ctx.setPositions(xyz_nm)
    st = ctx.getState(getEnergy=True)
    return st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def pick_platform():
    """默认 CPU（本机 GB 内核已优化，CUDA 缺 cl.exe 无法 JIT）。

    环境变量 MMGBSA_PLATFORM=CUDA 可强制走 GPU（Colab 上必须这样设）。
    """
    import os
    want = os.environ.get("MMGBSA_PLATFORM")
    if want:
        order = [want, "CUDA", "OpenCL", "CPU"]
    else:
        order = ["CPU", "OpenCL", "CUDA"]
    for nm in order:
        try:
            p = Platform.getPlatformByName(nm)
            if p:
                print(f"OpenMM 平台: {nm}")
                return p
        except Exception:
            continue
    raise SystemExit("no OpenMM platform available")


def run(name, stride=5):
    work = ROOT / "md_runs" / name
    top = mdlib.read_top(work / "processed.top")
    frames, _ = mdlib.load_frames(work)
    nf = frames.shape[0]
    sel = list(range(0, nf, stride))
    print(f"{name}: {nf} 帧，取 {len(sel)} 帧做 MM-GBSA")

    elem = top["elem"]
    n_lig = _n_lig(top)
    lig = np.arange(len(elem) - n_lig, len(elem))
    prot = np.arange(0, len(elem) - n_lig)
    rad = mbondi2_radii(elem, top["bonds"])                     # A（mbondi2，配 igb=5）
    sig = top["sigma"] * 10.0                                   # nm -> A
    eps = top["eps"] * KJ2KCAL                                  # kJ -> kcal
    q = top["charge"]

    platform = pick_platform()
    rec = []
    d_sasa = 0.0
    t0 = time.time()
    for k, s in enumerate(sel):
        x = frames[s] * 10.0                                    # nm -> A
        xl, xp = x[lig], x[prot]
        e_vdw, e_ele = mm_interaction(xl, xp, sig[lig], eps[lig], q[lig],
                                      sig[prot], eps[prot], q[prot])
        x_nm = x / 10.0
        g_all = gb_energy(x_nm, q, rad / 10.0, top["mass"], platform)
        g_p = gb_energy(x_nm[prot], q[prot], rad[prot] / 10.0, top["mass"][prot], platform)
        g_l = gb_energy(x_nm[lig], q[lig], rad[lig] / 10.0, top["mass"][lig], platform)
        d_gb = (g_all - g_p - g_l) * KJ2KCAL

        if True:                                                # SASA 每帧都算（很便宜）
            dl = np.linalg.norm(xl[:, None, :] - xp[None, :, :], axis=2)
            iface_p = np.where((dl < IFACE_CUT).any(0))[0]
            iface_l = np.where((dl < IFACE_CUT).any(1))[0]
            s_com = sasa_subset(x, rad, np.concatenate([iface_p, iface_l + len(prot)]))
            s_rec = sasa_subset(xp, rad[prot], iface_p)
            s_lig = sasa_subset(xl, rad[lig], iface_l)
            d_sasa = s_rec + s_lig - s_com
        rec.append(dict(frame=int(s), e_vdw=e_vdw, e_ele=e_ele, g_gb=d_gb,
                        d_sasa=float(d_sasa), g_np=float(-GAMMA * d_sasa)))
        if k % 10 == 0:
            print(f"  {k+1}/{len(sel)}  dE_vdw={e_vdw:8.2f} dE_ele={e_ele:8.2f} "
                  f"dG_gb={d_gb:8.2f} dG_np={-GAMMA*d_sasa:6.2f}  "
                  f"({time.time()-t0:.0f}s)")

    arr = {k: np.array([r[k] for r in rec], dtype=float)
           for k in ("e_vdw", "e_ele", "g_gb", "d_sasa", "g_np")}
    arr["total"] = arr["e_vdw"] + arr["e_ele"] + arr["g_gb"] + arr["g_np"]
    summ = {k: dict(mean=float(v.mean()), sem=float(v.std(ddof=1) / np.sqrt(len(v))),
                    std=float(v.std(ddof=1)))
            for k, v in arr.items()}
    # 分块均值：帧间强相关 + 位点漂移时 SEM 会严重低估误差，必须看 block 一致性
    nb = 5
    cut = np.linspace(0, len(arr["total"]), nb + 1).astype(int)
    blocks = [float(arr["total"][cut[i]:cut[i + 1]].mean()) for i in range(nb)]
    summ["total"]["blocks"] = blocks
    summ["total"]["block_range"] = [float(min(blocks)), float(max(blocks))]
    out = dict(system=name, n_frames=len(sel), stride=stride,
               frames=rec, summary=summ)
    (ROOT / "results" / f"mmgbsa_{name}.json").write_text(json.dumps(out, indent=1))
    print(f"\n{name}  MM-GBSA (kcal/mol, n={len(sel)})")
    for k in ("e_vdw", "e_ele", "g_gb", "g_np", "total"):
        print(f"   {k:8s} {summ[k]['mean']:9.2f} ± {summ[k]['sem']:.2f}")
    print("   分块均值(每 %.1f ns): %s" % (len(sel) * stride * 10 / 1000 / nb,
                                        "  ".join(f"{b:.2f}" for b in blocks)))
    return out


def _n_lig(top):
    """配体原子数：resname 与蛋白不同的最后一段。"""
    rn = top["resname"]
    i = len(rn) - 1
    while i > 0 and rn[i] == rn[-1]:
        i -= 1
    return len(rn) - i - 1


if __name__ == "__main__":
    for nm in (sys.argv[1:] or ["T3", "PFNA"]):
        run(nm)
