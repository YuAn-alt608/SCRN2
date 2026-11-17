"""08 - MM-GBSA (single-trajectory, OBC2 implicit solvent, no entropy term)

Must be run with the Python bundled with ChimeraX (the only interpreter here with openmm / parmed):
    D:\\Tools\\ChimeraX\\bin\\python.exe scripts\\08_mmgbsa.py PFNA

Formulas:
    dG_bind = <dE_vdW> + <dE_elec> + <dG_GB> + <dG_np>
  * dE_vdW / dE_elec: the **intermolecular** energy between the complex and the two monomers,
    summed directly over ligand-protein atom pairs with no cutoff (the usual MM-GBSA gas-phase
    convention). Intramolecular terms cancel exactly within a frame at the same parameters, so
    there is no need to actually split the system three ways.
  * dG_GB = G_pol(complex) - G_pol(receptor) - G_pol(ligand) from OpenMM GBSAOBCForce
    (OBC-II, igb=5), with eps_sol=78.3 and eps_solute=1. Radii are **mbondi2** (the set that
    pairs with igb=5) and the screening factor is 0.8.
  * dG_np = -gamma * dSASA_buried with gamma=0.0072 kcal/mol/A^2; SASA by Shrake-Rupley
    (probe 1.4 A). Buried area only involves interface atoms, so the sum runs over interface
    atoms within 10 A, which is equivalent and faster.
  * -TdS (conformational entropy) is not computed, so the result is an effective binding free
    energy, useful for **relative comparison** rather than absolute prediction.
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
IFACE_CUT = 10.0       # A, interface-atom cutoff
N_SPHERE = 92


def mbondi2_radii(elem, bonds):
    """The Amber mbondi2 radius set: the radii that pair with igb=5 (OBC-II).

    It has to be used rather than plain Bondi: mbondi2 shrinks hydrogens bonded to oxygen to
    0.80 A and those bonded to nitrogen to 1.30 A. The two carboxylate oxygens of PFNA enter
    the electrostatic term directly, so Bondi radii would make the desolvation penalty
    inconsistent with the OBC-II model and contaminate dG_GB.
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
    """Intermolecular ligand-protein LJ plus Coulomb energy (kcal/mol), no cutoff. Coordinates in A."""
    r = np.linalg.norm(xyz_a[:, None, :] - xyz_b[None, :, :], axis=2)
    inv = 1.0 / np.maximum(r, 1e-6)
    sig = 0.5 * (sig_a[:, None] + sig_b[None, :])
    eps = np.sqrt(eps_a[:, None] * eps_b[None, :])
    s6 = (sig * inv) ** 6
    e_vdw = float((4.0 * eps * (s6 * s6 - s6)).sum())
    e_ele = float((KCOUL * np.outer(q_a, q_b) * inv).sum())
    return e_vdw, e_ele


def sasa_subset(xyz_all, rad_all, subset, probe=1.4):
    """SASA for the atoms in subset only, but every atom contributes occlusion. Coordinates in A."""
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
    """CPU by default (the GB kernel is optimised here, and CUDA cannot JIT without cl.exe).

    Setting MMGBSA_PLATFORM=CUDA forces the GPU (required on Colab).
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
                print(f"OpenMM platform: {nm}")
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
    print(f"{name}: {nf} frames, {len(sel)} used for MM-GBSA")

    elem = top["elem"]
    n_lig = _n_lig(top)
    lig = np.arange(len(elem) - n_lig, len(elem))
    prot = np.arange(0, len(elem) - n_lig)
    rad = mbondi2_radii(elem, top["bonds"])                     # A (mbondi2, as required by igb=5)
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

        if True:                                                # SASA on every frame (very cheap)
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
    # block means: frames are strongly correlated and the SEM badly understates the error when
    # the site drifts, so block-to-block agreement has to be inspected
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
    print("   block means (every %.1f ns): %s" % (len(sel) * stride * 10 / 1000 / nb,
                                        "  ".join(f"{b:.2f}" for b in blocks)))
    return out


def _n_lig(top):
    """Number of ligand atoms: the last block whose resname differs from the protein."""
    rn = top["resname"]
    i = len(rn) - 1
    while i > 0 and rn[i] == rn[-1]:
        i -= 1
    return len(rn) - i - 1


if __name__ == "__main__":
    for nm in (sys.argv[1:] or ["PFNA"]):
        run(nm)
