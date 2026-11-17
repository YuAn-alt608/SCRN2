"""07 · 六项指标分析（配体 RMSD / 接触数 / 氢键 / RMSF / 蛋白 RMSD / Rg）

指标定义（报告里可直接引用）：
  1. 配体 RMSD      ：蛋白主链叠合后，配体重原子相对首帧的 RMSD（只反映配体自身漂移）
  2. 结合位点接触数 ：配体重原子与蛋白重原子距离 < 0.6 nm 的原子对数
  3. 氢键           ：供体-受体距离 <= 0.35 nm 且 D-H...A 角 >= 150 度（即 H-D-A <= 30 度）
                      寿命 = 间歇自相关 C(tau)=<h(0)h(tau)>/<h(0)> 积分到 100 ps
  4. 配体 RMSF      ：同一叠合下每原子涨落
  5. 蛋白 RMSD      ：主链 (N,CA,C,O) RMSD
  6. 回转半径 Rg    ：蛋白质量加权
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import mdlib  # noqa: E402

GMX = "gmx"
CUT_CONTACT = 0.6      # nm
HB_DA = 0.35           # nm
HB_ANGLE = 30.0        # deg, H-D-A 角上限
LIFE_LAG = 100.0       # ps, 自相关积分上限


def sh(cmd, cwd, inp=None):
    r = subprocess.run(cmd, cwd=str(cwd), input=inp, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2500:]); print(r.stderr[-1500:])
        raise SystemExit(f"FAILED: {' '.join(map(str, cmd))}")
    return r


def write_ndx(groups, path):
    out = []
    for tag, ids in groups.items():
        out.append(f"[ {tag} ]")
        for i in range(0, len(ids), 15):
            out.append(" ".join(str(v) for v in ids[i:i + 15]))
    Path(path).write_text("\n".join(out) + "\n", encoding="ascii")


def read_gro_meta(path):
    """从 complex.gro 读每个原子的残基号 / 残基名（拓扑里没有残基编号）。"""
    lines = Path(path).read_text(errors="ignore").splitlines()
    n = int(lines[1])
    resid, resname = [], []
    for ln in lines[2:2 + n]:
        resid.append(int(ln[:5]))
        resname.append(ln[5:10].strip())
    return np.array(resid), np.array(resname)


def dump_traj(work, name, dt_ps=10):
    sh([GMX, "make_ndx", "-f", "md.tpr", "-o", "std.ndx"], work, inp="q\n")
    groups, cur = {}, None
    for ln in (work / "std.ndx").read_text().splitlines():
        if ln.startswith("["):
            cur = ln.strip("[] ").strip()
            groups[cur] = []
        elif cur and ln.split():
            groups[cur] += [int(v) for v in ln.split()]
    prot1 = groups["Protein"]
    lig1 = groups[name]
    write_ndx({"Complex": sorted(set(prot1) | set(lig1))}, work / "full.ndx")
    sh([GMX, "trjconv", "-s", "md.tpr", "-f", "md.xtc", "-n", "full.ndx",
        "-o", "traj.pdb", "-pbc", "mol", "-dt", str(dt_ps)], work, inp="Complex\n")
    return np.array(prot1, dtype=int) - 1, np.array(lig1, dtype=int) - 1


def hbond_sets(y, pairs):
    """单帧氢键判定，返回存在的 (donor, acceptor) 集合。

    判据与 gmx hbond 一致：D-A <= 0.35 nm，且 H 处的 D-H...A 角 >= 180-30 = 150 度
    （线性氢键为 180 度）。一个供体可能带多个氢，取最优的那个。
    """
    out = set()
    for dset, aset in pairs:
        if len(dset) == 0 or len(aset) == 0:
            continue
        d = np.linalg.norm(y[dset][:, None, :] - y[aset][None, :, :], axis=2)
        ii, jj = np.where(d <= HB_DA)
        for a_i, b_i in zip(ii, jj):
            dd, ac = dset[a_i], aset[b_i]
            best = 0.0
            for h in HYD.get(dd, []):
                v1 = y[dd] - y[h]        # H -> D
                v2 = y[ac] - y[h]        # H -> A
                c = v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
                best = max(best, np.degrees(np.arccos(np.clip(c, -1, 1))))
            if best >= 180.0 - HB_ANGLE:
                out.add((int(dd), int(ac)))
    return out


HYD = {}


def anal(name, dt_ps=10):
    work = ROOT / "md_runs" / name
    top = mdlib.read_top(work / "processed.top")
    resid, resname = read_gro_meta(work / "complex.gro")
    pidx, lidx = dump_traj(work, name, dt_ps)
    frames, boxes = mdlib.load_frames(work)
    nf = frames.shape[0]
    print(f"  {name}: {nf} 帧 × {frames.shape[1]} 原子")

    elem, aname = top["elem"], top["aname"]
    heavy_l = lidx[elem[lidx] != "H"]
    heavy_p = pidx[elem[pidx] != "H"]
    bb = pidx[np.isin(aname[pidx], ["N", "CA", "C", "O"])]
    ca = pidx[aname[pidx] == "CA"]

    # 供体 / 受体
    Hmask = elem == "H"
    for a, b in top["bonds"]:
        if Hmask[a]:
            HYD.setdefault(b, []).append(a)
        if Hmask[b]:
            HYD.setdefault(a, []).append(b)
    polar = np.isin(elem, ["N", "O", "S"])
    donor = np.array([i for i in np.where(polar)[0] if i in HYD], dtype=int)
    accept = np.array([i for i in np.where(polar)[0]
                       if elem[i] in ("O", "S") or (elem[i] == "N" and i not in HYD)], dtype=int)
    pairs = [
        (donor[np.isin(donor, pidx)], accept[np.isin(accept, lidx)]),
        (donor[np.isin(donor, lidx)], accept[np.isin(accept, pidx)]),
    ]

    ref = frames[0]
    prot_rmsd, lig_rmsd, rg, contacts, mindist = [], [], [], [], []
    fitted_p = np.empty((nf, len(pidx), 3), dtype=np.float32)
    lig_fit = np.empty((nf, len(lidx), 3), dtype=np.float32)
    hb_frames = []
    site_sets = []

    reslab = np.array([f"{resname[i]}{resid[i]}" for i in pidx])
    res_of_heavy = reslab[np.isin(pidx, heavy_p)]
    contact_freq = {}

    for s in range(nf):
        x = frames[s].copy()
        x[lidx] = mdlib.unwrap_ligand(x[lidx], x[pidx].mean(0), boxes[s])
        R, t = mdlib.fit_transform(x[bb], ref[bb])
        y = x @ R + t
        prot_rmsd.append(float(np.sqrt(((y[bb] - ref[bb]) ** 2).sum(1).mean())))
        lig_rmsd.append(float(np.sqrt(((y[heavy_l] - ref[heavy_l]) ** 2).sum(1).mean())))
        rg.append(float(np.sqrt((top["mass"][pidx] *
                                 ((y[pidx] - y[pidx].mean(0)) ** 2).sum(1)).sum()
                                / top["mass"][pidx].sum())))
        d = np.linalg.norm(y[heavy_l][:, None, :] - y[heavy_p][None, :, :], axis=2)
        contacts.append(int((d < CUT_CONTACT).sum()))
        mindist.append(float(d.min()))
        # 逐帧接触残基集合（用于位点身份漂移分析，Jaccard）
        hit = set(res_of_heavy[(d < CUT_CONTACT).any(0)].tolist())
        site_sets.append(hit)
        for r in hit:
            contact_freq[r] = contact_freq.get(r, 0) + 1
        fitted_p[s] = y[pidx]
        lig_fit[s] = y[lidx]
        hb_frames.append(hbond_sets(y, pairs))

    contact_freq = {k: v / nf for k, v in sorted(contact_freq.items(), key=lambda kv: -kv[1])}

    # 位点身份保持度：与首帧接触残基集合的 Jaccard 相似度
    s0 = site_sets[0]
    site_jac = [len(s & s0) / max(len(s | s0), 1) for s in site_sets]
    site_nres = [len(s) for s in site_sets]

    # 氢键寿命
    hb_keys = sorted({k for st in hb_frames for k in st})
    hb_detail = []
    for k in hb_keys:
        v = np.array([1.0 if k in st else 0.0 for st in hb_frames])
        occ = v.mean()
        if occ < 0.05:
            continue
        maxlag = min(int(LIFE_LAG / dt_ps), nf // 2)
        tau = 0.0
        for lag in range(maxlag):
            c = (v[:nf - lag] * v[lag:]).sum() / max(v.sum(), 1)
            tau += c * dt_ps
        dd, ac = k
        l_side = dd if dd in set(lidx.tolist()) else ac
        p_side = ac if l_side == dd else dd
        hb_detail.append(dict(
            ligand_atom=f"{aname[l_side]}({elem[l_side]})",
            protein_atom=f"{aname[p_side]}({elem[p_side]})",
            protein_res=f"{resname[p_side]}{resid[p_side]}",
            occupancy=float(occ), lifetime_ps=float(tau)))
    hb_detail.sort(key=lambda d: (-d["occupancy"], -d["lifetime_ps"]))

    # RMSF
    rmsf_p = np.sqrt(((fitted_p - fitted_p.mean(0)) ** 2).sum(-1).mean(0))
    rmsf_l = np.sqrt(((lig_fit - lig_fit.mean(0)) ** 2).sum(-1).mean(0))

    # 表面二维扩散：配体重原子质心的 MSD（叠合后坐标，消除整体平动/转动）
    com = lig_fit[:, np.isin(lidx, heavy_l)].mean(1).astype(np.float64)
    nl = max(nf // 2, 2)
    msd = np.array([np.mean(np.sum((com[lag:] - com[:-lag]) ** 2, 1))
                    for lag in range(1, nl)])
    lag_ns = np.arange(1, nl) * dt_ps / 1000.0
    fit = (lag_ns >= 0.1) & (lag_ns <= 1.0)
    d_surf = (np.polyfit(lag_ns[fit], msd[fit], 1)[0] / 4.0 * 1e-5
              if fit.sum() > 5 else float("nan"))     # nm^2/ns -> cm^2/s (2D: MSD=4Dt)

    t_ps = np.arange(nf) * dt_ps
    c0 = contacts[0]
    half = 0.5 * max(c0, 1)
    below = np.array(contacts) < half
    t_half = None
    for i in range(nf):
        if below[i:min(i + 10, nf)].all() and len(below[i:min(i + 10, nf)]) >= 5:
            t_half = float(t_ps[i])
            break
    res = dict(
        system=name, n_frames=int(nf), dt_ps=float(dt_ps), total_ns=float(nf * dt_ps / 1000),
        prot_rmsd=prot_rmsd, lig_rmsd=lig_rmsd, rg=rg,
        contacts=contacts, mindist=mindist,
        hbonds=[len(s) for s in hb_frames],
        rmsf_protein=rmsf_p.tolist(), rmsf_ligand=rmsf_l.tolist(),
        lig_atoms=[str(a) for a in aname[lidx]],
        lig_elem=[str(e) for e in elem[lidx]],
        ca_res=[str(f"{resname[i]}{resid[i]}") for i in ca],
        contact_freq=contact_freq,
        hbonds_detail=hb_detail[:10],
        site_jaccard=site_jac,
        site_nres=site_nres,
        msd_nm2=msd.tolist(),
        msd_lag_ns=lag_ns.tolist(),
        d_surface_cm2s=float(d_surf),
        c0=int(c0),
        bound_frac=float((np.array(contacts) >= half).mean()),
        t_half_ps=t_half,
    )
    (ROOT / "results" / f"metrics_{name}.json").write_text(json.dumps(res, indent=1))

    # 逐帧时间序列（便于外部核查，不必解 JSON）
    import csv
    with open(ROOT / "results" / f"timeseries_{name}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time_ns", "prot_rmsd_nm", "lig_rmsd_nm", "rg_nm",
                    "contacts", "mindist_nm", "hbonds", "site_jaccard", "n_site_res"])
        for i in range(nf):
            w.writerow([f"{i * dt_ps / 1000:.2f}", f"{prot_rmsd[i]:.4f}",
                        f"{lig_rmsd[i]:.4f}", f"{rg[i]:.4f}", contacts[i],
                        f"{mindist[i]:.4f}", len(hb_frames[i]),
                        f"{site_jac[i]:.3f}", site_nres[i]])

    def stat(k):
        a = np.array(res[k])
        return f"{a.mean():.3f}±{a.std():.3f}"
    print(f"  蛋白主链 RMSD {stat('prot_rmsd')} nm | 配体 RMSD {stat('lig_rmsd')} nm")
    print(f"  接触数 {stat('contacts')} | 氢键 {stat('hbonds')} | Rg {stat('rg')} nm")
    print(f"  最小原子间距 {np.mean(mindist)*10:.2f} Å (>4 Å 的帧占 "
          f"{100 * np.mean(np.array(mindist) > 0.4):.1f}%)")
    print(f"  位点 Jaccard {np.mean(site_jac):.3f}±{np.std(site_jac):.3f}"
          f" | 每帧接触残基 {np.mean(site_nres):.1f} 个 | 全轨迹累计 {len(contact_freq)} 个")
    print(f"  质心 MSD(1 ns) {msd[min(int(1.0 / (dt_ps / 1000)) - 1, len(msd) - 1)]:.3f} nm^2"
          f" | 表面扩散 D {d_surf:.2e} cm^2/s")
    print(f"  氢键条目(occ>5%) {len(hb_detail)}；首要位点 {list(contact_freq)[:6]}")
    return res


if __name__ == "__main__":
    for nm in (sys.argv[1:] or ["PFNA"]):
        print(f"=== {nm} ===")
        anal(nm)
