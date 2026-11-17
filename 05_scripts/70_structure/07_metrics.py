"""07 - the six trajectory metrics (ligand RMSD / contacts / hydrogen bonds / RMSF / protein RMSD / Rg)

Definitions (quoted directly in the manuscript):
  1. ligand RMSD     : after superposition on the protein backbone, RMSD of ligand heavy atoms
                       to frame 1 (ligand drift alone)
  2. site contacts   : number of ligand-heavy to protein-heavy atom pairs closer than 0.6 nm
  3. hydrogen bonds  : donor-acceptor distance <= 0.35 nm and the D-H...A angle >= 150 degrees
                       (that is, H-D-A <= 30 degrees); lifetime = intermittent autocorrelation
                       C(tau)=<h(0)h(tau)>/<h(0)> integrated to 100 ps
  4. ligand RMSF     : per-atom fluctuation under the same superposition
  5. protein RMSD    : backbone (N,CA,C,O) RMSD
  6. radius of gyration Rg : mass-weighted over the protein
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
HB_ANGLE = 30.0        # deg, upper bound on the H-D-A angle
LIFE_LAG = 100.0       # ps, upper limit of the autocorrelation integral


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
    """Read residue number and name for every atom from complex.gro (the topology has no residue numbering)."""
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
    """Hydrogen bonds in one frame; returns the set of (donor, acceptor) pairs present.

    Criteria match gmx hbond: D-A <= 0.35 nm and the D-H...A angle at H >= 180-30 = 150 degrees
    (a linear hydrogen bond is 180 degrees). A donor may carry several hydrogens; the best is used.
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
    print(f"  {name}: {nf} frames x {frames.shape[1]} atoms")

    elem, aname = top["elem"], top["aname"]
    heavy_l = lidx[elem[lidx] != "H"]
    heavy_p = pidx[elem[pidx] != "H"]
    bb = pidx[np.isin(aname[pidx], ["N", "CA", "C", "O"])]
    ca = pidx[aname[pidx] == "CA"]

    # donors / acceptors
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
        # per-frame set of contacting residues (for the site-identity drift, Jaccard)
        hit = set(res_of_heavy[(d < CUT_CONTACT).any(0)].tolist())
        site_sets.append(hit)
        for r in hit:
            contact_freq[r] = contact_freq.get(r, 0) + 1
        fitted_p[s] = y[pidx]
        lig_fit[s] = y[lidx]
        hb_frames.append(hbond_sets(y, pairs))

    contact_freq = {k: v / nf for k, v in sorted(contact_freq.items(), key=lambda kv: -kv[1])}

    # site identity: Jaccard similarity with the frame-1 contact-residue set
    s0 = site_sets[0]
    site_jac = [len(s & s0) / max(len(s | s0), 1) for s in site_sets]
    site_nres = [len(s) for s in site_sets]

    # hydrogen-bond lifetime
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

    # lateral surface diffusion: MSD of the ligand heavy-atom centroid, in superposed coordinates,
    # so bulk translation and rotation are removed
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

    # per-frame time series (so a reader can check without parsing the JSON)
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
    print(f"  protein backbone RMSD {stat('prot_rmsd')} nm | ligand RMSD {stat('lig_rmsd')} nm")
    print(f"  contacts {stat('contacts')} | hydrogen bonds {stat('hbonds')} | Rg {stat('rg')} nm")
    print(f"  minimum atom distance {np.mean(mindist)*10:.2f} A (frames above 4 A: "
          f"{100 * np.mean(np.array(mindist) > 0.4):.1f}%)")
    print(f"  site Jaccard {np.mean(site_jac):.3f}+/-{np.std(site_jac):.3f}"
          f" | contacting residues per frame {np.mean(site_nres):.1f} | cumulative {len(contact_freq)}")
    print(f"  centroid MSD(1 ns) {msd[min(int(1.0 / (dt_ps / 1000)) - 1, len(msd) - 1)]:.3f} nm^2"
          f" | surface diffusion D {d_surf:.2e} cm^2/s")
    print(f"  hydrogen-bond entries (occ>5%) {len(hb_detail)}; leading residues {list(contact_freq)[:6]}")
    return res


if __name__ == "__main__":
    for nm in (sys.argv[1:] or ["PFNA"]):
        print(f"=== {nm} ===")
        anal(nm)
