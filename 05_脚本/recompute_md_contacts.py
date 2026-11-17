#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
recompute_md_contacts.py — 从 MD 轨迹重算 PFNA 与三个基序/沟槽残基的接触
=====================================================================
输入：DATA14/md_runs/PFNA/traj.pdb（拓扑）+ md.xtc（1,001 帧，10 ns）
输出：07_结构/MD接触重算.csv、07_结构/MD逐残基接触.csv

口径（稿件必须与此一致）
  - 只算重原子（去氢）
  - 截断：<5 Å 记为接触（"contact"）；另报 <4 Å、<6 Å
  - "最小距离" = 配体全部重原子 × 目标残基全部重原子 的距离矩阵最小值
  - 帧：全部 1,001 帧，不抽稀

用法：python recompute_md_contacts.py [--stride N]
"""
import os, sys, csv
import numpy as np

try:
    import MDAnalysis as mda
except ImportError:
    sys.exit("需要 MDAnalysis：pip install MDAnalysis")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
TOP = os.path.join(ROOT, "DATA14", "md_runs", "PFNA", "traj.pdb")
XTC = os.path.join(ROOT, "DATA14", "md_runs", "PFNA", "md.xtc")
OUT = os.path.join(HERE, "..", "04_数据", "07_结构")

stride = 1
if "--stride" in sys.argv:
    stride = int(sys.argv[sys.argv.index("--stride") + 1])

u = mda.Universe(TOP, XTC)
lig = u.select_atoms("resname PFNA and not name H*")
print("配体重原子 %d；总帧数 %d；dt %.1f ps" % (lig.n_atoms, len(u.trajectory), u.trajectory.dt))

TARGETS = {
    "QDPVR (339-343)": (339, 343),
    "QVDRR M2 (351-355)": (351, 355),
    "ELLRQ M1 (251-255)": (251, 255),
}
sels = {k: u.select_atoms("protein and resid %d-%d and not name H*" % v) for k, v in TARGETS.items()}
for k, s in sels.items():
    print("  %-20s %d 个重原子  %s" % (k, s.n_atoms, [(r.resid, r.resname) for r in s.residues]))

rows, dmin = [], {k: [] for k in sels}
for ts in u.trajectory[::stride]:
    for k, s in sels.items():
        d = np.linalg.norm(lig.positions[:, None, :] - s.positions[None, :, :], axis=-1)
        dmin[k].append(d.min())

for k, d in dmin.items():
    d = np.array(d)
    rows.append({
        "target": k, "n_frames": len(d), "stride": stride,
        "median_min_dist_A": round(float(np.median(d)), 2),
        "p05_A": round(float(np.percentile(d, 5)), 2),
        "p95_A": round(float(np.percentile(d, 95)), 2),
        "frac_lt4A": round(float((d < 4).mean()), 3),
        "frac_lt5A": round(float((d < 5).mean()), 3),
        "frac_lt6A": round(float((d < 6).mean()), 3),
        "frac_lt8A": round(float((d < 8).mean()), 3),
    })

os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "MD接触重算.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for r in rows:
    print("%-20s 中位 %5.2f Å | <4 Å %5.1f%% | <6 Å %5.1f%% | <8 Å %5.1f%%"
          % (r["target"], r["median_min_dist_A"], 100 * r["frac_lt4A"], 100 * r["frac_lt6A"], 100 * r["frac_lt8A"]))

# 沟槽与基序逐残基接触频率（<5 Å）
RES = [44, 46, 77, 81, 82, 119, 251, 252, 253, 254, 255,
       339, 340, 341, 342, 343, 351, 352, 353, 354, 355]
per = []
for r in RES:
    at = u.select_atoms("protein and resid %d and not name H*" % r)
    if at.n_residues == 0:
        continue
    hits, ds = 0, []
    for ts in u.trajectory[::stride]:
        d = np.linalg.norm(lig.positions[:, None, :] - at.positions[None, :, :], axis=-1).min()
        ds.append(d); hits += d < 5.0
    per.append({"resid": r, "resname": at.residues[0].resname,
                "frac_contact_lt5A": round(hits / len(ds), 3),
                "median_min_dist_A": round(float(np.median(ds)), 2)})
with open(os.path.join(OUT, "MD逐残基接触.csv"), "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(per[0].keys())); w.writeheader(); w.writerows(per)
print("\n逐残基接触（<5 Å 频率）已写入 04_数据/07_结构/MD逐残基接触.csv")
for x in per:
    print("  %3d %-4s %5.1f%%  中位 %5.2f Å" % (x["resid"], x["resname"], 100 * x["frac_contact_lt5A"], x["median_min_dist_A"]))
