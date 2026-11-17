"""06 · 平衡与生产 MD：NVT 100 ps -> NPT 200 ps -> 生产 MD

用法：
    python scripts/06_run_md.py T3            # 全套
    python scripts/06_run_md.py T3 nvt        # 只跑某一段
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MDP = ROOT / "mdp"
GMX = "gmx"

STAGES = [
    ("nvt", "em.gro", "em.gro"),
    ("npt", "nvt.gro", "nvt.gro"),
    ("md", "npt.gro", "npt.gro"),
]


def sh(cmd, cwd, check=True):
    r = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    if check and r.returncode != 0:
        print(r.stdout[-3000:]); print(r.stderr[-2000:])
        raise SystemExit(f"FAILED: {' '.join(map(str, cmd))}")
    return r


def run(name, stage):
    work = ROOT / "md_runs" / name
    _, cin, ref = [s for s in STAGES if s[0] == stage][0]
    cmd = [GMX, "grompp", "-f", str(MDP / f"{stage}.mdp"), "-c", cin,
           "-r", ref, "-p", "topol.top", "-n", "index.ndx",
           "-o", f"{stage}.tpr", "-maxwarn", "2"]
    sh(cmd, work)
    mcmd = [GMX, "mdrun", "-deffnm", stage, "-ntmpi", "1"]
    if (work / f"{stage}.cpt").exists():        # 断点续跑
        mcmd += ["-cpi", f"{stage}.cpt", "-append"]
    r = sh(mcmd, work, check=False)
    if r.returncode != 0:
        print(r.stdout[-2500:]); print(r.stderr[-1500:])
        raise SystemExit(f"mdrun FAILED: {stage}")
    log = (work / f"{stage}.log").read_text(errors="ignore").splitlines()
    perf = [ln for ln in log if "Performance:" in ln or "ns/day" in ln]
    print(f"  {name} {stage}: " + (perf[-1].strip() if perf else "见日志"))


if __name__ == "__main__":
    nm = sys.argv[1]
    want = sys.argv[2:] or [s[0] for s in STAGES]
    for st in want:
        print(f"=== {nm} {st} ===")
        run(nm, st)
