import csv, json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
from scipy.ndimage import distance_transform_edt

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

EVID = "D:/NewBeginning/MR1/DATA18_manuscript/SCRN2_TSH_story_evidence"
def load(p): return list(csv.DictReader(open(f"{EVID}/{p}", encoding="utf-8-sig")))

mr = load("01_discovery/MR_significant_FDR05.csv")
gw = load("02_primary_replication/gw_thyroid_mr.csv")
co = load("03_locus_gene_assignment/coloc_results.csv")
a2 = load("04_expression_context/a2_summary.csv")
ph = load("03_locus_gene_assignment/phewas_eqtl_targets.csv")

TOOLS = {"rs6503875", "rs62076102", "rs4794224"}
S = [
    {r["gene"] for r in mr if r["outcome"] == "TSH"},
    {r["gene"] for r in gw if r["outcome"] == "TSH" and float(r["fdr"]) < 0.05},
    {r["gene"] for r in co if r["outcome"] == "TSH" and float(r["PP.H4"]) >= 0.8},
    {r["gene"] for r in a2 if r["pct_same_sign"] not in ("", "nan") and float(r["pct_same_sign"]) >= 99},
    {r["gene"] for r in ph if r["snp"] in TOOLS and int(r["rank"]) <= 3},
]

members = {}
for g in set().union(*S):
    m = sum(1 << i for i in range(5) if g in S[i])
    members.setdefault(m, []).append(g)
cnt = {m: len(v) for m, v in members.items()}

A, B, D, ROT, DELTA = 1.3623110744260034, 0.5780690630376263, 0.49319202093268455, -59.336879058595386, -11.69853818521527
N = 1100
xs = np.linspace(-2.2, 2.2, N); ys = np.linspace(-2.2, 2.2, N)
X, Y = np.meshgrid(xs, ys)
M = np.zeros(X.shape, dtype=np.uint8)
centers = []
for k in range(5):
    th = np.deg2rad(90 + 72 * k + ROT)
    cx, cy = D * np.cos(th), D * np.sin(th)
    centers.append((cx, cy, th))
    ang = th + np.deg2rad(DELTA)
    dx, dy = X - cx, Y - cy
    u = dx * np.cos(ang) + dy * np.sin(ang)
    v = -dx * np.sin(ang) + dy * np.cos(ang)
    M |= (((u / A) ** 2 + (v / B) ** 2) <= 1).astype(np.uint8) << k

# ---------- label placement: greedy max-min distance inside each region ----------
EDT = distance_transform_edt(M > 0)
placed = []          # (x, y, r)
cands = {}
for m in range(1, 32):
    if cnt.get(m, 0) == 0:
        continue
    sub = (M == m)
    d = distance_transform_edt(sub)
    if d.max() <= 0:
        continue
    idx = np.argwhere(d >= 0.80 * d.max())
    pts = np.column_stack((xs[idx[:, 1]], ys[idx[:, 0]]))
    if len(pts) > 900:
        pts = pts[np.random.default_rng(0).choice(len(pts), 900, replace=False)]
    cands[m] = pts

locs = {}
order = sorted(cands, key=lambda m: (m != 31, len(cands[m])))
for m in order:
    pts = cands[m]
    if not placed:
        best = pts[np.argmax([distance_transform_edt(M == m)[int(np.argmin(np.abs(ys - p[1]))), int(np.argmin(np.abs(xs - p[0])))] for p in pts[:0]] or [0])] if False else pts[len(pts)//2]
        best = pts.mean(axis=0)
    dm = distance_transform_edt(M == m)
    def depth(p):
        i = int(np.argmin(np.abs(xs - p[0]))); j = int(np.argmin(np.abs(ys - p[1])))
        return dm[j, i] / max(dm.max(), 1e-9)
    bestd, best = -9, None
    for p in pts:
        sep = min((np.hypot(p[0]-q[0], p[1]-q[1]) for q in placed), default=9.9)
        sc = sep + 0.55 * depth(p)
        if sc > bestd:
            bestd, best = sc, p
    locs[m] = tuple(best)
    placed.append((best[0], best[1], 0.13))

# ---------- draw ----------
COLS = ["#7FB77E", "#E76F8E", "#F2C14E", "#6E9BD1", "#A98BC0"]

CFG = {
 "en": dict(
   font="Arial", out="Fig_SCRN2_five_layer_venn",
   title="Five evidence layers converge on a single gene: SCRN2",
   labels=["Whole-blood cis-eQTL MR\neQTLGen; TSH FDR<0.05  (n=337)",
           "Thyroid cis-eQTL MR\nGTEx v8; TSH FDR<0.05  (n=180)",
           "Colocalization\ncoloc.abf PP.H4\u22650.80  (n=26)",
           "Cross-tissue coherence\n\u226599% of 50 GTEx tissues  (n=3)",
           "Target-gene assignment\neQTL top-3 for instrument SNPs  (n=4)"],
   center="SCRN2", center_sub="1 gene",
   foot="Numbers are genes shared by that exact combination of layers; blank regions contain none. "
        "See caption for layer definitions (all layers start from the same TSH GWAS).",
 ),
 "zh": dict(
   font="Microsoft YaHei", out="Fig_SCRN2_five_layer_venn_zh",
   title="\u4e94\u5c42\u8bc1\u636e\u6536\u655b\u5230\u4e00\u4e2a\u57fa\u56e0\uff1aSCRN2",
   labels=["\u5168\u8840 cis-eQTL MR\neQTLGen\uff1bTSH FDR<0.05\uff08n=337\uff09",
           "\u7532\u72b6\u817a cis-eQTL MR\nGTEx v8\uff1bTSH FDR<0.05\uff08n=180\uff09",
           "\u5171\u5b9a\u4f4d\ncoloc.abf PP.H4\u22650.80\uff08n=26\uff09",
           "\u8de8\u7ec4\u7ec7\u4e00\u81f4\n50 \u4e2a GTEx \u7ec4\u7ec7\u4e2d\u226599%\u540c\u5411\uff08n=3\uff09",
           "\u9776\u57fa\u56e0\u5f52\u5c5e\n\u5de5\u5177 SNP \u7684 eQTL \u524d\u4e09\uff08n=4\uff09"],
   center="SCRN2", center_sub="\u552f\u4e00 1 \u4e2a\u57fa\u56e0",
   foot="\u6570\u5b57\u4e3a\u540c\u65f6\u6ee1\u8db3\u8be5\u7ec4\u5408\u7684\u57fa\u56e0\u6570\uff1b\u7a7a\u767d\u533a\u57df\u4e3a 0\u3002\u4e94\u5c42\u8bc1\u636e\u5747\u57fa\u4e8e\u540c\u4e00\u5957 TSH GWAS \u4e0e eQTL \u8d44\u6e90\uff0c\u4e92\u4e3a\u8865\u5145\u800c\u975e\u7edf\u8ba1\u72ec\u7acb\u3002",
 ),
}

def draw(cfg):
    plt.rcParams["font.family"] = cfg["font"]
    fig, ax = plt.subplots(figsize=(9.4, 8.2))
    for k in range(5):
        cx, cy, th = centers[k]
        ax.add_patch(Ellipse((cx, cy), 2*A, 2*B, angle=np.rad2deg(th + np.deg2rad(DELTA)),
                             facecolor=COLS[k], alpha=0.26, edgecolor=COLS[k], lw=2.0, zorder=2))
    for m in range(1, 32):
        if cnt.get(m, 0) == 0:
            continue
        x, y = locs[m]
        if m == 31:
            ax.text(x, y + 0.05, cfg["center"], ha="center", va="center", fontsize=15.5, fontweight="bold",
                    color="#1a1a1a", zorder=6,
                    path_effects=[pe.withStroke(linewidth=4, foreground="white", alpha=0.95)])
            ax.text(x, y - 0.17, cfg["center_sub"], ha="center", va="center", fontsize=9.5, color="#333333",
                    zorder=6, path_effects=[pe.withStroke(linewidth=3, foreground="white", alpha=0.95)])
        else:
            ax.text(x, y, str(cnt[m]), ha="center", va="center", fontsize=11.5, color="#2b2b2b", zorder=6,
                    path_effects=[pe.withStroke(linewidth=3.4, foreground="white", alpha=0.92)])
    for k in range(5):
        cx, cy, th = centers[k]
        ang = th + np.deg2rad(DELTA)
        tip = (cx + (A + 0.03)*np.cos(ang), cy + (A + 0.03)*np.sin(ang))
        lab = (cx + (A + 0.42)*np.cos(ang), cy + (A + 0.42)*np.sin(ang))
        ha = "left" if np.cos(ang) > 0.3 else ("right" if np.cos(ang) < -0.3 else "center")
        va = "bottom" if np.sin(ang) > 0.4 else ("top" if np.sin(ang) < -0.4 else "center")
        ax.annotate(cfg["labels"][k], xy=tip, xytext=lab, ha=ha, va=va, fontsize=9.6, color="#1a1a1a",
                    linespacing=1.5, zorder=7,
                    arrowprops=dict(arrowstyle="-", color=COLS[k], lw=1.1, shrinkA=0, shrinkB=3))
    ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.05, 2.55)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(cfg["title"], fontsize=13.5, fontweight="bold", pad=30)
    fig.text(0.5, 0.03, cfg["foot"], ha="center", fontsize=8.4, color="#666666")
    fig.tight_layout(rect=[0, 0.035, 1, 1])
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{cfg['out']}.{ext}", dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)

for cfg in CFG.values():
    draw(cfg)

json.dump({str(m): {"n": cnt.get(m, 0), "genes": sorted(members.get(m, []))} for m in range(1, 32)},
          open("region_counts.json", "w"), ensure_ascii=False, indent=1)
print("ok", {m: cnt[m] for m in sorted(cnt)})
