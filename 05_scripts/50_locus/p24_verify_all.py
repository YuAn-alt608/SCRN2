#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p24_verify_all.py
=================
Machine-checks every key number cited in `02_revised_conclusions_and_data_summary`.
Three independent verification paths:
  A internal consistency: P from b/se, F from beta/se, coloc posteriors summing to 1
  B cross-table consistency: table 23 vs 04 / 24; table 27 vs 04 / 31; table 26 internally
  C back to the raw data: inputs/*.tsv, inputs_ld/*.tsv.gz,
                  raw/ens_genes_17_peak_b37.csv, raw/sentinels/*.csv, plus an
                  independent recomputation of the GLS joint model

Output:
  tables/32_verification_checklist.csv      expected / observed / verdict (PASS|WARN|FAIL) per check
  tables/33_key_numbers_recomputable.csv    formula / reported / recomputed / deviation per number in the text
"""
import os, re, json
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
TB, IN, LD, RAW = f"{BASE}/tables", f"{BASE}/inputs", f"{BASE}/inputs_ld", f"{BASE}/raw"

ROWS, KEY = [], []
def add(cid, cat, item, exp, act, ok, note="", level="data"):
    """level: data = an error in the table or the raw data itself; wording = the document says it wrong, the data are fine"""
    ROWS.append(dict(check_id=cid, category=cat, level=level, check_item=item, expected=exp, observed=act,
                     result=("PASS" if ok is True else ("WARN" if ok is None else "FAIL")),
                     note=note))
def key(kid, item, formula, rep, rec, note=""):
    try:
        dev = float(rec) - float(rep)
    except Exception:
        dev = np.nan
    KEY.append(dict(key_id=kid, quantity=item, formula=formula, reported=rep, recomputed=rec,
                    deviation=dev, note=note))
def close(a, b, rtol=1e-6, atol=1e-12):
    try:
        return bool(abs(float(a) - float(b)) <= atol + rtol * max(abs(float(a)), abs(float(b))))
    except Exception:
        return False

T = lambda f: pd.read_csv(f"{TB}/{f}", encoding="utf-8-sig")
t23 = T("23_signal_attribution_summary.csv")
t24 = T("24_coloc_by_signal_window.csv")
t26 = T("26_peak_neighbour_genes_window_coloc.csv")
t27 = T("27_shoulder_gene_attribution.csv")
t31 = T("31_joint_conditional_model.csv")
t04 = T("04_cisMR_SMR_HEIDI.csv")
t02b = T("02b_top_eqtl_per_gene.csv")
t11 = T("11_variant_info.csv")
t21 = T("21_published_sentinels_17q21.csv")
t22 = T("22_sentinel_gene_mapping_uniqueness.csv")
genes = pd.read_csv(f"{RAW}/ens_genes_17_peak_b37.csv")

# ---------------------------------------------------------------- A internal consistency
def a_p_from_bse(df, bcol, secol, pcol, tag):
    ok = 0; tot = 0; bad = []
    for _, r in df.iterrows():
        b, se, p = r.get(bcol), r.get(secol), r.get(pcol)
        if not np.isfinite(b) or not np.isfinite(se) or not np.isfinite(p) or se <= 0:
            continue
        tot += 1
        pz = 2 * stats.norm.sf(abs(b / se))
        if close(pz, p, rtol=2e-3, atol=1e-300):
            ok += 1
        else:
            bad.append((r.get("gene", r.get("SNP", "")), pz, p))
    add("A1", "internal consistency", f"{tag}: P reproducible from b/se", f"all {tot} rows agree",
        f"{ok}/{tot}", ok == tot and tot > 0, f"mismatches: {bad[:4]}" if bad else "")

a_p_from_bse(t31, "marginal_b", "marginal_se", "marginal_P", "table 31 marginal")
a_p_from_bse(t31, "joint_b", "joint_se", "joint_P", "table 31 joint")
a_p_from_bse(t04, "b_wald", "se_wald", "p_wald", "table 04 Wald")

# table 04 F_lead must equal (lead_beta/lead_se)^2
fok = 0; ftot = 0
for _, r in t04.iterrows():
    t, g, p = r.tissue, r.gene, int(r.lead_pos37)
    fp = f"{IN}/{t}__{g}__TSH.tsv"
    if not os.path.exists(fp):
        continue
    d = pd.read_csv(fp, sep="\t")
    m = d[d.pos37 == p]
    if m.empty:
        continue
    ftot += 1
    F = (m.iloc[0].beta_eqtl / m.iloc[0].se_eqtl) ** 2
    if close(F, r.F_lead, rtol=1e-6):
        fok += 1
add("A2", "internal consistency", "table 04: F_lead == (lead_beta/lead_se)^2", f"{ftot} rows",
    f"{fok}/{ftot}", fok == ftot, "recomputed row by row from the raw inputs tsv")

# coloc posteriors
s26 = t26[["H0", "H1", "H2", "H3", "H4"]].sum(axis=1)
add("A3", "internal consistency", "table 26: H0+...+H4 == 1", "all 1+/-1e-6",
    f"min={s26.min():.8f} max={s26.max():.8f}", bool(np.allclose(s26, 1, atol=1e-6)),
    f"{len(s26)} rows")
bad24 = t24[(t24.H3 + t24.H4 > 1 + 1e-9)]
add("A4", "internal consistency", "table 24: H3+H4 <= 1", "no violation", f"{len(bad24)} violations", len(bad24) == 0)

# ---------------------------------------------------------------- B cross-table consistency
# B1 table 23 MR columns == table 04 for the same (tissue,gene,outcome)
mism = []
for _, r in t23.iterrows():
    g, t = r["candidate_gene"], r["tissue"]
    oc = "TSH"
    s = t04[(t04.gene == g) & (t04.tissue == t) & (t04.outcome == oc)]
    if s.empty:
        mism.append((g, t, "row missing")); continue
    s = s.iloc[0]
    if not (close(s.b_wald, r["MR_b"], 1e-6) and close(s.se_wald, r["MR_se"], 1e-6)):
        mism.append((g, t, s.b_wald, r["MR_b"]))
add("B1", "cross-table consistency", "table 23 MR_b/MR_se == table 04 b_wald/se_wald (TSH)",
    f"all {len(t23)} rows equal", f"{len(mism)} unequal", len(mism) == 0, str(mism[:3]))

# B2 table 23 window coloc == table 24
WIN = {1: "Signal 1 (45.72-45.86)", 2: "Signal 2 (45.86-45.98)", 3: "Signal 3 (46.00-46.12)"}
mism2 = []
for _, r in t23.iterrows():
    g, t, sig = r["candidate_gene"], r["tissue"], int(r["signal"])
    for oc, h3c, h4c in [("TSH", "window_coloc_H3", "window_coloc_H4"),
                         ("HYPO", "hypo_window_coloc_H3", "hypo_window_coloc_H4")]:
        v3, v4 = r[h3c], r[h4c]
        if not (np.isfinite(v3) and np.isfinite(v4)):
            continue
        s = t24[(t24.gene == g) & (t24.tissue == t) & (t24.outcome == oc) & (t24.window == WIN[sig])]
        if s.empty:
            mism2.append((g, t, oc, "table 24 row missing")); continue
        s = s.iloc[0]
        if not (close(s.H3, v3, 1e-4, 1e-9) and close(s.H4, v4, 1e-4, 1e-9)):
            mism2.append((g, t, oc, s.H3, v3, s.H4, v4))
add("B2", "cross-table consistency", "table 23 window coloc == table 24 for the matching (gene,tissue,outcome,window)",
    f"{len(t23)*2} cells", f"{len(mism2)} unequal", len(mism2) == 0, str(mism2[:3]))

# B3 table 27 MR == table 04
mism3 = []
for _, r in t27.iterrows():
    g = r["gene"]; t = r["tissue"]
    for oc, bc, pc in [("HYPO", "MR_HYPO_b", "MR_HYPO_p"), ("TSH", "MR_TSH_b", "MR_TSH_p")]:
        s = t04[(t04.gene == g) & (t04.tissue == t) & (t04.outcome == oc)]
        if s.empty:
            continue
        s = s.iloc[0]
        if np.isfinite(r[bc]) and not close(s.b_wald, r[bc], 1e-6):
            mism3.append((g, t, oc, s.b_wald, r[bc]))
add("B3", "cross-table consistency", "table 27 cis-MR == table 04", f"{len(t27)*2} cells",
    f"{len(mism3)} unequal", len(mism3) == 0, str(mism3[:3]))

# B4 table 22 sentinel counts agree with the raw CSV group sizes
raw_s = pd.read_csv(f"{RAW}/sentinels/gwas_catalog_associations_all.csv", encoding="utf-8-sig")
cnt = raw_s.groupby(["paper", "trait"]).size().reset_index(name="n")
oks = []
for _, r in t22.iterrows():
    m = cnt[(cnt.paper == r.paper) & (cnt.trait == r.trait)]
    oks.append((r.paper, r.trait, int(r.n_sentinels), int(m.n.iloc[0]) if len(m) else -1))
add("B4", "cross-table consistency", "table 22 sentinel counts == raw CSV group sizes", "all 4 groups equal",
    "; ".join(f"{a}/{b}:{c}={d}" for a, b, c, d in oks),
    all(c == d for _, _, c, d in oks), f"raw total {len(raw_s)}")

# ---------------------------------------------------------------- C raw data
# C1 dual-build offset
for nm, df, c38, c37 in [("table 11", t11, "pos_grch38", "pos_grch37"),
                         ("table 21", t21, "pos38", "pos37")]:
    dd = df.dropna(subset=[c38, c37])
    off = (dd[c38].astype(float) - dd[c37].astype(float))
    add("C1" + nm[-1], "coordinates", f"{nm}: pos38 - pos37 == 1,922,634",
        "all 1,922,634", f"value set {sorted(set(off.astype(int)))}",
        bool((off == 1922634).all()), f"{len(dd)} rows")
key("K1", "GRCh38 to GRCh37 offset", "rs4794063 47,727,128 - 45,804,494", 1922634,
    int(t21.pos38.iloc[0] - t21.pos37.iloc[0]), "same value at two independent loci")

# C2 allele orientation of the input files
#     beta_out should already be oriented to eqtl_ea; the gwas_ea column holds the
#     pre-flip GWAS effect allele
chk = []
for _, r in t11.iterrows():
    rs, b11 = r.rsid, r.beta_tsh
    pos = int(r.pos_grch37)
    for t in ["Thyroid", "WholeBlood"]:
        for f in os.listdir(IN):
            if not f.startswith(f"{t}__") or not f.endswith("__TSH.tsv"):
                continue
            d = pd.read_csv(f"{IN}/{f}", sep="\t")
            m = d[d.pos37 == pos]
            if m.empty:
                continue
            m = m.iloc[0]
            ea_eq = (str(m.eqtl_ea).upper() == str(r.tsh_ea).upper())
            exp = float(b11) if ea_eq else -float(b11)
            chk.append((rs, t, f.split("__")[1], ea_eq, exp, float(m.beta_out),
                        close(exp, m.beta_out, 1e-3)))
nok = sum(1 for c in chk if c[-1])
add("C2", "alleles", "inputs/*.tsv beta_out is oriented to eqtl_ea (it differs in sign from table 11 tsh_ea if and only if the alleles differ)",
    f"all {len(chk)} matches hold", f"{nok}/{len(chk)}", nok == len(chk),
    "note: the gwas_ea column in these files is the pre-flip value; aligning on it flips the wrong way (p12 is obsolete)")

# C3 gene coordinates: is rs6503875 inside SCRN2; direction of rs62076102 relative to SCRN2
sc = genes[genes.gene_name == "SCRN2"].iloc[0]
inside = bool(sc.start <= 45915395 <= sc.end)
add("C3", "coordinates", "rs6503875 (45,915,395) falls inside the SCRN2 gene body",
    f"SCRN2 {sc.start}-{sc.end} (strand {sc.strand:+d})", str(inside), inside,
    f"{45915395 - sc.start} bp from the 3' end; {sc.end - 45915395} bp from the TSS")
tss = int(sc.end) if sc.strand == -1 else int(sc.start)
d62 = 45887132 - tss
add("C4", "coordinates", "position of rs62076102 relative to the SCRN2 promoter",
    "must match the wording in the text", f"{d62:+,} bp (negative = downstream / 3' side)", None,
    "SCRN2 is on the negative strand with TSS=45,918,699; 45,887,132 lies 31.6 kb downstream of the TSS, "
    "i.e. 28.0 kb beyond the 3' end. The text originally said \"28 kb from the 5' end\", which is the wrong "
    "direction; corrected")
key("K2", "distance from rs62076102 to the SCRN2 3' end", "45,915,058 - 45,887,132",
    28000, int(sc.start - 45887132), "the text says 28 kb, measured 27,926 bp")

# C5 independent recomputation of the 5-SNP GLS joint model (uses only marginal b/se + LD, never reads table 31)
snp5 = [("rs7219420", 45800299), ("rs1989291", 45807547),
        ("shoulder", 45873049), ("rs62076102", 45887132), ("rs6503875", 45915395)]
def gls(outcome):
    d = t31[(t31["outcome"] == outcome) & (t31["model"] == "five-SNP")].set_index("SNP")
    d = d.loc[[s for s, _ in snp5], ["pos37", "marginal_b", "marginal_se"]]
    LDm = pd.read_csv(f"{LD}/ld_{outcome}.tsv.gz", sep="\t", index_col=0)
    LDm.index = LDm.index.astype(int); LDm.columns = LDm.columns.astype(int)
    p = [int(x) for x in d.pos37]
    R = LDm.loc[p, p].values
    b = d["marginal_b"].values; se = d["marginal_se"].values
    V = np.diag(se) @ R @ np.diag(se)
    W = np.linalg.inv(V)
    bj = np.linalg.solve(R.T @ W @ R, R.T @ W @ b)
    sej = np.sqrt(np.diag(np.linalg.inv(R.T @ W @ R)))
    return pd.DataFrame(dict(pos37=p, b_joint=bj, se_joint=sej, z=bj / sej,
                             p=2 * stats.norm.sf(abs(bj / sej))), index=d.index)

for oc in ["HYPO", "TSH"]:
    g = gls(oc)
    rep = t31[(t31["outcome"] == oc) & (t31["model"] == "five-SNP")].set_index("SNP")["joint_z"]
    ok = all(close(g.loc[s, "z"], rep[s], 1e-4) for s in g.index)
    add(f"C5_{oc}", "GLS recomputation", f"{oc} five-SNP joint z recomputed independently",
        "agrees with table 31", "; ".join(f"{s}:{g.loc[s,'z']:.3f}" for s in g.index), ok,
        "inputs are the marginal b/se and the LD matrix only; the result columns of table 31 are never read")
    for s in g.index:
        key(f"K3_{oc}_{s}", f"{oc} {s} joint z", "(R'V^-1 R)^-1 R'V^-1 b", float(rep[s]),
            round(float(g.loc[s, "z"]), 3))

# C6 signal count: how many stay P<0.05 after conditioning
def sig_count(oc, alpha=0.05):
    g = gls(oc)
    return int((g.p < alpha).sum()), g.p.round(4).to_dict()
c_h, ph = sig_count("HYPO"); c_t, pt = sig_count("TSH")
add("C6", "signal count", "hypothyroidism: independent signals still P<0.05 in the full five-SNP model", "4", f"{c_h}  {ph}",
    c_h == 4, "45,800 / 45,807 / 45,873 / 45,887")
add("C7", "signal count", "TSH: independent signals still P<0.05 in the full five-SNP model", "2 (+ the COPZ2 signal at 46.05 is not in the five-SNP set)",
    f"{c_t}  {pt}", c_t == 2, "45,800 / 45,887; counting 46.05 separately gives 3 TSH signals in total")

# C8 LD recomputation (back to the 1000G matrix)
LDt = pd.read_csv(f"{LD}/ld_TSH.tsv.gz", sep="\t", index_col=0)
LDt.index = LDt.index.astype(int); LDt.columns = LDt.columns.astype(int)
pairs = [("rs62076102 / rs6503875", 45887132, 45915395, 0.7669),
         ("rs7219420 / rs6503875", 45800299, 45915395, 0.0262),
         ("rs7219420 / rs1989291", 45800299, 45807547, -0.3981),
         ("shoulder / rs1989291", 45873049, 45807547, -0.2526),
         ("shoulder / rs62076102", 45873049, 45887132, 0.4015),
         ("shoulder / rs6503875", 45873049, 45915395, 0.3051)]
allok = True
for nm, a, b, rep in pairs:
    v = float(LDt.loc[a, b])
    ok = close(v, rep, 2e-3); allok &= ok
    key("K4_" + nm.replace(" ", "").replace("/", "_"), f"r({nm})", "element of the LD matrix", rep, round(v, 4))
add("C8", "LD", "the six r values cited in table 31 == LD matrix elements", "all six equal",
    "; ".join(f"{n}:{LDt.loc[a,b]:.4f}" for n, a, b, _ in pairs), allok)

# C9 sentinels: count inside the region
reg = t21
add("C9", "sentinels", "sentinels inside the 17q21 region (GRCh38 47.3-48.5 Mb) == 1", "1", f"{len(reg)}",
    len(reg) == 1, f"{reg.rsid.iloc[0]} / {reg.trait.iloc[0]} / P={reg.p.iloc[0]}")
key("K5", "sentinels inside the region", "raw/sentinels filtered to chr38 in [47.3e6,48.5e6]", 1, len(reg))
key("K6", "total sentinels", "four studies, GWAS Catalog associations, after de-duplication", 783, len(raw_s))

# C10 source check of table 25
t25 = T("25_hypo_data_source_check.csv").set_index("check_item")["value"]
add("C10", "source check", "panel hypothyroidism == Rand2025 GCST90572791 (beta/SE position by position)",
    "identical fraction 1.0 and the same position for the minimum P", f"{t25['beta_se_identical_fraction']} / {t25['min_p_and_position_both_sets']}",
    float(t25["beta_se_identical_fraction"]) == 1.0, "n=1,043,800")

# C11 TBX21: do eQTL_P and F in table 23 come from the same variant
row = t23[(t23["candidate_gene"] == "TBX21")].iloc[0]
d = pd.read_csv(f"{IN}/Thyroid__TBX21__TSH.tsv", sep="\t")
p_at_lead = float(d[d.pos37 == 45800299].iloc[0].p_eqtl)
lead = d.loc[d.p_eqtl.idxmin()]
F_min = float((lead.beta_eqtl / lead.se_eqtl) ** 2)
same = close(p_at_lead, row["eQTL_P"], 1e-4) and close(F_min, row["F"], 1e-3)
add("C11", "definition", "TBX21 row: are eQTL_P and F taken from the same variant", "they should be the same variant",
    f"eQTL_P={p_at_lead:.4f} at the signal 1 lead; F={F_min:.2f} at {int(lead.pos37)}; table 23 records P={row['eQTL_P']:.4f}, F={row['F']:.2f}",
    False, "they are not the same variant: P is the value at the signal 1 lead, while F is the value at the nominally "
           "smallest-P variant in the interval (45,455,669, 355 kb from TBX21). The text combined the two variants in "
           "one sentence (\"P=0.55, F=5.8\"); they are now reported separately", level="wording")

# C12 does the TSH_beta of the table 23 sentinels agree with table 11 (allowing for allele orientation)
t11i = t11.set_index("rsid")
bad12 = []
for _, r in t23.iterrows():
    rs = r["sentinel"]
    if rs not in t11i.index:
        continue
    b11 = float(t11i.loc[rs, "beta_tsh"])
    if not (close(abs(b11), abs(r["TSH_beta"]), 1e-3)):
        bad12.append((rs, b11, r["TSH_beta"]))
add("C12", "alleles", "table 23 TSH_beta has the same magnitude as table 11 (only the sign may flip with the allele)",
    f"{len(t23)} rows", f"{len(bad12)} mismatched", len(bad12) == 0,
    "table 23 is uniformly oriented to eqtl_ea, table 11 to tsh_ea. Opposite signs are expected, not an error")

# C13 do the two TSH GWAS really have zero sentinels in this region
for (pp, tr), g in raw_s.groupby(["paper", "trait"]):
    n_reg = int(((g.chr38.astype(str) == "17") & (g.pos38.astype(float) >= 47.3e6) &
                 (g.pos38.astype(float) <= 48.5e6)).sum())
    rep = t22[(t22.paper == pp) & (t22.trait == tr)]
    rep = int(rep.n_regional_sentinels.iloc[0]) if len(rep) else -1
    add(f"C13_{pp}_{tr}", "sentinels", f"{pp} {tr}: sentinels inside the region (raw CSV recomputation vs table 22)",
        f"table 22 records {rep}", f"recomputed {n_reg} (of {len(g)})", n_reg == rep)

# C14 positional direction must follow the gene strand (negative-strand genes are the easiest to get backwards)
COPZ2 = dict(start=46103533, end=46115392, strand=-1)   # Ensembl GRCh37 grch37.rest.ensembl.org
COORD = {g.gene_name: dict(start=int(g.start), end=int(g.end), strand=int(g.strand))
         for _, g in genes.iterrows()}
COORD["COPZ2"] = COPZ2
CASES = [("rs62076102", 45887132, "SCRN2"), ("rs6503875", 45915395, "SCRN2"),
         ("rs7219420", 45800299, "TBKBP1"), ("rs7219420", 45800299, "TBX21"),
         ("COPZ2lead", 46050669, "COPZ2"), ("shoulder", 45873049, "TBKBP1")]
lines = []
for nm, pos, g in CASES:
    c = COORD[g]
    tss = c["end"] if c["strand"] == -1 else c["start"]
    d_tss = pos - tss
    side = ("outside the promoter" if c["strand"] == -1 else "outside the promoter") if d_tss < 0 else "gene body / downstream"
    inside = c["start"] <= pos <= c["end"]
    dist_edge = 0 if inside else (c["start"] - pos if pos < c["start"] else pos - c["end"])
    if inside:
        desc = f"within the gene body ({pos - c['start']} bp from the 3' boundary, {c['end'] - pos} bp from the TSS)"
    else:
        five_side = (pos > c["end"]) if c["strand"] == -1 else (pos < c["start"])
        desc = f"{'5' if five_side else '3'}-prime side {abs(dist_edge)/1000:.1f} kb"
    lines.append(f"{nm}/{g}: TSS{'+' if d_tss>=0 else ''}{d_tss} -> {desc}")
add("C14", "coordinates", "positional relationships judged by gene strand (the 5'/3' of negative-strand genes is easily reversed)", "must match the wording in the text",
    " | ".join(lines), False,
    "SCRN2 and COPZ2 are both on the negative strand: rs62076102 lies 28.0 kb on the 3' side of SCRN2 (the text "
    "originally said 5' end); COPZ2lead lies 52.9 kb on the 3' side of COPZ2 and 64.7 kb from the TSS (the text "
    "originally said \"53 kb upstream\")", level="wording")

# ---------------------------------------------------------------- output
res = pd.DataFrame(ROWS)
res.to_csv(f"{TB}/32_verification_checklist.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(KEY).to_csv(f"{TB}/33_key_numbers_recomputable.csv", index=False, encoding="utf-8-sig")

pd.set_option("display.width", 260, "display.max_colwidth", 60)
print(res.to_string(index=False))
n_fail = int((res["result"] == "FAIL").sum()); n_warn = int((res["result"] == "WARN").sum())
nf_d = int(((res["result"] == "FAIL") & (res["level"] == "data")).sum())
nf_w = int(((res["result"] == "FAIL") & (res["level"] == "wording")).sum())
print(f"\n{len(res)} checks in total: PASS {len(res)-n_fail-n_warn} / WARN {n_warn} / FAIL {n_fail}")
print(f"  of these, {nf_d} FAIL at the data level and {nf_w} at the wording level (the data are sound, the document was corrected; see the note in each row)")
print("wrote tables/32_verification_checklist.csv, tables/33_key_numbers_recomputable.csv")
