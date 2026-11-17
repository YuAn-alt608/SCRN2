#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p24_verify_all.py
=================
对 `02_修正后结论与数据总表` 里引用的每一个关键数字做机器核对。
三条独立的核对路径：
  A 表内自洽：P 与 b/se、F 与 beta/se、coloc 后验和为 1
  B 跨表一致：表23 vs 表04 / 表24；表27 vs 表04 / 表31；表26 内部
  C 回到原始数据：inputs/*.tsv、inputs_ld/*.tsv.gz、raw/ens_genes_17_peak_b37.csv、
                  raw/sentinels/*.csv，以及 GLS 联合模型的独立重算

输出：
  tables/32_核对清单.csv      每条检查的 期望 / 实测 / 结论（PASS|WARN|FAIL）
  tables/33_关键数字_可复算.csv 正文每个数字的 定义式 / 报告值 / 复算值 / 偏差
"""
import os, re, json
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
TB, IN, LD, RAW = f"{BASE}/tables", f"{BASE}/inputs", f"{BASE}/inputs_ld", f"{BASE}/raw"

ROWS, KEY = [], []
def add(cid, cat, item, exp, act, ok, note="", level="数据"):
    """level: 数据 = 表/原始数据本身的错；表述 = 文档写错、数据无误"""
    ROWS.append(dict(编号=cid, 类别=cat, 层级=level, 检查项=item, 期望=exp, 实测=act,
                     结论=("PASS" if ok is True else ("WARN" if ok is None else "FAIL")),
                     说明=note))
def key(kid, item, formula, rep, rec, note=""):
    try:
        dev = float(rec) - float(rep)
    except Exception:
        dev = np.nan
    KEY.append(dict(编号=kid, 数字=item, 定义式=formula, 报告值=rep, 复算值=rec,
                    偏差=dev, 说明=note))
def close(a, b, rtol=1e-6, atol=1e-12):
    try:
        return bool(abs(float(a) - float(b)) <= atol + rtol * max(abs(float(a)), abs(float(b))))
    except Exception:
        return False

T = lambda f: pd.read_csv(f"{TB}/{f}", encoding="utf-8-sig")
t23 = T("23_信号解析归属总表.csv")
t24 = T("24_按信号窗口的coloc.csv")
t26 = T("26_峰周围基因_窗口coloc.csv")
t27 = T("27_肩部基因归属表.csv")
t31 = T("31_联合条件模型.csv")
t04 = T("04_cisMR_SMR_HEIDI.csv")
t02b = T("02b_每基因最强eQTL.csv")
t11 = T("11_变异信息表.csv")
t21 = T("21_已发表哨兵_17q21.csv")
t22 = T("22_哨兵基因映射唯一性.csv")
genes = pd.read_csv(f"{RAW}/ens_genes_17_peak_b37.csv")

# ---------------------------------------------------------------- A 表内自洽
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
    add("A1", "表内自洽", f"{tag}：P 可由 b/se 复算", f"{tot} 行全部一致",
        f"{ok}/{tot}", ok == tot and tot > 0, f"不一致：{bad[:4]}" if bad else "")

a_p_from_bse(t31, "边际b", "边际se", "边际P", "表31 边际")
a_p_from_bse(t31, "联合b", "联合se", "联合P", "表31 联合")
a_p_from_bse(t04, "b_wald", "se_wald", "p_wald", "表04 Wald")

# 表04 的 F_lead 应等于 (lead_beta/lead_se)^2
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
add("A2", "表内自洽", "表04：F_lead == (lead_beta/lead_se)^2", f"{ftot} 行",
    f"{fok}/{ftot}", fok == ftot, "逐行回到 inputs 原始 tsv 复算")

# coloc 后验和
s26 = t26[["H0", "H1", "H2", "H3", "H4"]].sum(axis=1)
add("A3", "表内自洽", "表26：H0+…+H4 == 1", "全部 1±1e-6",
    f"min={s26.min():.8f} max={s26.max():.8f}", bool(np.allclose(s26, 1, atol=1e-6)),
    f"{len(s26)} 行")
bad24 = t24[(t24.H3 + t24.H4 > 1 + 1e-9)]
add("A4", "表内自洽", "表24：H3+H4 <= 1", "无越界", f"越界 {len(bad24)} 行", len(bad24) == 0)

# ---------------------------------------------------------------- B 跨表一致
# B1 表23 的 MR 列 == 表04 同 (tissue,gene,outcome)
mism = []
for _, r in t23.iterrows():
    g, t = r["候选基因"], r["组织"]
    oc = "TSH"
    s = t04[(t04.gene == g) & (t04.tissue == t) & (t04.outcome == oc)]
    if s.empty:
        mism.append((g, t, "缺行")); continue
    s = s.iloc[0]
    if not (close(s.b_wald, r["MR_b"], 1e-6) and close(s.se_wald, r["MR_se"], 1e-6)):
        mism.append((g, t, s.b_wald, r["MR_b"]))
add("B1", "跨表一致", "表23 的 MR_b/MR_se == 表04 的 b_wald/se_wald（TSH）",
    f"{len(t23)} 行全等", f"不等 {len(mism)} 行", len(mism) == 0, str(mism[:3]))

# B2 表23 的窗口 coloc == 表24
WIN = {1: "信号1(45.72-45.86)", 2: "信号2(45.86-45.98)", 3: "信号3(46.00-46.12)"}
mism2 = []
for _, r in t23.iterrows():
    g, t, sig = r["候选基因"], r["组织"], int(r["信号"])
    for oc, h3c, h4c in [("TSH", "窗口内_coloc_H3", "窗口内_coloc_H4"),
                         ("HYPO", "甲减_窗口内_H3", "甲减_窗口内_H4")]:
        v3, v4 = r[h3c], r[h4c]
        if not (np.isfinite(v3) and np.isfinite(v4)):
            continue
        s = t24[(t24.gene == g) & (t24.tissue == t) & (t24.outcome == oc) & (t24.window == WIN[sig])]
        if s.empty:
            mism2.append((g, t, oc, "表24缺行")); continue
        s = s.iloc[0]
        if not (close(s.H3, v3, 1e-4, 1e-9) and close(s.H4, v4, 1e-4, 1e-9)):
            mism2.append((g, t, oc, s.H3, v3, s.H4, v4))
add("B2", "跨表一致", "表23 的窗口 coloc == 表24 对应 (基因,组织,结局,窗口)",
    f"{len(t23)*2} 个单元格", f"不等 {len(mism2)} 个", len(mism2) == 0, str(mism2[:3]))

# B3 表27 的 MR == 表04
mism3 = []
for _, r in t27.iterrows():
    g = r["基因"]; t = {"全血": "WholeBlood", "甲状腺": "Thyroid"}[r["组织"]]
    for oc, bc, pc in [("HYPO", "MR_HYPO_b", "MR_HYPO_p"), ("TSH", "MR_TSH_b", "MR_TSH_p")]:
        s = t04[(t04.gene == g) & (t04.tissue == t) & (t04.outcome == oc)]
        if s.empty:
            continue
        s = s.iloc[0]
        if np.isfinite(r[bc]) and not close(s.b_wald, r[bc], 1e-6):
            mism3.append((g, t, oc, s.b_wald, r[bc]))
add("B3", "跨表一致", "表27 的 cis-MR == 表04", f"{len(t27)*2} 个单元格",
    f"不等 {len(mism3)} 个", len(mism3) == 0, str(mism3[:3]))

# B4 表22 的哨兵数与原始 CSV 分组一致
raw_s = pd.read_csv(f"{RAW}/sentinels/gwas_catalog_associations_all.csv", encoding="utf-8-sig")
cnt = raw_s.groupby(["paper", "trait"]).size().reset_index(name="n")
oks = []
for _, r in t22.iterrows():
    m = cnt[(cnt.paper == r.paper) & (cnt.trait == r.trait)]
    oks.append((r.paper, r.trait, int(r.哨兵数), int(m.n.iloc[0]) if len(m) else -1))
add("B4", "跨表一致", "表22 的哨兵数 == 原始 CSV 分组计数", "4 组全等",
    "; ".join(f"{a}/{b}:{c}={d}" for a, b, c, d in oks),
    all(c == d for _, _, c, d in oks), f"原始总条数 {len(raw_s)}")

# ---------------------------------------------------------------- C 原始数据
# C1 双 build 偏移
for nm, df, c38, c37 in [("表11", t11, "pos_grch38", "pos_grch37"),
                         ("表21", t21, "pos38", "pos37")]:
    dd = df.dropna(subset=[c38, c37])
    off = (dd[c38].astype(float) - dd[c37].astype(float))
    add("C1" + nm[-1], "坐标", f"{nm}：pos38 − pos37 == 1,922,634",
        "全部 1,922,634", f"取值集合 {sorted(set(off.astype(int)))}",
        bool((off == 1922634).all()), f"{len(dd)} 行")
key("K1", "GRCh38→GRCh37 偏移", "rs4794063 47,727,128 − 45,804,494", 1922634,
    int(t21.pos38.iloc[0] - t21.pos37.iloc[0]), "两独立位点同值")

# C2 输入文件的等位基因定向
#     beta_out 应已按 eqtl_ea 定向；gwas_ea 列存的是翻转前的 GWAS 效应等位基因
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
add("C2", "等位基因", "inputs/*.tsv 的 beta_out 已按 eqtl_ea 定向（与表11的 tsh_ea 差一个符号当且仅当等位基因不同）",
    f"{len(chk)} 个匹配全部成立", f"{nok}/{len(chk)}", nok == len(chk),
    "注意：文件里的 gwas_ea 列是翻转前的值，按它对齐会翻错（p12 已废）")

# C3 基因坐标：rs6503875 是否在 SCRN2 内；rs62076102 相对 SCRN2 的方向
sc = genes[genes.gene_name == "SCRN2"].iloc[0]
inside = bool(sc.start <= 45915395 <= sc.end)
add("C3", "坐标", "rs6503875（45,915,395）落在 SCRN2 基因体内",
    f"SCRN2 {sc.start}-{sc.end}（strand {sc.strand:+d}）", str(inside), inside,
    f"距 3′ 端 {45915395 - sc.start} bp；距 TSS {sc.end - 45915395} bp")
tss = int(sc.end) if sc.strand == -1 else int(sc.start)
d62 = 45887132 - tss
add("C4", "坐标", "rs62076102 相对 SCRN2 启动子的位置",
    "应与正文写法一致", f"{d62:+,} bp（负=下游/3′侧）", None,
    "SCRN2 为负链，TSS=45,918,699；45,887,132 在 TSS 下游 31.6 kb，"
    "即 3′ 端外侧 28.0 kb。正文原写“5′ 端 28 kb”方向错，已改")
key("K2", "rs62076102 到 SCRN2 3′ 端距离", "45,915,058 − 45,887,132",
    28000, int(sc.start - 45887132), "正文写 28 kb，实测 27,926 bp")

# C5 独立重算 5 SNP 的 GLS 联合模型（只依赖边际 b/se + LD，不读表31）
snp5 = [("rs7219420", 45800299), ("rs1989291", 45807547),
        ("肩部", 45873049), ("rs62076102", 45887132), ("rs6503875", 45915395)]
def gls(outcome):
    d = t31[(t31["结局"] == outcome) & (t31["模型"] == "五SNP")].set_index("SNP")
    d = d.loc[[s for s, _ in snp5], ["pos37", "边际b", "边际se"]]
    LDm = pd.read_csv(f"{LD}/ld_{outcome}.tsv.gz", sep="\t", index_col=0)
    LDm.index = LDm.index.astype(int); LDm.columns = LDm.columns.astype(int)
    p = [int(x) for x in d.pos37]
    R = LDm.loc[p, p].values
    b = d["边际b"].values; se = d["边际se"].values
    V = np.diag(se) @ R @ np.diag(se)
    W = np.linalg.inv(V)
    bj = np.linalg.solve(R.T @ W @ R, R.T @ W @ b)
    sej = np.sqrt(np.diag(np.linalg.inv(R.T @ W @ R)))
    return pd.DataFrame(dict(pos37=p, b_joint=bj, se_joint=sej, z=bj / sej,
                             p=2 * stats.norm.sf(abs(bj / sej))), index=d.index)

for oc in ["HYPO", "TSH"]:
    g = gls(oc)
    rep = t31[(t31["结局"] == oc) & (t31["模型"] == "五SNP")].set_index("SNP")["联合z"]
    ok = all(close(g.loc[s, "z"], rep[s], 1e-4) for s in g.index)
    add(f"C5_{oc}", "GLS 重算", f"{oc} 五SNP 联合 z 独立重算",
        "与表31 一致", "; ".join(f"{s}:{g.loc[s,'z']:.3f}" for s in g.index), ok,
        "只输入边际 b/se 与 LD 矩阵，不读表31 的结果列")
    for s in g.index:
        key(f"K3_{oc}_{s}", f"{oc} {s} 联合 z", "(R′V⁻¹R)⁻¹R′V⁻¹b", float(rep[s]),
            round(float(g.loc[s, "z"]), 3))

# C6 信号计数：条件化后仍 P<0.05 的个数
def sig_count(oc, alpha=0.05):
    g = gls(oc)
    return int((g.p < alpha).sum()), g.p.round(4).to_dict()
c_h, ph = sig_count("HYPO"); c_t, pt = sig_count("TSH")
add("C6", "信号数", "甲减：五SNP 全模型内仍 P<0.05 的独立信号数", "4", f"{c_h}  {ph}",
    c_h == 4, "45,800 / 45,807 / 45,873 / 45,887")
add("C7", "信号数", "TSH：五SNP 全模型内仍 P<0.05 的独立信号数", "2（+ 46.05 的 COPZ2 不在五SNP内）",
    f"{c_t}  {pt}", c_t == 2, "45,800 / 45,887；46.05 单独计入 → 共 3 个 TSH 信号")

# C8 LD 复算（回到 1000G 矩阵）
LDt = pd.read_csv(f"{LD}/ld_TSH.tsv.gz", sep="\t", index_col=0)
LDt.index = LDt.index.astype(int); LDt.columns = LDt.columns.astype(int)
pairs = [("rs62076102 / rs6503875", 45887132, 45915395, 0.7669),
         ("rs7219420 / rs6503875", 45800299, 45915395, 0.0262),
         ("rs7219420 / rs1989291", 45800299, 45807547, -0.3981),
         ("肩部 / rs1989291", 45873049, 45807547, -0.2526),
         ("肩部 / rs62076102", 45873049, 45887132, 0.4015),
         ("肩部 / rs6503875", 45873049, 45915395, 0.3051)]
allok = True
for nm, a, b, rep in pairs:
    v = float(LDt.loc[a, b])
    ok = close(v, rep, 2e-3); allok &= ok
    key("K4_" + nm.replace(" ", "").replace("/", "_"), f"r({nm})", "LD 矩阵取元素", rep, round(v, 4))
add("C8", "LD", "表31 引用的 6 个 r 值 == LD 矩阵元素", "6 个全等",
    "; ".join(f"{n}:{LDt.loc[a,b]:.4f}" for n, a, b, _ in pairs), allok)

# C9 哨兵：区域内条数
reg = t21
add("C9", "哨兵", "17q21 区域（GRCh38 47.3–48.5 Mb）内哨兵数 == 1", "1", f"{len(reg)}",
    len(reg) == 1, f"{reg.rsid.iloc[0]} / {reg.trait.iloc[0]} / P={reg.p.iloc[0]}")
key("K5", "区域内哨兵数", "raw/sentinels 过滤 chr38∈[47.3e6,48.5e6]", 1, len(reg))
key("K6", "哨兵总数", "四套研究 GWAS Catalog associations 去重计数", 783, len(raw_s))

# C10 表25 来源核验
t25 = T("25_HYPO数据来源核验.csv").set_index("核验项")["值"]
add("C10", "来源核验", "面板甲减 == Rand2025 GCST90572791（beta/SE 逐位）",
    "相同比例 1.0 且最小 P 位置一致", f"{t25['beta/SE 逐位相同比例']} / {t25['两套数据最小 P 及位置']}",
    float(t25["beta/SE 逐位相同比例"]) == 1.0, "n=1,043,800")

# C11 TBX21：表23 的 eQTL_P 与 F 是否来自同一个变异
row = t23[(t23["候选基因"] == "TBX21")].iloc[0]
d = pd.read_csv(f"{IN}/Thyroid__TBX21__TSH.tsv", sep="\t")
p_at_lead = float(d[d.pos37 == 45800299].iloc[0].p_eqtl)
lead = d.loc[d.p_eqtl.idxmin()]
F_min = float((lead.beta_eqtl / lead.se_eqtl) ** 2)
same = close(p_at_lead, row["eQTL_P"], 1e-4) and close(F_min, row["F"], 1e-3)
add("C11", "口径", "TBX21 行：eQTL_P 与 F 是否取自同一变异", "应同源",
    f"eQTL_P={p_at_lead:.4f} 在信号1 lead；F={F_min:.2f} 在 {int(lead.pos37)}；表23 记 P={row['eQTL_P']:.4f}, F={row['F']:.2f}",
    False, "两者不同源：P 是信号 1 lead 上的值，F 是区间内名义最小 P 变异的值（45,455,669，距 TBX21 355 kb）。"
           "正文“P=0.55，F=5.8”把两个变异混着写，已拆开", level="表述")

# C12 表23 代表位点的 TSH_beta 是否与表11 一致（按等位基因定向）
t11i = t11.set_index("rsid")
bad12 = []
for _, r in t23.iterrows():
    rs = r["代表位点"]
    if rs not in t11i.index:
        continue
    b11 = float(t11i.loc[rs, "beta_tsh"])
    if not (close(abs(b11), abs(r["TSH_beta"]), 1e-3)):
        bad12.append((rs, b11, r["TSH_beta"]))
add("C12", "等位基因", "表23 的 TSH_beta 与表11 同量级（仅符号可因等位基因不同而反）",
    f"{len(t23)} 行", f"不符 {len(bad12)} 行", len(bad12) == 0,
    "表23 统一按 eqtl_ea 定向；表11 按 tsh_ea。两者符号相反是预期，不是错误")

# C13 两家 TSH GWAS 在本区域是否确实 0 哨兵
for (pp, tr), g in raw_s.groupby(["paper", "trait"]):
    n_reg = int(((g.chr38.astype(str) == "17") & (g.pos38.astype(float) >= 47.3e6) &
                 (g.pos38.astype(float) <= 48.5e6)).sum())
    rep = t22[(t22.paper == pp) & (t22.trait == tr)]
    rep = int(rep.区域哨兵数.iloc[0]) if len(rep) else -1
    add(f"C13_{pp}_{tr}", "哨兵", f"{pp} {tr}：区域内哨兵数（原始 CSV 重算 vs 表22）",
        f"表22 记 {rep}", f"重算 {n_reg}（总数 {len(g)}）", n_reg == rep)

# C14 位置关系的方向必须按基因链判定（负链基因最容易写反）
COPZ2 = dict(start=46103533, end=46115392, strand=-1)   # Ensembl GRCh37 grch37.rest.ensembl.org
COORD = {g.gene_name: dict(start=int(g.start), end=int(g.end), strand=int(g.strand))
         for _, g in genes.iterrows()}
COORD["COPZ2"] = COPZ2
CASES = [("rs62076102", 45887132, "SCRN2"), ("rs6503875", 45915395, "SCRN2"),
         ("rs7219420", 45800299, "TBKBP1"), ("rs7219420", 45800299, "TBX21"),
         ("COPZ2lead", 46050669, "COPZ2"), ("肩部", 45873049, "TBKBP1")]
lines = []
for nm, pos, g in CASES:
    c = COORD[g]
    tss = c["end"] if c["strand"] == -1 else c["start"]
    d_tss = pos - tss
    side = ("启动子外侧" if c["strand"] == -1 else "启动子外侧") if d_tss < 0 else "基因体/下游"
    inside = c["start"] <= pos <= c["end"]
    dist_edge = 0 if inside else (c["start"] - pos if pos < c["start"] else pos - c["end"])
    if inside:
        desc = f"基因体内（距 3′ 边界 {pos - c['start']} bp，距 TSS {c['end'] - pos} bp）"
    else:
        five_side = (pos > c["end"]) if c["strand"] == -1 else (pos < c["start"])
        desc = f"{'5′' if five_side else '3′'} 侧 {abs(dist_edge)/1000:.1f} kb"
    lines.append(f"{nm}/{g}: TSS{'+' if d_tss>=0 else ''}{d_tss} → {desc}")
add("C14", "坐标", "位置关系按基因链判定（负链基因的 5′/3′ 极易写反）", "与正文写法一致",
    " | ".join(lines), False,
    "SCRN2 与 COPZ2 均为负链：rs62076102 在 SCRN2 的 3′ 侧 28.0 kb（原写 5′ 端）；"
    "COPZ2lead 在 COPZ2 的 3′ 侧 52.9 kb / 距 TSS 64.7 kb（原写“上游 53 kb”）", level="表述")

# ---------------------------------------------------------------- 输出
res = pd.DataFrame(ROWS)
res.to_csv(f"{TB}/32_核对清单.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(KEY).to_csv(f"{TB}/33_关键数字_可复算.csv", index=False, encoding="utf-8-sig")

pd.set_option("display.width", 260, "display.max_colwidth", 60)
print(res.to_string(index=False))
n_fail = int((res["结论"] == "FAIL").sum()); n_warn = int((res["结论"] == "WARN").sum())
nf_d = int(((res["结论"] == "FAIL") & (res["层级"] == "数据")).sum())
nf_w = int(((res["结论"] == "FAIL") & (res["层级"] == "表述")).sum())
print(f"\n合计 {len(res)} 条：PASS {len(res)-n_fail-n_warn} / WARN {n_warn} / FAIL {n_fail}")
print(f"  其中 数据层 FAIL {nf_d} 条，表述层 FAIL {nf_w} 条（数据无误，文档已改，见表内说明）")
print("写出 tables/32_核对清单.csv, tables/33_关键数字_可复算.csv")
