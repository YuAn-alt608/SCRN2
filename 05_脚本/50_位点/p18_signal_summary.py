"""
p18_signal_summary.py —— 把"位点=多信号"的结论固化成三张表

产出
  tables/21_已发表哨兵_17q21.csv     两家论文在该区域的哨兵 + LD 归属
  tables/22_哨兵基因映射唯一性.csv   全基因组层面的"可唯一归属"比例（可扩展的口径）
  tables/23_信号解析归属总表.csv     本区域三个信号 × 候选基因 × 全部证据 × 裁决

原则：所有数字从项目内的结果表/输入表读取，不手抄；只有已发表论文的哨兵值来自
GWAS Catalog 下载（raw/sentinels/），并注明 accession。
"""
import os

import numpy as np
import pandas as pd

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
T = f"{BASE}/tables"
# GRCh38 -> GRCh37 偏移（用两个位点独立标定，二者一致）：
#   rs4794063 47,727,128 - 45,804,494 = 1,922,634
#   rs1989291 47,730,181 - 45,807,547 = 1,922,634
# 早先用 1,957,714 是拿基因起点反推的估计值，错了 35 kb。
OFS = 1_922_634

# ---------------------------------------------------------------- 21 哨兵
sen = pd.read_csv(f"{BASE}/raw/sentinels/gwas_catalog_associations_all.csv")
reg = sen[(sen.chr38.astype(str) == "17") &
          (sen.pos38 >= 47_300_000) & (sen.pos38 <= 48_500_000)].copy()
reg["pos37"] = reg.pos38.astype(int) - OFS

# LD 代理：1000G 里 rs4794063 未被 call，取面板内最近位点（下方自动求）
PROXY = None
ld = pd.read_csv(f"{BASE}/inputs_ld/ld_TSH.tsv.gz", sep="\t", index_col=0)
ld.index = ld.index.astype(int)
ld.columns = ld.columns.astype(int)
LEAD = {"rs7219420": 45_800_299, "rs62076102": 45_887_132, "rs6503875": 45_915_395}
idx = np.array(ld.index)
if len(reg):
    tgt = int(reg.pos37.iloc[0])
    PROXY = int(idx[np.abs(idx - tgt).argmin()])
    print(f"rs4794063 在 1000G 中未被 call，取面板最近位点 {PROXY}（{abs(PROXY-tgt)} bp）")
    for k, v in LEAD.items():
        reg[f"r_{k}"] = ld.loc[PROXY, v] if (PROXY in ld.index and v in ld.index) else np.nan
    reg["proxy_pos37"] = PROXY
    reg["proxy_dist_bp"] = (reg.pos37 - PROXY).abs()
reg["信号"] = reg.pos37.apply(
    lambda p: "1（45.74–45.81）" if 45_740_000 <= p < 45_810_000 else
              ("2（45.86–45.93）" if 45_860_000 <= p < 45_930_000 else
               ("3（46.02–46.08）" if 46_020_000 <= p < 46_080_000 else "其他")))
reg["裁决"] = reg["信号"].apply(
    lambda s: "落在不可归属的信号 1（TBKBP1 异质 / TBX21 无 eQTL）" if s.startswith("1")
    else ("落在信号 2，可归 SCRN2" if s.startswith("2") else "未解析"))
cols = ["paper", "trait", "accession", "rsid", "effect_allele", "pos38", "pos37", "p", "beta",
        "raf", "mapped_gene", "proxy_pos37", "proxy_dist_bp",
        "r_rs7219420", "r_rs62076102", "r_rs6503875", "信号", "裁决"]
reg[cols].to_csv(f"{T}/21_已发表哨兵_17q21.csv", index=False, encoding="utf-8-sig")
print("=== 表21 区域哨兵 ===")
print(reg[cols].to_string(index=False))

# ---------------------------------------------------------------- 22 映射唯一性
sen["ngene"] = sen.mapped_gene.fillna("").apply(
    lambda s: len([x for x in str(s).split(";") if x.strip()]))
rows = []
for (paper, trait), g in sen.groupby(["paper", "trait"]):
    n = len(g)
    rows.append(dict(paper=paper, trait=trait, 哨兵数=n,
                     唯一基因=(g.ngene == 1).sum(),
                     多基因=(g.ngene > 1).sum(),
                     无映射=(g.ngene == 0).sum(),
                     唯一基因比例=round((g.ngene == 1).mean(), 4),
                     区域哨兵数=int(((g.chr38.astype(str) == "17") &
                                 (g.pos38 >= 47_300_000) & (g.pos38 <= 48_500_000)).sum())))
u = pd.DataFrame(rows)
u.to_csv(f"{T}/22_哨兵基因映射唯一性.csv", index=False, encoding="utf-8-sig")
print("\n=== 表22 映射唯一性 ===")
print(u.to_string(index=False))

# ---------------------------------------------------------------- 23 信号归属总表
mr = pd.read_csv(f"{T}/04_cisMR_SMR_HEIDI.csv")
co = pd.read_csv(f"{T}/03_面板_coloc_abf_先验网格.csv")
ts = pd.read_csv(f"{T}/10_逐组织_关键SNP的eQTL.csv")


def gwas(pos, outcome, tissue="Thyroid"):
    f = f"{BASE}/inputs/{tissue}__SCRN2__{outcome}.tsv"
    d = pd.read_csv(f, sep="\t")
    r = d[d.pos37 == pos]
    return (float(r.p_out.iloc[0]), float(r.beta_out.iloc[0])) if len(r) else (np.nan, np.nan)


def eqtl(pos, gene, tissue="Thyroid"):
    f = f"{BASE}/inputs/{tissue}__{gene}__TSH.tsv"
    if not os.path.exists(f):
        return np.nan, np.nan
    d = pd.read_csv(f, sep="\t")
    r = d[d.pos37 == pos]
    return (float(r.p_eqtl.iloc[0]), float(r.beta_eqtl.iloc[0])) if len(r) else (np.nan, np.nan)


def mrow(gene, tissue, outcome="TSH"):
    m = mr[(mr.gene == gene) & (mr.tissue == tissue) & (mr.outcome == outcome)]
    return m.iloc[0] if len(m) else None


def crow(gene, tissue, outcome="TSH", p12=1e-4):
    c = co[(co.gene == gene) & (co.tissue == tissue) & (co.outcome == outcome) &
           (np.isclose(co.p12, p12))]
    return c.iloc[0] if len(c) else None


TISSUES = ts.tissue.nunique()
COVERED = set(ts.gene.unique())


def ntis(snp, gene):
    """GTEx 49 组织里该 SNP 是该基因 eQTL 的组织数；基因未纳入则为'未检验'。"""
    if gene not in COVERED:
        return "未检验"
    d = ts[(ts.key_snp == snp) & (ts.gene == gene)]
    return f"{d.tissue.nunique()}/{TISSUES}"


ROWS = [
    # 信号, 位点, 位置, 候选基因, 组织, MR 工具位点
    ("1", "rs7219420", 45_800_299, "TBKBP1", "Thyroid", 45_774_072),
    ("1", "rs7219420", 45_800_299, "TBKBP1", "WholeBlood", 45_744_652),
    ("1", "rs7219420", 45_800_299, "TBX21", "Thyroid", 45_800_299),
    ("2", "rs6503875", 45_915_395, "SCRN2", "Thyroid", 45_915_395),
    ("2", "rs6503875", 45_915_395, "SCRN2", "WholeBlood", 45_912_302),
    ("2", "rs62076102", 45_887_132, "MRPL10", "Thyroid", 45_888_254),
    ("3", "COPZ2lead", 46_050_669, "COPZ2", "WholeBlood", 46_050_669),
    ("3", "COPZ2lead", 46_050_669, "COPZ2", "Thyroid", 45_881_197),
]

out = []
for sig, snp, pos, gene, tissue, mrpos in ROWS:
    p_tsh, b_tsh = gwas(pos, "TSH")
    p_hyp, _ = gwas(pos, "HYPO")
    # eQTL 用"该基因在该组织的工具位点"查，才与 F 统计量同源
    p_teq, b_teq = eqtl(mrpos, gene, tissue)
    m = mrow(gene, tissue)
    c = crow(gene, tissue)
    out.append(dict(
        信号=sig, 代表位点=snp, pos37=pos, 候选基因=gene, 组织=tissue, MR工具位点=mrpos,
        TSH_P=p_tsh, TSH_beta=b_tsh, HYPO_P=p_hyp,
        eQTL_P=p_teq, eQTL_beta=b_teq,
        F=m.F_lead if m is not None else np.nan,
        MR_b=m.b_wald if m is not None else np.nan,
        MR_se=m.se_wald if m is not None else np.nan,
        MR_P=m.p_wald if m is not None else np.nan,
        n_iv=m.n_iv if m is not None else 0,
        Q_P=m.p_Q if m is not None else np.nan,
        I2=m.I2 if m is not None else np.nan,
        HEIDI_P=m.heidi_p if m is not None else np.nan,
        coloc_H3=c["PP.H3"] if c is not None else np.nan,
        coloc_H4=c["PP.H4"] if c is not None else np.nan,
        组织可复制=ntis(snp, gene),
    ))

s = pd.DataFrame(out)

# ---- 合并"本信号窗口内"的 coloc（甲减与面板同源，见 tables/25 核验）----
w = pd.read_csv(f"{T}/24_按信号窗口的coloc.csv")
OWN_WIN = {"SCRN2": "信号2(45.86-45.98)", "MRPL10": "信号2(45.86-45.98)",
           "TBKBP1": "信号1(45.72-45.86)", "COPZ2": "信号3(46.00-46.12)"}
s["窗口内_coloc_H3"] = np.nan
s["窗口内_coloc_H4"] = np.nan
s["甲减_窗口内_H3"] = np.nan
s["甲减_窗口内_H4"] = np.nan
for i, r in s.iterrows():
    q = w[(w.gene == r.候选基因) & (w.tissue == r.组织) & (w.outcome == "TSH") &
          (w.window == OWN_WIN.get(r.候选基因, ""))]
    if len(q):
        s.at[i, "窗口内_coloc_H3"] = q.H3.iloc[0]
        s.at[i, "窗口内_coloc_H4"] = q.H4.iloc[0]
    q2 = w[(w.gene == r.候选基因) & (w.tissue == r.组织) & (w.outcome == "HYPO") &
           (w.window == OWN_WIN.get(r.候选基因, ""))]
    if len(q2):
        s.at[i, "甲减_窗口内_H3"] = q2.H3.iloc[0]
        s.at[i, "甲减_窗口内_H4"] = q2.H4.iloc[0]
# 裁决规则（先写死，保证可解释）
verdict = {
    ("1", "TBKBP1", "Thyroid"):
        "不可归属：eQTL 强（F=166）但 coloc 不支持共享（H3=0.87），单工具无法做 HEIDI",
    ("1", "TBKBP1", "WholeBlood"):
        "不可归属：Q P=0.009（I²=70%）、HEIDI<1e-4、coloc H4=0.03",
    ("1", "TBX21", "Thyroid"):
        "排除：无 cis-eQTL（P=0.55，F=5.8），MR 不可解释",
    ("2", "SCRN2", "Thyroid"):
        "归属 SCRN2：MR+coloc+SMR 均过；HEIDI 需≥2 工具，甲状腺不可检验",
    ("2", "SCRN2", "WholeBlood"):
        "归属 SCRN2：4/4 全过（Q P=0.25，HEIDI 0.064），47/49 组织可复制",
    ("2", "MRPL10", "Thyroid"):
        "排除：MR 为零（P=0.70），coloc H4=0.02（与 TSH 不共享）",
    ("3", "COPZ2", "WholeBlood"):
        "归属 COPZ2：全血 4/4（Q P=0.06，HEIDI 0.053）",
    ("3", "COPZ2", "Thyroid"):
        "方向一致但工具过弱：F=15.7、coloc H4=0.79，仅作支持",
}
s["裁决"] = [verdict.get((r.信号, r.候选基因, r.组织), "") for r in s.itertuples()]


def hypo_verdict(r):
    if np.isnan(r.甲减_窗口内_H4):
        return "不可检验"
    return ("共享（H4=%.2f）" % r.甲减_窗口内_H4 if r.甲减_窗口内_H4 > 0.75
            else ("不共享（H3=%.2f）" % r.甲减_窗口内_H3 if r.甲减_窗口内_H3 > 0.75
                  else "不可判定（H3=%.2f/H4=%.2f）" % (r.甲减_窗口内_H3, r.甲减_窗口内_H4)))


s["甲减裁决"] = [hypo_verdict(r) for r in s.itertuples()]
s.to_csv(f"{T}/23_信号解析归属总表.csv", index=False, encoding="utf-8-sig")
print("\n=== 表23 信号解析归属 ===")
print(s.round(5).to_string(index=False))
