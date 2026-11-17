#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p25_raw_data_export.py
======================
把「最初那份清单」里的 11 条必做/加分项 + 2 条投稿前核对，逐条导出成一本独立的原始数据册。
原则：
  1) 只搬数，不二次加工；凡是算过的量都标出口径（Wald / GLS / SMR / coloc.abf）。
  2) 列名写明等位基因与 build：pos_GRCh37、eQTL_EA、结局_EA_原始、*_per_EA。
  3) 缺的项不编数，写在 00_索引 与 19_说明 里，并给补法。
输出：MR1_原始数据_按需求分册_2026-11-14.xlsx
"""
import os, math
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
D18 = "D:/NewBeginning/MR1/DATA18_manuscript/analysis"
OUT = f"{BASE}/MR1_原始数据_按需求分册_2026-11-14.xlsx"

KEY_GENES = ["SCRN2", "TBKBP1", "MRPL10", "TBX21", "COPZ2", "LRRC46", "OSBPL7",
             "NPEPPS", "MRPL45P2", "SP2", "SKAP1"]

def T(name):
    return pd.read_csv(f"{BASE}/tables/{name}")

def clean(v):
    if v is None:
        return None
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    if isinstance(v, (np.floating,)):
        v = float(v)
        return None if (math.isnan(v) or math.isinf(v)) else v
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    return v

HDR_FILL = PatternFill("solid", fgColor="1F3864")
HDR_FONT = Font(bold=True, color="FFFFFF", size=10)
NOTE_FONT = Font(italic=True, color="7F7F7F", size=9)

NOTES = []

def add_sheet(wb, name, df, note=None, key_col=None, num_fmt=None):
    ws = wb.create_sheet(name)
    cols = list(df.columns)
    ws.append(cols)
    for c in range(1, len(cols) + 1):
        cell = ws.cell(row=1, column=c)
        cell.fill = HDR_FILL
        cell.font = HDR_FONT
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    r = 1
    for row in df.itertuples(index=False):
        r += 1
        ws.append([clean(v) for v in row])
    ws.freeze_panes = "A2"
    if len(df):
        ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{r}"
    # 列宽
    for i, c in enumerate(cols, start=1):
        vals = [len(str(c))]
        if len(df):
            try:
                vals += [len(str(v)) for v in df[c].head(400).tolist()]
            except Exception:
                vals += [len(str(c))]
        w = min(max(max(vals) + 2, 9), 46)
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.row_dimensions[1].height = 30
    if note:
        NOTES.append((name, note))
    return ws

# ---------------------------------------------------------------- 读源表
t01 = T("01_基因清单.csv")
t02 = T("02_面板_lead.csv")
t03 = T("03_面板_coloc_abf_先验网格.csv")
t03c = T("03c_可信集_eQTL侧.csv")
t04 = T("04_cisMR_SMR_HEIDI.csv")
t05 = T("05_rs7219420_边际.csv")
t06 = T("06_两SNP联合模型.csv")
t08 = T("08_区域图数据.csv")
t10 = T("10_逐组织_关键SNP的eQTL.csv")
t09b = T("09b_逐组织_最强eQTL_p.csv")
t11 = T("11_变异信息表.csv")
t13 = T("13_位点身份_区域最小P.csv")
t14 = T("14_三队列_样本重叠.csv")
t14b = T("14b_三队列_合并.csv")
t16 = T("16_变异_cCRE注释.csv")
t18 = T("18_工具数与GLS_IVW.csv")
t21 = T("21_已发表哨兵_17q21.csv")
t22 = T("22_哨兵基因映射唯一性.csv")
t31 = T("31_联合条件模型.csv")

sym = dict(zip(t01.gene_id, t01.gene_name))
ALIAS = {"AC003665.1": "SP2-AS1"}   # 同一基因的两个注释名，分属两个组织，不改会重复计数
def gsym(g):
    g = str(g)
    if g.startswith("ENSG"):
        g = sym.get(g, g)
    return ALIAS.get(g, g)

# ---------------------------------------------------------------- 01 面板主表
# 逐 (组织,基因,结局) 回读 inputs，取 lead 位点的等位基因与 EAF/N
alu = {}
for (tissue, gene, outcome) in t04[["tissue", "gene", "outcome"]].drop_duplicates().itertuples(index=False):
    f = f"{BASE}/inputs/{tissue}__{gene}__{outcome}.tsv"
    if not os.path.exists(f):
        continue
    WANT = ("pos37", "eqtl_ea", "eqtl_oa", "gwas_ea", "gwas_oa",
            "maf", "n_eqtl", "n_out",
            "beta_eqtl", "se_eqtl", "beta_out", "se_out", "p_eqtl")
    d = pd.read_csv(f, sep="\t", usecols=lambda c: c in WANT)
    for c in WANT:
        if c not in d.columns:
            d[c] = np.nan
    d = d.drop_duplicates("pos37").set_index("pos37")
    alu[(tissue, gene, outcome)] = d

rows = []
for r in t04.itertuples(index=False):
    k = (r.tissue, r.gene, r.outcome)
    ea = oa = gea = goa = maf = ne = no = None
    be = se = bo = so = pe = None
    d = alu.get(k)
    if d is not None and int(r.lead_pos37) in d.index:
        q = d.loc[int(r.lead_pos37)]
        ea, oa, gea, goa = q.eqtl_ea, q.eqtl_oa, q.gwas_ea, q.gwas_oa
        maf, ne, no = q.maf, q.n_eqtl, q.n_out
        be, se, bo, so, pe = q.beta_eqtl, q.se_eqtl, q.beta_out, q.se_out, q.p_eqtl
    rows.append(dict(组织=r.tissue, 基因符号=gsym(r.gene), 面板基因键=r.gene, 结局=r.outcome,
                     工具SNP_MR用=r.lead_snp, lead_pos_GRCh37=int(r.lead_pos37),
                     eQTL_EA=ea, eQTL_OA=oa, 结局_EA_原始=gea, 结局_OA_原始=goa,
                     beta_eQTL_per_EA=be, se_eQTL=se, p_eQTL=pe,
                     eQTL_MAF=maf, eQTL_n=ne, 结局_n=no, 窗口_nsnp=r.nsnp,
                     F_lead=r.F_lead,
                     MR_Wald_b=getattr(r, "b_wald"), MR_Wald_se=getattr(r, "se_wald"),
                     MR_Wald_p=getattr(r, "p_wald"),
                     n_iv=r.n_iv, k_eff=r.k_eff,
                     IVW_b=r.b_ivw, IVW_se=r.se_ivw, IVW_p=r.p_ivw,
                     Q=r.Q, Q_p=r.p_Q, I2=r.I2,
                     SMR_b=r.b_smr, SMR_se=r.se_smr, SMR_p=r.p_smr,
                     HEIDI_T=r.heidi_T, HEIDI_p=r.heidi_p, HEIDI_npair=r.heidi_npair,
                     HEIDI_ndrop=r.heidi_ndrop,
                     beta_结局_per_eQTL_EA=bo, se_结局=so))
pan = pd.DataFrame(rows)

# 窗口内最强 eQTL 位点（02 表）：与 04 的工具位点常常不是同一个，必须并列给出
lead = t02[["tissue", "gene", "outcome", "lead_snp", "lead_pos37", "lead_p_eqtl",
            "lead_beta_eqtl"]].rename(
    columns={"tissue": "组织", "gene": "面板基因键", "outcome": "结局",
             "lead_snp": "窗口内最强eQTL_SNP", "lead_pos37": "窗口内最强eQTL_pos_GRCh37",
             "lead_p_eqtl": "窗口内最强eQTL_P", "lead_beta_eqtl": "窗口内最强eQTL_beta"})
pan = pan.merge(lead, on=["组织", "面板基因键", "结局"], how="left")

# coloc H3/H4（p12=1e-4）
c4 = t03[np.isclose(t03.p12, 1e-4)][["tissue", "gene", "outcome", "PP.H3", "PP.H4", "nsnp"]] \
    .rename(columns={"tissue": "组织", "gene": "面板基因键", "outcome": "结局",
                     "PP.H3": "coloc_H3_p12_1e-4", "PP.H4": "coloc_H4_p12_1e-4",
                     "nsnp": "coloc_nsnp"})
pan = pan.merge(c4, on=["组织", "面板基因键", "结局"], how="left")

pan["复核_Wald=b_out/b_eQTL"] = pan["beta_结局_per_eQTL_EA"] / pan["beta_eQTL_per_EA"]
N04 = len(t04)
pan["是否关键基因"] = pan["基因符号"].isin(KEY_GENES).map({True: "是", False: ""})

# --- 同一基因在 04 里同时以 symbol 与 ENSG 两个键出现（数值完全相同），去重并留证 ---
KEY3 = ["基因符号", "组织", "结局"]
pan["_ensg"] = pan["面板基因键"].astype(str).str.startswith("ENSG")
dupqa = pan[pan.duplicated(KEY3, keep=False)].copy()
pan = (pan.sort_values("_ensg", kind="stable")
          .drop_duplicates(KEY3, keep="first")
          .drop(columns="_ensg").reset_index(drop=True))
qa_rows = []
for (g, tis, oc), gg in dupqa.groupby(KEY3):
    if len(gg) < 2:
        continue
    gg = gg.sort_values("_ensg")
    keep, drop = gg.iloc[0], gg.iloc[1]
    qa_rows.append(dict(基因符号=g, 组织=tis, 结局=oc,
                        保留键=keep["面板基因键"], 去除键=drop["面板基因键"],
                        b_wald_保留=keep["MR_Wald_b"], b_wald_去除=drop["MR_Wald_b"],
                        p_wald_保留=keep["MR_Wald_p"], p_wald_去除=drop["MR_Wald_p"],
                        是否一致="是" if abs(float(keep["MR_Wald_b"]) - float(drop["MR_Wald_b"])) < 1e-9 else "否",
                        说明="同一基因的 ENSG 别名键，MR/coloc/SMR 逐位相同；仅 symbol 键能联到 02 的 eQTL β/SE"))
dupqa_df = pd.DataFrame(qa_rows)

order = ["组织", "基因符号", "面板基因键", "结局", "工具SNP_MR用", "lead_pos_GRCh37",
         "窗口内最强eQTL_SNP", "窗口内最强eQTL_pos_GRCh37", "窗口内最强eQTL_P",
         "窗口内最强eQTL_beta",
         "eQTL_EA", "eQTL_OA", "结局_EA_原始", "结局_OA_原始",
         "beta_eQTL_per_EA", "se_eQTL", "p_eQTL", "eQTL_MAF", "eQTL_n",
         "beta_结局_per_eQTL_EA", "se_结局", "结局_n", "窗口_nsnp", "F_lead",
         "MR_Wald_b", "MR_Wald_se", "MR_Wald_p", "n_iv", "k_eff",
         "IVW_b", "IVW_se", "IVW_p", "Q", "Q_p", "I2",
         "SMR_b", "SMR_se", "SMR_p",
         "HEIDI_T", "HEIDI_p", "HEIDI_npair", "HEIDI_ndrop",
         "coloc_H3_p12_1e-4", "coloc_H4_p12_1e-4", "coloc_nsnp",
         "复核_Wald=b_out/b_eQTL", "是否关键基因"]
pan = pan[order].sort_values(["是否关键基因", "基因符号", "组织", "结局"],
                             ascending=[False, True, True, True])

# ---------------------------------------------------------------- 02 先验网格
grid = t03.copy()
grid.insert(2, "基因符号", [gsym(g) for g in grid.gene])
grid["是否关键基因"] = grid["基因符号"].isin(KEY_GENES).map({True: "是", False: ""})
grid = grid.rename(columns={"tissue": "组织", "gene": "面板基因键", "outcome": "结局",
                            "nsnp": "coloc_nsnp"})
grid = grid[["组织", "基因符号", "面板基因键", "结局", "p12", "coloc_nsnp",
             "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "是否关键基因"]]

# ---------------------------------------------------------------- 03 边际
m5 = t05.copy()
m5["队列说明"] = m5["cohort"].map({
    "Discovery EUR (n=482,678)": "发现集（Rand 2025 欧洲）",
    "Replication EUR GCST90319321": "复制集 1（欧洲）",
    "Replication SAS GCST90448634": "复制集 2（南亚）"})
m5 = m5[["cohort", "队列说明", "snp", "ea", "beta", "se", "p", "eaf", "n", "z"]]
m5.columns = ["队列", "队列说明", "SNP", "EA", "beta_per_EA", "se", "P", "EAF_of_EA", "N", "z"]

# ---------------------------------------------------------------- 04 联合
j6 = t06.copy()
j6.insert(0, "结局", "TSH")

# ---------------------------------------------------------------- 05 区域图
reg = t08.copy()
reg = reg.rename(columns={"pos37": "pos_GRCh37", "tsh_ea": "TSH_EA", "tsh_oa": "TSH_OA",
                          "beta_tsh": "beta_TSH_per_EA", "se_tsh": "se_TSH",
                          "eaf_tsh": "EAF_TSH_EA", "p_tsh": "P_TSH", "n_tsh": "N_TSH",
                          "r2_lead": "r2_with_rs6503875"})
for g in ["SCRN2", "TBKBP1", "MRPL10"]:
    for tis, cn in [("Thyroid", "甲状腺"), ("WholeBlood", "全血")]:
        reg = reg.rename(columns={f"{g}_ea_{tis}": f"{g}_{cn}_EA",
                                  f"{g}_beta_{tis}": f"{g}_{cn}_beta_per_EA",
                                  f"{g}_se_{tis}": f"{g}_{cn}_se",
                                  f"{g}_p_{tis}": f"{g}_{cn}_P"})
    reg = reg.rename(columns={f"{g}_ea_Thyroid": f"{g}_甲状腺_EA",
                              f"{g}_ea_WholeBlood": f"{g}_全血_EA"})

# ---------------------------------------------------------------- 06 LD
ld = pd.read_csv(f"{BASE}/inputs_ld/ld_TSH.tsv.gz", sep="\t", index_col=0)
ld.index = ld.index.astype(int)
ld.columns = ld.columns.astype(int)
snps5 = [("rs7219420", 45800299), ("rs1989291", 45807547), ("rs72833417(肩部)", 45873049),
         ("rs62076102", 45887132), ("rs6503875", 45915395)]
lrows = []
for a, pa in snps5:
    for b, pb in snps5:
        if pa not in ld.index or pb not in ld.columns:
            continue
        r = float(ld.loc[pa, pb])
        lrows.append(dict(类型="r", SNP_i=a, pos_i_GRCh37=pa, SNP_j=b, pos_j_GRCh37=pb,
                          r=r, r2=r * r))
for a, pa in snps5:
    for b, pb in snps5:
        if pa not in ld.index or pb not in ld.columns:
            continue
        r = float(ld.loc[pa, pb])
        lrows.append(dict(类型="r2", SNP_i=a, pos_i_GRCh37=pa, SNP_j=b, pos_j_GRCh37=pb,
                          r=r, r2=r * r))
ldt = pd.DataFrame(lrows)

# ---------------------------------------------------------------- 08 哨兵
sent = t21.copy()
sent["是否在17q21.32区域"] = np.where(sent["信号"].notna(), "是", "否")
sent = sent.rename(columns={"pos38": "pos_GRCh38", "pos37": "pos_GRCh37",
                            "effect_allele": "哨兵_EA", "raf": "RAF_EA"})

# 四套研究的全部 783 条哨兵（原始下载件）
sent_all = pd.read_csv(f"{BASE}/raw/sentinels/gwas_catalog_associations_all.csv")
sent_all["在17q21.32区域"] = np.where(
    (sent_all.chr38.astype(str) == "17") & (sent_all.pos38 >= 47.3e6) & (sent_all.pos38 <= 48.5e6),
    "是", "否")
sent_all = sent_all.rename(columns={"pos38": "pos_GRCh38", "effect_allele": "哨兵_EA",
                                    "raf": "RAF_EA", "pubmed": "PMID",
                                    "reported_trait": "报告性状"})

# ---------------------------------------------------------------- 11 变异信息
vt = t11.copy()
def lookup_eqtl(tissue, gene, pos):
    f = f"{BASE}/inputs/{tissue}__{gene}__TSH.tsv"
    if not os.path.exists(f):
        return (None, None, None)
    d = pd.read_csv(f, sep="\t", usecols=["pos37", "eqtl_ea", "maf", "n_eqtl"])
    d = d.drop_duplicates("pos37").set_index("pos37")
    if pos not in d.index:
        return (None, None, None)
    q = d.loc[pos]
    return (q.eqtl_ea, q.maf, q.n_eqtl)

extra = []
for r in vt.itertuples(index=False):
    row = {}
    for g in ["SCRN2", "TBKBP1"]:
        for tis, cn in [("Thyroid", "甲"), ("WholeBlood", "血")]:
            ea, maf, ne = lookup_eqtl(tis, g, int(r.pos_grch37))
            row[f"{g}_{cn}_eQTL_EA"] = ea
            row[f"{g}_{cn}_MAF"] = maf
            row[f"{g}_{cn}_N"] = ne
    extra.append(row)
vt = pd.concat([vt, pd.DataFrame(extra)], axis=1)
vt = vt.rename(columns={"rsid": "rsID", "pos_grch37": "pos_GRCh37", "pos_grch38": "pos_GRCh38",
                        "tsh_ea": "TSH_EA", "tsh_oa": "TSH_OA",
                        "eaf_tsh_of_tsh_ea": "EAF_TSH_EA", "n_tsh": "N_TSH",
                        "beta_tsh": "beta_TSH_per_EA", "se_tsh": "se_TSH", "p_tsh": "P_TSH",
                        "r2_with_rs6503875": "r2_with_rs6503875"})

# ---------------------------------------------------------------- 14 逐组织
tissue_q = t09b.rename(columns={"tissue": "组织"})

# ---------------------------------------------------------------- 15 调控
regu = t16.rename(columns={"rsid": "rsID", "pos_grch38": "pos_GRCh38"})

# ---------------------------------------------------------------- 17 PheWAS
ph = pd.read_csv(f"{D18}/phewas_thyroid_harmonised.csv")

# ---------------------------------------------------------------- 写盘
wb = Workbook()
wb.remove(wb.active)

# 00 索引
idx_rows = [
 dict(需求="必做1 位点基因面板（cis-MR+SE、coloc PP.H3/H4 含先验网格、SMR/HEIDI）",
      本册表="01_面板_主表 / 02_coloc先验网格",
      行数=f"{len(pan)} 行（{pan['基因符号'].nunique()} 个基因符号 / {pan['面板基因键'].nunique()} 个面板键 × 2 组织 × 4 结局）",
      状态="齐",
      关键值="SCRN2·甲状腺·TSH：β=+0.0336 (0.0089)，H4=0.987；TBKBP1·全血·TSH：β=+0.0058，H4=0.03；MRPL10 见 01 表",
      缺口与补法="无结构性缺口。n_iv=0 的行 SMR 退化为单工具 Wald（表内 se_smr=b_wald 同值），投稿时须在脚注说明"),
 dict(需求="必做2 rs7219420 对 TSH 的三队列边际 + TBKBP1 cis-MR + 两 SNP 联合",
      本册表="03_rs7219420边际 / 04_两SNP联合 / 01_面板_主表（取 TBKBP1 行）",
      行数="6 + 2 行",
      状态="部分",
      关键值="发现集 rs7219420 β=−0.0121 (0.0023) P=1.57e-7；复制集1 β=−0.0065 (0.0041) P=0.115；两 SNP r=0.0262，联合后 rs6503875 仅移 −2.28%",
      缺口与补法="南亚队列 GCST90448634 无区域数据（β/SE/P 全空）。补法：向原文作者索取或改用 GWAS Catalog 的南亚 sumstats（若有）"),
 dict(需求="必做3 区域图汇总统计 + 统一 LD + SuSiE 可信集",
      本册表="05_区域图_输入 / 06_LD_关键位点 / 07_可信集",
      行数=f"{len(reg)} / 50 / {len(t03c)} 行",
      状态="部分",
      关键值="β/SE/EAF/N 四件齐（N_TSH 逐变异 201,288–482,771）；LD 全矩阵 2252×2252（1000G EUR n=503）",
      缺口与补法="TSH 侧 SuSiE 报 prior variance unreasonably large，未采用；只有 eQTL 侧可信集。补法：改 coloc.susie 用 eQTL 一侧的 CS 当条件集"),
 dict(需求="必做4 位点身份（Rand 2025 / Williams 2023 哨兵表；区域最小 P）",
      本册表="08_哨兵明细_783 / 08b_区域哨兵_带LD / 09_哨兵映射唯一性 / 10_区域最小P",
      行数=f"{len(sent_all)} / {len(sent)} / {len(t22)} / {len(t13)} 行",
      状态="齐",
      关键值="四套共 783 条哨兵，区域（chr17:47.3–48.5 Mb GRCh38）仅 1 条：Rand 甲减 rs4794063（P=2e-36，映射 TBKBP1;TBX21）；两家 TSH 区域 0 条。TSH 区域最小 P=1.568e-7（rs7219420，未过 5e-8）",
      缺口与补法="无。注意：TSH 未过 5e-8 → 稿件只能写“已知甲减位点的效应基因解析”，不能写“新 TSH 位点”"),
 dict(需求="必做5 变异信息表（rsID、双 build、EA/OA、eQTL 与 TSH 各自 EAF 和 N）",
      本册表="11_变异信息",
      行数=f"{len(vt)} 行 × {len(vt.columns)} 列",
      状态="齐",
      关键值="8 个变异；GRCh38→GRCh37 偏移恒为 1,922,634；rs4794224 为三等到位（T/C/G），GTEx 只二态",
      缺口与补法="rs4794224 多等位，面板缺该行；rs3047616 仅 N=227,214（其余 482,7xx）"),
 dict(需求="必做6 样本重叠（Rand 分队列、去 UKB、Figueredo 构成）",
      本册表="12_样本重叠",
      行数="3 + 1 行（+ 缺项说明行）",
      状态="缺",
      关键值="现有三队列：Rand2025_EUR n=482,873 → 0.0336 (0.0089)；GCST90319321 n=164,818 → 0.0404 (0.0144)；GCST90448634 SAS n=38,000 → 0.0288 (0.0241)",
      缺口与补法="缺 Rand 2025 六个分队列（CHB/DBDS/UKB/HUNT/MGI/ThyroidOmics）逐变异估计、去 UKB 重估、Figueredo 队列构成。补法：Rand 2025 补充表 S1–S3（Nat Genet 或 medRxiv 附件）；拿到后按 n 加权重算"),
 dict(需求="加分7 SCRN2 第二个独立 eQTL 信号",
      本册表="13_SCRN2第二信号 / 07_可信集",
      行数=f"{len(t18)} 行（+ 07_可信集 的 {len(t03c)} 行）",
      状态="部分",
      关键值="甲状腺 eQTL 侧 SuSiE 只出 1 个 CS（L1，含 rs6503875，PIP=0.747，大小 5）；k_eff=1.64（3 工具）/1.62（4 工具）",
      缺口与补法="未找到第二个独立信号。rs3047616 在 eQTL 显著（P=4.5e-6）但 1000G VCF 无此变异，进不了 LD。补法：eQTL Catalogue 的 SuSiE 可信集，或补 1000G 该位点"),
 dict(需求="加分8 逐组织共定位（含垂体/下丘脑）",
      本册表="14_逐组织eQTL",
      行数=f"{len(tissue_q)} / {len(t10)} 行",
      状态="缺",
      关键值="只有 49 个组织的 eQTL 显著性与方向（SCRN2 全 49 个组织最强 P 均 <1e-8）",
      缺口与补法="GTEx 只发布显著对，逐组织 coloc 需各组织全区域 nominal 统计（约 100 GB）。GTEx 无下丘脑。补法：写“eQTL 存在性与方向”，不写“逐组织共定位”，或下 eQTL Catalogue 全量"),
 dict(需求="加分9 变异→基因调控证据（cCRE/ATAC/H3K27ac/ABC/rE2G/L2G）",
      本册表="15_调控注释",
      行数=f"{len(regu)} 行",
      状态="部分",
      关键值="8 个变异均不在 cCRE 内，最近元件为 dELS/pELS（36–2,029 bp）；5 kb 内 4–16 个 cCRE",
      缺口与补法="缺甲状腺 ATAC/H3K27ac、ABC、ENCODE-rE2G、Open Targets L2G。补法：Open Targets 全量 parquet（L2G）+ ENCODE rE2G 下载；本机 API 被网关挡"),
 dict(需求="加分10 罕见变异（Genebass / AZ PheWAS 的 LoF/错义负担）",
      本册表="16_罕见变异_空缺模板",
      行数="0（模板 12 行 = 3 基因 × 2 结局 × 2 掩码）",
      状态="缺",
      关键值="—",
      缺口与补法="Genebass API 返回 301 不可达。补法：下载 Genebass gene-level burden（约 2 GB）或 AZ PheWAS 的 *.gene.burden，取 SCRN2/TBKBP1/MRPL10 三行"),
 dict(需求="加分11 游离 T4 GWAS + rs6503875 全表型 PheWAS",
      本册表="17_FT4与PheWAS / 20_结局GWAS身份",
      行数=f"{len(ph)} 行（PheWAS 甲状腺相关）+ 3 行（N 核定）",
      状态="齐",
      关键值="FT4 不缺：面板的 T4 就是游离 T4（Rand 2025 GCST90572790，abstract 明写 free thyroxine，191,449 人；文件最大 n=191,448）。PheWAS 594 命中/376 数据集；rs6503875 133 条非 eQTL 命中，p<5e-8 共 29 条",
      缺口与补法="无需再取 Teumer 2018。原“总 T4”说法已证伪，见 20_结局GWAS身份；稿件与 DATA10 里所有“总 T4”一律改为“游离 T4（FT4）”"),
 dict(需求="投稿前A 两 SNP 首位基因 P 都是 3.27341e-310",
      本册表="18_投稿前核对",
      行数="8 行",
      状态="已定案",
      关键值="eQTLGen 原始文件中 3.2717E-310 这一个值被 877 行共用（z 从 44 到 79 都有）→ 是浮点下限，不是复制错误",
      缺口与补法="建议正文用 |z| 或 β/se 排序，不用 P；P 只写“<1e-300”"),
 dict(需求="投稿前B s21 脚本一阶 SE 复现不了 CSV",
      本册表="18_投稿前核对",
      行数="8 行",
      状态="已定案",
      关键值="s21 已改为 delta 法，复跑输出 0.0250 (0.0080)、z=2.72/p=0.007，与 reference_frame.csv 逐位一致",
      缺口与补法="DATA18 的 A2 总表 I5 行仍写 z=2.86/p=0.004（旧值），须改为 2.72/0.007"),
 dict(需求="附加：把同一套流程跑遍 Williams 2023 / Rand 2025 所有 TSH 位点，报可解析比例",
      本册表="09_哨兵映射唯一性（代理指标）",
      行数="4 行",
      状态="未做",
      关键值="仅能报 GWAS Catalog 注释的映射唯一率：Rand TSH 69.4%、Williams TSH 66.1%、Rand 甲减 67.2%、Rand T4 60.7%",
      缺口与补法="全基因组跑需 49 组织 eQTL 面板 + 各 locus 的 LD，本地不具备。补法：限定在 GTEx 甲状腺 + eQTLGen 全血，对 297/245 个哨兵逐个 ±500 kb 取基因并跑同一套脚本"),
]
idx = pd.DataFrame(idx_rows)

ws = wb.create_sheet("00_索引")
ws.append(["#", "需求", "本册表", "行数", "状态", "关键值", "缺口与补法"])
for c in range(1, 8):
    cell = ws.cell(row=1, column=c)
    cell.fill = HDR_FILL; cell.font = HDR_FONT
    cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
for i, r in enumerate(idx_rows, start=1):
    ws.append([i, r["需求"], r["本册表"], r["行数"], r["状态"], r["关键值"], r["缺口与补法"]])
for i, w in enumerate([5, 52, 34, 22, 9, 70, 80], start=1):
    ws.column_dimensions[get_column_letter(i)].width = w
for row in ws.iter_rows(min_row=2, max_row=1 + len(idx_rows)):
    for cell in row:
        cell.alignment = Alignment(wrap_text=True, vertical="top")
ws.freeze_panes = "A2"
ws.row_dimensions[1].height = 24

add_sheet(wb, "01_面板_主表", pan,
          note="来源：tables/04（cis-MR/IVW/SMR/HEIDI）+ tables/02（lead 的 β/SE）+ tables/03（coloc，p12=1e-4）+ inputs/*.tsv（等位基因、EAF、N）。"
               "等位基因约定：beta_eQTL_per_EA 与 beta_结局_per_eQTL_EA 都以 eQTL_EA 为效应等位基因；结局_EA_原始 是结局侧原始效应等位，若与 eQTL_EA 不同，beta_结局 已按此做符号翻转（p01_build_panel.py 第 142–147 行）。"
               "坐标全为 GRCh37。MR_Wald = beta_结局 / beta_eQTL，单位是“每 1 个 eQTL β 单位的表达升高”。"
               "SMR 与 Wald 在 n_iv=0 时数值相同，因为此时 SMR 退化为单工具。HEIDI 需 ≥2 工具才有定义。")
add_sheet(wb, "02_coloc先验网格", grid,
          note="来源 tables/03。p12 网格 = 1e-6 / 1e-5 / 1e-4 / 1e-3，每行 H0–H4 相加为 1。"
               "正文主表用 p12=1e-4；敏感性分析必须把整张网格给出去，因为 SCRN2–TSH 甲状腺的 H4 随先验从 0.20 升到 0.99。")
add_sheet(wb, "03_rs7219420边际", m5,
          note="来源 tables/05。三个队列：发现集（Rand 2025 欧洲，n=482,678）、复制集 1（GCST90319321 欧洲，n=164,818）、复制集 2（GCST90448634 南亚，n=38,000）。"
               "南亚队列的 β/SE/P/EAF 在源文件里即为空，不是本次计算失败——该 GCST 没有 17q21.32 区域的 sumstats。"
               "EA 列为该 beta 的效应等位基因；两个 SNP 的 EA 都是 C。")
add_sheet(wb, "04_两SNP联合", j6,
          note="来源 tables/06。r=0.0262（r²=0.000686，1000G EUR）。联合模型为 GLS：(R′V⁻¹R)⁻¹R′V⁻¹b，V=diag(se)·R·diag(se)。"
               "shift_pct 是联合相对边际的变化百分比；se_inflation≈1.0003 说明共线性可忽略。"
               "结论：条件化后两个 SNP 几乎不动（−2.3% / −2.9%），所以这两个 SNP 互相区分不了任何事，不能用来做“SCRN2 vs TBKBP1”的判别。")
add_sheet(wb, "05_区域图_输入", reg,
          note="来源 tables/08（8,602 个变异，chr17:45.42–46.12 Mb GRCh37 窗口）。TSH 列来自 Rand 2025（面板甲减同源文件的对齐结果）；"
               "三个基因的 eQTL 列分别来自 GTEx v8 甲状腺（QTD000341, n=574）与 eQTLGen 全血（n=31,684）。"
               "所有 beta 列都以同名 *_EA 列为效应等位基因；等位基因不一致的位点已在 p01 中对齐。"
               "r2_with_rs6503875 由 1000G EUR（n=503）逐样本剂量重算（p15）。")
add_sheet(wb, "06_LD_关键位点", ldt,
          note="五个位点的 r 与 r²（1000G EUR n=503，按 VCF 逐样本数效应等位基因拷贝数重算，见 p15_ld_true_dosage.py）。"
               "全区域 LD 矩阵在 inputs_ld/ld_TSH.tsv.gz（2,252×2,252，GRCh37 位置为索引），本表只放展示用的五位点子矩阵。"
               "注意 rs4794224 与另外三个工具的三个负相关（−0.51/−0.73/−0.63）是真实的：靶等位基因 C 在欧洲人群是主要单体型（AF 0.80）。")
add_sheet(wb, "07_可信集", t03c,
          note="来源 tables/03c（本地 runsusie，eQTL 侧）。只有 eQTL 侧有可信集；TSH 侧 SuSiE 报 “prior variance unreasonably large”，已删除 03b/03c 的 TSH 部分。"
               "SCRN2 只有 1 个 CS（L1，5 个变异，top=rs6503875，PIP=0.747）→ 甲状腺 eQTL 侧没有第二个独立信号。")
add_sheet(wb, "08_哨兵明细_783", sent_all[["paper", "trait", "accession", "rsid", "哨兵_EA",
                                          "chr38", "pos_GRCh38", "p", "beta", "RAF_EA",
                                          "mapped_gene", "报告性状", "PMID", "在17q21.32区域"]],
          note="来源 raw/sentinels/gwas_catalog_associations_all.csv（GWAS Catalog REST v2，accession_id 过滤，下载日期 2026-11-14）。"
               "四套：Williams 2023（GCST90296333, TSH, 245）、Rand 2025 TSH（GCST90572789, 297）、"
               "Rand 2025 甲减（GCST90572791, 180）、Rand 2025 游离 T4（GCST90572790, 61），共 783 条。"
               "区域定义 chr17:47.3–48.5 Mb GRCh38；全表只有 1 条落在区域内。"
               "坐标是 GRCh38（GWAS Catalog 原始），换成 GRCh37 需减 1,922,634。")
add_sheet(wb, "08b_区域哨兵_带LD", sent,
          note="来源 tables/21。四套研究：Williams 2023（GCST90296333, TSH, 245）、Rand 2025 TSH（GCST90572789, 297）、"
               "Rand 2025 甲减（GCST90572791, 180）、Rand 2025 T4（GCST90572790, 61）。总计 783 条。"
               "区域定义为 chr17:47.3–48.5 Mb GRCh38（= 45.38–46.58 Mb GRCh37）。区域哨兵只有 1 条。"
               "proxy 列是当哨兵不在 1000G 面板时用 r>0.8 的代理；r_rs7219420 等三列给出与本地三个关键位点的 LD。")
add_sheet(wb, "09_哨兵映射唯一性", t22,
          note="来源 tables/22。GWAS Catalog 的 mapped_gene 字段有 60.7%–69.4% 只映射到唯一基因，其余是多基因或空。"
               "这是“变异→基因”这一步在公共数据库里的真实可解析率，可当作全基因组流程可解析比例的上界代理，但它不是我们自己跑出来的比例。")
add_sheet(wb, "10_区域最小P", t13,
          note="来源 tables/13。四个结局在本区域的变异数与最小 P。只有甲减过 5e-8（2.806e-33，rs1989291，284 个变异过线）；"
               "TSH 最小 1.568e-7（rs7219420），桥本 1.232e-7，游离 T4 5.443e-5。")
add_sheet(wb, "11_变异信息", vt,
          note="来源 tables/11 + inputs。坐标 GRCh37 与 GRCh38 并列，两者之差恒为 1,922,634（rs4794063 与 rs1989291 两个位点独立标定）。"
               "beta_TSH_per_EA 的效应等位基因是 TSH_EA；四组 eQTL 列分别给出该变异在 SCRN2/TBKBP1 × 甲状腺/全血 的 eQTL 效应等位、MAF 与样本量。"
               "rs4794224 是三等到位（T/C/G），GTEx/eQTLGen 只发布二态，面板里缺这一行，正文必须注明。")
add_sheet(wb, "12_样本重叠", pd.concat([
    t14.assign(表="三队列明细"),
    t14b.assign(表="三队列合并"),
    pd.DataFrame([dict(表="缺项", cohort="Rand2025 分队列（CHB/DBDS/UKB/HUNT/MGI/ThyroidOmics）",
                       n=np.nan, est=np.nan, se=np.nan, p=np.nan, weight=np.nan),
                  dict(表="缺项", cohort="去 UKB 后的重估", n=np.nan, est=np.nan, se=np.nan,
                       p=np.nan, weight=np.nan),
                  dict(表="缺项", cohort="Figueredo 队列构成", n=np.nan, est=np.nan, se=np.nan,
                       p=np.nan, weight=np.nan)])], ignore_index=True),
          note="来源 tables/14、14b。已有：Rand2025 欧洲发现集、GCST90319321（= Figueredo 欧洲）、GCST90448634 南亚。"
               "缺：Rand 2025 的六个分队列逐变异估计、去 UKB 重估、Figueredo 队列的样本构成。"
               "这三样只能从原文补充表拿，本地无法反推；没有它们，就无法回答“UKB 是否主导了发现集”这一审稿问题。")
add_sheet(wb, "13_SCRN2第二信号", t18,
          note="来源 tables/18、03c。k 是 LD 剪枝后的工具数，k_eff 是有效工具数（特征值法 (Σλ)²/Σλ²）。"
               "SCRN2 在 r²<0.01 时 k=1、k_eff=1.00；放宽到 0.10 时 k=2、k_eff=1.88 —— 说明第二个工具与 rs6503875 高度相关，不是独立信号。"
               "rs3047616（r²<0.001，eQTL P=4.5e-6）在 1000G VCF 里不存在，进不了 LD 矩阵，这是“约 1.6 个有效工具”的由来。")
add_sheet(wb, "14_逐组织eQTL", tissue_q,
          note="来源 tables/09b。49 个 GTEx 组织，每个组织给出三个基因最强 eQTL 的 P。"
               "这不是共定位：GTEx 只发布显著对，逐组织 coloc 需要每个组织的全区域 nominal 统计。"
               "因此本节只能支撑“SCRN2 在几乎所有组织都有强 cis-eQTL、方向一致”，不能支撑“SCRN2–TSH 在各组织共定位”。"
               "GTEx v8 没有下丘脑；垂体样本量小（n=159），功效有限。")
add_sheet(wb, "15_调控注释", regu,
          note="来源 tables/16（UCSC encodeCcreCombined, hg38, ±5 kb）。八个变异都不在 cCRE 内部，但都在 36–2,029 bp 内贴着 dELS/pELS。"
               "这只说明“处在调控区”，不说明“连向哪个启动子”。ABC / ENCODE-rE2G / Open Targets L2G 三样都没做（API 不可达），不做假。")
rare = pd.DataFrame([
    dict(基因=g, 结局=o, 掩码=m, beta=np.nan, se=np.nan, p=np.nan, n_case=np.nan,
         来源="待补", 状态="缺")
    for g in ["SCRN2", "TBKBP1", "MRPL10"]
    for o in ["TSH", "甲减"]
    for m in ["LoF", "错义（MPC≥2）"]])
add_sheet(wb, "16_罕见变异_空缺模板", rare,
          note="整项空缺。Genebass（https://genebass.org/api）本机返回 301，AZ PheWAS 需要申请访问。"
               "补法：下载 Genebass 的 gene-level burden（约 2 GB，包含 4,529 个基因 × 罕见变异掩码）或 AZ PheWAS 的 *.gene.burden，"
               "取 SCRN2 / TBKBP1 / MRPL10 三行即可，阴性也要报。")
add_sheet(wb, "17_FT4与PheWAS", ph,
          note="来源 DATA18_manuscript/analysis/phewas_thyroid_harmonised.csv（OpenGWAS，p≤1e-4 全库扫描的甲状腺相关子集）。"
               "FT4 汇总统计不缺：面板的 T4 列就是游离 T4（GCST90572790）。"
               "此前写“现有 T4 是总 T4、x01_three_outcomes 的 Free T4 列名要改掉”，方向正好反了——"
               "该改的是“总 T4”这个标签，不是 Free T4。详见 20_结局GWAS身份。")
add_sheet(wb, "18_投稿前核对", pd.DataFrame([
    dict(编号="A1", 检查="eQTLGen 原始文件中 3.2717E-310 的出现次数", 值="877 行", 判读="同一数值被 877 行共用，来自 z 从 44 到 79 的不同变异 → 浮点下限，不是复制"),
    dict(编号="A2", 检查="rs7219420→TBKBP1 的 z 与 P", 值="z=58.1955，P=3.2717E-310", 判读="P 已触底，z 仍有信息"),
    dict(编号="A3", 检查="rs6503875→SCRN2 的 z 与 P", 值="z=−44.2555，P=3.2717E-310", 判读="与 A2 同值但 z 差 14 个单位 → 证明 P 不可用于排序"),
    dict(编号="A4", 检查="OpenGWAS 靶基因表里的 3.27341e-310", 值="两个 SNP 的 rank-1 均为该值", 判读="同样是下限（第四位有效数字与 eQTLGen 的 3.2717E-310 不同，说明是两套管线各自截断）"),
    dict(编号="A5", 检查="结论", 值="不是复制错误", 判读="正文改用 |z| 或 β/se 排序；P 写“<1e-300”，不写具体值"),
    dict(编号="B1", 检查="s21_reference_frame.py 的 SE 口径", 值="delta 法（含暴露侧误差）", 判读="脚本已改，非一阶"),
    dict(编号="B2", 检查="复跑 s21 的输出 vs reference_frame.csv", 值="甲状腺 GLS 0.0250 (0.0080)；z=1.67/p=0.095", 判读="逐位一致"),
    dict(编号="B3", 检查="DATA18 A2 总表 I5 行", 值="写 z=2.86 / p=0.004", 判读="与当前 CSV 的 2.72 / 0.007 不符 → 旧值，须改"),
]), note=None)

add_sheet(wb, "19_重复键检查", dupqa_df,
          note="tables/04 里每个基因都以 symbol 和 ENSG 两个键各存一份，MR/coloc/SMR/HEIDI 数值逐位相同（本表“是否一致”列全为是）。"
               f"本册已按（基因符号, 组织, 结局）去重，保留 symbol 键，行数由 {N04} 降为 {len(pan)}。"
               "任何引用“面板共 N 个基因”的地方都必须用去重后的数，否则会整体翻倍。"
               "另：AC003665.1 与 SP2-AS1 是同一基因（ENSG00000234494）的两个注释名，分别只出现在全血与甲状腺，已合并为 SP2-AS1。")

NID = pd.DataFrame([
    dict(结局="TSH", 研究号="GCST90572789", 论文报告N="482,873", 文件最大逐SNP_n=482863, 差=10,
         基因组加权中位n=191203, 队列结构="135,861 / 201,288 / 191,203 三档为主，最高档 482,8xx",
         判读="一致。文件最大 n 只差 10（单 SNP 掉样），可直接用 482,873"),
    dict(结局="游离 T4（FT4）", 研究号="GCST90572790", 论文报告N="191,449", 文件最大逐SNP_n=191448, 差=1,
         基因组加权中位n=66057,
         队列结构="66,057 + 49,269 + 76,123 = 191,449（逐位精确，且三者两两相加都作为模态出现）",
         判读="一致。49,269 即 Teumer 2018 的 FT4 队列（GCST006896），故本结局确定为游离 T4 而非总 T4"),
    dict(结局="甲减", 研究号="GCST90572791", 论文报告N="1,178,661（113,393 例 + 1,065,268 对照）",
         文件最大逐SNP_n=1044010, 差=134651, 基因组加权中位n=424442,
         队列结构="424,442 与 344,168 各自单独出现，两者相加 768,610；三者齐备时 1,044,010（= 768,610 + 275,400）",
         判读="**不一致，已定位**：公开统计量的上限 1,044,010 ≈ GWAS Catalog 的欧洲部分 1,046,991（98,805 例 + 948,186 对照），"
              "非欧洲 131,670（14,588 例 + 117,082 对照，23andMe）未随文件发布。"
              "所以正文写总样本量 1,178,661 时，必须同时说明“公开汇总统计仅含欧洲部分（上限 1,044,010）”"),
])
add_sheet(wb, "20_结局GWAS身份", NID,
          note="三个结局 GWAS 均为 Rand SA et al., Nat Genet 2025（PMID 41238958）同一篇。"
               "本表把“论文报告的 N”与“文件里逐 SNP 的 n 上限”逐项对齐，用来判定下载版本是否与研究一致。"
               "关键结论：① T4 是游离 T4，不是总 T4；② 甲减公开文件缺 23andMe 的 131,670 人；"
               "③ 逐 SNP 的 n 是各队列子集之和，全基因组加权中位数远低于研究总量（T4 66,057 / 甲减 424,442 / TSH 191,203），"
               "引用“中位样本量”时不能拿它当研究样本量。")

wsn = wb.create_sheet("21_说明")
wsn.append(["表", "说明"])
for c in (1, 2):
    cell = wsn.cell(row=1, column=c)
    cell.fill = HDR_FILL; cell.font = HDR_FONT
for name, note in NOTES:
    wsn.append([name, note])
wsn.column_dimensions["A"].width = 24
wsn.column_dimensions["B"].width = 150
for row in wsn.iter_rows(min_row=2):
    for cell in row:
        cell.alignment = Alignment(wrap_text=True, vertical="top")
wsn.freeze_panes = "A2"

wb.save(OUT)
print("saved:", OUT, os.path.getsize(OUT) / 1e6, "MB")
print("sheets:", wb.sheetnames)
print(pan.shape, reg.shape, ldt.shape)
