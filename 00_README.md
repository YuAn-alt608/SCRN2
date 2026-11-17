# DATA1 · 数据、脚本与证据归档

**建立日期**：2026-11-14
**对应稿件**：`03_稿件/MR1_manuscript_v5_2026-11-14.docx`
**叙事来源**：`01_故事线_中文_总纲.md`
**数字校验**：`05_脚本/verify_numbers.py` → `05_脚本/核对报告.txt`

本目录只放**支撑故事线所需的最小充分集**：每一个数字都有对应的结果表，每一张结果表都有生成它的脚本。
未复制的原始大文件（eQTL 区域汇总、GWAS 全量、MD 轨迹）在本文件第 4 节列明位置与再生方式。

---

## 1. 目录结构

```
DATA1/
├─ 00_README.md                  ← 本文件
├─ 01_故事线_中文_总纲.md         ← 写作思路与七层证据链（唯一叙事来源）
├─ 02_证据对照表.md               ← 每一条主张 → 数据文件 → 行号/字段 → 校验输出
├─ 03_稿件/
│   ├─ MR1_manuscript_v5_2026-11-14.docx
│   └─ 稿件结构说明.md
├─ 04_数据/                       ← 55 个结果表，5.5 MB
├─ 05_脚本/                       ← 70 个分析脚本 + 2 个自写校验脚本 + 核对报告.txt
└─ 06_图表规格.md                 ← 9 图 6 表的逐面板规格（中文，含数据源与已核对值）
```

---

## 2. 数据清单（按故事线的层组织）

### 2.1 `01_筛查/` — 第 1 层：名单从哪来

| 文件 | 行数 | 内容 | 用在哪 |
|---|---|---|---|
| `edc_chemicals.csv` | 1,580 | EDC 化学物（MeSH ID、名称、是否 MeSH 标注、类别匹配） | 4,526 基因的上游 |
| `gene_degree_all.csv` | 24,264 | 每个基因被多少种 EDC 调控（`n_edc_chemicals`）、互作数、是否甲状腺表达 | 第 1 层的"剂量无关"分析 |
| `moduleA_stats.txt` | 8 | 构建链条的关键计数：1,580 / 18,306 / 256,899 / 24,264 / 5,272 / 4,526 | 表 1 |
| `edc_thyroid_genes.txt` | 4,526 | 最终候选基因集（每行一个 symbol） | 全篇入口 |
| `ivs.tsv` | 3,758 | 每个基因的 cis 工具（SNP、β、SE、F、r²） | 表 1、方法 |
| `iv_summary.csv` | 3,559 | 每基因的工具数与最强工具 | 方法 |
| `MR_all_results.csv` | 12,685 | 全部 MR 估计；**唯一 gene–outcome 对 12,636**（49 行为第二方法） | 表 2、图 2 |
| `evidence_tiers.csv` | 733 | FDR<0.05 的 733 对 + A/B/C/D 分级 + 共定位字段 | 表 2、图 2 |
| `MR_significant_FDR05.csv` | 733 | 同上，精简列 | 备用 |

> 一致性：`evidence_tiers` 的 733 对与 `MR_all_results` 中 `fdr<0.05` 的 733 对**逐对相同**（校验脚本 §1 输出 `True`）。

### 2.2 `02_共定位/` — 第 4 层

| 文件 | 行数 | 内容 |
|---|---|---|
| `coloc_results.csv` | 2,464 | 每对基因–结局的 PP.H0–H4；733 个显著对 + 1,731 个对照对 |
| `h4_circularity.csv` | 71 | PP.H4>0.8 的 71 对，逐对标注"峰 SNP 是否就是工具"等三个循环性指标 |
| `h4_independent_loci.csv` | 10 | 把窗口缩到各自独立位点后仍 H4>0.8 的 10 个 |
| `coloc_p12_summary.txt` | — | 先验敏感性：p12 = 1e-4/1e-5/1e-6 时 H4>0.8 为 196/71/18，排序不变 |

### 2.3 `03_复现/` — 第 2 层

| 文件 | 行数 | 内容 |
|---|---|---|
| `finngen_replication.csv` | 733 | 每条命中在 FinnGen R12 的 β、P、是否复现 |
| `finngen_summary.txt` | 6 | 分结局复现率：甲减 94.3% 同向 / 66.9% 复现 |
| `finngen_null_control.txt` | — | 阴性对照（FDR>0.5 抽 388 条）6.7%；对照组也用它自己的工具做 IVW |

### 2.4 `04_组织/` — 第 3 层

| 文件 | 行数 | 内容 |
|---|---|---|
| `thyroid_3source_validation.csv` | 70 | 70 个工具在三个甲状腺来源的 P 与方向 |
| `thyroid_3source_summary.txt` | 7 | 可比 23/62/13；显著 13/27/13；同向 6/45/9 |
| `gtex_tissue_rates.csv` | 49 | 同一批 70 个工具在 49 个 GTEx 组织的检出率（全血 41.4%、甲状腺 31.4%、中位 12.9%） |
| `gtex_tissue_summary.txt` | 23 | 上述的文本汇总与阈值说明 |

### 2.5 `05_位点/` — 第 6 层

| 文件 | 行数 | 内容 |
|---|---|---|
| `23_信号解析归属总表.csv` | 8 | **核心**：3 个信号 × 组织 × 候选基因，含 eQTL、cis-MR、Q/HEIDI、窗口内 coloc、裁决 |
| `31_联合条件模型.csv` | 84 | GLS 联合模型，五种模型 × 四个结局；条件化位移的唯一来源 |
| `24_按信号窗口的coloc.csv` | 96 | 把窗口缩到信号内部后的 coloc（SCRN2 对 TSH 0.988、对甲减 H3=1.00） |
| `04_cisMR_SMR_HEIDI.csv` | 587 | 面板内各基因的 cis-MR/SMR/HEIDI |
| `02b_每基因最强eQTL.csv` | 120 | 每基因的 lead eQTL 与位置 |
| `11_变异信息表.csv` | 8 | 七条信号的 rsID、GRCh37/b38 位置、等位基因 |
| `12_GWAS目录_位点身份.csv` | 6 | GWAS Catalog 在该位点给出的基因映射（"TBKBP1;TBX21"这类） |
| `21_已发表哨兵_17q21.csv` | 1 | 已发表哨兵的落点 |
| `00_LD面板位点.csv` | 2,279 | 1000G EUR（n=503）面板内的位点与 MAF |

### 2.6 `06_表达与单细胞/` — 第 9 层与暴露端

| 文件 | 行数 | 内容 |
|---|---|---|
| `genes_expression_profile.csv` | 10 | 54 组织 TPM：SCRN2 甲状腺 80.5 / 全血 8.2 / 甲状腺排名 1 / τ=0.594 |
| `a2_summary.csv` | 5 | 50 组织 eQTL：50/50 显著、49 同向 |
| `SCRN2_thyroid_celltype.csv` | 22 | 单细胞逐类型：滤泡细胞 2.96% 检出（最高） |
| `s07_panel_pc10.csv` | 19 | 共表达：NKX2-1 ρ=0.36、PAX8 0.19（校正 10 PC） |
| `SCRN2_CTD_EDC_human.csv` | 23 | **PFNA 选择的直接依据**：23 条人类 EDC 记录，PFAS 3 条同源研究 |
| `SCRN2_CTD_EDC_summary.txt` | — | 上述的解读与方向性矛盾的原始记录 |
| `SCRN2_geo_expression_validation.txt` | — | 两个 GEO bulk 数据集的阴性结果 + 功效界限（≥0.24 log2FC） |

### 2.7 `07_结构/` — 第 10 层

| 文件 | 内容 |
|---|---|
| `dock_T3.csv` / `dock_PFNA.csv` | 盲对接前 10 位姿：T3 −10.27、PFNA −5.78 kcal/mol |
| `summary_table.csv` | 六项动力学 + MM-GBSA 分量（T3 vs PFNA 并排） |
| `metrics_T3.json` / `metrics_PFNA.json` | 逐帧原始指标 |
| `mmgbsa_T3.json` / `mmgbsa_PFNA.json` | MM-GBSA 逐帧分解 |
| `motif_occlusion.json` | 基序遮挡分析原始结果 |
| `pockets_ligsite.csv` | 空腔普查：最大内部封闭空腔 80.4 Å³（T3 需 ≈330 Å³） |
| `timeseries_PFNA.csv` | 1,001 帧时间序列 |
| `MD接触重算.csv` | **本目录自算**：三个基序的最小距离与各截断下帧比例 |
| `MD逐残基接触.csv` | **本目录自算**：21 个残基的接触频率 |

### 2.8 `08_平面与跨族裔/` — 第 7–8 层

| 文件 | 内容 |
|---|---|
| `gw_thyroid_mr.csv` | 甲状腺 eQTL 工具（GTEx v8，n=574）对四结局的 3,697 个估计 |
| `s08_plane_genes.csv` | 915 个双结局基因的 b_TSH、b_T4、比率 r 与置信区间 |
| `s08_plane.txt` / `s08_ratio_summary.csv` | 平面分析的文本输出与参照值（TPO −1.07、SCRN2 +0.15） |
| `s08b_power.txt` | 功效匹配后的平面复核 |
| `three_cohort_meta.csv` | 三队列：Rand2025 EUR 0.0336 / GCST90319321 EUR 0.0404 / GCST90448634 SAS 0.0288 |
| `x02_transancestry.txt` | 固定效应合并（0.0330，95%CI 0.0167–0.0493）、Q=0.034、P=0.853、I²=0% |

---

## 3. 脚本清单

### 3.1 自写校验脚本（本目录新增，可直接运行）

| 脚本 | 作用 | 依赖 |
|---|---|---|
| `verify_numbers.py` | 从结果表重算故事线与稿件里的**每一个数字**，输出 `核对报告.txt` | pandas / numpy / scipy |
| `recompute_md_contacts.py` | 从 MD 轨迹重算三个基序的最小距离与逐残基接触频率 | MDAnalysis |

### 3.2 原始分析脚本（按流程编号，均已复制）

| 子目录 | 脚本 | 对应步骤 |
|---|---|---|
| `10_筛查/` | `01a_build_chemical_set.py` … `38_scrn2_ctd_edc.py`（18 个） | EDC 集 → 基因集 → 工具 → MR → 共定位 → 分级 |
| `20_共定位/` | `09_coloc.R`、`12a_coloc_triage.py`、`12c_coloc_susie.R`、`21_coloc_p12_sensitivity.R`、`22_h4_independent.py` | abf / susie / 先验敏感性 / 缩窗 |
| `30_复现/` | `15_fetch_finngen.py`、`16_finngen_replication.py` | FinnGen R12 |
| `40_组织/` | `17_…`、`24_/26_gtex_tissue_control.py`、`28–30_thyroid_3source_*` | 组织可迁移性 |
| `50_位点/` | `p01`–`p25`（16 个） | 面板构建、GLS 联合模型、逐信号窗口 coloc、归属裁决 |
| `60_表达/` | `31_gtex_expression_profile.py`、`32_scrn2_celltype.py`、`37_scrn2_geo_expression.py`、`38_scrn2_ctd_edc.py`、`s07_coexpression_hardened.py` | 组织表达、单细胞、GEO、CTD、共表达 |
| `70_结构/` | `01_prep_ligands.py` … `17_hspa8_insertion_test.py`（10 个） | 参数化 → 空腔 → 盲对接 → MD → 指标 → MM-GBSA → 遮挡 |
| `80_图表/` | `make_venn.py`、`s08_tsh_t4_plane.py`、`s08b_power_matched_plane.py` | 图 5、图 7 |

> **注意**：原始脚本里的路径多为绝对路径（`D:/NewBeginning/MR1/...`）。迁移到别的机器时需要改 `BASE`。
> 本目录的两个自写脚本已改为相对路径（按 `__file__` 上溯两级），可直接运行。

---

## 4. 未复制的原始大文件（及再生方式）

| 未复制的内容 | 原位置 | 体量 | 为什么没复制 | 怎么再生 |
|---|---|---|---|---|
| eQTL 区域汇总 | `DATA10/data/eqtl_regions_all.txt.gz`、`region_*.tsv.gz` | 195 MB | 可由 `08_extract_outcome_regions.py` 从 eQTLGen 与四个 GWAS 重建 | 跑 `08_extract_outcome_regions.py` |
| GWAS 全量 | `DATA10/raw/{TSH,T4,HYPO}.tsv.gz`、`HASH_full.tsv.gz` | ~1 GB | 公开数据，按第 5 节链接重新下载 | `02b_extract_outcomes.py`、`02d_fetch_opengwas.py` |
| eQTLGen 全量 | `DATA10/raw/eQTLGen_full.txt.gz` | 大 | 同上 | eQTLGen 官网 |
| MD 轨迹 | `DATA14/md_runs/PFNA/md.xtc`（1,001 帧） | 数百 MB | 可由 `05_build_systems.py` + `06_run_md.py` 重跑（约 10 ns × 2 体系） | 见 `70_结构/` |
| GTEx/单细胞中间件 | `DATA17_repro_checks/raw/` | 1.5 GB | 外部下载 | 见原目录 `scripts/` |

---

## 5. 外部数据源（类型、版本、规模）

| 数据 | 版本 / 标识 | 规模 | 用途 |
|---|---|---|---|
| CTD | Chemical–Gene Interactions，人源 | 256,899 条互作 / 1,580 EDC | 基因集与 EDC 度 |
| eQTLGen | cis-eQTL，全血 | n = 31,684 | 主工具 |
| GTEx | v8，甲状腺 | n = 574 | 独立工具集、组织基准 |
| eQTL Catalogue | QTD000341（=GTEx v8 甲状腺） | — | 三来源之一 |
| Zenodo 23068645 | 甲状腺 eQTL（salmon 定量） | — | 三来源之一 |
| GTEx 官方 API | v8/v10 逐组织 | 49–50 组织 | 跨组织 |
| TSH / T4 / 甲减 | 甲状腺功能 GWAS（见下方待确认项） | TSH n≈482,800；T4 n≈191,300；甲减 n≈1.04 M | 结局 |
| 桥本甲状腺炎 | `ebi-a-GCST90018855` | n = 395,640 | 结局 |
| FinnGen | R12，`E4_HYTHY_AI_STRICT` | 54,752 例 / 362,639 对照 | 复现 |
| 南亚队列 | `GCST90448634` | n = 38,000 | 跨族裔 |
| AlphaFold DB | Q96FV2（SCRN2），425 aa，pLDDT 94.7 | — | 结构 |

### ✅ 两个待确认项已定案（2026-11-14 11:3x，`DATA21/06_结局GWAS身份核定`）

三项结局全部来自同一篇：**Rand SA et al., _Nat Genet_ 2025;57:3007–3015（PMID 41238958）**，同一批汇总统计。

**1. T4 = 游离 T4（FT4），不是总 T4。** 三条独立证据：

- 摘要原文：*"…hypothyroidism (113,393 cases and 1,065,268 controls), **free thyroxine (191,449 individuals)** and thyroid-stimulating hormone (482,873 individuals)"*；
- GWAS Catalog `GCST90572790` 的 `initial_sample_size` = **191,449**，而 `T4.tsv.gz` 的逐 SNP n 上限是 **191,448**（差 1，单样本掉样），两者是同一批人；
- 文件里的 n 是**队列子集之和**，且可精确分解：**66,057 + 49,269 + 76,123 = 191,449**。其中 **49,269 正是 Teumer 2018 游离 T4 队列**（`GCST006896`，trait "Free thyroxine concentration"，49,269 欧洲人）。总 T4 的 meta 不可能装进一个 FT4 队列。

→ 内部文档里所有"总 T4"都是错的，已改。`x01_three_outcomes.csv` 的 "Free T4" 列名**本来就是对的**，当初要改的方向反了。
→ 生化含义也变了：TSH↑ + FT4 正常 = **亚临床甲减**的定义，比"总 T4 不动"贴切，可进摘要。

**2. 甲减 GWAS 出处已确认，n 的差异是结构性的，不是下载错版本。**

| 口径 | 数值 | 来源 |
|---|---|---|
| 论文总量 | **1,178,661** = 113,393 例 + 1,065,268 对照 | 摘要；GWAS Catalog 分族裔加总亦为此数 |
| 欧洲部分 | **1,046,991** = 98,805 例 + 948,186 对照 | GWAS Catalog `GCST90572791` |
| 非欧洲部分 | **131,670** = 14,588 例 + 117,082 对照（23andMe） | 同上，未随文件发布 |
| 文件逐 SNP n 上限 | **1,044,010** | `HYPO.tsv.gz` 实测 |

1,044,010 ≈ 1,046,991 − 2,981（0.28%，单 SNP 掉样）；距 1,178,661 差的 **134,651 ≈ 23andMe 的 131,670**。
即**公开汇总统计只含欧洲队列**。正文写 1,178,661 时须加一句"公开统计仅含欧洲部分（上限 1,044,010）"。

另：`HYPO.tsv.gz` 的**逐 SNP n 中位数不是 1,043,640，而是 424,442**（34,556,406 个变异里 40.9% 只有 424,442 这一个队列）。
1,043,640 是"n>1e6 的那 814 万个变异"的中位数（1,043,670），被误当成了全基因组中位数。
同样，`T4.tsv.gz` 的逐 SNP 中位数是 **66,057**，不是 191,300。

---

## 6. 一键复现

```bash
cd D:/NewBeginning/MR1/DATA1/05_脚本

# 1) 所有数字回源（秒级）
python verify_numbers.py
#    输出：核对报告.txt

# 2) 连 MD 接触一起重算（需 MDAnalysis，约 2 分钟）
python verify_numbers.py --md
python recompute_md_contacts.py
#    输出：04_数据/07_结构/MD接触重算.csv、MD逐残基接触.csv

# 3) 重建 6 张表（从 04_数据 实算，不硬编码）+ 9 张图（约 1 分钟）
python 90_图表/make_figures.py
#    输出：03_稿件/figs/fig1-9.png（300 dpi）与 .pdf（矢量文字）

# 4) 重建稿件（图在文末 Figures 节，表就地，字数自动回填标题页）
cd ../03_稿件
python build_v5.py
python export_pdf.py          # 可选：导出 PDF 并把每页切成 PNG 供核对
```

稿件里的每一个百分比，都能在 `核对报告.txt` 里找到对应的分子与分母。**这是本目录的验收标准。**

三条自检口子（都应通过）：

```bash
cd D:/NewBeginning/MR1/DATA1/03_稿件
python -c "from docx import Document; d=Document('MR1_manuscript_v5_2026-11-14.docx'); \
print('占位', len([p for p in d.paragraphs if p.text.startswith('【')]), '真表格', len(d.tables))"
# 期望：占位 0 真表格 6
```

图表状态（2026-11-14）：正文不再保留任何灰底占位框；6 张表全部由 `03_稿件/tables_data.py`
从 `04_数据` 实算后渲染为三线表，9 张图由 `05_脚本/90_图表/make_figures.py` 生成后置于文末
`Figures` 节（每图一页，图注随图）。正文已补 `(Figure n)` / `(Table n)` 引用。
