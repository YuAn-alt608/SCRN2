# -*- coding: utf-8 -*-
"""
Step 01a: 构建 EDC 化学物清单 与 甲状腺表达基因集
输出:
  data/edc_chemicals.csv       EDC 化学物（MeSH D052244 子类 + 六大类模式匹配）
  data/thyroid_expressed.txt   GTEx 甲状腺 TPM>=1 的基因符号
"""
import gzip, re, csv, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, DATA = os.path.join(BASE, "raw"), os.path.join(BASE, "data")
os.makedirs(DATA, exist_ok=True)

# ---------------------------------------------------------------- GTEx 甲状腺表达
def load_thyroid_genes(path, tpm_cut=1.0):
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        f.readline(); f.readline()          # #1.2 和维度行
        header = f.readline().rstrip("\n").split("\t")
        idx = header.index("Thyroid")
        genes = []
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= idx:
                continue
            try:
                v = float(p[idx])
            except ValueError:
                continue
            if v >= tpm_cut:
                genes.append(p[1])          # Description = gene symbol
    return set(genes)

thyroid = load_thyroid_genes(os.path.join(RAW, "GTEx_median_tpm.gct.gz"))
with open(os.path.join(DATA, "thyroid_expressed.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(sorted(thyroid)))
print(f"[GTEx] 甲状腺 TPM>=1 基因数: {len(thyroid)}")

# ---------------------------------------------------------------- EDC 化学物清单
# 依据: (1) MeSH 主题词 Endocrine Disruptors (D052244) 的下位化学物
#       (2) 六大类环境内分泌干扰物的名称模式（PFAS / 邻苯二甲酸酯 / 双酚 /
#           PCB / PBDE / 有机氯农药 / 对羟基苯甲酸酯 / 三氯生 / 二噁英 / 重金属）
MESH_EDC = "MESH:D052244"

# 环境内分泌干扰物的操作性定义：按化学类别界定（与 TEDX 分类体系一致的六大类，
# 并补充多环芳烃、烷基酚、二苯甲酮、高氯酸盐与四种经典干扰金属）
CLASS_PATTERNS = [
    r"perfluor|polyfluor|\bPFOA\b|\bPFOS\b|\bPFNA\b|\bPFHxS\b|\bPFDA\b|\bPFBA\b|\bPFHpA\b|\bPFDoA\b",
    r"phthalat",
    r"bisphenol",
    r"polychlorinated biphenyl|\bPCB\s?\d|\bPCB\b",
    r"polybrominated diphenyl ether|\bPBDE\b|\btetrabromobisphenol\b",
    r"\bDDT\b|dichlorodiphenyltrichloroethane|\bDDE\b|\bDDD\b",
    r"hexachlorocyclohexane|lindane|chlordane|heptachlor|\baldrin\b|dieldrin|endosulfan|toxaphene|mirex|hexachlorobenzene",
    r"paraben",
    r"triclosan|triclocarban",
    r"\bdioxin\b|dibenzodioxin|dibenzofuran|tetrachlorodibenzodioxin|\bTCDD\b",
    r"\bcadmium\b|\blead\b|\bmercury\b|\bmethylmercury\b|\barsenic\b|arsenite|arsenate",
    r"perchlorate|thiocyanate",
    r"octylphenol|nonylphenol|alkylphenol|benzophenone",
    r"glyphosate|atrazine|chlorpyrifos|malathion|permethrin|cypermethrin|imidacloprid|fipronil|diazinon|parathion",
    r"benzo\[?a\]?pyrene|fluoranthene|polycyclic aromatic",
]

pat = re.compile("|".join(CLASS_PATTERNS), re.I)
chemicals = {}
with gzip.open(os.path.join(RAW, "CTD_chemicals.tsv.gz"), "rt",
               encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 6:
            continue
        name, cid = c[0], c[1]
        parents = c[8] if len(c) > 8 else ""
        synonyms = "|".join(c[11:13]) if len(c) > 11 else ""
        is_mesh = MESH_EDC in parents
        is_class = bool(pat.search(name)) or bool(pat.search(synonyms))
        if is_mesh or is_class:
            chemicals[cid] = {
                "chemical_id": cid,
                "chemical_name": name,
                "mesh_edc": int(is_mesh),
                "class_matched": int(is_class),
            }

with open(os.path.join(DATA, "edc_chemicals.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["chemical_id", "chemical_name", "mesh_edc", "class_matched"])
    w.writeheader()
    for v in chemicals.values():
        w.writerow(v)

n_mesh = sum(v["mesh_edc"] for v in chemicals.values())
print(f"[EDC] 化学物总数: {len(chemicals)} (MeSH 子类 {n_mesh}, 模式补充 {len(chemicals)-n_mesh})")
