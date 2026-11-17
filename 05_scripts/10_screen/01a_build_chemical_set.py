# -*- coding: utf-8 -*-
"""
Step 01a: build the EDC chemical list and the thyroid-expressed gene set
Output:
  data/edc_chemicals.csv       EDC chemicals (MeSH D052244 descendants plus pattern matching)
  data/thyroid_expressed.txt   genes with GTEx thyroid TPM>=1
"""
import gzip, re, csv, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, DATA = os.path.join(BASE, "raw"), os.path.join(BASE, "data")
os.makedirs(DATA, exist_ok=True)

# ---------------------------------------------------------------- GTEx thyroid expression
def load_thyroid_genes(path, tpm_cut=1.0):
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        f.readline(); f.readline()          # the #1.2 line and the dimension line
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
print(f"[GTEx] genes with thyroid TPM>=1: {len(thyroid)}")

# ---------------------------------------------------------------- EDC chemical list
# Sources: (1) chemicals below the MeSH heading Endocrine Disruptors (D052244)
#          (2) name patterns for the main classes of environmental endocrine disruptors
#              (PFAS / phthalates / bisphenols / PCBs / PBDEs / organochlorine pesticides /
#               parabens / triclosan / dioxins / metals)
MESH_EDC = "MESH:D052244"

# Operational definition: endocrine disruptors are defined by chemical class, following the six
# main classes of the TEDX scheme, and extended with polycyclic aromatic hydrocarbons,
# alkylphenols, benzophenones, perchlorate and four classical interfering metals.
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
print(f"[EDC] chemicals in total: {len(chemicals)} (MeSH descendants {n_mesh}, added by pattern {len(chemicals)-n_mesh})")
