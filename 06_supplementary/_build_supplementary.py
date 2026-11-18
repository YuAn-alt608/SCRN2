#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the standalone Supplementary Material file.

Three tables, one document, numbered in the order the manuscript cites them:

  S1  five evidence layers: layer sizes and threshold sensitivity (portrait)
  S2  PFNA at the SCRN2 surface — docking, MD, MM-GBSA (portrait)
  S3  two-sample MR of the nine instrumentable EDC exposures (landscape, 32 rows)

S3 is copied element-for-element out of the manuscript so that the typesetting
checked in v15 is preserved; S1 and S2 are re-set from the CSV files that
`verify_numbers.py` reads, in the same booktabs style used for the manuscript
tables. S2 is transposed — quantities as rows, configurations as columns —
because eleven columns do not fit a portrait page and the reader compares
configurations anyway.

The document is built on the v15 file so that the manuscript's stylesheet, page
size and page-number footer carry over unchanged.
"""
import copy
import csv
import os

from docx import Document
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.table import WD_ALIGN_VERTICAL
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt
from docx.table import Table
from docx.text.paragraph import Paragraph

BASE = "03_manuscript/MR1_manuscript_v15_2026-11-18.docx"
OUT = "06_supplementary/Supplementary_Material_SCRN2.docx"

EDC_CSV = "06_supplementary/Table_S3_EDC_exposure_arm_MR.csv"
SENS_CSV = "06_supplementary/Table_S1_sensitivity.csv"

TITLE = ("Colocalisation-resolved expression, thyroid function and a "
         "perfluoroalkyl congener at SCRN2")

# ------------------------------------------------------------------ template
doc = Document(BASE)

src_blocks = []
for ch in doc.element.body.iterchildren():
    if ch.tag == qn("w:p"):
        src_blocks.append(Paragraph(ch, doc))
    elif ch.tag == qn("w:tbl"):
        src_blocks.append(Table(ch, doc))

i_cap3 = next(k for k, b in enumerate(src_blocks)
              if isinstance(b, Paragraph) and b.text.strip().startswith("Supplementary Table S3 |"))
i_note3 = next(k for k, b in enumerate(src_blocks)
               if isinstance(b, Paragraph) and b.text.strip().startswith("A single ligand was simulated"))
s3_caption = copy.deepcopy(src_blocks[i_cap3]._element)
s3_table = copy.deepcopy(src_blocks[i_cap3 + 1]._element)
s3_note = copy.deepcopy(src_blocks[i_note3]._element)

body = doc.element.body
for ch in list(body.iterchildren()):
    if ch.tag in (qn("w:p"), qn("w:tbl")):
        body.remove(ch)

# the manuscript's figure images are now unreferenced; drop the relationships so
# that the package does not carry 7 MB of unused media
for rid, rel in list(doc.part.rels.items()):
    if rel.reltype.endswith("/image"):
        doc.part.drop_rel(rid)


# ------------------------------------------------------------------ helpers
def para(text="", size=12, bold=False, italic=False, align=WD_ALIGN_PARAGRAPH.LEFT,
         before=0, after=6, line=1.5, keep=False, indent=None, hanging=None):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.alignment = align
    pf.space_before = Pt(before)
    pf.space_after = Pt(after)
    pf.line_spacing = line
    pf.keep_with_next = keep
    if indent is not None:
        pf.left_indent = Cm(indent)
    if hanging is not None:
        pf.first_line_indent = Cm(-hanging)
    if text:
        r = p.add_run(text)
        r.font.name = "Times New Roman"
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.italic = italic
        r._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    return p


def cell_text(cell, text, size=7.5, bold=False, align=WD_ALIGN_PARAGRAPH.LEFT,
              before=0.5, after=0.5, line=0.95):
    p = cell.paragraphs[0]
    pf = p.paragraph_format
    pf.alignment = align
    pf.space_before = Pt(before)
    pf.space_after = Pt(after)
    pf.line_spacing = line
    r = p.add_run(text)
    r.font.name = "Times New Roman"
    r.font.size = Pt(size)
    r.font.bold = bold
    r._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER


def table_rules(table, header_row=0):
    tblPr = table._tbl.tblPr
    for el in tblPr.findall(qn("w:tblBorders")):
        tblPr.remove(el)
    borders = OxmlElement("w:tblBorders")
    spec = [("top", "single", 10, "333333"), ("left", "nil", 0, "auto"),
            ("bottom", "single", 10, "333333"), ("right", "nil", 0, "auto"),
            ("insideH", "single", 2, "BFBFBF"), ("insideV", "nil", 0, "auto")]
    for name, val, sz, color in spec:
        el = OxmlElement("w:" + name)
        el.set(qn("w:val"), val)
        el.set(qn("w:sz"), str(sz))
        el.set(qn("w:space"), "0")
        el.set(qn("w:color"), color)
        borders.append(el)
    tblPr.append(borders)
    for cell in table.rows[header_row].cells:
        tcPr = cell._tc.get_or_add_tcPr()
        tcb = OxmlElement("w:tcBorders")
        el = OxmlElement("w:bottom")
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), "6")
        el.set(qn("w:space"), "0")
        el.set(qn("w:color"), "333333")
        tcb.append(el)
        tcPr.append(tcb)


def set_widths(table, widths_twips):
    table.autofit = False
    tblPr = table._tbl.tblPr
    for el in tblPr.findall(qn("w:tblLayout")):
        tblPr.remove(el)
    layout = OxmlElement("w:tblLayout")
    layout.set(qn("w:type"), "fixed")
    tblPr.append(layout)
    tblW = tblPr.find(qn("w:tblW"))
    if tblW is not None:
        tblW.set(qn("w:w"), str(sum(widths_twips)))
        tblW.set(qn("w:type"), "dxa")
    grid = table._tbl.tblGrid
    for col, w in zip(grid.findall(qn("w:gridCol")), widths_twips):
        col.set(qn("w:w"), str(w))
    for row in table.rows:
        for cell, w in zip(row.cells, widths_twips):
            cell.width = Pt(w / 20.0)


def vmerge(table, r0, r1, c):
    """Vertical merge that leaves the grid regular (no cell is removed)."""
    for r, val in ((r0, "restart"), (r1, "continue")):
        tcPr = table.cell(r, c)._tc.get_or_add_tcPr()
        for el in tcPr.findall(qn("w:vMerge")):
            tcPr.remove(el)
        el = OxmlElement("w:vMerge")
        el.set(qn("w:val"), val)
        tcPr.append(el)


def repeat_header(table, n=1):
    for k in range(n):
        trPr = table.rows[k]._tr.get_or_add_trPr()
        trPr.append(OxmlElement("w:tblHeader"))


def no_split(table):
    for row in table.rows:
        trPr = row._tr.get_or_add_trPr()
        trPr.append(OxmlElement("w:cantSplit"))


def landscape(section):
    section.orientation = WD_ORIENT.LANDSCAPE
    w, h = section.page_width, section.page_height
    section.page_width, section.page_height = max(w, h), min(w, h)
    section.left_margin = Cm(2.0)
    section.right_margin = Cm(2.0)
    section.top_margin = Cm(2.0)
    section.bottom_margin = Cm(2.0)


def portrait(section):
    section.orientation = WD_ORIENT.PORTRAIT
    w, h = section.page_width, section.page_height
    section.page_width, section.page_height = min(w, h), max(w, h)
    section.left_margin = Cm(2.0)
    section.right_margin = Cm(2.0)
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)


def fmt(v, nd=3, dash="\u2014"):
    if v is None or str(v).strip() == "":
        return dash
    return f"{float(v):.{nd}f}"


# ------------------------------------------------------------------ title page
para("Supplementary Material", size=13, bold=True, after=4)
para(TITLE, size=15, bold=True, line=1.3, after=18)
para("Contents", size=12, bold=True, before=6, after=8)

contents = [
    ("Supplementary Table S1",
     "Five evidence layers on the TSH branch: layer sizes and the five-layer intersection "
     "under each threshold configuration."),
    ("Supplementary Table S2",
     "PFNA at the SCRN2 surface: docking, 10 ns all-atom MD and MM-GBSA."),
    ("Supplementary Table S3",
     "Two-sample Mendelian randomisation of the nine instrumentable EDC exposures on four "
     "thyroid outcomes (32 exposure\u2013outcome pairs)."),
]
for head, text in contents:
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.space_before = Pt(4)
    pf.space_after = Pt(4)
    pf.line_spacing = 1.25
    pf.left_indent = Cm(3.5)
    pf.first_line_indent = Cm(-3.5)
    pf.tab_stops.add_tab_stop(Cm(3.5))
    r = p.add_run(head + "\t")
    r.font.name = "Times New Roman"
    r.font.size = Pt(11)
    r.font.bold = True
    r2 = p.add_run(text)
    r2.font.name = "Times New Roman"
    r2.font.size = Pt(11)

para("", size=10, after=12)
para("Machine-readable versions are deposited as Table_S1_five_layer_membership_and_sensitivity.xlsx "
     "(with Table_S1_sensitivity.csv) and Table_S3_EDC_exposure_arm_MR.xlsx. "
     "The per-gene layer membership of all 1,398 screened genes in each configuration, and the "
     "71-gene cross-tissue coherence panel, are carried by the latter and are not reproduced here.",
     size=9, italic=True, line=1.15, after=0)

# ------------------------------------------------------------------ Table S1

sec = doc.add_section(WD_SECTION.NEW_PAGE)
portrait(sec)

para("Supplementary Table S1 | Five evidence layers: layer sizes and the five-layer "
     "intersection under each threshold configuration.", size=10, bold=True, after=4,
     line=1.15, keep=True)
para("Layers are defined on the TSH branch of the screen. Layer 1, whole-blood eQTL MR at "
     "FDR < 0.05; layer 2, thyroid eQTL MR at FDR < 0.05; layer 3, TSH colocalisation at the "
     "stated PP.H4; layer 4, cross-tissue coherence at the stated concordance; layer 5, "
     "assignment as one of the top-ranked eQTL targets of the instrument. Each column relaxes "
     "one or more thresholds and reports the resulting layer sizes. Layer 4 requires the same "
     "direction of effect in every testable tissue (concordance 1.00): at most 48 tissues are "
     "testable for a given instrument, so 0.99 and 1.00 select the same genes and the relaxed "
     "configuration uses 0.90.",
     size=8.5, line=1.1, after=6)

with open(SENS_CSV, encoding="utf-8-sig") as fh:
    s2 = list(csv.DictReader(fh))
assert len(s2) == 5, len(s2)

label = {
    "baseline": "Baseline",
    "PP.H4>=0.70": "PP.H4 \u2265 0.70",
    "top-5 targets": "Top-5 targets",
    "concordance>=90%": "Concordance \u2265 90%",
    "all three relaxed": "All three relaxed",
}
rows_s2 = [
    ("PP.H4 minimum", "PP_H4_min", 1),
    ("Target rank maximum", "target_rank_max", 0),
    ("Concordance minimum", "concordance_min", 2),
    ("Layer 1 \u00b7 whole-blood MR", "layer1_whole_blood_MR", 0),
    ("Layer 2 \u00b7 thyroid MR", "layer2_thyroid_MR", 0),
    ("Layer 3 \u00b7 colocalisation", "layer3_colocalisation", 0),
    ("Layer 4 \u00b7 cross-tissue", "layer4_cross_tissue", 0),
    ("Layer 5 \u00b7 target assignment", "layer5_target", 0),
    ("Five-layer genes, n", "n_five_layer", 0),
    ("Gene passing all five", "five_layer_genes", None),
]

W2 = [2000, 1520, 1520, 1520, 1520, 1560]      # sum 9,640 twips = 17.0 cm
t2 = doc.add_table(rows=1 + len(rows_s2), cols=6)
table_rules(t2)

hdr = ["Threshold or layer"] + [label.get(r["configuration"], r["configuration"]) for r in s2]
for j, h in enumerate(hdr):
    cell_text(t2.cell(0, j), h, bold=True,
              align=WD_ALIGN_PARAGRAPH.LEFT if j == 0 else WD_ALIGN_PARAGRAPH.CENTER,
              before=1.2, after=1.2, line=1.0)

for i, (name, key, nd) in enumerate(rows_s2, start=1):
    if nd is None:
        vals = [r[key] for r in s2]
    elif nd == 0:
        vals = [str(int(round(float(r[key])))) for r in s2]
    else:
        vals = [f"{float(r[key]):.{nd}f}" for r in s2]
    cell_text(t2.cell(i, 0), name)
    for j, v in enumerate(vals, start=1):
        cell_text(t2.cell(i, j), v, align=WD_ALIGN_PARAGRAPH.CENTER)

set_widths(t2, W2)
repeat_header(t2)
no_split(t2)

para("Notes. The five-layer intersection is one gene in every configuration. Counts are genes, "
     "not gene\u2013outcome pairs. The three thresholds relaxed here are the ones a reader would "
     "reasonably have set differently; the two significance thresholds are fixed at FDR < 0.05. "
     "Per-gene membership of all 1,398 screened genes in each configuration, as a layer bitmask, "
     "and the 71-gene cross-tissue coherence panel, are in "
     "Table_S1_five_layer_membership_and_sensitivity.xlsx, generated by "
     "05_scripts/90_figures/42_five_layer_sensitivity.py and "
     "05_scripts/60_expression/41_cross_tissue_coherence_panel.py.",
     size=8, line=1.05, before=6, after=0)

# ------------------------------------------------------------------ Table S2

# the caption is copied from the manuscript, where it is numbered S3; renumber
for _t in s3_caption.iter(qn("w:t")):
    if _t.text and "Supplementary Table S3" in _t.text:
        _t.text = _t.text.replace("Supplementary Table S3", "Supplementary Table S2")

# this table now sits between S1 and S3, so it is placed right after the layer table
_cap = Paragraph(s3_caption, doc)
_cap.paragraph_format.page_break_before = True
_prev = [ch for ch in doc.element.body.iterchildren()
         if ch.tag in (qn('w:p'), qn('w:tbl'))][-1]
for _el in (s3_caption, s3_table, s3_note):
    _prev.addnext(_el)
    _prev = _el

# ------------------------------------------------------------------ Table S3

sec = doc.add_section(WD_SECTION.NEW_PAGE)
landscape(sec)

para("Supplementary Table S3 | Two-sample MR of EDC exposures on thyroid outcomes.",
     size=10, bold=True, before=0, after=4, line=1.15, keep=True)
para("Nine EDC exposures could be instrumented, giving 32 testable exposure\u2013outcome pairs "
     "(47 instrument entries, 13 distinct variants). Estimates are inverse-variance weighted "
     "where two or more instruments were available and a Wald ratio where one was.",
     size=8.5, line=1.1, after=6)

with open(EDC_CSV, encoding="utf-8-sig") as fh:
    rows = list(csv.DictReader(fh))
assert len(rows) == 32, len(rows)

W1 = [1400, 1650, 620, 950, 880, 820, 880, 820, 850, 850, 700, 850, 880, 820, 880]
t1 = doc.add_table(rows=34, cols=15)          # 2 header rows + 32 pairs
table_rules(t1, header_row=1)

top = ["EDC exposure", "Outcome", "n SNP", "Method", "\u03b2", "SE", "P", "FDR q",
       "Cochran Q", "P (Q)", "I\u00b2 (%)", "MDE80 (SD)", "Weighted median", "", ""]
low = [""] * 12 + ["\u03b2", "SE", "P"]
for j, h in enumerate(top):
    cell_text(t1.cell(0, j), h, bold=True,
              align=WD_ALIGN_PARAGRAPH.LEFT if j < 2 else WD_ALIGN_PARAGRAPH.CENTER)
for j, h in enumerate(low):
    cell_text(t1.cell(1, j), h, bold=True, align=WD_ALIGN_PARAGRAPH.CENTER)

for i, r in enumerate(rows, start=2):
    ivw = r["Method"].startswith("IVW")
    vals = [r["EDC exposure"], r["Outcome"], r["n SNP"], "IVW" if ivw else "Wald ratio",
            fmt(r["beta"]), fmt(r["SE"]), fmt(r["P"]), fmt(r["FDR q"]),
            fmt(r["Cochran Q"], 2) if ivw else "\u2014",
            fmt(r["P (Q)"]) if ivw else "\u2014",
            fmt(r["I\u00b2 (%)"], 1) if ivw else "\u2014",
            fmt(r["MDE80 (SD)"], 2),
            fmt(r["beta (WM)"]) if ivw else "\u2014",
            fmt(r["SE (WM)"]) if ivw else "\u2014",
            fmt(r["P (WM)"]) if ivw else "\u2014"]
    for j, v in enumerate(vals):
        cell_text(t1.cell(i, j), v,
                  align=WD_ALIGN_PARAGRAPH.LEFT if j < 2 else WD_ALIGN_PARAGRAPH.CENTER)

set_widths(t1, W1)
for j in range(12):
    vmerge(t1, 0, 1, j)
repeat_header(t1, 2)
no_split(t1)

para("Notes. Outcomes: TSH, thyroid-stimulating hormone; free T4, free thyroxine; "
     "hypothyroidism and Hashimoto thyroiditis as dichotomous endpoints. \u03b2 is in SD units of "
     "the outcome per SD of genetically predicted exposure. FDR q is Benjamini\u2013Hochberg across "
     "the 32 pairs; Cochran Q and I\u00b2 describe heterogeneity between instruments; MDE80 is the "
     "smallest effect detectable with 80% power at the observed instrument strength. Three pairs "
     "reached nominal P < 0.05 (MECPP\u2013hypothyroidism, MEHHP\u2013hypothyroidism, "
     "PFOS\u2013free T4); none survived FDR correction (minimum q = 0.21). Weighted-median "
     "estimates are a pleiotropy-robust sensitivity analysis and are undefined for "
     "single-instrument pairs. Source: b1_edc_thyroid_mr.csv, reproduced by "
     "05_scripts/verify_numbers.py, section 10.",
     size=8, line=1.05, before=6, after=0)

doc.save(OUT)
print("saved", OUT, os.path.getsize(OUT))
