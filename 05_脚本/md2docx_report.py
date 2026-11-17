# -*- coding: utf-8 -*-
"""Markdown -> docx 渲染器（专用于本项目的修改指导类报告）

用法：
    python md2docx_report.py <input.md> <output.docx>

设计要点
  * run() 用标记栈做分词，支持 **加粗** 与 `行内代码` 任意嵌套（旧版 re.split 在
    "加粗里含代码" / "代码里含加粗" 时会漏出原始标记字符）
  * 表格：按列数指定 cm 列宽、行不跨页(cantSplit)、跨页重复表头(tblHeader)、
    表头行 keep_with_next（避免表头孤悬页底）
  * 标题：keep_with_next，避免标题落在页底而正文翻页
"""
import re
import sys
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

LAT = "Times New Roman"
CJK = "SimSun"
MONO = "Consolas"

WIDTHS = {
    3: [2.4, 7.4, 6.8],
    4: [1.2, 4.0, 8.4, 3.0],
    5: [1.25, 3.45, 5.45, 3.35, 2.9],
}
TOTAL_CM = 16.5  # A4 21.0 − 左右各 2.2 边距


def auto_widths(rows, total=TOTAL_CM, mn=1.15, mx=9.4):
    """按各列最长内容分配列宽：中文按 2 个字符宽计，开 0.55 次方压平极值。"""
    ncol = max(len(r) for r in rows)

    def vlen(s):
        return sum(2 if ord(ch) > 127 else 1 for ch in s)

    raw = []
    for c in range(ncol):
        L = max([vlen(r[c]) for r in rows if c < len(r)] + [4])
        raw.append(min(L, 110) ** 0.55)
    s = sum(raw)
    w = [total * x / s for x in raw]
    for _ in range(6):
        w = [min(max(x, mn), mx) for x in w]
        s = sum(w)
        if abs(s - total) < 0.02:
            break
        w = [x * total / s for x in w]
    return w

ACCENT = RGBColor(0x1F, 0x3B, 0x5C)
MUTED = RGBColor(0x44, 0x4A, 0x52)


# ---------------------------------------------------------------- tokenizer
def tokenize(text):
    """把一段 markdown 文本切成 [(text, bold, mono, italic), ...]，支持四种标记任意嵌套。"""
    out, buf = [], []
    bold = mono = ital = False

    def flush():
        if buf:
            out.append(("".join(buf), bold, mono, ital))
            buf.clear()

    i, n = 0, len(text)
    while i < n:
        if text.startswith("**", i):
            flush()
            bold = not bold
            i += 2
            continue
        if text[i] == "*":
            flush()
            ital = not ital
            i += 1
            continue
        if text[i] == "`":
            flush()
            mono = not mono
            i += 1
            continue
        buf.append(text[i])
        i += 1
    flush()
    return [t for t in out if t[0]]


def run(p, text, bold=False, size=None, color=None, italic=False):
    for seg, b, m, it in tokenize(text):
        r = p.add_run(seg)
        r.bold = bold or b
        r.italic = italic or it
        r.font.name = MONO if m else LAT
        r.font.size = Pt(size) if size else Pt(8.6 if m else 10.5)
        if color is not None:
            r.font.color.rgb = color
        rPr = r._element.get_or_add_rPr()
        rFonts = rPr.find(qn("w:rFonts"))
        if rFonts is None:
            rFonts = OxmlElement("w:rFonts")
            rPr.append(rFonts)
        rFonts.set(qn("w:eastAsia"), CJK)


# ---------------------------------------------------------------- helpers
def shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), hexcolor)
    tcPr.append(shd)


def bottom_border(p, size=6, color="1F3B5C"):
    pPr = p._p.get_or_add_pPr()
    pBdr = OxmlElement("w:pBdr")
    b = OxmlElement("w:bottom")
    b.set(qn("w:val"), "single")
    b.set(qn("w:sz"), str(size))
    b.set(qn("w:space"), "2")
    b.set(qn("w:color"), color)
    pBdr.append(b)
    pPr.append(pBdr)


def keep_next(p):
    p.paragraph_format.keep_with_next = True
    return p


def mk_table(doc, rows):
    n = max(len(r) for r in rows)
    t = doc.add_table(rows=0, cols=n)
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    w = auto_widths(rows)
    for ri, r in enumerate(rows):
        row = t.add_row()
        trPr = row._tr.get_or_add_trPr()
        trPr.append(OxmlElement("w:cantSplit"))
        if ri == 0:
            trPr.append(OxmlElement("w:tblHeader"))
        for ci in range(n):
            cell = row.cells[ci]
            if w:
                cell.width = Cm(w[ci])
            p = cell.paragraphs[0]
            pf = p.paragraph_format
            pf.space_before = Pt(0)
            pf.space_after = Pt(0.5)
            pf.line_spacing = 1.06
            if ri == 0:
                pf.keep_with_next = True
            run(p, r[ci] if ci < len(r) else "", bold=(ri == 0), size=8.6)
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            if ri == 0:
                shade(cell, "E8EDF3")
    sp = doc.add_paragraph()
    sp.paragraph_format.space_after = Pt(2)
    sp.paragraph_format.line_spacing = 1.0
    for r in sp.runs:
        r.font.size = Pt(4)
    return t


# ---------------------------------------------------------------- main
def build(md_path, out_path):
    lines = open(md_path, encoding="utf-8").read().split("\n")
    doc = Document()

    sec = doc.sections[0]
    sec.top_margin = Cm(2.0)
    sec.bottom_margin = Cm(2.0)
    sec.left_margin = Cm(2.2)
    sec.right_margin = Cm(2.2)

    st = doc.styles["Normal"]
    st.font.name = LAT
    st.font.size = Pt(10.5)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), CJK)
    st.paragraph_format.space_after = Pt(3.5)
    st.paragraph_format.line_spacing = 1.21

    i = 0
    N = len(lines)
    while i < N:
        raw = lines[i]
        s = raw.strip()
        if not s:
            i += 1
            continue

        # 分隔线：章节已由标题底边框分隔，跳过以免线条堆叠
        if re.fullmatch(r"(-{3,}|\*{3,}|_{3,})", s):
            i += 1
            continue

        # 表格：连续的 | 行
        if s.startswith("|"):
            block = []
            while i < N and lines[i].strip().startswith("|"):
                block.append(lines[i].strip())
                i += 1
            rows = []
            for b in block:
                if re.fullmatch(r"\|[\s:\-|]+\|", b):
                    continue  # 分隔行
                cells = [c.strip() for c in b.strip("|").split("|")]
                rows.append(cells)
            if rows:
                mk_table(doc, rows)
            continue

        # 标题
        m = re.match(r"^(#{1,4})\s+(.*)$", s)
        if m:
            lvl, txt = len(m.group(1)), m.group(2)
            p = doc.add_paragraph()
            pf = p.paragraph_format
            pf.keep_with_next = True
            if lvl == 1:
                pf.space_before = Pt(0)
                pf.space_after = Pt(8)
                run(p, txt, bold=True, size=19, color=ACCENT)
                bottom_border(p, size=12)
            elif lvl == 2:
                pf.space_before = Pt(15)
                pf.space_after = Pt(5)
                run(p, txt, bold=True, size=13.5, color=ACCENT)
                bottom_border(p, size=5, color="A9B6C6")
            elif lvl == 3:
                pf.space_before = Pt(11)
                pf.space_after = Pt(3.5)
                run(p, txt, bold=True, size=11, color=ACCENT)
            else:
                pf.space_before = Pt(8)
                pf.space_after = Pt(3)
                run(p, txt, bold=True, size=10.5, color=MUTED)
            i += 1
            continue

        # 引用
        if s.startswith(">"):
            p = doc.add_paragraph()
            pf = p.paragraph_format
            pf.left_indent = Cm(0.5)
            pf.space_before = Pt(4)
            pf.space_after = Pt(4)
            run(p, s.lstrip("> ").strip(), size=9.8, color=MUTED, italic=True)
            i += 1
            continue

        # 无序列表
        if s.startswith(("- ", "* ")):
            p = doc.add_paragraph()
            pf = p.paragraph_format
            pf.left_indent = Cm(0.62)
            pf.first_line_indent = Cm(-0.62)
            pf.space_after = Pt(2.6)
            run(p, "•  ")
            run(p, s[2:].strip())
            i += 1
            continue

        # 有序列表
        m = re.match(r"^(\d+)\.\s+(.*)$", s)
        if m:
            p = doc.add_paragraph()
            pf = p.paragraph_format
            pf.left_indent = Cm(0.75)
            pf.first_line_indent = Cm(-0.75)
            pf.space_after = Pt(2.6)
            run(p, m.group(1) + ".  ")
            run(p, m.group(2).strip())
            i += 1
            continue

        # 缩进续行
        if raw.startswith("  ") and not s.startswith("- "):
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.75)
            run(p, s, size=9.8, color=MUTED)
            i += 1
            continue

        # 普通段落
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(4)
        run(p, s)
        i += 1

    doc.save(out_path)
    return doc


if __name__ == "__main__":
    a = sys.argv[1:]
    d = build(a[0], a[1])
    print("paras %d tables %d" % (len(d.paragraphs), len(d.tables)))
