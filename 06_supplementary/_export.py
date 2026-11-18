# -*- coding: utf-8 -*-
"""Export the supplement to PDF and render its pages for inspection."""
import glob
import os

import fitz
import win32com.client

HERE = os.path.dirname(os.path.abspath(__file__))
DOCX = os.path.join(HERE, "Supplementary_Material_SCRN2.docx")
PDF = os.path.join(HERE, "Supplementary_Material_SCRN2.pdf")

word = win32com.client.DispatchEx("Word.Application")
word.Visible = False
word.DisplayAlerts = 0
try:
    d = word.Documents.Open(DOCX, ReadOnly=False)
    d.ExportAsFixedFormat(PDF, 17)
    print("pages:", d.ComputeStatistics(2))
    d.Close(False)
finally:
    word.Quit()

doc = fitz.open(PDF)
print("pdf pages:", len(doc))
for pno, page in enumerate(doc, 1):
    blocks = [b for b in page.get_text("blocks") if b[4].strip()]
    first = blocks[0][4].replace("\n", " ")[:70] if blocks else ""
    last = max((b[3] for b in blocks), default=0)
    print("p%d  bottom=%.0f  %s" % (pno, last, first))
    print("     orientation: %.0f x %.0f" % (page.rect.width, page.rect.height))

out = os.path.join(HERE, "_pages")
os.makedirs(out, exist_ok=True)
for f in glob.glob(os.path.join(out, "*.png")):
    os.remove(f)
for i, page in enumerate(doc, 1):
    page.get_pixmap(dpi=110).save(os.path.join(out, "p%02d.png" % i))
print("rendered", len(doc))
