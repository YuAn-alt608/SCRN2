#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p06_tissue.py  (v2, parsed by column header)
===============================
Per-tissue eQTLs: extract the significant cis-eQTLs of SCRN2 / TBKBP1 / MRPL10 in
17q21.32 (chr17:47.2-48.5 Mb, GRCh38) from the local GTEx v8 bundle (49 tissues).

GTEx only releases significant pairs (permutation-based), so what can be reported
here is whether a given SNP is a significant eQTL in a tissue, its direction,
effect size and TSS distance -- not per-tissue colocalization.
"""
import tarfile, io, gzip, os
import pandas as pd

BASE = "D:/NewBeginning/MR1"
OUT = f"{BASE}/DATA21_complete_2026-11-14"
TAR = f"{BASE}/DATA10/tmp/GTEx_Analysis_v8_eQTL.tar"
GENES = {"ENSG00000141295": "SCRN2", "ENSG00000198933": "TBKBP1", "ENSG00000159111": "MRPL10"}
LO, HI = 47_200_000, 48_500_000
KEYPOS = {"rs6503875": 47838029, "rs7219420": 47722933,
          "rs62076102": 47809766, "rs4794224": 47880097}

rows, egrows = [], []
with tarfile.open(TAR, "r") as tf:
    members = [m for m in tf.getmembers() if m.name.endswith(".txt.gz")]
    print("files in the tar:", len(members), flush=True)
    for m in members:
        base = os.path.basename(m.name)
        tissue = base.split(".v8.")[0]
        fh = tf.extractfile(m)
        if fh is None:
            continue
        gz = gzip.GzipFile(fileobj=io.BufferedReader(fh))
        txt = io.TextIOWrapper(gz, encoding="utf-8", errors="replace")
        head = txt.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(head)}
        if "egenes" in base:
            for line in txt:
                p = line.rstrip("\n").split("\t")
                if len(p) < len(head):
                    continue
                gid = p[ix["gene_id"]].split(".")[0]
                if gid not in GENES:
                    continue
                egrows.append(dict(tissue=tissue, gene=GENES[gid], gene_id=gid,
                                   num_var=p[ix["num_var"]],
                                   pval_nominal=float(p[ix["pval_nominal"]]),
                                   slope=float(p[ix["slope"]]),
                                   variant_id=p[ix["variant_id"]],
                                   rs_id=p[ix["rs_id_dbSNP151_GRCh38p7"]],
                                   maf=float(p[ix["maf"]])))
        else:
            for line in txt:
                p = line.rstrip("\n").split("\t")
                if len(p) < len(head):
                    continue
                vid = p[ix["variant_id"]]
                if not vid.startswith("chr17_"):
                    continue
                try:
                    pos = int(vid.split("_")[1])
                except Exception:
                    continue
                if pos < LO or pos > HI:
                    continue
                gid = p[ix["gene_id"]].split(".")[0]
                if gid not in GENES:
                    continue
                rows.append(dict(tissue=tissue, gene=GENES[gid], variant_id=vid,
                                 pos38=pos, ref=vid.split("_")[2], alt=vid.split("_")[3],
                                 tss_distance=int(p[ix["tss_distance"]]),
                                 maf=float(p[ix["maf"]]),
                                 pval_nominal=float(p[ix["pval_nominal"]]),
                                 pval_threshold=float(p[ix["pval_nominal_threshold"]]),
                                 slope=float(p[ix["slope"]]),
                                 slope_se=float(p[ix["slope_se"]])))
        fh.close()

sp = pd.DataFrame(rows)
sp.to_csv(f"{OUT}/raw/gtex_v8_per_tissue_signif_pairs.csv.gz", index=False, compression="gzip")
eg = pd.DataFrame(egrows)
eg.to_csv(f"{OUT}/tables/09_per_tissue_egenes.csv", index=False, encoding="utf-8-sig")
print("significant pairs:", sp.shape, "| egenes:", eg.shape, flush=True)

sp["key_snp"] = sp.pos38.map({v: k for k, v in KEYPOS.items()})
key = sp[sp.key_snp.notna()].copy()
key.to_csv(f"{OUT}/tables/10_per_tissue_key_snp_eqtl.csv", index=False, encoding="utf-8-sig")

piv = sp.groupby(["tissue", "gene"]).pval_nominal.min().unstack()
piv.to_csv(f"{OUT}/tables/09b_per_tissue_top_eqtl_p.csv", encoding="utf-8-sig")
print("\nmost significant eQTL per tissue (-log10 P):")
print((-piv.apply(lambda c: c.map(lambda v: __import__("math").log10(v) if v and v > 0 else None)))
      .round(1).to_string())
print("\nnumber of tissues hit by each key SNP:")
print(key.groupby(["gene", "key_snp"]).size().to_string())
