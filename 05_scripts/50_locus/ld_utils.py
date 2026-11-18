#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ld_utils.py  (v2, position-matched)
===================================
The ID column of the 1000G phase3 VCFs is almost entirely ".", with no rsIDs, so LD
cannot be looked up by rsID. Matching is therefore done on **GRCh37 position**: r2 is
insensitive to the direction in which alleles are coded (swapping them leaves r2
unchanged), so a position match is sufficient and ref/alt need not be reconciled.

Only EUR samples are parsed (CEU/FIN/GBR/IBS/TSI, 503 individuals), and variants with
MAF<0.05 are dropped to keep the scan fast.
"""
import os, gzip, pickle
import numpy as np

BASE = "D:/NewBeginning/MR1/DATA15_thyroid_eQTL"
VCF = os.path.join(BASE, "tmp/1000g")
CACHE = os.path.join(BASE, "tmp/ld_cache")
os.makedirs(CACHE, exist_ok=True)
PANEL = "D:/NewBeginning/MR1/DATA10/raw/1000g/integrated_call_samples_v3.20130502.ALL.panel"
EUR = {"CEU", "FIN", "GBR", "IBS", "TSI"}
MAF_MIN = 0.05


def eur_indices(n_samples):
    if not os.path.exists(PANEL):
        return np.arange(n_samples)
    idx, i = [], 0
    with open(PANEL) as f:
        next(f, None)
        for line in f:
            p = line.split()
            if len(p) >= 2:
                if p[1] in EUR:
                    idx.append(i)
                i += 1
    return np.array(idx) if idx else np.arange(n_samples)


def get_region(chrom, start, end, force=False):
    tag = f"v2_chr{chrom}_{start}_{end}"
    fn = os.path.join(CACHE, tag + ".pkl")
    if os.path.exists(fn) and not force:
        with open(fn, "rb") as f:
            return pickle.load(f)

    path = os.path.join(VCF, f"chr{chrom}.vcf.gz")
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    pos_l, gt_l = [], []
    eur = None
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                if line.startswith("#CHROM"):
                    n = len(line.rstrip("\n").split("\t")) - 9
                    eur = eur_indices(n)
                continue
            parts = line.rstrip("\n").split("\t")
            pp = int(parts[1])
            if pp < start:
                continue
            if pp > end:
                break
            if "," in parts[4]:
                continue
            g = np.empty(len(eur), dtype=np.int8)
            for k, j in enumerate(eur):
                gt = parts[9 + j]
                c = gt.find(":")
                if c >= 0:
                    gt = gt[:c]
                if gt in ("0|0", "0/0"):
                    g[k] = 0
                elif gt in ("0|1", "1|0", "0/1", "1/0"):
                    g[k] = 1
                elif gt in ("1|1", "1/1"):
                    g[k] = 2
                else:
                    g[k] = -1
            miss = (g < 0).mean()
            if miss > 0.05:
                continue
            af = g[g >= 0].mean() / 2.0
            maf = min(af, 1 - af)
            if maf < MAF_MIN:
                continue
            pos_l.append(pp)
            gt_l.append(g)

    reg = dict(chrom=chrom, start=start, end=end,
               pos=np.array(pos_l, dtype=np.int64),
               G=np.array(gt_l, dtype=np.int8) if gt_l else np.zeros((0, 0), np.int8))
    with open(fn, "wb") as f:
        pickle.dump(reg, f)
    return reg


def ld_by_position(reg, positions):
    """Look up r2 by GRCh37 position. Returns (matched positions, r2 matrix)."""
    idx = {p: i for i, p in enumerate(reg["pos"])}
    keep = [(p, idx[p]) for p in positions if p in idx]
    if not keep:
        return np.array([], dtype=np.int64), np.zeros((0, 0))
    sel = [k[1] for k in keep]
    G = reg["G"][sel].astype(np.float64)
    mu = np.array([g[g >= 0].mean() for g in G])
    G = np.where(G < 0, mu[:, None], G)
    Z = G - G.mean(axis=1, keepdims=True)
    sd = Z.std(axis=1, keepdims=True)
    sd[sd == 0] = 1
    Z = Z / sd
    R = (Z @ Z.T) / Z.shape[1]
    return np.array([k[0] for k in keep]), R ** 2
