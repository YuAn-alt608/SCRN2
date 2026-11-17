#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
p03_mr_smr.py
=============
For every gene x tissue x outcome:
  1) single-instrument Wald ratio (lead eQTL SNP, delta-method SE)
  2) multi-instrument IVW (instruments with eQTL p<5e-8 and F>10, clumped at r2<0.01), with Q, I2 and the effective number of instruments
  3) SMR (= the lead-SNP Wald ratio) + HEIDI (heterogeneity among instruments clumped at r2<0.1)
Output: tables/04_cisMR_SMR_HEIDI.csv
"""
import os, glob
import numpy as np
import pandas as pd
from scipy import stats

BASE = "D:/NewBeginning/MR1/DATA21_complete_2026-11-14"
IN, OUT = f"{BASE}/inputs", f"{BASE}/tables"
LD_DIR = f"{BASE}/inputs_ld"
_ld = {}


def get_ld(oc):
    if oc in _ld:
        return _ld[oc]
    f = f"{LD_DIR}/ld_{oc}.tsv.gz"
    if not os.path.exists(f):
        _ld[oc] = None
        return None
    d = pd.read_csv(f, sep="\t", index_col=0)
    d.index = d.index.astype(int)
    d.columns = d.columns.astype(int)
    _ld[oc] = d
    return d


def wald(bg, seg, bo, seo):
    if not np.isfinite(bg) or not np.isfinite(bo) or bg == 0 or bo == 0:
        return np.nan, np.nan
    b = bo / bg
    se = abs(b) * np.sqrt((seo / bo) ** 2 + (seg / bg) ** 2)
    return b, se


def clump(pos, z, R, r2_thr):
    order = np.argsort(-np.abs(z))
    sel = []
    for i in order:
        if all(R[pos[i], pos[j]] ** 2 < r2_thr for j in sel):
            sel.append(i)
    return sel


def ivw(b, se):
    w = 1.0 / se ** 2
    bb = float(np.sum(w * b) / np.sum(w))
    ss = float(np.sqrt(1.0 / np.sum(w)))
    Q = float(np.sum(w * (b - bb) ** 2))
    df = len(b) - 1
    pq = float(stats.chi2.sf(Q, df)) if df > 0 else np.nan
    i2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 and df > 0 else 0.0
    return bb, ss, Q, pq, i2


def eff_k(Rsub):
    if Rsub.shape[0] < 2:
        return 1.0
    ev = np.clip(np.linalg.eigvalsh(Rsub), 1e-8, None)
    return float(ev.sum() ** 2 / (ev ** 2).sum())


def heidi(bs, ses):
    k = len(bs)
    if k < 2:
        return np.nan, np.nan, 0
    T, npair = 0.0, 0
    for i in range(k):
        for j in range(i + 1, k):
            d = bs[i] - bs[j]
            T += d ** 2 / (ses[i] ** 2 + ses[j] ** 2)
            npair += 1
    return T, float(stats.chi2.sf(T, df=npair)), npair


rows = []
for f in sorted(glob.glob(f"{IN}/*.tsv")):
    tag = os.path.basename(f)[:-4]
    tissue, gene, oc = tag.split("__")
    d = pd.read_csv(f, sep="\t")
    if len(d) < 30:
        continue
    LD = get_ld(oc)
    d["pos37"] = d.pos37.astype(int)
    if LD is not None:
        d = d[d.pos37.isin(LD.index)].copy()
    if len(d) < 30:
        continue
    d = d.drop_duplicates("pos37").reset_index(drop=True)
    ii = d.p_eqtl.astype(float).idxmin()
    df_lead = d.loc[ii]

    b_lead, se_lead = wald(df_lead.beta_eqtl, df_lead.se_eqtl, df_lead.beta_out, df_lead.se_out)
    z_lead = b_lead / se_lead
    p_lead = 2 * stats.norm.sf(abs(z_lead))

    # instruments
    inst = d[(d.p_eqtl < 5e-8) & ((d.beta_eqtl / d.se_eqtl) ** 2 > 10)].copy()
    n_iv, b_ivw, se_ivw, Q, pq, i2, k_eff = 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    bs = ses = None
    if LD is not None and len(inst) >= 1:
        inst = inst.reset_index(drop=True)
        p_ok = [p for p in inst.pos37 if p in LD.index]
        inst = inst[inst.pos37.isin(p_ok)].reset_index(drop=True)
        if len(inst):
            posmap = {p: i for i, p in enumerate(LD.index)}
            R = LD.values if LD.shape[0] == LD.shape[1] else None
            pidx = [posmap[p] for p in inst.pos37]
            sel = clump(pidx, (inst.beta_eqtl / inst.se_eqtl).values, R, 0.01)
            inst = inst.iloc[sel].reset_index(drop=True)
            n_iv = len(inst)
            if n_iv >= 1:
                bs, ses = [], []
                for _, r in inst.iterrows():
                    b, s = wald(r.beta_eqtl, r.se_eqtl, r.beta_out, r.se_out)
                    bs.append(b); ses.append(s)
                bs, ses = np.array(bs), np.array(ses)
                ok = np.isfinite(bs) & np.isfinite(ses) & (ses > 0)
                bs, ses = bs[ok], ses[ok]
                inst = inst.iloc[np.where(ok)[0]].reset_index(drop=True)
                n_iv = len(inst)
                if n_iv == 0:
                    continue
                if n_iv >= 2:
                    b_ivw, se_ivw, Q, pq, i2 = ivw(bs, ses)
                    sub = R[np.ix_([posmap[p] for p in inst.pos37],
                                   [posmap[p] for p in inst.pos37])]
                    k_eff = eff_k(sub)
                else:
                    b_ivw, se_ivw = bs[0], ses[0]
                    k_eff = 1.0

    # HEIDI: clumped at r2<0.1
    he_T = he_p = np.nan; he_n = 0; he_k = 0; he_ndrop = 0
    if LD is not None and len(inst) >= 2:
        posmap = {p: i for i, p in enumerate(LD.index)}
        R = LD.values
        pidx = [posmap[p] for p in inst.pos37]
        sel = clump(pidx, (inst.beta_eqtl / inst.se_eqtl).values, R, 0.1)
        bb, ss_ = [], []
        for i in sel:
            r = inst.iloc[i]
            b, s = wald(r.beta_eqtl, r.se_eqtl, r.beta_out, r.se_out)
            if np.isfinite(b) and np.isfinite(s) and s > 0:
                bb.append(b); ss_.append(s)
        he_k = len(bb)
        he_T, he_p, he_n = heidi(bb, ss_)

    rows.append(dict(tissue=tissue, gene=gene, outcome=oc, nsnp=len(d),
                     lead_snp=df_lead.snp, lead_pos37=int(df_lead.pos37),
                     lead_p_eqtl=float(df_lead.p_eqtl), F_lead=float((df_lead.beta_eqtl / df_lead.se_eqtl) ** 2),
                     b_wald=b_lead, se_wald=se_lead, p_wald=p_lead,
                     n_iv=n_iv, k_eff=k_eff, b_ivw=b_ivw, se_ivw=se_ivw,
                     p_ivw=(2 * stats.norm.sf(abs(b_ivw / se_ivw)) if n_iv >= 2 else np.nan),
                     Q=Q, p_Q=pq, I2=i2,
                     b_smr=b_lead, se_smr=se_lead, p_smr=p_lead,
                     heidi_T=he_T, heidi_p=he_p, heidi_npair=he_n,
                     heidi_k=he_k, heidi_ndrop=n_iv - he_k))

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/04_cisMR_SMR_HEIDI.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250, "display.max_columns", 60)
show = res[res.outcome == "TSH"].sort_values("p_wald")
print(show[["tissue", "gene", "lead_snp", "lead_p_eqtl", "n_iv", "k_eff", "b_wald", "se_wald",
            "p_wald", "p_ivw", "p_Q", "heidi_p"]].to_string(index=False))
print("\nwrote", len(res), "rows")
