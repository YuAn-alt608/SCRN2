# -*- coding: utf-8 -*-
"""
Step 24: GTEx multi-tissue control - how specific is the "tissue validation" in thyroid

Question: section 06 looked at thyroid alone, so a 34% detection rate cannot be judged high or
      low. If the same instruments were also detected in 30% of liver, muscle or pancreas, the
      "thyroid validation" would not hold.

Approach: for all 70 A/B-tier (gene, instrument) pairs, query GTEx v8 for the eQTL records of
      each variant across all tissues (singleTissueEqtl, filtered by geneSymbol to the target
      gene), then compute the detection rate (p<0.05) and directional agreement per tissue and
      see where thyroid ranks.

Source: GTEx Portal API v2 with datasetId=gtex_v8. Checked: the thyroid record for RERE/rs301802
      is p=7.04189e-09, nes=0.148, identical to the local Thyroid.v8.signif_variant_gene_pairs
      file, so the API and v8 are the same source. Note that v2 returns gtex_v10 by default, so
      datasetId=gtex_v8 must be passed explicitly.

Variant ID: the API accepts only variantId (chr_pos_ref_alt_b38), built from the NCBI cache; on
        failure the REF/ALT pair is swapped, since a multi-allelic site may exist in either form.

Output:
  results/gtex_tissue_control.csv    long table (gene, snp, tissue, p, nes)
  results/gtex_tissue_summary.txt    detection rate per tissue and the thyroid rank
"""
import os, json, time, subprocess
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
DATA = os.path.join(BASE, "data")
TMP = os.path.join(BASE, "tmp")
CACHE = os.path.join(TMP, "gtex_api_cache.json")
VAR_CACHE = os.path.join(TMP, "ensembl_var_cache.json")
API = "https://gtexportal.org/api/v2/association/singleTissueEqtl"


def query(vid, tries=3):
    url = f"{API}?variantId={vid}&datasetId=gtex_v8&itemsPerPage=1000"
    for _ in range(tries):
        try:
            p = subprocess.run(["curl", "-sS", "--max-time", "30", url],
                               capture_output=True, text=True, timeout=40)
            if p.returncode == 0 and p.stdout.strip():
                return json.loads(p.stdout).get("data", [])
        except Exception:
            pass
        time.sleep(0.5)
    return None


def main():
    tier = pd.read_csv(os.path.join(RES, "evidence_tiers_v2.csv"), encoding="utf-8-sig")
    tier.columns = [c.strip() for c in tier.columns]
    genes = sorted(tier.loc[tier["tier"].isin(["A_multi", "B_strong"]), "gene"].unique())

    iv = pd.read_csv(os.path.join(DATA, "ivs.tsv"), sep="\t")
    pairs = iv[iv["gene"].isin(genes)][["gene", "snp"]].drop_duplicates()
    print(f"{len(genes)} A/B-tier genes, {len(pairs)} (gene, SNP) pairs", flush=True)

    vcache = json.load(open(VAR_CACHE, encoding="utf-8"))
    cache = json.load(open(CACHE, encoding="utf-8")) if os.path.exists(CACHE) else {}

    rows, fail = [], []
    for i, r in enumerate(pairs.itertuples(), 1):
        key = f"{r.gene}|{r.snp}"
        if key in cache:
            recs = cache[key]
        else:
            m = vcache.get(r.snp)
            recs = None
            if m:
                for ref, alt in [(m["ref"], m["alt"]), (m["alt"], m["ref"])]:
                    vid = f"chr{m['chr']}_{m['pos']}_{ref}_{alt}_b38"
                    recs = query(vid)
                    if recs:
                        break
            cache[key] = recs
            if (i % 10) == 0:
                json.dump(cache, open(CACHE, "w", encoding="utf-8"))
                print(f"  {i}/{len(pairs)}", flush=True)
        if recs is None:
            fail.append(key)
            continue
        for x in recs:
            if x.get("geneSymbol") != r.gene:
                continue
            rows.append(dict(gene=r.gene, snp=r.snp,
                             tissue=x["tissueSiteDetailId"],
                             p=x["pValue"], nes=x["nes"]))
    json.dump(cache, open(CACHE, "w", encoding="utf-8"))

    d = pd.DataFrame(rows).drop_duplicates(["gene", "snp", "tissue"])
    d.to_csv(os.path.join(RES, "gtex_tissue_control.csv"), index=False,
             encoding="utf-8-sig")

    iv2 = iv.drop_duplicates("snp").set_index("snp")
    d["ea"] = d["snp"].map(iv2["ea"])
    d["z"] = d["snp"].map(iv2["z"])
    d["same"] = (d["nes"] * d["z"]) > 0

    n_pairs = len(pairs)
    sig = d[d["p"] < 0.05]
    per = sig.groupby("tissue").agg(n_sig=("snp", "size"),
                                    n_same=("same", "sum"))
    per["detect_pct"] = 100 * per["n_sig"] / n_pairs
    per["same_pct"] = 100 * per["n_same"] / per["n_sig"]
    per = per.sort_values("n_sig", ascending=False)
    per.to_csv(os.path.join(RES, "gtex_tissue_rates.csv"), encoding="utf-8-sig")

    thy = per.loc["Thyroid"] if "Thyroid" in per.index else None
    rank = (per["n_sig"] > thy["n_sig"]).sum() + 1 if thy is not None else None
    n_tis = d["tissue"].nunique()
    med = per["n_sig"].median()

    lines = [
        f"A/B-tier (gene, SNP) pairs: {n_pairs}; failed queries: {len(fail)}" + (f" ({fail[:5]})" if fail else ""),
        f"tissues covered by GTEx v8: {n_tis}; tissues with a significant eQTL in the table: {len(per)}",
        "",
    ]
    if thy is not None:
        lines += [
            f"thyroid: detected {int(thy['n_sig'])}/{n_pairs} ({thy['detect_pct']:.1f}%),"
            f" same direction {int(thy['n_same'])}/{int(thy['n_sig'])} ({thy['same_pct']:.1f}%)",
            f"thyroid detection rank: {rank} of {len(per)} (ties counted at the upper bound)",
            f"median detection across tissues: {med:.0f}/{n_pairs} ({100*med/n_pairs:.1f}%)",
            "",
            "8 tissues with the highest detection count:",
        ]
        for t, r in per.head(8).iterrows():
            lines.append(f"  {t:<38s} {int(r['n_sig']):>3d}/{n_pairs} "
                         f"({r['detect_pct']:.0f}%)  same dir. {r['same_pct']:.0f}%")
        lines += ["", "5 tissues with the lowest detection count:"]
        for t, r in per.tail(5).iterrows():
            lines.append(f"  {t:<38s} {int(r['n_sig']):>3d}/{n_pairs} "
                         f"({r['detect_pct']:.0f}%)  same dir. {r['same_pct']:.0f}%")
        above = per[per["n_sig"] > thy["n_sig"]]
        if len(above):
            lines += ["", f"tissues detecting more than thyroid ({len(above)}): " + ", ".join(above.index[:12])]
    txt = "\n".join(lines)
    print("\n" + txt)
    open(os.path.join(RES, "gtex_tissue_summary.txt"), "w", encoding="utf-8").write(txt + "\n")


if __name__ == "__main__":
    main()
