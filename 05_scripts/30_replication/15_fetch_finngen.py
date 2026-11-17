"""
15_fetch_finngen.py - stream the association results for target SNPs out of the FinnGen R12 bucket.

The 800 MB release is not downloaded; it is decompressed and filtered on the fly, keeping only the
instruments of the genes of interest.
Usage: python 05_scripts/30_replication/15_fetch_finngen.py E4_HYTHY_AI_STRICT
"""
import os
import sys
import gzip
import subprocess

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP = os.path.join(BASE, "tmp")
URL = ("https://storage.googleapis.com/finngen-public-data-r12/"
       "summary_stats/release/finngen_R12_{}.gz")


def load_targets(p=None):
    p = p or os.path.join(TMP, "fg_snp_list.txt")
    return {l.strip() for l in open(p) if l.strip()}


def main():
    pheno = sys.argv[1]
    # optional: the second argument names a SNP list, the third an output suffix
    # (the negative control needs every instrument pulled down, not just those of the significant genes)
    lst = sys.argv[2] if len(sys.argv) > 2 else None
    suffix = sys.argv[3] if len(sys.argv) > 3 else ""
    targets = load_targets(lst)
    url = URL.format(pheno)
    print(f"[{pheno}] fetching {len(targets)} target SNPs <- {url}", flush=True)

    proc = subprocess.Popen(["curl", "-sS", "--max-time", "2400", url],
                            stdout=subprocess.PIPE)
    gz = gzip.GzipFile(fileobj=proc.stdout)
    header = None
    rows = []
    n = 0
    for line in gz:
        n += 1
        s = line.decode("utf-8", "replace").rstrip("\n")
        if s.startswith("#"):
            header = s.lstrip("#").split("\t")
            continue
        f = s.split("\t")
        if len(f) < 11:
            continue
        # the rsids column may hold several comma-separated IDs
        ids = {x.strip() for x in f[4].split(",")}
        if ids & targets:
            rows.append(s)
    proc.wait()

    print(f"[{pheno}] scanned {n:,} rows, {len(rows)} matched", flush=True)
    out = os.path.join(TMP, f"finngen_{pheno}{suffix}.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("#" + "\t".join(header) + "\n")
        fh.write("\n".join(rows) + ("\n" if rows else ""))
    print(f"[{pheno}] written {out}", flush=True)


if __name__ == "__main__":
    main()
