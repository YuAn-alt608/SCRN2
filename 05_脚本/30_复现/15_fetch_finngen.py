"""
15_fetch_finngen.py — 从 FinnGen R12 公共 bucket 流式抓取目标 SNP 的关联结果。

不下载 800 MB 全量文件，而是边解压边过滤，只保留命中基因的工具变量。
用法: python scripts/15_fetch_finngen.py E4_HYTHY_AI_STRICT
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
    # 可选：第二个参数指定 SNP 清单，第三个参数指定输出后缀
    # （做阴性对照时需要把"全部工具变量"都抓下来，而不只是显著基因的）
    lst = sys.argv[2] if len(sys.argv) > 2 else None
    suffix = sys.argv[3] if len(sys.argv) > 3 else ""
    targets = load_targets(lst)
    url = URL.format(pheno)
    print(f"[{pheno}] 抓取 {len(targets)} 个目标 SNP <- {url}", flush=True)

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
        # rsids 列可能是逗号分隔的多个 ID
        ids = {x.strip() for x in f[4].split(",")}
        if ids & targets:
            rows.append(s)
    proc.wait()

    print(f"[{pheno}] 扫描 {n:,} 行，命中 {len(rows)} 行", flush=True)
    out = os.path.join(TMP, f"finngen_{pheno}{suffix}.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("#" + "\t".join(header) + "\n")
        fh.write("\n".join(rows) + ("\n" if rows else ""))
    print(f"[{pheno}] 已写 {out}", flush=True)


if __name__ == "__main__":
    main()
