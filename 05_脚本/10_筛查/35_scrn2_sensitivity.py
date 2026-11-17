#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SCRN2 敏感性分析：Steiger 方向性 + MR-Egger 截距 + 留一法 + 多工具变量 MR
=====================================================================
目标：强化 "SCRN2 表达 → 甲状腺功能" 的方向性，排除反向因果。

方法：
1. Steiger 方向性检验：比较 eQTL 解释的方差 (R²_exposure) 与结局解释的方差
   (R²_outcome)，判断因果方向（暴露应解释更多方差）。
2. 多工具变量 MR（clump 后独立 SNP）：IVW + MR-Egger（含截距）+ 加权中位数。
   Egger 截距用于检测水平多效性。
3. 留一法：逐 SNP 剔除后重估，检验是否由单点驱动。

数据：tmp/susie/SCRN2__<outcome>.tsv（harmonized 的 eQTL beta/varbeta + 结局 beta/se）
      tmp/susie/SCRN2__<outcome>.ld.csv.gz（LD 矩阵）
"""
import pandas as pd
import numpy as np
import os

os.makedirs("results", exist_ok=True)

# ============ 配置 ============
OUTCOMES = {
    "TSH":  dict(n_out=482873, type="quant"),
    "HYPO": dict(n_out=1046991, s=98805/1046991, type="cc"),
    "HASH": dict(n_out=395640, s=15654/395640, type="cc"),
}
IV = "rs11079797"

# ============ 工具函数 ============
def r2_from_z(z, n, maf):
    """从 z-score 反推方差解释比例 R² = z²/(z² + n - 2)。"""
    return z**2 / (z**2 + n - 2)

def steiger_test(beta_eqtl, se_eqtl, n_eqtl, maf, beta_out, se_out, n_out):
    """Steiger 方向性检验：暴露应解释更多方差。"""
    z_eqtl = beta_eqtl / se_eqtl
    z_out = beta_out / se_out
    r2_eqtl = r2_from_z(z_eqtl, n_eqtl, maf)
    r2_out = r2_from_z(z_out, n_out, maf)
    # Steiger Z 检验
    # 用标准公式：比较两个 R²（经校正）
    return r2_eqtl, r2_out, r2_eqtl > r2_out

def clump_independent(snps, R, pval, r2_cut=0.1):
    """贪心 clump：按 p 值排序，保留 LD r² < cut 的 SNP。"""
    order = np.argsort(pval)
    keep = []
    for i in order:
        if all(abs(R[i, j]) < np.sqrt(r2_cut) for j in keep):
            keep.append(i)
    return sorted(keep)

def ivw(beta, se):
    """逆方差加权。"""
    w = 1 / se**2
    b = np.sum(beta * w) / np.sum(w)
    se_b = np.sqrt(1 / np.sum(w))
    return b, se_b

def mr_egger(beta_x, se_x, beta_y, se_y):
    """MR-Egger 回归（含截距）。"""
    # 以 se_y 为权重做加权回归：beta_y ~ beta_x
    w = 1 / se_y**2
    x = beta_x
    y = beta_y
    # 加权线性回归
    W = np.diag(w)
    X = np.column_stack([np.ones(len(x)), x])
    XtW = X.T @ W
    beta_hat = np.linalg.solve(XtW @ X, XtW @ y)
    intercept, slope = beta_hat
    # 标准误
    resid = y - X @ beta_hat
    dof = len(x) - 2
    if dof > 0:
        sigma2 = (resid @ W @ resid) / dof
        cov = np.linalg.inv(XtW @ X) * sigma2
        se_intercept = np.sqrt(cov[0, 0])
        se_slope = np.sqrt(cov[1, 1])
    else:
        se_intercept = se_slope = np.nan
    return intercept, se_intercept, slope, se_slope

def weighted_median(beta, se, weights=None):
    """加权中位数估计。"""
    if weights is None:
        weights = np.ones(len(beta))
    order = np.argsort(beta)
    beta_s = beta[order]
    w_s = weights[order]
    cum = np.cumsum(w_s)
    total = cum[-1]
    # 找到累计权重过半的区间
    idx = np.searchsorted(cum, total/2)
    return beta_s[idx]

# ============ 主分析 ============
print("=" * 72)
print("SCRN2 敏感性分析（Steiger + MR-Egger + 留一法）")
print("=" * 72)

all_rows = []

for oc, cfg in OUTCOMES.items():
    f = f"tmp/susie/SCRN2__{oc}.tsv"
    ld_f = f"tmp/susie/SCRN2__{oc}.ld.csv.gz"
    if not os.path.exists(f) or not os.path.exists(ld_f):
        print(f"\n[跳过] {oc}: 无数据")
        continue

    d = pd.read_csv(f, sep="\t")
    R = pd.read_csv(ld_f, header=None).values
    n_eqtl = int(d["n"].median())

    print(f"\n{'='*72}\n【{oc}】 n_eQTL={n_eqtl}, n_outcome={cfg['n_out']}\n{'='*72}")

    # ---- 1. Steiger 方向性检验（工具变量本身）----
    iv = d[d.snp == IV].iloc[0]
    z_eqtl_iv = iv["z"]
    z_out_iv = iv["beta"] / iv["se"]
    r2_e, r2_o = r2_from_z(z_eqtl_iv, n_eqtl, iv["maf"]), r2_from_z(z_out_iv, cfg["n_out"], iv["maf"])
    steiger_ok = r2_e > r2_o
    print(f"\n[1] Steiger 方向性检验（工具变量 {IV}）")
    print(f"    eQTL: z={z_eqtl_iv:.1f}, R²={r2_e:.3e}")
    print(f"    结局: z={z_out_iv:.1f}, R²={r2_o:.3e}")
    print(f"    R²(exposure) > R²(outcome): {steiger_ok}  → {'方向正确(eQTL→结局)' if steiger_ok else '方向存疑(可能反向因果)'}")

    # ---- 2. clump 得到独立工具变量 ----
    pval = d["p"].values
    keep = clump_independent(d["snp"].values, R, pval, r2_cut=0.1)
    # 确保工具变量在内（或至少保留它）
    iv_pos = list(d["snp"]).index(IV)
    if iv_pos not in keep:
        keep = sorted(set(keep) | {iv_pos})
    d_ind = d.iloc[keep].reset_index(drop=True)

    # 结局 beta 方向对齐：确保效应等位基因一致（susie 数据已 harmonized，ea 列是 eQTL 的）
    # susie 里 beta_eqtl 和 beta 都是对同一 ea 的效应
    beta_x = d_ind["beta_eqtl"].values
    se_x = np.sqrt(d_ind["varbeta_eqtl"].values)
    beta_y = d_ind["beta"].values
    se_y = d_ind["se"].values
    maf = d_ind["maf"].values
    snps = d_ind["snp"].values

    n_iv = len(snps)
    print(f"\n[2] 多工具变量 MR（clump 后 {n_iv} 个独立 SNP，含 {IV}）")

    # IVW
    b_ivw, se_ivw = ivw(beta_y, se_y)
    z_ivw = b_ivw / se_ivw
    p_ivw = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(z_ivw)))
    print(f"    IVW:          β={b_ivw:.4f}, se={se_ivw:.4f}, z={z_ivw:.2f}, p={p_ivw:.2e}")

    # MR-Egger
    inter, se_inter, slope, se_slope = mr_egger(beta_x, se_x, beta_y, se_y)
    z_inter = inter / se_inter if se_inter and not np.isnan(se_inter) else np.nan
    p_inter = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(z_inter))) if not np.isnan(z_inter) else np.nan
    print(f"    MR-Egger:     截距={inter:.4f}, se={se_inter:.4f}, p(截距)={p_inter:.3f}")
    print(f"                 斜率={slope:.4f}, se={se_slope:.4f}")
    egger_interp = "无水平多效性证据" if not np.isnan(p_inter) and p_inter > 0.05 else ("有水平多效性信号" if not np.isnan(p_inter) else "不可估计")
    print(f"                 → {egger_interp}")

    # 加权中位数
    w = 1 / se_y**2
    b_wm = weighted_median(beta_y, se_y, weights=w)
    print(f"    加权中位数:   β={b_wm:.4f}")

    # ---- 3. 留一法 ----
    print(f"\n[3] 留一法（leave-one-out，{n_iv} 个 SNP）")
    loo_bs = []
    for i in range(n_iv):
        mask = np.ones(n_iv, dtype=bool)
        mask[i] = False
        b, _ = ivw(beta_y[mask], se_y[mask])
        loo_bs.append(b)
    loo_bs = np.array(loo_bs)
    # 找出剔除后偏离最大的 SNP
    b_all = b_ivw
    dev = np.abs(loo_bs - b_all)
    worst_i = np.argmax(dev)
    print(f"    IVW 全量 β={b_all:.4f}")
    print(f"    留一法 β 范围: [{loo_bs.min():.4f}, {loo_bs.max():.4f}]")
    print(f"    剔除 {snps[worst_i]} 后 β={loo_bs[worst_i]:.4f} (偏离最大)")
    is_iv_worst = snps[worst_i] == IV
    print(f"    是否工具变量 {IV} 驱动: {'是（单点驱动，需警惕）' if is_iv_worst else '否（稳健，非单点驱动）'}")

    all_rows.append(dict(
        outcome=oc, n_iv=n_iv,
        r2_exposure=r2_e, r2_outcome=r2_o, steiger_ok=steiger_ok,
        b_ivw=b_ivw, se_ivw=se_ivw, p_ivw=p_ivw,
        egger_intercept=inter, egger_intercept_p=p_inter,
        b_egger_slope=slope, b_wm=b_wm,
        loo_min=loo_bs.min(), loo_max=loo_bs.max(),
        loo_worst_snp=snps[worst_i], iv_driven=is_iv_worst,
    ))

res = pd.DataFrame(all_rows)
res.to_csv("results/SCRN2_sensitivity.csv", index=False, encoding="utf-8-sig")
print("\n\n已写入 results/SCRN2_sensitivity.csv")
print(res.to_string(index=False))
