# -*- coding: utf-8 -*-
"""
x02_transancestry.py
====================
Replaces "cross-ancestry replication" with two quantities that can go in an abstract:
 1) fixed-effect IVW pooling of the two cohorts (discovery = thyroid eQTL/EUR;
    check = South Asian, n = 38,000), with Cochran Q and I2 - "no detectable
    heterogeneity" is both more accurate and a stronger statement than "replication"
 2) the Bayesian posterior of directional agreement: how far the South Asian data
    push down the probability that the true effect points the other way
Output: x02_transancestry.txt / x02_transancestry.csv
"""
import numpy as np
from scipy import stats

# ---- inputs (from DATA15 mr_thyroid_lead.csv; DATA17 s03_replication_lead.csv) ----
disc = dict(name="Discovery (EUR, GTEx v8 thyroid)", b=0.0335876875001555, se=0.008851453010536043)
rep  = dict(name="Replication (South Asian, n=38,000)", b=0.028822465435781586, se=0.02413250467117654)
lead_rep = dict(b=-0.0115847, se=0.00969965)   # rs6503875 -> TSH, SAS
lead_dis = dict(b=-0.0135000, se=0.00320000)   # rs6503875 -> TSH, EUR

lines = []
P = lines.append
P("=" * 72)
P("Cross-ancestry concordance: fixed-effect pooling and a heterogeneity test")
P("=" * 72)
w1, w2 = 1/disc["se"]**2, 1/rep["se"]**2
b_pool = (w1*disc["b"] + w2*rep["b"]) / (w1+w2)
se_pool = 1/np.sqrt(w1+w2)
z = b_pool/se_pool
p = 2*stats.norm.sf(abs(z))
Q = w1*(disc["b"]-b_pool)**2 + w2*(rep["b"]-b_pool)**2
pQ = stats.chi2.sf(Q, df=1)
I2 = max(0, (Q-1)/Q*100) if Q > 0 else 0.0

P(f"Discovery set  b = {disc['b']:+.4f} (se {disc['se']:.4f})  z = {disc['b']/disc['se']:.2f}")
P(f"Check set      b = {rep['b']:+.4f} (se {rep['se']:.4f})  z = {rep['b']/rep['se']:.2f}  p = {2*stats.norm.sf(abs(rep['b']/rep['se'])):.3f}")
P("")
P(f"Fixed-effect pooled  b = {b_pool:+.4f} (se {se_pool:.4f})")
P(f"              95% CI [{b_pool-1.96*se_pool:+.4f}, {b_pool+1.96*se_pool:+.4f}]")
P(f"              z = {z:.2f}, p = {p:.2e}")
P(f"Cochran Q = {Q:.3f} (df = 1), p = {pQ:.3f}, I2 = {I2:.1f}%")
P(f"Ratio of the two cohort point estimates = {rep['b']/disc['b']:.3f}")
P("")
P("Wording of the conclusion: there is no significant heterogeneity between the two cohorts")
P("          (Q p = %.2f), and the pooled effect is almost identical to the discovery set." % pQ)
P("          This supports 'direction-concordant without detectable heterogeneity',")
P("          but it cannot be called 'replication' (the check set alone is not significant and its")
P("          95%% CI is as wide as +/-%.3f)." % (1.96*rep['se']))

# ---- Bayesian: strength of evidence for directional agreement ----
P("")
P("=" * 72)
P("Bayesian weight of the directional agreement")
P("=" * 72)
P("Model: the true effect is +d or -d (50% prior each), observation normal. d = discovery effect.")
for d, tag in [(disc["b"], "d = discovery effect (0.0336)"),
               (disc["b"]/2, "d = half the discovery effect (0.0168)"),
               (0.05, "d = 0.05 (a larger effect)")]:
    p_pos = stats.norm.sf(0, loc=d, scale=rep["se"])       # P(X>0 | true effect = +d)
    p_neg = stats.norm.sf(0, loc=-d, scale=rep["se"])      # P(X>0 | true effect = -d)
    LR = p_pos/p_neg
    post = LR/(1+LR)
    P(f"{tag}")
    P(f"  P(observed positive | true effect +d) = {p_pos:.3f};  P(observed positive | true effect -d) = {p_neg:.3f}")
    P(f"  likelihood ratio = {LR:.2f}  ->  given a positive observation, posterior probability that the")
    P(f"  true effect is positive = {post*100:.1f}%")
    P(f"  i.e. given the observation, the posterior probability that 'the true effect points the other")
    P(f"  way' = {(1-post)*100:.1f}% (prior 50%)")
    P("")

# ---- magnitude agreement at the lead SNP ----
ratio = lead_rep["b"]/lead_dis["b"]
se_ratio = abs(ratio)*np.sqrt((lead_rep["se"]/lead_rep["b"])**2 + (lead_dis["se"]/lead_dis["b"])**2)
P("=" * 72)
P("Lead-SNP level (rs6503875 -> TSH)")
P("=" * 72)
P(f"EUR b = {lead_dis['b']:+.5f} (se {lead_dis['se']:.5f})")
P(f"SAS b = {lead_rep['b']:+.5f} (se {lead_rep['se']:.5f})")
P(f"magnitude ratio = {ratio:.3f} (95% CI {ratio-1.96*se_ratio:.2f}-{ratio+1.96*se_ratio:.2f})  -- CI is wide, descriptive only")
wa, wb = 1/lead_dis["se"]**2, 1/lead_rep["se"]**2
bb = (wa*lead_dis["b"] + wb*lead_rep["b"])/(wa+wb)
sb = 1/np.sqrt(wa+wb)
P(f"  b = {bb:+.5f} (se {sb:.5f}), p = {2*stats.norm.sf(abs(bb/sb)):.2e}")

txt = "\n".join(lines)
open("x02_transancestry.txt", "w", encoding="utf-8").write(txt + "\n")
print(txt)

import pandas as pd
pd.DataFrame([
    dict(cohort="Discovery (EUR)", b=disc["b"], se=disc["se"], role="discovery"),
    dict(cohort="Replication (SAS)", b=rep["b"], se=rep["se"], role="replication"),
    dict(cohort="Pooled (fixed effect)", b=b_pool, se=se_pool, role="pooled"),
]).to_csv("x02_transancestry.csv", index=False)
