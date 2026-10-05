r"""Discriminator for the step-3 pair-split gain: is it r12 (cusp/pair) content or a one-body /
radial reshaping that needs no r12 at all?  gamma = 0 class columns are r12-FREE:
G_cc = sum_{i<j} w_i w_j, G_cv, G_vv span {1, e1 = sum_i w_i, e2 = sum_{i<j} w_i w_j} x Psi_CI
(a linear electron-nucleus "Jastrow" in the core switch).  Same correlated sampling as the
gamma = 1, 2 columns.  Usage: python debug/pair_geminal_step3_analyze_g0.py
"""
from __future__ import annotations

import os
import sys
import glob
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ci_relax_step2_analyze import solve  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def main():
    files = sorted(glob.glob(os.path.join(DATA, "pairgem_90detg0_lm_s*.npz")))
    Sb, Hb = [], []
    for f in files:
        z = np.load(f, allow_pickle=True)
        Sb.append(z["S"]); Hb.append(z["H"]); cvec = z["cvec"]; N = int(z["N"])
        psi_cols = [tuple(x) for x in z["psi_cols"]]
    Sb = np.concatenate(Sb); Hb = np.concatenate(Hb); nb = len(Sb); K = Sb.shape[1]
    ng = min(50, nb); per = nb // ng
    Sg = Sb[:ng * per].reshape(ng, per, K, K).mean(1)
    Hg = Hb[:ng * per].reshape(ng, per, K, K).mean(1)
    S = Sg.mean(0); H = Hg.mean(0)
    print(f"files={len(files)} blocks={nb} K={K} N={N}")
    e = lambda i: np.eye(K, dtype=complex)[:, i]
    def bv(b):
        v = np.zeros(K, complex); v[b * N:(b + 1) * N] = cvec; return v
    PHI = [e(i) for i in range(N)]; F15 = [e(N + i) for i in range(N)]
    pidx = {k: 2 * N + j for j, k in enumerate(psi_cols)}
    P = lambda r0, c, g: e(pidx[(r0, c, float(g))])
    A = PHI + [bv(1)]; B = PHI + F15
    Q = np.column_stack
    sub = {"C": Q([bv(0), bv(1)]), "A": Q(A), "B": Q(B)}
    for r0 in (0.6, 1.0):
        for g in (0.0, 1.0, 2.0):
            sub[f"A+cc_g{g:g}_r{r0}"] = Q(A + [P(r0, "cc", g)])
            sub[f"A+3cls_g{g:g}_r{r0}"] = Q(A + [P(r0, c, g) for c in ("cc", "cv", "vv")])
        sub[f"A+3cls_g0+g1_r{r0}"] = Q(A + [P(r0, c, g) for c in ("cc", "cv", "vv") for g in (0.0, 1.0)])
        sub[f"A+3cls_g1+g2_r{r0}"] = Q(A + [P(r0, c, g) for c in ("cc", "cv", "vv") for g in (1.0, 2.0)])
        sub[f"A+3cls_all_r{r0}"] = Q(A + [P(r0, c, g) for c in ("cc", "cv", "vv") for g in (0.0, 1.0, 2.0)])
    sub["A+g0_both_radii"] = Q(A + [P(r0, c, 0.0) for r0 in (0.6, 1.0) for c in ("cc", "cv", "vv")])
    sub["A+allcols"] = Q(A + [e(i) for i in range(2 * N, K)])
    sub["B+g0_both_radii"] = Q(B + [P(r0, c, 0.0) for r0 in (0.6, 1.0) for c in ("cc", "cv", "vv")])
    sub["B+allcols"] = np.eye(K, dtype=complex)
    res = {k: solve(S, H, q, 1e-10)[0] for k, q in sub.items()}
    jk = {k: [] for k in sub}
    for g in range(ng):
        Sj = (S * ng - Sg[g]) / (ng - 1); Hj = (H * ng - Hg[g]) / (ng - 1)
        for k, q in sub.items():
            jk[k].append(solve(Sj, Hj, q, 1e-10)[0])
    jk = {k: np.array(v) for k, v in jk.items()}
    jerr = lambda a: np.sqrt((ng - 1) / ng * np.sum((a - a.mean()) ** 2))
    print(f"{'ansatz':26s} {'E':>11s} {'err':>8s}  {'-E(A) mHa':>16s}  nparam")
    for k in sub:
        d = (res[k] - res["A"]) * 1e3; ed = jerr(jk[k] - jk["A"]) * 1e3
        print(f"{k:26s} {res[k]:11.6f} {jerr(jk[k]):8.6f}  {d:+7.3f}+/-{ed:5.3f}  {sub[k].shape[1]}")


if __name__ == "__main__":
    main()
