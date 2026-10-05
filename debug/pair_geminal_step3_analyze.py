r"""Combine LM blocks from debug/data/pairgem_<ref>_lm_s*.npz and solve the generalized eigenproblem
on every step-3 ansatz subspace (exact correlated sampling: one S,H for all).  Jackknife over 50
groups for E and for differences vs (A) and (B).  Writes optimum vectors (full K basis) to
debug/data/pairgem_<ref>_vecs.npz for the re-sampling VMC.
Usage: python debug/pair_geminal_step3_analyze.py 90det [tol]
"""
from __future__ import annotations

import os
import sys
import glob
import itertools
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ci_relax_step2_analyze import solve  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def main():
    ref = sys.argv[1] if len(sys.argv) > 1 else "90det"
    tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-10
    files = sorted(glob.glob(os.path.join(DATA, f"pairgem_{ref}_lm_s*.npz")))
    Sb, Hb = [], []
    for f in files:
        z = np.load(f, allow_pickle=True)
        Sb.append(z["S"]); Hb.append(z["H"]); cvec = z["cvec"]; N = int(z["N"])
        labels = list(z["labels"]); psi_cols = [tuple(x) for x in z["psi_cols"]]
        dress = [tuple(x) for x in z["dress"]]
    Sb = np.concatenate(Sb); Hb = np.concatenate(Hb); nb = len(Sb); K = Sb.shape[1]
    ng = min(50, nb); per = nb // ng
    Sg = Sb[:ng * per].reshape(ng, per, K, K).mean(1)
    Hg = Hb[:ng * per].reshape(ng, per, K, K).mean(1)
    S = Sg.mean(0); H = Hg.mean(0)
    print(f"ref={ref} files={len(files)} blocks={nb} groups={ng} K={K} N={N}")
    if len(dress) == 0 and len(psi_cols) == 0:          # lite run: step-2 ansatz set only
        from ci_relax_step2_analyze import projectors
        Qs = projectors(cvec)
        res = {k: solve(S, H, Qs[k], tol)[0] for k in Qs}
        jk = {k: [] for k in Qs}
        for g in range(ng):
            Sj = (S * ng - Sg[g]) / (ng - 1); Hj = (H * ng - Hg[g]) / (ng - 1)
            for k in Qs:
                jk[k].append(solve(Sj, Hj, Qs[k], tol)[0])
        jk = {k: np.array(v) for k, v in jk.items()}
        jerr = lambda a: np.sqrt((ng - 1) / ng * np.sum((a - a.mean()) ** 2))
        for k in Qs:
            print(f"  E({k}) = {res[k]:.6f} +/- {jerr(jk[k]):.6f}")
        for a, b in (("C", "0"), ("A", "C"), ("B", "C"), ("B", "A")):
            print(f"  E({a})-E({b}) = {(res[a]-res[b])*1e3:+.3f} +/- {jerr(jk[a]-jk[b])*1e3:.3f} mHa")
        np.savez(os.path.join(DATA, f"pairgem_{ref}_vecs.npz"),
                 **{f"v_{k}": solve(S, H, Qs[k], tol)[1] for k in Qs},
                 **{f"E_{k}": res[k] for k in Qs})
        return

    def e(i):
        v = np.zeros(K, complex); v[i] = 1.0; return v
    def blockvec(b):            # cvec placed on block b (block 0 = Phi, 1 = F1.5 Phi, 2.. dressed)
        v = np.zeros(K, complex); v[b * N:(b + 1) * N] = cvec; return v
    PHI = [e(i) for i in range(N)]
    F15 = [e(N + i) for i in range(N)]
    off = (2 + len(dress)) * N
    pidx = {k: off + j for j, k in enumerate(psi_cols)}
    dblk = {d: 2 + j for j, d in enumerate(dress)}
    def P(r0, c, g):
        return e(pidx[(r0, c, float(g))])
    def Fglob(g):               # global F_g Psi = sum of the classes (r0-independent)
        return P(1.0, "cc", g) + P(1.0, "cv", g) + P(1.0, "vv", g)
    def D(c, g):
        b = dblk[(c, g)]; return [e(b * N + i) for i in range(N)]
    def Q(cols):
        return np.column_stack(cols)

    A = PHI + [blockvec(1)]
    B = PHI + F15
    sub = {
        "C": Q([blockvec(0), blockvec(1)]),
        "A": Q(A), "B": Q(B),
    }
    # consistency: F1.5 Psi column vs the class-sum at 1.5 (should coincide in S,H up to noise)
    for r0 in (0.6, 1.0, 1.5):
        sub[f"A_split1.5_r{r0}"] = Q(PHI + [P(r0, c, 1.5) for c in ("cc", "cv", "vv")])
    for g in (1.0, 2.0, 3.0, 4.0, 5.0):
        sub[f"A3_glob{g}"] = Q(A + [Fglob(g)])
        for r0 in (0.6, 1.0, 1.5):
            sub[f"A2cc_r{r0}_g{g}"] = Q(A + [P(r0, "cc", g)])
    # full per-class scan (A2): Phi + one geminal per class with its own gamma
    scan = {}
    for r0 in (0.6, 1.0, 1.5):
        for gc in (2.0, 3.0, 4.0, 5.0):
            for gcv in (1.0, 1.5, 2.0, 3.0):
                for gv in (1.0, 1.5):
                    E, _, _ = solve(S, H, Q(PHI + [P(r0, "cc", gc), P(r0, "cv", gcv), P(r0, "vv", gv)]), tol)
                    scan[(r0, gc, gcv, gv)] = E
    for r0 in (0.6, 1.0, 1.5):
        best = min((k for k in scan if k[0] == r0), key=lambda k: scan[k])
        _, gc, gcv, gv = best
        sub[f"A2best_r{r0}(gc{gc},gcv{gcv},gv{gv})"] = Q(PHI + [P(r0, "cc", gc), P(r0, "cv", gcv),
                                                                 P(r0, "vv", gv)])
        sub[f"A+A2best_r{r0}"] = Q(A + [P(r0, "cc", gc), P(r0, "cv", gcv), P(r0, "vv", gv)])
        sub[f"A_allclass_r{r0}"] = Q(PHI + [P(r0, c, g) for c in ("cc", "cv", "vv")
                                            for g in (1.0, 1.5, 2.0, 3.0, 4.0, 5.0)])
    sub["A_glob_all"] = Q(PHI + [Fglob(g) for g in (1.0, 1.5, 2.0, 3.0, 4.0, 5.0)])
    sub["A_allP"] = Q(PHI + [e(i) for i in range(off, K)])
    # dressed (B-type)
    sub["B2_cc3"] = Q(B + D("cc", 3.0))
    sub["B2_cc5"] = Q(B + D("cc", 5.0))
    sub["B3_glob3"] = Q(B + [a + b for a, b in zip(D("cc", 3.0), D("ncc", 3.0))])
    sub["B4_split3"] = Q(B + D("cc", 3.0) + D("ncc", 3.0))
    sub["B5_cc3cc5"] = Q(B + D("cc", 3.0) + D("cc", 5.0))
    sub["B_all_dressed"] = Q(B + D("cc", 3.0) + D("cc", 5.0) + D("ncc", 3.0))
    sub["B+A_allclass_r1.0"] = Q(B + [P(1.0, c, g) for c in ("cc", "cv", "vv")
                                      for g in (1.0, 1.5, 2.0, 3.0, 4.0, 5.0)])
    sub["B+glob3Psi"] = Q(B + [Fglob(3.0)])
    sub["B+cc3Psi_r1.0"] = Q(B + [P(1.0, "cc", 3.0)])
    sub["EVERYTHING"] = np.eye(K, dtype=complex)

    res = {}; vecs = {}
    for k, Qk in sub.items():
        E, v, nk = solve(S, H, Qk, tol)
        res[k] = E; vecs[k] = v
    # jackknife
    jk = {k: [] for k in sub}
    for g in range(ng):
        Sj = (S * ng - Sg[g]) / (ng - 1); Hj = (H * ng - Hg[g]) / (ng - 1)
        for k, Qk in sub.items():
            jk[k].append(solve(Sj, Hj, Qk, tol)[0])
    jk = {k: np.array(v) for k, v in jk.items()}
    def jerr(a):
        return np.sqrt((ng - 1) / ng * np.sum((a - a.mean()) ** 2))
    print(f"\n{'ansatz':42s} {'E':>11s} {'err':>8s}  {'-E(A) mHa':>16s}  {'-E(B) mHa':>16s}  nparam")
    for k in sub:
        dA = (res[k] - res["A"]) * 1e3; eA = jerr(jk[k] - jk["A"]) * 1e3
        dB = (res[k] - res["B"]) * 1e3; eB = jerr(jk[k] - jk["B"]) * 1e3
        print(f"{k:42s} {res[k]:11.6f} {jerr(jk[k]):8.6f}  {dA:+7.3f}+/-{eA:5.3f}  {dB:+7.3f}+/-{eB:5.3f}  {sub[k].shape[1]}")
    print("\nA2 scan (E rel. to A, mHa, no jackknife), sorted, top 12 per radius:")
    for r0 in (0.6, 1.0, 1.5):
        ks = sorted((k for k in scan if k[0] == r0), key=lambda k: scan[k])
        print(f"  r0={r0}: " + "; ".join(f"(gc{k[1]:g},gcv{k[2]:g},gv{k[3]:g}) {1e3*(scan[k]-res['A']):+.2f}"
                                        for k in ks[:12]))
        print(f"    worst: {1e3*(scan[ks[-1]]-res['A']):+.2f}")
    out = {f"v_{k}": vecs[k] for k in vecs}
    out.update({f"E_{k}": res[k] for k in res}); out.update({f"err_{k}": jerr(jk[k]) for k in res})
    np.savez(os.path.join(DATA, f"pairgem_{ref}_vecs.npz"), **out)


if __name__ == "__main__":
    main()
