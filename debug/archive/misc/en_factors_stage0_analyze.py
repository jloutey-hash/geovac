r"""Stage-0 analysis: LM energies with the parameter-free core switch vs the step-3 radial switches.

Reads debug/data/pairgem_<tag>_lm_s*.npz (tag en0 = 90det, en09det = 9det) from
debug/en_factors_stage0_lm.py.  Correlated sub-ansatzes of one sampled (S, H):
    C, A, B                              (step-2 baselines)
    A+W   = A + {G_cc, G_cv, G_vv}(gamma 0)  == A + {W1 psi, W2 psi}            (phase 1)
    A+G   = A + class geminals (gamma 1.5)   == A + {F_cc psi, F_c1 psi}
    A+W+G = both                                                           (phase 2)
for every switch key (0.6, 1.0 radial; -1 = DENS ruled; -2 = FRAC variant).  50-group jackknife.
--project : also print the LM matrices projected onto the ANALYTIC basis of en_factors
            [chi_I, xi, W1 psi, W2 psi, F_cc psi, F_c1 psi] (DENS switch; electronic H) and save them
            to debug/data/en_factors_stage0_<tag>_proj.npz for gate G2/G3 comparison.
Option b (PI 2026-10-09; tags enab = 90det, enab9det = 9det, keys 0.6, -1 DENS, -3 AMP):
    A+W[-1]+W[-3]  = A + {W1, W2}(dens) + {A1, A2}(amp)   (the combined phase 1)
    A+W+G[-1]+W[-3] etc.;  --project then projects onto [chi_I, xi, W1, W2, A1, A2].
Usage: python debug/en_factors_stage0_analyze.py en0 [--project]
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ci_relax_step2_analyze import solve  # noqa: E402

DATA = os.path.join(HERE, "data")
V_NN = 3.0 / 3.015
FBAR = {"en0": None, "en09det": None}


def load(tag):
    files = sorted(glob.glob(os.path.join(DATA, f"pairgem_{tag}_lm_s*.npz")))
    Sb, Hb = [], []
    for f in files:
        z = np.load(f, allow_pickle=True)
        Sb.append(z["S"]); Hb.append(z["H"]); cvec = z["cvec"]; N = int(z["N"])
        pc = [tuple(x) for x in z["psi_cols"]]
    return files, np.concatenate(Sb), np.concatenate(Hb), cvec, N, pc


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "en0"
    files, Sb, Hb, cvec, N, pc = load(tag)
    nb, K = len(Sb), Sb.shape[1]
    ng = min(50, nb); per = nb // ng
    Sg = Sb[:ng * per].reshape(ng, per, K, K).mean(1); Hg = Hb[:ng * per].reshape(ng, per, K, K).mean(1)
    S = Sg.mean(0); H = Hg.mean(0)
    print(f"tag={tag} files={len(files)} blocks={nb} K={K} N={N}")
    I = np.eye(K, dtype=complex); e = lambda i: I[:, i]

    def bv(b):
        v = np.zeros(K, complex); v[b * N:(b + 1) * N] = cvec; return v
    pidx = {(float(k[0]), k[1], float(k[2])): 2 * N + j for j, k in enumerate(pc)}
    P = lambda r0, c, g: e(pidx[(float(r0), c, float(g))])
    PHI = [e(i) for i in range(N)]; F15 = [e(N + i) for i in range(N)]
    A = PHI + [bv(1)]; B = PHI + F15
    Q = np.column_stack
    sub = {"C": Q([bv(0), bv(1)]), "A": Q(A), "B": Q(B)}
    keys = sorted(set(k[0] for k in pidx), key=lambda x: (x < 0, abs(x)))
    for r0 in keys:
        W = [P(r0, c, 0.0) for c in ("cc", "cv", "vv")]
        G = [P(r0, c, 1.5) for c in ("cc", "cv", "vv")]
        sub[f"A+W[{r0:g}]"] = Q(A + W)
        sub[f"A+G[{r0:g}]"] = Q(A + G)
        sub[f"A+W+G[{r0:g}]"] = Q(A + W + G)
        sub[f"C+W[{r0:g}]"] = Q([bv(0), bv(1)] + W)
    if -1.0 in keys and -3.0 in keys:                       # option b: density + amplitude switch
        Wd = [P(-1.0, c, 0.0) for c in ("cc", "cv", "vv")]; Wa = [P(-3.0, c, 0.0) for c in ("cc", "cv", "vv")]
        Gd = [P(-1.0, c, 1.5) for c in ("cc", "cv", "vv")]; Ga = [P(-3.0, c, 1.5) for c in ("cc", "cv", "vv")]
        sub["A+W[-1]+W[-3]"] = Q(A + Wd + Wa)
        sub["A+W+G[-1]+W[-3]"] = Q(A + Wd + Wa + Gd)
        sub["A+W+G[-1]+W+G[-3]"] = Q(A + Wd + Wa + Gd + Ga)
        sub["C+W[-1]+W[-3]"] = Q([bv(0), bv(1)] + Wd + Wa)
    res = {k: solve(S, H, q, 1e-10)[0] for k, q in sub.items()}
    jk = {k: [] for k in sub}
    for g in range(ng):
        Sj = (S * ng - Sg[g]) / (ng - 1); Hj = (H * ng - Hg[g]) / (ng - 1)
        for k, q in sub.items():
            jk[k].append(solve(Sj, Hj, q, 1e-10)[0])
    jk = {k: np.array(v) for k, v in jk.items()}
    je = lambda a: np.sqrt((ng - 1) / ng * np.sum((a - a.mean()) ** 2))
    print(f"{'ansatz':22s} {'E':>11s} {'err':>9s}  {'vs A (mHa)':>18s} {'vs B (mHa)':>18s} np")
    for k in sub:
        dA = (res[k] - res["A"]) * 1e3; eA = je(jk[k] - jk["A"]) * 1e3
        dB = (res[k] - res["B"]) * 1e3; eB = je(jk[k] - jk["B"]) * 1e3
        print(f"{k:22s} {res[k]:11.6f} {je(jk[k]):9.6f}  {dA:+8.3f} +/- {eA:5.3f}  {dB:+8.3f} +/- {eB:5.3f} {sub[k].shape[1]}")
    np.savez(os.path.join(DATA, f"en_factors_stage0_{tag}_E.npz"),
             names=np.array(list(sub)), E=np.array([res[k] for k in sub]),
             jk=np.array([jk[k] for k in sub]))
    if "--project" in sys.argv:
        # analytic basis [chi_I (N), xi, W1 psi, W2 psi, Fcc psi, Fc1 psi]; chi = Phi/2, psi = Psi_LM/2
        import json
        gold = json.load(open(os.path.join(os.path.dirname(HERE), "tests", "data", "lih_golden", "golden.json"),
                              encoding="utf-8"))
        items = gold.get("items", gold)
        it = items["pi_marriage_exp1.5_1e-3" if N > 20 else "pi_marriage_exp1.5_1e-2"]["values"]
        Fbar = float(it["pieces"]["value"]["Fbar"]["repr"])
        k = -1.0
        cols = [0.5 * e(i) for i in range(N)]
        cols.append(0.5 * (bv(1) - Fbar * bv(0)))
        cols.append(0.5 * (P(k, "cv", 0.0) + 2 * P(k, "cc", 0.0)) / 3.0)
        cols.append(0.5 * P(k, "cc", 0.0))
        if -3.0 in keys:                                    # option b basis: + A1 psi, A2 psi
            cols.append(0.5 * (P(-3.0, "cv", 0.0) + 2 * P(-3.0, "cc", 0.0)) / 3.0)
            cols.append(0.5 * P(-3.0, "cc", 0.0))
        else:
            cols.append(0.5 * P(k, "cc", 1.5))
            cols.append(0.5 * (P(k, "cv", 1.5) + 2 * P(k, "cc", 1.5)))
        Qa = np.column_stack(cols)
        Sp = np.array([Qa.conj().T @ s @ Qa for s in Sg]).real
        Hp = np.array([Qa.conj().T @ h @ Qa for h in Hg]).real - V_NN * Sp
        np.savez(os.path.join(DATA, f"en_factors_stage0_{tag}_proj.npz"), S=Sp.mean(0), H=Hp.mean(0),
                 S_err=Sp.std(0, ddof=1) / np.sqrt(ng), H_err=Hp.std(0, ddof=1) / np.sqrt(ng), Fbar=Fbar)
        print(f"projected analytic-basis matrices saved (Fbar={Fbar:.6f})")


if __name__ == "__main__":
    main()
