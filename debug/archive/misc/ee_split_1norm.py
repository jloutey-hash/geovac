"""THE DECIDING MEASUREMENT: does the e-e split lower the LCU 1-norm, and does it
beat plain density fitting at matched accuracy?

Context.  The split g = g_sep - W has one property generic low-rank factorization does
NOT: g_sep is EXACTLY a one-body operator on the N-electron space (bit-exact, verified),
so folding it into h1 moves weight from the expensive two-body block to the cheap
one-body block BY CONSTRUCTION rather than by approximation.  Whether that helps is a
measurement, and the corpus's prior is NEGATIVE: "analytic momentum-space factorization
as a 1-norm lever" (2026-08-21) found the lambda wins were the generic density-fitting
advantage, matched by eigen-Cholesky DF within 1-13%.  Also measured: ||W||_F = 1.08
||g||_F, so the remainder is NOT smaller than the original in Frobenius norm.

Four variants, same system, same JW convention (tracked geovac lcu_lambda):
  A  baseline           h1,      g                       (exact)
  B  split, exact       h_eff,   W          (h_eff = h1 + (N-1)/2 V)   (exact)
  C  split, truncated   h_eff,   W at rank m                (approximate)
  D  CONTROL: plain DF  h1,      g at rank m                (approximate)

DECISION RULE (pre-registered).
  If B < A materially            -> the exact one-body fold is a real 1-norm lever.
  If C beats D at MATCHED accuracy -> the split adds something beyond density fitting.
  If C ~ D                        -> it IS density fitting; corpus prior confirmed;
                                     the arc is structural/pedagogical only.
"""
import importlib.util
import io
import json
import os
import sys

import numpy as np
from scipy.linalg import eigh

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, ROOT)

from geovac import transcorrelated_sturmian as TC  # noqa: E402


def build(ns, k=2.0, Z=2.0, Ng=700):
    r, wr = TC.make_grid(k, Ng=Ng)
    S, h1, Rtab, W2 = TC.build_one_body(ns, r, wr, k, Z)
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    D = {(i, j): Rtab[i][0] * Rtab[j][0] * W2
         for i in range(1, ns + 1) for j in range(1, ns + 1)}

    def eri(K):
        g = np.zeros((ns,) * 4)
        for a, i in enumerate(range(1, ns + 1)):
            for b, j in enumerate(range(1, ns + 1)):
                for c, kk in enumerate(range(1, ns + 1)):
                    for d, ll in enumerate(range(1, ns + 1)):
                        g[a, b, c, d] = D[(i, kk)] @ K @ D[(j, ll)]
        return g

    g = eri(1.0 / np.maximum(R1, R2))
    gs = eri(0.5 * (1.0 / R1 + 1.0 / R2))
    gw = eri(0.5 * np.abs(1.0 / R1 - 1.0 / R2))
    V = np.array([[np.sum(Rtab[i][0] * Rtab[j][0] * r * wr) for j in range(1, ns + 1)]
                  for i in range(1, ns + 1)])
    return S, h1, V, g, gs, gw


def pairmat(t, n):
    return t.transpose(0, 2, 1, 3).reshape(n * n, n * n)


def unpair(M, n):
    return M.reshape(n, n, n, n).transpose(0, 2, 1, 3)


def sym4(t):
    return 0.25 * (t + t.transpose(2, 1, 0, 3) + t.transpose(0, 3, 2, 1)
                   + t.transpose(2, 3, 0, 1))


def trunc(t, n, m):
    lam, U = np.linalg.eigh(pairmat(t, n))
    o = np.argsort(-np.abs(lam))[:m]
    return sym4(unpair((U[:, o] * lam[o]) @ U[:, o].T, n))


def lam_and_E(S, h1x, gx, n_elec, ns):
    X = TC.lowdin(S)
    h1o, go = TC.transform_1(h1x, X), TC.transform_2(gx, X)
    nso = 2 * ns
    hso, asym = TC.h_spin(h1o, nso), TC.asym_from_phys(go, nso)
    lam = TC.lcu_lambda(hso, asym, nso)
    dets, didx = TC.make_dets(nso, n_elec)
    E = float(eigh(TC.build_H(dets, didx, hso, asym, nso), eigvals_only=True)[0])
    return lam, E


out = {}
N_ELEC = 2
print("=" * 92)
print("A/B -- does folding the EXACT one-body head lower the 1-norm?  (exact variants)")
print("=" * 92)
print(f"{'ns':>3}{'lam_A (full)':>14}{'lam_B (split)':>15}{'ratio B/A':>11}"
      f"{'E_A':>14}{'E_B':>14}{'dE (Ha)':>11}")
sysd = {}
for ns in (3, 4, 5, 6):
    S, h1, V, g, gs, gw = build(ns)
    sysd[ns] = (S, h1, V, g, gs, gw)
    h_eff = h1 + 0.5 * (N_ELEC - 1) * V
    lamA, EA = lam_and_E(S, h1, g, N_ELEC, ns)
    lamB, EB = lam_and_E(S, h_eff, -gw, N_ELEC, ns)
    kA, kB = lamA["lam"], lamB["lam"]
    print(f"{ns:>3}{kA:>14.4f}{kB:>15.4f}{kB/kA:>11.3f}{EA:>14.8f}{EB:>14.8f}"
          f"{abs(EA-EB):>11.1e}")
    out.setdefault("AB", {})[str(ns)] = dict(lamA=kA, lamB=kB, ratio=kB / kA,
                                             EA=EA, EB=EB, dE=abs(EA - EB))

print()
print("=" * 92)
print("C vs D -- split-truncated vs PLAIN density fitting, at matched rank")
print("=" * 92)
print(f"{'ns':>3}{'rank':>6}{'lam_C':>11}{'err_C (mHa)':>13}{'lam_D':>11}"
      f"{'err_D (mHa)':>13}{'verdict':>26}")
for ns in (4, 6):
    S, h1, V, g, gs, gw = sysd[ns]
    h_eff = h1 + 0.5 * (N_ELEC - 1) * V
    _, E_ex = lam_and_E(S, h1, g, N_ELEC, ns)
    for m in (2, 3, 4, 6, 8):
        lamC, EC = lam_and_E(S, h_eff, -trunc(gw, ns, m), N_ELEC, ns)
        lamD, ED = lam_and_E(S, h1, trunc(g, ns, m), N_ELEC, ns)
        kC, kD = lamC["lam"], lamD["lam"]
        eC, eD = 1000 * (EC - E_ex), 1000 * (ED - E_ex)
        # which is better: lower lambda AND lower |error|?
        if abs(eC) < abs(eD) and kC < kD:
            v = "C wins both"
        elif abs(eC) < abs(eD):
            v = "C: better acc, worse lam"
        elif kC < kD:
            v = "C: better lam, worse acc"
        else:
            v = "D wins both"
        print(f"{ns:>3}{m:>6}{kC:>11.4f}{eC:>13.4f}{kD:>11.4f}{eD:>13.4f}{v:>26}")
        out.setdefault("CD", {}).setdefault(str(ns), []).append(
            dict(m=m, lamC=kC, errC=eC, lamD=kD, errD=eD, verdict=v))

print()
print("=" * 92)
print("norm diagnostics (why): entrywise 1-norms of the tensors themselves")
print("=" * 92)
print(f"{'ns':>3}{'sum|g|':>12}{'sum|W|':>12}{'ratio':>9}{'||g||_F':>11}"
      f"{'||W||_F':>11}{'ratio':>9}")
for ns in (3, 4, 5, 6):
    S, h1, V, g, gs, gw = sysd[ns]
    s1, s2 = np.abs(g).sum(), np.abs(gw).sum()
    f1, f2 = np.linalg.norm(g), np.linalg.norm(gw)
    print(f"{ns:>3}{s1:>12.4f}{s2:>12.4f}{s2/s1:>9.3f}{f1:>11.4f}{f2:>11.4f}{f2/f1:>9.3f}")
    out.setdefault("norms", {})[str(ns)] = dict(sum_g=float(s1), sum_W=float(s2),
                                                ratio_1=float(s2 / s1),
                                                fro_g=float(f1), fro_W=float(f2),
                                                ratio_F=float(f2 / f1))

os.makedirs("debug/data", exist_ok=True)
with io.open("debug/data/ee_split_1norm.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2, default=float)
print("\nwrote debug/data/ee_split_1norm.json")
