r"""Analyse debug/data/mg_lm_s*.npz (multi-geminal linear method) + mg_pade_*.npz.

Pooled jackknife (blocks within each seed, all seeds pooled) on the generalized eigenproblem
H c = E S c for every subset.  Per-seed energies + chi2/dof as a consistency check.
"""
from __future__ import annotations

import os
import sys
import glob
from itertools import combinations

import numpy as np
import scipy.linalg as sl

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
E_ANALYTIC_K1 = -8.042181268689944
C_ANALYTIC = 0.22158757687188235
FBAR = 3.497901461468177
NBLK = 25   # blocks per seed


def lowest(S, H, sub):
    s = S[np.ix_(sub, sub)]; h = H[np.ix_(sub, sub)]
    w, v = sl.eig(h, s)
    w = w.real
    k = int(np.argmin(np.where(np.isfinite(w), w, np.inf)))
    return w[k], v[:, k].real


def jack(Sblk, Hblk, sub):
    """Sblk,Hblk: (nb, n, n) block means (equal weight). Returns E, err, and per-block
    leave-one-out energies (for correlated differences)."""
    nb = len(Sblk)
    S = Sblk.mean(0); H = Hblk.mean(0)
    E, v = lowest(S, H, sub)
    tot_S = Sblk.sum(0); tot_H = Hblk.sum(0)
    Ej = np.array([lowest((tot_S - Sblk[b]) / (nb - 1), (tot_H - Hblk[b]) / (nb - 1), sub)[0]
                   for b in range(nb)])
    err = np.sqrt((nb - 1) / nb * np.sum((Ej - Ej.mean()) ** 2))
    return E, err, Ej, v


def jdiff(Ej1, Ej2):
    d = Ej1 - Ej2; nb = len(d)
    return np.sqrt((nb - 1) / nb * np.sum((d - d.mean()) ** 2))


def main() -> None:
    files = sorted(glob.glob(os.path.join(DATA, "mg_lm_s*.npz")))
    Sb, Hb, perseed = [], [], []
    for f in files:
        z = np.load(f)
        S = z['S']; H = z['H']; labels = list(z['labels'])
        n = len(S) // NBLK * NBLK
        Sbl = S[:n].reshape(NBLK, -1, *S.shape[1:]).mean(1)
        Hbl = H[:n].reshape(NBLK, -1, *H.shape[1:]).mean(1)
        Sb.append(Sbl); Hb.append(Hbl); perseed.append((os.path.basename(f), Sbl, Hbl, len(S)))
    Sb = np.concatenate(Sb); Hb = np.concatenate(Hb)
    print(f"files: {len(files)}  total blocks {len(Sb)}  labels {labels}")
    iLX = labels.index('linexp0.5'); iF2 = labels.index('F2(linexp0.5)')
    gem = [i for i in range(1, len(labels)) if i != iF2]

    # ---- V0
    E0, e0, Ej0, _ = jack(Sb, Hb, [0])
    E1, e1, Ej1, v1 = jack(Sb, Hb, [0, iLX])
    c = v1[1] / (v1[0] + v1[1] * FBAR)
    # jackknife c
    tS = Sb.sum(0); tH = Hb.sum(0); nb = len(Sb); cj = []
    for b in range(nb):
        _, vv = lowest((tS - Sb[b]) / (nb - 1), (tH - Hb[b]) / (nb - 1), [0, iLX])
        cj.append(vv[1] / (vv[0] + vv[1] * FBAR))
    cj = np.array(cj); ce = np.sqrt((nb - 1) / nb * np.sum((cj - cj.mean()) ** 2))
    print(f"\nV0: E(Psi0) = {E0:.6f} +/- {e0:.6f}  (E_trunc -8.025184)")
    print(f"V0: E_k1(linexp0.5) = {E1:.6f} +/- {e1:.6f}  analytic {E_ANALYTIC_K1:.6f}  "
          f"d = {(E1-E_ANALYTIC_K1)*1e3:+.3f} mHa = {(E1-E_ANALYTIC_K1)/e1:+.2f} sigma")
    print(f"V0: c = {c:+.5f} +/- {ce:.5f}  analytic {C_ANALYTIC:+.5f}")
    print("per-seed k1:")
    Es, es = [], []
    for name, Sbl, Hbl, ns in perseed:
        E, e, _, _ = jack(Sbl, Hbl, [0, iLX]); Es.append(E); es.append(e)
        print(f"   {name}: {E:.6f} +/- {e:.6f}  (nsamp {ns})")
    Es = np.array(Es); es = np.array(es); w = 1 / es ** 2
    Ew = np.sum(w * Es) / w.sum()
    chi2 = np.sum(((Es - Ew) / es) ** 2)
    print(f"   inv-var {Ew:.6f} +/- {1/np.sqrt(w.sum()):.6f}  chi2/dof = {chi2:.2f}/{len(Es)-1}")

    # ---- subset scan
    print("\nBest subsets per k (geminal pool, 8 functions):")
    results = {}
    for k in range(1, 5):
        best = None
        for comb in combinations(gem, k):
            sub = [0] + list(comb)
            E, _ = lowest(Sb.mean(0), Hb.mean(0), sub)
            if best is None or E < best[0]:
                best = (E, sub)
        E, e, Ej, _ = jack(Sb, Hb, best[1])
        results[k] = (E, e, Ej, best[1])
        print(f"  k={k}: {[labels[i] for i in best[1][1:]]}  E = {E:.6f} +/- {e:.6f}  "
              f"dE vs k1(lx0.5) = {(E-E1)*1e3:+.3f} +/- {jdiff(Ej, Ej1)*1e3:.3f} mHa")
    # nested sequences
    print("\nNested sequences:")
    seqs = {
        "linexp ladder": [labels.index(x) for x in ('linexp0.5', 'linexp2.0', 'linexp0.25', 'linexp1.0')],
        "mixed": [labels.index(x) for x in ('linexp0.5', 'exp2.0', 'linexp0.25', 'exp1.0')],
    }
    for nm, seq in seqs.items():
        for k in range(1, 5):
            E, e, Ej, _ = jack(Sb, Hb, [0] + seq[:k])
            print(f"  {nm} k={k} {[labels[i] for i in seq[:k]]}: E = {E:.6f} +/- {e:.6f}  "
                  f"dE = {(E-E1)*1e3:+.3f} +/- {jdiff(Ej, Ej1)*1e3:.3f} mHa")
    # all 8
    E8, e8, Ej8, _ = jack(Sb, Hb, [0] + gem)
    print(f"  all 8 geminals: E = {E8:.6f} +/- {e8:.6f}  dE = {(E8-E1)*1e3:+.3f} +/- "
          f"{jdiff(Ej8, Ej1)*1e3:.3f} mHa")
    # quadratic probe
    Eq, eq, Ejq, _ = jack(Sb, Hb, [0, iLX, iF2])
    print(f"\nQuadratic probe {{Psi0, F, F^2}} (linexp0.5): E = {Eq:.6f} +/- {eq:.6f}  dE = "
          f"{(Eq-E1)*1e3:+.3f} +/- {jdiff(Ejq, Ej1)*1e3:.3f} mHa")
    Eq8, eq8, Ejq8, _ = jack(Sb, Hb, [0] + gem + [iF2])
    print(f"All 8 + F^2: E = {Eq8:.6f} +/- {eq8:.6f}  dE = {(Eq8-E1)*1e3:+.3f} +/- "
          f"{jdiff(Ejq8, Ej1)*1e3:.3f} mHa")

    # ---- Pade
    pf = sorted(glob.glob(os.path.join(DATA, "mg_pade_*.npz")))
    byg = {}
    for f in pf:
        z = np.load(f); byg.setdefault(float(z['gamma']), []).append((float(z['E']), float(z['err'])))
    print("\nPade (same trunc-1e-2 reference):")
    bestP = None
    for g in sorted(byg):
        a = np.array(byg[g]); w = 1 / a[:, 1] ** 2
        Ew = np.sum(w * a[:, 0]) / w.sum(); ew = 1 / np.sqrt(w.sum())
        chi2 = np.sum(((a[:, 0] - Ew) / a[:, 1]) ** 2)
        print(f"  gamma={g}: {Ew:.6f} +/- {ew:.6f}  (n={len(a)}, chi2/dof {chi2:.2f}/{len(a)-1})")
        if bestP is None or Ew < bestP[0]:
            bestP = (Ew, ew, g)
    if bestP:
        EP, eP, gP = bestP
        gap = E1 - EP
        print(f"\nGap k1 -> Pade(gamma={gP}): {gap*1e3:.3f} +/- {np.hypot(e1, eP)*1e3:.3f} mHa")
        for k in range(1, 5):
            E, e, Ej, _ = results[k]
            dk = E1 - E; edk = jdiff(Ej, Ej1)
            frac = dk / gap
            efrac = abs(frac) * np.sqrt((edk / max(abs(dk), 1e-12)) ** 2 + (np.hypot(e1, eP) / gap) ** 2)
            print(f"  k={k}: closed {100*frac:.1f} +/- {100*efrac:.1f} %")
        for nm, Ex, Ejx in (("F,F^2", Eq, Ejq), ("all8", E8, Ej8), ("all8+F^2", Eq8, Ejq8)):
            print(f"  {nm}: closed {100*(E1-Ex)/gap:.1f} %")


if __name__ == "__main__":
    main()
