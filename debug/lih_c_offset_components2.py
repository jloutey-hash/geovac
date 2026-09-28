r"""Analyse debug/data/cvmc2_s*.npz (lih_c_offset_vmc2.py): per-component VMC for linexp0.25 and
linexp0.5 on the trunc-1e-2 pi reference, with the kinetic estimated TWO independent ways
(non-symmetric -1/2 lap, symmetric 1/2|grad|^2), vs the analytic sparse pieces.  Block jackknife.
"""
from __future__ import annotations

import os
import glob
import numpy as np
import scipy.linalg as sl

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
NBLK = 25
VNN = 0.9950248756218905
AN_PAT = {1: "lih_exp_geminal_linexp_g0.25_72x44x28_p28-48_1e-02.npz",
          2: "lih_exp_geminal_linexp_g0.5_72x44x28_1e-02.npz"}


def analytic(i):
    a = np.load(os.path.join(DATA, AN_PAT[i]))
    d = {k: float(a[k]) for k in ('Fbar', 'sigma2', 'h_T', 'h_Vne', 'h_Vee', 'g_T', 'g_Vne', 'g_Vee', 'c')}
    d['E'] = float(a['E_R12'])
    d['g_tot'] = d['g_T'] + d['g_Vne'] + d['g_Vee'] + VNN * d['sigma2']
    d['D'] = d['g_tot'] - d['E'] * d['sigma2']
    d['h'] = d['h_T'] + d['h_Vne'] + d['h_Vee']
    return d


def q(B, i, Fb_an, tkind):
    S, Tns, Tsym, Vne, Vee = B
    T = Tns if tkind == 'ns' else Tsym
    sub = [0, i]
    S2 = S[np.ix_(sub, sub)]
    comps = {'T': T[np.ix_(sub, sub)], 'Vne': Vne[np.ix_(sub, sub)], 'Vee': Vee[np.ix_(sub, sub)]}
    H2 = comps['T'] + comps['Vne'] + comps['Vee'] + VNN * S2
    S00 = S2[0, 0]; Fb = S2[0, 1] / S00; s2 = S2[1, 1] / S00 - Fb ** 2
    o = dict(Fbar=Fb, sigma2=s2)
    for X, HX in comps.items():
        o[f'g_{X}'] = (HX[1, 1] - Fb * (HX[0, 1] + HX[1, 0]) + Fb ** 2 * HX[0, 0]) / S00
        o[f'h_{X}'] = (HX[1, 0] - Fb * HX[0, 0]) / S00          # col form
    o['h'] = (H2[1, 0] - Fb * H2[0, 0]) / S00
    o['g_tot'] = (H2[1, 1] - Fb * (H2[0, 1] + H2[1, 0]) + Fb ** 2 * H2[0, 0]) / S00
    w, v = sl.eig(H2, S2); k = int(np.argmin(w.real)); v = v[:, k].real
    o['E'] = w.real[k]; o['c'] = v[1] / (v[0] + v[1] * Fb_an); o['D'] = o['g_tot'] - o['E'] * s2
    return o


def main():
    bl = []
    files = sorted(glob.glob(os.path.join(DATA, "cvmc2_s*.npz")))
    for f in files:
        X = np.load(f)['X']; n = len(X) // NBLK * NBLK
        bl.append(X[:n].reshape(NBLK, -1, *X.shape[1:]).mean(1))
    Bb = np.concatenate(bl); nb = len(Bb); tot = Bb.sum(0)
    print(f"v2 component VMC: {len(files)} seeds, {nb} blocks")
    for i, name in ((1, 'linexp0.25'), (2, 'linexp0.5')):
        an = analytic(i)
        for tk in ('ns', 'sym'):
            full = q(Bb.mean(0), i, an['Fbar'], tk)
            J = [q((tot - Bb[b]) / (nb - 1), i, an['Fbar'], tk) for b in range(nb)]
            print(f"\n== {name}   kinetic estimator: {tk}")
            print(f"{'':>8s} {'VMC':>12s} {'err':>10s} {'analytic':>12s} {'diff':>11s} {'sigma':>7s}")
            for k in full:
                jk = np.array([j[k] for j in J]); e = np.sqrt((nb - 1) / nb * np.sum((jk - jk.mean()) ** 2))
                a = an.get(k, np.nan)
                print(f"{k:>8s} {full[k]:+12.7f} {e:10.7f} {a:+12.7f} {full[k]-a:+11.7f} {(full[k]-a)/e:+7.2f}")
        # ns-vs-sym kinetic difference (correlated)
        dg = lambda B: q(B, i, an['Fbar'], 'ns')['g_T'] - q(B, i, an['Fbar'], 'sym')['g_T']
        v = dg(Bb.mean(0)); jk = np.array([dg((tot - Bb[b]) / (nb - 1)) for b in range(nb)])
        e = np.sqrt((nb - 1) / nb * np.sum((jk - jk.mean()) ** 2))
        print(f"   g_T(ns) - g_T(sym) = {v:+.6f} +/- {e:.6f}   ({v/e:+.2f} sigma)")


if __name__ == "__main__":
    main()
