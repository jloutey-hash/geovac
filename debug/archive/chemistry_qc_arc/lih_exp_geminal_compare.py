r"""Job 1 gate table: analytic exp-geminal marriage (debug/data/lih_exp_geminal_*.npz) vs the
saved multi-geminal linear-method VMC (debug/data/mg_lm_s*.npz, same trunc-1e-2 reference),
ingredient by ingredient (Fbar, sigma2, h_col, D = g_tot - E sigma2, E, c), block jackknife.
Plus the Kato-cusp estimate for the multiplier 1 + c (F - Fbar).
"""
from __future__ import annotations

import os
import sys
import glob
import numpy as np
import scipy.linalg as sl

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
NBLK = 25
VNN = 0.9950248756218905


def load_vmc():
    Sb, Hb = [], []
    for f in sorted(glob.glob(os.path.join(DATA, "mg_lm_s*.npz"))):
        d = np.load(f); S, H = d['S'], d['H']; n = len(S) // NBLK * NBLK
        Sb.append(S[:n].reshape(NBLK, -1, *S.shape[1:]).mean(1))
        Hb.append(H[:n].reshape(NBLK, -1, *H.shape[1:]).mean(1)); lab = list(d['labels'])
    return np.concatenate(Sb), np.concatenate(Hb), [str(x) for x in lab]


def ingredients(S, H, i, Fb_an):
    S00, S01, S11 = S[0, 0], S[0, i], S[i, i]
    H00, H01, H10, H11 = H[0, 0], H[0, i], H[i, 0], H[i, i]
    Fb = S01 / S00; sig2 = S11 / S00 - Fb ** 2
    hcol = (H10 - Fb * H00) / S00
    g = (H11 - Fb * (H01 + H10) + Fb ** 2 * H00) / S00
    w, v = sl.eig(H[np.ix_([0, i], [0, i])], S[np.ix_([0, i], [0, i])])
    k = int(np.argmin(w.real)); v = v[:, k].real; E = w.real[k]
    c = v[1] / (v[0] + v[1] * Fb_an)
    return dict(Fbar=Fb, sigma2=sig2, h_col=hcol, g_tot=g, D=g - E * sig2, E=E, c=c)


def main():
    Sb, Hb, lab = load_vmc(); nb = len(Sb); tS, tH = Sb.sum(0), Hb.sum(0)
    files = sorted(glob.glob(os.path.join(DATA, "lih_exp_geminal_*_1e-02*.npz")))
    for f in files:
        a = np.load(f)
        gem = str(a['geminal']); g = float(a['gam'])
        tag = os.path.basename(f)
        print(f"\n=== {tag}")
        an = dict(Fbar=float(a['Fbar']), sigma2=float(a['sigma2']), h_col=float(a['h']),
                  E=float(a['E_R12']), c=float(a['c']))
        an['g_tot'] = float(a['g_T'] + a['g_Vne'] + a['g_Vee']) + VNN * an['sigma2']
        an['D'] = an['g_tot'] - an['E'] * an['sigma2']
        print(f"  analytic: E_R12 {an['E']:.6f}  corr {(an['E']-float(a['E_trunc']))*1e3:+.3f} mHa  c {an['c']:+.6f}"
              f"  Fbar {an['Fbar']:.6f}  sigma2 {an['sigma2']:.6f}  h {an['h_col']:+.6f}")
        print("  pieces: " + "  ".join(f"{k} {float(a[k]):+.6f}" for k in
                                     ('h_T', 'h_Vne', 'h_Vee', 'g_T', 'g_Vne', 'g_Vee')))
        print(f"  grid {int(a['NXI'])}x{int(a['NETA'])}x{int(a['NPHI'])}  nphi gen/mode "
              f"{int(a['nphi_gen']) if 'nphi_gen' in a.files else 28}/{int(a['nphi_mode']) if 'nphi_mode' in a.files else 48}  wall {float(a['wall_s']):.0f}s  "
              f"rss {float(a['peak_rss_gb']):.2f} GB")
        # cusp estimate: F_rest ~ (5/6) Fbar (five spectator pairs at their mean)
        if gem == 'exp':
            f0, fp0 = 1.0, -g
        else:
            f0, fp0 = 0.0, 1.0
        c = an['c']
        cusp = c * fp0 / (1 + c * (f0 + 5.0 / 6.0 * an['Fbar'] - an['Fbar']))
        print(f"  cusp at coalescence ~ {cusp:+.4f}  (Kato antiparallel 1/2; f(0)={f0}, f'(0)={fp0})")
        lbl = f"{gem}{g:g}" if f"{gem}{g:g}" in lab else f"{gem}{g}"
        if lbl not in lab:
            print(f"  (no VMC basis function {lbl} in the saved pool)")
            continue
        i = lab.index(lbl)
        full = ingredients(Sb.mean(0), Hb.mean(0), i, an['Fbar'])
        J = [ingredients((tS - Sb[b]) / (nb - 1), (tH - Hb[b]) / (nb - 1), i, an['Fbar']) for b in range(nb)]
        print(f"  {'':>8s} {'VMC':>12s} {'err':>10s} {'analytic':>12s} {'diff':>11s} {'sigma':>6s}")
        for k in ('Fbar', 'sigma2', 'h_col', 'g_tot', 'D', 'E', 'c'):
            jk = np.array([j[k] for j in J]); err = np.sqrt((nb - 1) / nb * np.sum((jk - jk.mean()) ** 2))
            print(f"  {k:>8s} {full[k]:+12.7f} {err:10.7f} {an[k]:+12.7f} {full[k]-an[k]:+11.7f} "
                  f"{(full[k]-an[k])/err:+6.2f}")


if __name__ == "__main__":
    main()
