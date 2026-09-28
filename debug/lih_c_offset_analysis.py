r"""Job 2 (c-offset): decompose the saved multi-geminal linear-method VMC samples
(debug/data/mg_lm_s*.npz, trunc-1e-2 pi reference, guide 1+c(F-Fbar)) into the analytic 2x2's
scalar ingredients for the linexp0.5 geminal and compare with the analytic sparse pieces.

From the basis {phi0 = Psi_CI, phi1 = F Psi_CI} (weighted estimators S_ab, H_ab, H non-symmetric):
  Fbar = S01/S00,  E0 = H00/S00,  sigma2 = S11/S00 - Fbar^2,
  h(row) = (H01 - Fbar H00)/S00,  h(col) = (H10 - Fbar H00)/S00,
  g_tot  = (H11 - Fbar (H01+H10) + Fbar^2 H00)/S00      (includes V_NN * sigma2)
and c = h / (E sigma2 - g) etc.  Errors: pooled block jackknife (25 blocks/seed x 8 seeds) on the
FULL nonlinear functional (ratio estimators included), plus a bootstrap over blocks, plus
the 'plug-in' c from the 2x2 built from VMC matrix elements vs the generalized-eigen c.
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

z = np.load(os.path.join(DATA, "lih_pi_m16_energy_trunc1e-2.npz"))
AN = dict(Fbar=float(z['linexp_Fbar']), sigma2=float(z['linexp_sigma2']), h=float(z['linexp_h']),
          g_elec=float(z['linexp_g_T'] + z['linexp_g_Vne'] + z['linexp_g_Vee']),
          E0=float(z['linexp_E0_grid']), E_R12=float(z['linexp_E_R12']), c=float(z['linexp_c']))
AN['g_tot'] = AN['g_elec'] + VNN * AN['sigma2']


def scalars(S, H, i=2):
    sub = [0, i]
    S00, S01, S11 = S[0, 0], S[0, i], S[i, i]
    H00, H01, H10, H11 = H[0, 0], H[0, i], H[i, 0], H[i, i]
    Fb = S01 / S00
    E0 = H00 / S00
    sig2 = S11 / S00 - Fb ** 2
    h_row = (H01 - Fb * H00) / S00
    h_col = (H10 - Fb * H00) / S00
    g = (H11 - Fb * (H01 + H10) + Fb ** 2 * H00) / S00
    # generalized eigen (the production LM estimator), c with the ANALYTIC Fbar (as analyze.py)
    w, v = sl.eig(H[np.ix_(sub, sub)], S[np.ix_(sub, sub)])
    k = int(np.argmin(w.real)); v = v[:, k].real; E = w.real[k]
    c_eig = v[1] / (v[0] + v[1] * AN['Fbar'])
    c_eig_selfF = v[1] / (v[0] + v[1] * Fb)
    # plug-in symmetric 2x2 from the scalars (h symmetrized)
    hs = 0.5 * (h_row + h_col)
    Mt = np.array([[E0, hs / np.sqrt(sig2)], [hs / np.sqrt(sig2), g / sig2]])
    ww, vv = np.linalg.eigh(Mt); v0 = vv[:, 0]
    c_plug = (v0[1] / v0[0]) / np.sqrt(sig2)
    return dict(Fbar=Fb, E0=E0, sigma2=sig2, h_row=h_row, h_col=h_col, h_sym=hs, g_tot=g,
                g_over_s2=g / sig2, E_eig=E, c_eig=c_eig, c_eig_selfF=c_eig_selfF,
                E_plug=ww[0], c_plug=c_plug, dE_eig=E - E0,
                c_formula=hs / (ww[0] * sig2 - g))


def main():
    files = sorted(glob.glob(os.path.join(DATA, "mg_lm_s*.npz")))
    Sb, Hb = [], []
    for f in files:
        d = np.load(f)
        S, H = d['S'], d['H']
        n = len(S) // NBLK * NBLK
        Sb.append(S[:n].reshape(NBLK, -1, *S.shape[1:]).mean(1))
        Hb.append(H[:n].reshape(NBLK, -1, *H.shape[1:]).mean(1))
        lab = list(d['labels'])
    assert lab[2] == 'linexp0.5'
    Sb = np.concatenate(Sb); Hb = np.concatenate(Hb); nb = len(Sb)
    full = scalars(Sb.mean(0), Hb.mean(0))
    tS, tH = Sb.sum(0), Hb.sum(0)
    J = [scalars((tS - Sb[b]) / (nb - 1), (tH - Hb[b]) / (nb - 1)) for b in range(nb)]
    rng = np.random.default_rng(0)
    B = []
    for _ in range(2000):
        idx = rng.integers(0, nb, nb)
        B.append(scalars(Sb[idx].mean(0), Hb[idx].mean(0)))
    print(f"blocks {nb} ({len(files)} seeds)")
    print(f"{'quantity':>12s} {'VMC':>13s} {'jk err':>10s} {'jk-bias-corr':>13s} {'boot mean':>13s}"
          f" {'boot err':>10s} {'analytic':>13s} {'diff/err':>9s}")
    anmap = dict(Fbar='Fbar', sigma2='sigma2', h_sym='h', h_row='h', h_col='h', g_tot='g_tot',
                 E0='E0', E_eig='E_R12', c_eig='c', c_plug='c', c_formula='c', E_plug='E_R12')
    for k in full:
        jk = np.array([j[k] for j in J]); bs = np.array([b[k] for b in B])
        err = np.sqrt((nb - 1) / nb * np.sum((jk - jk.mean()) ** 2))
        bc = nb * full[k] - (nb - 1) * jk.mean()          # jackknife bias-corrected
        an = AN.get(anmap.get(k, ''), np.nan)
        dd = (full[k] - an) / err if np.isfinite(an) else np.nan
        print(f"{k:>12s} {full[k]:+13.7f} {err:10.7f} {bc:+13.7f} {bs.mean():+13.7f} {bs.std():10.7f}"
              f" {an:+13.7f} {dd:+9.2f}")
    # implied-c sensitivity: c recomputed with analytic sigma2/g/h swapped in one at a time
    print("\nc sensitivity (analytic 2x2 with ONE ingredient replaced by its VMC value):")
    E0a = AN['E0'] - VNN

    def c_of(s2, h, ge, E0e):
        Mt = np.array([[E0e, h / np.sqrt(s2)], [h / np.sqrt(s2), ge / s2]])
        ww, vv = np.linalg.eigh(Mt); v0 = vv[:, 0]
        return (v0[1] / v0[0]) / np.sqrt(s2), ww[0] + VNN
    base = c_of(AN['sigma2'], AN['h'], AN['g_elec'], E0a)
    print(f"   analytic         c={base[0]:+.6f} E={base[1]:.6f}")
    vs = dict(sigma2=full['sigma2'], h=full['h_sym'], g_elec=full['g_tot'] - VNN * full['sigma2'],
              E0e=full['E0'] - VNN)
    for k in ['sigma2', 'h', 'g_elec', 'E0e']:
        args = dict(s2=AN['sigma2'], h=AN['h'], ge=AN['g_elec'], E0e=E0a)
        key = {'sigma2': 's2', 'h': 'h', 'g_elec': 'ge', 'E0e': 'E0e'}[k]
        args[key] = vs[k]
        if k == 'sigma2':   # keep g/sigma2 (the normalized diagonal) consistent: scale g too
            pass
        cc = c_of(**args)
        print(f"   swap {k:7s}     c={cc[0]:+.6f} E={cc[1]:.6f}")
    # normalized-diagonal comparison (the quantity the 2x2 actually uses)
    print(f"\n g/sigma2 (tot): VMC {full['g_over_s2']:+.6f}  analytic {AN['g_tot']/AN['sigma2']:+.6f}")
    print(f" h/sigma  : VMC {full['h_sym']/np.sqrt(full['sigma2']):+.6f}  analytic "
          f"{AN['h']/np.sqrt(AN['sigma2']):+.6f}")


if __name__ == "__main__":
    main()
