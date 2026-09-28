r"""Analyse debug/data/cvmc_s*.npz (lih_c_offset_vmc.py): per-component VMC of the linexp0.5 2x2
vs the analytic sparse pieces (trunc-1e-2 pi reference).  Block jackknife, 25 blocks/seed.
Also: c from the generalized eigenproblem on the TOTAL matrices (independent replicate of the
multigeminal c), c by plug-in, and the pooled (old + new) c.
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
AN = {k: float(z['linexp_' + k]) for k in ('Fbar', 'sigma2', 'h_T', 'h_Vne', 'h_Vee', 'g_T', 'g_Vne',
                                            'g_Vee', 'E0_grid', 'E_R12', 'c', 'h')}
AN['g_tot'] = AN['g_T'] + AN['g_Vne'] + AN['g_Vee'] + VNN * AN['sigma2']
AN['E'] = AN['E_R12']; AN['E0'] = AN['E0_grid']; AN['h_col'] = AN['h']; AN['h_row'] = AN['h']
AN['D'] = AN['g_tot'] - AN['E_R12'] * AN['sigma2']
for X in ('T', 'Vne', 'Vee'):
    AN[f'hrow_{X}'] = AN[f'h_{X}']; AN[f'hcol_{X}'] = AN[f'h_{X}']


def q(B):
    """B: (4,2,2) = S, H_T, H_ne, H_ee (walker/sample means)."""
    S, HT, Hne, Hee = B
    H = HT + Hne + Hee + VNN * S
    S00 = S[0, 0]; Fb = S[0, 1] / S00; s2 = S[1, 1] / S00 - Fb ** 2
    out = dict(Fbar=Fb, sigma2=s2, E0=H[0, 0] / S00)
    for X, HX in (('T', HT), ('Vne', Hne), ('Vee', Hee)):
        out[f'g_{X}'] = (HX[1, 1] - Fb * (HX[0, 1] + HX[1, 0]) + Fb ** 2 * HX[0, 0]) / S00
        out[f'hrow_{X}'] = (HX[0, 1] - Fb * HX[0, 0]) / S00
        out[f'hcol_{X}'] = (HX[1, 0] - Fb * HX[0, 0]) / S00
    out['h_col'] = (H[1, 0] - Fb * H[0, 0]) / S00
    out['h_row'] = (H[0, 1] - Fb * H[0, 0]) / S00
    out['g_tot'] = (H[1, 1] - Fb * (H[0, 1] + H[1, 0]) + Fb ** 2 * H[0, 0]) / S00
    w, v = sl.eig(H, S); k = int(np.argmin(w.real)); v = v[:, k].real
    out['E'] = w.real[k]
    out['c'] = v[1] / (v[0] + v[1] * AN['Fbar'])
    out['D'] = out['g_tot'] - out['E'] * s2
    return out


def main():
    blocks = []
    files = sorted(glob.glob(os.path.join(DATA, "cvmc_s*.npz")))
    for f in files:
        X = np.load(f)['X']; n = len(X) // NBLK * NBLK
        blocks.append(X[:n].reshape(NBLK, -1, 4, 2, 2).mean(1))
    Bb = np.concatenate(blocks); nb = len(Bb); tot = Bb.sum(0)
    full = q(Bb.mean(0)); J = [q((tot - Bb[b]) / (nb - 1)) for b in range(nb)]
    print(f"component VMC: {len(files)} seeds, {nb} blocks")
    print(f"{'':>10s} {'VMC':>13s} {'err':>10s} {'analytic':>13s} {'diff':>11s} {'sigma':>7s}")
    errs = {}
    for k in full:
        jk = np.array([j[k] for j in J]); e = np.sqrt((nb - 1) / nb * np.sum((jk - jk.mean()) ** 2))
        errs[k] = e
        an = AN.get(k, np.nan)
        print(f"{k:>10s} {full[k]:+13.7f} {e:10.7f} {an:+13.7f} {full[k]-an:+11.7f} {(full[k]-an)/e:+7.2f}")
    # per-seed c
    print("\nper-seed c (generalized eigen, analytic Fbar):")
    cs = []
    for i, bl in enumerate(blocks):
        o = q(bl.mean(0)); jj = [q((bl.sum(0) - bl[b]) / (len(bl) - 1)) for b in range(len(bl))]
        e = np.sqrt((len(bl) - 1) / len(bl) * np.sum((np.array([x['c'] for x in jj]) - np.mean([x['c'] for x in jj])) ** 2))
        cs.append((o['c'], e)); print(f"   {os.path.basename(files[i])}: c = {o['c']:+.5f} +/- {e:.5f}")
    c = np.array(cs); w = 1 / c[:, 1] ** 2
    cm = np.sum(w * c[:, 0]) / w.sum(); ce = 1 / np.sqrt(w.sum())
    chi2 = np.sum(w * (c[:, 0] - cm) ** 2)
    print(f"   inv-var {cm:+.5f} +/- {ce:.5f}   chi2/dof {chi2:.2f}/{len(c)-1}   analytic {AN['c']:+.5f}"
          f"  ({(cm-AN['c'])/ce:+.2f} sigma)")


if __name__ == "__main__":
    main()
