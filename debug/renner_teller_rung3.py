"""m-seam rung 3: one electron, central Z_c plus two bent protons, in a
MULTI-CENTRE Coulomb-Sturmian basis, Shibuya-Wulfman fixed-S kernel posing
K = D^1/2 S D^1/2 (eigenvalues k; only two-centre integrals), fixed R by the
1-D self-consistency S/k(S) = R as in geovac/sturmian_sigma_basis.kernel_energy.
SW metric from geovac/sturmian_general_sw.py.  Plan:
debug/m_seam_rung3_scoping_memo.md.

Usage:  python debug/renner_teller_rung3.py calib
        python debug/renner_teller_rung3.py sweep  NcLc,NhLh ...   (e.g. 3,2:2,1)
"""
import json
import os
import sys
import time

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from geovac import sturmian_general_sw as g  # noqa: E402
from geovac import sturmian_sigma_basis as sb  # noqa: E402

mp.mp.dps = 20


def geometry(sbend: float, scale: float):
    """Centre at origin, protons at scale*(s, 0, +-c)."""
    c = np.sqrt(1 - sbend ** 2)
    return [(0.0, 0.0, 0.0), (scale * sbend, 0.0, scale * c), (scale * sbend, 0.0, -scale * c)]


def kernel(sbend: float, S: float, cb, hb, Zc: float = 1.0):
    X = geometry(sbend, S)
    cen = [(X[0], Zc, cb), (X[1], 1.0, hb), (X[2], 1.0, hb)]
    Smat, Dg, index = g.assemble(cen)
    Sn = np.array(Smat.tolist(), dtype=float)
    d = np.sqrt(np.array([float(v) for v in Dg]))
    return d[:, None] * Sn * d[None, :], index


def blocks(index):
    c = [i for i, (_, lab) in enumerate(index) if lab[3] == 'c']
    s = [i for i, (_, lab) in enumerate(index) if lab[3] == 's']
    return c, s


def pi_branch(K, index, idx, target=None):
    """Eigenvalue of the y-parity block idx: at s = 0 (target None) the largest
    one whose eigenvector lives on |m| = 1 functions; else the one nearest target."""
    ev, V = np.linalg.eigh(K[np.ix_(idx, idx)])
    if target is not None:
        return ev[np.argmin(np.abs(ev - target))]
    m1 = np.array([index[i][1][2] == 1 for i in idx])
    w = np.sum(V[m1, :] ** 2, axis=0)
    cand = [t for t in range(len(ev)) if w[t] > 1 - 1e-9]
    return max(ev[t] for t in cand)


def fixed_R(sbend, R, cb, hb, kind, k0=None):
    """Self-consistent k of the pi branch in y-parity block `kind` at fixed R."""
    def kpi(S):
        K, index = kernel(sbend, S, cb, hb)
        c, s = blocks(index)
        return pi_branch(K, index, c if kind == 'c' else s, k0)
    lo, hi = 0.3 * R, 3.0 * R
    S = brentq(lambda SS: SS / kpi(SS) - R, lo, hi, xtol=1e-13, rtol=1e-14)
    return kpi(S)


def calib():
    """H2+ (two Z=1 centres) through the general engine == sigma module kernel posing."""
    for N, R in ((3, 2.0), (4, 2.0)):
        sig = [(n, l, 0, 'c') for (n, l) in sb.basis(N)]
        def kmax(S):
            Smat, Dg, _ = g.assemble([((0, 0, -S / 2), 1.0, sig), ((0, 0, S / 2), 1.0, sig)])
            d = np.sqrt(np.array([float(v) for v in Dg]))
            return np.linalg.eigvalsh(d[:, None] * np.array(Smat.tolist(), dtype=float) * d[None, :])[-1]
        S = brentq(lambda SS: SS / kmax(SS) - R, 1e-2, 40.0, xtol=1e-12)
        e_new = -0.5 * kmax(S) ** 2
        e_ref = sb.kernel_energy(N, R)
        print('H2+ N=%d R=%.1f  general engine %.10f  sigma module %.10f  diff %.1e'
              % (N, R, e_new, e_ref, e_new - e_ref))


def sweep(specs, R=2.0):
    out = os.path.join(ROOT, 'debug', 'data', 'renner_teller_rung3.jsonl')
    for spec in specs:
        cpart, hpart = spec.split(':')
        Nc, Lc = map(int, cpart.split(','))
        Nh, Lh = map(int, hpart.split(','))
        cb, hb = g.labels(Nc, Lc), g.labels(Nh, Lh)
        t0 = time.time()
        k0 = fixed_R(0.0, R, cb, hb, 'c')
        k0y = fixed_R(0.0, R, cb, hb, 's')
        def dE(sv):
            kx = fixed_R(sv, R, cb, hb, 'c', k0)
            ky = fixed_R(sv, R, cb, hb, 's', k0y)
            return -(kx ** 2 - ky ** 2) / 2
        h = [4e-3, 2e-3, 1e-3]
        q = [dE(v) / v ** 2 for v in h]
        r1, r2 = (4 * q[1] - q[0]) / 3, (4 * q[2] - q[1]) / 3
        row = {'centre': [Nc, Lc], 'proton': [Nh, Lh], 'N': len(cb) + 2 * len(hb),
               'E_pi': -k0 ** 2 / 2, 'linear_pair_gap': k0 - k0y, 'a2': r2,
               'a2_spread': abs(r1 - r2), 'seconds': round(time.time() - t0, 1)}
        print(json.dumps(row))
        sys.stdout.flush()
        with open(out, 'a') as f:
            f.write(json.dumps(row) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'calib':
        calib()
    else:
        sweep(sys.argv[2:])
