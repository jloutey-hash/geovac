"""Scoring driver (to be frozen by the PM) for the four-centre (XY|ZW) Legendre coefficients on the one-ray slice.

    python debug/fc_run.py G1            -> debug/data/fc_run_G1.json
    python debug/fc_run.py G2 [--lmax N] [--out PATH]

Computes F_L, L = 0..Lmax (default 260), with debug/fc_engine.py at levels A and B (B: more bits AND finer
theta, radial and sub-interval grids), and writes, for G_L = F_L / rho_dom^L (rho_dom = the larger of rho_Z, rho_W):
  - beta_eff(L) = log|G_{L+1}/G_{L-1}| / log((L+1)/(L-1)) at L = 20..min(240, Lmax-1), from level B, with the
    level A - B difference of beta_eff at each L (the precision floor);
  - the relative self-difference |F_A - F_B| / |F_B| at every L (no coefficient values are written);
  - sign pattern of F_L;
  - amplitude fits on level B: (i) beta fixed at the design value -4: G_L L^4 = sum_{k=0..K} c_k L^-k, least squares
    over L in [L0, Lmax], K = 2, 3, 4, L0 = Lmax/2 and 2 Lmax/3; (ii) beta free: log G_L = log A + beta log L +
    sum_{k=1..K} d_k L^-k, same windows.  A is the dominant point's amplitude in the F/rho^L normalisation (unnormalised
    densities, Legendre-coefficient multipoles; see fc_engine).  The subdominant point's term rho_sub^L is not fitted.
Every number written is derived from the engine; nothing is tuned.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import mpmath as mp
from flint import arb, ctx

import sys
sys.path.insert(0, str(Path(__file__).parent))
import fc_engine as FE


def to_mp(x):
    return mp.mpf(x.mid().str(mp.mp.dps + 10, radius=False))


def fits(G, Lmax):
    """amplitude fits on G_L (mpf list)."""
    out = {}
    for L0 in sorted({Lmax // 2, (2 * Lmax) // 3}):
        Ls = list(range(max(L0, 2), Lmax + 1))
        for K in (2, 3, 4):
            if len(Ls) <= K + 2:
                continue
            # (i) beta fixed at -4
            Am = mp.matrix([[mp.mpf(L) ** (-k) for k in range(K + 1)] for L in Ls])
            y = mp.matrix([G[L] * mp.mpf(L) ** 4 for L in Ls])
            c = mp.qr_solve(Am, y)[0]
            out[f'beta=-4,L0={L0},K={K}'] = dict(c0=float(c[0]), c1=float(c[1]))
            # (ii) beta free
            Am = mp.matrix([[1, mp.log(L)] + [mp.mpf(L) ** (-k) for k in range(1, K + 1)] for L in Ls])
            y = mp.matrix([mp.log(abs(G[L])) for L in Ls])
            c = mp.qr_solve(Am, y)[0]
            out[f'beta free,L0={L0},K={K}'] = dict(A=float(mp.exp(c[0]) * mp.sign(G[Lmax])), beta=float(c[1]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('case')
    ap.add_argument('--lmax', type=int, default=260)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    Lmax = a.lmax
    out_path = Path(a.out) if a.out else Path(__file__).parent / 'data' / f'fc_run_{a.case}.json'
    t0 = time.time()
    FA, ia = FE.coefficients(a.case, Lmax, 'A')
    FB, ib = FE.coefficients(a.case, Lmax, 'B')
    mp.mp.dps = 120
    lan = FE.landau(a.case)
    rho = lan['rhoZ'] if lan['dominant'] in ('Z', 'tie') else lan['rhoW']
    rm = mp.mpf(rho.numerator) / rho.denominator
    FAm = [to_mp(f) for f in FA]
    FBm = [to_mp(f) for f in FB]
    GA = [FAm[L] / rm ** L for L in range(Lmax + 1)]
    GB = [FBm[L] / rm ** L for L in range(Lmax + 1)]

    def beta(G, L):
        return mp.log(abs(G[L + 1] / G[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1))

    betas = {}
    for L in range(20, min(240, Lmax - 1) + 1):
        bB = beta(GB, L)
        betas[str(L)] = dict(beta_eff=float(bB), A_minus_B=float(beta(GA, L) - bB))
    selfdiff = [float(abs((FAm[L] - FBm[L]) / FBm[L])) for L in range(Lmax + 1)]
    signs = sorted({int(mp.sign(f)) for f in FBm})
    res = dict(case=a.case, params=FE.CASES[a.case], Lmax=Lmax,
               landau={k: (str(v) if not isinstance(v, str) else v) for k, v in lan.items()},
               rho_used=str(rho), beta_eff=betas, rel_selfdiff_A_B=selfdiff, signs=signs,
               fits_B=fits(GB, Lmax), fits_A=fits(GA, Lmax), timing_A=ia, timing_B=ib, t_total=time.time() - t0,
               conventions='unnormalised Slater densities; F_L Legendre coefficients of I(x); G_L = F_L/rho^L')
    out_path.write_text(json.dumps(res, indent=1))
    print('wrote', out_path, round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
