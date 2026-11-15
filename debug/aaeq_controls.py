"""Controls for a pre-registration record on (AA|BC) at equal arms, D1 = D2 (record AE; inventory row 8, the
D1 = D2 test of the |a|^7 h_l(a) mechanism).  Nothing here computes a scored quantity (no exact split, no decay
exponent of the integrals).

Geometry: B at the origin, C at (0, 0, D), A at D (sin g, 0, cos g), x = cos g.  At x = 1, A sits on C.
Integrals: S = (1s_A 1s_A | 1s_B 1s_C) (AB's), PC = (1s_A 1s_A | 1s_B 2pz_C) (AP's), PD = (1s_A 1s_A | 1s_B 3ds_C)
(AD's).  Cases E1 (D = 2, zB = 13/10, zC = 1, zA = 9/10), E2 (D = 12/5, zB = 1, zC = 6/5, zA = 4/5).

    python debug/aaeq_controls.py referee : R2 (AB's Gaussian fits, n = 30 / 34) at two angles per case, R1 (prolate)
    python debug/aaeq_controls.py toy     : beta_eff of the exact Legendre coefficients of (1 - x)^alpha,
                                            alpha = 7/2, 9/2 (positive control for the algebraic-decay instrument)
Output debug/data/aaeq_controls.json.
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parents[1]))
import aabc_controls as AB  # noqa: E402
import aapc_controls as AP  # noqa: E402
import aadz_controls as AD  # noqa: E402

OUT = Path(__file__).parent / "data" / "aaeq_controls.json"
CASES_EQ = {'E1': dict(D='2', zB='13/10', zC='1', zA='9/10'), 'E2': dict(D='12/5', zB='1', zC='6/5', zA='4/5')}
POINTS = {'E1': (0.5, 0.0), 'E2': (-1 / 3, 0.8)}


def load():
    return json.loads(OUT.read_text()) if OUT.exists() else {}


def referee():
    d = load()
    out = {}
    t0 = time.time()
    fits = {nn: AB.fit(nn) for nn in (30, 34)}
    for case, g0 in CASES_EQ.items():
        g = {k: float(Fraction(v)) for k, v in g0.items()}
        D = g['D']
        B, C = np.zeros(3), np.array([0.0, 0.0, D])
        for x in POINTS[case]:
            A = D * np.array([np.sqrt(1 - x * x), 0.0, x])
            rec = {}
            for name, fn in (('S', AB.gauss_aabc), ('PC', AP.gauss_PC), ('PD', AD.gauss_PD)):
                rec[name] = {str(nn): fn(al, c, g['zA'], g['zB'], g['zC'], A, B, C) for nn, (al, c, _) in fits.items()}
            gr = (120, 120, 144)
            rec['R1'] = dict(S=AB.prolate_aabc(g['zA'], g['zB'], g['zC'], D, A, *gr),
                             PC=AP.prolate('PC', g['zA'], g['zB'], g['zC'], D, A, *gr),
                             PD=AD.prolate_PD(g['zA'], g['zB'], g['zC'], D, A, *gr))
            out[f'{case}_x={x:.6g}'] = dict(x=x, **rec)
            print(case, x, json.dumps(rec), round(time.time() - t0), 's', flush=True)
    d['referee'] = out
    OUT.write_text(json.dumps(d, indent=1))


def toy():
    """c_L of (1 - x)^a: |int (1-x)^a P_L| = 2^{a+1} Gamma(a+1)^2 / (Gamma(a+L+2) |Gamma(a+1-L)|) (magnitude), times
    (2L+1)/2; beta_eff(L) = log(c_{L+1}/c_{L-1}) / log((L+1)/(L-1)); expected limit -(2a + 1)."""
    import mpmath as mp
    mp.mp.dps = 60
    d = load()
    res = {}
    for a in (mp.mpf(7) / 2, mp.mpf(9) / 2):
        def c(L):
            return (2 * L + 1) / mp.mpf(2) * 2 ** (a + 1) * mp.gamma(a + 1) ** 2 / (
                mp.gamma(a + L + 2) * abs(mp.gamma(a + 1 - L)))
        beta = {L: float(mp.log(c(L + 1) / c(L - 1)) / mp.log(mp.mpf(L + 1) / (L - 1))) for L in (40, 100, 160, 220)}
        res[str(float(a))] = beta
        print('alpha', float(a), beta, flush=True)
    # sign check of the formula's magnitude against a direct quadrature at small L
    chk = {}
    for L in (1, 2, 3, 6):
        q = mp.quad(lambda x: (1 - x) ** mp.mpf(3.5) * mp.legendre(L, x), [-1, 1])
        mag = 2 ** mp.mpf(4.5) * mp.gamma(4.5) ** 2 / (mp.gamma(5.5 + L) * abs(mp.gamma(4.5 - L)))
        chk[L] = [float(q), float(mag)]
    print('check |I_L| vs quad (alpha 7/2):', chk)
    d['toy'] = dict(beta=res, check=chk)
    OUT.write_text(json.dumps(d, indent=1))


def amplitude():
    """Predicted leading Legendre coefficients at equal arms, from the local analysis at a = A - C -> 0 (no engine):
    the non-analytic part of the BC density at C convolved with V_A's s^3 term, v3 = (2/3) zA^4.
    K(m, n) = c_m c_n / c_{m+n+3}, c_q = int |r|^q e^{-ik.r} d^3r / k^{-q-3}; h2(u)|u| = h2(grad)|u|^5 / 15.
      S : N_B N_C e^{-zB D} (-zC) v3 K(1,3) |a|^7
      PC: N_B N_C' e^{-zB D} zC v3 [3 K33/(2D) + zB (K33 - (99/15) K53)/3] |a|^9   (a_z = -|a|^2/(2D) on the sphere;
          terms a_z|a|^7, and zB u_z^2 |u| from e^{-zB r} about C; vanishes at zB D = 9)
      PD: N_B N_d e^{-zB D} zC v3 (99/30) K53 |a|^9   (h2(a) = -|a|^2/2 + O(|a|^4))
    |a|^2 = 2 D^2 (1 - x); F_L = amplitude (2 D^2)^alpha c_L(alpha), c_L = (2L+1)/2 int (1-x)^alpha P_L, sign (-1)^L."""
    import mpmath as mp
    mp.mp.dps = 50
    d = load()

    def c(q):
        return 2 ** (q + 3) * mp.pi ** 1.5 * mp.gamma((q + 3) / mp.mpf(2)) / mp.gamma(-q / mp.mpf(2))

    def K(m, n):
        return c(m) * c(n) / c(m + n + 3)

    def cL(a, L):
        mag = 2 ** (a + 1) * mp.gamma(a + 1) ** 2 / (mp.gamma(a + L + 2) * abs(mp.gamma(a + 1 - L)))
        return (-1) ** L * (2 * L + 1) / mp.mpf(2) * mag
    out = {}
    for case, g0 in CASES_EQ.items():
        D, zB, zC, zA = (mp.mpf(Fraction(g0[k]).numerator) / Fraction(g0[k]).denominator for k in ('D', 'zB', 'zC', 'zA'))
        v3 = mp.mpf(2) / 3 * zA ** 4
        NB = mp.sqrt(zB ** 3 / mp.pi)
        f0 = mp.exp(-zB * D)
        amp = dict(S=(NB * mp.sqrt(zC ** 3 / mp.pi) * f0 * (-zC) * v3 * K(1, 3), mp.mpf(7) / 2),
                   PC=(NB * mp.sqrt(zC ** 5 / mp.pi) * f0 * zC * v3 *
                       (3 * K(3, 3) / (2 * D) + zB * (K(3, 3) - mp.mpf(99) / 15 * K(5, 3)) / 3), mp.mpf(9) / 2),
                   PD=(NB * mp.sqrt(5 * (2 * zC) ** 7 / (4 * mp.pi * 720)) * f0 * zC * v3 * mp.mpf(99) / 30 * K(5, 3),
                       mp.mpf(9) / 2))
        rec = {}
        for name, (A, a) in amp.items():
            rec[name] = dict(amplitude_a_power=mp.nstr(A, 20), alpha=float(a),
                             pred_F={L: mp.nstr(A * (2 * D * D) ** a * cL(a, L), 15) for L in (180, 200, 219, 220, 221)})
        rec['PC_cancellation_zB_D'] = float(zB * D)
        out[case] = rec
        print(case, json.dumps(rec, indent=1), flush=True)
    d['amplitude'] = out
    OUT.write_text(json.dumps(d, indent=1))


if __name__ == '__main__':
    {'referee': referee, 'toy': toy, 'amplitude': amplitude}[sys.argv[1]]()
