"""Controls for record AE2 (follow-up to AE, PI-chosen 2026-11-15): equal arms at fresh geometries E3 (zB D = 9, where
the closed form says PC's (1 - x)^{9/2} coefficient vanishes) and E4 (generic).  No scored quantity is computed.

    python debug/aaeq2_controls.py referee   : R2 (n = 30 / 34) and R1 at two angles per case, S, PC, PD
    python debug/aaeq2_controls.py toy       : beta_eff of (1 - x)^alpha Legendre coefficients, alpha = 7/2, 9/2, 11/2,
                                               at L = 220, 400, 580; the sign rule checked by quadrature to L = 21
    python debug/aaeq2_controls.py amplitude : AE's closed-form leading amplitudes with the CORRECTED sign rule
                                               sign c_L(alpha) = sign(Gamma(alpha+1)/Gamma(alpha+1-L)) (-1)^L
Output debug/data/aaeq2_controls.json.
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

OUT = Path(__file__).parent / "data" / "aaeq2_controls.json"
CASES = {'E3': dict(D='3', zB='3', zC='1', zA='9/10'), 'E4': dict(D='5/2', zB='7/10', zC='13/10', zA='11/10')}
POINTS = {'E3': (0.5, 0.0), 'E4': (-1 / 3, 0.8)}
LS = (220, 400, 579, 580, 581)


def load():
    return json.loads(OUT.read_text()) if OUT.exists() else {}


def referee():
    d = load()
    out = {}
    t0 = time.time()
    fits = {nn: AB.fit(nn) for nn in (30, 34)}
    for case, g0 in CASES.items():
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


def _cL(mp, a, L):
    """c_L(alpha) = (2L+1)/2 int (1-x)^alpha P_L, exact magnitude and the corrected sign."""
    mag = 2 ** (a + 1) * mp.gamma(a + 1) ** 2 / (mp.gamma(a + L + 2) * abs(mp.gamma(a + 1 - L)))
    sgn = mp.sign(mp.gamma(a + 1) / mp.gamma(a + 1 - L)) * (-1) ** L
    return sgn * (2 * L + 1) / mp.mpf(2) * mag


def toy():
    import mpmath as mp
    mp.mp.dps = 60
    d = load()
    res, chk = {}, {}
    for a in (mp.mpf(7) / 2, mp.mpf(9) / 2, mp.mpf(11) / 2):
        res[str(float(a))] = {L: float(mp.log(abs(_cL(mp, a, L + 1) / _cL(mp, a, L - 1))) /
                                       mp.log(mp.mpf(L + 1) / (L - 1))) for L in (220, 400, 580)}
        rows = []
        for L in (1, 2, 5, 6, 7, 10, 21):
            q = (2 * L + 1) / mp.mpf(2) * mp.quad(lambda x: (1 - x) ** a * mp.legendre(L, x),
                                                   [-1, 0, 0.5, 0.9, 0.99, 1], maxdegree=10)
            rows.append((L, float(q), float(_cL(mp, a, L))))
        chk[str(float(a))] = rows
        print('alpha', float(a), res[str(float(a))], '\n  sign/magnitude check (L, quad, formula):', rows, flush=True)
    d['toy'] = dict(beta=res, check=chk)
    OUT.write_text(json.dumps(d, indent=1))


def amplitude():
    import mpmath as mp
    mp.mp.dps = 50
    d = load()

    def c(q):
        return 2 ** (q + 3) * mp.pi ** 1.5 * mp.gamma((q + 3) / mp.mpf(2)) / mp.gamma(-q / mp.mpf(2))

    def K(m, n):
        return c(m) * c(n) / c(m + n + 3)
    out = {}
    for case, g0 in CASES.items():
        D, zB, zC, zA = (mp.mpf(Fraction(g0[k]).numerator) / Fraction(g0[k]).denominator for k in ('D', 'zB', 'zC', 'zA'))
        v3 = mp.mpf(2) / 3 * zA ** 4
        NB = mp.sqrt(zB ** 3 / mp.pi)
        f0 = mp.exp(-zB * D)
        bracket = 3 * K(3, 3) / (2 * D) + zB * (K(3, 3) - mp.mpf(99) / 15 * K(5, 3)) / 3
        amp = dict(S=(NB * mp.sqrt(zC ** 3 / mp.pi) * f0 * (-zC) * v3 * K(1, 3), mp.mpf(7) / 2),
                   PC=(NB * mp.sqrt(zC ** 5 / mp.pi) * f0 * zC * v3 * bracket, mp.mpf(9) / 2),
                   PD=(NB * mp.sqrt(5 * (2 * zC) ** 7 / (4 * mp.pi * 720)) * f0 * zC * v3 * mp.mpf(99) / 30 * K(5, 3),
                       mp.mpf(9) / 2))
        rec = {name: dict(amplitude=mp.nstr(A, 20), alpha=float(a),
                          pred_F={L: mp.nstr(A * (2 * D * D) ** a * _cL(mp, a, L), 15) for L in LS})
               for name, (A, a) in amp.items()}
        rec['zB_D'] = float(zB * D)
        rec['PC_bracket'] = mp.nstr(bracket, 10)
        out[case] = rec
        print(case, json.dumps(rec, indent=1), flush=True)
    d['amplitude'] = out
    OUT.write_text(json.dumps(d, indent=1))


def amp11():
    """PC at E3: leading (1 - x)^{11/2} amplitude from the general local expansion (aaeq2_local.py), and the predicted
    F_L = A_{11/2} c_L(11/2); also every case's t^{7/2}, t^{9/2}, t^{11/2} coefficients (next-order ratios)."""
    import mpmath as mp
    import sympy as sp
    import aaeq2_local as LOC
    mp.mp.dps = 50
    d = load()
    out = {}
    for case, g0 in CASES.items():
        D, zB, zC, zA = (sp.Rational(g0[k]) for k in ('D', 'zB', 'zC', 'zA'))
        NB = mp.sqrt(mp.mpf(sp.N(zB, 50)) ** 3 / mp.pi)
        z = mp.mpf(sp.N(zC, 50))
        N = {'S': mp.sqrt(z ** 3 / mp.pi), 'PC': mp.sqrt(z ** 5 / mp.pi),
             'PD': mp.sqrt(5 * (2 * z) ** 7 / (4 * mp.pi * 720))}
        out[case] = {kind: {str(k): mp.nstr(v * NB * N[kind], 15) for k, v in sorted(LOC.expansion(kind, D, zB, zC, zA).items())}
                     for kind in ('S', 'PC', 'PD')}
    A = mp.mpf(out['E3']['PC']['11/2'])
    out['E3_PC_pred_F'] = {L: mp.nstr(A * _cL(mp, mp.mpf(11) / 2, L), 15) for L in LS}
    print(json.dumps(out, indent=1))
    d['local_expansion'] = out
    OUT.write_text(json.dumps(d, indent=1))


if __name__ == '__main__':
    {'referee': referee, 'toy': toy, 'amplitude': amplitude, 'amp11': amp11}[sys.argv[1]]()
