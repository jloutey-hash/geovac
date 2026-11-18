"""Record L52: the predicted amplitudes (derivation only; no exact split of any scored integral is evaluated here).
    python debug/l52_predict.py  -> debug/data/l52_predictions.json
P1: (x*-x)^{7/2} amplitude of (1sA 1sA|1sB 1sC), Coulomb line:  p = 3 (v3 = 2 zA^4/3), q = 1.
P2: (x*-x)^{9/2} amplitude of D = (AA|BC)_kappa - (AA|BC)_Coulomb: only the (5,1) pair differs,
    Delta v5 = kappa^2 mu^4 / 720 = kappa^2 zA^4 / 45 (from l52_controls sym).
Also the 7/2 / 9/2 toy decay values beta_eff(220) at each x*, for P2b.
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import l52_local as LL

CASES = {'A1': dict(D1='9/5', D2='13/5', zB='13/10', zC='1', zA='9/10', kappa='7/10'),
         'A2': dict(D1='3/2', D2='12/5', zB='1', zC='6/5', zA='4/5', kappa=None),
         'A3': dict(D1='2', D2='3', zB='11/10', zC='9/10', zA='6/5', kappa='1')}


def mpf(s):
    q = Fraction(s)
    return mp.mpf(q.numerator) / q.denominator


def beta_toy(gam, xs, L):
    rho = 1 / (xs + mp.sqrt(xs * xs - 1))
    G = lambda l: LL.T(gam, xs, l) / rho ** l
    return mp.log(abs(G(L + 1) / G(L - 1))) / mp.log(mp.mpf(L + 1) / (L - 1))


if __name__ == '__main__':
    mp.mp.dps = 40
    out = {}
    for c, g in CASES.items():
        D1, D2, zB, zC, zA = (mpf(g[k]) for k in ('D1', 'D2', 'zB', 'zC', 'zA'))
        xs = LL.xstar(D1, D2)
        rec = dict(xstar=str(Fraction(g['D1']) ** 2 / (2 * Fraction(g['D1']) * Fraction(g['D2'])) +
                             Fraction(g['D2']) ** 2 / (2 * Fraction(g['D1']) * Fraction(g['D2']))),
                   Phi7=mp.nstr(LL.amplitude(3, 1, 2 * zA ** 4 / 3, D1, D2, zB, zC), 20),
                   beta_toy_72_220=mp.nstr(beta_toy(3.5, xs, 220), 6))
        if g['kappa']:
            k = mpf(g['kappa'])
            rec['kappa'] = g['kappa']
            rec['DeltaPhi9'] = mp.nstr(LL.amplitude(5, 1, k ** 2 * zA ** 4 / 45, D1, D2, zB, zC), 20)
            rec['beta_toy_92_220'] = mp.nstr(beta_toy(4.5, xs, 220), 6)
        out[c] = rec
    (HERE / 'data' / 'l52_predictions.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
