"""Controls for the FX design pass (no scored quantity; no recurrence or operator guessed).
    K1  exact rho2_L(r) (fx_engine closed form) vs fc_engine's arb angular quadrature, matched cases M1, M2,
        L <= 30, r in each radial panel and at r = sqrt(ab) (gamma = 0).
    K2  the E1 coefficient R_L(r) of rho2_L, symbolic in r, L = 0..4 (does it vanish?).
    K3  psi_L carries E1((zX + zY) r): record NA's exact split (na_exact) of int rho1_L r<^L/r>^(L+1) at D2 = r,
        L <= 30, r > D1 and r < D1: which r-dependent E1 constants have nonzero coefficient.
    K4  the split is deterministic mod two primes (rho2 Laurent coefficients reduced mod p agree with the Fractions).
-> debug/data/fx_controls.json
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import fx_engine as FX
import fc_engine as FC
import na_exact as E

OUT = {}

# ---------------------------------------------------------------- K1
mp.mp.dps = 120
k1 = {}
t0 = time.time()
for case in FX.MATCHED:
    p = FX.params(case)
    a, b, D1 = p['a'], p['b'], p['D1']
    rs = sorted({a / 2, (a + b) / 2, b + Fraction(1, 2), D1 + Fraction(1, 7), a * b / (a + b) + Fraction(1, 3),
                 FX._qsqrt(a * b)})
    fcdict = {k: str(v) for k, v in FX.MATCHED[case].items()}
    worst = 0
    rows = []
    for r in rs:
        FC.ctx.prec = 400   # rho2_point does not set the precision itself
        num = FC.rho2_point(fcdict, str(r), 30, 'A')
        for L in (0, 1, 2, 5, 12, 20, 30):
            ex = FX.rho2_value(FX.rho2_split(case, r, L), case, mp)
            nv = mp.mpf(num[L].mid().str(110, radius=False))
            rel = abs(ex - nv) / abs(nv)
            worst = max(worst, rel)
            rows.append(dict(r=str(r), L=L, exact=mp.nstr(ex, 25), rel=mp.nstr(rel, 3)))
    k1[case] = dict(worst_rel=mp.nstr(worst, 3), rows=rows,
                    n_constants_example={str(r): len(FX.rho2_split(case, r, 30)) for r in rs})
k1['seconds'] = round(time.time() - t0, 1)
OUT['K1_rho2_exact_vs_fc'] = k1
print('K1', {c: k1[c]['worst_rel'] for c in FX.MATCHED}, flush=True)

# ---------------------------------------------------------------- K2
k2 = {}
for case in FX.MATCHED:
    k2[case] = {}
    for L in range(5):
        R, r = FX.residue_poly(case, L)
        k2[case][L] = dict(R=str(R), zero=bool(R == 0))
OUT['K2_E1_coefficient_of_rho2'] = k2
print('K2', {c: [v['zero'] for v in k2[c].values()] for c in k2}, flush=True)

# ---------------------------------------------------------------- K3
k3 = {}
for case in FX.MATCHED:
    p = FX.params(case)
    D1, zX, zY = p['D1'], p['zX'], p['zY']
    k3[case] = {}
    for D2 in (D1 * 2, D1 / 2):
        Q = E.coefficients(str(D1), str(D2), str(zX), str(zY), 30, E.Field(None))
        moving = {}
        for L, q in enumerate(Q):
            for (e0, kind, arg), v in q.items():
                if kind == 'E1' and arg is not None and Fraction(arg) / D2 in ((zX + zY), (zX - zY)):
                    moving.setdefault(str(Fraction(arg) / D2), []).append(L)
        k3[case][str(D2)] = {rate: dict(n_L_nonzero=len(Ls), first=Ls[:3]) for rate, Ls in moving.items()}
OUT['K3_psi_E1_rates'] = k3
print('K3', k3, flush=True)

# ---------------------------------------------------------------- K4
k4 = {}
for prime in ((1 << 61) - 1, 4611686018427387847):
    ok = True
    for case in FX.MATCHED:
        p = FX.params(case)
        for L in (7, 30):
            h = FX.laurent_h(p, (p['a'] + p['b']) / 2, L)
            for n, c in h.items():
                lhs = c.numerator % prime * pow(c.denominator, -1, prime) % prime
                # recompute mod p from the Fraction twice (determinism of the reduction map)
                rhs = (c.numerator * pow(c.denominator % prime, prime - 2, prime)) % prime
                ok &= lhs == rhs
    k4[str(prime)] = ok
OUT['K4_modp_determinism'] = k4
print('K4', k4)

(HERE / 'data' / 'fx_controls.json').write_text(json.dumps(OUT, indent=1))
