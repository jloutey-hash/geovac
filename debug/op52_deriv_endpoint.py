"""OP52 derivation memo: does the endpoint part of the Landau-pinch germ predict the pieces' local branch series?

Test (exact, then mod p): within the D1 block of f' = NA - Y_mu, the normalized branch series of the NA piece (B_N)
and of a Y half (B_Y) differ first at t^{5/2} (t = x* - x).  If each D1 piece's germ is the D1-endpoint part of the
pinch germ of the potential it sees (linear in V), the t^{5/2} coefficient of B_Y - B_N is

    (mu^2/2) * [C_{1,1} (2 D1 D2)^{5/2} E_1] / [C_{-1,1} (2 D1 D2)^{3/2} E_{-1}] ,

E_p = the e^{-zB D1}-endpoint part (prefactor removed) of int_0^1 u^{(p+1)/2} (1-u) e^{-zB R_B(u)} du,
R_B(u)^2 = (1-u) D2^2 + u D1^2  (L52 segment formula, q = 1 cusp of C).  The p = -1 parts cancel in the difference
because 1/s and e^{-mu s}/s share v_{-1} = 1.  Also computed: the D2-endpoint version (for the record of what the
same reading would give at D2).  Compared with the mod-p value from debug/op52_deriv_local.py's local series.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
import importlib.util
spec = importlib.util.spec_from_file_location('loc', ROOT / 'debug' / 'op52_deriv_local.py')
loc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(loc)

r, u = sp.symbols('r u', positive=True)


def C(p, q):
    return sp.pi ** sp.Rational(3, 2) * sp.gamma(-sp.Rational(p + q + 3, 2)) / (
        sp.gamma(-sp.Rational(p, 2)) * sp.gamma(-sp.Rational(q, 2)))


def endpoint_parts(p, D1, D2, zB):
    """int_0^1 u^{(p+1)/2} (1-u) e^{-zB R} du = [F]_{D1}^{D2}-type split: returns (part at D1, part at D2), each
    with its exponential e^{-zB D} removed."""
    Dl = D2 ** 2 - D1 ** 2
    uu = (D2 ** 2 - r ** 2) / Dl
    integrand = sp.expand(uu ** sp.Rational(p + 1, 2) * (1 - uu) * 2 * r / Dl)     # times e^{-zB r}, r: D1..D2
    F = sp.integrate(integrand * sp.exp(-zB * r), r)                               # elementary for integer powers
    P = sp.simplify(F * sp.exp(zB * r))                                            # polynomial part
    return -P.subs(r, D1), P.subs(r, D2)                                           # int_{D1}^{D2} = F(D2) - F(D1)


def pred_check(x):
    return sp.simplify(x).is_Rational


def main():
    out = {}
    local = json.loads((ROOT / 'debug' / 'data' / 'op52_deriv_local.json').read_text())
    emb = json.loads(loc.EMBED.read_bytes().decode())['entries']
    import flint
    for case, vals in loc.CASES.items():
        D1, D2, zB, zC, zA = (sp.Rational(v) for v in vals)
        mu = 2 * zA
        e1m, e2m = endpoint_parts(-1, D1, D2, zB)
        e11, e21 = endpoint_parts(1, D1, D2, zB)
        ratioC = sp.simplify(C(1, 1) / C(-1, 1))
        assert pred_check(ratioC)
        pred_D1 = (sp.simplify(mu ** 2 / 2 * ratioC * (2 * D1 * D2) * e11 / e1m))
        pred_D2 = (sp.simplify(mu ** 2 / 2 * ratioC * (2 * D1 * D2) * e21 / e2m))
        row = {'C11/C-11': str(ratioC), 'E-1_D1': str(e1m), 'E1_D1': str(e11), 'pred_D1': str(pred_D1),
               'pred_D2': str(pred_D2)}
        for p in loc.PRIMES:
            e = emb[f'{case}|{p}']
            P = [flint.nmod_poly(cc, p) for cc in e['piece_ops'][loc.endpoint_keys(case)[0]['NA']]]
            BN = loc.local_halfint_space(P, loc.XSTAR[case], p)[0][0]
            P = [flint.nmod_poly(cc, p) for cc in e['piece_ops'][loc.endpoint_keys(case)[0]['Y-']]]
            BY = loc.local_halfint_space(P, loc.XSTAR[case], p)[0][0]
            meas = (BY[2] - BN[2]) % p
            q = Fraction(int(sp.numer(pred_D1)), int(sp.denom(pred_D1)))
            row[f'meas_{p}'] = meas
            row[f'pred_D1_mod_{p}'] = loc.fp(q, p)
            row[f'match_D1_{p}'] = meas == loc.fp(q, p)
            q2 = Fraction(int(sp.numer(pred_D2)), int(sp.denom(pred_D2)))
            row[f'match_D2_{p}'] = meas == loc.fp(q2, p)
            row[f'meas_ratrec_{p}'] = str(loc.ratrec(meas, p, 10 ** 5))
        out[case] = row
        print(case, row)
    (ROOT / 'debug' / 'data' / 'op52_deriv_endpoint.json').write_text(json.dumps(out, indent=1, default=str))


if __name__ == '__main__':
    main()
