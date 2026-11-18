"""OP52 derivation memo: next-order endpoint germs and the between-endpoint collision (M2, v = 1), f' = NA - Y_mu.

Pinch germ of a potential term s^p against C's cusp term |r - C|^q (density factor g_q = N (-zC)^q / q! e^{-zB R}),
to relative order |a|^2 (derived here; L52's leading term plus the t-Taylor and R(x) corrections):

  I_pq ~ C_pq |a|^{p+q+3} int_0^1 w_pq(u) [ g(R0) + u(1-u) |a|^2 ( zB/(2 R0) g(R0) - lap g / (2 (p+q+5)) ) ] du,
  w_pq = u^{(p+1)/2} (1-u)^{(q+1)/2},  R0^2 = (1-u) D2^2 + u D1^2,  lap g = (zB^2 - 2 zB / R0) g,  g = e^{-zB R0},
  |a|^2 = 2 D1 D2 t,  t = x* - x.
(R(u; x)^2 = R0^2 - u(1-u)|a|^2 gives the zB/(2R0) term; the isotropic t-Hessian gives -lap g / (2(p+q+5)).)

Endpoint part J_E: the u-integral in r = R0 (u = (D2^2 - r^2)/(D2^2 - D1^2)) is elementary, F(D2) - F(D1); the
D1 part is -F(D1) e^{zB D1}, the D2 part F(D2) e^{zB D2}.

Predictions for the V' germ (v_1 = -mu^2/2, v_3 = -mu^4/24; no p = -1 term) at endpoint E:
  rho_E = [t^{7/2}] / [t^{5/2}]  from (1,1) next order, (1,3) leading, (3,1) leading.
Measured (mod p, debug/op52_deriv_local.py): rho_D2 = t^{7/2} coeff of the normalized D2 branch series; rho_D1 = that
of the normalized D1 V'-germ (B_Y - B_N)/lead; the collision residual is rho_D2 - rho_D1 (v = 1 iff nonzero).
Also the D1-block t^{7/2} coefficient of B_Y - B_N at t^{3/2}-normalization is checked.
"""
from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path

import flint
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location('loc', ROOT / 'debug' / 'op52_deriv_local.py')
loc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(loc)

r = sp.symbols('r', positive=True)


def Cpq(p, q):
    return sp.simplify(sp.pi ** sp.Rational(3, 2) * sp.gamma(-sp.Rational(p + q + 3, 2)) /
                       (sp.gamma(-sp.Rational(p, 2)) * sp.gamma(-sp.Rational(q, 2))))


def J(expr_u_r, D1, D2, zB):
    """endpoint parts of int_0^1 expr(u, r) e^{-zB r} du (expr polynomial in u, Laurent in r)."""
    Dl = D2 ** 2 - D1 ** 2
    uu = (D2 ** 2 - r ** 2) / Dl
    integrand = sp.expand(expr_u_r(uu) * 2 * r / Dl)
    F = sp.integrate(integrand * sp.exp(-zB * r), r)
    P = sp.simplify(sp.expand(F * sp.exp(zB * r)))
    assert not P.has(sp.Ei) and not P.has(sp.expint), P
    return {'D1': sp.simplify(-P.subs(r, D1)), 'D2': sp.simplify(P.subs(r, D2))}


def germ_coeffs(V, D1, D2, zB, zC, ps=(-1, 1, 3), qs=(1, 3)):
    """[t^{3/2}], [t^{5/2}], [t^{7/2}] coefficients (common N dropped) of the germ of V = {p: v_p} at both
    endpoints, through relative order |a|^2 for the (p, q) with p + q + 3 <= 7."""
    K = 2 * D1 * D2
    out = {E: {3: 0, 5: 0, 7: 0} for E in ('D1', 'D2')}
    for p in ps:
        vp = V.get(p, 0)
        if vp == 0:
            continue
        for q in qs:
            n = p + q + 3
            if n > 7:
                continue
            gq = (-zC) ** q / sp.factorial(q)
            a, b = sp.Rational(p + 1, 2), sp.Rational(q + 1, 2)
            lead = J(lambda u: u ** a * (1 - u) ** b, D1, D2, zB)
            for E in out:
                out[E][n] += vp * Cpq(p, q) * gq * K ** sp.Rational(n, 2) * lead[E]
            if n + 2 <= 7:
                corr = J(lambda u: u ** a * (1 - u) ** b * u * (1 - u) *
                         (zB / (2 * r) - (zB ** 2 - 2 * zB / r) / (2 * (p + q + 5))), D1, D2, zB)
                for E in out:
                    out[E][n + 2] += vp * Cpq(p, q) * gq * K ** sp.Rational(n + 2, 2) * corr[E]
    return {E: {k: sp.simplify(v) for k, v in d.items()} for E, d in out.items()}


def frac(x):
    x = sp.simplify(x)
    assert x.is_Rational, x
    return Fraction(int(sp.numer(x)), int(sp.denom(x)))


def main():
    emb = json.loads(loc.EMBED.read_bytes().decode())['entries']
    out = {}
    for case, vals in loc.CASES.items():
        D1, D2, zB, zC, zA = (sp.Rational(v) for v in vals)
        mu = 2 * zA
        Vp = {1: -mu ** 2 / 2, 3: -mu ** 4 / 24}                    # V' = (1 - e^{-mu s})/s, odd terms
        Pt = {-1: 1}                                                # 1/s
        gV = germ_coeffs(Vp, D1, D2, zB, zC)
        gN = germ_coeffs(Pt, D1, D2, zB, zC)
        rho = {E: gV[E][7] / gV[E][5] for E in gV}
        # D1-block quantities normalized by the point-charge lead: (B_Y - B_N) = -G_V'/G_N.lead
        d5 = -gV['D1'][5] / gN['D1'][3]
        d7 = -gV['D1'][7] / gN['D1'][3]
        pred = {'rho_D1': frac(rho['D1']), 'rho_D2': frac(rho['D2']), 'resid': frac(rho['D2'] - rho['D1']),
                'BYmBN_5': frac(d5), 'BYmBN_7': frac(d7)}
        row = {k: str(v) for k, v in pred.items()}
        for p in loc.PRIMES:
            e = emb[f'{case}|{p}']
            d1, d2 = loc.endpoint_keys(case)
            ser = {}
            for nm, key in {**d1, **d2}.items():
                P = [flint.nmod_poly(cc, p) for cc in e['piece_ops'][key]]
                ser[nm] = loc.local_halfint_space(P, loc.XSTAR[case], p, N=6)[0][0]
            diff = [(a - b) % p for a, b in zip(ser['Y-'], ser['NA'])]
            meas = {'BYmBN_5': diff[2], 'BYmBN_7': diff[3],
                    'rho_D1': diff[3] * pow(diff[2], p - 2, p) % p,
                    'rho_D2': ser['D2+'][3], 'resid': (ser['D2+'][3] - diff[3] * pow(diff[2], p - 2, p)) % p}
            for k, v in meas.items():
                row[f'{k}_match_{p}'] = v == loc.fp(pred[k], p)
        out[case] = row
        print(case, json.dumps(row, indent=1))
    (ROOT / 'debug' / 'data' / 'op52_deriv_nextorder.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
