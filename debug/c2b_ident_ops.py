"""EC1-EC3 (docs/predictions/2026-10-22_c2b_identification.md): the order-two operators tested,
in d/dx form a2 y'' + a1 y' + a0 y = 0 with integer coefficient lists (ascending in x), plus the
base point and the loops used for the monodromy.

  C2b              the order-two composition factor of L9 at gamma = 3/2 (theta_t = 2 x d/dx)
  heun             negative control: Heun operator, singular points 0, 1/4, 25/4, infinity,
                   exponent difference 0 at each, accessory parameter 7/3
  legendre_gauged  positive control: Legendre's K(l) operator transported by A = 1 + (l^2 + 3) d/dl
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G  # noqa: E402

x = sp.Symbol("x")
DATA = Path(__file__).parent / "data"


def roots_of(polys):
    rs = []
    for f in polys:
        pf = sp.Poly(f, x)
        if pf.degree() >= 1:
            rs += [mp.mpc(r) for r in mp.polyroots([int(c) for c in pf.all_coeffs()],
                                                    maxsteps=800, extraprec=600)]
    return rs


def int_lists(exprs):
    """common-denominator integer coefficient lists of polynomial expressions"""
    polys = [sp.Poly(sp.expand(e), x) for e in exprs]
    den = sp.ilcm(*[sp.Rational(c).q for p in polys for c in p.all_coeffs()])
    return [[int(c * den) for c in p.all_coeffs()[::-1]] for p in polys]


def c2b():
    D = json.load(open(DATA / "correlated_triangle_g32_series_final.json"))
    L = G.parse_op(D["factor_C2b"])
    a = G.dform(L)
    a = [[int(v) for v in aj.coeffs()] for aj in a]
    a2 = sum(c * x**i for i, c in enumerate(a[2]))
    return a, [f for f, m in sp.factor_list(a2)[1]]


def heun(q=sp.Rational(7, 3)):
    a2 = sp.expand(16 * x * (x - sp.Rational(1, 4)) * (x - sp.Rational(25, 4)))
    a1 = sp.diff(a2, x)              # gamma = delta = epsilon = 1
    a0 = 16 * (x - q)                # alpha = beta = 1
    return int_lists([a0, a1, a2]), [a2]


def transported(a0, a1, a2, r0, r1):
    """operator annihilating z = r0 y + r1 y' for y in Sol(a2 y'' + a1 y' + a0 y)"""
    p, q = a1 / a2, a0 / a2
    s0 = sp.diff(r0, x) - r1 * q
    s1 = r0 + sp.diff(r1, x) - r1 * p
    u0 = sp.diff(s0, x) - s1 * q
    u1 = s0 + sp.diff(s1, x) - s1 * p
    P, Q = sp.symbols("P Q")
    sol = sp.solve([u0 + P * s0 + Q * r0, u1 + P * s1 + Q * r1], [P, Q], dict=True)[0]
    Pf, Qf = sp.together(sol[P]), sp.together(sol[Q])
    den = sp.lcm(sp.fraction(Pf)[1], sp.fraction(Qf)[1])
    A2 = sp.factor(den)
    A1 = sp.cancel(Pf * den)
    A0 = sp.cancel(Qf * den)
    return [A0, A1, A2]


def legendre_gauged():
    ops = transported(-sp.Rational(1, 4), 1 - 2 * x, x * (1 - x), sp.Integer(1), x**2 + 3)
    a = int_lists(ops)
    a2 = sum(c * x**i for i, c in enumerate(a[2]))
    return a, [f for f, m in sp.factor_list(a2)[1]]


def get(which):
    from c2b_ident_B_monodromy import Op2
    if which == "C2b":
        a, facs = c2b()
        cusps = [(sp.Integer(0), "0.05"), (sp.Rational(1, 4), "0.05"), (sp.Rational(25, 4), "0.6")]
        base = mp.mpc("2.0", "0.3")
    elif which == "heun":
        a, facs = heun()
        cusps = [(sp.Integer(0), "0.05"), (sp.Rational(1, 4), "0.05"), (sp.Rational(25, 4), "0.6")]
        base = mp.mpc("2.0", "0.3")
    elif which == "legendre_gauged":
        a, facs = legendre_gauged()
        cusps = [(sp.Integer(0), "0.2"), (sp.Integer(1), "0.2")]
        base = mp.mpc("0.5", "0.1")
    else:
        raise ValueError(which)
    roots = roots_of(facs + [x])
    loops = [(mp.mpf(s.p) / s.q, mp.mpf(r)) for s, r in cusps]
    return Op2(a, roots), base, loops
