"""Positive-geometry chain check (2026-11-18).

One-centre s-type Yukawa/Coulomb tree integrals vs flat-space wavefunction
coefficients of Arkani-Hamed, Benincasa, Postnikov (arXiv:1709.02813).

Identities tested (all with the full d^3r measures, no extra 1/r factors):

  I2(a,b;mu)      = int int d^3r1 d^3r2 e^{-a r1} e^{-b r2} e^{-mu r12}/r12
                  = 16 pi^2  d_a d_b psi2(a,b;mu)
  I3(a,b,c;mu,nu) = int int int d^3r1 d^3r2 d^3r3 e^{-a r1 - b r2 - c r3}
                    e^{-mu r12}/r12 e^{-nu r23}/r23
                  = 64 pi^3  d_a d_c psi3(a,b,c;mu,nu)      (no derivative at the middle vertex)

  psi2 = 1/[(a+b)(a+mu)(b+mu)]
  psi3 = (a+2b+c+mu+nu)/[(a+b+c)(a+mu)(c+nu)(b+mu+nu)(a+b+nu)(b+c+mu)]

Mechanism: angular average of e^{-mu r12}/r12 over the direction of r1 is
G_mu(r1,r2)/(r1 r2), G_mu = (e^{-mu|r1-r2|} - e^{-mu(r1+r2)})/(2 mu), the
half-line Dirichlet Green's function.  A vertex of degree d with density
r^k e^{-x r} carries radial weight r^{k+2-d}.

Output: debug/data/posgeom_chain_check.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction

import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

OUT = os.path.join(os.path.dirname(__file__), "data", "posgeom_chain_check.json")
res: dict = {}

# ---------------------------------------------------------------- symbolic
a, b, c, mu, nu = sp.symbols("a b c mu nu", positive=True)
r1, r2, r3 = sp.symbols("r1 r2 r3", positive=True)

psi2 = 1 / ((a + b) * (a + mu) * (b + mu))
psi3 = (a + 2 * b + c + mu + nu) / (
    (a + b + c) * (a + mu) * (c + nu) * (b + mu + nu) * (a + b + nu) * (b + c + mu))


def say(msg: str) -> None:
    print(msg, flush=True)


def leaf(x: sp.Symbol, m: sp.Symbol) -> dict:
    """int_0^inf e^{-x r} G_m(r, s) dr = (e^{-x s} - e^{-m s})/(m^2 - x^2), as {rate: coeff}.

    By hand (no sp.integrate): the r<s and r>s pieces sum to
    (2m/(m^2-x^2)) (e^{-x s} - e^{-m s}); divided by 2m.  Checked by quadrature below.
    """
    return {x: 1 / (m**2 - x**2), m: -1 / (m**2 - x**2)}


def outer(y: sp.Symbol, *leaves: dict) -> sp.Expr:
    """int_0^inf e^{-y s} prod(leaves)(s) ds, termwise: int e^{-k s} ds = 1/k."""
    terms = [(sp.Integer(0), sp.Integer(1))]
    for lf in leaves:
        terms = [(k + r, cf * c2) for k, cf in terms for r, c2 in lf.items()]
    return sum(cf / (y + k) for k, cf in terms)


half2 = outer(b, leaf(a, mu))
half3 = outer(b, leaf(a, mu), leaf(c, nu))
res["sym_halfline_2site_equals_psi2"] = bool(sp.cancel(sp.together(half2 - psi2)) == 0)
res["sym_halfline_3site_equals_psi3"] = bool(sp.cancel(sp.together(half3 - psi3)) == 0)
say(f"symbolic half-line: {res['sym_halfline_2site_equals_psi2']} {res['sym_halfline_3site_equals_psi3']}")

psi2c = 1 / (a * b * (a + b))
psi3c = (a + 2*b + c) / ((a + b + c) * a * b * c * (a + b) * (b + c))
I2_sym = sp.factor(16 * sp.pi**2 * sp.diff(psi2, a, b))
I3_sym = 64 * sp.pi**3 * sp.together(sp.diff(psi3, a, c))
I2_c = sp.factor(16 * sp.pi**2 * sp.diff(psi2c, a, b))
I3_c = sp.factor(64 * sp.pi**3 * sp.diff(psi3c, a, c))
res["I2_closed"] = str(I2_sym)
res["I2_coulomb"] = str(I2_c)
res["I3_coulomb"] = str(I3_c)
res["I3_denominator_factors"] = str(sp.factor_list(sp.denom(sp.cancel(I3_sym)))[1])
# memo's (pi/4) normalisation
memo = sp.pi * (a**2 + 3*a*b + b**2 + 2*a*mu + 2*b*mu + mu**2) / (2 * (a+b)**3 * (a+mu)**2 * (b+mu)**2)
res["ratio_I2_over_memo_expression"] = str(sp.cancel(I2_sym / memo))
res["ratio_memo_expression_over_dadb_psi2"] = str(sp.cancel(memo / sp.diff(psi2, a, b)))
say("closed forms done")

# ---------------------------------------------------------------- numeric
mp.mp.dps = 20
A, B, C, MU, NU = mp.mpf("1.3"), mp.mpf("0.7"), mp.mpf("2.1"), mp.mpf("0.45"), mp.mpf("0.9")


def G(m, r, s):
    if m == 0:
        return min(r, s)
    return (mp.exp(-m * abs(r - s)) - mp.exp(-m * (r + s))) / (2 * m)


# (i) the angular average itself, by quadrature over cos(theta)
def ang_avg(m, r, s):
    f = lambda t: mp.exp(-m * mp.sqrt(r*r + s*s - 2*r*s*t)) / mp.sqrt(r*r + s*s - 2*r*s*t)
    return mp.quad(f, [-1, 1]) / 2


pts = [(MU, mp.mpf("0.8"), mp.mpf("1.9")), (NU, mp.mpf("2.5"), mp.mpf("0.3")), (mp.mpf(0), mp.mpf("1.1"), mp.mpf("0.6"))]
res["angular_average_max_abs_err"] = mp.nstr(max(abs(ang_avg(m, r, s) - G(m, r, s) / (r * s)) for m, r, s in pts), 5)


def leaf_num(x, m, s, power=1):
    """int_0^inf r^power e^{-x r} G_m(r,s) dr by quadrature (kink at r=s)."""
    return mp.quad(lambda r: r**power * mp.exp(-x * r) * G(m, r, s), [0, s, s + 5, mp.inf])


def I2_num(x, y, m):
    return 16 * mp.pi**2 * mp.quad(lambda s: s * mp.exp(-y * s) * leaf_num(x, m, s), [0, 1, 5, 20, mp.inf])


def I3_num(x, y, z, m, n):
    return 64 * mp.pi**3 * mp.quad(
        lambda s: mp.exp(-y * s) * leaf_num(x, m, s) * leaf_num(z, n, s), [0, 1, 5, 20, mp.inf])


def ev(expr, **kw):
    sub = {a: sp.Float(str(kw.get("a", A)), 40), b: sp.Float(str(kw.get("b", B)), 40),
           c: sp.Float(str(kw.get("c", C)), 40), mu: sp.Float(str(kw.get("mu", MU)), 40),
           nu: sp.Float(str(kw.get("nu", NU)), 40)}
    return mp.mpf(str(sp.N(expr.subs(sub), 40)))


lf_err = abs(leaf_num(A, MU, mp.mpf("1.7"), power=0) - (mp.exp(-A*mp.mpf("1.7")) - mp.exp(-MU*mp.mpf("1.7")))/(MU**2 - A**2))
res["leaf_closed_form_abs_err"] = mp.nstr(lf_err, 5)
say("angular + leaf checks done")
n2, p2 = I2_num(A, B, MU), ev(I2_sym)
res["two_site"] = {"params": "a=1.3 b=0.7 mu=0.45", "quadrature": mp.nstr(n2, 25),
                   "16pi^2 dadb psi2": mp.nstr(p2, 25), "rel_err": mp.nstr(abs(n2 / p2 - 1), 5)}
say("two-site done")
n3, p3 = I3_num(A, B, C, MU, NU), ev(I3_sym)
res["three_site"] = {"params": "a=1.3 b=0.7 c=2.1 mu=0.45 nu=0.9", "quadrature": mp.nstr(n3, 25),
                     "64pi^3 dadc psi3": mp.nstr(p3, 25), "rel_err": mp.nstr(abs(n3 / p3 - 1), 5)}
# wrong-answer controls: derivative at the middle vertex / missing total-energy pole
w1 = ev(64 * sp.pi**3 * sp.diff(psi3, a, b, c))
w2 = ev(64 * sp.pi**3 * sp.diff(psi3 * (a + b + c), a, c))
res["three_site_controls"] = {"dadbdc_rel_err": mp.nstr(abs(n3 / abs(w1) - 1), 5),
                              "no_total_pole_rel_err": mp.nstr(abs(n3 / w2 - 1), 5)}
say("three-site done")
# Coulomb limits by quadrature
n2c, p2c = I2_num(A, B, mp.mpf(0)), ev(I2_c)
n3c, p3c = I3_num(A, B, C, mp.mpf(0), mp.mpf(0)), ev(I3_c)
res["two_site_coulomb_rel_err"] = mp.nstr(abs(n2c / p2c - 1), 5)
res["three_site_coulomb_rel_err"] = mp.nstr(abs(n3c / p3c - 1), 5)

say("coulomb quadrature done")
# ---------------------------------------------------------------- exact corpus values
from geovac.hypergeometric_slater import compute_rk_exact  # noqa: E402

# F0(1s,1s) = 5/8: densities 4 e^{-2r} r^2 each; R^0 = int int r1 r2 P1 P2 G_0
F0_1s1s = sp.Rational(16) * sp.diff(psi2c, a, b).subs({a: 2, b: 2})
# F0(1s,2s) = 17/81: 2s density (1/8)(4 - 4r + r^2) e^{-r}; r^k -> (-d_b)^k
d = lambda e, n: sp.diff(e, b, n) if n else e
F0_1s2s = sp.Rational(4, 8) * sp.diff(-(4 * (-d(psi2c, 1)) - 4 * d(psi2c, 2) + (-d(psi2c, 3))), a).subs({a: 2, b: 1})
q1, s1 = compute_rk_exact(1, 0, 1, 0, 1, 0, 1, 0, 0)
q2, s2 = compute_rk_exact(1, 0, 1, 0, 2, 0, 2, 0, 0)
res["geovac_F0_1s1s"] = {"geovac": str(q1), "sqrt": s1, "from_psi2": str(F0_1s1s), "equal": Fraction(str(F0_1s1s)) == q1 and s1 == 1}
res["geovac_F0_1s2s"] = {"geovac": str(q2), "sqrt": s2, "from_psi2": str(F0_1s2s), "equal": Fraction(str(F0_1s2s)) == q2 and s2 == 1}
# textbook: int int e^{-2Z r1} e^{-2Z r2}/r12 = 5 pi^2/(8 Z^5)
res["textbook_5pi2_over_8Z5"] = bool(sp.simplify(I2_c.subs({a: 2, b: 2}) - sp.Rational(5, 8) * sp.pi**2) == 0)

# Paper 63 one-centre chain (tests/test_paper63_three_electron_k3.py): central zeta first
chain = {}
for zs, want in {(1, 1, 1): sp.Rational(23, 54), (1, 2, 3): sp.Rational(133, 144), (2, 1, 3): sp.Rational(289, 225)}.items():
    zc, zl, zr = zs
    val = sp.Rational(zc * zl * zr)**3 * 64 * sp.diff(psi3c, a, c).subs({b: 2 * zc, a: 2 * zl, c: 2 * zr})
    chain[str(zs)] = {"from_psi3": str(val), "paper63": str(want), "equal": bool(val == want)}
res["paper63_one_centre_chain"] = chain

say("corpus values done")
# ---------------------------------------------------------------- scope: degree-3 vertex gives a log
mp.mp.dps = 30
pot = lambda r: 1 / r - mp.exp(-2 * r) * (1 + 1 / r)
star = mp.quad(lambda r: 4 * r**2 * mp.exp(-2 * r) * pot(r)**3, [0, 1, 4, 12, mp.inf])
rel = mp.pslq([star, 1, mp.log(2), mp.log(3)], maxcoeff=10**5, maxsteps=10**4)
rel_rat = mp.pslq([star, 1], maxcoeff=10**5, maxsteps=10**4)
res["star4_1s_unit_zeta"] = {"value": mp.nstr(star, 30), "pslq_[star,1,log2,log3]": str(rel), "pslq_[star,1]": str(rel_rat),
                             "note": "central 1s vertex of degree 3 has radial weight r^{-1}: not rational"}

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    json.dump(res, fh, indent=1)
print(json.dumps(res, indent=1))
