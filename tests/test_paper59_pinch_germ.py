"""Guard for a sentence PROPOSED for Paper 59 sec:scope (not in the paper when this file was written, 2026-11-18):

    "Both steps are accounted for by one mechanism, checked exactly but not proved: the half-integer local solution
    of each endpoint piece at x* agrees with the endpoint part, at the sphere through C or through A, of the
    Landau-pinch germ of the potential that piece sees, computed to one order beyond the leading term. ..."

Written as a separate step (CLAUDE.md sec. 9, guard-writing rule).  Two sides, kept apart in the code:

(A) PREDICTION (section A below): exact rationals from the geometry (D1, D2, zB, zC) and the potential's odd Taylor
    coefficients only.  No operator, no prime, no embedded file is read by this side.
        I_pq ~ C_pq |a|^{p+q+3} int_0^1 u^{(p+1)/2} (1-u)^{(q+1)/2}
                 [ g(R0) + u(1-u)|a|^2 ( zB/(2 R0) g(R0) - lap g / (2 (p+q+5)) ) ] du,     g = e^{-zB R},
    R0^2 = (1-u) D2^2 + u D1^2, |a|^2 = 2 D1 D2 (x* - x); u = 1 is C (sphere of radius D1 about B, "D1"), u = 0 is A
    ("D2").  With r = R0 the u-integral is F(D2) - F(D1), F = e^{-zB r} P(r); the D1 part is -P(D1), the D2 part P(D2).
(B) MEASUREMENT (section B): for each endpoint piece operator mod p, the local solution sum_n a_n t^{1/2+n},
    t = x* - x, from the null space of the truncated (lower-triangular) linear system L y = 0, by this file's own
    elimination; normalised to leading coefficient 1.  Not the memo's recursion, not s1_controls.points.

COVERAGE, exactly (object, geometry, prime):
    f' = NA - Y_mu   A1  1073741441, 1073741381      A4  1073741441, 1073741381      (d5, d7, rho_D2; + d3 of NA)
    (AA|BC)          A1  1073741441, 1073741381      A2  1073741441, 1073741381      (d5 = 0, d7, d9, rho_D2)
The f' piece operators are those of tests/data/paper59_op52_2026_11_18.json.  Record MO's embedded file holds only the
FULL (AA|BC) operators, so the five endpoint piece operators of (AA|BC) are carried in
tests/data/paper59_pinch_germ_2026_11_18.json, built 2026-11-18 from the scored MO run's cached sequences with the
frozen guesser copy; BOTH primes at both geometries (the derivation memo ran one prime per geometry).  That file is
not trusted on its own: each piece operator must annihilate its sequence (T-form route), record MO's embedded full
operator must annihilate each of the five sequences, and the first terms are recomputed from the frozen engine copy.
A1 is the same geometry for both objects: three geometries, four (object, geometry) cases.

Claims guarded, each with the wrong answer it rejects:
(1) the next-order germ formula: Schwinger weights, the segment-point shift, the isotropic Hessian term, the
    canonical endpoint split (sympy; rejects 2(p+q+3), a shift of the other sign, u <-> 1-u);
(2) every endpoint piece has exactly ONE half-integer indicial root at x*, simple, and a one-dimensional log-free
    half-integer solution space (rejects two branches per piece);
(3) the two halves of one family are two non-proportional sequences with ONE minimal operator, hence one branch
    series (rejects a branch that depends on the exponential constant).  Found while writing: 'the halves have
    identical series' is this operator equality, not a separate coincidence of two series;
(4) the predicted rationals equal the measured coefficients mod p (rejects the leading-only germ, the D2 endpoint
    part used at D1, another geometry's germ, the wrong potential, a wrong prime);
(5) rho_D2 - rho_D1 != 0, predicted exactly and measured, at the four cases (rejects a zero residual);
(6) block and full operators: log-free at x* (Frobenius: every resonant right-hand side vanishes; a detector
    control with a logarithm is read as such) and their half-integer space is the span of the pieces' series;
(7) one branch per piece is not automatic: (x* - x)^{1/2} log(x - x1) has exponents {1/2, 3/2} at x*.
NOT guarded, the sentence's own cap: the operators are guessed and read mod p; the identification is checked to one
order beyond leading only; one-dimensionality per piece is measured, not proved.

Fire-test record: the comment block at the end of this file.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import flint
import numpy as np
import pytest
import sympy as sp

HERE = Path(__file__).parent
SUP = HERE / 'paper59_intclass_support'
EMBED = HERE / 'data' / 'paper59_op52_2026_11_18.json'          # f': sequences, piece / block / full operators
MO_EMBED = HERE / 'data' / 'paper59_mo_2026_11_17.json'         # (AA|BC): full operators only
PG_EMBED = HERE / 'data' / 'paper59_pinch_germ_2026_11_18.json'  # (AA|BC): the five endpoint pieces

PRIMES = (1073741441, 1073741381)
SEEDS = (2026, 7)
# D1 = |BC|, D2 = |BA|, zB, zC, zA
GEOMETRY = {'A1': ('9/5', '13/5', '13/10', '1', '9/10'),
            'A4': ('8/5', '5/2', '13/10', '4/5', '7/10'),
            'A2': ('3/2', '12/5', '1', '6/5', '4/5')}
XSTAR = {'A1': Fraction(125, 117), 'A4': Fraction(881, 800), 'A2': Fraction(89, 80)}
# what each endpoint's pieces see: D1 = (the lone piece, the two-halves family), D2 = the pair
OBJECTS = {'FPRIME': {'cases': ('A1', 'A4'), 'D1': ('point', 'yukawa'), 'D2': 'vprime'},
           'AABC': {'cases': ('A1', 'A2'), 'D1': ('point', 'screened'), 'D2': 'bounded'}}
OCP = [(o, c, p) for o in OBJECTS for c in OBJECTS[o]['cases'] for p in PRIMES]
OC = [(o, c) for o in OBJECTS for c in OBJECTS[o]['cases']]
NTERMS = 12


def params(case):
    return tuple(Fraction(v) for v in GEOMETRY[case])


# ======================================================================================================================
# (A) PREDICTION: geometry + potential data, exact rationals.  Nothing below this banner, down to section (B), reads
#     an operator or a prime.
# ======================================================================================================================
_s, _mu, _z = sp.symbols('s mu zeta', positive=True)
POT = {'point': 1 / _s,                                           # NA's kernel
       'yukawa': sp.exp(-_mu * _s) / _s,                          # Y_mu, mu = 2 zA
       'vprime': (1 - sp.exp(-_mu * _s)) / _s,                    # potential of mu^2 e^{-mu s}/(4 pi s)
       'screened': (1 / _s + _z) * sp.exp(-2 * _z * _s),          # S, zeta = zA
       'bounded': 1 / _s - (1 / _s + _z) * sp.exp(-2 * _z * _s)}  # potential of the 1s^2 density
ODD = (-1, 1, 3, 5)


@lru_cache(maxsize=None)
def _odd_symbolic(name):
    ser = sp.series(POT[name], _s, 0, ODD[-1] + 2).removeO()
    return tuple(sp.simplify(ser.coeff(_s, n)) for n in ODD)


def odd_coefficients(name, zA):
    """{p: v_p} for the odd powers s^p, p = -1..5, of the potential at mu = 2 zA, zeta = zA."""
    zA = sp.Rational(zA.numerator, zA.denominator)
    out = {}
    for n, c in zip(ODD, _odd_symbolic(name)):
        v = sp.Rational(c.subs({_mu: 2 * zA, _z: zA}))
        out[n] = Fraction(int(v.p), int(v.q))
    return out


def gamma_half(k):
    """Gamma(k/2) / sqrt(pi) for odd k, exact."""
    assert k % 2
    if k == 1:
        return Fraction(1)
    if k > 1:
        return Fraction(k - 2, 2) * gamma_half(k - 2)
    return gamma_half(k + 2) / Fraction(k, 2)


def pinch_constant(p, q):
    """C_pq / pi,  C_pq = pi^{3/2} Gamma(-(p+q+3)/2) / (Gamma(-p/2) Gamma(-q/2)),  p, q odd."""
    return gamma_half(-(p + q + 3)) / (gamma_half(-p) * gamma_half(-q))


def correction_bracket(p, q, zB):
    """The relative-order |a|^2 bracket multiplying u (1-u) |a|^2 g(R0), as (constant, coefficient of 1/R0):
         zB / (2 R0)                       segment point: R^2 = R0^2 - u (1-u) |a|^2,  g = e^{-zB R}
       - (zB^2 - 2 zB / R0) / (2 (p+q+5))  isotropic Hessian: - lap g / (2 (p+q+5))."""
    shift = (Fraction(0), zB / 2)
    hess = (-zB ** 2 / (2 * (p + q + 5)), 2 * zB / (2 * (p + q + 5)))
    return shift[0] + hess[0], shift[1] + hess[1]


# ---- polynomials in r, ascending lists of Fractions
def pmul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def ppow(a, n):
    out = [Fraction(1)]
    for _ in range(n):
        out = pmul(out, a)
    return out


def pval(a, x):
    return sum((c * x ** i for i, c in enumerate(a)), Fraction(0))


def pdiff(a):
    return [i * c for i, c in enumerate(a)][1:]


def segment_u(D1, D2):
    """u as a polynomial in r = R0 on the segment: u = (D2^2 - r^2) / (D2^2 - D1^2); u = 1 at r = D1 (the point C),
    u = 0 at r = D2 (the point A)."""
    Dl = D2 ** 2 - D1 ** 2
    return [D2 ** 2 / Dl, Fraction(0), -1 / Dl]


def antiderivative(poly, zB):
    """The polynomial P with d/dr [e^{-zB r} P] = e^{-zB r} poly, i.e. P' - zB P = poly (unique)."""
    P, d, k = [Fraction(0)] * len(poly), list(poly), 0
    while d:
        for i, c in enumerate(d):
            P[i] -= c / zB ** (k + 1)
        d, k = pdiff(d), k + 1
    return P


def endpoint_parts(a, b, extra, D1, D2, zB):
    """int_0^1 u^a (1-u)^b extra(r) e^{-zB R0} du with r = R0, du = -2 r dr / (D2^2 - D1^2):
    = [e^{-zB r} P]_{D1}^{D2}.  Returns the two endpoint parts with their exponentials removed.
    extra = (c0, c1) stands for c0 + c1 / r."""
    Dl = D2 ** 2 - D1 ** 2
    u = segment_u(D1, D2)
    w = pmul(ppow(u, a), ppow([1 - u[0], -u[1], -u[2]], b))
    integrand = pmul(w, [2 * extra[1] / Dl, 2 * extra[0] / Dl])          # (c0 + c1/r) * 2 r / Dl
    P = antiderivative(integrand, zB)
    return {'D1': -pval(P, D1), 'D2': pval(P, D2)}


def germ(v, D1, D2, zB, zC, correction=True):
    """Endpoint parts of sum_{p,q} v_p (-zC)^q/q! I_pq, as {E: {n: coefficient of t^{n/2}}} for the leading odd n0
    and n0 + 2 (the common factor N pi (2 D1 D2)^{3/2} dropped).  Returns (n0, parts)."""
    K = 2 * D1 * D2
    n0 = min(p for p, vp in v.items() if vp) + 1 + 3
    out = {E: {n0: Fraction(0), n0 + 2: Fraction(0)} for E in ('D1', 'D2')}
    for p, vp in v.items():
        if not vp:
            continue
        for q in (1, 3, 5):
            n = p + q + 3
            if n > n0 + 2:
                continue
            c = vp * pinch_constant(p, q) * (-zC) ** q / sp.factorial(q)
            a, b = (p + 1) // 2, (q + 1) // 2
            lead = endpoint_parts(a, b, (Fraction(1), Fraction(0)), D1, D2, zB)
            for E in out:
                out[E][n] += c * K ** ((n - 3) // 2) * lead[E]
            if correction and n + 2 <= n0 + 2:
                corr = endpoint_parts(a + 1, b + 1, correction_bracket(p, q, zB), D1, D2, zB)
                for E in out:
                    out[E][n + 2] += c * K ** ((n + 2 - 3) // 2) * corr[E]
    return n0, out


def predictions(obj, case, correction=True, d1_endpoint='D1', force_zero_residual=False):
    """The parameter-free predictions for one (object, geometry):
      lead   leading exponents of the D1 pieces and of the D2 pieces;
      d3     [t^{5/2}] of the normalised series of the lone point-charge piece at D1 (B_N = t^{3/2} (1 + d3 t ...));
      d      {n: [t^{n/2}] of B_family - B_N} for the two orders the difference germ is known to;
      rho_D1, rho_D2   next / leading coefficient of the germ of the D2 potential at each endpoint;
      resid  rho_D2 - rho_D1."""
    D1, D2, zB, zC, zA = params(case)
    lone, fam = OBJECTS[obj]['D1']
    vl, vf, v2 = (odd_coefficients(nm, zA) for nm in (lone, fam, OBJECTS[obj]['D2']))
    diff = {p: vf[p] - vl[p] for p in ODD}
    nl, gl = germ(vl, D1, D2, zB, zC, correction)
    nf, gf = germ(vf, D1, D2, zB, zC, correction)
    nd, gd = germ(diff, D1, D2, zB, zC, correction)
    n2, g2 = germ(v2, D1, D2, zB, zC, correction)
    E = d1_endpoint
    lead = gl[E][nl]
    out = {'lead': {'D1': (Fraction(nl, 2), Fraction(nf, 2)), 'D2': Fraction(n2, 2)},
           'd3': gl[E][nl + 2] / lead,
           'd': {n: (gd[E][n] / lead if n in gd[E] else Fraction(0)) for n in range(5, nd + 3, 2)},
           'nd': nd, 'n2': n2,
           'rho_D1': gd[E][nd + 2] / gd[E][nd],
           'rho_D2': g2['D2'][n2 + 2] / g2['D2'][n2]}
    if force_zero_residual:
        out['rho_D1'] = out['rho_D2']
    out['resid'] = out['rho_D2'] - out['rho_D1']
    return out


# typed from the derivation memo's data files (debug/data/op52_deriv_{nextorder,aabc_*}.json), produced by another
# author's sympy route (sp.integrate); not read from any file here
TYPED = {
    ('FPRIME', 'A1'): {'d5': '-55222849332/5770496875', 'd7': '-8363055453905000997/79423676363281250',
                       'rho_D1': '11471955355150893/1042624818235000', 'rho_D2': '325288437828681/30021182916725',
                       'resid': '-133315426199032841829/795157299960054955000'},
    ('FPRIME', 'A4'): {'d5': '-294034080676/66094864875', 'd7': '-473335076843648079592/12674356245544865625',
                       'rho_D1': '16904824172987431414/2013714191535207525',
                       'rho_D2': '28528997043499024/3467861873368425',
                       'resid': '-32734421418205156806/194655690989567650325'},
    ('AABC', 'A1'): {'d5': '0', 'd7': '-1220992024002715071/14440668429687500',
                     'd9': '-26799622745914812579083307/28393964299873046875000',
                     'rho_D1': '50428218270541438427/4517482782087453750',
                     'rho_D2': '342361942152264847/30637470212302525',
                     'resid': '288281497059030345777797657/24611871515529793041349646250'},
    ('AABC', 'A2'): {'d5': '0', 'd7': '-300547039232/3604453125', 'd9': '-5090388295233437696/6167669853515625',
                     'rho_D1': '621385289945488/62777533264875', 'rho_D2': '119400025449416/12059969644125',
                     'resid': '194396849147304/84036963070764275'}}


def test_next_order_germ_formula_is_derived():
    """Claim (1).  Schwinger form |y|^p = Gamma(-p/2)^{-1} int alpha^{-p/2-1} e^{-alpha y^2} d alpha for the potential
    term (at A) and the cusp term (at C); alpha = lam (1-u), beta = lam u.  Steps asserted:
      a. alpha |r-A|^2 + beta |r-C|^2 = lam |r - r_u|^2 + lam u (1-u) |A-C|^2,  r_u = (1-u) A + u C;
      b. |r_u|^2 = (1-u) D2^2 + u D1^2 - u (1-u) |a|^2 and |a|^2 = 2 D1 D2 (x* - x)   (B at the origin);
      c. Jacobian lam; with the Gaussian (pi/lam)^{3/2} the lam power is -(nu+1), nu = (p+q+3)/2, and after
         int lam^{-nu-1} e^{-lam M} = Gamma(-nu) M^nu the weights are u^{(p+1)/2} (1-u)^{(q+1)/2};
      d. the Gaussian average of g is g + lap g / (4 lam): one more 1/lam, so the factor is
         Gamma(-nu-1) / (4 Gamma(-nu)) = -1 / (2 (p+q+5)) per u (1-u) |a|^2;
      e. lap e^{-zB R} = (zB^2 - 2 zB / R) e^{-zB R};  d/d eps e^{-zB sqrt(R0^2 - eps)} at 0 = zB / (2 R0) g;
      f. C_pq / pi is the rational used.
    Wrong answers rejected: 2 (p+q+3) in the Hessian term; the segment shift with the other sign; weights with p and
    q exchanged; a pinch constant off by a rational."""
    lam, u, M, R, eps, zB, y = sp.symbols('lam u M R eps zB y', positive=True)
    p, q = sp.symbols('p q')
    # a
    r, A, C = (sp.Matrix(sp.symbols(f'{n}1:4', real=True)) for n in 'rAC')
    ru = (1 - u) * A + u * C
    sq = lambda v: (v.T * v)[0]                                                    # noqa: E731
    assert sp.expand(lam * (1 - u) * sq(r - A) + lam * u * sq(r - C) - lam * sq(r - ru)
                     - lam * u * (1 - u) * sq(A - C)) == 0
    # b
    D1, D2, x = sp.symbols('D1 D2 x', positive=True)
    AA, CC, AC = D2 ** 2, D1 ** 2, D1 * D2 * x                                     # |A|^2, |C|^2, A.C
    a2 = AA + CC - 2 * AC
    assert sp.expand((1 - u) ** 2 * AA + u ** 2 * CC + 2 * u * (1 - u) * AC
                     - ((1 - u) * AA + u * CC - u * (1 - u) * a2)) == 0
    assert sp.simplify(a2 - 2 * D1 * D2 * ((D1 ** 2 + D2 ** 2) / (2 * D1 * D2) - x)) == 0
    # c
    al, be = lam * (1 - u), lam * u
    jac = sp.Matrix([[sp.diff(al, lam), sp.diff(al, u)], [sp.diff(be, lam), sp.diff(be, u)]]).det()
    assert sp.simplify(jac - lam) == 0
    nu = (p + q + 3) / 2
    lam_power = (-p / 2 - 1) + (-q / 2 - 1) + 1 - sp.Rational(3, 2)
    assert sp.simplify(lam_power + nu + 1) == 0
    assert sp.simplify(nu + (-q / 2 - 1) - (p + 1) / 2) == 0                       # exponent of u
    assert sp.simplify(nu + (-p / 2 - 1) - (q + 1) / 2) == 0                       # exponent of 1 - u
    sgm = sp.Symbol('sigma', positive=True)
    mellin = sp.integrate(lam ** (sgm - 1) * sp.exp(-lam * M), (lam, 0, sp.oo), conds='none')
    assert sp.simplify(mellin - sp.gamma(sgm) * M ** (-sgm)) == 0
    assert sp.simplify(sp.integrate(lam ** sp.Rational(-1, 2) * sp.exp(-lam * y ** 2), (lam, 0, sp.oo))
                       / sp.gamma(sp.Rational(1, 2)) - 1 / y) == 0                 # the Schwinger form at p = -1
    # d
    second = sp.integrate(y ** 2 * sp.exp(-lam * y ** 2), (y, -sp.oo, sp.oo)) / sp.integrate(
        sp.exp(-lam * y ** 2), (y, -sp.oo, sp.oo))
    assert sp.simplify(second - 1 / (2 * lam)) == 0                                # <y_i y_j> = delta_ij / (2 lam)
    w = sp.Symbol('w')
    ratio = sp.gammasimp(sp.gamma(w) / sp.gamma(w + 1))                            # Gamma(-nu-1) / Gamma(-nu)
    assert sp.simplify(ratio - 1 / w) == 0
    hess_factor = ratio.subs(w, -nu - 1) * sp.Rational(1, 2) * sp.Rational(1, 2)   # (1/2) d_i d_j g <y_i y_j>
    assert sp.simplify(hess_factor + 1 / (2 * (p + q + 5))) == 0
    # e
    g = sp.exp(-zB * R)
    lap = sp.simplify((sp.diff(R ** 2 * sp.diff(g, R), R) / R ** 2) / g)
    assert sp.simplify(lap - (zB ** 2 - 2 * zB / R)) == 0
    shift = sp.simplify(sp.diff(sp.exp(-zB * sp.sqrt(R ** 2 - eps)), eps).subs(eps, 0) / g)
    assert sp.simplify(shift - zB / (2 * R)) == 0          # R^2 = R0^2 - eps, eps = u (1-u) |a|^2: g grows
    # the literals of correction_bracket against d and e, R^2 = R0^2 - u (1-u) |a|^2
    for pp, qq in ((-1, 1), (1, 1), (3, 1), (1, 3), (5, 3)):
        zb = Fraction(13, 10)
        c0, c1 = correction_bracket(pp, qq, zb)
        want = (shift + hess_factor * lap).subs({p: pp, q: qq, zB: sp.Rational(13, 10)})
        got = sp.Rational(c0.numerator, c0.denominator) + sp.Rational(c1.numerator, c1.denominator) / R
        assert sp.simplify(want - got) == 0, (pp, qq)
    # f
    for pp, qq in ((-1, 1), (1, 1), (3, 1), (1, 3), (-1, 3), (5, 1), (3, 3)):
        want = sp.pi ** sp.Rational(3, 2) * sp.gamma(-sp.Rational(pp + qq + 3, 2)) / (
            sp.gamma(-sp.Rational(pp, 2)) * sp.gamma(-sp.Rational(qq, 2))) / sp.pi
        c = pinch_constant(pp, qq)
        assert sp.simplify(want - sp.Rational(c.numerator, c.denominator)) == 0, (pp, qq)
    assert pinch_constant(1, 1) / pinch_constant(-1, 1) == Fraction(1, 5)


def test_endpoint_split_orientation_and_uniqueness():
    """Claim (1), the split.  u = 1 at r = D1 (C's sphere) and u = 0 at r = D2 (A's sphere); the polynomial P is the
    one sympy's integral gives (e^{-zB r} P is the antiderivative; a polynomial solution of P' - zB P = poly is
    unique because the homogeneous solutions are e^{zB r}); the two parts add up to the segment integral
    (quadrature); nothing but e^{-zB r} poly arises.  Wrong answers rejected: u <-> 1-u; a split that does not sum
    to the integral; an antiderivative off by a polynomial."""
    D1, D2, zB = Fraction(9, 5), Fraction(13, 5), Fraction(13, 10)
    u = segment_u(D1, D2)
    assert pval(u, D1) == 1 and pval(u, D2) == 0
    for r0 in (D1, Fraction(2), D2):                                 # R0^2 = (1-u) D2^2 + u D1^2
        uu = pval(u, r0)
        assert (1 - uu) * D2 ** 2 + uu * D1 ** 2 == r0 ** 2
    r = sp.Symbol('r', positive=True)
    poly = [Fraction(3, 7), Fraction(-2), Fraction(0), Fraction(5, 3), Fraction(1, 11)]
    P = antiderivative(poly, zB)
    as_sym = lambda a: sum(sp.Rational(c.numerator, c.denominator) * r ** i for i, c in enumerate(a))   # noqa: E731
    zb = sp.Rational(13, 10)
    assert sp.simplify(sp.diff(sp.exp(-zb * r) * as_sym(P), r) - sp.exp(-zb * r) * as_sym(poly)) == 0
    F = sp.integrate(as_sym(poly) * sp.exp(-zb * r), r)
    assert sp.simplify(sp.expand(F * sp.exp(zb * r)) - as_sym(P)) == 0
    import mpmath as mp
    with mp.workdps(30):
        for a, b, extra in ((0, 1, (Fraction(1), Fraction(0))), (2, 2, correction_bracket(1, 1, zB)),
                            (3, 2, correction_bracket(3, 1, zB))):
            parts = endpoint_parts(a, b, extra, D1, D2, zB)
            f = lambda t: (t ** a * (1 - t) ** b * (float(extra[0]) + float(extra[1]) / mp.sqrt(   # noqa: E731
                (1 - t) * float(D2) ** 2 + t * float(D1) ** 2)) * mp.exp(-float(zB) * mp.sqrt(
                    (1 - t) * float(D2) ** 2 + t * float(D1) ** 2)))
            quad = mp.quad(f, [0, 0.5, 1])
            total = float(parts['D1']) * mp.exp(-float(zB * D1)) + float(parts['D2']) * mp.exp(-float(zB * D2))
            assert abs(total / quad - 1) < 1e-10, (a, b)
            assert parts['D1'] != parts['D2'] and parts['D1'] != 0


def test_potentials_seen_by_the_pieces():
    """The potentials: the two D1 potentials of each object share the s^-1 coefficient 1, and their difference is
    minus the D2 potential (so 'the two endpoints' germs of one potential' is meant literally); V' leads with s^1
    (-mu^2/2), the bounded potential with s^3 (2 zeta^4 / 3) and has no s^1 term.  Wrong answers rejected: a D2
    potential that is not the difference; an s^1 term in the bounded potential."""
    for obj, o in OBJECTS.items():
        lone, fam = o['D1']
        assert sp.simplify(POT[fam] - POT[lone] + POT[o['D2']]) == 0
        zA = Fraction(9, 10)
        vl, vf, v2 = (odd_coefficients(nm, zA) for nm in (lone, fam, o['D2']))
        assert vl[-1] == vf[-1] == 1 and v2[-1] == 0
        assert all(vf[p] - vl[p] == -v2[p] for p in ODD)
    v = odd_coefficients('vprime', Fraction(9, 10))
    assert v[1] == -Fraction(9, 5) ** 2 / 2 and v[3] == -Fraction(9, 5) ** 4 / 24
    b = odd_coefficients('bounded', Fraction(4, 5))
    assert b[1] == 0 and b[3] == Fraction(2, 3) * Fraction(4, 5) ** 4 and b[5] != 0
    assert odd_coefficients('point', Fraction(1)) == {-1: 1, 1: 0, 3: 0, 5: 0}


@pytest.mark.parametrize('obj,case', OC)
def test_predicted_rationals(obj, case):
    """Side (A) alone: the predictions equal the rationals typed from the memo's data (another author's route:
    sp.integrate for the endpoint parts, sympy gamma for the constants); the leading exponents are 3/2 for the D1
    pieces and 5/2 (f') or 7/2 ((AA|BC)) for the D2 pieces; the D1 difference starts at 5/2 for f' and at 7/2 for
    (AA|BC), where [t^{5/2}] is exactly 0; rho_D2 - rho_D1 is not 0 (the '+1' condition).  Wrong answers rejected:
    a prediction code that drifts; a zero residual; a (AA|BC) difference that starts at 5/2."""
    pr, ty = predictions(obj, case), TYPED[(obj, case)]
    assert pr['lead']['D1'] == (Fraction(3, 2), Fraction(3, 2))
    assert pr['lead']['D2'] == (Fraction(5, 2) if obj == 'FPRIME' else Fraction(7, 2))
    assert pr['nd'] == pr['n2'] == (5 if obj == 'FPRIME' else 7)
    for key in ('d5', 'd7', 'd9'):
        if key in ty:
            assert pr['d'][int(key[1])] == Fraction(ty[key]), key
    assert sorted(pr['d']) == ([5, 7] if obj == 'FPRIME' else [5, 7, 9])
    for key in ('rho_D1', 'rho_D2', 'resid'):
        assert pr[key] == Fraction(ty[key]), key
    assert pr['resid'] != 0 and pr['rho_D1'] != pr['rho_D2']
    assert pr['d'][pr['nd']] != 0
    if obj == 'AABC':
        assert pr['d'][5] == 0


def test_prediction_controls_differ_from_the_prediction():
    """The wrong answers used as negative controls below are different numbers, exactly (so 'fails to match' cannot
    be a control that equals the prediction).  At A1, f': leading-only rho_D2, the D2 part at D1, and the
    forced-zero residual."""
    pr = predictions('FPRIME', 'A1')
    lo = predictions('FPRIME', 'A1', correction=False)
    sw = predictions('FPRIME', 'A1', d1_endpoint='D2')
    assert lo['rho_D2'] != pr['rho_D2'] and lo['d'][7] != pr['d'][7] and lo['d'][5] == pr['d'][5]
    assert sw['d'][5] == Fraction(-477023536746, 41997616375) != pr['d'][5]       # typed from the memo's data
    assert predictions('FPRIME', 'A1', force_zero_residual=True)['resid'] == 0


# ======================================================================================================================
# (B) MEASUREMENT: local solutions of the embedded operators mod p.  Nothing below reads a potential or calls germ().
# ======================================================================================================================
def fp(q, p):
    q = Fraction(q)
    assert q.denominator % p, 'denominator vanishes mod p'
    return q.numerator % p * pow(q.denominator % p, p - 2, p) % p


def small_rational(a, p, dmax=64, nmax=4096):
    for d in range(1, dmax + 1):
        n = a * d % p
        if n > p // 2:
            n -= p
        if abs(n) <= nmax:
            return Fraction(n, d)
    return None


@lru_cache(maxsize=None)
def _json(path):
    return json.loads(Path(path).read_bytes().replace(b'\r\n', b'\n').decode())


def entry(obj, case, p):
    return _json(str(EMBED if obj == 'FPRIME' else PG_EMBED))['entries'][f'{case}|{p}']


def endpoint_keys(case):
    """repr of the engine's keys of the endpoint constants: D1 = e^{-zB D1} {1, e^{-mu (D2 -+ D1)}},
    D2 = e^{+-zC D1 - (zB + zC) D2}, mu = 2 zA (both objects)."""
    D1, D2, zB, zC, zA = params(case)
    mu = 2 * zA
    one = lambda e: repr((Fraction(e), '1', None))                                # noqa: E731
    return {'lone': one(-zB * D1), 'fam-': one(-zB * D1 - mu * (D2 - D1)), 'fam+': one(-zB * D1 - mu * (D2 + D1)),
            'D2+': one(zC * D1 - (zB + zC) * D2), 'D2-': one(-zC * D1 - (zB + zC) * D2)}


def expand_at(c_lists, xs, p):
    """L = sum_k c_k(x) (d/dx)^k in t = x* - x: rows C[k] = coefficients in t of (-1)^k c_k(x* - t)."""
    sub = flint.nmod_poly([fp(xs, p), p - 1], p)
    C = []
    for k, co in enumerate(c_lists):
        P = flint.nmod_poly([int(v) for v in co], p)
        e = [int(v) for v in P(sub).coeffs()] if P.degree() >= 0 else []
        C.append([(-v) % p for v in e] if k % 2 else e)
    return C


def _ff(m, k, p):
    out = 1
    for j in range(k):
        out = out * ((m - j) % p) % p
    return out


def lowest(C):
    """sigma = min_k (valuation of C[k] - k): L t^m = t^{m + sigma} (chi(m) + O(t))."""
    return min(next(j for j, v in enumerate(row) if v) - k for k, row in enumerate(C) if any(row))


def indicial_roots(C, p):
    """All roots of the indicial polynomial at x* as Fractions (asserts regular singular and that every root is a
    small rational)."""
    sig, R = lowest(C), len(C) - 1
    chi = flint.nmod_poly([], p)
    for k, row in enumerate(C):
        j = sig + k
        if 0 <= j < len(row) and row[j]:
            fall = flint.nmod_poly([1], p)
            for i in range(k):
                fall = fall * flint.nmod_poly([(-i) % p, 1], p)
            chi = chi + fall * row[j]
    assert chi.degree() == R, 'irregular point'
    out = []
    for r, mult in chi.roots():
        q = small_rational(int(r), p)
        assert q is not None, int(r)
        out += [q] * int(mult)
    assert len(out) == R, 'roots unaccounted for'
    return sorted(out)


def local_matrix(C, p, N):
    """M[n][i] = coefficient of t^{1/2 + n + sigma} in L t^{1/2 + i}; lower triangular; rows n < N are the complete
    conditions on a_0..a_{N-1}."""
    sig, half = lowest(C), fp(Fraction(1, 2), p)
    M = [[0] * N for _ in range(N)]
    for i in range(N):
        m = (half + i) % p
        for k, row in enumerate(C):
            f = _ff(m, k, p)
            if not f:
                continue
            for n in range(i, N):
                j = n - i + sig + k
                if 0 <= j < len(row) and row[j]:
                    M[n][i] = (M[n][i] + row[j] * f) % p
    assert all(M[n][i] == 0 for n in range(N) for i in range(n + 1, N))
    return M


def rref(rows, p):
    """Reduced row echelon form mod p (pivot = first nonzero column, normalised to 1); returns (rows, pivots)."""
    A, piv, r = [list(v) for v in rows], [], 0
    ncol = len(A[0]) if A else 0
    for col in range(ncol):
        pr = next((i for i in range(r, len(A)) if A[i][col]), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        inv = pow(A[r][col], p - 2, p)
        A[r] = [v * inv % p for v in A[r]]
        for i in range(len(A)):
            if i != r and A[i][col]:
                f = A[i][col]
                A[i] = [(v - f * w) % p for v, w in zip(A[i], A[r])]
        piv.append(col)
        r += 1
    return A[:r], piv


def half_integer_space(c_lists, xs, p, N=NTERMS):
    """Reduced echelon basis (by valuation) of {a : M a = 0}: the log-free solutions sum a_n t^{1/2+n} to N terms.
    Row pivots n are the exponents 1/2 + n."""
    C = expand_at(c_lists, xs, p)
    M = local_matrix(C, p, N)
    R, piv = rref(M, p)
    free = [c for c in range(N) if c not in piv]
    null = []
    for fc in free:
        v = [0] * N
        v[fc] = 1
        for row, pc in zip(R, piv):
            v[pc] = (-row[fc]) % p
        null.append(v)
    for v in null:
        assert all(sum(a * b for a, b in zip(row, v)) % p == 0 for row in M)
    basis, pivots = rref(null, p) if null else ([], [])
    return basis, [Fraction(1, 2) + n for n in pivots]


def frobenius_obstructions(c_lists, xs, p, N=NTERMS):
    """The Frobenius criterion, explicitly.  For each half-integer root e_a = 1/2 + n_a run the recursion from
    a_{n_a} = 1 with every later free coefficient set to 0, and at each later root n_b record the right-hand side
    sum_{i < n_b} M[n_b][i] a_i that a_{n_b} cannot absorb (its diagonal entry is 0).  No logarithm iff all vanish.
    Returns {(e_a, e_b): rhs}."""
    C = expand_at(c_lists, xs, p)
    M = local_matrix(C, p, N)
    roots = [n for n in range(N) if M[n][n] == 0]
    out = {}
    for na in roots:
        a = [0] * N
        a[na] = 1
        for n in range(na + 1, N):
            rhs = sum(M[n][i] * a[i] for i in range(n)) % p
            if M[n][n]:
                a[n] = (-rhs) * pow(M[n][n], p - 2, p) % p
            else:
                out[(Fraction(1, 2) + na, Fraction(1, 2) + n)] = rhs
    return roots, out


def halves_of(roots):
    return [r for r in roots if r.denominator == 2]


@lru_cache(maxsize=None)
def branch_series(obj, case, p, N=NTERMS):
    """{piece: (exponent, normalised series from that exponent on)} for the five endpoint pieces; asserts each has
    a one-dimensional half-integer space."""
    e, keys = entry(obj, case, p), endpoint_keys(case)
    out = {}
    for nm, key in keys.items():
        basis, exps = half_integer_space(e['piece_ops'][key], XSTAR[case], p, N)
        assert len(basis) == 1, (nm, len(basis))
        i0 = next(i for i, v in enumerate(basis[0]) if v)
        assert basis[0][i0] == 1 and exps == [Fraction(1, 2) + i0]
        out[nm] = (exps[0], tuple(basis[0][i0:]))
    return out


def measured(obj, case, p, e=None):
    """The measured counterparts of predictions(): d3, d = {n: [t^{n/2}] (B_family - B_lone)}, rho_D1, rho_D2,
    resid, all mod p."""
    ser = branch_series(obj, case, p) if e is None else e
    el, bl = ser['lone']
    ef, bf = ser['fam-']
    e2, b2 = ser['D2+']
    assert el == ef == Fraction(3, 2)
    d = {3 + 2 * j: (bf[j] - bl[j]) % p for j in range(1, 4)}                 # index j <-> t^{3/2 + j}
    nd = next(n for n in sorted(d) if d[n])
    out = {'lead': {'D1': (el, ef), 'D2': e2}, 'd3': bl[1], 'd': d, 'nd': nd,
           'rho_D1': d[nd + 2] * pow(d[nd], p - 2, p) % p, 'rho_D2': b2[1]}
    out['resid'] = (out['rho_D2'] - out['rho_D1']) % p
    return out


# ------------------------------------------------------------------ data integrity (T-form; not the search's code)
def T_form(c, p):
    """(1 - x^2)^R L = sum_i a_i(x) T^i, T = (1 - x^2) d/dx (as tests/test_paper59_operator_exponents.py)."""
    R = len(c) - 1
    w = flint.nmod_poly([1, 0, p - 1], p)
    zero = flint.nmod_poly([], p)
    t = [[flint.nmod_poly([1], p)]]
    for k in range(R):
        t.append([w * ((t[k][m].derivative() if m <= k else zero) + (t[k][m - 1] if m >= 1 else zero))
                  for m in range(k + 2)])
    rem = [w ** R * ck for ck in c]
    a = [None] * (R + 1)
    for i in range(R, -1, -1):
        q, r = divmod(rem[i], w ** i)
        if r.degree() >= 0:
            return None
        a[i] = q
        for m in range(i + 1):
            rem[m] = rem[m] - q * t[i][m]
    return a if all(r.degree() < 0 for r in rem) else None


def apply_T(a, seq, p):
    n = len(seq) + 2
    L = np.arange(n, dtype=np.int64)
    inv = lambda v: np.array([pow(int(z) % p, p - 2, p) for z in v], dtype=np.int64)     # noqa: E731
    xm, xp = L % p * inv(2 * L - 1) % p, (L + 1) % p * inv(2 * L + 3) % p
    tm, tp = (L - 1) * L % p * inv(2 * L - 1) % p, (L + 1) * (L + 2) % p * inv(2 * L + 3) % p

    def xop(s):
        m = len(s) - 1
        out = xp[:m] * s[1:] % p
        out[1:] = (out[1:] + xm[1:m] * s[:m - 1]) % p
        return out

    def top(s):
        m = len(s) - 1
        out = tp[:m] * s[1:] % p
        out[1:] = (out[1:] - tm[1:m] * s[:m - 1]) % p
        return out
    S = [np.array([int(v) % p for v in seq], dtype=np.int64)]
    for _ in range(len(a) - 1):
        S.append(top(S[-1]))
    total = None
    for i, ai in enumerate(a):
        co = [int(v) for v in ai.coeffs()]
        if not co:
            continue
        r = co[-1] * S[i] % p
        for cj in reversed(co[:-1]):
            r = xop(r)
            r = (r + cj * S[i][:len(r)]) % p
        m = len(r) if total is None else min(len(total), len(r))
        total = r[:m] if total is None else (total[:m] + r[:m]) % p
    return total


def annihilates(c_lists, seqs, p, min_len=400):
    a = T_form([flint.nmod_poly([int(v) for v in co], p) for co in c_lists], p)
    assert a is not None
    for sq in seqs:
        r = apply_T(a, sq, p)
        assert len(r) >= min_len, len(r)
        if np.any(r):
            return False
    return True


@pytest.mark.parametrize('obj,case,p', OCP)
def test_endpoint_piece_operators_annihilate_their_sequences(obj, case, p):
    """The five endpoint piece operators read here annihilate their embedded sequences (L <= 1000), by the T-form
    route; for (AA|BC), whose pieces are carried by this guard's own data file, record MO's embedded FULL operator
    (both weight seeds) annihilates each of the five sequences too, and a D1 piece operator does not annihilate a D2
    sequence.  Wrong answers rejected: a corrupted or mislabelled piece operator; (AA|BC) sequences that are not
    the pieces of MO's function."""
    e, keys = entry(obj, case, p), endpoint_keys(case)
    assert e['case'] == case and e['prime'] == p and set(keys.values()) <= set(e['piece_ops'])
    for nm, key in keys.items():
        assert len(e['seqs'][key]) == 1001
        assert annihilates(e['piece_ops'][key], [e['seqs'][key]], p), nm
    assert not annihilates(e['piece_ops'][keys['lone']], [e['seqs'][keys['D2+']]], p)
    assert [len(e['piece_ops'][keys[nm]]) - 1 for nm in ('lone', 'fam-', 'fam+', 'D2+', 'D2-')] == [5, 6, 6, 10, 10]
    if obj == 'AABC':
        mo = _json(str(MO_EMBED))['entries']
        for sd in SEEDS:
            m = mo[f'{case}|{p}|base|{sd}']
            assert set(keys.values()) <= set(m['keys'])
            assert annihilates(m['op'], [e['seqs'][k] for k in keys.values()], p, min_len=600), sd


# -------------------------------------------------------------------------------------- claims (2), (3)
@pytest.mark.parametrize('obj,case,p', OCP)
def test_each_endpoint_piece_has_exactly_one_half_integer_root(obj, case, p):
    """Claim (2): at x* each of the five endpoint piece operators is regular singular, all its indicial roots are
    accounted for, exactly one is a half-integer and it is simple; the other roots are non-negative integers; the
    log-free half-integer solution space is one-dimensional and starts at that root: 3/2 for the three D1 pieces;
    5/2 (f') or 7/2 ((AA|BC)) for the two D2 pieces.  Wrong answers rejected: two branches for one piece (the
    counter-example below has them); a half-integer root with a logarithm; a D2 piece at 3/2."""
    e, keys = entry(obj, case, p), endpoint_keys(case)
    want = {'lone': Fraction(3, 2), 'fam-': Fraction(3, 2), 'fam+': Fraction(3, 2)}
    want['D2+'] = want['D2-'] = Fraction(5, 2) if obj == 'FPRIME' else Fraction(7, 2)
    for nm, key in keys.items():
        C = expand_at(e['piece_ops'][key], XSTAR[case], p)
        roots = indicial_roots(C, p)
        h = halves_of(roots)
        assert h == [want[nm]], (nm, h)
        assert all(r.denominator == 1 and r >= 0 for r in roots if r not in h), (nm, roots)
        basis, exps = half_integer_space(e['piece_ops'][key], XSTAR[case], p)
        assert len(basis) == 1 and exps == h, (nm, exps)
        hroots, obs = frobenius_obstructions(e['piece_ops'][key], XSTAR[case], p)
        assert [Fraction(1, 2) + n for n in hroots] == h and obs == {}


@pytest.mark.parametrize('obj,case,p', OCP)
def test_family_halves_carry_one_branch(obj, case, p):
    """Claim (3).  The two halves of the D1 family (constants e^{-mu (D2 - D1)}, e^{-mu (D2 + D1)}) are two different,
    non-proportional sequences with ONE minimal operator, and so are the two D2 halves (e^{+zC D1}, e^{-zC D1}).
    'Identical branch series' is therefore not a separate measurement: it is this operator equality together with
    claim (2) (one half-integer root, one-dimensional space).  The lone D1 piece has another operator and another
    series with the same leading exponent.  For f' the embedded D2 block operator (the generic sum of the two
    halves) is the piece operator again: the family adds no dimension.  Wrong answers rejected: halves with
    different operators (a branch that depends on the constant); halves that are one sequence up to a factor (then
    'a family' would be one piece); a lone piece with the family's series."""
    ser = branch_series(obj, case, p)
    e, keys = entry(obj, case, p), endpoint_keys(case)
    for a, b in (('fam-', 'fam+'), ('D2+', 'D2-')):
        assert e['piece_ops'][keys[a]] == e['piece_ops'][keys[b]], (a, b)
        sa, sb = e['seqs'][keys[a]], e['seqs'][keys[b]]
        i = next(j for j, v in enumerate(sa) if v)
        ratio = sb[i] * pow(sa[i], p - 2, p) % p
        assert any((x * ratio - y) % p for x, y in zip(sa, sb)), (a, b)          # not proportional
        assert ser[a] == ser[b] and len(ser[a][1]) == NTERMS - int(ser[a][0] - Fraction(1, 2))
    assert e['piece_ops'][keys['lone']] != e['piece_ops'][keys['fam-']]
    assert ser['lone'][0] == ser['fam-'][0] and ser['lone'][1] != ser['fam-'][1]
    assert ser['lone'][1][0] == ser['fam-'][1][0] == ser['D2+'][1][0] == 1
    assert any(ser['fam-'][1][1:]) and any(ser['D2+'][1][1:])
    if obj == 'FPRIME':
        for sd in SEEDS:
            assert e['ops'][f'D2|{sd}'] == e['piece_ops'][keys['D2+']]


# ------------------------------------------------------------------------------------------ claims (4), (5)
def reduce_prediction(q, p):
    """A predicted rational as a residue mod p (the only place side A meets a prime).  The denominator must be a
    unit; kept apart from fp(), which the measurement uses for x* and 1/2."""
    q = Fraction(q)
    num, den = q.numerator % p, q.denominator % p
    assert den, 'predicted denominator vanishes mod p'
    return num * pow(den, p - 2, p) % p


def _compare(pr, me, p):
    """Predicted rationals against measured residues.  Returns {name: bool}."""
    out = {'lead_D1': pr['lead']['D1'] == me['lead']['D1'], 'lead_D2': pr['lead']['D2'] == me['lead']['D2'],
           'd3': reduce_prediction(pr['d3'], p) == me['d3'],
           'rho_D1': reduce_prediction(pr['rho_D1'], p) == me['rho_D1'],
           'rho_D2': reduce_prediction(pr['rho_D2'], p) == me['rho_D2'],
           'resid': reduce_prediction(pr['resid'], p) == me['resid']}
    for n in pr['d']:
        out[f'd{n}'] = reduce_prediction(pr['d'][n], p) == me['d'][n]
    return out


@pytest.mark.parametrize('obj,case,p', OCP)
def test_pinch_germ_matches_the_piece_operators(obj, case, p):
    """Claim (4), the sentence's check.  Predicted (side A) = measured (side B) mod p:
      f':       d5, d7 (D1: Y half minus NA piece), rho_D2 (D2 piece);                3 independent coefficients
      (AA|BC):  d5 = 0, d7, d9 (D1: screened half minus point piece), rho_D2;        4 (one of them a zero)
    and, not in the memo, d3 = [t^{5/2}] of the lone point-charge piece itself (the sentence says 'each endpoint
    piece'); rho_D1 = d_{n+2}/d_n follows from the d's.  The leading exponents agree with the germ's leading power.
    Wrong answers rejected: see the controls below; a prediction reduced with the wrong prime; a perturbed
    operator; a series normalised at the wrong coefficient."""
    pr, me = predictions(obj, case), measured(obj, case, p)
    got = _compare(pr, me, p)
    assert all(got.values()), got
    assert me['nd'] == pr['nd']
    assert sorted(k for k in got if k.startswith('d') and k != 'd3') == (
        ['d5', 'd7'] if obj == 'FPRIME' else ['d5', 'd7', 'd9'])
    if obj == 'AABC':
        assert me['d'][5] == 0 and me['d'][7] != 0
    else:
        assert me['d'][5] != 0


@pytest.mark.parametrize('obj,case,p', OCP)
def test_between_endpoint_residual_is_not_zero(obj, case, p):
    """Claim (5), the '+1' condition: rho_D2 - rho_D1 != 0 as an exact rational and mod p as measured, and the two
    agree.  Wrong answers rejected: a zero residual (then the collision would not resolve at the next order); a
    residual that is nonzero only on one side."""
    pr, me = predictions(obj, case), measured(obj, case, p)
    assert pr['resid'] != 0 and reduce_prediction(pr['resid'], p) != 0
    assert me['resid'] != 0 and me['rho_D1'] != me['rho_D2']
    assert reduce_prediction(pr['resid'], p) == me['resid']
    assert reduce_prediction(pr['rho_D1'], p) == me['rho_D1']
    assert reduce_prediction(pr['rho_D2'], p) == me['rho_D2']


@pytest.mark.parametrize('obj,case,p', OCP)
def test_negative_controls_fail_to_match(obj, case, p):
    """Controls that must NOT match, on every case and prime:
      leading-only germ (no |a|^2 correction):  rho_D2 fails, the last d fails, d3 fails (the first d still holds);
      D2 endpoint part used for the D1 pieces:  every nonzero d fails, d3 fails;
      another geometry's prediction:            rho_D2 and the nonzero d's fail;
      the residual forced to 0:                 fails.
    Wrong answer rejected: a comparison that passes whatever is predicted."""
    me = measured(obj, case, p)
    last = max(predictions(obj, case)['d'])
    first = predictions(obj, case)['nd']
    lo = _compare(predictions(obj, case, correction=False), me, p)
    assert not lo['rho_D2'] and not lo[f'd{last}'] and not lo['d3'] and not lo['resid']
    assert lo[f'd{first}'] and lo['lead_D2']
    sw = _compare(predictions(obj, case, d1_endpoint='D2'), me, p)
    assert not sw[f'd{first}'] and not sw[f'd{last}'] and not sw['d3'] and not sw['rho_D1']
    assert sw['rho_D2']                                              # the D2 side is untouched by this control
    other = next(c for c in OBJECTS[obj]['cases'] if c != case)
    ot = _compare(predictions(obj, other), me, p)
    assert not ot['rho_D2'] and not ot[f'd{first}'] and not ot[f'd{last}'] and not ot['d3']
    assert not _compare(predictions(obj, case, force_zero_residual=True), me, p)['resid']


def test_reviewer_residues_at_fprime_A1():
    """An independent reviewer's from-scratch residues at f', A1, p = 1073741441 (typed): d5 = 121155036,
    rho_D1 = 1025160734, rho_D2 = 87906117; the reviewer's failing controls: leading-only rho_D2 = 952265355, the
    D2 value in D1's place d5 = 17011549.  Wrong answer rejected: a measurement or reduction that differs from a
    second implementation's."""
    p = 1073741441
    me, pr = measured('FPRIME', 'A1', p), predictions('FPRIME', 'A1')
    assert (me['d'][5], me['rho_D1'], me['rho_D2']) == (121155036, 1025160734, 87906117)
    assert tuple(reduce_prediction(v, p) for v in (pr['d'][5], pr['rho_D1'], pr['rho_D2'])) == (
        121155036, 1025160734, 87906117)
    assert reduce_prediction(predictions('FPRIME', 'A1', correction=False)['rho_D2'], p) == 952265355
    assert reduce_prediction(predictions('FPRIME', 'A1', d1_endpoint='D2')['d'][5], p) == 17011549


# ---------------------------------------------------------------------------------------------- claim (6)
def _span(rows, p):
    return rref([list(r) for r in rows], p)[0]


@pytest.mark.parametrize('case,p', [(c, p) for c in OBJECTS['FPRIME']['cases'] for p in PRIMES])
def test_fprime_blocks_and_full_operator_are_log_free_and_span_the_pieces(case, p):
    """Claim (6), f'.  D1 block: half-integer roots {3/2, 5/2}; D2 block {5/2}; full operator {3/2, 5/2, 7/2}; at
    every resonance the Frobenius right-hand side vanishes (no logarithm), so the log-free space has the full
    dimension; and it equals the span of the pieces' series (D1 block: lone + family; D2 block: the D2 series;
    full: all three), both weight seeds, to NTERMS terms.  Wrong answers rejected: a logarithm at 5/2 or 7/2; a
    block branch that is not a combination of the pieces' germs; a family contributing two dimensions."""
    e, xs = entry('FPRIME', case, p), XSTAR[case]
    ser = branch_series('FPRIME', case, p)
    pad = lambda nm: [0] * int(ser[nm][0] - Fraction(1, 2)) + list(ser[nm][1])   # noqa: E731
    want = {'D1': ([Fraction(3, 2), Fraction(5, 2)], ['lone', 'fam-', 'fam+']),
            'D2': ([Fraction(5, 2)], ['D2+', 'D2-']),
            'full': ([Fraction(3, 2), Fraction(5, 2), Fraction(7, 2)], ['lone', 'fam-', 'fam+', 'D2+', 'D2-'])}
    for name, (exps, pieces) in want.items():
        for sd in SEEDS:
            op = e['ops'][f'{name}|{sd}']
            assert halves_of(indicial_roots(expand_at(op, xs, p), p)) == exps
            hroots, obs = frobenius_obstructions(op, xs, p)
            assert [Fraction(1, 2) + n for n in hroots] == exps
            assert len(obs) == len(exps) * (len(exps) - 1) // 2 and all(v == 0 for v in obs.values()), obs
            basis, got = half_integer_space(op, xs, p)
            assert got == exps and len(basis) == len(exps)
            assert _span([pad(nm) for nm in pieces], p) == basis, (name, sd)


@pytest.mark.parametrize('case,p', [(c, p) for c in OBJECTS['AABC']['cases'] for p in PRIMES])
def test_aabc_full_operator_is_log_free_and_spans_the_pieces(case, p):
    """Claim (6), (AA|BC): record MO's embedded full operator (both seeds) has half-integer roots {3/2, 7/2, 9/2}
    at x*, no logarithm, and its log-free half-integer space is the span of this guard's five piece series.  This
    also ties the piece data carried here to the committed MO operators.  Wrong answers rejected: pieces of another
    function; a 5/2; a logarithm."""
    mo, xs = _json(str(MO_EMBED))['entries'], XSTAR[case]
    ser = branch_series('AABC', case, p)
    pad = lambda nm: [0] * int(ser[nm][0] - Fraction(1, 2)) + list(ser[nm][1])   # noqa: E731
    exps = [Fraction(3, 2), Fraction(7, 2), Fraction(9, 2)]
    for sd in SEEDS:
        op = mo[f'{case}|{p}|base|{sd}']['op']
        assert halves_of(indicial_roots(expand_at(op, xs, p), p)) == exps
        hroots, obs = frobenius_obstructions(op, xs, p)
        assert [Fraction(1, 2) + n for n in hroots] == exps and len(obs) == 3 and not any(obs.values()), obs
        basis, got = half_integer_space(op, xs, p)
        assert got == exps
        assert _span([pad(nm) for nm in ser], p) == basis, sd


def _operator_mod_p(y1, y2, x, p):
    """Order-2 operator with solutions y1, y2 (Wronskian construction, checked), denominators cleared, mod p."""
    W = sp.simplify(y1 * sp.diff(y2, x) - y2 * sp.diff(y1, x))
    p1 = sp.simplify(-sp.diff(W, x) / W)
    p0 = sp.simplify(-(sp.diff(y1, x, 2) + p1 * sp.diff(y1, x)) / y1)
    for y in (y1, y2):
        assert sp.simplify(sp.diff(y, x, 2) + p1 * sp.diff(y, x) + p0 * y) == 0
    p0, p1 = sp.cancel(sp.together(p0)), sp.cancel(sp.together(p1))
    assert p0.is_rational_function(x) and p1.is_rational_function(x)
    den = sp.lcm(sp.denom(p0), sp.denom(p1))
    co = [sp.Poly(sp.cancel(v * den), x) for v in (p0, p1, sp.Integer(1))]
    return [[fp(Fraction(int(sp.Rational(q).p), int(sp.Rational(q).q)), p) for q in reversed(P.all_coeffs())]
            for P in co]


@pytest.mark.parametrize('p', PRIMES)
def test_one_branch_per_piece_is_not_automatic_and_the_log_detector_sees_a_log(p):
    """Claim (7) and the instrument's controls, exact (sympy), read by this file's code at x* = 5/3:
      (i)   f = (x* - x)^{1/2} log(x - x1): the operator of f and its continuation (x* - x)^{1/2} has TWO
            half-integer exponents {1/2, 3/2}, no logarithm at x*, and a two-dimensional half-integer space;
      (ii)  contrast: solutions (x* - x)^{3/2} (1 + x^2) and 1: one half-integer exponent {3/2};
      (iii) solutions t^{3/2} and t^{1/2} + t^{3/2} log t: roots {1/2, 3/2} but the resonant right-hand side is
            NOT zero and the log-free space is one-dimensional: the Frobenius check can fail.
    Wrong answers rejected: a reader for which every piece has one branch; a log check that cannot see a log."""
    x = sp.Symbol('x')
    xs = sp.Rational(5, 3)
    s = sp.sqrt(xs - x)
    h12, h32 = Fraction(1, 2), Fraction(3, 2)
    op = _operator_mod_p(s * sp.log(x + 2), s, x, p)
    assert halves_of(indicial_roots(expand_at(op, Fraction(5, 3), p), p)) == [h12, h32]
    basis, exps = half_integer_space(op, Fraction(5, 3), p)
    assert exps == [h12, h32] and len(basis) == 2
    assert frobenius_obstructions(op, Fraction(5, 3), p)[1] == {(h12, h32): 0}
    op = _operator_mod_p(s ** 3 * (1 + x ** 2), sp.Integer(1), x, p)
    assert halves_of(indicial_roots(expand_at(op, Fraction(5, 3), p), p)) == [h32]
    assert half_integer_space(op, Fraction(5, 3), p)[1] == [h32]
    t = xs - x
    op = _operator_mod_p(t ** sp.Rational(3, 2), sp.sqrt(t) + t ** sp.Rational(3, 2) * sp.log(t), x, p)
    assert halves_of(indicial_roots(expand_at(op, Fraction(5, 3), p), p)) == [h12, h32]
    obs = frobenius_obstructions(op, Fraction(5, 3), p)[1]
    assert list(obs) == [(h12, h32)] and obs[(h12, h32)] != 0
    basis, exps = half_integer_space(op, Fraction(5, 3), p)
    assert exps == [h32] and len(basis) == 1


# ------------------------------------------------------------------------------ provenance of the (AA|BC) pieces
def _load_engine():
    sys.path.insert(0, str(SUP))
    try:
        mods = []
        for name in ('s1_controls', 'aabc_exact'):
            spec = importlib.util.spec_from_file_location(f'p59pg_{name}', SUP / f'{name}.py')
            m = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(m)
            mods.append(m)
    finally:
        sys.path.remove(str(SUP))
    return mods


def _aabc_prefix(case, p, lmax):
    S1, X = _load_engine()
    D1, D2, zB, zC, zA = params(case)
    K = X.E.Field(p)
    f = X.combine(X.pieces(case, lmax, K, params_override=(D1, D2, zB, zC, zA)), zA, K)
    return {repr(k): [int(q.get(k, flint.nmod(0, p))) for q in f] for k in set(k for q in f for k in q)}


@pytest.mark.parametrize('case', OBJECTS['AABC']['cases'])
def test_aabc_piece_sequences_match_the_engine_prefix(case):
    """The five (AA|BC) endpoint sequences carried by this guard are the frozen engine's (aabc_exact.pieces +
    combine) at this geometry: terms L <= 12, one prime in the fast tier.  Wrong answers rejected: sequences of
    another geometry, of f', or of NA alone."""
    p = PRIMES[0]
    e, keys = entry('AABC', case, p), endpoint_keys(case)
    seqs = _aabc_prefix(case, p, 12)
    for nm, key in keys.items():
        assert seqs[key] == e['seqs'][key][:13], nm
        assert any(seqs[key][2:])
    fpr = entry('FPRIME', 'A1', p)
    assert e['seqs'][keys['fam-']][:13] != fpr['seqs'][endpoint_keys('A1')['fam-']][:13]


# ============================================================================================================ slow tier
@pytest.mark.slow
@pytest.mark.parametrize('case,p', [(c, p) for c in OBJECTS['AABC']['cases'] for p in PRIMES])
def test_aabc_piece_operators_recomputed_by_gcrd(case, p):
    """The (AA|BC) piece operators carried here are what the frozen guesser returns on the carried sequences (about
    1 min per case and prime), and the engine prefix L <= 40 on both primes."""
    S1, _ = _load_engine()
    e, keys = entry('AABC', case, p), endpoint_keys(case)
    for nm, key in keys.items():
        o = S1.piece_operator(e['seqs'][key], p)
        assert [[int(a) for a in ck.coeffs()] for ck in o['c']] == e['piece_ops'][key], nm
    seqs = _aabc_prefix(case, p, 40)
    for nm, key in keys.items():
        assert seqs[key] == e['seqs'][key][:41], nm


# Fire-test record, 2026-11-18.  Every plant restored; no paper, record, debug/ file or existing tests/data file was
# written.  "match" = test_pinch_germ_matches_the_piece_operators.
# (a) in this file, python debug/qa/fire_test.py tests/test_paper59_pinch_germ.py --plant "O=>N" -k <sel>:
#   F01 Hessian term 2 (p+q+5) => 2 (p+q+3)                   match, formula_is_derived, predicted_rationals   FIRED x3
#   F02 segment-shift term sign (zB/2 => -zB/2)               match, formula_is_derived, predicted_rationals   FIRED x3
#   F03 u <-> 1-u in the endpoint split                       match, orientation, predicted_rationals          FIRED x3
#   F04 f' D2 pieces see the point charge instead of V'       match, potentials_seen                           FIRED x2
#   F04b (AA|BC) D2 pieces see V'                             match                                            FIRED
#   F04c (AA|BC) D1 family sees the Yukawa potential          match                                            FIRED
#   F05 swapped primes (entry of the other prime)             match, between_endpoint                          FIRED x2
#   F06 reduction of the prediction: numerator only           match, between_endpoint, reviewer                FIRED x3
#   F07 residual forced to 0                                  between_endpoint, match, predicted_rationals     FIRED x3
#   F08 leading-only germ as the prediction                   match, negative_controls                         FIRED x2
#   F09 D2 endpoint part used for the D1 pieces               match, negative_controls                         FIRED x2
#   F10 weights with p and q exchanged                        match                                            FIRED
#   F11 cusp coefficient without 1/q!                         match                                            FIRED
#   F13 measured rho_D2 one index off                         match                                            FIRED
#   F14 Frobenius right-hand side recorded as 0               log_detector                                     FIRED
#   F15 typed rho_D2, last digit                              predicted_rationals                              FIRED
#   F16 measured difference => sum                            match                                            FIRED
#   F17 pinch constant Gamma(-(p+q+1)/2)                      match, formula_is_derived                        FIRED x2
#   F18 du = 2 r dr / (D2^2 - D1^2) with the r dropped        match, orientation                               FIRED x2
#   F22 x*(A1) 125/117 => 125/118                             match, exactly_one_half                          FIRED x2
#   F23 |a|^2 power of the correction term                    match                                            FIRED
#   F24 d3 read one coefficient later                         match                                            FIRED
#   F25 half-integer space truncated to one vector            not_automatic, log_free                          FIRED x2
#   F12 cusp sign (-zC)^q => zC^q                             match          DID NOT FIRE
#   F19 antiderivative with the opposite sign                 match          DID NOT FIRE;   orientation FIRED
#   F21 D1 endpoint part with the opposite sign               match          DID NOT FIRE;   orientation FIRED
#        F12, F19, F21: not wrong answers for the match.  Each multiplies every coefficient of one endpoint's germ by
#        -1, and only normalised series (ratios within one endpoint) are compared.  The guard does not test the
#        overall amplitude or sign of a germ, and the sentence does not claim it.  F19 and F21 are caught by the
#        segment-integral check of the split; F12 by nothing here.
# (b) temp copies of the embedded files in a scratch directory, module paths pointed at the copy, test functions
#     called directly:
#   T01 f' A1: lone piece operator, leading coefficient entry +1       all five (annihilation, one root, family,
#        match, residual) FIRED
#   T02 f' A1: lone piece operator, c_0 constant term +1               annihilation FIRED; the rest DID NOT FIRE
#        (expected: c_0 enters the local system 5 or more orders above the window; annihilation is its guard)
#   T03 f' A1: family half operator, one c_3 entry +1                  annihilation, family, match, residual FIRED;
#        one root DID NOT FIRE (expected: the entry is not in the indicial polynomial)
#   T04 f' A1: D2 piece operator, one c_10 entry +1                    all five FIRED
#   T05 f' A1: the two primes' entries swapped                         all five FIRED
#   T06 f' A1 and A4 entries swapped                                   all five FIRED
#   T07 (AA|BC) A1: D2 piece operator, one c_10 entry +1               all five FIRED, span against MO's operator FIRED
#   T08 (AA|BC) A1, second prime: lone piece operator, leading +1      all five FIRED
#   T09 (AA|BC) A1 and A2 entries swapped                              all five FIRED
#   T10 (AA|BC) A1: primes swapped                                     all five FIRED
#   T11 (AA|BC) A1: family sequence term L = 5 +1                      annihilation FIRED, engine prefix FIRED; match
#        DID NOT FIRE (expected: the match reads operators, not sequences)
#   T12 (AA|BC) A1: family sequence term L = 500 +1                    annihilation FIRED; prefix DID NOT FIRE (L <= 12)
#   T13 f' A1: full operator, leading coefficient entry +1             block span FIRED; match DID NOT FIRE (expected:
#        the match reads piece operators only)
# Not planted: an error in the shared engine or guesser (the frozen copies; hash-pinned by
# tests/test_paper59_operator_exponents.py).
