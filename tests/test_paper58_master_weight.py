"""Paper 58 -- closed form of the two-centre master integral f0 at the exchange point, and why it is not weight one.

Guard for this result (derived and cold-reviewed 2026-11-19; this file is the guard only, written as a separate
pass, CLAUDE.md Sec. 9):

    For 0 < a < b, Pachucki's master integral f0(R; a, b) = int_R^oo F^(-1)(s) ds/s (arXiv:0905.1635 labels 58-60;
    f0' = -F^(-1)/R) has the explicit closed form
        f0 = c sum_i k_i T(lam_i, mu_i; R) + c LC [E1(2bR) - E1(2aR)],   c = 1/(2(a^2 - b^2)),  LC = ln((a+b)/(b-a)),
        T(lam, mu; R) = K - g ln R - (1/2) ln^2 R - g Ein(lam R) - ln R Ein(lam R) + E2(lam R) - Phi(lam R, mu R),
    (definitions in tests/_master_weight.py); the six-term Li2 constant C(a, b) = f0(0) vanishes; and (Theorem A)
    f0 is not a C(R)-linear combination of e^{wR}, e^{wR} ln R, e^{wR} Ei(nu R) for any finite rate sets, because
    its ln R coefficient at R = 0 is -Lam, Lam = 2c [Cinh(2bR) - Cinh(2aR)] = -R^2 + O(R^4), whose derivative has
    e^{2bR} component -c/R, and rho' + 2b rho = -c/R has no rational solution.  At a -> b the weight-two part
    W2 = -c sum k_i Phi tends to a weight-one limit at rate O(eps).

WHAT THIS FILE TESTS.
  A1  the six (k, lam, mu) and the LC terms equal Calc.FM1 of tests/_pachucki_extension.py exactly (over Q at four
      pairs and over Q(t)).
  A2  f0' = -F^(-1)/R: in the helper calculus exactly, and for the closed form numerically.
  B1  each T(lam, mu; R) against direct quadrature, both signs of mu (mu > 0 uses Re Li2(p/mu), p/mu > 1).
  B2  the double series Phi against its integral form int_0^1 e^{xu} Ein(yu) du/u.
  B3  the assembled closed form against the helper's quadrature f0_quadrature at (1, 2) (a derivation pair) and
      (3/5, 7/4) (used by neither the derivation nor the review), >= 45 digits at R <= 1.6.
  C1  C(a, b) = 0 to 1e-45 at five pairs, two of them irrational ratios (1, sqrt 2), (1, sqrt 3).
  D1  -Lam = R^2 + (a^2 + b^2) R^4/6 + O(R^6) from its sympy definition; numeric Lam equals the sympy one.
  D2  the ln R coefficient of the TRUE f0 (helper quadrature, not the closed form) is -Lam: f0 + Lam ln R has no
      R^2 ln R term (a second-difference test at R ~ 1e-6 that reads 8 ln 2 x the residual log coefficient).
  D3  symbolically, d/dR(-Lam) = -L/R with L the ln-coefficient of F^(-1) built from the six terms, and its
      e^{2bR} component is exactly -c/R.
  D4  the right side found in D3 admits no rational rho (ansatz P(R)/(R^m0 (R-1)^m1), deg P <= 5) with a positive
      control that the ansatz solver does find solutions when they exist; and the pole-order bookkeeping
      (a pole of order m of rho gives order m + 1 of rho' + 2b rho; a polynomial gives a polynomial).
  E1  a -> b at a = 1, R = 0.7: (W2 - limit)/eps bounded and ~ -0.287 at eps = 1e-3, 1e-4; d_x Phi formula.

WHAT THIS FILE DOES NOT TEST (and must not be read as asserting).
  * Theorem A's deduction (log-splitting, independence of exponentials over C(R)) is an argument; its inputs are
    pinned (D1-D4), the deduction is not a computation.
  * Theorem B (no dilogarithm / product / algebraic form) is an UNPROVED sketch: nothing here asserts it.
  * Novelty of the closed form relative to the classical literature; anything about values f0(R0) at fixed R0.
  * "C = 0" is checked numerically at five pairs, not proved as a Li2 identity here.
  * The a -> b statement for the exchange integral K itself (that rests on the symbolic Sugiura-limit guard in
    tests/test_paper58_pachucki_extension.py); E1 is about W2 only.

FAST TIER / SLOW TIER.  Fast tier (default, ~30-40 s): A1, A2, B1, B2, B3 at two points, C1, D1-D4 at (1, 2),
E1 at R = 0.7.  Slow tier adds: B3 at R = 0.35/1.6 for both pairs and at R = 3.1 (tolerance relaxed to 1e-38 there:
the entire power series cancel at size e^{2(a+b)R}, losing ~8 digits), D2 at (3/5, 7/4), D4 at (3/5, 7/4),
E1 at R = 1.9.  Fast-tier gaps: the large-R regime of B3 and the second pair of D2/D4.

Each test's docstring names the wrong answer it rejects; fire-test plants are listed at the end of the file.
Nothing here reads or imports debug/.
"""
from __future__ import annotations

from fractions import Fraction

import mpmath as mp
import pytest
import sympy as sp

import _master_weight as MW
import _pachucki_extension as X

F = Fraction
DPS = 50

# C1 pairs: rational and irrational ratios (an exponent-snapping derivation would pass the rational ones only).
C_PAIRS = (("1", "2"), ("3/4", "2"), ("2", "3"), ("1", "sqrt2"), ("1", "sqrt3"))


def _num(s: str):
    if s.startswith("sqrt"):
        return mp.sqrt(int(s[4:]))
    return MW.mpq(F(s))


# ================================================================ A: transcription and the sign of f0'
@pytest.mark.parametrize("a, b", [(F(1), F(2)), (F(1), F(3, 2)), (F(1), F(3)), (F(3, 5), F(7, 4))])
def test_terms_equal_helper_FM1(a, b):
    """Rejects: a mistranscribed (k, lam, mu) table, or the LC terms with the wrong sign, against Calc.FM1."""
    c = 1 / (2 * (a * a - b * b))
    mine = {}
    for k, lam, mu in MW.terms(a, b):
        mine.setdefault(("Ei", lam, mu), {0: 0})[0] += k * c
    mine[("E", -2 * b, 1)] = {0: c}
    mine[("E", -2 * a, 1)] = {0: -c}
    assert X.clean(mine) == X.Calc(X.Ctx("Q", a, b)).FM1


def test_terms_equal_helper_FM1_over_Qt():
    """Rejects: a table right at the tested pairs only (identity over Q(t), a = 1, b = t)."""
    ctx = X.Ctx("Qt")
    a, b = ctx.a, ctx.b
    c = 1 / (2 * (a * a - b * b))
    mine = {}
    for k, lam, mu in MW.terms(a, b):
        mine.setdefault(("Ei", lam, mu), {0: ctx.zero})[0] += ctx.c(k) * c
    mine[("E", -2 * b, 1)] = {0: c}
    mine[("E", -2 * a, 1)] = {0: -c}
    assert X.clean(mine) == X.Calc(ctx).FM1


def test_f0_derivative_is_minus_FM1_over_R():
    """Rejects: f0' = +F^(-1)/R (in the helper calculus, or a closed form built on that sign)."""
    a, b = F(1), F(2)
    calc = X.Calc(X.Ctx("Q", a, b))
    assert calc.diff({("f",): {0: F(1)}}) == X.fscale(calc.FM1, F(-1), -1)
    with mp.workdps(40):
        for (pa, pb, Rv) in ((F(1), F(2), "0.8"), (F(3, 5), F(7, 4), "1.3")):
            fm1 = X.Calc(X.Ctx("Q", pa, pb)).FM1
            for R in (mp.mpf(Rv),):
                d = mp.diff(lambda x: MW.f0_closed(pa, pb, x), R)
                target = -X.evaluate(pa, pb, fm1, R) / R
                assert abs(d - target) < mp.mpf(10) ** -25 * max(1, abs(target)), (pa, pb, R, d, target)


# ================================================================ B: the closed form
# (lam, mu): mu < 0 with lam < 0 and lam > 0; mu > 0 (principal-value Ei, p/mu > 1).
T_CASES = ((-5, -1), (3, -9), (-7, 3), (-4, 1))


def test_T_closed_vs_quadrature():
    """Rejects: +pi^2/12 (for -5 pi^2/12) when mu > 0; complex Li2(p/mu) instead of its real part; + sign on the
    ln^2 R term; Phi without 1/(j+k); a wrong finite part K for mu < 0."""
    with mp.workdps(DPS):
        for lam, mu in T_CASES:
            for R in (mp.mpf("0.3"), mp.mpf("1.8")):
                q = MW.T_quad(lam, mu, R)
                cf = MW.T_closed(lam, mu, R)
                assert abs(cf - q) < mp.mpf(10) ** -45 * max(1, abs(q)), (lam, mu, R, cf, q)


def test_Phi_series_vs_integral():
    """Rejects: the double series without its 1/(j+k) (it would be a different function)."""
    with mp.workdps(DPS):
        for x, y in ((mp.mpf("1.4"), mp.mpf("-2.8")), (mp.mpf("-3.2"), mp.mpf("1.1")), (mp.mpf("4"), mp.mpf("-6"))):
            ser = MW.Phi(x, y)
            itg = mp.quad(lambda u: mp.exp(x * u) * MW.Ein(y * u) / u, [0, 0.5, 1])
            assert abs(ser - itg) < mp.mpf(10) ** -45 * max(1, abs(itg)), (x, y)


def _check_f0(a, b, R, tol):
    with mp.workdps(DPS):
        R = mp.mpf(R)
        h = X.f0_quadrature(a, b, R)
        cf = MW.f0_closed(a, b, R)
        assert abs(cf - h) < tol * abs(h), (a, b, R, cf, h)


@pytest.mark.parametrize("a, b, R", [(F(1), F(2), "0.35"), (F(3, 5), F(7, 4), "1.6")])
def test_f0_closed_vs_helper_quadrature(a, b, R):
    """Rejects: any error in the assembled closed form (LC sign, a term, the overall sign implied by f0' = +F/R),
    against an independent quadrature of label 59; (3/5, 7/4) was used by neither the derivation nor the review."""
    _check_f0(a, b, R, mp.mpf(10) ** -45)


@pytest.mark.slow
@pytest.mark.parametrize("a, b, R, tol", [(F(1), F(2), "1.6", -45), (F(3, 5), F(7, 4), "0.35", -45),
                                          (F(1), F(2), "3.1", -38), (F(3, 5), F(7, 4), "3.1", -38)])
def test_f0_closed_vs_helper_quadrature_slow(a, b, R, tol):
    """As above; at R = 3.1 the tolerance is 1e-38 (entire series cancel at size e^{2(a+b)R}, ~8 digits lost)."""
    _check_f0(a, b, R, mp.mpf(10) ** tol)


def test_C_vanishes():
    """Rejects: an LC sign flip (C(3/4, 2) = -0.225); a constant correct only at rational exponent ratios
    (checked at (1, sqrt 2), (1, sqrt 3)); a wrong K in any of the six terms."""
    with mp.workdps(DPS):
        for sa, sb in C_PAIRS:
            C = MW.C_ab(_num(sa), _num(sb))
            assert abs(C) < mp.mpf(10) ** -45, (sa, sb, C)


# ================================================================ D: the log coefficient and Theorem A's inputs
Rs, As, Bs = sp.symbols("R a b", positive=True)


def test_log_coefficient_series():
    """Rejects: Lam = 0 or sign-flipped: -Lam must be R^2 + (a^2 + b^2) R^4/6 + O(R^6); the mpmath Lam must be the
    same function as the sympy one."""
    lam = MW.Lam_sym(As, Bs, Rs)
    ser = sp.series(-lam, Rs, 0, 6).removeO()
    assert sp.simplify(ser - (Rs ** 2 + (As ** 2 + Bs ** 2) * Rs ** 4 / 6)) == 0, ser
    with mp.workdps(30):
        for (a, b, R) in ((F(1), F(2), "0.7"), (F(3, 5), F(7, 4), "1.9")):
            v_sym = sp.N(MW.Lam_sym(sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator),
                                    sp.Rational(R)), 30)
            assert abs(MW.Lam(a, b, mp.mpf(R)) - mp.mpf(str(v_sym))) < mp.mpf(10) ** -25


def _second_difference(a, b):
    """Q/r^2 for g = f0 + Lam ln R at r, 2r, 4r: O(r) if the ln R coefficient of f0 is -Lam, 8 ln 2 x the residual
    R^2 ln R coefficient otherwise.  f0 from the helper's quadrature (independent of the closed form)."""
    with mp.workdps(30):
        r = mp.mpf("1e-6")
        g = {}
        for m in (1, 2, 4):
            R = m * r
            g[m] = X.f0_quadrature(a, b, R) + MW.Lam(a, b, R) * mp.log(R)
        Q = (g[4] - 2 * g[2]) - 4 * (g[2] - 2 * g[1])
        return Q / r ** 2


def test_log_coefficient_of_f0_numeric():
    """Rejects: f0 without ln R at R = 0 (Lam = 0: reads 8 ln 2 = 5.5) or with Lam sign-flipped (reads 11.1)."""
    v = _second_difference(F(1), F(2))
    assert abs(v) < mp.mpf("1e-3"), v


@pytest.mark.slow
def test_log_coefficient_of_f0_numeric_second_pair():
    v = _second_difference(F(3, 5), F(7, 4))
    assert abs(v) < mp.mpf("1e-3"), v


def _e2b_component():
    """(d/dR of the log coefficient -Lam, its e^{2bR} coefficient), symbolic in a, b, R."""
    dA = sp.diff(-MW.Lam_sym(As, Bs, Rs), Rs)
    dA = sp.expand(dA.rewrite(sp.exp))
    comp = sp.simplify(dA.coeff(sp.exp(2 * Bs * Rs)))
    return dA, comp


def test_log_coefficient_derivative_symbolic():
    """Rejects: a log coefficient whose derivative is not -L/R (L = ln-coefficient of F^(-1) from the six terms),
    and any e^{2bR} component other than exactly -c/R."""
    c = 1 / (2 * (As ** 2 - Bs ** 2))
    L = c * sum(k * sp.exp(lam * Rs) for k, lam, _ in MW.terms(As, Bs))
    dA, comp = _e2b_component()
    assert sp.simplify(dA - sp.expand(-L / Rs)) == 0
    assert sp.simplify(comp - (-c / Rs)) == 0, comp


@pytest.mark.parametrize("a, b", [(F(1), F(2)), pytest.param(F(3, 5), F(7, 4), marks=pytest.mark.slow)])
def test_no_rational_solution(a, b):
    """Rejects: a rational rho with rho' + 2b rho = (e^{2bR} component from D3).  The right side is taken from the
    symbolic extraction, so a log coefficient with a different component (0, or no 1/R) would make it solvable."""
    _, comp = _e2b_component()
    av, bv = sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator)
    rhs = sp.simplify(comp.subs({As: av, Bs: bv}))
    for m0 in range(4):
        for m1 in range(3):
            sols = MW.rational_ansatz_solutions(rhs, bv, Rs, 5, m0, m1)
            assert sols == [], (m0, m1, sols)
    # positive controls: the ansatz solver finds solutions when they exist
    cv = 1 / (2 * (av ** 2 - bv ** 2))
    assert MW.rational_ansatz_solutions(-cv, bv, Rs, 2, 0, 0) != []
    target = sp.diff(3 / Rs ** 2 + Rs, Rs) + 2 * bv * (3 / Rs ** 2 + Rs)
    assert MW.rational_ansatz_solutions(target, bv, Rs, 3, 2, 0) != []


def test_pole_order_bookkeeping():
    """Rejects: the bookkeeping step 'a pole of order m of rho gives order m + 1 of rho' + 2b rho' (and a polynomial
    gives a polynomial), on which 'no rational solution' rests; the right side -c/R has a pole of order 1 at 0."""
    bv = sp.Rational(2)
    cv = 1 / (2 * (1 - bv ** 2))
    assert MW.pole_order(-cv / Rs, Rs, 0) == 1
    cands = [(1 / Rs, 0, 1), (5 / Rs ** 3 + Rs, 0, 3), ((Rs ** 2 + 1) / Rs ** 2, 0, 2),
             (3 / (Rs - 1) ** 2, 1, 2), ((2 * Rs + 7) / (Rs * (Rs - 1)), 1, 1), (Rs ** 3 - 2, 0, 0)]
    for rho, p, m in cands:
        assert MW.pole_order(rho, Rs, p) == m, (rho, p)
        lhs = sp.diff(rho, Rs) + 2 * bv * rho
        assert MW.pole_order(lhs, Rs, p) == (m + 1 if m >= 1 else 0), (rho, p)
        if m == 0 and sp.denom(sp.together(rho)) == 1:
            assert sp.denom(sp.together(lhs)) == 1


# ================================================================ E: a -> b
def _w2_ratio(R, eps):
    with mp.workdps(DPS):
        a = mp.mpf(1)
        e = mp.mpf(eps)
        return (MW.W2(a, a + e, mp.mpf(R)) - MW.W2_limit(a, mp.mpf(R))) / e


def test_dxPhi_formula():
    """Rejects: a wrong weight-one limit function (d_x Phi by parts vs its integral)."""
    with mp.workdps(DPS):
        for x, y in ((mp.mpf("1.4"), mp.mpf("-2.8")), (mp.mpf("3.8"), mp.mpf("-7.6"))):
            itg = mp.quad(lambda u: mp.exp(x * u) * MW.Ein(y * u), [0, 1])
            assert abs(MW.dxPhi(x, y) - itg) < mp.mpf(10) ** -45 * max(1, abs(itg))


def test_W2_limit_rate_R07():
    """Rejects: W2 diverging like 1/eps as a -> b (Phi sum not cancelling c ~ 1/eps), or a wrong limit:
    (W2 - limit)/eps must be bounded, ~ -0.287 at a = 1, R = 0.7, and settle as eps shrinks."""
    r3, r4 = _w2_ratio("0.7", "1e-3"), _w2_ratio("0.7", "1e-4")
    assert abs(r3 - r4) < mp.mpf("5e-4"), (r3, r4)
    assert abs(r4 - mp.mpf("-0.28704")) < mp.mpf("1e-4"), r4


@pytest.mark.slow
def test_W2_limit_rate_R19():
    r3, r4 = _w2_ratio("1.9", "1e-3"), _w2_ratio("1.9", "1e-4")
    assert abs(r3 - r4) < mp.mpf("5e-2"), (r3, r4)
    assert abs(r4 - mp.mpf("-30.483")) < mp.mpf("1e-2"), r4


# ================================================================ FIRE TESTS (2026-11-19)
# Each plant was applied in memory by a scratch pytest plugin that replaces one function of the imported module
# after collection (tests/_master_weight.py = (M), tests/_pachucki_extension.py = (H), this file's C_PAIRS = (T));
# no repo file was edited.  The named selection was run, then the file was run clean (18 passed, 7 skipped).
# "item" = the numbered list "Wrong answers a regression guard should reject" of the cold review.
#
#   item  plant                                                    -k selector                         result
#   1     Calc.diff: f0' = +F^(-1)/R (H)                           f0_derivative                       FIRED
#   1     closed form negated (f0 built on f0' = +F/R) (M)         f0_derivative / f0_closed_vs_helper FIRED x2
#   1     helper quadrature negated (H)                            f0_closed_vs_helper                 FIRED
#   1     LC sign flipped (M)                                      C_vanishes (C(1,2) = -0.254)        FIRED
#           same plant                                             f0_closed_vs_helper                 FIRED
#   2     +pi^2/12 for mu > 0 (M)                                  T_closed_vs_quadrature / C_vanishes FIRED x2
#   2     complex Li2(p/mu), no real part (M)                      T_closed_vs_quadrature              FIRED
#   2     + sign on the (1/2) ln^2 R term (M)                      T_closed_vs_quadrature              FIRED
#   3     Phi without 1/(j+k) (= e^x Ein(y)) (M)                   T_closed_vs_quadrature / Phi_series FIRED x2
#   4     LC from b snapped to a rational (denominator <= 1e6) (M) C_vanishes (fails at (1, sqrt 2) only, 5.5e-13)  FIRED
#   4     same plant AND C_PAIRS reduced to rational pairs (M, T) C_vanishes                          NOT FIRED (expected: why the irrational pairs are there)
#   5     Lam = 0 (numeric and sympy) (M)                          log_coefficient_series / _of_f0_numeric (5.545 = 8 ln 2) / _derivative / no_rational_solution   FIRED x4
#   5     Lam sign flipped (both) (M)                              log_coefficient_series / _of_f0_numeric (11.09) / _derivative   FIRED x3
#   5     numeric Lam flipped only (M)                             log_coefficient_of_f0_numeric / log_coefficient_series   FIRED x2
#   6     W2 term 1 x (1 + 1e-6): Phi sum not O(eps) (M)           W2_limit_rate_R07 / _R19 --slow     FIRED x2
#   6     W2 limit doubled (M)                                     W2_limit_rate_R07                   FIRED
#   7     log coefficient with cosh for Cinh: e^{2bR} part -2bc, no 1/R; rational rho = -c solves (M)
#                                                                  log_coefficient_derivative / no_rational_solution   FIRED x2
#   7     ansatz solver returns [] always (M)                      no_rational_solution (positive control)   FIRED
#   7     pole bookkeeping: order m instead of m + 1 (M)           pole_order_bookkeeping              FIRED
#   --    k_4 = +2 -> -2 in the six-term table (M)                 terms_equal_helper (x5) / C_vanishes  FIRED x2
#
# NOT PINNED, and why.
#   * Theorem A's deduction itself (log split, C(R)-independence of exponentials): an argument, not a computation.
#   * Theorem B, novelty, values at fixed R: not claims of this guard (see header).
#   * C = 0 as an identity: five pairs numerically, two irrational.
#   * The large-R closed form (R = 3.1) and the second pair of D2/D4 are slow tier only.
