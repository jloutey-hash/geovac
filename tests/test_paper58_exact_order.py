"""Paper 58 -- the unequal-exponent exchange integrals satisfy no ODE in R of order < 9.

Guard for this result (derived and cold-reviewed 2026-11-18; this file is the guard only, written as a separate
pass, CLAUDE.md Sec. 9):

    For every real 0 < a < b, the two-centre exchange integral K(R; a, b) -- 1s, and the four sigma = 0 integrals
    2psig = (2p0_A 1s_B | 2p0_A 1s_B), 3dsig = (3d0_A 1s_B | 3d0_A 1s_B), ppsig = (2p0_A 2p0_B | 2p0_A 2p0_B),
    pppi0 = the complex pair 2p+1_A, 2p+1_B -- satisfies no homogeneous linear ODE in R with polynomial
    coefficients of order < 9.

The argument has two legs, given the closed form K = sum_i c_i(t, R) g_i on nine functions (t = b/a; guarded by
tests/test_paper58_pachucki_extension.py):

    (L) the nine functions are linearly independent over C(R).  A HAND PROOF: each g_i = A_i ln R + B_i with A_i,
        B_i entire; the log coefficients A_i carry the five distinct rates 0, +-2a, +-2b, and A_1 (of f0) has
        A_1' = -L/R with L the log coefficient of F^(-1), which no combination of exponentials with rational
        coefficients can have.
    (D) the 9 x 9 matrix M(t, R) of coordinates of K, K', ..., K^(8) on the nine functions has determinant not
        identically zero in R at every real t > 1.  EXACT COMPUTATION: a bad t0 would be a common root of the
        numerators of det M(t, R_j) for every rational R_j, since the entries have poles in t only at 0, +-1.

WHAT THIS FILE TESTS.
  * The inputs of (L), not (L) itself: the derivative table of the nine functions used in the proof equals the
    helper calculus identically over Q(t) (A1); F^(-1) = c(g7 - g6 - g9 + g8) exactly (A2); f0, defined by
    f0' = -F^(-1)/R, against an independent integral (A3); the log coefficients A_2..A_9 and A_1 numerically, at
    rational pairs and at b/a = sqrt 2 (A4-A6).
  * (D) completely, for all five cases: entry denominators (B1), the common factor of the determinant numerators
    at three values of R (C1), and that the polynomials it is computed from are the determinants of independent
    fresh solves (C2).
  * The stored closed form against fresh exact solves at ratios other than t = 2, including t = 3 where two Ei
    arguments coincide (B2).

WHAT THIS FILE DOES NOT TEST.
  * The hand proof of (L) is not machine-checked.  Steps 2-4 (splitting a relation by the logarithm, the
    dominant-rate argument, "rho' + 2b rho = -c/R has no rational solution") are an argument, not a test.
  * The closed form itself, identically in t: that is tests/test_paper58_pachucki_extension.py.  If a stored
    coefficient were wrong the determinant here would be that of the wrong vector (B2 compares at two or three
    ratios, not identically).
  * The bound order <= 9 (same file).  "Exactly 9" = this file + that bound.
  * Nothing is asserted for a > b in the mixed cases 2psig, 3dsig (the l > 0 orbital carries the SMALLER
    exponent here), for m != 0 integrals, or at a = b (order 4): the proof does not claim them.

FAST TIER / SLOW TIER.  Fast (default): everything for 1s and 2psig; legs (L) inputs; leg (D) and the denominator
audit for all five cases.  The fast tier does NOT cover, for 3dsig, ppsig, pppi0: the stored coefficients against
a fresh solve at t = 3 and t = 7/3, and the determinant polynomials against fresh-solve determinants (both need
the 576-index solve, about half a minute per ratio; slow tier).  In the fast tier those three cases rest on the
stored coefficients as guarded by tests/test_paper58_pachucki_extension.py.

Why leg (D) is tested through a gcd over several R and never through "det != 0 at one (t, R)": for the four l > 0
cases det M(t, R0) at fixed R0 has real roots in t > 1 (they move with R), so a single-R test fails at an unlucky
t and cannot see a ratio that is bad at every R.

Each test's docstring names the wrong answer it rejects; fire-test plants are listed at the end of the file.
Nothing here reads or imports debug/.
"""
from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

import _exact_order as EO
import _pachucki_extension as X

DATA = json.loads((Path(__file__).resolve().parent / "data" / "paper58_pachucki_ext_2026_11_18.json").read_text())

F = Fraction
CASES = X.CASES
FAST_CASES = ("1s", "2psig")                    # live solve over Q: 0.1 s, 3 s; the other three ~30 s per ratio

# Values of R for leg (D).  None was used by the derivation (1, 2, 7/5, 1/3, 5) or by its review (3/2, 11/7, 9).
R_VALUES = (F(2, 3), F(13, 4), F(6))
# Ratios for the fresh solves: t = 3 (Ei(-(a+b)R) and the first Ei of psi_a coincide), and one more.
FRESH_T = {"1s": (F(3), F(17, 10)), "2psig": (F(3), F(17, 10)),
           "3dsig": (F(3), F(7, 3)), "ppsig": (F(3), F(7, 3)), "pppi0": (F(3), F(7, 3))}
# Exponent pairs for the numerics of Step 1: the coincidence ratio, a pair with a != 1, an irrational ratio.
PAIRS = {"t=3": lambda: (mp.mpf(1), mp.mpf(3)),
         "a=4/5,b=3/2": lambda: (mp.mpf(4) / 5, mp.mpf(3) / 2),
         "t=sqrt2": lambda: (mp.mpf(1), mp.sqrt(2))}
DPS = 40
EPS = "1e-4"                                    # node spacing of the fifth difference
TOL_SMOOTH = mp.mpf("1e-14")                    # measured <= 8e-17 for g2..g9 at the three pairs
TOL_SMOOTH_F0 = mp.mpf("1e-8")                  # measured <= 7e-11 (in units of e^2); wrong A_1: >= 0.059


def _split(cases=CASES, fast=FAST_CASES):
    return [c if c in fast else pytest.param(c, marks=pytest.mark.slow) for c in cases]


# ---------------------------------------------------------------- shared exact objects
def stored_coefficients(case: str):
    """{name: {power of R: rational function of t, string}} -- the frozen closed form (a = 1, b = t)."""
    return DATA["sym"][case]


@lru_cache(maxsize=None)
def leg_d(case: str):
    """(ctx, rows, per-row multiplicities of t, t-1, t+1, other denominator cofactors) for the stored K."""
    ctx = X.Ctx("Qt")
    rows = EO.derivative_rows(ctx, EO.stored_K(ctx, stored_coefficients(case)), 8)
    per_row, other = EO.denominator_audit(rows)
    return ctx, rows, per_row, other


@lru_cache(maxsize=None)
def numerators(case: str):
    """[(R0, N, scale)]: N(t) = scale * prod_i D_i(t) * det M(t, R0), integer polynomial; D_i clears row i."""
    ctx, rows, per_row, other = leg_d(case)
    cleared = EO.clear_rows(ctx, rows, per_row)
    return [(R0,) + EO.determinant_polynomial(cleared, R0) for R0 in R_VALUES]


def common_root_report(polys):
    """(number of real roots in (1, oo) of the gcd of the given integer polynomials, gcd degree, stripped
    multiplicities).  The powers of t, t - 1, t + 1 are divided out first: they have no root in (1, oo), and
    what is left of each polynomial is coprime to them, so the gcd's roots above 1 are those of the gcd of the
    stripped polynomials."""
    stripped = [EO.strip_linear(list(N)) for N in polys]
    g = EO.common_factor([r for _, r in stripped])
    return EO.real_roots_above_one(g), g.degree(), [m for m, _ in stripped]


@lru_cache(maxsize=None)
def _fresh_solve(Dkey: tuple, tv: Fraction):
    ctx = X.Ctx("Q", 1, tv)
    calc, fS, Fs, P2, nfree = X.solve_all(ctx, list(Dkey))
    return ctx, calc, fS


@lru_cache(maxsize=None)
def fresh(case: str, tv: Fraction):
    """(ctx, calc, K) from a fresh exact solve of the derivative equations over Q at (a, b) = (1, tv)."""
    ctx, calc, fS = _fresh_solve(tuple(X.case_downset(case)), tv)
    return ctx, calc, X.assemble_K(ctx, case, fS)


def _nz(d):
    return {p: v for p, v in d.items() if v != 0}


def _mpq(v: Fraction):
    return mp.mpf(v.numerator) / v.denominator


# ================================================================ A. the inputs of the independence lemma (L)
@pytest.mark.parametrize("field", ["Qt", "Q(4/5,3/2)"])
def test_derivative_table_of_the_proof_equals_the_helper_calculus_identically(field):
    """A1.  The table "g_i' on the nine functions" of the proof, typed by hand in tests/_exact_order.py, equals
    in_nine(Calc.diff(g_i)) of the reviewed helper: as rational functions of t (a = 1, b = t), and over Q at a
    pair with a != 1.  This is closure of the nine-function span under d/dR and the rule f0' = -F^(-1)/R.

    Rejects: one sign flipped in a row (e.g. -2 g3/R -> +2 g3/R in psi_a'); f0' = +F^(-1)/R; g6' with e^{-2aR}/R
    for e^{-2bR}/R.  Rejects them IDENTICALLY in t: entries carry the denominator t^2 - 1, so agreement at three
    sample ratios would not determine them.
    """
    ctx = X.Ctx("Qt") if field == "Qt" else X.Ctx("Q", F(4, 5), F(3, 2))
    mine, theirs = EO.hand_table(ctx), EO.helper_table(ctx)
    assert len(mine) == len(theirs) == 9
    for i in range(9):
        assert theirs[i] is not None, f"g_{i + 1}' leaves the span of the nine"
        m = {j: _nz(d) for j, d in mine[i].items() if _nz(d)}
        h = {j: _nz(d) for j, d in theirs[i].items() if _nz(d)}
        assert m == h, (X.NAMES9[i], m, h)
    # the five rates of the log coefficients and the three plain exponentials are what Steps 3-4 separate
    assert mine[0] and set(mine[0]) == {5, 6, 7, 8}


@pytest.mark.parametrize("field", ["Qt", "Q(4/5,3/2)"])
def test_F_minus_one_is_c_times_g7_minus_g6_minus_g9_plus_g8(field):
    """A2.  The helper's F^(-1) (label 60, the function whose integral defines f0) is exactly
    c (g7 - g6 - g9 + g8), c = 1/(2(a^2 - b^2)), on the nine functions.  The proof reads the log coefficient L
    of F^(-1) off this identity.

    Rejects: the g6/g7 pair or the psi_a/psi_b pair with swapped signs (which would flip or cancel L), a wrong
    prefactor, an LC term that does not combine into psi.
    """
    ctx = X.Ctx("Qt") if field == "Qt" else X.Ctx("Q", F(4, 5), F(3, 2))
    c = ctx.one / (2 * (ctx.a * ctx.a - ctx.b * ctx.b))
    co = EO.fm1_on_nine(ctx)
    assert co is not None, "F^(-1) is not in the span of the nine functions"
    got = {j: _nz(d) for j, d in enumerate(co) if _nz(d)}
    assert got == {5: {0: -c}, 6: {0: c}, 7: {0: c}, 8: {0: -c}}, got


def test_f0_quadrature_equals_the_master_integral_at_ratio_three():
    """A3.  f0(R; a, b) = int_R^oo F^(-1)(s)/s ds (i.e. f0' = -F^(-1)/R with f0 -> 0) against Pachucki's master
    integral at the exchange point, integrated directly in prolate spheroidal coordinates (shares no code with
    the closed form), at (a, b, R) = (1, 3, 6/5): the ratio where two Ei arguments coincide, not used by the
    existing guard (4/5, 3/2, 23/10).

    Rejects: f0' = +F^(-1)/R (f0 changes sign); F^(-1) with the g6/g7 signs swapped (about 10% off).
    """
    a, b, R = F(1), F(3), F(6, 5)
    with mp.workdps(30):
        f0 = X.f0_quadrature(a, b, _mpq(R))
        direct, per = X.neumann_integral("f0", a, b, R, dps=9, guard=8)
        assert len(per) >= 5 and f0 > 0 and abs(direct / f0 - 1) < mp.mpf(10) ** -9, (mp.nstr(direct, 15),
                                                                                      mp.nstr(f0, 15))


@pytest.mark.parametrize("pair", [(F(1), F(3)), (F(4, 5), F(3, 2))])
def test_mp_nine_functions_are_the_helper_functions(pair):
    """A4.  The mpmath functions g2..g9 and F^(-1) of tests/_exact_order.py (used below at an irrational ratio,
    where the helper's exact Func cannot go) equal the helper's nine functions and its label-60 F^(-1) at two
    rational pairs, to 1e-35.

    Rejects: a psi typed with Ei(2x) - Ei(x) in place of Ei(2x) - 2 Ei(x), a missing LC, e^{2aR} for e^{2bR} in
    g6/g7 -- in either copy.
    """
    a, b = pair
    ctx = X.Ctx("Q", a, b)
    with mp.workdps(DPS):
        g, fm1 = EO.nine_mp(_mpq(a), _mpq(b))
        for x in (mp.mpf("0.37"), mp.mpf("2.9")):
            for i, (gf, _) in enumerate(X.nine_functions(ctx)[1:]):
                ref = X.evaluate(a, b, gf, x)
                assert abs(g[i](x) - ref) < mp.mpf(10) ** -35 * (1 + abs(ref)), (X.NAMES9[i + 1], x)
            ref = X.evaluate(a, b, X.Calc(ctx).FM1, x)
            assert abs(fm1(x) - ref) < mp.mpf(10) ** -35 * (1 + abs(ref))


@pytest.mark.parametrize("name", list(PAIRS))
def test_log_coefficients_of_g2_to_g9(name):
    """A5.  Step 1 of the proof for g2..g9: g_i - A_i ln R is analytic at R = 0 with
        A_2 = 1;  A_3 = A_4 = A_5 = 0;  A_6, A_7 = e^{2aR}, e^{2bR};  A_8, A_9 = -e^{-2aR}, -e^{-2bR}.
    Method: the fifth difference on the nodes e, ..., 6e (e = 1e-4) of an analytic function is O(e^5) (measured
    <= 8e-17); a leftover delta(R) ln R contributes 0.0875 delta(0) - 0.060 e delta'(0) + ..., so the A_i are
    pinned through their R^3 terms, using only real values of Ei (no series for Ei is typed in).

    Rejects: psi with a log count of 0 (bracket Ei(2x) - Ei(x): A_8 = 0, psi_a log-free, Step 4 fails) -- the
    test also computes that candidate and requires it to be visibly wrong; A_8 = +e^{-2aR}; A_6 and A_7 swapped;
    A_8 and A_9 swapped; A_2 = 0.  At a rational ratio, at a != 1 and at b/a = sqrt 2.
    """
    with mp.workdps(DPS):
        a, b = PAIRS[name]()
        e = mp.mpf(EPS)
        g, _ = EO.nine_mp(a, b)
        A, _, _ = EO.log_coefficients(a, b)

        def leftover(i, Ai):
            return abs(EO.fifth_difference(lambda R: g[i](R) - Ai(R) * mp.log(R), e))

        for i in range(8):
            assert leftover(i, A[i]) < TOL_SMOOTH, (X.NAMES9[i + 1], mp.nstr(leftover(i, A[i]), 5))
        # the detector is not blind: the wrong candidates named above are far outside the tolerance
        zero = lambda R: mp.mpf(0)                                   # noqa: E731
        assert leftover(0, zero) > mp.mpf("0.05")                    # A_2 = 0
        assert leftover(6, zero) > mp.mpf("0.05") and leftover(7, zero) > mp.mpf("0.05")     # psi log-free
        assert leftover(6, lambda R: -A[6](R)) > mp.mpf("0.1")       # A_8 = +e^{-2aR}
        assert leftover(4, A[5]) > 1000 * TOL_SMOOTH and leftover(5, A[4]) > 1000 * TOL_SMOOTH    # A_6 <-> A_7
        assert leftover(6, A[7]) > 1000 * TOL_SMOOTH and leftover(7, A[6]) > 1000 * TOL_SMOOTH    # A_8 <-> A_9


@pytest.mark.parametrize("name", list(PAIRS))
def test_log_coefficient_of_f0(name):
    """A6.  Step 1 for f0, the function that carries the proof's key step:
      (i)  A_1 = -2c [Cinh(2bR) - Cinh(2aR)] satisfies A_1' = -L/R and A_1(0) = 0, with L = c(A_7 - A_6 - A_9 +
           A_8) the log coefficient of F^(-1) built from the A's pinned in A5 (identity A2);
      (ii) f0 - A_1 ln R is analytic at 0: its fifth difference on e..6e, computed from f0' = -F^(-1)/R alone
           (five short quadratures; the integration constant of f0 drops out), is O(e^5); in units of e^2 it is
           measured <= 7e-11.  Since A_1 = R^2 + O(R^4), a wrong multiple k A_1 leaves 0.1185 (1 - k) e^2.

    Rejects: A_1 with the sign flipped (0.237), A_1 = 0 (0.119), A_1 halved (0.059), A_1 plus a constant; the
    test computes the first three and requires them to be visibly wrong.  At a rational ratio, at a != 1 and at
    b/a = sqrt 2.
    NOT pinned by (ii): the R^4 and higher terms of A_1 (they sit at e^4); those rest on (i).
    """
    with mp.workdps(DPS):
        a, b = PAIRS[name]()
        e = mp.mpf(EPS)
        _, fm1 = EO.nine_mp(a, b)
        _, A1, L = EO.log_coefficients(a, b)
        assert A1(mp.mpf(0)) == 0
        for R in (mp.mpf("0.3"), mp.mpf("1.7")):
            assert abs(mp.diff(A1, R) + L(R) / R) < mp.mpf(10) ** -25 * (1 + abs(L(R) / R)), R
        assert abs(L(mp.mpf("0.9"))) > mp.mpf("0.1"), "L vanishes: the key step has nothing to work with"
        d5 = EO.f0_fifth_difference(fm1, e)

        def leftover(cand):
            return abs(d5 - EO.fifth_difference(lambda R: cand(R) * mp.log(R), e)) / e ** 2

        assert leftover(A1) < TOL_SMOOTH_F0, mp.nstr(leftover(A1), 5)
        assert leftover(lambda R: -A1(R)) > mp.mpf("0.2")
        assert leftover(lambda R: mp.mpf(0)) > mp.mpf("0.1")
        assert leftover(lambda R: A1(R) / 2) > mp.mpf("0.05")
        assert leftover(lambda R: A1(R) + mp.mpf("1e-6")) > mp.mpf(1)


# ================================================================ B. the stored closed form
@pytest.mark.parametrize("case", CASES)
def test_entry_denominators_in_t_are_only_t_tm1_tp1(case):
    """B1.  Every Laurent coefficient of every entry of the 9 x 9 coordinate matrix of (K, ..., K^(8)), as a
    rational function of t, has denominator a constant times t^i (t - 1)^j (t + 1)^k.  The gcd argument of leg (D)
    needs this: at a pole t0 of an entry, "det M(t0, R) = 0 for all R" would not make t0 a root of the cleared
    numerators.

    Rejects: a new denominator factor (t - 3, 2t - 3, t^2 + 1, ...) in a stored coefficient or produced by the
    derivative rows.
    """
    ctx, rows, per_row, other = leg_d(case)
    assert len(rows) == 9 and all(len(rw) == 9 for rw in rows)
    assert other == [], f"{case}: denominator factors other than t, t - 1, t + 1: {sorted(set(other))[:5]}"
    assert all(m[1] > 0 and m[2] > 0 for m in per_row), per_row       # the poles at t = +-1 are really there


@pytest.mark.parametrize("case", _split())
def test_stored_coefficients_equal_fresh_solves_at_other_ratios(case):
    """B2.  The stored Q(t) coefficients of K on the nine functions, specialised at t = 3 and at one more ratio
    (17/10 for 1s, 2psig; 7/3 for the three slow cases), equal a fresh exact solve of the derivative equations
    over Q at (1, t) -- coefficient by coefficient.  (t = 2 is the existing guard's
    test_live_closed_form_at_1_2_is_the_stored_one; t = 15/8 with a != 1 its scaling-law test.)

    Rejects: a stored coefficient that is right at t = 2 only (e.g. one multiplied by 1 + (t - 2) q(t)); a
    closed form that breaks at t = 3, where Ei(-(a+b)R) and the first Ei of psi_a are the same function.
    """
    for tv in FRESH_T[case]:
        ctx, calc, K = fresh(case, tv)
        co = X.in_nine(ctx, K)
        assert co is not None, f"{case}, t = {tv}: fresh K is not in the span of the nine functions"
        stored = X.coeff_strings_at(stored_coefficients(case), tv)
        for i in range(9):
            assert _nz(co[i]) == stored[i], (case, str(tv), X.NAMES9[i])
        assert _nz(co[0]), "no f0 component"


# ================================================================ C. leg (D)
@pytest.mark.parametrize("case", CASES)
def test_determinant_numerators_have_no_common_root_above_one(case):
    """C1.  Leg (D).  N_j(t) = cleared numerator of det M(t, R_j), an integer polynomial computed exactly
    (integer t, Bareiss, Newton interpolation, checked on six extra nodes) for R_j = 2/3, 13/4, 6.  A ratio
    t0 > 1 at which det M(t0, R) vanished identically in R would be a root of every N_j (B1: no pole there), hence
    of their gcd.  The gcd has no real root in (1, oo).

    Rejects: a bad exponent ratio t0 > 1 (rank < 9 for all R at that ratio; e.g. the f0 coordinate of K carrying
    a factor t - 2, which leaves K in the d/dR-closed span of g2..g9 at t = 2); rank 8 for all t (K without its
    f0 or psi component: every N_j is the zero polynomial).
    Does NOT assert det != 0 at any single (t, R): for l > 0 that is false at some t > 1 for each fixed R.
    """
    polys = [N for _, N, _ in numerators(case)]
    assert len(polys) >= 2
    for (R0, N, _) in numerators(case):
        assert any(N), f"{case}: det M(t, R = {R0}) is identically zero in t -- rank < 9"
    nroots, gdeg, mults = common_root_report(polys)
    assert nroots == 0, f"{case}: the determinant numerators at R = {R_VALUES} share {nroots} real root(s) > 1"
    # every pair of R values separately: the conclusion does not hang on one lucky third value
    for i in range(len(polys)):
        for j in range(i + 1, len(polys)):
            assert common_root_report([polys[i], polys[j]])[0] == 0, (case, R_VALUES[i], R_VALUES[j])


def test_common_root_detector_sees_a_planted_bad_ratio_and_a_rank_drop():
    """C0.  Positive controls for C1's machinery, on 1s (not a claim of the proof; a check that C1 can fail).
    (i) K with its f0 coordinate multiplied by (t - 2): at t = 2 the function lies in the span of g2..g9, which
    is closed under d/dR, so the rank is <= 8 for all R and t = 2 must appear as a common root above 1.
    (ii) K with its f0 coordinate removed: every determinant numerator must be the zero polynomial.

    Rejects: a C1 whose gcd / root count / interpolation would report "no common root" whatever it is given.
    """
    ctx = X.Ctx("Qt")
    K = EO.stored_K(ctx, stored_coefficients("1s"))
    rest = {k: v for k, v in K.items() if k != ("f",)}
    for factor, expect_zero in ((ctx.b - 2, False), (ctx.zero, True)):
        Kp = X.fadd(rest, X.fscale({("f",): K[("f",)]}, factor)) if not expect_zero else rest
        rows = EO.derivative_rows(ctx, Kp, 8)
        per_row, other = EO.denominator_audit(rows)
        assert other == []
        cleared = EO.clear_rows(ctx, rows, per_row)
        polys = [EO.determinant_polynomial(cleared, R0)[0] for R0 in R_VALUES[:2]]
        if expect_zero:
            assert all(not any(N) for N in polys)
        else:
            assert all(any(N) for N in polys)
            nroots, gdeg, _ = common_root_report(polys)
            assert nroots == 1 and gdeg >= 1, (nroots, gdeg)
            g = EO.common_factor([EO.strip_linear(list(N))[1] for N in polys])
            assert g.eval(2) == 0


@pytest.mark.parametrize("case", _split())
def test_determinant_polynomials_equal_fresh_solve_determinants(case):
    """C2.  The polynomials N_j of C1 are what they are said to be.  At t = 3 and one more ratio the derivative
    rows are rebuilt from a FRESH exact solve over Q (no stored coefficient, no Q(t) arithmetic), the determinant
    is taken over Q by a different elimination, and N_j(t) = scale * prod_i D_i(t) * det must hold exactly for
    every R_j.  Also: at each of these ratios the determinant is nonzero for at least one of the three R (rank 9
    there; not required at each R).

    Rejects: an interpolation with too low a degree bound or wrong clearing exponents (C1 would then examine
    polynomials unrelated to the determinant); rows built from a closed form that is right at t = 2 only.
    """
    qt, rows_t, per_row, other = leg_d(case)
    for tv in FRESH_T[case]:
        ctx, calc, K = fresh(case, tv)
        rows = EO.derivative_rows(ctx, K, 8)
        clearing = F(1)
        for (e0, e1, e2) in per_row:
            clearing *= tv ** e0 * (tv - 1) ** e1 * (tv + 1) ** e2
        nonzero = 0
        for (R0, N, scale) in numerators(case):
            Mq = [[sum((F(x) * R0 ** p for p, x in cd.items()), F(0)) for cd in rw] for rw in rows]
            det = EO.fraction_det(Mq)
            value = sum((F(c) * tv ** k for k, c in enumerate(N)), F(0))
            assert value == scale * clearing * det, (case, str(tv), str(R0))
            nonzero += det != 0
        assert nonzero >= 1, f"{case}: det M({tv}, R) = 0 at all of R = {R_VALUES}"


# ================================================================ FIRE TESTS (2026-11-18)
# Each plant was applied in memory by a pytest plugin that patches the imported modules after collection (no
# committed file was edited), the named selection was run, and the run was repeated clean.  "item" = the numbered
# list "Wrong answers a regression test should reject" of the cold review (item 9 is prose scope, not a test).
# (H) = tests/_pachucki_extension.py, (E) = tests/_exact_order.py, (J) = the frozen coefficients as loaded.
#
#   item  plant                                                         -k selector                           result
#   1     f0 coordinate of K x (t - 2), 1s (J)                          no_common_root and 1s                 FIRED (1 common root > 1)
#   1     same, 2psig / 3dsig / ppsig / pppi0 (J)                       no_common_root and <case>             FIRED x4
#   1     all nine coordinates x (t - 2), 1s and 3dsig (J)              no_common_root and <case>             FIRED x2
#   1     all nine coordinates x (3t - 5), 2psig (J)                    no_common_root and 2psig              FIRED
#   2     f0 coordinate removed, 1s and pppi0 (J)                       no_common_root and <case>             FIRED x2 (det identically 0)
#   2     f0 and psi_a coordinates removed, 1s (J)                      no_common_root and 1s                 FIRED
#   2     psi_a coordinate ALONE removed, 1s (J)                        no_common_root and 1s                 NOT FIRED (correct: see note 1)
#           same plant                                                  fresh_solves and 1s                   FIRED
#   3     table: -2 g3/R -> +2 g3/R in psi_a' (E)                       derivative_table                      FIRED (both fields)
#   3     table row f0' x (1 + (t-2)(t-3)(2t-5)): right at 3 ratios (E) derivative_table                      FIRED (both fields)
#   3     helper: d Ei((a-b)R) = -e^{..}/R (H)                          derivative_table                      FIRED
#   4     psi_a bracket Ei(2x) - Ei(x), mpmath copy (E)                 log_coefficients_of_g2                FIRED (3 pairs)
#           same plant                                                  mp_nine_functions                     FIRED
#   4     psi_a bracket Ei(2x) - Ei(x), helper (H)                      mp_nine_functions / derivative_table / no_common_root and 1s   FIRED x3
#   4     typed A_8 = 0; A_8 = +e^{-2aR}; A_6 <-> A_7; A_2 = 0 (E)      log_coefficients_of_g2                FIRED x4 (3 pairs each)
#   5     A_1 sign flipped; A_1 = 0; A_1 halved; A_1 + 1e-9 (E)         log_coefficient_of_f0                 FIRED x4 (3 pairs each)
#   5     A_1 AND L both flipped (passes A_1' = -L/R) (E)               log_coefficient_of_f0                 FIRED (0.237 against 1e-8)
#   5     A_1 AND L both halved (E)                                     log_coefficient_of_f0                 FIRED (0.059 against 1e-8)
#   6     f0' = +F^(-1)/R in Calc.diff (H)                              derivative_table                      FIRED
#           same plant                                                  f0_quadrature                         NOT FIRED (expected: the quadrature does not use Calc.diff)
#   6     f0 quadrature sign flipped (H)                                f0_quadrature                         FIRED
#   6     label 60: g6/g7 pair sign swapped (H)                         F_minus_one / f0_quadrature / mp_nine_functions   FIRED x3
#   6     label 60: psi_a/psi_b pair sign swapped (H)                   F_minus_one / f0_quadrature           FIRED x2
#   6     label 60: prefactor x 2 (H)                                   F_minus_one                           FIRED
#   6     mpmath F^(-1): g6/g7 swapped (E)                              mp_nine_functions / log_coefficient_of_f0   FIRED x2
#   7     one psi_a coefficient x (1 + (t - 2)), 1s, 2psig (J)          fresh_solves and <case>               FIRED x2
#           same plant, 1s                                              fresh_solve_determinants and 1s       FIRED
#   7     one coefficient x (1 + (t-2)(10t-17)(3t-7)): wrong at t = 3   fresh_solves and 1s                   FIRED
#         only (J)                                                      fresh_solve_determinants and 2psig    FIRED
#   7     one psi_a coefficient x (1 + (t - 2)), 3dsig (J), slow        fresh_solves and 3dsig --slow         FIRED
#           same plant, fast tier                                       fresh_solves and 3dsig                NOT FIRED (skipped: see header, fast-tier gap)
#   7     wrong-at-3-only coefficient, pppi0 (J), slow                  fresh_solve_determinants and pppi0 --slow   FIRED
#   8     one coefficient / (t - 3), 1s (J)                             entry_denominators and 1s             FIRED
#           same plant                                                  no_common_root and 1s                 FIRED (clearing refuses)
#   8     one coefficient / (t^2 + 1), 2psig (J)                        entry_denominators and 2psig          FIRED
#   8     all coefficients / (2t - 3), ppsig (J)                        entry_denominators / no_common_root and ppsig   FIRED x2
#   --    machinery: root counter returns 0 always (E)                  common_root_detector                  FIRED
#           same plant                                                  no_common_root and 1s                 NOT FIRED (expected: why C0 exists)
#   --    machinery: interpolation degree bound 30 too low (E)          no_common_root and 1s                 FIRED (extra nodes)
#   --    machinery: clearing exponent of t off by one (E)              fresh_solve_determinants and 1s       FIRED
#   --    design: a single value of R                                   no_common_root and 2psig              FIRED (needs >= 2)
#
# Note 1.  Removing the psi_a coordinate alone does not lower the rank: psi_a returns through f0' = -F^(-1)/R, and
#   the remaining eight functions are not closed under d/dR.  The review's "(or one psi)" is too strong as worded;
#   leg (D) is right to pass, and the corrupted closed form is rejected by B2.  Only the f0 coordinate (g2..g9
#   are closed) or f0 together with others gives rank < 9.
#
# NOT PINNED, and why.
#   * The hand proof of the independence lemma (Steps 2-4).  Its inputs are pinned (A1-A6); the deduction is not
#     a computation.
#   * A_1 beyond its R^2 term by the smoothness test (the R^4 term sits at e^4 = 1e-16 of the signal); the full
#     function rests on A_1' = -L/R, checked at two R per pair, and that is all the proof uses.
#   * The stored coefficients identically in t (B2 compares at two ratios per case, the existing guard at two
#     more); and, in the fast tier, B2 and C2 for 3dsig, ppsig, pppi0.
#   * Real roots of the common factor are excluded on (1, oo) only; nothing is asserted on (0, 1) (the mirrored
#     ordering of the mixed cases is outside the claim).
