"""Paper 58 -- the four mixed exchange integrals at a > b: the same closed forms, real branch, order exactly 9.
Plus a regression guard for geovac.two_center_eri.radial_norm (exact normalisation, no nsimplify).

Guard for this result (derived and cold-reviewed 2026-11-19; this file is the guard only, written as a separate
pass, CLAUDE.md Sec. 9):

    For the four mixed two-centre exchange integrals
        2psig = (2p0_A 1s_B | 2p0_A 1s_B),     3dsig = (3d0_A 1s_B | 3d0_A 1s_B),
        2ppi  = (2p+1_A 1s_B | 2p-1_A 1s_B),   3ddel = (3d+2_A 1s_B | 3d-2_A 1s_B),
    with a the exponent of the l > 0 orbital and t = b/a, the stored Q(t) closed forms (frozen in
    tests/data/paper58_pachucki_ext_2026_11_18.json and tests/data/paper58_exchange_m_2026_11_18.json, derived for
    t > 1) hold UNCHANGED for t in (0, 1), on the real branch: Ei of a positive argument = principal value and
    LC = ln((a+b)/|b-a|).  This is the analytic continuation in t around t = 1 (every i pi picked up by the Ei and
    Log terms cancels inside each psi bracket).  The order is exactly 9 for every a != b: the same nine functions,
    the independence lemma is symmetric in a <-> b, and the gcd over several R of the determinant numerators has
    no real root in (0, 1) as well as in (1, oo).  Nothing is claimed at a = b.

WHAT THIS FILE TESTS.
  A. The continuation itself: the stored coefficients with the nine functions continued from t > 1 to
     t = s +- i eps (s = 1/2, 1/3), both half-planes, against the real-branch value: imaginary part O(eps),
     conjugate in the two half-planes, real part equal to O(eps^2); and the continued form equals the real form
     at t = 2.  The real-branch evaluator refuses any complex intermediate (no real part is ever taken).
  B. Fresh exact solves over Q at (a, b) = (2, 1) and (3, 1) (t = 1/2, and t = 1/3 where -(a+b) = 2(b-a)): every
     solve re-checked as an exact identity, zero free unknowns, every key decaying at that (a, b), all nine
     functions present, coordinates EQUAL to the stored Q(t) form specialised and rescaled.  The frozen audit:
     every intermediate key decays on all of (0, 1), Ei arguments vanish only at t = +-1, poles only at 0, +-1.
  C. Values: the real-branch closed form against exchange_hp (Neumann sum, never sees the closed form) at a > b,
     frozen hp numbers in the fast tier, live exchange_hp in the slow tier; near t = 1 ((1, 9/10), (1, 97/100),
     (1, 99/100)) with the working precision raised; the LC-flipped and a <-> b relabelled values rejected.
  D. Order exactly 9 on (0, 1): gcd over three R of the cleared determinant numerators has no real root in (0, 1)
     (nor in (1, oo)); at some single R the numerator does have roots in (0, 1); a planted bad ratio t = 1/2 is seen
     by the (0, 1) count and NOT by the old (1, oo) count; exact determinants nonzero at four t in (0, 1).
  E. The inputs of the independence lemma at t < 1 (Step 1, log coefficients, real branch), as in
     tests/test_paper58_exact_order.py for t > 1.
  F. radial_norm: exact (1/sqrt of an exact Rational), equal to the textbook normalisation for ~400 rational Z,
     n <= 3, l < n; at Z = 99/100 the norm squared times the exact integral is exactly 1; the nsimplify-based
     version is detected (it fails 27 of the 2370 p/q cases, plus (99/100, 1, 0) and three 97/100 cases).

WHAT IT DOES NOT TEST.
  * The hand argument of the continuation (holomorphy of the parameter integral in t on Re t > 0, identity
    theorem).  A tests its numerical consequence at two s, one R; B tests the exact consequence at two ratios.
  * Steps 2-4 of the independence lemma (an argument; Step 1's inputs only, at three ratios).
  * Pachucki's source formula at a > b as a separate statement: it is tested end to end by B and C.
  * Anything at a = b (order 4); a > b for other (l, m); other R for the continuation.
  * Whether the other nsimplify calls in geovac (two_center_eri.radial_product, gaunt_LM, qfd_core.s_radial)
    can misfire the same way; only radial_norm is guarded here.

FAST TIER / SLOW TIER.  Fast (default): A for 2ppi at s = 1/2; B for 2psig and 2ppi at both pairs (3dsig, 3ddel
fresh solves take 40-60 s per pair: slow); the frozen audit for all four; C against frozen hp numbers for 2ppi
and 3ddel (2psig, 3dsig and the near-t = 1 rows slow); D for 2ppi (2psig about 10-20 s, 3dsig about 85 s, 3ddel about
17 s: slow) and its positive control; E at one ratio; F entire.  Slow tier (--slow): the rest of A, B, C, D, E, and
live exchange_hp at (3, 2), R = 6/5 for all four cases, at (1, 9/10) and at (1, 99/100) (the radial_norm bug
point, end to end).  Measured wall times: in the trailing comment.

Each test's docstring names the wrong answer it rejects; fire-test plants are listed at the end of the file.
Nothing here reads or imports debug/.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import factorial

import mpmath as mp
import pytest
import sympy as sp

import _exact_order as EO
import _exchange_reverse as XR
import _pachucki_extension as X

F = Fraction
CASES = XR.CASES
FRESH_PAIRS = ((F(2), F(1)), (F(3), F(1)))           # t = 1/2 and t = 1/3 (the coincidence -(a+b) = 2(b-a))
R_VALUES = (F(2, 3), F(13, 4), F(6))                 # leg (D); none used by the a > b derivation (1, 2, 7/5, 1/3, 5,
                                                     # 3/2, 11/7, 9); 13/4 and 2/3 were used by its review for 2ppi
CONT_S = (F(1, 2), F(1, 3))
CONT_R = F(7, 5)
CONT_EPS = mp.mpf("1e-10")
CONT_DPS = 30

# exchange_hp numbers at a > b (Neumann sum over tau; independent of the closed form).  (case, a, b, R):
# (value, tolerance, closed-form working precision).  Sources: the derivation's hp runs (tau_max 24, dps 50;
# 3dsig tau_max 18, dps 48) and the review's (tau_max 14-20).  Tolerances: about 1000 x the measured difference,
# which is the size of the dropped tau tail (the closed form is exact).
HP = {
    ("2psig", F(2), F(1), F(7, 5)): ("0.213670604546883717941975415966665000588809023787", "1e-46", 60),
    ("2psig", F(3, 2), F(4, 5), F(23, 10)): ("0.140104392480534939053220628164355801123341223005", "1e-45", 60),
    ("2psig", F(5, 2), F(1), F(2)): ("0.0886391294226717081608053379725115541654046289494", "1e-33", 60),
    ("3dsig", F(2), F(1), F(7, 5)): ("0.1234259169318707540234076268615385750666804454", "1e-41", 60),
    ("3dsig", F(5, 2), F(1), F(7, 5)): ("0.1129540421886717799462961057246702863800523", "1e-35", 60),
    ("2ppi", F(5, 2), F(1), F(7, 5)): ("-0.046451008899994523282568029702237691842034793104", "1e-44", 60),
    ("2ppi", F(2), F(1), F(3)): ("-0.00522187854922980064591911044875213969838367502349", "1e-47", 60),
    ("3ddel", F(2), F(1), F(7, 5)): ("0.0215387891733121922834187780117461666527309192218", "1e-47", 60),
    ("3ddel", F(3, 2), F(4, 5), F(23, 10)): ("0.0105315786132040018249691296541742403757050817002", "1e-46", 60),
    # near t = 1: the coefficients carry poles up to (t - 1)^-9, about 18 digits cancel at t = 0.9
    ("2ppi", F(1), F(9, 10), F(2)): ("-0.0312241487371145378885618610562088107723777701229", "1e-46", 80),
    ("3dsig", F(1), F(9, 10), F(2)): ("0.09324359559429483146362937256228233826440141", "1e-41", 80),
    ("3ddel", F(1), F(9, 10), F(7, 5)): ("0.0112031473599387359900110160786098244635809286472", "1e-46", 90),
    ("2ppi", F(1), F(97, 100), F(3, 2)): ("-0.042030697632537461970338539988273018771636057388585590227", "1e-54",
                                          110),
}
HP_FAST = {("2ppi", F(5, 2), F(1), F(7, 5)), ("3ddel", F(2), F(1), F(7, 5))}

LM = {"2psig": ((1, 0), (1, 0)), "3dsig": ((2, 0), (2, 0)), "2ppi": ((1, 1), (1, -1)), "3ddel": ((2, 2), (2, -2))}


def _slow(*vals):
    return pytest.param(*vals, marks=pytest.mark.slow)


def _nz(d):
    return {p: v for p, v in d.items() if v != 0}


# ================================================================ A. the continuation around t = 1
@pytest.mark.parametrize("case,s", [("2ppi", F(1, 2))] + [_slow(c, s) for c in CASES for s in CONT_S
                                                          if (c, s) != ("2ppi", F(1, 2))])
def test_continuation_from_t_above_one_is_the_real_branch(case, s):
    """A1.  F(t) = sum c_i(t, R) g_i(t, R) with every g_i continued holomorphically from t > 1 (LC = Log(1+t) -
    Log(t-1); Ei of an argument positive at t > 1 as gamma + Log w - Ein(-w); of one negative at t > 1 as -E1(-w))
    at t = s + i eps and t = s - i eps, eps = 1e-10, R = 7/5, dps 30, against the real-branch K(s):
        |Im F| <= 10 eps, Im F(s + i eps) = -Im F(s - i eps) to 1e-6 relative, |Re F - K(s)| <= 10 eps^2.
    (Measured: |Im F|/eps = dK/dt, 0.008 to 0.29; |Re F - K| <= 2.5e-21.)  Also F = K_real at t = 2 + i 0.

    Rejects: a continuation with an i pi slip (for instance LC continued as Log(1+t) - Log(1-t), or a psi_b Ei
    continued from the wrong side): it leaves a +-i pi e^{-2bR} term, |Im F| ~ 1e-2, which a guard comparing real
    parts only would accept -- the imaginary part is asserted before the real part is used.
    """
    with mp.workdps(CONT_DPS):
        sv, R, e = XR.mpq(s), XR.mpq(CONT_R), CONT_EPS
        vr = XR.K1_typed(case, sv, R, "real")
        assert isinstance(vr, mp.mpf)
        up = XR.K1_typed(case, sv + 1j * e, R, "cont")
        lo = XR.K1_typed(case, sv - 1j * e, R, "cont")
        for v in (up, lo):
            assert abs(mp.im(v)) <= 10 * e, f"Im F = {mp.nstr(mp.im(v), 5)}: a branch slip"
            assert abs(mp.re(v) - vr) <= 10 * e ** 2, mp.nstr(mp.re(v) - vr, 5)
        assert abs(mp.im(up) + mp.im(lo)) <= mp.mpf("1e-6") * abs(mp.im(up)) + mp.mpf(10) ** -(CONT_DPS - 5)
        t2 = XR.mpq(F(2))
        assert abs(XR.K1_typed(case, t2 + 1j * mp.mpf(10) ** -25, R, "cont") - XR.K1_typed(case, t2, R, "real")) \
            < mp.mpf(10) ** -20


@pytest.mark.parametrize("pair", [(F(1), F(1, 2)), (F(1), F(1, 3))])
def test_typed_nine_functions_are_the_helper_functions_at_t_below_one(pair):
    """A2.  The typed real-branch nine functions and F^(-1) used by A1 and E equal the helper's own nine
    functions and its label-60 F^(-1) (X.nine_functions, X.Calc.FM1), evaluated on the real branch, at t = 1/2
    and 1/3, two R, to 1e-35.

    Rejects: a typed psi with +LC, a missing factor 2 on Ei((t-1)R), e^{2aR} and e^{2bR} swapped -- in either copy.
    """
    a, b = pair
    ctx = X.Ctx("Q", a, b)
    with mp.workdps(45):
        for x in (mp.mpf("0.37"), mp.mpf("2.9")):
            g, fm1 = XR.nine_typed(XR.mpq(b), x, "real")
            for i, (gf, _) in enumerate(X.nine_functions(ctx)[1:], start=1):
                ref = XR.evaluate_real(a, b, gf, x)
                assert abs(g[i] - ref) < mp.mpf(10) ** -35 * (1 + abs(ref)), (X.NAMES9[i], x)
            ref = XR.evaluate_real(a, b, X.Calc(ctx).FM1, x)
            assert abs(fm1 - ref) < mp.mpf(10) ** -35 * (1 + abs(ref))


def test_real_branch_evaluator_refuses_a_complex_intermediate(monkeypatch):
    """A3.  The real-branch evaluator asserts that LC and every base value are real mpf before using them.  With
    LC taken as the principal log of the negative number (a+b)/(b-a) (the t > 1 formula used blindly at a > b, as
    X.evaluate does), it raises; it never takes a real part.

    Rejects: an evaluator that silently returns Re() of a mis-branched expression.
    """
    a, b = F(2), F(1)
    K = XR.coords_to_func(a, b, XR.stored_coords("2ppi", a, b))
    with mp.workdps(20):
        assert isinstance(XR.evaluate_real(a, b, K, mp.mpf(3), mp.mpf(0)), mp.mpf)
        monkeypatch.setattr(XR, "lc_real", lambda a_, b_: mp.log(XR.mpq((a_ + b_) / (b_ - a_))))
        with pytest.raises(AssertionError, match="not real"):
            XR.evaluate_real(a, b, K, mp.mpf(3), mp.mpf(0))


# ================================================================ B. fresh exact solves at a > b
@lru_cache(maxsize=None)
def fresh(case: str, a: Fraction, b: Fraction):
    return XR.fresh_solve(case, a, b)


@pytest.mark.parametrize("case,a,b", [(c, a, b) if c in ("2psig", "2ppi") else _slow(c, a, b)
                                      for c in CASES for a, b in FRESH_PAIRS])
def test_fresh_solve_at_a_above_b_equals_the_stored_form(case, a, b):
    """B1.  A fresh exact solve over Q at (a, b) with a > b (t = 1/2, 1/3), with the reviewed solvers: every
    solve_L0 / solve_D2 output re-checked as an exact identity (no failures), zero free unknowns (per solve and in
    total), every key of every solved jet decaying at this (a, b) (rate lam, resp. lam + mu, < 0), K in the span
    of the nine functions with all nine present, and its coordinates EQUAL, coordinate by coordinate over Q, to
    the frozen Q(t) closed form specialised at t = b/a and rescaled (R^p g_i: a^(p+1) c_ip; f0: a^(p+3) c_0p).
    At (3, 1) the Ei arguments -(a+b) and 2(b-a) coincide.

    Rejects: a closed form valid only for t > 1 (any coefficient different at t < 1); the K(R; b, a) relabelling;
    a wrong rescaling; a solve with a free unknown or a growing key at t < 1; a wrong solver output.
    """
    K, co, st = fresh(case, a, b)
    assert st["solves"] > 50
    assert st["residual_failures"] == 0
    assert max(st["free"]) == 0 and st["returned_free"] == 0
    rates = [XR.key_rate(k) for k in st["keys"]]
    bad = [k for k, r in zip(st["keys"], rates) if r is not None and not r < 0]
    assert bad == [], bad
    assert co is not None, "K left the span of the nine functions"
    assert all(co), [X.NAMES9[i] for i, d in enumerate(co) if not d]
    want = XR.stored_coords(case, a, b)
    for i, (d, e) in enumerate(zip(co, want)):
        assert _nz(d) == _nz(e), (case, X.NAMES9[i])
    if (a, b) == (F(3), F(1)):
        assert -(a + b) == 2 * (b - a)
    # the relabelled form is a different vector (the coordinates are not symmetric in a <-> b)
    assert any(_nz(d) != _nz(e) for d, e in zip(co, XR.stored_coords(case, b, a)))


def test_residual_and_free_counters_can_fire(monkeypatch):
    """B0.  Positive control for B1's wrapper: a solver whose output is doubled is caught as a residual failure.

    Rejects: a B1 whose residual check is blind.
    """
    orig = X.Calc.solve_L0
    calls = []

    def first_doubled(self, src):
        sol = orig(self, src)
        calls.append(1)
        return X.fscale(sol, 2) if len(calls) == 1 else sol

    monkeypatch.setattr(X.Calc, "solve_L0", first_doubled)
    ctx = X.Ctx("Q", F(2), F(1))
    with XR.checked_solves() as st:
        try:
            X.solve_all(ctx, X.downset([(1, 1, 1, 1)]))
        except RuntimeError:            # later solves may see the corrupted jet as an inconsistent source
            pass
    assert st["solves"] > 0 and st["residual_failures"] > 0


@pytest.mark.parametrize("case", CASES)
def test_frozen_audit_decays_on_zero_one(case):
    """B2.  The frozen audit over Q(t) of every intermediate jet of the derivation: every function key has a rate
    linear in t that is < 0 on the whole of (0, 1); every Ei argument mu(t) is nonzero on (0, 1) (zeros only at
    t = +-1); every denominator factor is t, t - 1 or t + 1.  So the residual identities over Q(t) specialise at
    every t0 in (0, 1), and the specialised jets are the decaying solutions.

    Rejects: a key that grows somewhere in (0, 1) (for instance e^{2(1-t)R} without its Ei), an Ei argument
    vanishing inside (0, 1), a pole at some t0 in (0, 1).
    """
    t = sp.Symbol("t")
    au = XR.audit(case)
    dens = au["denominator_factors_max_multiplicity"]
    for grp in dens.values():
        for f in grp:
            assert sp.expand(sp.sympify(f, locals={"t": t})) in (t, t - 1, t + 1), (grp, f)
    nkeys = 0
    for row in au["keys"]:
        if row["kind"] == "f":
            continue
        lam = sp.sympify(row["lam"], locals={"t": t})
        mu = sp.sympify(row["mu"], locals={"t": t}) if row["kind"] == "Ei" else sp.Integer(0)
        rate = sp.expand(lam + mu)
        assert sp.Poly(rate, t).degree() <= 1
        r0, r1 = rate.subs(t, 0), rate.subs(t, 1)
        assert r0 <= 0 and r1 <= 0 and not (r0 == 0 and r1 == 0), (row, rate)
        if row["kind"] == "Ei":
            assert mu != 0
            assert all(not (z.is_real and 0 < z < 1) for z in sp.solve(mu, t)), row
        nkeys += 1
    assert nkeys >= 12


# ================================================================ C. values at a > b against exchange_hp
@pytest.mark.parametrize("key", [k if k in HP_FAST else _slow(k) for k in HP])
def test_closed_form_at_a_above_b_against_frozen_exchange_hp(key):
    """C1.  The stored closed form, real branch (principal-value Ei, LC = ln((a+b)/|b-a|)), at a > b against
    exchange_hp numbers (Neumann sum over tau; never sees the closed form), at the derivation's and the review's
    points, including t = 0.9 and 0.97 with the working precision raised (80-110 digits).

    Rejects: the closed form wrong at t < 1; a wrong rescaling; too little working precision near t = 1 (at dps 60
    the t = 0.97 row misses its tolerance by 7.6e-51 against 1e-54; the t = 0.9 rows still pass at dps 60).
    """
    case, a, b, R = key
    ref, tol, dps = HP[key]
    v = XR.value_real(case, a, b, R, dps)
    with mp.workdps(dps):
        assert abs(v - mp.mpf(ref)) < mp.mpf(tol), mp.nstr(v - mp.mpf(ref), 5)


@pytest.mark.parametrize("case,a,b,R", [("2ppi", F(5, 2), F(1), F(7, 5)), ("3ddel", F(2), F(1), F(7, 5))])
def test_wrong_continuations_are_far_off(case, a, b, R, monkeypatch):
    """C2.  The two wrong answers named by the review, computed and required to miss the hp number by more than
    10%: (i) LC with the wrong sign inside psi (and F^(-1)): off by a factor ~99 for 2ppi at (5/2, 1);
    (ii) the a <-> b relabelled value K(R; b, a) (the "symmetry" shortcut): relative error 0.84 for 2ppi.

    Rejects: a C1 tolerance or evaluator loose enough to accept either.
    """
    ref = mp.mpf(HP[(case, a, b, R)][0])
    with mp.workdps(20):
        swapped = XR.value_real(case, b, a, R, 20)
        assert abs(swapped / ref - 1) > mp.mpf("0.1"), mp.nstr(swapped, 10)
        good = XR.lc_real
        monkeypatch.setattr(XR, "lc_real", lambda a_, b_: -good(a_, b_))
        flipped = XR.value_real(case, a, b, R, 20)
        assert abs(flipped / ref - 1) > mp.mpf("0.1"), mp.nstr(flipped, 10)


def _hp_live(case, a, b, R, tau, dps):
    from geovac import qfd_core as Q
    (l1, m1), (l2, m2) = LM[case]
    A, B = sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator)
    with mp.workdps(dps):
        v, per = Q.exchange_hp((l1 + 1) * A, (l1 + 1, l1, m1), (1, 0, 0), B, (l2 + 1, l2, m2), (1, 0, 0),
                               sp.Rational(R.numerator, R.denominator), tau_max=tau, dps=dps)
        return +v, [+x for x in per]


@pytest.mark.slow
@pytest.mark.parametrize("case", CASES)
def test_closed_form_against_live_exchange_hp_at_a_new_point(case):
    """C3 (slow, 30-160 s each).  Live exchange_hp at (a, b) = (3, 2), R = 6/5 (t = 2/3; used by neither the
    derivation nor the review), tau_max 8, dps 40, against the real-branch closed form (dps 50).  Tolerance
    0.05 x |last tau term| + 1e-35: the tau tail ratio is ~3e-3 per step, so the dropped tail is below 0.01 x the
    last term (measured relative agreement 1e-14 to 1e-17).

    Rejects: the closed form wrong at t = 2/3; the LC flip (x 99) or the a <-> b relabelling (0.84) at a fresh
    point, against a route that never sees the closed form.
    """
    a, b, R = F(3), F(2), F(6, 5)
    hp, per = _hp_live(case, a, b, R, 8, 40)
    v = XR.value_real(case, a, b, R, 50)
    with mp.workdps(50):
        tol = mp.mpf("0.05") * abs(per[-1]) + mp.mpf("1e-35")
        assert abs(per[-1]) < mp.mpf("1e-10") * abs(hp)
        assert abs(v - hp) < tol, (mp.nstr(v, 30), mp.nstr(hp, 30), mp.nstr(per[-1], 3))


@pytest.mark.slow
@pytest.mark.parametrize("a,b,R,tau,dps,cdps,tol", [
    (F(1), F(9, 10), F(2), 10, 40, 80, "1e-33"),
    (F(1), F(99, 100), F(1), 14, 60, 110, "1e-50"),
])
def test_near_t_one_against_live_exchange_hp(a, b, R, tau, dps, cdps, tol):
    """C4 (slow, 1-2 min each).  2ppi near t = 1, live exchange_hp against the closed form at raised working
    precision.  The (1, 99/100) row is the radial_norm regression end to end: with the old nsimplify-based
    radial_norm, the 1s normalisation at Z = 99/100 is off by 1.76e-18 relative and this row fails (the review
    saw a constant 1.76e-18 relative offset there); with the fix it agrees to the hp precision.

    Rejects: radial_norm(99/100, 1, 0) not exact; the closed form at t near 1 evaluated with too little precision.
    """
    hp, per = _hp_live("2ppi", a, b, R, tau, dps)
    v = XR.value_real("2ppi", a, b, R, cdps)
    with mp.workdps(cdps):
        assert abs(v - hp) < mp.mpf(tol), mp.nstr(v - hp, 5)


# ================================================================ D. exactly 9 on (0, 1): leg (D)
@lru_cache(maxsize=None)
def numerators(case: str):
    ctx = X.Ctx("Qt")
    rows = EO.derivative_rows(ctx, EO.stored_K(ctx, XR.stored(case)), 8)
    per_row, other = EO.denominator_audit(rows)
    assert other == [], sorted(set(other))[:5]
    cleared = EO.clear_rows(ctx, rows, per_row)
    return [(R0,) + EO.determinant_polynomial(cleared, R0) for R0 in R_VALUES]


def gcd_report(polys):
    """(roots of the gcd in (0, 1), roots in (1, oo), [per-R residual roots in (0, 1)]) after dividing out t, t-1,
    t+1 (which have no root in either open interval)."""
    res = [EO.strip_linear(list(N))[1] for N in polys]
    g = EO.common_factor(res)
    single = [XR.roots_in(EO.common_factor([r]), 0, 1) for r in res]
    return XR.roots_in(g, 0, 1), XR.roots_in(g, 1, sp.oo), single


@pytest.mark.parametrize("case", ["2ppi", _slow("2psig"), _slow("3dsig"), _slow("3ddel")])
def test_determinant_numerators_have_no_common_root_in_zero_one(case):
    """D1.  Leg (D) on (0, 1): N_j(t) = cleared numerator of det M(t, R_j) (Bareiss at integer t, Newton
    interpolation, extra nodes), R_j = 2/3, 13/4, 6.  The gcd, and every pairwise gcd, has no real root in (0, 1)
    and none in (1, oo).  Every entry's denominator is t^i (t-1)^j (t+1)^k (asserted while clearing), so a ratio
    t0 in (0, 1) with det M(t0, R) = 0 for all R would be a common root.  At some single R_j the numerator HAS real
    roots in (0, 1) (asserted), so a single-R check is meaningless.

    Rejects: a common root in (0, 1) (a bad ratio a > b); a leg (D) that counts roots only in (1, oo) (the old
    criterion, see D2); a single-R or two-R check (>= 3 distinct R asserted).
    """
    assert len(R_VALUES) >= 3 and len(set(R_VALUES)) == len(R_VALUES)
    polys = [N for _, N, _ in numerators(case)]
    for (R0, N, _) in numerators(case):
        assert any(N), f"det M(t, {R0}) identically 0"
    n01, n1oo, single = gcd_report(polys)
    assert n01 == 0 and n1oo == 0, (n01, n1oo)
    for i in range(len(polys)):
        for j in range(i + 1, len(polys)):
            assert gcd_report([polys[i], polys[j]])[0] == 0, (case, R_VALUES[i], R_VALUES[j])
    assert max(single) >= 1, single


def test_zero_one_detector_sees_a_planted_bad_ratio_that_the_old_criterion_misses():
    """D2.  Positive control for D1's detector on 2ppi: every numerator N_j multiplied by (2t - 1), i.e. t = 1/2
    made a common root.  The (0, 1) count sees it (1 root); the old (1, oo) count (EO.real_roots_above_one) reports
    0, so a guard using only the old count would pass a bad ratio below 1.  (The same bad ratio planted upstream,
    f0 coordinate x (2t - 1) in the frozen data, is fire-tested against D1; rebuilding the rows for it takes about
    6 min, too long for this tier.)

    Rejects: a D1 whose (0, 1) count is blind; D1 rewritten with the old count.
    """
    polys = []
    for _, N, _ in numerators("2ppi")[:2]:
        polys.append([-c for c in N] + [0])                 # (2t - 1) N = 2 t N - N, ascending coefficients
        for k, c in enumerate(N):
            polys[-1][k + 1] += 2 * c
    n01, n1oo, _ = gcd_report(polys)
    assert n01 == 1 and n1oo == 0
    g = EO.common_factor([EO.strip_linear(list(N))[1] for N in polys])
    assert g.eval(sp.Rational(1, 2)) == 0
    assert EO.real_roots_above_one(g) == 0


def test_exact_determinants_nonzero_at_four_ratios_below_one():
    """D3.  For 2ppi, the 9 x 9 coordinate matrix of (K, ..., K^(8)) evaluated exactly at t in {1/10, 1/3, 1/2,
    9/10}, R = 13/4 and 2/3 (Fraction elimination, no interpolation): every determinant is nonzero, and equals
    N_j(t) / (scale x clearing), tying D1's polynomials to the matrix.

    Rejects: wrong clearing exponents or an interpolation degree bound too low in D1's polynomials.
    """
    ctx = X.Ctx("Qt")
    rows = EO.derivative_rows(ctx, EO.stored_K(ctx, XR.stored("2ppi")), 8)
    per_row, _ = EO.denominator_audit(rows)
    nums = {R0: (N, sc) for R0, N, sc in numerators("2ppi")}
    for tv in (F(1, 10), F(1, 3), F(1, 2), F(9, 10)):
        clearing = F(1)
        for (e0, e1, e2) in per_row:
            clearing *= tv ** e0 * (tv - 1) ** e1 * (tv + 1) ** e2
        for R0 in (F(13, 4), F(2, 3)):
            det = EO.fraction_det(EO.rows_at(ctx, rows, tv, R0))
            assert det != 0, (tv, R0)
            N, sc = nums[R0]
            assert sum((F(c) * tv ** k for k, c in enumerate(N)), F(0)) == sc * clearing * det, (tv, R0)


# ================================================================ E. the independence lemma's inputs at t < 1
PAIRS_E = {"t=1/2": lambda: (mp.mpf(1), mp.mpf(1) / 2), "t=1/3": lambda: (mp.mpf(1), mp.mpf(1) / 3),
           "t=1/sqrt5": lambda: (mp.mpf(1), 1 / mp.sqrt(5))}


@pytest.mark.parametrize("name", ["t=1/2", _slow("t=1/3"), _slow("t=1/sqrt5")])
def test_log_coefficients_at_t_below_one(name):
    """E1.  Step 1 of the lemma on the real branch at a > b: g_i - A_i ln R is analytic at R = 0 with the SAME
    A_i as for a < b (A_2 = 1, A_3..5 = 0, A_6, A_7 = e^{2aR}, e^{2bR}, A_8, A_9 = -e^{-2aR}, -e^{-2bR}): the fifth
    difference on e..6e (e = 1e-4) is <= 1e-14; and f0 - A_1 ln R is analytic (fifth difference from
    f0' = -F^(-1)/R alone, in units of e^2, <= 1e-8).  ln|x| has the same log coefficient for either sign of x, and
    the real branch is what enters.

    Rejects: a psi bracket that is log-free at a > b (A_8 = 0); A_8 = +e^{-2aR}; A_1 = 0 or sign-flipped (each
    computed here and required to be visibly wrong).
    """
    with mp.workdps(40):
        a, b = PAIRS_E[name]()
        e = mp.mpf("1e-4")
        g = [lambda R, i=i: XR.nine_typed(b, R, "real")[0][i] for i in range(1, 9)]
        fm1 = lambda R: XR.nine_typed(b, R, "real")[1]                         # noqa: E731
        A, A1, L = EO.log_coefficients(a, b)

        def leftover(i, Ai):
            return abs(EO.fifth_difference(lambda R: g[i](R) - Ai(R) * mp.log(R), e))

        for i in range(8):
            assert leftover(i, A[i]) < mp.mpf("1e-14"), (X.NAMES9[i + 1], mp.nstr(leftover(i, A[i]), 5))
        zero = lambda R: mp.mpf(0)                                              # noqa: E731
        assert leftover(0, zero) > mp.mpf("0.05")
        assert leftover(6, zero) > mp.mpf("0.05") and leftover(7, zero) > mp.mpf("0.05")
        assert leftover(6, lambda R: -A[6](R)) > mp.mpf("0.1")
        d5 = EO.f0_fifth_difference(fm1, e)

        def left_f0(cand):
            return abs(d5 - EO.fifth_difference(lambda R: cand(R) * mp.log(R), e)) / e ** 2

        assert left_f0(A1) < mp.mpf("1e-8"), mp.nstr(left_f0(A1), 5)
        assert left_f0(lambda R: -A1(R)) > mp.mpf("0.2")
        assert left_f0(zero) > mp.mpf("0.1")


# ================================================================ F. radial_norm regression (exact, no nsimplify)
def _textbook_norm2(Z: Fraction, n: int, l: int):
    """N^2 for r^l L_{n-l-1}^{2l+1}(lam r) e^{-Z r/n}, lam = 2Z/n: the textbook (2Z/n)^3 (n-l-1)!/(2n (n+l)!) times
    lam^(2l) (radial_poly omits the lam^l of rho^l)."""
    lam = sp.Rational(2 * Z.numerator, Z.denominator * n)
    return lam ** (2 * l + 3) * sp.Rational(factorial(n - l - 1), 2 * n * factorial(n + l))


def _exact_integral(Z: Fraction, n: int, l: int):
    """int_0^oo (sum_k c_k r^k)^2 e^{-2ar} r^2 dr by sympy integration (no shared code with radial_norm)."""
    from geovac import two_center_eri as TCE
    r = sp.Symbol("r", positive=True)
    c, a = TCE.radial_poly(Z, n, l)
    return sp.integrate(sum(v * r ** k for k, v in c.items()) ** 2 * r ** 2 * sp.exp(-2 * a * r), (r, 0, sp.oo))


SWEEP_Z = sorted({F(p, q) for q in range(1, 17) for p in range(1, 41)} | {F(99, 100), F(199, 200), F(97, 100)})


def test_radial_norm_is_exact_at_z_99_over_100():
    """F1.  radial_norm(99/100, 1, 0)^2 x the exact integral of the unnormalised radial function = 1 exactly
    (sympy, no floats), at the point where the old nsimplify-based version returned a wrong surd (relative error
    1.76e-18); also at (99/100, 2, 1), (19/10, 1, 0), (31/15, 2, 1), (40/11, 3, 0).  qfd_core uses the same function
    object (so the fix reaches exchange_hp).

    Rejects: radial_norm = 1/sqrt(nsimplify(s)) (the pre-2026-11-19 code); any float in the normalisation.
    """
    from geovac import qfd_core as Q
    from geovac import two_center_eri as TCE
    assert Q.radial_norm is TCE.radial_norm
    for Z, n, l in ((F(99, 100), 1, 0), (F(99, 100), 2, 1), (F(19, 10), 1, 0), (F(31, 15), 2, 1), (F(40, 11), 3, 0)):
        N = TCE.radial_norm(Z, n, l)
        I = _exact_integral(Z, n, l)
        assert I.is_Rational
        assert N ** 2 * I == 1, (Z, n, l, N)
        assert (N ** 2).is_Rational


def test_radial_norm_sweep_is_exact():
    """F2.  For 2388 (Z, n, l): Z = p/q (p <= 40, q <= 16, plus 99/100, 199/200, 97/100), n <= 3, l < n:
    radial_norm(Z, n, l)^2 is a Rational equal to the textbook normalisation exactly.  (Measured: the old
    nsimplify version fails 27 of the 2370 sweep cases -- Z = 19/14, 19/10, 19/2, 31/15, 40/11, 37/3, ... --
    with relative errors 2e-18 to 1.5e-17; plus (99/100, 1, 0) and 97/100 at (n, l) = (2, 1), (3, 1), (3, 2),
    up to 1.6e-16.)

    Rejects: the nsimplify-based radial_norm; a normalisation off by a factor (e.g. the textbook N^2 without
    lam^(2l)).
    """
    from geovac import two_center_eri as TCE
    bad = []
    for Z in SWEEP_Z:
        for n in (1, 2, 3):
            for l in range(n):
                N2 = TCE.radial_norm(Z, n, l) ** 2
                if not (N2.is_Rational and N2 == _textbook_norm2(Z, n, l)):
                    bad.append((str(Z), n, l))
    assert bad == [], bad[:10]
    assert len(SWEEP_Z) * 6 > 2000


def test_nsimplify_version_would_be_caught():
    """F3.  Positive control for F1/F2: the pre-fix code (1/sqrt(nsimplify(s))) rebuilt here fails the F1
    identity at Z = 99/100.  If a future sympy fixes nsimplify on this input the control is skipped (it would then
    prove nothing), not failed.

    Rejects: an F1/F2 that cannot tell the two versions apart.
    """
    from geovac import two_center_eri as TCE

    def old(Z, n, l):
        c, a = TCE.radial_poly(Z, n, l)
        s = sum(c1 * c2 * sp.factorial(k1 + k2 + 2) / (2 * a) ** (k1 + k2 + 3)
                for k1, c1 in c.items() for k2, c2 in c.items())
        return 1 / sp.sqrt(sp.nsimplify(s))

    s = sp.Rational(250000, 970299)
    if sp.nsimplify(s) == s:
        pytest.skip("this sympy's nsimplify returns the exact Rational here")
    Z = F(99, 100)
    N = old(Z, 1, 0)
    assert not (N ** 2 * _exact_integral(Z, 1, 0) == 1)
    assert not ((N ** 2).is_Rational and N ** 2 == _textbook_norm2(Z, 1, 0))
    assert abs(sp.N(N ** 2 / _textbook_norm2(Z, 1, 0) - 1, 30)) > 1e-19


# ================================================================ FIRE TESTS (2026-11-19)
# Written and fired as a separate step after the guard.  Each plant was applied in memory by a pytest plugin that
# patches the imported modules (or this module's constants) at configure / collection time; no committed file was
# edited.  (XR) = tests/_exchange_reverse.py, (EO) = tests/_exact_order.py, (J) = the frozen data as loaded,
# (T) = a constant of this file, (G) = geovac.  FIRED = the selection failed, for the reason named.
#
#   plant                                                         selection                         result
#   LC continued as Log(1+t) - Log(1-t) (i pi slip) (XR)          A1 2ppi                           FIRED (Im F = 43.8)
#   psi_b Ei continued from the negative side (XR)                A1 2ppi                           FIRED (Im F = -57.1)
#   same LC slip, A1 copy with the Im assertions removed          A1 2ppi                           FIRED anyway, by the
#        (scratch copy)                                           Re part: at finite eps the slip x Im c(s + i eps) leaks
#                                                                 an O(eps) real term (2.6e-8 > 10 eps^2).  A real-axis
#                                                                 Re() check (eps = 0) would NOT fire; A3 covers that.
#   typed psi_a with +LC (XR)                                     A2, A1, E1                        FIRED x4
#   typed psi_a bracket Ei(2x) - Ei(x) (XR)                       A2, E1                            FIRED x3
#   real-ness assertion disabled in evaluate_real (XR)            A3                                FIRED (no raise)
#   psi_a coefficient x 2t: right at t = 1/2 only (J)             B1 fast                           FIRED at (3, 1) x2,
#                                                                                                   passes at (2, 1)
#        same plant                                               B1 3dsig --slow                   FIRED at (3, 1)
#   rescaling a^p in place of a^(p+1) (XR)                        B1, C1 fast                       FIRED x6
#   a <-> b relabelled coordinates (XR)                           B1, C1 fast                       FIRED x6
#   residual check of the wrapper disabled (XR)                   B0                                FIRED
#   f0 coordinate x (2t - 1), 2ppi (J)                            D1 2ppi, B1 2ppi                  FIRED x3 (D1: one
#                                                                 common root in (0, 1), none above 1; ~9 min run)
#   LC sign flipped in the real-branch evaluator (XR)             C1 fast                           FIRED x2
#        same plant                                               C3 2ppi --slow (live hp)          FIRED
#   value_real(a, b) -> value_real(b, a) (XR)                     C1 fast                           FIRED x2
#   value_real symmetrised in a, b (XR)                           C2                                FIRED x2 (swap unseen)
#   working precision 60 for the near-t = 1 rows (T)              C1 near rows --slow               FIRED at t = 0.97
#                                                                 (7.6e-51 vs 1e-54); NOT FIRED at the three t = 0.9
#                                                                 rows (dps 60 suffices there; the docstring says so)
#   root count restricted to (1, oo) (the old criterion) (XR)     D1, D2                            FIRED (D2); D1 NOT
#                                                                 FIRED (correct: the true gcd has no root anywhere)
#   one coefficient of a determinant polynomial + 1 (EO)          D3                                FIRED
#   R_VALUES = one R; = two R (T)                                 D1                                FIRED x2 (design)
#   A_8 = +e^{-2aR} (EO)                                          E1                                FIRED (0.174)
#   A_1 = 0 (EO)                                                  E1                                FIRED
#   pre-fix radial_norm (1/sqrt(nsimplify(s))) in (G) two_center_eri and qfd_core
#                                                                 F1, F2                            FIRED x2 (F3 passes,
#                                                                 as it should); C4 --slow: FIRED at (1, 99/100)
#                                                                 (9.9e-20 vs 1e-50), passes at (1, 9/10)
#   pre-fix radial_norm in two_center_eri only (G)                F1                                FIRED (identity and
#                                                                 the exactness check)
#
# NOT PINNED, and why.
#   * The holomorphy / identity-theorem argument (an argument; A1 checks its consequence at s = 1/2, 1/3, R = 7/5).
#   * Fresh solves at t = 1/2 and 1/3 only; other t in (0, 1) rest on the Q(t) audit (B2) and the values (C).
#   * The other nsimplify sites in geovac (radial_product, gaunt_LM, qfd_core.s_radial): measured, the first and
#     last also return values off by 1e-18..1e-17 at some rational Z (e.g. s_radial at Z = 19/2, 17/4;
#     radial_product at Z pairs with 19/10, 99/100); not guarded here.  None is hit at the points this file uses
#     (C3, C4 pass at 1e-35 / 1e-50).
#
# MEASURED (2026-11-19, Windows, this machine).  Fast: 24 passed, 33 skipped, 42.5 s.  --slow: 57 passed, 873.6 s
# (3dsig leg (D) 205 s; live exchange_hp 26-110 s per case at (3, 2), 47 s and 104 s near t = 1).
