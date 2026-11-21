"""Paper 58 -- the m != 0 exchange integrals 2p_pi and 3d_delta are in the nine-function span, order exactly 9.

Guard for this result (derived and cold-reviewed 2026-11-18; this file is the guard only, written as a separate
pass, CLAUDE.md Sec. 9):

    For every real 0 < a < b, the two-centre exchange integrals
        2ppi  = (2p+1_A 1s_B | 2p-1_A 1s_B)    and    3ddel = (3d+2_A 1s_B | 3d-2_A 1s_B)
    lie in the Q(a, b)(R)-span of the same nine functions as the 1s exchange integral, hence satisfy an ODE in R of
    order <= 9.  The route is an equation derived here from Pachucki's label 39 (his Eq. (42) in the PDF, the u1 = t form of his Eq. (16), which the paper cites) (arXiv:0905.1635), first order in
    the r12 exponent, by the inverse Laplace transform:
        q f_{m+1}'' = c f_{m+1} + m^2 (d^4 - Sig d^2 + e) f_{m-1} + m (m-1)^2 (m-2) f_{m-3}'' + S_m(r)      (**)
    At (a, b) = (1, 2) the derived operators equal the committed ones, (9, 25) and (9, 35).  The sigma = 2 part of
    the real 2p_x 2p_x exchange integral is ppdel/2 with ppdel = 3 (b/a)^2 x 3ddel.  The order is exactly 9 for
    every 0 < a < b: the nine functions are independent over C(R) (proved and guarded in
    tests/test_paper58_exact_order.py) and the 9 x 9 coordinate determinant of K..K^(8) is not identically zero in
    R at any t = b/a > 1 (the gcd over several R of its numerators has no real root t > 1; tested here).

WHAT THIS FILE TESTS.
  A. The equation (**): label 13 transcription; S_1 = -(label 63)''; f_2 = label 61 (n = 1); f_4 against label 61
     (n = 3, coefficients 18, -9, -12 typed from the label); the residual of (**) with q, c, Sig, e typed HERE;
     the Laplacian controls nabla^2 r12 = 2/r12 (f_2 rung) and nabla^2 r12^3 = 12 r12 (f_4 rung) as exact
     identities; zero free unknowns in every solve, and a counter that can see one.
  B. 2ppi and 3ddel over Q: nine-function span, closure under d/dR, the unique operator at (1, 2) equal to the
     committed one, the cross-family operator NOT annihilating, values against independent numbers.
  C. Every 0 < a < b: the frozen Q(t) closed form (tests/data/paper58_exchange_m_2026_11_18.json) against live
     solves at t = 2, 3, 21/8 (and 15/8 via a = 4/5); denominators t, t - 1, t + 1 only (closed form, derivative
     rows, every intermediate jet); every intermediate key decays for every t > 1.
  D. Exactly 9 (leg (D)): gcd over three R of the determinant numerators, no common root t > 1.
  E. sigma = 2 part of real 2p_x 2p_x: normalisations and the azimuthal average from first principles; slow:
     against geovac.qfd_core.exchange_hp (Neumann expansion, no Pachucki equation).

WHAT IT DOES NOT TEST.
  * (**) as a theorem: it is proved in Paper 58 (given Pachucki's Eq. (16), for u2 != u3, w2 != w3; proof memo
    debug/pachucki_m_proof_memo.md, cold-reviewed 2026-11-20). This file checks its consequences (labels 61/63, the
    exact controls, the independent numbers), not the proof.
  * Operator equality away from (1, 2); a > b for these mixed cases; anything at a = b beyond the existing guards;
    other (l, m).
  * The independence lemma (tests/test_paper58_exact_order.py).
  * 3ddel over Q(t) live (the derivation needed ~2 h): the frozen record rests on the live checks at four ratios
    (slow tier) and on the frozen audit.

FAST TIER / SLOW TIER.  Fast (default, about 70 s; 37 tests): everything for 2ppi over Q and over the frozen Q(t) data;
stage-1 jets (36 multi-indices) with f_2 AND f_4 at two pairs, so the m = 3 machinery (m^2, the f_{m-3}'' term, its
12) is exercised in the fast tier through the f_4 control and label 61 (n = 3); leg (D) for both cases.
FAST-TIER GAPS: no 3ddel closed form is built live (576-index solve, about 25 s per pair), so the 3ddel values,
the 3ddel operator, OP_P_1_2 on 3ddel, the frozen 3ddel coefficients against live solves, the 3ddel residuals of
(**) on the full downsets and the determinant cross-check for 3ddel are slow tier; no Q(t) solve is run live
(2ppi about 3 min; the stage-1 controls over Q(t) about 4 min); no exchange_hp evaluation (20-90 s each).
Slow tier (--slow): 21 more tests, about 16 min in all.

Each test's docstring names the wrong answer it rejects; fire-test plants are listed at the end of the file.
Nothing here reads or imports debug/.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

import _exact_order as EO
import _exchange_m_nonzero as M
import _pachucki_exchange as PX
import _pachucki_extension as X

DATA = json.loads((Path(__file__).resolve().parent / "data" / "paper58_exchange_m_2026_11_18.json").read_text())

F = Fraction
Z0 = (0, 0, 0, 0)
CASES = M.CASES_M
PAIRS_STAGE1 = [(F(1), F(2)), (F(2, 3), F(7, 4))]
UNUSED = (F(2, 3), F(7, 4), F(19, 10))              # used by the review only, never by the derivation
N_INDICES = {"2ppi": {0: 144, 1: 16}, "3ddel": {0: 576, 1: 144, 2: 16}}
OP_DEGREE = {"2ppi": 25, "3ddel": 35}
MAX_POLE = {"2ppi": 7, "3ddel": 9}
R_VALUES = (F(2, 3), F(13, 4), F(6))                 # leg (D); none used by the m != 0 derivation or its review

# Controls, typed here: nabla^2 r12 = 2/r12 gives 16 f_{1111} = 8 [..f_2..]; nabla^2 r12^3 = 12 r12.
CTRL_A = (16, 8)
CTRL_A_PRIME = 12
LABEL61_N3 = (18, 9, 12)                            # label 61 at n = 3: 18 (u^2+w^2) f_2, -9 f_2'', -12 f_0

# Values.  (1, 2) R = 1.4 and (4/5, 3/2) R = 2.3: closed form of the derivation (dps 60), equal to exchange_hp to
# 47-48 digits there.  UNUSED point: exchange_hp of the review (tau_max 16, dps 32, tail ~3e-34), an independent route.
VAL = {
    ("2ppi", F(1), F(2), F(7, 5)): "-0.013353601891216017682930492464386025458452519549795",
    ("3ddel", F(1), F(2), F(7, 5)): "0.0013461743233320120177136736257341354157956824513622",
    ("2ppi", F(4, 5), F(3, 2), F(23, 10)): "-0.0072708870648988060628723129702682786806276558388417",
    ("3ddel", F(4, 5), F(3, 2), F(23, 10)): "0.00085209572345054007851997145804439878674346695593029",
}
VAL_UNUSED = {"2ppi": "-0.00515968878885911303770792463991", "3ddel": "0.000355893627169501579264221099214"}
TOL_MEMO = "1e-45"                                  # compared inside mp.workdps(60)
TOL_UNUSED = "1e-29"                     # 30 printed digits; measured agreement 4e-33 and 6e-34


def _slow(*vals):
    return pytest.param(*vals, marks=pytest.mark.slow)


def _mpq(v: Fraction):
    return mp.mpf(v.numerator) / v.denominator


# ---------------------------------------------------------------- shared exact objects
@lru_cache(maxsize=None)
def stage1(kind: str, a=None, b=None):
    """Jets on the downset of (2, 2, 1, 1): f_0 (engine), f_2 (m = 1) and f_4 (m = 3) by (**)."""
    ctx = X.Ctx("Q", a, b) if kind == "Q" else X.Ctx("Qt")
    if kind == "Qt":
        M.install_gcd_fallback()
    D = X.downset([(2, 2, 1, 1)])
    calc0, f0, Fs, P2, nfree0 = X.solve_all(ctx, D)
    calc = M.MCalc(ctx)
    RJ = M.RJets(ctx, D, calc)
    f2, fr2 = M.solve_fn(calc, RJ, D, 1, f0)
    f4, fr4 = M.solve_fn(calc, RJ, D, 3, f2, f0)
    return ctx, D, calc, RJ, {0: f0, 1: f2, 2: f4}, nfree0, fr2 + fr4


@lru_cache(maxsize=None)
def solved(case: str, a: Fraction, b: Fraction):
    """(ctx, calc, jets {k: f_{2k}}, free unknowns of f_0, free unknowns of each solve_D2, K) over Q."""
    ctx = X.Ctx("Q", a, b)
    calc, fl, nf0, free = M.solve_case(ctx, case)
    return ctx, calc, fl, nf0, free, M.assemble_K(ctx, case, fl)


@lru_cache(maxsize=None)
def f0_at(a: Fraction, b: Fraction, x: Fraction, dps: int = 60):
    with mp.workdps(dps):
        return +X.f0_quadrature(a, b, _mpq(x))


def value(case, a, b, R, dps=60):
    ctx, calc, fl, nf0, free, K = solved(case, a, b)
    with mp.workdps(dps):
        return +X.evaluate(a, b, K, _mpq(R), f0_at(a, b, R, dps))


def own_scalar_jets(ctx, D):
    """Derivative values at the exchange point of q, c, Sig, e of (**), typed HERE (not the helper's)."""
    u2, u3, w2, w3 = sp.symbols("u2 u3 w2 w3")
    ex = {"q": (u2**2 - u3**2) * (w2**2 - w3**2),
          "c": (u2**2 - u3**2 + w2**2 - w3**2) * (u2**2 * w2**2 - u3**2 * w3**2),
          "Sig": u2**2 + u3**2 + w2**2 + w3**2,
          "e": (u2**2 - w3**2) * (u3**2 - w2**2)}
    sa, sb = sp.Rational(ctx.a.numerator, ctx.a.denominator), sp.Rational(ctx.b.numerator, ctx.b.denominator)
    sub = {u3: sa, u2: sb, w2: sa, w3: sb}
    out = {}
    for nm, e in ex.items():
        d = {}
        for al in D:
            v = sp.diff(e, u3, al[0], u2, al[1], w2, al[2], w3, al[3]).subs(sub)
            if v != 0:
                d[al] = F(int(v.p), int(v.q))
        out[nm] = d
    return out


def residual_of_star(ctx, calc, RJ, D, m, fn, flow, flow3, sc):
    """The multi-indices of D at which (**) differentiated by Leibniz is NOT satisfied (written here):
    sum_U binom [q_U f''_T - c_U f_T - m^2 (delta_U0 f''''_low - Sig_U f''_low + e_U f_low) - delta_U0 m(m-1)^2(m-2)
    f''_low3] - S_m = 0,  T = alpha - U."""
    S = RJ.source(m)
    bad = []
    for al in D:
        res = X.fscale(S[al], -1)
        for U in itertools.product(*[range(x + 1) for x in al]):
            T = tuple(s - u for s, u in zip(al, U))
            bn = M.binom(al, U)
            if U in sc["q"]:
                res = X.fadd(res, X.fscale(calc.D2(fn[T]), bn * sc["q"][U]))
            if U in sc["c"]:
                res = X.fadd(res, X.fscale(fn[T], -bn * sc["c"][U]))
            if U == Z0:
                res = X.fadd(res, X.fscale(calc.D2(calc.D2(flow[T])), -m * m))
                if m >= 3:
                    res = X.fadd(res, X.fscale(calc.D2(flow3[T]), -m * (m - 1) ** 2 * (m - 2)))
            if U in sc["Sig"]:
                res = X.fadd(res, X.fscale(calc.D2(flow[T]), m * m * bn * sc["Sig"][U]))
            if U in sc["e"]:
                res = X.fadd(res, X.fscale(flow[T], -m * m * bn * sc["e"][U]))
        if res != {}:
            bad.append(al)
    return bad


def laplacian_combination(ctx, jets):
    """(a^2+b^2) j_{1111} + 2a j_{0111} + 2b j_{1011} + ab (j_{2011} + j_{0211}) - ab R^2 j_{0011}: the jets of
    nabla_1^2 psi = [a^2 + b^2 - 2a/r_A - 2b/r_B + ab (r_A/r_B + r_B/r_A - R^2/(r_A r_B))] psi (typed here)."""
    a, b = ctx.a, ctx.b
    tot = {}
    for al, co, sh in (((1, 1, 1, 1), a * a + b * b, 0), ((0, 1, 1, 1), 2 * a, 0), ((1, 0, 1, 1), 2 * b, 0),
                       ((2, 0, 1, 1), a * b, 0), ((0, 2, 1, 1), a * b, 0), ((0, 0, 1, 1), -a * b, 2)):
        tot = X.fadd(tot, X.fscale(jets[al], co, sh))
    return tot


def label63_W1(ctx):
    """W_1^{(-2)} of label 63 at u = a, w = b (typed from the label)."""
    u, w, C_ = ctx.a, ctx.b, ctx.c
    return {("E", -2 * u, 0): {-1: -(-1 / (w - u) - 1 / (u + w))},
            ("E", -2 * w, 0): {-1: -(-1 / (u - w) - 1 / (u + w))},
            ("E", -(u + w), 0): {0: -C_(2), -1: -2 / (u + w)},
            ("Ei", ctx.zero, -(u + w)): {0: C_(2)}}


def committed_operator(case: str):
    import test_paper58_exchange_d_pp as DPP
    import test_paper58_exchange_l1 as L1
    s = {"2ppi": L1.OP_P_1_2, "3ddel": DPP.OP_3DDEL_1_2}[case]
    return X.primitive(X.parse_operator(s))


def apply_operator(calc, K, op):
    ders = calc.derivatives(K, len(op) - 1)
    tot = {}
    for i, c in enumerate(op):
        for j, x in enumerate(c):
            if x:
                tot = X.fadd(tot, X.fscale(ders[i], F(int(x)), j))
    return tot


def stored_K_q(case, tv: Fraction):
    ctx = X.Ctx("Q", 1, tv)
    K = {}
    for (g, _), c in zip(X.nine_functions(ctx), X.coeff_strings_at(DATA["sym"][case], tv)):
        for p, x in c.items():
            K = X.fadd(K, X.fscale(g, x, p))
    return ctx, K


def scaled_stored(case, a, b):
    """Coefficients predicted from the frozen Q(t) ones by K(R; a, b) = a K(aR; 1, b/a) (f0: a^{p+3})."""
    return [{p: a ** (p + (3 if i == 0 else 1)) * x for p, x in c.items()}
            for i, c in enumerate(X.coeff_strings_at(DATA["sym"][case], b / a))]


def _nz(d):
    return {p: v for p, v in d.items() if v != 0}


# ================================================================ A. the equation (**)
def test_label13_log_decomposition_is_label13_verbatim():
    """A1.  The source of (**) uses label 13 split into four logarithms ln(t + c + eps w1); label 13 verbatim minus
    that split is a polynomial in t (its third t-derivative vanishes at a generic point, while P's own does not).

    Rejects: a mistranscribed label-13 term in the split (a sign in A3..A6, D3..D6, a wrong log argument) -- which
    the exact residual tests cannot see, since they solve whatever source they are given.
    """
    rem = M.P_full() - sum(C * sp.log(M.ts + c + e * M.w1s) for (C, c, e) in M.P_logs())
    pt = {M.u2s: sp.Rational(17, 10), M.u3s: sp.Rational(9, 10), M.w2s: sp.Rational(11, 10), M.w3s: sp.Rational(23, 10),
          M.w1s: sp.Rational(3, 10), M.ts: sp.Rational(7, 5)}
    val = abs(sp.N(sp.diff(rem, M.ts, 3).subs(pt), 40))
    ctl = abs(sp.N(sp.diff(M.P_full(), M.ts, 3).subs(pt), 40))
    assert ctl > 1e-3 and val < 1e-30, (val, ctl)


@pytest.mark.parametrize("a,b", PAIRS_STAGE1)
def test_source_and_f2_match_labels_61_63_at_the_exchange_point(a, b):
    """A2.  At the exchange point (base jet): S_1 = -(W_1^{(-2)})'' with W_1^{(-2)} typed from label 63, and
    f_2 = [2(a^2+b^2) f_0 - f_0'' + W_1^{(-2)}] / (a^2-b^2)^2 (label 61, n = 1).  Exact, at (1, 2) and at the
    unused pair (2/3, 7/4).

    Rejects: S_m -> -S_m; Sig -> -Sig; the d^4 term dropped; a wrong label-13 log term.  (Blind to c and e, which
    vanish at the base jet: see A5.)
    """
    ctx, D, calc, RJ, fl, nfree0, free = stage1("Q", a, b)
    W = label63_W1(ctx)
    assert X.fadd(calc.D2(W), RJ.source(1)[Z0]) == {}
    u, w = ctx.a, ctx.b
    lab61 = X.fscale(X.fadd(X.fadd(X.fscale(fl[0][Z0], 2 * (u * u + w * w)), calc.D2(fl[0][Z0]), -1), W),
                     1 / (u * u - w * w) ** 2)
    assert X.fadd(lab61, fl[1][Z0], -1) == {}


@pytest.mark.parametrize("a,b", PAIRS_STAGE1)
def test_f4_matches_label61_at_n3(a, b):
    """A3.  (a^2-b^2)^2 f_4'' = 18 (a^2+b^2) f_2'' - 9 f_2'''' - 12 f_0'' - S_3 at the exchange point, the
    coefficients 18, 9, 12 typed from label 61 (n = 3), not from (**).

    Rejects: m^2 -> m (9 -> 3); the f_{m-3}'' term dropped; its coefficient m(m-1)^2(m-2) = 12 -> m(m-1)(m-2) = 6.
    """
    ctx, D, calc, RJ, fl, nfree0, free = stage1("Q", a, b)
    u, w = ctx.a, ctx.b
    c18, c9, c12 = LABEL61_N3
    lhs = X.fscale(calc.D2(fl[2][Z0]), (u * u - w * w) ** 2)
    rhs = X.fadd(X.fscale(calc.D2(fl[1][Z0]), c18 * (u * u + w * w)), X.fscale(calc.D2(calc.D2(fl[1][Z0])), -c9))
    rhs = X.fadd(X.fadd(rhs, X.fscale(calc.D2(fl[0][Z0]), -c12)), RJ.source(3)[Z0], -1)
    assert X.fadd(lhs, rhs, -1) == {}


@pytest.mark.parametrize("a,b", PAIRS_STAGE1)
def test_star_residuals_vanish_with_coefficients_typed_here(a, b):
    """A4.  Every jet of f_2 and f_4 on the stage-1 downset (36 multi-indices) satisfies (**) differentiated by
    Leibniz, with q, c, Sig, e and the m-coefficients typed in this file; every solve left 0 free unknowns.

    Rejects: c -> -c, e dropped, high-|U| c-jets dropped, m^2 -> m, the f_{m-3}'' term dropped or 12 -> 6 -- in
    the helper's solver or in its scalar jets; a non-unique solve.
    """
    ctx, D, calc, RJ, fl, nfree0, free = stage1("Q", a, b)
    assert nfree0 == 0 and len(free) == 2 * len(D) and set(free) == {0}, free
    sc = own_scalar_jets(ctx, D)
    assert max(sum(U) for U in sc["c"]) >= 4                      # high-|U| c-jets are present and exercised
    assert residual_of_star(ctx, calc, RJ, D, 1, fl[1], fl[0], None, sc) == []
    assert residual_of_star(ctx, calc, RJ, D, 3, fl[2], fl[1], fl[0], sc) == []


@pytest.mark.parametrize("a,b", PAIRS_STAGE1)
def test_laplacian_control_f2_rung(a, b):
    """A5.  Control (a), exact: nabla^2 r12 = 2/r12 gives 16 f_{1111} = 8 [(a^2+b^2) f2_{1111} + 2a f2_{0111} +
    2b f2_{1011} + ab (f2_{2011} + f2_{0211}) - ab R^2 f2_{0011}] as an identity in the function algebra; f_0 from
    the reviewed engine, f_2 from (**), including transverse derivatives and alpha_i = 0, 2.

    Rejects: c -> -c; e dropped; c-jets with |U| >= 4 dropped; S_1 -> -S_1.  Does NOT reject m^2 -> m or any
    f_{m-3} error (m = 1 only): that is A3 and A6.
    """
    ctx, D, calc, RJ, fl, nfree0, free = stage1("Q", a, b)
    s16, s8 = CTRL_A
    assert X.fadd(X.fscale(fl[0][(1, 1, 1, 1)], ctx.c(s16)),
                  X.fscale(laplacian_combination(ctx, fl[1]), ctx.c(s8)), -1) == {}


@pytest.mark.parametrize("a,b", PAIRS_STAGE1)
def test_laplacian_control_f4_rung(a, b):
    """A6.  Control (a'), exact: nabla^2 r12^3 = 12 r12 gives 12 f2_{1111} = the same Laplacian combination of f_4
    jets.  The f_4 jets come from (**) at m = 3, so this is the control that sees the m-dependence.

    Rejects: m^2 -> m; the f_{m-3}'' term dropped; its 12 -> 6; c -> -c; e dropped.  Also a control written with
    6 for 12.
    """
    ctx, D, calc, RJ, fl, nfree0, free = stage1("Q", a, b)
    assert X.fadd(X.fscale(fl[1][(1, 1, 1, 1)], ctx.c(CTRL_A_PRIME)), laplacian_combination(ctx, fl[2]), -1) == {}
    assert X.in_nine(ctx, fl[1][(1, 1, 1, 1)]) is not None and X.in_nine(ctx, fl[2][(1, 1, 1, 1)]) is not None


def test_free_unknown_counter_can_see_a_free_unknown():
    """A7.  Positive control for "0 free unknowns" (A4, B1): the counter the solver reports is live.  (i) An
    under-determined condition system reports its free directions; (ii) a solve_D2 whose termination conditions are
    thinned out reports at least one free unknown instead of silently setting it to 0.

    Rejects: a solver that zeroes free unknowns without reporting them (a residual check would still pass with an
    arbitrary member of the solution family).
    """
    calc = M.MCalc(X.Ctx("Q", 1, 2))
    calc._solve_conds([[F(1), F(1), F(0), F(0)]], 3)
    assert calc.last_free == 2
    calc._solve_conds([], 2)
    assert calc.last_free == 2
    ctx, D, calc1, RJ, fl, nfree0, free = stage1("Q", F(1), F(2))
    src = calc1.D2(fl[1][(1, 1, 1, 1)])
    thin = M.MCalc(ctx)
    orig = thin._solve_conds
    thin._solve_conds = lambda conds, n: orig(conds[:1], n)
    thin.solve_D2(src)
    assert thin.last_free >= 1
    jets, fr = M.solve_fn(thin, RJ, [Z0], 1, fl[0])              # (iii) solve_fn passes the count on
    assert fr and fr[0] >= 1, fr
    full = M.MCalc(ctx)
    assert X.fadd(full.solve_D2(src), fl[1][(1, 1, 1, 1)], -1) == {} and full.last_free == 0


def test_gcd_fallback_is_exact():
    """A8.  The python-flint fallback for sympy's heuristic gcd (installed only for Q(t) runs): on structured and
    random integer polynomials it returns sympy's gcd up to sign and exact cofactors; installed, a forced heugcd
    failure is routed to it.

    Rejects: a fallback with wrong cofactors or a non-maximal gcd taken as exact; a patch that changes the result
    when heugcd succeeds.
    """
    import random

    import sympy.polys.heuristicgcd as H
    from sympy.polys.polyerrors import HeuristicGCDFailed
    from sympy.polys.rings import ring
    Rr, t = ring("t", sp.ZZ)
    rnd = random.Random(20261118)
    pairs = [((t - 1) ** 5 * (t + 1) * 6, (t - 1) ** 3 * (t ** 2 + 1) * 4), (t ** 7 - t, t ** 4 - 1),
             (t ** 3 + 2, t ** 5 + 3 * t + 1)]
    for _ in range(25):
        g = Rr.from_dense([rnd.randint(-9, 9) for _ in range(rnd.randint(1, 4))] or [1])
        f1 = Rr.from_dense([rnd.randint(-20, 20) for _ in range(rnd.randint(1, 8))])
        f2 = Rr.from_dense([rnd.randint(-20, 20) for _ in range(rnd.randint(1, 8))])
        if g and f1 and f2:
            pairs.append((g * f1, g * f2))
    for f, g in pairs:
        h, cf, cg = M._gcd_ZZ_flint(f, g)
        ref = sp.Poly(f.as_expr(), sp.Symbol("t")).gcd(sp.Poly(g.as_expr(), sp.Symbol("t")))
        hp = sp.Poly(h.as_expr(), sp.Symbol("t"))
        assert hp == ref or hp == -ref, (f, g, h, ref)
        assert h * cf == f and h * cg == g
    orig = H.heugcd

    def boom(f, g):
        raise HeuristicGCDFailed("forced")
    try:
        H.heugcd = boom
        f, g = pairs[0]
        assert M._gcd_ZZ_safe(f, g) == M._gcd_ZZ_flint(f, g)
    finally:
        H.heugcd = orig
    f, g = pairs[1]
    assert M._gcd_ZZ_safe(f, g) == orig(f, g)


@pytest.mark.slow
def test_controls_identically_in_t():
    """A9 (slow, about 4 min).  A2, A3, A5, A6 and the residuals of A4's form with the helper's scalar jets, over
    Q(t) (a = 1, b = t): the controls are identities of rational functions of t, not checks at sampled ratios.

    Rejects: an error in (**) that vanishes at the two sampled pairs (a coefficient wrong by a factor (t - 2)(3t - 7)).
    """
    ctx, D, calc, RJ, fl, nfree0, free = stage1("Qt")
    assert nfree0 == 0 and set(free) == {0}
    u, w = ctx.a, ctx.b
    W = label63_W1(ctx)
    assert X.fadd(calc.D2(W), RJ.source(1)[Z0]) == {}
    lab61 = X.fscale(X.fadd(X.fadd(X.fscale(fl[0][Z0], 2 * (u * u + w * w)), calc.D2(fl[0][Z0]), -1), W),
                     1 / (u * u - w * w) ** 2)
    assert X.fadd(lab61, fl[1][Z0], -1) == {}
    s16, s8 = CTRL_A
    assert X.fadd(X.fscale(fl[0][(1, 1, 1, 1)], ctx.c(s16)),
                  X.fscale(laplacian_combination(ctx, fl[1]), ctx.c(s8)), -1) == {}
    assert X.fadd(X.fscale(fl[1][(1, 1, 1, 1)], ctx.c(CTRL_A_PRIME)), laplacian_combination(ctx, fl[2]), -1) == {}


# ================================================================ B. the two cases over Q
@pytest.mark.parametrize("case,a,b", [("2ppi", F(1), F(2)), ("2ppi", F(2, 3), F(7, 4)),
                                      _slow("3ddel", F(1), F(2)), _slow("3ddel", F(2, 3), F(7, 4))])
def test_closed_form_in_nine_function_span_closed_under_d_dR(case, a, b):
    """B1.  Over Q: every solve of (**) leaves 0 free unknowns, every f_{2k} jet the case needs satisfies (**) with
    the coefficients typed here, K is in the span of the nine functions with all nine present, K..K^(10) stay in it,
    and their coordinate rank at R = 7/5 is 9.

    Rejects: a non-unique solve; an m != 0 integral outside the nine-function span (a tenth function); a K missing
    one of the nine (the order bound would then be about a smaller span).
    """
    ctx, calc, fl, nf0, free, K = solved(case, a, b)
    Dk = M.downsets(case)
    assert {k: len(d) for k, d in Dk.items()} == N_INDICES[case]
    assert nf0 == 0 and free and set(free) == {0}, free
    for k in range(1, max(Dk) + 1):
        RJ = M.RJets(ctx, Dk[k], calc)
        sc = own_scalar_jets(ctx, Dk[k])
        assert residual_of_star(ctx, calc, RJ, Dk[k], 2 * k - 1, fl[k], fl[k - 1], fl[k - 2] if k >= 2 else None,
                                sc) == [], (case, k)
    co = X.in_nine(ctx, K)
    assert co is not None and all(co), "K outside the span, or one of the nine absent"
    ders = calc.derivatives(K, 10)
    assert all(X.in_nine(ctx, A) is not None for A in ders)
    Mx = PX.coordinate_matrix_at(ders, PX.all_keys(ders), F(7, 5))
    assert PX.rank_q(Mx[:9]) == 9 and PX.rank_q(Mx) == 9


@pytest.mark.parametrize("case", ["2ppi", _slow("3ddel")])
def test_operator_at_1_2_is_unique_and_equals_the_committed_one(case):
    """B2.  At (a, b) = (1, 2): exactly one operator of order <= 9 and degree <= 25 (2ppi) / 35 (3ddel) annihilates
    K, none of degree one lower, none of order <= 8 up to degree + 12; it equals the committed operator OP_P_1_2
    (tests/test_paper58_exchange_l1.py) / OP_3DDEL_1_2 (tests/test_paper58_exchange_d_pp.py), every coefficient.

    Rejects: a closed form of the wrong integral (the committed operators were guessed from exact Neumann-series
    data and never saw this route); a sigma/pi mix-up; the wrong sign of the r12^2 term in X; a shifted alpha.
    """
    ctx, calc, fl, nf0, free, K = solved(case, F(1), F(2))
    d = OP_DEGREE[case]
    ders = calc.derivatives(K, 9)
    ops = PX.operators_annihilating(ders, 9, d)
    assert len(ops) == 1, len(ops)
    assert PX.operators_annihilating(ders, 9, d - 1) == []
    assert PX.operators_annihilating(ders[:9], 8, d + 12) == []
    assert X.primitive(ops[0]) == committed_operator(case)


@pytest.mark.parametrize("case", ["2ppi", _slow("3ddel")])
def test_committed_operators_annihilate_their_own_family_only(case):
    """B3.  Applied as exact identities in the function algebra at (1, 2): the committed operator of the case
    annihilates K; the committed operator of the OTHER family leaves a nonzero residual.

    Rejects: a wrong-family operator (OP_3DDEL_1_2 passed for 2ppi or the reverse); a closed form that every
    order-9 operator of this shape would annihilate.
    """
    ctx, calc, fl, nf0, free, K = solved(case, F(1), F(2))
    other = {"2ppi": "3ddel", "3ddel": "2ppi"}[case]
    assert apply_operator(calc, K, committed_operator(case)) == {}
    res = apply_operator(calc, K, committed_operator(other))
    assert len(res) >= 5, len(res)


@pytest.mark.parametrize("case,a,b,R", [("2ppi", F(1), F(2), F(7, 5)), _slow("3ddel", F(1), F(2), F(7, 5)),
                                        ("2ppi", F(4, 5), F(3, 2), F(23, 10)),
                                        _slow("3ddel", F(4, 5), F(3, 2), F(23, 10))])
def test_values_at_the_derivation_points(case, a, b, R):
    """B4.  K(R) from the live closed form (f0 by its one quadrature, dps 60) equals the derivation's value to
    1e-45 (that value agreed with exchange_hp, an independent Neumann expansion, to 47-48 digits).

    Rejects: a wrong normalisation (a^7/12 -> a^7/18); 2X^2 - rho^2 rho^2 -> X^2; the m-errors of A3 for 3ddel
    (value moves to -0.24, 0.11, 0.05); a wrong power of a (visible at a = 4/5 only).
    """
    v = value(case, a, b, R)
    with mp.workdps(60):
        ref = mp.mpf(VAL[(case, a, b, R)])
        assert abs(v / ref - 1) < mp.mpf(TOL_MEMO), mp.nstr(v, 50)


@pytest.mark.parametrize("case", ["2ppi", _slow("3ddel")])
def test_values_at_the_unused_pair_against_exchange_hp_digits(case):
    """B5.  At (a, b, R) = (2/3, 7/4, 19/10), a point the derivation never used: the live closed form against the
    review's exchange_hp values (tau_max 16, dps 32), typed here: 2ppi -0.00515968878885911303770792463991, 3ddel
    0.000355893627169501579264221099214; agreement to 1e-29 relative.

    Rejects: an error tuned away at the derivation's points; a wrong b-dependence or a-power together.
    """
    a, b, R = UNUSED
    v = value(case, a, b, R)
    with mp.workdps(60):
        ref = mp.mpf(VAL_UNUSED[case])
        assert abs(v / ref - 1) < mp.mpf(TOL_UNUSED), mp.nstr(v, 40)


@pytest.mark.slow
@pytest.mark.parametrize("case", CASES)
def test_values_at_the_unused_pair_against_live_exchange_hp(case):
    """B6 (slow; 2ppi about 50 s, 3ddel about 90 s).  The same point against geovac.qfd_core.exchange_hp computed
    here (tau_max 12, dps 25; measured truncation error <= 2e-21 relative), not typed digits.

    Rejects: typed reference digits that were copied from the closed form rather than computed independently.
    """
    from geovac import qfd_core as Q
    a, b, R = UNUSED
    A, B = sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator)
    with mp.workdps(25):
        if case == "2ppi":
            hp, per = Q.exchange_hp(2 * A, (2, 1, 1), (1, 0, 0), B, (2, 1, -1), (1, 0, 0), sp.Rational(19, 10),
                                    tau_max=12, dps=25)
        else:
            hp, per = Q.exchange_hp(3 * A, (3, 2, 2), (1, 0, 0), B, (3, 2, -2), (1, 0, 0), sp.Rational(19, 10),
                                    tau_max=12, dps=25)
        hp = +hp
    v = value(case, a, b, R)
    with mp.workdps(40):
        assert abs(v / hp - 1) < mp.mpf(10) ** -19, (mp.nstr(v, 30), mp.nstr(hp, 30))


@pytest.mark.parametrize("a,b", [(F(1), F(2)), (F(2, 3), F(7, 4))])
def test_2ppi_f0_coefficient_closed_form(a, b):
    """B7.  The coefficient of f0 in 2ppi is 2 a^7 b^5 [225 - 180(a^2+b^2)R^2 - (102a^4+336a^2b^2+102b^4)R^4
    - (24a^6+176a^4b^2+176a^2b^4+24b^6)R^6 + (a^8+16a^6b^2+36a^4b^4+16a^2b^6+b^8)R^8] / (R^3 (b^2-a^2)^6), exactly.

    Rejects: a wrong normalisation or scaling in a or b (two pairs, a != 1 in the second).
    """
    ctx, calc, fl, nf0, free, K = solved("2ppi", a, b)
    R_, a_, b_ = sp.symbols("R a b")
    form = 2 * a_**7 * b_**5 * (225 - 180 * (a_**2 + b_**2) * R_**2 - (102 * a_**4 + 336 * a_**2 * b_**2 + 102 * b_**4) * R_**4
                                - (24 * a_**6 + 176 * a_**4 * b_**2 + 176 * a_**2 * b_**4 + 24 * b_**6) * R_**6
                                + (a_**8 + 16 * a_**6 * b_**2 + 36 * a_**4 * b_**4 + 16 * a_**2 * b_**6 + b_**8) * R_**8) \
        / (R_**3 * (b_**2 - a_**2)**6)
    mine = sum(sp.Rational(c.numerator, c.denominator) * R_**p for p, c in K[("f",)].items())
    ref = form.subs({a_: sp.Rational(a.numerator, a.denominator), b_: sp.Rational(b.numerator, b.denominator)})
    assert sp.expand(mine - ref) == 0


@pytest.mark.slow
def test_3ddel_f0_coefficient_leading_terms():
    """B8 (slow).  The coefficient of f0 in 3ddel starts a^9 b^5 [99225 - 28350 (a^2+b^2) R^2 + ...] /
    (12 R^5 (b^2-a^2)^8): its R^-5 and R^-3 Laurent coefficients at (2/3, 7/4).

    Rejects: a wrong 3ddel normalisation or power of a, b.
    """
    a, b = UNUSED[:2]
    ctx, calc, fl, nf0, free, K = solved("3ddel", a, b)
    den = 12 * (b * b - a * a) ** 8
    assert K[("f",)][-5] == a ** 9 * b ** 5 * 99225 / den
    assert K[("f",)][-3] == -a ** 9 * b ** 5 * 28350 * (a * a + b * b) / den


# ================================================================ C. every 0 < a < b
@pytest.mark.parametrize("case", CASES)
def test_frozen_closed_form_poles_only_at_t_0_1_minus1(case):
    """C1.  The frozen Q(t) coefficients of K on the nine functions: denominators t, t - 1, t + 1 only, all nine
    present; pole order at t = 1: 7 (2ppi), 9 (3ddel).

    Rejects: a closed form with a pole at some b/a > 1 (e.g. a coefficient over t - 3 or 2t - 3).
    """
    t = sp.Symbol("t")
    mult = {}
    for nm in X.NAMES9:
        assert DATA["sym"][case][nm], nm
        for p, s in DATA["sym"][case][nm].items():
            for f, m in sp.factor_list(sp.denom(sp.together(sp.sympify(s, locals={"t": t}))))[1]:
                mult[str(f)] = max(mult.get(str(f), 0), int(m))
    assert set(mult) == {"t", "t - 1", "t + 1"}, mult
    assert mult["t - 1"] == MAX_POLE[case]


@pytest.mark.parametrize("case", CASES)
def test_frozen_audit_of_every_intermediate_jet(case):
    """C2.  The frozen audit of the Q(t) runs (every jet of f_0, f_2, f_4 the case needs): coefficient denominators
    factor into t, t - 1, t + 1 only; every function key decays for every t > 1 with no vanishing Ei argument,
    re-judged here by the reviewed rule X.key_verdict.  Uniqueness by decay is used at every multi-index, so the
    identity over Q(t) holds at every real t > 1 only if this holds.

    Rejects: an intermediate with a pole or a growing key at some t > 1 that cancels in K.  (Does not reject a
    wrong record: C3, C4 compare it with live runs.)
    """
    rec = DATA["audit"][case]
    assert rec["multi_indices"] == {f"f_{2 * k}": n for k, n in N_INDICES[case].items()}
    for lvl, dens in rec["denominator_factors_max_multiplicity"].items():
        assert set(dens) <= {"t", "t - 1", "t + 1"}, (lvl, dens)
    for r in rec["keys"]:
        v = X.key_verdict(r)
        assert v["decays"] and v["ei_zero"] == [], (r, v)


@pytest.mark.parametrize("case", ["2ppi", _slow("3ddel")])
def test_frozen_keys_are_the_live_keys_at_a_specialisation(case):
    """C3.  The frozen key list at t = 21/8, rates scaled by a = 2/3, is exactly the set of keys of the live f_0,
    f_2 (and f_4) jets at (a, b) = (2/3, 7/4).

    Rejects: a frozen audit that omits a key of an intermediate (C2 would pass by omission).
    """
    a, b = UNUSED[:2]
    ctx, calc, fl, nf0, free, K = solved(case, a, b)
    live = {k for d in fl.values() for A in d.values() for k in A}
    t = sp.Symbol("t")

    def val(s):
        v = sp.Rational(sp.sympify(s, locals={"t": t}).subs(t, sp.Rational(21, 8)))
        return a * F(int(v.p), int(v.q))
    stored = set()
    for r in DATA["audit"][case]["keys"]:
        if r["kind"] == "f":
            stored.add(("f",))
        elif r["kind"] == "E":
            stored.add(("E", val(r["lam"]), r["log_power"]))
        else:
            stored.add(("Ei", val(r["lam"]), val(r["mu"])))
    assert live == stored and len(stored) == len(DATA["audit"][case]["keys"])


@pytest.mark.parametrize("case,tv", [("2ppi", F(2)), ("2ppi", F(3)), _slow("3ddel", F(2)), _slow("3ddel", F(3))])
def test_frozen_coefficients_equal_live_solves(case, tv):
    """C4.  The frozen Q(t) coefficients at t = 2 and t = 3 equal a live exact solve over Q at (1, t), coefficient
    by coefficient (t = 3: two Ei arguments coincide).

    Rejects: frozen coefficients that are not the engine's, or right at t = 2 only.
    """
    ctx, calc, fl, nf0, free, K = solved(case, F(1), tv)
    co = X.in_nine(ctx, K)
    assert co is not None
    assert [_nz(c) for c in co] == X.coeff_strings_at(DATA["sym"][case], tv)


@pytest.mark.parametrize("case,a,b", [("2ppi", F(2, 3), F(7, 4)), _slow("3ddel", F(2, 3), F(7, 4)),
                                      _slow("3ddel", F(4, 5), F(3, 2))])
def test_scaling_law_and_frozen_coefficients_at_a_not_1(case, a, b):
    """C5.  At a != 1 the live coefficients on the nine functions are a^{p+1} c_p(b/a) (f0: a^{p+3} c_p(b/a)) with
    c_p the frozen Q(t) ones, i.e. K(R; a, b) = a K(aR; 1, b/a).

    Rejects: a wrong power of a in a normalisation (invisible at a = 1); frozen coefficients wrong at t = 21/8.
    """
    ctx, calc, fl, nf0, free, K = solved(case, a, b)
    co = X.in_nine(ctx, K)
    assert [_nz(c) for c in co] == scaled_stored(case, a, b)
    assert [{p: a * x for p, x in c.items()} for c in co] != scaled_stored(case, a, b)


@pytest.mark.slow
def test_2ppi_Qt_run_recomputed():
    """C6 (slow, about 3 min).  The 2ppi run over Q(t) from scratch: 0 free unknowns, (**) residuals of f_2
    identically 0 in t with the scalar jets of the helper (A9 pins them), the audit of every intermediate equals the
    frozen one, and K's coefficients are the frozen ones.

    Rejects: frozen 2ppi data the engine does not reproduce identically in t.
    """
    M.install_gcd_fallback()
    ctx = X.Ctx("Qt")
    calc, fl, nf0, free = M.solve_case(ctx, "2ppi")
    assert nf0 == 0 and set(free) == {0}
    K = M.assemble_K(ctx, "2ppi", fl)
    co = X.in_nine(ctx, K)
    assert co is not None
    t = ctx.t
    for nm, c in zip(X.NAMES9, co):
        st = DATA["sym"]["2ppi"][nm]
        assert {str(p) for p in c} == set(st), nm
        for p, x in c.items():
            assert x == ctx.K.from_sympy(sp.sympify(st[str(p)], locals={"t": t})), (nm, p)
    dens, keys = {}, set()
    for k, d in fl.items():
        acc = dens.setdefault(f"f_{2 * k}", {})
        for A in d.values():
            for key, dd in A.items():
                keys.add(key)
                for cc in dd.values():
                    den = sp.Poly(ctx.K.to_sympy(cc.denom), t)
                    if den.degree() > 0:
                        for f, m in sp.factor_list(den.as_expr())[1]:
                            acc[str(f)] = max(acc.get(str(f), 0), int(m))
    assert dens == DATA["audit"]["2ppi"]["denominator_factors_max_multiplicity"]
    assert len(keys) == len(DATA["audit"]["2ppi"]["keys"])


# ================================================================ D. exactly 9: leg (D)
@lru_cache(maxsize=None)
def leg_d(case: str):
    ctx = X.Ctx("Qt")
    rows = EO.derivative_rows(ctx, EO.stored_K(ctx, DATA["sym"][case]), 8)
    per_row, other = EO.denominator_audit(rows)
    return ctx, rows, per_row, other


@lru_cache(maxsize=None)
def numerators(case: str):
    ctx, rows, per_row, other = leg_d(case)
    cleared = EO.clear_rows(ctx, rows, per_row)
    return [(R0,) + EO.determinant_polynomial(cleared, R0) for R0 in R_VALUES]


def common_root_report(polys):
    stripped = [EO.strip_linear(list(N)) for N in polys]
    g = EO.common_factor([r for _, r in stripped])
    return EO.real_roots_above_one(g), g.degree(), [m for m, _ in stripped]


@pytest.mark.parametrize("case", CASES)
def test_derivative_rows_have_denominators_t_tm1_tp1_only(case):
    """D1.  Every Laurent coefficient of every entry of the 9 x 9 coordinate matrix of (K, ..., K^(8)) over Q(t)
    has denominator a constant times t^i (t - 1)^j (t + 1)^k (the gcd argument needs no pole at t0 > 1).

    Rejects: a denominator factor other than t, t - 1, t + 1 in a frozen coefficient or produced by the rows.
    """
    ctx, rows, per_row, other = leg_d(case)
    assert len(rows) == 9 and all(len(rw) == 9 for rw in rows)
    assert other == [], sorted(set(other))[:5]
    assert all(m[1] > 0 and m[2] > 0 for m in per_row)


@pytest.mark.parametrize("case", CASES)
def test_determinant_numerators_have_no_common_root_above_one(case):
    """D2.  Leg (D): N_j(t) = cleared numerator of det M(t, R_j), exact integer polynomials (Bareiss at integer t,
    Newton interpolation, six extra nodes) for R_j = 2/3, 13/4, 6.  A ratio t0 > 1 where det M(t0, R) vanished
    identically in R would be a root of every N_j; the gcd (and each pairwise gcd) has no real root in (1, oo).
    The design needs at least three R: at each single R the numerator HAS real roots t > 1 (asserted here for one R
    at least), so a single-R test is meaningless.

    Rejects: a bad exponent ratio (rank < 9 for all R at some t0 > 1, e.g. the f0 coordinate times t - 2); rank 8
    for all t (f0 coordinate removed: every N_j is 0); a leg (D) run at one value of R.
    """
    assert len(R_VALUES) >= 3 and len(set(R_VALUES)) == len(R_VALUES)
    polys = [N for _, N, _ in numerators(case)]
    for (R0, N, _) in numerators(case):
        assert any(N), f"det M(t, {R0}) identically 0"
    nroots, gdeg, mults = common_root_report(polys)
    assert nroots == 0, nroots
    for i in range(len(polys)):
        for j in range(i + 1, len(polys)):
            assert common_root_report([polys[i], polys[j]])[0] == 0
    single = [EO.real_roots_above_one(EO.common_factor([EO.strip_linear(list(N))[1]])) for N in polys]
    assert max(single) >= 1, single


def test_determinant_detector_sees_a_planted_bad_ratio_and_a_rank_drop():
    """D0.  Positive controls for D2's machinery on 2ppi: the f0 coordinate times (t - 2) must give a common root
    at t = 2; the f0 coordinate removed must give zero determinants.

    Rejects: a D2 that reports "no common root" whatever it is given.
    """
    ctx = X.Ctx("Qt")
    K = EO.stored_K(ctx, DATA["sym"]["2ppi"])
    rest = {k: v for k, v in K.items() if k != ("f",)}
    for factor in (ctx.b - 2, None):
        Kp = rest if factor is None else X.fadd(rest, X.fscale({("f",): K[("f",)]}, factor))
        rows = EO.derivative_rows(ctx, Kp, 8)
        per_row, other = EO.denominator_audit(rows)
        cleared = EO.clear_rows(ctx, rows, per_row)
        polys = [EO.determinant_polynomial(cleared, R0)[0] for R0 in R_VALUES[:2]]
        if factor is None:
            assert all(not any(N) for N in polys)
        else:
            nroots, gdeg, _ = common_root_report(polys)
            assert nroots == 1
            assert EO.common_factor([EO.strip_linear(list(N))[1] for N in polys]).eval(2) == 0


@pytest.mark.parametrize("case,tv", [("2ppi", F(3)), _slow("3ddel", F(3)), _slow("3ddel", F(7, 3))])
def test_determinant_polynomials_equal_fresh_solve_determinants(case, tv):
    """D3.  The N_j of D2 are the determinants they are said to be: rows rebuilt from a FRESH exact solve over Q at
    (1, t) (no frozen data, no Q(t)), determinant by a different elimination, N_j(t) = scale x clearing x det for
    every R_j; det != 0 for at least one R_j.

    Rejects: an interpolation degree bound too low, wrong clearing exponents, rows from a closed form right at
    t = 2 only.
    """
    qt, rows_t, per_row, other = leg_d(case)
    ctx, calc, fl, nf0, free, K = solved(case, F(1), tv)
    rows = EO.derivative_rows(ctx, K, 8)
    clearing = F(1)
    for (e0, e1, e2) in per_row:
        clearing *= tv ** e0 * (tv - 1) ** e1 * (tv + 1) ** e2
    nonzero = 0
    for (R0, N, scale) in numerators(case):
        det = EO.fraction_det([[sum((F(x) * R0 ** p for p, x in cd.items()), F(0)) for cd in rw] for rw in rows])
        assert sum((F(c) * tv ** k for k, c in enumerate(N)), F(0)) == scale * clearing * det, (case, str(R0))
        nonzero += det != 0
    assert nonzero >= 1


# ================================================================ E. the sigma = 2 part of real 2p_x 2p_x exchange
def _norm2(s: int, z):
    """N^2 for N (x + i y)^s e^{-z r}: 1 / int |x + i y|^{2s} e^{-2 z r} d^3r."""
    r, th = sp.symbols("r theta", positive=True)
    I = 2 * sp.pi * sp.integrate(r ** (2 * s + 2) * sp.exp(-2 * z * r), (r, 0, sp.oo)) \
        * sp.integrate(sp.sin(th) ** (2 * s + 1), (th, 0, sp.pi))
    return sp.simplify(1 / I)


def test_sigma2_identities_from_first_principles():
    """E1.  Fast, from definitions only.  (i) ppdel = (2p+1_A 2p-1_B | 2p-1_A 2p+1_B) and 3ddel = (3d+2_A 1s_B |
    3d-2_A 1s_B) have the same spatial integrand (x1+iy1)^2 (x2-iy2)^2 e^{..}/r12, so their ratio is the ratio of
    normalisations N_2p(a)^2 N_2p(b)^2 / (N_3d(a)^2 N_1s(b)^2), computed here by integration: 3 (b/a)^2, the
    helper's SIGMA2_RATIO.  (ii) (2p_x 2p_x | 2p_x 2p_x) = pppi0 + ppdel/2: for every Fourier mode e^{ik(phi1-phi2)}
    of 1/r12 (k = -4..4), the real weight N_x^4 cos^2 phi1 cos^2 phi2 integrates like N_+^4 (1 + cos 2(phi1-phi2)/2),
    with N_x^2 = z^5/pi, N_+^2 = z^5/(2 pi) from (i)'s integrals.

    Rejects: the identity factor 3 -> 2 or 4; (b/a)^2 -> (b/a); ppdel instead of ppdel/2; a 2p_x normalisation
    taken equal to the complex one.
    """
    a, b = sp.symbols("a b", positive=True)
    n2p = lambda z: _norm2(1, z)                       # noqa: E731
    n3d = lambda z: _norm2(2, z)                       # noqa: E731
    n1s = lambda z: _norm2(0, z)                       # noqa: E731
    ratio = sp.simplify(n2p(a) * n2p(b) / (n3d(a) * n1s(b)))
    assert sp.simplify(ratio - M.SIGMA2_RATIO(a, b)) == 0, ratio
    assert sp.simplify(ratio - 3 * b / a) != 0 and sp.simplify(ratio - 2 * (b / a) ** 2) != 0
    r, th = sp.symbols("r theta", positive=True)
    nx = lambda z: 1 / (sp.pi * sp.integrate(r ** 4 * sp.exp(-2 * z * r), (r, 0, sp.oo))      # noqa: E731
                        * sp.integrate(sp.sin(th) ** 3, (th, 0, sp.pi)))
    p1, p2 = sp.symbols("phi1 phi2", real=True)
    for k in range(-4, 5):
        mode = sp.exp(sp.I * k * (p1 - p2))
        lhs = nx(a) * nx(b) * sp.integrate(sp.integrate(sp.cos(p1) ** 2 * sp.cos(p2) ** 2 * mode, (p1, 0, 2 * sp.pi)),
                                           (p2, 0, 2 * sp.pi))
        rhs = n2p(a) * n2p(b) * sp.integrate(sp.integrate((1 + sp.cos(2 * (p1 - p2)) / 2) * mode, (p1, 0, 2 * sp.pi)),
                                             (p2, 0, 2 * sp.pi))
        assert sp.simplify(lhs - rhs) == 0, k


@pytest.mark.slow
def test_ppdel_against_exchange_hp():
    """E2 (slow, about 90 s).  ppdel from exchange_hp (tau_max 12, dps 25) at the unused point (2/3, 7/4, 19/10)
    equals SIGMA2_RATIO(a, b) x the live 3ddel closed form, to 1e-18 relative.

    Rejects: the factor 3 -> anything else, (b/a)^2 -> (b/a), against an independent evaluator.
    """
    from geovac import qfd_core as Q
    a, b, R = UNUSED
    A, B = sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator)
    with mp.workdps(25):
        hp, per = Q.exchange_hp(2 * A, (2, 1, 1), (2, 1, -1), 2 * B, (2, 1, -1), (2, 1, 1), sp.Rational(19, 10),
                                tau_max=12, dps=25)
        hp = +hp
    v = value("3ddel", a, b, R)
    with mp.workdps(40):
        pred = _mpq(M.SIGMA2_RATIO(a, b)) * v
        assert abs(pred / hp - 1) < mp.mpf(10) ** -18, (mp.nstr(pred, 25), mp.nstr(hp, 25))


@pytest.mark.slow
def test_real_2px_exchange_decomposes_into_pppi0_plus_half_ppdel():
    """E3 (slow, several minutes).  Real (2p_x,A 2p_x,B | 2p_x,A 2p_x,B) at (1, 2), R = 7/5 from the 16 complex
    exchange_hp terms of 2p_x = (2p-1 - 2p+1)/sqrt 2, against pppi0 + ppdel/2 from exchange_hp (tau_max 8, dps 22;
    measured difference 1.7e-19 on 0.15).

    Rejects: (xx|xx) = pppi0 + ppdel (off by 0.008), or a sign error in the sigma = 2 part.
    """
    from geovac import qfd_core as Q
    A, B, R = sp.Integer(1), sp.Integer(2), sp.Rational(7, 5)
    with mp.workdps(22):
        c = {-1: 1 / mp.sqrt(2), 1: -1 / mp.sqrt(2)}
        tot = mp.mpf(0)
        for ma, mb, mc, md in itertools.product((-1, 1), repeat=4):
            v, _ = Q.exchange_hp(2 * A, (2, 1, ma), (2, 1, mb), 2 * B, (2, 1, mc), (2, 1, md), R, tau_max=8, dps=22)
            tot += c[ma] * c[mb] * c[mc] * c[md] * v
        pppi0, _ = Q.exchange_hp(2 * A, (2, 1, 1), (2, 1, 1), 2 * B, (2, 1, 1), (2, 1, 1), R, tau_max=8, dps=22)
        ppdel, _ = Q.exchange_hp(2 * A, (2, 1, 1), (2, 1, -1), 2 * B, (2, 1, -1), (2, 1, 1), R, tau_max=8, dps=22)
        assert abs(ppdel) > mp.mpf("0.01")
        assert abs(tot - pppi0 - ppdel / 2) < mp.mpf(10) ** -15 * abs(tot), mp.nstr(tot - pppi0 - ppdel / 2, 5)


# ================================================================ FIRE TESTS (2026-11-18)
# Written and fired as a separate step after the guard.  Each plant was applied in memory by a pytest plugin that
# patches the imported modules after collection (no committed file was edited), the named selection was run with -x,
# and the file was then run clean.  (H) = tests/_exchange_m_nonzero.py, (E) = tests/_pachucki_extension.py,
# (J) = the frozen data as loaded, (T) = a constant of this file.  FIRED = the selection failed.
#
#   plant                                              -k selector                                    result
#   m^2 -> m (H)                                       f4_matches_label61 / star_residuals / laplacian_control_f4   FIRED x3
#     same                                             laplacian_control_f2                           NOT FIRED (expected: m = 1 only)
#     same                                             values_at_the_derivation_points and 2ppi       NOT FIRED (expected: 2ppi has no f_4)
#     same, slow                                       values_at_the_derivation_points and 3ddel and R1 --slow     FIRED
#   f_{m-3}'' term dropped (H)                         f4_matches_label61 / star_residuals / laplacian_control_f4   FIRED x3
#   f_{m-3}'' coefficient 12 -> 6 (H)                  f4_matches_label61 / star_residuals / laplacian_control_f4   FIRED x3
#     same                                             laplacian_control_f2                           NOT FIRED (expected)
#   control written with 6 for 12 (T)                  laplacian_control_f4                           FIRED
#   label 61 n = 3 typed 18, 9, 6 (T)                  f4_matches_label61                             FIRED
#   c -> -c in the solver (H)                          laplacian_control_f2 / star_residuals          FIRED x2
#     same                                             source_and_f2                                  NOT FIRED (expected: c = 0 at the base jet)
#     same, over Q(t), slow                            controls_identically_in_t --slow               FIRED
#   c -> -c in the helper's scalar jets (H)            star_residuals / laplacian_control_f2          FIRED x2
#   e term dropped (H)                                 laplacian_control_f2 / star_residuals          FIRED x2
#     same                                             source_and_f2                                  NOT FIRED (expected: e = 0 at the base jet)
#   c-jets with |U| >= 4 dropped (H)                   laplacian_control_f2 / star_residuals          FIRED x2
#     same                                             values_at_the_derivation_points and 2ppi       NOT FIRED (expected, as the review found)
#     same, slow                                       values_at_the_derivation_points and 3ddel and R1 --slow     FIRED
#   S_m -> -S_m; Sig -> -Sig; d^4 dropped (H)          source_and_f2                                  FIRED x3
#   _solve_conds silently zeroes free unknowns (E)     free_unknown_counter                           FIRED
#   solve_fn reports 0 free whatever the solver says   free_unknown_counter                           FIRED
#   label 13 split: A5 sign flipped (H)                label13                                        FIRED
#   gcd fallback returns a wrong cofactor (H)          gcd_fallback                                   FIRED
#   r12^2 sign in X (H)                                operator_at_1_2 and 2ppi / values_at_the_derivation_points and 2ppi   FIRED x2
#   2ppi normalisation a^5 -> a^4 (H)                  values .. 2ppi and R2 / unused_pair_against_exchange_hp_digits and 2ppi /
#                                                      scaling_law and 2ppi / f0_coefficient_closed_form                     FIRED x4
#     same                                             values_at_the_derivation_points and 2ppi and R0 NOT FIRED (expected: a = 1)
#   3ddel normalisation a^7/12 -> a^7/18 (H), slow     values_at_the_derivation_points and 3ddel and R1 --slow     FIRED
#   3ddel density 2X^2 - rho^2 rho^2 -> X^2 (H), slow  values_at_the_derivation_points and 3ddel and R1 --slow     FIRED
#   wrong family: OP_3DDEL_1_2 for 2ppi and reverse (T) committed_operators_annihilate and 2ppi / operator_at_1_2 and 2ppi   FIRED x2
#     same, slow                                       committed_operators_annihilate and 3ddel --slow             FIRED
#   SIGMA2_RATIO 3(b/a)^2 -> 3(b/a); -> 2(b/a)^2 (H)   sigma2_identities                              FIRED x2
#     3(b/a), slow                                     ppdel_against --slow                           FIRED
#   frozen coefficient / (t - 3), 2ppi (J)             frozen_closed_form_poles and 2ppi / derivative_rows and 2ppi   FIRED x2
#   frozen coefficient / (t - 3), 3ddel (J)            derivative_rows and 3ddel                      FIRED
#   frozen f0 coordinate x (t - 2), 2ppi, 3ddel (J)    no_common_root and <case>                      FIRED x2
#   frozen coefficient x (1 + (t-2)(t-3)), 2ppi (J)    frozen_coefficients_equal_live_solves and 2ppi NOT FIRED (expected: right at t = 2, 3)
#     same                                             scaling_law and 2ppi                           FIRED (t = 21/8)
#     same, 3ddel, slow                                scaling_law and 3ddel --slow                   FIRED
#   frozen audit: growing key e^{(t-3)R}, both (J)     frozen_audit and <case>                        FIRED x2
#   frozen audit: denominator t - 3, 3ddel (J)         frozen_audit and 3ddel                         FIRED
#   frozen audit: a key omitted, 2ppi (J)              frozen_keys_are_the_live_keys and 2ppi         FIRED
#   leg (D) at a single R (T)                          no_common_root                                 FIRED
#
# NOT PINNED, and why.
#   * A non-unique solve can only be detected through the reported count (A7 shows the count is live); no physical
#     case has produced one.
#   * In the fast tier the 3ddel-only wrong answers (normalisation, density, wrong-family on 3ddel) are slow-only;
#     the m = 3 coefficients ARE covered fast through A3, A4, A6.
#   * The frozen 3ddel Q(t) record is checked against live solves at t = 2, 3, 7/3 (determinants), 21/8 and 15/8,
#     not recomputed over Q(t).
