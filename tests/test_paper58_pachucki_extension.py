"""Paper 58 -- the Pachucki exchange bound at general exponents, in the limit b -> a, and for l > 0.

Guard for the three sentences proposed for Paper 58 (papers/group2_quantum_chemistry/paper_58_abelian_residue.tex):

  1. 1s exchange K(R; a, b): with the exponent ratio t = b/a kept as an indeterminate, the fifteen derivative
     equations of Pachucki's second-order equation at the exchange point are solved, as identities of rational
     functions of t, in the span of decaying functions; K lies in the span of nine functions with coefficients in
     Q(a, b)(R) whose denominators are powers of b - a and b + a.  Hence order <= 9 for every 0 < a < b.
  2. As b -> a the poles in b - a cancel and the closed form tends to Sugiura's equal-exponent formula, in the
     span of e^{-2R}, e^{-2R}(gamma + ln R), E1(2R), e^{2R} E1(4R).
  3. The same derivation, with the l > 0 orbital products written as polynomials in r_A^2, r_B^2 and mapped to
     higher exponent derivatives of the master integral, places four m = 0 exchange integrals in the same span for
     all 0 < a < b:  2psig = (2p0_A 1s_B | 2p0_A 1s_B),  3dsig = (3d0_A 1s_B | 3d0_A 1s_B),
     ppsig = (2p0_A 2p0_B | 2p0_A 2p0_B),  pppi0 = the exchange integral of the COMPLEX pair 2p+1_A, 2p+1_B
     (orbital product (1/2) N_A N_B (x^2 + y^2) e^{..}; record EQ's `pppi0`; NOT the real 2p_x-2p_x exchange).
     At (a, b) = (1, 2) the derived operators equal the guessed ones over Q: 2psig (9, 34), ppsig (9, 53),
     pppi0 (9, 39); 3dsig, unrecorded there, gives a unique (9, 53).

Written as a SEPARATE pass from the derivation (CLAUDE.md Sec. 9).  Prior art: K. Pachucki, Phys. Rev. A 80, 032520
(2009), arXiv:0905.1635; Y. Sugiura, Z. Phys. 45, 484 (1927) for the equal-exponent formula.  The engine is
tests/_pachucki_extension.py (jets in the four exponents, direct solver, field Q or Q(t)); the reviewed helper
tests/_pachucki_exchange.py is the reference it is tested against.  Nothing here reads or imports debug/.  Exact
outputs of the derivation's Q(t) runs are frozen in tests/data/paper58_pachucki_ext_2026_11_18.json; the fast tier
tests them, the slow tier recomputes them.

WHAT IS PROVED, AND WHAT IS NOT.
  * Exact: residuals of all derivative equations vanish identically (1s over Q(t), fast; l > 0 over Q at two
    pairs, fast or slow, and over Q(t), slow); K in the nine-function span; that span closed under d/dR.  Order
    <= 9 follows: ten derivatives in a nine-dimensional module are dependent.
  * The zero residual accepts ANY source: the solver solves whatever F it is handed.  Source correctness rests on
    (i) the jets being the derivatives of the reviewed transcription (test A1), (ii) the numerical tests.
  * Uniqueness by decay is used at every multi-index, so "for every t > 1" needs the denominators and the decay
    of every intermediate, not only of K (tests C).
  * NOT proved: order exactly 9.  That needs the nine functions independent over C(R).  The rank 9 computed here
    is a rank of coordinates in the formal span (test_formal_rank_nine_is_measured_in_the_formal_span_only).
  * NOT covered: m != 0 on the axis (2p+-1 1s exchange, 3d+-2), real 2p_x-2p_x exchange.

Each test's docstring names the wrong answer it rejects; fire-test plants are listed at the end of the file.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

import _pachucki_exchange as PX
import _pachucki_extension as X

DATA = json.loads((Path(__file__).resolve().parent / "data" / "paper58_pachucki_ext_2026_11_18.json").read_text())

F = Fraction
CASES = X.CASES
N_INDICES = {"1s": 16, "2psig": 144, "3dsig": 576, "ppsig": 576, "pppi0": 576}
FAST_CASES = ("1s", "2psig")                       # solved live in the fast tier (0.1 s, 3 s over Q)


def _split(cases=CASES, fast=FAST_CASES):
    return [c if c in fast else pytest.param(c, marks=pytest.mark.slow) for c in cases]


# The point the reviewer chose: a != 1 (every check in the derivation memo had a = 1), R not used elsewhere.
PT_A, PT_B, PT_R = F(4, 5), F(3, 2), F(23, 10)
# Reviewer's independent integrals at that point (Neumann expansion, dps 40), typed from the review memo.
ANCHOR = {
    "1s": "0.0955959302201830564664864485839835",
    "f0": "0.0022558102214482983799805712969021",
    "2psig": "0.171366050015632245442926777969856",
    "ppsig": "0.0427216693803001655774249080965235",
    "pppi0": "0.0797334869400471005133449208692647",
    "3dsig": "0.0850202908833205952575807791934367",
}
# K_1s at a = b = 1, R = 7/5 (derivation memo: limit = Sugiura = exchange_hp to 35 digits)
ANCHOR_EQUAL = "0.32329114155307318238571155893367021"

TOL_ANCHOR = mp.mpf(10) ** -30                    # the anchors are typed to 32-35 significant digits
NEU_FAST = dict(dps=9, guard=8)                   # fast-tier setting of the in-test integrator (2-3 s per integral)
TOL_NEU_FAST = mp.mpf(10) ** -9                   # measured: 1e-14 or better at that setting

ALLOWED_FACTORS = {"t", "t - 1", "t + 1"}
MAX_POLE = {"1s": 5, "2psig": 7, "3dsig": 9, "ppsig": 9, "pppi0": 9}      # order of the pole of K at t = 1
OP12_DEGREE = {"1s": 16, "2psig": 34, "ppsig": 53, "pppi0": 39, "3dsig": 53}
OP11_DEGREE = {"1s": 8, "2psig": 16, "ppsig": 24, "3dsig": 24, "pppi0": 16}


def committed_operator(case: str, pair: str):
    """The guessed operators, frozen over Q in the EX / EP / EQ guards (CRT + rational reconstruction of the
    records' mod-p operators).  pair '12' = (a, b) = (1, 2), '11' = (1, 1).  None where no record exists."""
    import test_paper58_exchange_d_pp as DPP
    import test_paper58_exchange_l1 as L1
    import test_paper58_exchange_unequal as U
    table = {
        ("1s", "12"): U.OP_EXCH_1_2, ("1s", "11"): U.OP_EXCH_1_EQUAL,
        ("2psig", "12"): L1.OP_S_1_2, ("2psig", "11"): L1.OP_S_1_1,
        ("ppsig", "12"): DPP.OP_PPSIG_1_2, ("ppsig", "11"): DPP.OP_PPSIG_1_1,
        ("pppi0", "12"): DPP.OP_PPPI0_1_2, ("pppi0", "11"): DPP.OP_PPPI0_1_1,
        ("3dsig", "11"): DPP.OP_3DSIG_1_1,
    }
    s = table.get((case, pair))
    return None if s is None else X.primitive(X.parse_operator(s))


# ---------------------------------------------------------------- shared exact objects
# 3dsig, ppsig and pppi0 need the SAME downset of 576 multi-indices (their densities have the same monomials in
# r_A^2, r_B^2), so they share the equations, the solutions f_alpha and every intermediate; only the final
# combination K differs.  The solve, the residual check and the audit are therefore cached per downset.
@lru_cache(maxsize=None)
def _solve(kind: str, a, b, Dkey: tuple):
    ctx = X.Ctx("Q", a, b) if kind == "Q" else X.Ctx("Qt")
    D = list(Dkey)
    calc, fS, Fs, P2, nfree = X.solve_all(ctx, D)
    return ctx, D, calc, fS, Fs, P2, nfree


@lru_cache(maxsize=None)
def solved_q(case: str, a: Fraction, b: Fraction):
    ctx, D, calc, fS, Fs, P2, nfree = _solve("Q", a, b, tuple(X.case_downset(case)))
    return ctx, D, calc, fS, Fs, P2, nfree, X.assemble_K(ctx, case, fS)


@lru_cache(maxsize=None)
def solved_qt(case: str):
    ctx, D, calc, fS, Fs, P2, nfree = _solve("Qt", None, None, tuple(X.case_downset(case)))
    return ctx, D, calc, fS, Fs, P2, nfree, X.assemble_K(ctx, case, fS)


@lru_cache(maxsize=None)
def nonzero_residuals(kind: str, a, b, Dkey: tuple):
    """The multi-indices of the downset whose equation (*) is NOT satisfied identically (expected: none)."""
    ctx, D, calc, fS, Fs, P2, nfree = _solve(kind, a, b, Dkey)
    return [al for al in D if sum(al) and residual(ctx, calc, al, fS, Fs, P2) != {}]


@lru_cache(maxsize=None)
def qt_audit(Dkey: tuple):
    ctx, D, calc, fS, Fs, P2, nfree = _solve("Qt", None, None, Dkey)
    return X.audit_intermediates(ctx, D, fS, Fs, P2)


@lru_cache(maxsize=None)
def px_solved(a: Fraction, b: Fraction):
    """The reviewed helper's own solution at (a, b): (alg, fS, Fs, P2, K)."""
    alg = PX.Algebra(a, b)
    fS, free, Fs, P2 = alg.solve_all()
    assert set(free.values()) == {0}
    return alg, fS, Fs, P2, alg.K(fS)


def residual(ctx, calc, al, fS, Fs, P2):
    """(r f_al')' + F_al - r sum_{0 != U <= al} binom(al, U) (d_U p^2) f_{al-U}, written here, not in the helper."""
    res = X.fadd(calc.L0(fS[al]), Fs[al])
    for U in itertools.product(*[range(x + 1) for x in al]):
        if sum(U) == 0:
            continue
        T = tuple(s - u for s, u in zip(al, U))
        bn = 1
        for s, u in zip(al, U):
            bn *= comb(s, u)
        res = X.fadd(res, X.fscale(fS[T], -ctx.c(bn) * P2[U], 1))
    return res


@lru_cache(maxsize=None)
def stored_K(case: str, tv: Fraction):
    """K(x; 1, tv) as a Func over Q, rebuilt from the STORED Q(t) coefficients on the nine functions."""
    ctx = X.Ctx("Q", 1, tv)
    K = {}
    for (g, _), c in zip(X.nine_functions(ctx), X.coeff_strings_at(DATA["sym"][case], tv)):
        for p, x in c.items():
            K = X.fadd(K, X.fscale(g, x, p))
    return ctx, K


def nine_coefficients_scaled(case: str, a: Fraction, b: Fraction):
    """Coefficients of R^p g_i(R; a, b) predicted from the stored Q(t) coefficients c_p(t) by
    K(R; a, b) = a K(aR; 1, b/a):  a^{p+1} c_p(b/a)  (f0: a^{p+3} c_p(b/a), since f0(R; a, b) = a^-2 f0(aR; 1, b/a))."""
    out = []
    for i, c in enumerate(X.coeff_strings_at(DATA["sym"][case], b / a)):
        out.append({p: a ** (p + (3 if i == 0 else 1)) * x for p, x in c.items()})
    return out


def unique_operator(calc, K, order: int, degree: int):
    """(the operator, number at degree - 1): operators of order <= `order` annihilating K in the function algebra."""
    ders = calc.derivatives(K, order)
    ops = PX.operators_annihilating(ders, order, degree)
    assert len(ops) == 1, f"{len(ops)} independent operators of order <= {order}, degree <= {degree}"
    return X.primitive(ops[0]), len(PX.operators_annihilating(ders, order, degree - 1))


def _mpq(v: Fraction):
    return mp.mpf(v.numerator) / v.denominator


@lru_cache(maxsize=None)
def f0_at(a: Fraction, b: Fraction, x: Fraction, dps: int = 60):
    """f0(x; a, b) by the label-59 quadrature of label 60 (the closed form's one quadrature)."""
    with mp.workdps(dps):
        return +X.f0_quadrature(a, b, _mpq(x))


# ================================================================ A. the sources
ALS_FAST = [(2, 0, 0, 0), (0, 2, 0, 0), (2, 1, 0, 0), (1, 2, 1, 0), (2, 0, 2, 0)]
ALS_SLOW = [(0, 0, 0, 2), (3, 0, 0, 0), (0, 0, 2, 0), (0, 3, 0, 0), (0, 2, 0, 2), (2, 1, 1, 1), (1, 1, 1, 2),
            (3, 1, 1, 1), (1, 3, 1, 1),
            (2, 2, 1, 1), (3, 1, 3, 1)]
JET_PAIRS = [(F(1), F(2)), (F(4, 5), F(3, 2))]
HELPER_PAIRS = [(F(1), F(2)), pytest.param(F(4, 5), F(3, 2), marks=pytest.mark.slow)]


def _jet_check(a, b, als):
    ctx = X.Ctx("Q", a, b)
    Fs, P2 = X.sources(ctx, X.downset(als))
    alg = PX.Algebra(a, b)
    sub = PX.exchange_point(sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator))
    for al in als:
        e, q = PX._dsym(al)                                    # plain sympy derivatives of labels 30-34 and 20
        assert Fs[al], al
        assert Fs[al] == alg.parse(e.subs(sub)), f"jet of F at {al} differs from the sympy derivative"
        assert P2[al] == PX._frac(sp.nsimplify(q.subs(sub))), f"jet of p^2 at {al}"


@pytest.mark.parametrize("a,b", JET_PAIRS)
def test_jet_sources_equal_sympy_derivatives_at_repeated_exponents(a, b):
    """A1.  d^alpha F and d^alpha p^2 from the jet engine equal the plain sympy derivatives of the reviewed
    helper's transcription (labels 30-34, 20), exactly, at multi-indices with a REPEATED exponent (order 2 in one
    or two variables, mixed with first orders; order 3 in the slow tier), at (1, 2) and (4/5, 3/2).

    Rejects: a jet rule right at order 1 and wrong at order 2 (the 1/k! of the exponential or Ei jet, the
    factorial that turns a Taylor coefficient into a derivative, the inverse or logarithm series truncated one
    term early).  The zero residual cannot see any of these: the solver solves whatever source it is given.
    """
    assert all(max(al) >= 2 for al in ALS_FAST)
    _jet_check(a, b, ALS_FAST)


@pytest.mark.slow
@pytest.mark.parametrize("a,b", JET_PAIRS)
def test_jet_sources_equal_sympy_derivatives_up_to_3131(a, b):
    """A1, the rest of the reviewer's list: eleven more multi-indices up to (3, 1, 3, 1), the top index of 2psig.

    Rejects: the same as the fast test, at total order up to 8 and with two exponents repeated.
    """
    _jet_check(a, b, ALS_SLOW)


@pytest.mark.parametrize("a,b", HELPER_PAIRS)
def test_engine_agrees_with_the_reviewed_helper_on_subsets(a, b):
    """A2.  On the 16 subsets (all orders <= 1) the engine reproduces the reviewed helper exactly: same sources,
    same label-60 function, same f_S and the same K; and d/dr of label 60 is the order-0 source.

    Rejects: drift between the copy of the engine in tests/_pachucki_extension.py and the reviewed transcription
    (a changed sign in FM1, a wrong a+-, a solver that returns a different element of the span).
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_q("1s", a, b)
    alg, fS0, Fs0, P20, K0 = px_solved(a, b)
    assert calc.FM1 == alg.FM1 and calc.diff(calc.FM1) == Fs[(0, 0, 0, 0)]
    assert len(D) == 16 and nfree == 0
    for al in D:
        assert Fs[al] == Fs0[al] and P2[al] == P20[al], al
        assert fS[al] == fS0[al], al
    assert K == K0 and K


# ================================================================ B. 1s over Q(t)
def test_1s_fifteen_residuals_vanish_identically_over_Qt():
    """B1.  With a = 1, b = t an indeterminate: each of the 15 equations (*) has residual identically 0 as an
    element of Q(t)[r, 1/r] x keys, no unknown of the solver is left free, and the p^2 coupling is not vacuous.

    Rejects: a solution valid only at sampled exponent pairs, a dropped or mis-signed coupling term, a system
    solved with a free direction.  Does NOT reject a wrong source (see A1, and the numerical tests).
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_qt("1s")
    assert len(D) == 16 and nfree == 0
    n = 0
    for al in D:
        if sum(al) == 0:
            continue
        assert fS[al]
        assert residual(ctx, calc, al, fS, Fs, P2) == {}, al
        if sum(al) >= 2:
            assert X.fadd(calc.L0(fS[al]), Fs[al]) != {}, al        # without the coupling the residual is not 0
        n += 1
    assert n == 15
    t = ctx.K.from_sympy(ctx.t)
    assert P2[(1, 1, 0, 0)] == 4 * t * (t * t + 1) / (1 - t * t) ** 2         # Hessian entry 4ab(a^2+b^2)/(a^2-b^2)^2


@pytest.mark.parametrize("tv", [F(2), pytest.param(F(3, 2), marks=pytest.mark.slow)])
def test_1s_Qt_closed_form_specialises_to_the_reviewed_helper(tv):
    """B2.  K over Q(t), specialised at t = 2 (and t = 3/2, slow tier), is the reviewed helper's K at (1, 2) and
    (1, 3/2), coefficient for coefficient (so it inherits that guard's agreement with the integral); and its
    coefficients on the nine functions are the frozen ones of tests/data, as rational functions of t.

    Rejects: a Q(t) closed form that solves its equations but is not the function the reviewed guard pins.
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_qt("1s")
    assert X.specialise(ctx, K, tv) == px_solved(F(1), tv)[4]
    stored = DATA["sym"]["1s"]
    for nm, c in zip(X.NAMES9, X.in_nine(ctx, K)):
        assert {str(p) for p in c} == set(stored[nm]), nm
        for p, x in c.items():
            assert x == ctx.K.from_sympy(sp.sympify(stored[nm][str(p)], locals={"t": ctx.t})), (nm, p)


def test_1s_K_in_nine_function_span_closed_under_differentiation_over_Qt():
    """B3.  Over Q(t): K is a Q(t)[R, 1/R]-combination of the nine functions (so the two ties, coordinate of
    Ei((a-b)R) = -2 x coordinate of Ei(2(a-b)R) and coordinate of LC e^{-2aR} = - coordinate of Ei(2(a-b)R), and
    their b-twins, hold as identities in t); each of the nine occurs; and d/dR maps each of the nine into the span.

    Rejects: ties that hold only at sampled pairs, a tenth function, a span that differentiation leaves (then
    "ten derivatives in a nine-dimensional module" would not follow).
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_qt("1s")
    co = X.in_nine(ctx, K)
    assert co is not None and all(co), "K not in the span, or one of the nine does not occur"
    assert len(K) == 13                                    # 11 functions + LC e^{-2aR}, LC e^{-2bR}
    a, b = ctx.a, ctx.b
    for x, y in ((a, b), (b, a)):
        lead = K[("Ei", -2 * x, 2 * (x - y))]
        assert K[("Ei", -2 * x, (x - y))] == {p: -2 * c for p, c in lead.items()}
        assert K[("E", -2 * x, 1)] == {p: -c for p, c in lead.items()}
    for g, _ in X.nine_functions(ctx):
        assert X.in_nine(ctx, calc.diff(g)) is not None
    # control: a function outside the span is reported as outside
    assert X.in_nine(ctx, {("E", -3 * b, 0): {0: ctx.one}}) is None
    assert X.in_nine(ctx, {("Ei", -2 * a, (a - b)): {0: ctx.one}}) is None       # half a bracket


def _raw_rank(ctx, K, tv, Rv):
    """Ranks of the first n rows, n = 1..11, of K, ..., K^(10) (K over Q(t) specialised at t = tv, then
    differentiated over Q) on the raw keys, coordinates evaluated at R = Rv."""
    ders = X.Calc(X.Ctx("Q", 1, tv)).derivatives(X.specialise(ctx, K, tv), 10)
    keys = PX.all_keys(ders)
    M = PX.coordinate_matrix_at(ders, keys, Rv)
    return len(keys), [PX.rank_q(M[:n]) for n in range(1, len(M) + 1)]


def test_1s_order_at_most_nine_ten_derivatives_are_dependent():
    """B4.  K, K', ..., K^(10) written on the THIRTEEN raw keys (no use of the nine-function reduction): at five
    (t, R), including the key-collision values t = 2, 3, the first ten rows already have rank <= 9, and so do
    all eleven.  With B3 (identically in t) this is the statement "order <= 9 for every t > 1".

    Rejects: an order bound of 10 or more (a missed tie would give rank 10 on the raw keys).  Says nothing about
    order >= 9: see the next test.
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_qt("1s")
    for tv, Rv in ((F(3), F(7, 5)), (F(2), F(-3, 11)), (F(15, 8), F(46, 25)), (F(11, 10), F(5, 3)), (F(40), F(1, 7))):
        nkeys, ranks = _raw_rank(ctx, K, tv, Rv)
        assert nkeys == 13
        assert ranks[9] <= 9 and ranks[10] <= 9, (tv, Rv, ranks)


def test_formal_rank_nine_is_measured_in_the_formal_span_only():
    """B5.  The rank of K, ..., K^(8) in the formal coordinates is 9 at the same five points, and at the collision
    pair (1, 3) run directly over Q (where e^{-(a+b)R} and e^{-2bR} Ei(2(b-a)R) share a rate).

    This is NOT a proof that the order is exactly 9, and the test does not claim it: the coordinates are formal,
    and order >= 9 for the function needs the nine functions to be linearly independent over C(R), which nothing
    here establishes.  "Exactly 9" stays measured (the guessed operators), as before.

    Rejects: a formal rank 8 (a tie the nine-function list missed, i.e. a span stated larger than needed), and
    rank 9 failing at a collision value of t.
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_qt("1s")
    for tv, Rv in ((F(3), F(7, 5)), (F(2), F(-3, 11)), (F(15, 8), F(46, 25)), (F(11, 10), F(5, 3)), (F(40), F(1, 7))):
        _, ranks = _raw_rank(ctx, K, tv, Rv)
        assert ranks[:9] == list(range(1, 10)) and ranks[9:] == [9, 9], (tv, Rv, ranks)
    ctx3, _, calc3, _, _, _, nfree3, K3 = solved_q("1s", F(1), F(3))
    ders3 = calc3.derivatives(K3, 10)
    assert nfree3 == 0 and all(X.in_nine(ctx3, A) is not None for A in ders3)
    M = PX.coordinate_matrix_at(ders3, PX.all_keys(ders3), F(7, 5))
    assert PX.rank_q(M[:9]) == 9 and PX.rank_q(M) == 9


def test_1s_K_coefficient_denominators_are_powers_of_b_minus_a_and_b_plus_a():
    """B6.  In (a, b), by K(R; a, b) = a K(aR; 1, b/a): every coefficient of K on the nine functions is a
    polynomial in a, b, R, 1/R over (b - a)^i (a + b)^j with i, j <= 5 (both reached), nothing else in the
    denominator; and the master integral enters with
        16 a^5 b^5 [9 + 6 (a^2 + b^2) R^2 + (a^4 + 4 a^2 b^2 + b^4) R^4] / (R (b^2 - a^2)^4).

    Rejects: a hidden pole at some b/a > 1 (a factor t - 3, 2t - 3, t^2 - 2), a surviving power of a or b in the
    denominator, a pole order 6, the normalisation 16 a^3 b^3 or the scaling a^{p+1} in place of a^{p+3} for f0.
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_qt("1s")
    co = X.in_nine(ctx, K)
    a, b, R, t = sp.symbols("a b R t", positive=True)
    mult = {}
    exprs = []
    for i, c in enumerate(co):
        e = sp.Integer(0)
        for p, x in c.items():
            ct = ctx.K.to_sympy(x)
            for f, m in sp.factor_list(sp.denom(sp.together(ct)))[1]:
                assert str(f) in ("t - 1", "t + 1"), (X.NAMES9[i], p, f)      # K itself: not even t
            e += a ** (p + (3 if i == 0 else 1)) * ct.subs(ctx.t, b / a) * R ** p
        e = sp.factor(sp.together(e))
        exprs.append(e)
        for f, m in sp.factor_list(sp.denom(e))[1]:
            mult[str(f)] = max(mult.get(str(f), 0), int(m))
    assert mult == {"R": 1, "a - b": 5, "a + b": 5}, mult
    f0_coefficient = 16 * a**5 * b**5 * (9 + 6 * (a**2 + b**2) * R**2 + (a**4 + 4 * a**2 * b**2 + b**4) * R**4) \
        / (R * (b**2 - a**2) ** 4)
    assert sp.simplify(exprs[0] - f0_coefficient) == 0
    assert sp.simplify(exprs[0] - f0_coefficient * a**2 / b**2) != 0


# ================================================================ C. every intermediate, every t > 1
def _check_audit(au: dict):
    for name, d in au["denominator_factors_max_multiplicity"].items():
        assert set(d) <= ALLOWED_FACTORS, (name, d)
    assert any(r["kind"] == "f" for r in au["keys"])
    assert sum(r["kind"] == "Ei" for r in au["keys"]) >= 7
    for r in au["keys"]:
        v = X.key_verdict(r)
        assert v["decays"], ("key does not decay for every t > 1", r, v)
        assert v["ei_zero"] == [], ("Ei argument vanishes at some t > 1", r, v)


def test_key_verdict_rejects_exceptional_values_of_t():
    """C0.  The decay rule itself: a key decays "for every t > 1" iff its rate lam (+ mu), linear in t, is < 0 on
    the whole half-line; an Ei argument must not vanish there.

    Rejects: a rule that samples t (t = 2 passes e^{(t-3)R}, which grows for t > 3), accepts rate 0, or misses
    Ei((t-2)R), undefined at t = 2.
    """
    ok = X.key_verdict({"kind": "Ei", "lam": "-2", "mu": "1 - t"})
    assert ok == {"rate": "-t - 1", "decays": True, "ei_zero": []}
    assert X.key_verdict({"kind": "Ei", "lam": "-2*t", "mu": "2*t - 2"})["decays"]            # rate -2, constant
    assert not X.key_verdict({"kind": "E", "lam": "t - 3", "log_power": 0})["decays"]         # grows for t > 3
    assert X.key_verdict({"kind": "E", "lam": "3 - 3*t", "log_power": 0})["decays"]           # rate 0 only at t = 1
    assert not X.key_verdict({"kind": "E", "lam": "0", "log_power": 1})["decays"]             # a constant
    assert not X.key_verdict({"kind": "Ei", "lam": "2*t", "mu": "-t - 1"})["decays"]          # rate t - 1 > 0
    assert X.key_verdict({"kind": "Ei", "lam": "-2*t", "mu": "t - 2"})["ei_zero"] == ["2"]
    assert X.key_verdict({"kind": "Ei", "lam": "-2", "mu": "0"})["ei_zero"] == ["identically"]
    assert X.key_verdict({"kind": "Ei", "lam": "-2", "mu": "2 - 2*t"})["ei_zero"] == []       # zero at t = 1 only


def test_1s_every_intermediate_has_denominators_t_tm1_tp1_and_decays_for_all_t():
    """C1 (live, over Q(t)).  For every one of the 16 multi-indices: the coefficient denominators of f_alpha,
    F_alpha and d_alpha p^2 factor into t, t - 1, t + 1 only, and every function key decays for every t > 1 with
    no vanishing Ei argument.  This is what turns the identity over Q(t) into a statement at every real t > 1
    (uniqueness by decay is used at every multi-index).  The frozen record of the run is reproduced exactly.

    Rejects: a hidden exceptional ratio in an intermediate that cancels in K (a pole at t = 3, a key e^{(b-3a)R}).
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_qt("1s")
    au = X.audit_intermediates(ctx, D, fS, Fs, P2)
    _check_audit(au)
    assert au["denominator_factors_max_multiplicity"]["f_alpha"] == {"t": 2, "t - 1": 5, "t + 1": 5}
    assert au == {k: DATA["denoms"]["1s"][k] for k in au}


@pytest.mark.parametrize("case", CASES)
def test_stored_audit_of_every_intermediate(case):
    """C2 (frozen record of the Q(t) runs; 2psig 144 multi-indices; 3dsig, ppsig, pppi0 share one downset of 576,
    so their records coincide; about 40 min for that solve).
    Same statement as C1 for the l > 0 cases, read from tests/data; the keys are re-judged here by the rule of
    C0, not read as stored verdicts.  The slow tier recomputes the record (C3).

    Rejects: a record containing a factor other than t, t - 1, t + 1 or a key that fails for some t > 1.  It does
    not by itself reject a record that is wrong: that is C3's job.
    """
    rec = DATA["denoms"][case]
    assert rec["multi_indices"] == N_INDICES[case] == len(X.case_downset(case))
    _check_audit(rec)
    assert set(rec["denominator_factors_max_multiplicity"]) == {"f_alpha", "F_alpha", "d_alpha_p2"}
    assert all(rec["denominator_factors_max_multiplicity"][k] for k in ("f_alpha", "F_alpha", "d_alpha_p2"))


@pytest.mark.parametrize("case", _split())
def test_live_keys_at_a_specialisation_are_the_stored_keys(case):
    """C2'.  The stored key list at t = 15/8, rates scaled by a = 4/5, is exactly the set of function keys that
    occur in the f_alpha and F_alpha of the live run over Q at (a, b) = (4/5, 3/2).

    Rejects: a stored key list that omits a key of an intermediate (the audit would then pass by omission).
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_q(case, PT_A, PT_B)
    live = {k for al in D for A in (fS[al], Fs[al]) for k in A}
    t = sp.Symbol("t")

    def val(s):
        v = sp.Rational(sp.sympify(s, locals={"t": t}).subs(t, sp.Rational(15, 8)))
        return PT_A * F(int(v.p), int(v.q))

    stored = set()
    for r in DATA["denoms"][case]["keys"]:
        if r["kind"] == "f":
            stored.add(("f",))
        elif r["kind"] == "E":
            stored.add(("E", val(r["lam"]), r["log_power"]))
        else:
            stored.add(("Ei", val(r["lam"]), val(r["mu"])))
    assert live == stored and len(stored) == len(DATA["denoms"][case]["keys"])


@pytest.mark.slow
@pytest.mark.parametrize("case", ["2psig", "3dsig", "ppsig", "pppi0"])
def test_Qt_run_recomputed_residuals_audit_and_closed_form(case):
    """C3 (slow: 2psig about 4.5 min; the three 576-index cases share one downset, hence one solve and one audit:
    about 65 min for the first of them, seconds for the other two).  The Q(t) run from scratch: every residual
    identically 0 with the multinomial binomials, no free unknown; the audit of every intermediate equals the
    frozen record and passes C1's conditions; K is in the nine-function span and its coefficients are the
    frozen ones.

    Rejects: frozen data that the engine does not reproduce (so every fast test on the frozen data is void), an
    l > 0 solution valid only at sampled pairs, an exceptional t in an intermediate.
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_qt(case)
    assert len(D) == N_INDICES[case] and nfree == 0
    assert nonzero_residuals("Qt", None, None, tuple(D)) == []
    au = qt_audit(tuple(D))
    _check_audit(au)
    assert au == {k: DATA["denoms"][case][k] for k in au}
    co = X.in_nine(ctx, K)
    assert co is not None
    stored = DATA["sym"][case]
    t = ctx.t
    for nm, c in zip(X.NAMES9, co):
        assert {str(p) for p in c} == set(stored[nm]), nm
        for p, x in c.items():
            assert x == ctx.K.from_sympy(sp.sympify(stored[nm][str(p)], locals={"t": t})), (nm, p)


@pytest.mark.parametrize("case", CASES)
def test_stored_closed_form_poles_only_at_t_0_1_minus1(case):
    """C4.  The frozen coefficients of K on the nine functions: denominators t - 1, t + 1 (and t for l > 0) only;
    the pole at t = 1 has order 5 (1s), 7 (2psig), 9 (3dsig, ppsig, pppi0).

    Rejects: a closed form with a pole at some b/a > 1; a pole order other than the one the limit test expands to.
    """
    t = sp.Symbol("t")
    mult = {}
    for nm in X.NAMES9:
        assert DATA["sym"][case][nm], nm
        for p, s in DATA["sym"][case][nm].items():
            for f, m in sp.factor_list(sp.denom(sp.together(sp.sympify(s, locals={"t": t}))))[1]:
                mult[str(f)] = max(mult.get(str(f), 0), int(m))
    assert set(mult) == ({"t - 1", "t + 1"} if case == "1s" else ALLOWED_FACTORS), mult
    assert mult["t - 1"] == MAX_POLE[case]


# ================================================================ D. l > 0 over Q, and the scaling in a
@pytest.mark.parametrize("case", _split())
def test_derivative_equations_solve_exactly_at_the_reviewer_point(case):
    """D1.  Over Q at (a, b) = (4/5, 3/2): every equation (*) of the downset (16, 144 or 576 multi-indices) has
    residual identically 0, with the multinomial binomials binom(alpha, U); no free unknown; K is in the
    nine-function span and so are its first ten derivatives.

    Rejects: a dropped binomial (invisible for 1s, where every binomial is 1), a solution with a free direction,
    an l > 0 integral outside the nine-function span.
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_q(case, PT_A, PT_B)
    assert len(D) == N_INDICES[case] and nfree == 0
    assert nonzero_residuals("Q", PT_A, PT_B, tuple(D)) == []
    nontrivial = sum(any(comb(s, u) > 1 for U in itertools.product(*[range(x + 1) for x in al])
                         for s, u in zip(al, U) if sum(U) and P2[U] != 0) for al in D)
    assert (nontrivial > 0) == (case != "1s")
    assert all(X.in_nine(ctx, A) is not None for A in calc.derivatives(K, 10))


@pytest.mark.parametrize("case", _split())
def test_scaling_law_and_stored_coefficients_exact(case):
    """D2.  The live closed form over Q at (a, b) = (4/5, 3/2) has, on the nine functions, exactly the
    coefficients a^{p+1} c_p(b/a) (f0: a^{p+3} c_p(b/a)) with c_p the frozen Q(t) coefficients: the scaling
    K(R; a, b) = a K(aR; 1, b/a), f0(R; a, b) = a^-2 f0(aR; 1, b/a), and the frozen data, at a != 1.

    Rejects: a wrong power of a in a normalisation (a^5 b^3 typed a^3 b^5 or a^4 b^4: invisible at a = 1), a
    scaling a^2 K or K, frozen coefficients that are not the engine's.
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_q(case, PT_A, PT_B)
    co = X.in_nine(ctx, K)
    assert co is not None and all(co)
    assert co == nine_coefficients_scaled(case, PT_A, PT_B)
    bad = [{p: PT_A * x for p, x in c.items()} for c in co]
    assert bad != nine_coefficients_scaled(case, PT_A, PT_B)


# ================================================================ E. the closed form is the integral
def closed_form_live(case: str, dps: int):
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_q(case, PT_A, PT_B)
    with mp.workdps(dps):
        x = _mpq(PT_R)
        return +X.evaluate(PT_A, PT_B, K, x, f0_at(PT_A, PT_B, PT_R))


def closed_form_stored(case: str, dps: int):
    """a K(aR; 1, b/a) from the frozen Q(t) coefficients."""
    tv = PT_B / PT_A
    ctx, K = stored_K(case, tv)
    with mp.workdps(dps):
        x = _mpq(PT_A * PT_R)
        return +(_mpq(PT_A) * X.evaluate(F(1), tv, K, x, f0_at(F(1), tv, PT_A * PT_R)))


def test_master_integral_equals_independent_integral():
    """E0.  f0(23/10; 4/5, 3/2): the label-59 quadrature of label 60 against the master integral (label 01)
    integrated directly in prolate spheroidal coordinates inside this test, and against the reviewer's digits.

    Rejects: the sign of the e^{2bR} E1 - e^{2aR} E1 pair of label 60 flipped (11% off; the one bug met in the
    derivation), a wrong prefactor 1/(2(a^2 - b^2)), a wrong LC.
    """
    with mp.workdps(45):
        f0 = f0_at(PT_A, PT_B, PT_R)
        ref = mp.mpf(ANCHOR["f0"])
        assert abs(f0 / ref - 1) < TOL_ANCHOR, mp.nstr(f0, 40)
        direct, per = X.neumann_integral("f0", PT_A, PT_B, PT_R, **NEU_FAST)
        assert len(per) >= 5 and abs(direct / f0 - 1) < TOL_NEU_FAST, (mp.nstr(direct, 20), mp.nstr(f0, 20))


@pytest.mark.parametrize("case", FAST_CASES)
def test_closed_form_equals_independent_integral_fast(case):
    """E1.  K(23/10; 4/5, 3/2) for 1s and 2psig, three routes: (i) the live closed form (engine over Q at the
    point, f0 by one quadrature), (ii) the integral computed inside this test by a Neumann expansion in prolate
    spheroidal coordinates, which never sees the closed form (to 1e-9 here; 30 digits in the slow tier),
    (iii) the reviewer's digits (to 1e-30).

    Rejects: a wrong power of a (a != 1 here; a^2 or a^0 for a is 20-25% off), and for 2psig a multi-index
    shifted by one, a sign (-1)^{|n|}, a wrong normalisation a^5 b^3.
    """
    with mp.workdps(45):
        v = closed_form_live(case, 45)
        ref = mp.mpf(ANCHOR[case])
        assert abs(v / ref - 1) < TOL_ANCHOR, (case, mp.nstr(v, 40))
        direct, per = X.neumann_integral(case, PT_A, PT_B, PT_R, **NEU_FAST)
        assert len(per) >= 5 and abs(direct / v - 1) < TOL_NEU_FAST, (mp.nstr(direct, 20), mp.nstr(v, 20))


@pytest.mark.parametrize("case", CASES)
def test_stored_closed_form_equals_the_reviewer_integrals(case):
    """E2.  All five cases from the FROZEN Q(t) coefficients through the scaling law a K(aR; 1, b/a), against the
    reviewer's independent integrals at (4/5, 3/2, 23/10), to 1e-30; and, for the three cases whose live closed
    form is slow-tier only (3dsig, ppsig, pppi0), against the integral computed inside this test, to 1e-9.

    Rejects (for the frozen data; the live engine is tested by E1 and, for the three 576-index cases, E3): the
    1/4 of the complex 2p+1 pair read as 1/2, the 1/18 of 3d0, a shifted multi-index, a wrong power of a in the
    scaling.
    """
    with mp.workdps(45):
        v = closed_form_stored(case, 45)
        assert abs(v / mp.mpf(ANCHOR[case]) - 1) < TOL_ANCHOR, (case, mp.nstr(v, 40))
        if case not in FAST_CASES:
            direct, per = X.neumann_integral(case, PT_A, PT_B, PT_R, **NEU_FAST)
            assert len(per) >= 5 and abs(direct / v - 1) < TOL_NEU_FAST, (mp.nstr(direct, 20), mp.nstr(v, 20))


@pytest.mark.slow
@pytest.mark.parametrize("case", CASES)
def test_closed_form_equals_independent_integral_to_30_digits(case):
    """E3 (slow).  All five cases live, three routes as in E1, the in-test integral at 30 digits.

    Rejects: the same as E1/E2, for the live engine, at 1e-30: a shifted alpha, a sign (-1)^{|n|}, 1/4 read as
    1/2, a wrong power of a.
    """
    with mp.workdps(50):
        v = closed_form_live(case, 50)
        assert abs(v / mp.mpf(ANCHOR[case]) - 1) < TOL_ANCHOR, (case, mp.nstr(v, 40))
        direct, per = X.neumann_integral(case, PT_A, PT_B, PT_R, dps=32)
        assert abs(direct / v - 1) < mp.mpf(10) ** -30, (mp.nstr(direct, 40), mp.nstr(v, 40))
        assert abs(direct / mp.mpf(ANCHOR[case]) - 1) < mp.mpf(10) ** -30


@pytest.mark.slow
def test_master_integral_equals_independent_integral_to_30_digits():
    """E0, slow: the in-test direct integral of label 01 at 30 digits.  Rejects: as E0, at 1e-30."""
    with mp.workdps(50):
        f0 = f0_at(PT_A, PT_B, PT_R)
        direct, _ = X.neumann_integral("f0", PT_A, PT_B, PT_R, dps=32)
        assert abs(direct / f0 - 1) < mp.mpf(10) ** -30
        assert abs(direct / mp.mpf(ANCHOR["f0"]) - 1) < mp.mpf(10) ** -30


# ================================================================ F. the limit b -> a
def sugiura():
    """Sugiura's equal-exponent 1s exchange integral at a = b = 1, in the symbols of X.LIM_SYMBOLS."""
    R, e2, Lg, X2, X4 = X.LIM_SYMBOLS[:5]
    S = 1 + R + R**2 / 3
    Sp = 1 - R + R**2 / 3
    return sp.Rational(1, 5) * (-e2 * (-sp.Rational(25, 8) + sp.Rational(23, 4) * R + 3 * R**2 + R**3 / 3)
                                + 6 / R * (e2 * S**2 * Lg - Sp**2 * X4 + 2 * S * Sp * X2))


def stored_limit(case: str):
    return sp.sympify(DATA["limit"][case], locals={str(s): s for s in X.LIM_SYMBOLS})


def limit_value(expr, Rv, dps: int):
    R, e2, Lg, X2, X4 = X.LIM_SYMBOLS[:5]
    with mp.workdps(dps):
        x = mp.mpf(Rv)
        vals = {R: x, e2: mp.exp(-2 * x), Lg: mp.euler + mp.log(x), X2: mp.e1(2 * x), X4: mp.exp(2 * x) * mp.e1(4 * x)}
        tot = mp.mpf(0)
        for c, lead in zip(X.four_coords(expr), (e2 * Lg, e2, X2, X4)):
            cv = sp.lambdify(R, c, "mpmath")(x)
            tot += cv * (vals[e2] * vals[Lg] if lead == e2 * Lg else vals[lead])
        return +tot


def test_stored_limit_is_sugiura_and_lies_in_the_four_function_span():
    """F1.  The frozen eps^0 term of K_1s(R; 1, 1 + eps) minus Sugiura's formula is 0 symbolically; it lies in
    the span of e^{-2R}(gamma + ln R), e^{-2R}, E1(2R), e^{2R} E1(4R) over Q[R, 1/R] (none of the three outside
    integrals J1, J2, J3); that span is closed under d/dR, so five derivatives are dependent (order <= 4); the
    one operator of order <= 4 and degree <= 8 is the committed equal-exponent operator over Q.

    Rejects: a limit that is not Sugiura's function (a wrong sign or factor in any of the four coordinates), a
    limit containing an integral outside the class, an order bound other than 4 at a = b.
    """
    lim = stored_limit("1s")
    R, e2, Lg, X2, X4, J1, J2, J3 = X.LIM_SYMBOLS
    assert sp.simplify(lim - sp.expand(sugiura())) == 0
    assert sp.simplify(lim - sp.expand(sugiura()) - e2 / R) != 0
    c4 = X.four_coords(lim)
    assert all(c != 0 for c in c4) and not lim.has(J1, J2, J3)
    with pytest.raises(ValueError):
        X.four_coords(lim + J3)
    # d/dR in coordinates is the derivative of the functions (checked on the limit itself, at a rational point)
    rep = {e2: sp.exp(-2 * R), Lg: sp.EulerGamma + sp.log(R), X2: sp.expint(1, 2 * R),
           X4: sp.exp(2 * R) * sp.expint(1, 4 * R)}
    d1 = sum(c * g for c, g in zip(X.four_dcoords(c4), (e2 * Lg, e2, X2, X4)))
    assert sp.simplify(sp.diff(lim.subs(rep), R) - d1.subs(rep)) == 0
    ops = X.four_operators(lim, 4, OP11_DEGREE["1s"])
    assert len(ops) == 1 and X.primitive(ops[0]) == committed_operator("1s", "11")
    assert X.four_operators(lim, 3, 30) == [] and X.four_operators(lim, 4, OP11_DEGREE["1s"] - 1) == []


def test_limit_formula_equals_independent_integral_at_equal_exponents():
    """F2.  Sugiura's formula (= the frozen limit) at a = b = 1, R = 7/5, against the equal-exponent exchange
    integral computed inside this test (Neumann expansion), and against the derivation memo's digits.

    Rejects: a limit that is a function of the right class but not the integral (e.g. the Coulomb integral, a
    factor 1/5 dropped), which the symbolic comparison with a mistyped Sugiura formula would let through.
    """
    with mp.workdps(40):
        v = limit_value(stored_limit("1s"), "1.4", 40)
        assert abs(v / mp.mpf(ANCHOR_EQUAL) - 1) < mp.mpf(10) ** -33, mp.nstr(v, 38)
        direct, per = X.neumann_integral("1s", 1, 1, F(7, 5), **NEU_FAST)
        assert len(per) >= 5 and abs(direct / v - 1) < TOL_NEU_FAST, (mp.nstr(direct, 20), mp.nstr(v, 20))


def test_closed_form_tends_to_the_limit_as_b_tends_to_a():
    """F3.  The live Q(t) closed form of K_1s, evaluated at t = 1 + eps for eps = 1/200, 1/400 (dps 60; the
    individual terms are of size eps^-5), approaches Sugiura's value at R = 7/5 linearly in eps: the difference
    halves with eps (to 5%) and the last one is below 1% of the value.

    Rejects: poles in b - a that do not cancel (the differences would grow like a negative power of eps; this
    is what the label-60 sign flip produced in the derivation), a limit that is not the eps^0 term.
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_qt("1s")
    lim = limit_value(sugiura(), "1.4", 60)
    with mp.workdps(60):
        x = mp.mpf("1.4")
        d = []
        for n in (200, 400):
            tv = 1 + F(1, n)
            Kq = X.specialise(ctx, K, tv)
            d.append(X.evaluate(F(1), tv, Kq, x, X.f0_quadrature(F(1), tv, x)) - lim)
        assert 0 < abs(d[1]) < abs(lim) / 100, [mp.nstr(v, 8) for v in d]
        assert abs(d[0] / d[1] - 2) < mp.mpf("0.05"), [mp.nstr(v, 8) for v in d]


@pytest.mark.parametrize("case", ["2psig", "3dsig", "ppsig", "pppi0"])
def test_stored_limits_for_l_gt_0_reproduce_the_equal_exponent_operators(case):
    """F4.  The frozen b -> a limits of the four l > 0 closed forms lie in the four-function span, and the one
    operator of order <= 4 at the recorded degree (2psig 16, 3dsig 24, ppsig 24, pppi0 16) that annihilates each
    is the committed equal-exponent operator of records EP / EQ over Q; none exists one degree lower.

    Rejects: a limit outside the equal-exponent class, a limit of the wrong integral (each operator is
    case-specific: the 2psig limit is not annihilated by the pppi0 operator of the same shape).
    """
    lim = stored_limit(case)
    deg = OP11_DEGREE[case]
    ops = X.four_operators(lim, 4, deg)
    assert len(ops) == 1
    assert X.primitive(ops[0]) == committed_operator(case, "11")
    assert X.four_operators(lim, 4, deg - 1) == []
    other = {"2psig": "pppi0", "pppi0": "2psig", "3dsig": "ppsig", "ppsig": "3dsig"}[case]
    assert X.primitive(ops[0]) != committed_operator(other, "11")


def _limit_recomputed(case: str, cft):
    out = X.limit_b_to_a(cft)
    assert out["max_pole"] == MAX_POLE[case]
    assert sorted(out["poles"]) == list(range(1, MAX_POLE[case] + 2))
    assert all(v == 0 for v in out["poles"].values()), {m: v for m, v in out["poles"].items() if v != 0}
    assert all(j == 0 for row in out["J_in_Ik"] for j in row)       # every eps-coefficient of f0 is in the class
    assert len(out["J_in_Ik"]) == MAX_POLE[case] + 1
    X.four_coords(out["limit"])                                     # raises if a J survives
    assert sp.expand(out["limit"] - stored_limit(case)) == 0
    return out["limit"]


def _live_1s_coefficient_strings():
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_qt("1s")
    return {nm: {str(p): str(sp.factor(ctx.K.to_sympy(x))) for p, x in c.items()}
            for nm, c in zip(X.NAMES9, X.in_nine(ctx, K))}


@pytest.mark.slow
def test_1s_limit_recomputed_poles_cancel():
    """F5 (slow, 10 s).  From the live Q(t) closed form: the expansion at t = 1 + eps (elementary pieces by Taylor,
    the psi brackets by the series of Ei, f0 term by term through label 60, each integral reduced by parts with rules
    checked by differentiation) has every pole eps^-1 .. eps^-6 equal to 0, no f0 coefficient outside the
    four-function class, and its eps^0 term is the frozen limit, hence Sugiura's formula.

    Rejects: uncancelled poles (the label-60 sign flip), a frozen limit that is not the limit of the closed form.
    """
    lim = _limit_recomputed("1s", _live_1s_coefficient_strings())
    assert sp.simplify(lim - sp.expand(sugiura())) == 0


@pytest.mark.slow
@pytest.mark.parametrize("case", ["2psig", "3dsig", "ppsig", "pppi0"])
def test_l_gt_0_limits_recomputed_poles_cancel(case):
    """F5 for l > 0 (slow), from the frozen Q(t) coefficients (themselves recomputed by C3): all poles up to
    order 8 (2psig) or 10 cancel, nothing outside the class, eps^0 term = the frozen limit.

    Rejects: as F5, per case.
    """
    _limit_recomputed(case, DATA["sym"][case])


# ================================================================ G. the operators at (1, 2)
OP12_SHAPE = {c: [d - 2, d - 1] + [d] * 8 for c, d in OP12_DEGREE.items()}


@pytest.mark.parametrize("case", _split(["1s", "2psig", "ppsig", "pppi0"]))
def test_operator_of_stored_closed_form_is_the_committed_one_at_1_2(case):
    """G1.  From the frozen Q(t) coefficients at t = 2: exactly one operator of order <= 9 and degree <= 16 (1s),
    34 (2psig), 53 (ppsig), 39 (pppi0) annihilates the closed form, none one degree lower, and it equals the guessed
    operator of records EX / EP / EQ over Q, every coefficient (both made primitive).  ppsig (8 s) and pppi0 (5 s)
    are in the slow tier.

    Rejects: a closed form of the wrong integral (the guessed operators come from exact Neumann-series data and
    never saw the closed form), an operator that only shares order and degree with the committed one.
    """
    ctx, K = stored_K(case, F(2))
    calc = X.Calc(ctx)
    op, n_lower = unique_operator(calc, K, 9, OP12_DEGREE[case])
    assert n_lower == 0
    assert [len(c) - 1 for c in op] == OP12_SHAPE[case]
    ref = committed_operator(case, "12")
    assert op == ref
    bad = [list(c) for c in ref]
    bad[3][5] += 1
    assert op != bad


@pytest.mark.slow
def test_3dsig_operator_at_1_2_is_unique_of_degree_53():
    """G2 (slow, 8 s).  3dsig has no record at (1, 2).  From the frozen coefficients at t = 2: exactly one operator
    of order <= 9 and degree <= 53, none of degree <= 52, coefficient degrees 51, 52, 53 (x 8); it is not the ppsig
    operator, which has the same shape.

    Rejects: a degree other than 53, a non-unique operator, the 3dsig closed form being the ppsig one.
    """
    ctx, K = stored_K("3dsig", F(2))
    op, n_lower = unique_operator(X.Calc(ctx), K, 9, 53)
    assert n_lower == 0
    assert [len(c) - 1 for c in op] == OP12_SHAPE["3dsig"]
    assert op != committed_operator("ppsig", "12")


@pytest.mark.parametrize("case", _split(fast=("1s",)))
def test_live_closed_form_at_1_2_is_the_stored_one(case):
    """G3.  The live engine over Q at (1, 2) gives, on the nine functions, exactly the frozen Q(t) coefficients at
    t = 2 (so G1/G2 are statements about the live closed form too); rank of K, ..., K^(10) in the formal
    coordinates is 9.

    Rejects: frozen coefficients that differ from the engine's at t = 2.
    """
    ctx, D, calc, fS, Fs, P2, nfree, K = solved_q(case, F(1), F(2))
    assert nfree == 0
    assert K == stored_K(case, F(2))[1]
    ders = calc.derivatives(K, 10)
    M = PX.coordinate_matrix_at(ders, PX.all_keys(ders), F(7, 5))
    assert PX.rank_q(M[:9]) == 9 and PX.rank_q(M) == 9


# ---------------------------------------------------------------------------------------------------------------
# TIMINGS (2026-11-18; machine shared with two other long jobs).
# Fast tier: 45 tests, 48-55 s (reviewed helper at (1, 2) 7.5 s; sympy derivatives for the jet test 6 s; the Q(t) solve
# of 1s 4.5 s; the in-test integrals 1-3 s each).
# Slow tier (--slow): 35 more tests.  Jets up to (3, 1, 3, 1): 7.6 min (sympy differentiation, 407 s + 49 s).
# Q(t) run of 2psig: 3.7 min.  Limits recomputed: 1s 8 s, 2psig 36 s, the three others 105 s each.  Live closed
# forms of the 576-index cases over Q and the 30-digit in-test integrals: 9-19 s each.  Q(t) run of the 576-index
# downset (shared by 3dsig, ppsig, pppi0): solve 32-36 min, residuals 25 min, i.e. about 65 min for the first of
# the three cases and seconds for the other two.  Whole slow tier: about 1 h 30 min.
# RUN RECORD.  pytest: fast tier 45 passed; slow tier 32 of 35 passed in two pytest runs (13 min 11 s + 9 min 8 s).
# The remaining three, test_Qt_run_recomputed_residuals_audit_and_closed_form[3dsig | ppsig | pppi0], were executed
# as function calls on in-process Q(t) solutions by the derivation-side driver that produced the frozen audit
# (one process per case, 60-65 min each, all passed; the last process also re-ran the other two bodies on its own
# solution, passed), not under the pytest runner.
#
# FIRE-TEST PLANTS (2026-11-18), written and fired as a separate step after the guard.  Every plant was applied
# to a TEMP COPY of tests/{conftest.py, _pachucki_exchange.py, _pachucki_extension.py, this file, the three operator
# guards, data/paper58_pachucki_ext_2026_11_18.json} and pytest run there; the tree was never mutated.  FIRED =
# pytest exit 1 (tests ran and failed); a control run of each selection without a plant passed first.
# H = tests/_pachucki_extension.py, T = this file, J = the data json, L1 = tests/test_paper58_exchange_l1.py.
# 82 runs: 77 FIRED; 5 did not fire and were expected not to (they show what a pin cannot see, and another pin
# covers each).  One of the 77 fired where the review memo predicted it could not.  No pin whose plant was expected
# to fire stayed green.
#
#   A. jets
#   Taylor coefficient not turned into a derivative, "fac *= factorial(x_)" -> "fac *= 1" (H): right at every
#     order <= 1, wrong from order 2                      -k jet_sources..._at_repeated               FIRED
#     same plant                                          -k engine_agrees (subsets only)             NOT FIRED (expected)
#     same plant                                          -k fifteen_residuals (1s, Q(t))             NOT FIRED (expected)
#     same plant                                          -k derivative_equations and 2psig           FIRED
#        (the review memo expected a zero residual for any source; in fact the solver's own consistency
#        conditions reject this one, "inconsistent: not in span".  Not relied upon.)
#     same plant                                          -k independent_integral_fast and 2psig      FIRED
#     same plant, slow list                               -k up_to_3131                               not run (7 min; same code path)
#   Ei jet: sign of the k-th derivative dropped (H)       -k jet_sources                              FIRED
#   log series: sign from the second term (H)             -k jet_sources                              FIRED
#   label 60 pair sign flipped in the engine (H)          -k engine_agrees                            FIRED
#   B. 1s over Q(t)
#   coupling sign in the solver (H)                       -k fifteen_residuals                        FIRED
#   coupling dropped in the solver (H)                    -k fifteen_residuals                        FIRED
#   1s normalisation b^3 -> b^4 (H)                       -k specialises_to_the_reviewed              FIRED
#   frozen 1s coefficient 144 -> 143 (J)                  -k specialises_to_the_reviewed              FIRED
#   psi tie -2 -> +2 in the nine functions (H)            -k nine_function_span_closed                FIRED
#   e^{-(a+b)R} removed from the nine (H)                 -k nine_function_span_closed                FIRED
#   label 60: -2 Ei -> -3 Ei (tie broken at the source, H) -k ten_derivatives / formal_rank_nine      FIRED x2 (in the solver)
#   a tenth function added to K (T)                       -k ten_derivatives_are_dependent            FIRED (rank 10)
#   K replaced by e^{-2bR}/R (T)                          -k formal_rank_nine                         FIRED
#   pole b = 3a planted in d_alpha p^2 (H)                -k coefficient_denominators / every_intermediate   FIRED x2 (in the solver)
#   factor t - 3 injected into one K coefficient (T)      -k coefficient_denominators                 FIRED
#   pole order 6 in one K coefficient (T)                 -k coefficient_denominators                 FIRED
#   1s normalisation b^3 -> b^4 (H)                       -k coefficient_denominators (pinned f0 term) FIRED
#   f0 scaling a^{p+3} -> a^{p+1} (T)                     -k coefficient_denominators                 FIRED
#   C. intermediates
#   decay rule samples t = 2 (H)                          -k key_verdict                              FIRED
#   Ei-argument zeros not reported (H)                    -k key_verdict                              FIRED
#   factor t - 3 injected into one f_alpha (T)            -k every_intermediate_has_denominators      FIRED
#   key e^{(t-3)R} injected into one f_alpha (T)          -k every_intermediate_has_denominators      FIRED
#   frozen audit 1s: factor t - 3 (J)                     -k stored_audit and 1s                      FIRED
#   frozen audit 1s: key rate t - 3 (J)                   -k stored_audit and 1s                      FIRED
#   frozen audit 1s: Ei argument t - 2 (J)                -k stored_audit and 1s                      FIRED
#   frozen audit pppi0: factor 2t - 3 (J)                 -k stored_audit and pppi0                   FIRED
#   frozen audit ppsig: key rate t - 5 (T)                -k stored_audit and ppsig                   FIRED
#   frozen audit 3dsig: Ei argument 2t - 5 (T)            -k stored_audit and 3dsig                   FIRED
#   frozen audit 1s: one key row removed (J)              -k live_keys and 1s                         FIRED
#   frozen 1s closed form: (t - 1)^4 -> (t - 3)^4 (J)     -k poles_only and 1s                        FIRED
#   frozen 3dsig closed form: pole at t = 3 (T)           -k poles_only and 3dsig                     FIRED
#   frozen 2psig coefficient x 2 (J), slow                -k Qt_run_recomputed and 2psig              FIRED
#   frozen 2psig audit multiplicity altered (J), slow     -k Qt_run_recomputed and 2psig              FIRED
#   D. l > 0 over Q, scaling
#   multinomial binomial dropped in the solver (H)        -k derivative_equations and 2psig           FIRED
#     same plant                                          -k derivative_equations and 1s              NOT FIRED (expected: all binomials 1)
#   1s normalisation a^3 -> a^2 (H)                       -k scaling_law and 1s                       FIRED
#     same plant                                          -k specialises_to_the_reviewed (a = 1)      NOT FIRED (expected: why a != 1)
#   2psig normalisation a^5 b^3 -> a^3 b^5 (H)            -k scaling_law and 2psig                    FIRED
#   scaling prediction a^{p+1} for f0 (T)                 -k scaling_law and 1s                       FIRED
#   pppi0 1/4 -> 1/2 (H), slow                            -k scaling_law and pppi0                    FIRED
#   E. numerics
#   label 60 pair sign flipped, engine (H)                -k test_master_integral                     FIRED (11% off)
#   label 60 pair sign flipped, f0 quadrature only (H)    -k test_master_integral / stored_closed_form and ppsig / tends_to_the_limit   FIRED x3
#   1s normalisation a^3 -> a^2 (H)                       -k independent_integral_fast and 1s         FIRED
#   2psig normalisation a^5 b^3 -> a^3 b^5 (H)            -k independent_integral_fast and 2psig      FIRED
#   2psig multi-index shifted by one (H)                  -k independent_integral_fast and 2psig      FIRED
#   per-term sign (-1)^{|alpha|/2} in the map (H)         -k independent_integral_fast and 2psig      FIRED
#   in-test integrator prefactor 2/R -> 1/R (H)           -k independent_integral_fast and 1s         FIRED
#   anchor K_1s, digit 20 altered (T)                     -k independent_integral_fast and 1s         FIRED
#   frozen pppi0 coefficient x 2 (J)                      -k stored_closed_form and pppi0             FIRED
#   scaling a K(aR) -> K(aR) for frozen forms (T)         -k stored_closed_form and 3dsig             FIRED
#   pppi0 1/4 -> 1/2 (H), slow                            -k to_30_digits and pppi0                   FIRED
#   3dsig density 3z^2 - r^2 -> 2z^2 - r^2 (H), slow      -k to_30_digits and 3dsig                   FIRED
#   F. limit
#   frozen limit 4/15 -> 4/17 (J)                         -k stored_limit_is_sugiura                  FIRED
#   frozen limit + J3/R (J)                               -k stored_limit_is_sugiura                  FIRED
#   Sugiura 23/4 -> 21/4 (T)                              -k stored_limit_is_sugiura                  FIRED
#   equal-exponent 1s operator -> the 2psig one (T)       -k stored_limit_is_sugiura                  FIRED
#   limit and Sugiura both x 2 (T)                        -k limit_formula_equals_independent         FIRED
#     same plant                                          -k stored_limit_is_sugiura                  NOT FIRED (expected: scale-free)
#   Sugiura with the E1(2R) term sign flipped (T)         -k limit_formula_equals_independent         FIRED
#   label 60 pair sign flipped, engine (H)                -k tends_to_the_limit                       FIRED (in the solver)
#   frozen 2psig limit altered (J)                        -k stored_limits_for_l_gt_0 and 2psig       FIRED
#   pppi0 (1,1) operator -> the 2psig one, same shape (T) -k stored_limits_for_l_gt_0 and pppi0       FIRED
#   label 60 pair sign flipped in the expansion (H), slow -k test_1s_limit_recomputed                 FIRED (poles survive)
#   frozen 2psig limit altered (J), slow                  -k l_gt_0_limits_recomputed and 2psig       FIRED
#   G. operators
#   committed 2psig (1,2) operator: 96 -> 97 in p_0 (L1)  -k operator_of_stored and 2psig             FIRED
#   table: 2psig (1,2) -> the pppi0 operator (T)          -k operator_of_stored and 2psig             FIRED
#   2psig degree 34 -> 33 (T)                             -k operator_of_stored and 2psig             FIRED (none)
#   1s degree 16 -> 17 (T)                                -k operator_of_stored and 1s                FIRED (two)
#   frozen pppi0 coefficient x 2 (J), slow                -k operator_of_stored and pppi0             FIRED
#   table: ppsig (1,2) -> the pppi0 operator (T), slow    -k operator_of_stored and ppsig             FIRED
#   frozen 3dsig coefficient x 2 (J), slow                -k 3dsig_operator                           FIRED
#   ppsig closed form passed off as 3dsig (T), slow       -k 3dsig_operator                           FIRED
#   1s normalisation b^3 -> b^4 (H)                       -k live_closed_form_at_1_2 and 1s           FIRED
#   2psig normalisation b^3 -> b^4 (H), slow              -k live_closed_form_at_1_2 and 2psig        FIRED
#
# NOT PINNED, and why.
#   * Jet sources against sympy at orders 4 and 5 in one exponent (needed by the 576-index downset, top index
#     (5, 1, 5, 1)): plain sympy differentiation at (3, 1, 3, 1) already costs 7 min.  Those orders are covered
#     only end to end (E2/E3 numerics, G operators).
#   * That the decaying solution found IS d^alpha f (decay of d^alpha f itself): an argument, not a test, as in
#     tests/test_paper58_pachucki_exchange.py.
#   * Order exactly 9 (independence of the nine functions over C(R)): not proved, not claimed.
#   * m != 0 on the axis (2p+-1 1s, 3d+-2) and real 2p_x-2p_x exchange: outside the derivation.
