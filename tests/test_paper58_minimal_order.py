"""Paper 58 -- minimal order in R of the two-centre integrals: what is PROVED, and nothing more.

Guard for the minimal-order results derived and cold-reviewed 2026-11-18 (this file is the guard only, written as
a separate pass, CLAUDE.md Sec. 9).  The proofs are by argument (Theorem A, Lemma GI); this file pins their inputs
and the exact computations they rest on.

    Theorem A.  f = sum_j r_j(R) e^{lam_j R}, lam_j distinct real, r_j in C(R) nonzero: minimal order exactly k
                (the number of distinct rates with a nonzero coefficient), at ANY exponents.
    Lemma GI.   The atoms e^{cR}, e^{cR} ln R, e^{cR} Ei(kappa R) (E1(kR) = -Ei(-kR)) with distinct labels
                (c, kappa), real rates, are linearly independent over C(R).
    Class A.    one-electron S, T, <a|1/r_A|b>, <a|1/r_B|a>; (AA|BB); hybrids with an s-type same-centre pair:
                order = number of distinct rates (Theorem A), every exponent choice.
    Class B.    hybrids with an L >= 1 same-centre pair: order PROVED AT THE COMPUTED EXPONENTS ONLY (GI + an exact
                minor + an exact generic rank).  Rates -a_d, -a_c, -(a_p + a_c), +a_d; the coincidences
                a_d = a_c AND a_d = a_p + a_c each drop the order by one.  At (2p0 2p0|1s 1s_B), Z = (3, 6) the closed
                form carries Ei and gamma + ln R together, three rates, order 4.
    Class C.    equal-exponent exchange, seven integrals (1s, 2psig, 2ppi, 3dsig, 3ddel, ppsig, pppi0): order
                exactly 4.  The d/dR-closed list e^{-2R}, e^{-2R} ln R, E1(2R), e^{2R} E1(4R) gives <= 4; a nonzero,
                gamma-free 4 x 4 determinant gives >= 4; every a = b by K(R; z) = z K(zR; 1).

HOW ORDERS ARE READ (tests/_minimal_order.py, own code path; nothing from debug/, nothing from the decompose /
d_atoms routines of tests/test_paper58_r_operators.py):
  * lower bound: an r x r minor of (f, ..., f^(r-1)) on the atoms at R = 37/59, exact as a polynomial in the
    constant symbols, evaluated at the TRUE constants in interval arithmetic; never "constants as free symbols";
  * upper bound: closure size, or the GENERIC rank over Q(R) of the split vectors of f, ..., f^(r); never a rank
    at a point;
  * logs of numbers reduced to logs of primes, and E1/Ei spellings merged onto one label, before anything is
    counted.

NOT ASSERTED (the proofs do not claim it): order of a class-B hybrid at exponents other than those listed
("rates + 1 at every ratio" is measured, not proved); any m != 0 unequal-exponent exchange integral (no closed form
on a finite d/dR-closed list is established); the unequal-exponent exchange (order 9: tests/test_paper58_exact_order
.py); the hand proofs of Theorem A and Lemma GI themselves.

FAST / SLOW.  Slow only: 3dsig of class C (its closed form takes ~8 s to build, ~11 s with the analysis), and the
independent prolate quadrature of the (3, 6) hybrid (~20 s).  The fast tier therefore covers six of the seven
class-C integrals; the 3d tau-cut plant is caught in the fast tier through 3ddel.

Each test's docstring names the wrong answer it rejects; fire-test plants are listed at the end of the file.
"""
from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

import _minimal_order as M
from geovac import qfd_core as Q
from geovac import two_center_one_electron as E
from geovac.two_center_eri import aabb_closed_form, hybrid_closed_form

R = M.R
F = Fraction
S1, P0, P1 = (1, 0, 0), (2, 1, 0), (2, 1, 1)
X0 = sp.Rational(37, 59)
DATA = json.loads((Path(__file__).parent / "data" / "paper58_minimal_order_2026_11_18.json").read_text())


def lab(c, kind, k=0):
    return (sp.Rational(c), kind, sp.Rational(k))


# ---------------------------------------------------------------- closed forms (live, from geovac)
def _o(s):
    return (int(s[0]), "spd".index(s[1]), int(s[2:]) if len(s) > 2 else 0)


def _qa(z):
    return ("A", F(z), 1)


def _qb(z):
    return ("B", F(z), 1)


ONE_CENTRE_1S = {
    "S": lambda a, b: Q.overlap(_qa(a), _qb(b), R),
    "T": lambda a, b: Q.kinetic(_qa(a), _qb(b), R),
    "VA(ab)": lambda a, b: Q._inv_r(_qa(a), _qb(b), "A", R),
    "VB(aa)": lambda a, b: Q._inv_r(_qa(a), _qa(a), "B", R),
    "(AA|BB)1s": lambda a, b: aabb_closed_form(F(a), S1, S1, F(b), S1, S1),
    "(AA|AB)1s": lambda a, b: hybrid_closed_form(F(a), S1, S1, S1, F(b), S1),
}

R8_A = {   # record R8 class-A cases (l <= 2), same definitions as the record's panel
    "E5VB": lambda: E.VB(("A", F(3)) + _o("3p1"), ("B", F(2)) + _o("3d1")),
    "E4T": lambda: E.T_lap(("A", F(3)) + _o("3d2"), ("B", F(3, 2)) + _o("3d2")),
    "F5": lambda: E.VB(("A", F(3)) + _o("3d2"), ("A", F(3)) + _o("3d2")),
    "C4": lambda: aabb_closed_form(F(3), _o("2p0"), _o("3d0"), F(2), _o("1s"), _o("2p0")),
    "H9": lambda: hybrid_closed_form(F(2), _o("2s"), _o("2s"), _o("3d0"), F(3), _o("3d0")),
}


@lru_cache(maxsize=None)
def fn_1s(name: str, a, b) -> M.Fn:
    return M.Fn(ONE_CENTRE_1S[name](a, b))


@lru_cache(maxsize=None)
def fn_r8(name: str) -> M.Fn:
    return M.Fn(R8_A[name]())


def hyb_expr(fam: str, za, zb) -> sp.Expr:
    if fam == "p0":      # (2p0 2p0 | 1s 1s_B): a_p = Z_A, a_c = Z_A, a_d = Z_B
        return hybrid_closed_form(F(za), P0, P0, S1, F(zb), S1)
    return hybrid_closed_form(F(za), P1, P1, P1, F(zb), P1)   # (2p1 2p1|2p1 2p1_B): a_p = Z_A, a_c = Z_A/2, a_d = Z_B/2


def hyb_rates(fam: str, za, zb):
    """(a_p, a_c, a_d) for the two families."""
    za, zb = sp.Rational(za), sp.Rational(zb)
    return (za, za, zb) if fam == "p0" else (za, za / 2, zb / 2)


@lru_cache(maxsize=None)
def fn_hyb(fam: str, za, zb) -> M.Fn:
    return M.Fn(hyb_expr(fam, za, zb))


XEQ = {
    "1s": (1, S1, S1, 1, S1, S1),
    "2psig": (2, (2, 1, 0), S1, 1, (2, 1, 0), S1),
    "2ppi": (2, (2, 1, 1), S1, 1, (2, 1, -1), S1),
    "3dsig": (3, (3, 2, 0), S1, 1, (3, 2, 0), S1),
    "3ddel": (3, (3, 2, 2), S1, 1, (3, 2, -2), S1),
    "ppsig": (2, (2, 1, 0), (2, 1, 0), 2, (2, 1, 0), (2, 1, 0)),
    "pppi0": (2, (2, 1, 1), (2, 1, 1), 2, (2, 1, 1), (2, 1, 1)),
}
TAU_MAX = 6


def exch_terms(case: str, scale=1, tau_max: int = TAU_MAX):
    ZA, oa, ob, ZB, oc, od = XEQ[case]
    return Q.exchange_closed_form(F(ZA) * scale, oa, ob, F(ZB) * scale, oc, od, R, tau_max=tau_max,
                                  return_terms=True)


@lru_cache(maxsize=None)
def exch_expr(case: str, scale=1) -> sp.Expr:
    return sp.Add(*exch_terms(case, scale))


@lru_cache(maxsize=None)
def fn_exch(case: str, scale=1) -> M.Fn:
    return M.Fn(exch_expr(case, scale))


def _xeq_param(case):
    return pytest.param(case, marks=pytest.mark.slow) if case == "3dsig" else case


XEQ_PARAMS = [_xeq_param(c) for c in XEQ]
CLASS_C_LIST = [lab(-2, "1"), lab(-2, "y", 0), lab(0, "y", -2), lab(2, "y", -4)]


# ================================================================ L. inputs of Theorem A / Lemma GI; the three traps
@pytest.mark.parametrize("which", ["exch1s", "hyb36", "hyb26", "aabb12", "VBaa32"])
def test_atom_calculus_is_d_dR(which):
    """Rejects: a wrong derivative rule for an atom (Ei' sign, the e^{(c+kappa)R}/R term, ln R' = 1/R), or a
    collector that loses or invents a term.  f, f', f'' from the helper's atom calculus against sympy.diff of the
    live closed form, numerically at R = 7/5 (40 digits)."""
    expr = {"exch1s": lambda: exch_expr("1s"), "hyb36": lambda: hyb_expr("p0", 3, 6),
            "hyb26": lambda: hyb_expr("p1", 2, 6), "aabb12": lambda: ONE_CENTRE_1S["(AA|BB)1s"](1, 2),
            "VBaa32": lambda: ONE_CENTRE_1S["VB(aa)"](F(3, 2), 1)}[which]()
    f = M.Fn(expr)
    x0 = sp.Rational(7, 5)
    for j in range(3):
        with mp.workdps(45):
            mine = M.numeric(f, j, x0, dps=45)
            ref = mp.mpf(str(sp.N(sp.diff(expr, R, j).subs(R, x0), 50)))
            assert abs(mine - ref) <= mp.mpf(10) ** -30 * max(1, abs(ref)), (which, j, mine, ref)


def test_redundant_spelling_is_one_atom():
    """Rejects: E1(kR) and Ei(-kR) treated as two atoms (review G3).  f = e^{-R}E1(2R) + e^{-R}Ei(-2R) + e^{-3R}
    IS e^{-3R}: order 1, and a claim of order 2 must be refused."""
    f3 = sp.exp(-R) * sp.expint(1, 2 * R) + sp.exp(-R) * sp.Ei(-2 * R) + sp.exp(-3 * R)
    f = M.Fn(f3)
    assert f.labels() == [lab(-3, "1")]
    assert M.order_is(f, 1)
    assert not M.order_at_least(f, 2)


@pytest.mark.parametrize("variant", ["R_exp", "E1"])
def test_dependent_constants_are_not_independent(variant):
    """Rejects: gamma / logs of numbers read as independent symbols (review G1, memo K2).  The coefficient
    ln 6 - ln 2 - ln 3 is ZERO: f has order 1, and order 2 must be refused -- by the prime-log reduction AND,
    independently, by the true-constant evaluation of the minor (checked with the reduction switched off)."""
    c = sp.log(6) - sp.log(2) - sp.log(3)
    f_expr = (c * R * sp.exp(-R) if variant == "R_exp" else c * sp.exp(-R) * sp.expint(1, R)) + sp.exp(-2 * R)
    f = M.Fn(f_expr)
    assert M.order_is(f, 1)
    assert not M.order_at_least(f, 2)
    g = M.Fn(f_expr, reduce_logs=False)           # logs kept unreduced: symbols k_i for ln 2, ln 3, ln 6
    assert len(g.syms) == 3
    assert not M.order_at_least(g, 2)             # the minor's interval at the true constants contains 0
    assert M.order_at_least(g, 1)


def test_upper_bound_is_generic_rank_not_point_rank():
    """Rejects: an order (upper bound) read from a rank at a single point (review G4).  The plant
    e^{-R} + e^{-2R} + (2R/3 + 103/177) e^{-3R} has true order 3 but its f, f', f'' have rank 2 at R = 37/59."""
    g = M.Fn(sp.exp(-R) + sp.exp(-2 * R) + (2 * R / 3 + sp.Rational(103, 177)) * sp.exp(-3 * R))
    assert M.point_rank(g, 3, X0) == 2            # the trap is real at 37/59 ...
    assert M.point_rank(g, 3, sp.Rational(11, 7)) == 3
    assert M.generic_rank(g, 3) == 3              # ... and the generic rank sees through it
    assert not M.order_at_most(g, 2)
    assert M.order_is(g, 3)


def test_vanishing_coefficient_is_not_a_rate():
    """Rejects: a rate whose coefficient vanishes in truth being counted (Theorem A counts rates with NONZERO
    coefficient).  S(1, 2) + (ln 6 - ln 2 - ln 3) R e^{-5R}: two rates, order 2, not 3 -- with the log reduction on,
    and with it off (then the true-constant evaluation of the coefficient drops the rate)."""
    base = ONE_CENTRE_1S["S"](1, 2)
    plant = base + (sp.log(6) - sp.log(2) - sp.log(3)) * R * sp.exp(-5 * R)
    for f in (M.Fn(plant), M.Fn(plant, reduce_logs=False)):
        assert len(f.rates()) == 2
        assert M.order_is(f, 2)
        assert not M.order_at_least(f, 3)


# ================================================================ A. Theorem A classes (exponential-polynomial)
_ORDER_1S = {   # (name) -> (order at a = b, order at a != b) = number of distinct rates (VB(aa): -2a and the 1/R tail 0)
    "S": (1, 2), "T": (1, 2), "VA(ab)": (1, 2), "VB(aa)": (2, 2), "(AA|BB)1s": (2, 3), "(AA|AB)1s": (2, 3),
}
CLASS_A_1S = [(nm, a, b, _ORDER_1S[nm][0 if a == b else 1])
              for nm in ONE_CENTRE_1S for (a, b) in ((1, 1), (1, 2), (2, 1), (F(3, 2), 1))]
CLASS_A_R8 = [("E5VB", 2), ("E4T", 2), ("F5", 2), ("C4", 3), ("H9", 3)]


def _check_class_a(f: M.Fn, order: int):
    labels = f.labels()
    assert all(k == "1" for (_c, k, _kap) in labels), labels      # exponential-polynomial: no ln / E1 / Ei atom
    assert len(f.rates()) == order                                   # rates with nonzero coefficient
    assert M.closure(labels) == labels                               # d/dR-closed: order <= #rates
    mn = M.lower_minor(f, order)                                     # n x n, true constants: order >= #rates
    assert mn.nonzero
    assert M.order_at_most(f, order)
    if order > 1:
        assert not M.order_at_most(f, order - 1)


@pytest.mark.parametrize("name,a,b,order", CLASS_A_1S, ids=[f"{n}-{a}-{b}" for n, a, b, _ in CLASS_A_1S])
def test_class_a_1s_order_is_number_of_rates(name, a, b, order):
    """Rejects: a claimed order one above or one below the number of distinct rates with nonzero coefficient
    ((AA|BB) 4, S(1,2) 1, ...); a planted ln / E1 atom passing as class A; a rate with a vanishing coefficient
    counted (the full n x n determinant at the true constants would vanish).  The 1/R tail of <a|1/r_B|a> counts
    as the rate 0."""
    _check_class_a(fn_1s(name, a, b), order)


@pytest.mark.parametrize("name,order", CLASS_A_R8, ids=[n for n, _ in CLASS_A_R8])
def test_class_a_record_r8_order_is_number_of_rates(name, order):
    """Same rejections, for five record-R8 cases with l <= 2 (two-centre one-electron, (AA|BB), s-type-pair hybrid)."""
    _check_class_a(fn_r8(name), order)


@pytest.mark.parametrize("name", ["S", "(AA|BB)1s", "(AA|AB)1s"])
def test_merging_two_rates_lowers_the_order_by_exactly_one(name):
    """Rejects: a merged pair of rates not lowering the order (Theorem A: coincident rates are counted once), and
    an 'exceptional ratio' near the coincidence: at b = a the order is one less than at b = 2a AND at b = 1001a/1000."""
    o_eq = len(fn_1s(name, 1, 1).rates())
    for b in (2, F(1001, 1000)):
        f = fn_1s(name, 1, b)
        assert M.order_is(f, o_eq + 1), (name, b)
    assert M.order_is(fn_1s(name, 1, 1), o_eq)
    assert not M.order_at_least(fn_1s(name, 1, 1), o_eq + 1)


# ================================================================ B. hybrids with an L >= 1 pair, computed exponents
CLASS_B = [   # (family, Z_A, Z_B, order, number of rates) -- each proved at THESE exponents only
    ("p0", 3, 2, 5, 4),     # a_d < a_c
    ("p0", 3, 3, 4, 3),     # a_d = a_c            (first coincidence)
    ("p0", 3, 4, 5, 4),     # a_c < a_d < a_p + a_c
    ("p0", 3, 5, 5, 4),
    ("p0", 3, 6, 4, 3),     # a_d = a_p + a_c      (second coincidence)
    ("p0", 3, 7, 5, 4),     # a_d > a_p + a_c
    ("p1", 2, 2, 4, 3),     # H2: a_d = a_c = 1
    ("p1", 2, 3, 5, 4),
    ("p1", 2, 6, 4, 3),     # a_d = 3 = a_p + a_c
]


@pytest.mark.parametrize("fam,za,zb,order,nrates", CLASS_B, ids=[f"{c[0]}-{c[1]}-{c[2]}" for c in CLASS_B])
def test_class_b_order_at_the_computed_exponents(fam, za, zb, order, nrates):
    """Rejects: order r - 1 or r + 1 accepted (order 5 at (3, 6): every 5 x 5 minor vanishes identically, generic
    rank 4); an upper bound read at a point (generic rank over Q(R) only); a minor that depends on gamma or a log
    (pinned gamma- and log-free, so no arithmetic property of gamma, ln 2, ... is used); 'order = rates + 1'
    failing at a computed point."""
    f = fn_hyb(fam, za, zb)
    assert len(f.rates()) == nrates
    assert order == nrates + 1
    mn = M.lower_minor(f, order)
    assert mn.nonzero, (fam, za, zb)
    assert mn.depends_on(f) == [], mn.depends_on(f)            # gamma-free and log-free
    assert len(M.closure(f.labels())) > order                  # so the upper bound needs the generic rank ...
    assert M.generic_rank(f, order + 1) == order               # ... over Q(R): an order-r operator exists
    assert M.generic_rank(f, order) == order
    assert not M.order_at_least(f, order + 1)


@pytest.mark.parametrize("fam,za,zb,order,nrates", CLASS_B, ids=[f"{c[0]}-{c[1]}-{c[2]}" for c in CLASS_B])
def test_class_b_function_content_names_both_coincidences(fam, za, zb, order, nrates):
    """Rejects: 'gamma + ln R iff a_d = a_c' (fails at (2p0 2p0|1s 1s_B) Z = (3, 6) and (2p1 2p1|2p1 2p1_B)
    Z = (2, 6), where a_d = a_p + a_c, a_d != a_c and the log comes WITH Ei); 'Ei iff a_d > a_c' broken; the order
    drop attributed to one coincidence only."""
    f = fn_hyb(fam, za, zb)
    a_p, a_c, a_d = hyb_rates(fam, za, zb)
    labels = f.labels()
    has_log = any(k == "y" and kap == 0 for (_c, k, kap) in labels)
    has_ei = any(k == "y" and kap > 0 for (_c, k, kap) in labels)
    coincidence = (a_d == a_c) or (a_d == a_p + a_c)
    assert has_log == coincidence
    assert has_ei == (a_d > a_c)
    assert (order == 4) == coincidence and (order == 5) == (not coincidence)
    assert sorted(f.rates()) == sorted({-a_d, -a_c, -(a_p + a_c), a_d})
    if (fam, za, zb) in (("p0", 3, 6), ("p1", 2, 6)):
        assert has_log and has_ei and a_d != a_c and a_d == a_p + a_c


def test_class_b_r_operators_frozen_point_agrees():
    """Rejects: this file and tests/test_paper58_r_operators.py disagreeing on the (3, 3) and (3, 4) orders the
    latter freezes (4 and 5)."""
    assert M.order_is(fn_hyb("p0", 3, 3), 4)
    assert M.order_is(fn_hyb("p0", 3, 4), 5)


# ---- independent value of the (3, 6) closed form: the cold review's from-scratch prolate quadrature (copied)
def _potential_2p0sq(Z, r, ct):
    ap = mp.mpf(Z)
    N2 = mp.mpf(Z) ** 5 / (32 * mp.pi)      # [2p0_Z]^2 = N2 r^2 e^{-Z r} (1/3 + 2/3 P2(cos))
    tot = mp.mpf(0)
    for L, cL in ((0, mp.mpf(1) / 3), (2, mp.mpf(2) / 3)):
        PL = mp.mpf(1) if L == 0 else (3 * ct * ct - 1) / 2
        inner = mp.gammainc(L + 5, 0, ap * r) / ap ** (L + 5)
        outer = mp.gammainc(4 - L, ap * r) / ap ** (4 - L)
        tot += 4 * mp.pi / (2 * L + 1) * cL * PL * N2 * (inner / r ** (L + 1) + outer * r ** L)
    return tot


def _hybrid_quadrature(ZA, ZB, Rv):
    Rv = mp.mpf(Rv)
    nc, nd = mp.sqrt(mp.mpf(ZA) ** 3 / mp.pi), mp.sqrt(mp.mpf(ZB) ** 3 / mp.pi)

    def integrand(xi, eta):
        rA, rB = Rv * (xi + eta) / 2, Rv * (xi - eta) / 2
        if rA == 0:
            return mp.mpf(0)
        ct = (1 + xi * eta) / (xi + eta)
        return (xi * xi - eta * eta) * _potential_2p0sq(ZA, rA, ct) * mp.exp(-ZA * rA - ZB * rB)

    scale = 2 / ((ZA + ZB) * Rv)
    cuts = [0, scale, 4 * scale, 12 * scale, 40 * scale, 120 * scale]
    val = mp.quad(lambda t, eta: integrand(1 + t, eta), cuts, [-1, -0.5, 0, 0.5, 1], maxdegree=6)
    return 2 * mp.pi * (Rv / 2) ** 3 * nc * nd * val


def test_class_b_second_coincidence_is_the_limit_of_the_generic_branch():
    """Rejects (fast tier): a wrong closed form AT the coincidence a_d = a_p + a_c, Z = (3, 6), whose order and
    atom content are right but whose weights are not (e.g. gamma + ln R with a wrong coefficient).  The
    coincidence form must be the limit of the generic-branch forms at Z_B = 6 -+ 1e-6 (which carry Ei(+-1e-6 R)
    instead of gamma + ln R): the symmetric second difference is O(h^2), here ~ 1e-13 relative, at R = 7/5."""
    x0, h = sp.Rational(7, 5), F(1, 10 ** 6)
    v = [sp.N(hyb_expr("p0", 3, F(6) + s * h).subs(R, x0), 40) for s in (-1, 0, 1)]
    assert abs((v[0] + v[2]) / 2 - v[1]) <= sp.Float(10) ** -12 * abs(v[1]), v


@pytest.mark.slow
def test_class_b_second_coincidence_closed_form_against_quadrature():
    """Rejects: a wrong closed form at the second coincidence a_d = a_p + a_c, Z = (3, 6) (the one point no
    existing guard covered), e.g. Ei and gamma + ln R mixed with a wrong weight.  Prolate quadrature with the
    potential of 2p0^2 from the Laplace expansion (no geovac integral code), R = 1.4; ~20 s."""
    with mp.workdps(18):
        q = _hybrid_quadrature(3, 6, mp.mpf("1.4"))
        cf = mp.mpf(str(sp.N(hyb_expr("p0", 3, 6).subs(R, sp.Rational(7, 5)), 30)))
        assert abs(q - cf) <= mp.mpf(10) ** -15 * abs(cf), (q, cf)


# ================================================================ C. equal-exponent exchange, order exactly 4
@pytest.mark.parametrize("case", XEQ_PARAMS)
def test_class_c_tau_sum_stops(case):
    """Rejects: a closed form whose tau sum was cut before it stops (the last two tau terms must vanish
    identically at tau_max = 6)."""
    terms = exch_terms(case)
    assert all(sp.simplify(t) == 0 for t in terms[-2:])


@pytest.mark.parametrize("case", XEQ_PARAMS)
def test_class_c_order_exactly_four(case):
    """Rejects: order 5 or order 3 accepted; a dropped ln R atom (the closure would be three functions); a
    tau sum cut early (the frozen exact determinant changes); a determinant that depends on gamma or a log
    (pinned gamma-free: gamma enters only as gamma x the ln R row); a decision made with gamma as a free symbol."""
    f = fn_exch(case)
    assert f.labels() == CLASS_C_LIST                         # e^{-2R}, e^{-2R} ln R, E1(2R), e^{2R} E1(4R)
    assert M.closure(f.labels()) == CLASS_C_LIST              # d/dR-closed, four functions: order <= 4
    for x0 in (X0, sp.Rational(11, 7)):
        mn = M.lower_minor(f, 4, x0, rows=CLASS_C_LIST)
        assert mn.nonzero
        assert mn.depends_on(f) == []
        assert mn.value.is_Rational
        assert mn.value == sp.Rational(DATA["exchange_equal_exponent_det4"][case][str(x0)])
    assert M.order_is(f, 4)
    assert not M.order_at_least(f, 5)


def test_class_c_1s_determinant_is_the_reviews_number():
    """Rejects: a recomputation that drifts from the cold review's own-code-path value
    |det| = 524742295998009665281408/1526641852514853284611875 at R = 37/59 (recomputed here, not pasted)."""
    mn = M.lower_minor(fn_exch("1s"), 4, X0, rows=CLASS_C_LIST)
    assert abs(mn.value) == sp.Rational(524742295998009665281408, 1526641852514853284611875)


def test_class_c_1s_against_sugiura():
    """Rejects: a wrong 1s exchange closed form (tau sum or E1 weights) passing because the guard only reads its
    own output: Sugiura's formula, typed independently, at R = 0.5, 1.4, 3 (30 digits)."""
    f = exch_expr("1s")
    with mp.workdps(40):
        for Rv in ("0.5", "1.4", "3"):
            x = mp.mpf(Rv)
            S = (1 + x + x ** 2 / 3) * mp.exp(-x)
            Sp = (1 - x + x ** 2 / 3) * mp.exp(x)
            sug = (-mp.exp(-2 * x) * (-mp.mpf(25) / 8 + 23 * x / 4 + 3 * x ** 2 + x ** 3 / 3)
                   + (6 / x) * (S ** 2 * (mp.euler + mp.log(x)) + Sp ** 2 * mp.ei(-4 * x)
                                - 2 * S * Sp * mp.ei(-2 * x))) / 5
            cf = mp.mpf(str(sp.N(f.subs(R, sp.Rational(Rv)), 45)))
            assert abs(cf - sug) <= mp.mpf(10) ** -30 * abs(cf), (Rv, cf, sug)


@pytest.mark.parametrize("case,z", [("1s", F(3, 2)), ("1s", 2), ("2psig", 2)])
def test_class_c_scaling_to_every_equal_exponent(case, z):
    """Rejects: order 4 failing at a = b != 1, or the scaling K(R; z) = z K(zR; 1) that carries order 4 to every
    a = b being false.  The scaled closure is the 1s list with every rate times z."""
    f = fn_exch(case, z)
    zz = sp.Rational(z)
    assert f.labels() == sorted([(c * zz, k, kap * zz) for (c, k, kap) in CLASS_C_LIST], key=str)
    assert M.order_is(f, 4)
    g = exch_expr(case)
    for Rv in (sp.Rational(7, 10), sp.Rational(9, 5)):
        a = sp.N(exch_expr(case, z).subs(R, Rv), 35)
        b = sp.N(zz * g.subs(R, zz * Rv), 35)
        assert abs(a - b) <= sp.Float(10) ** -28 * abs(a)


# ================================================================ FIRE TESTS (2026-11-18)
# Plants applied in memory only (a pytest plugin in the session scratchpad patched the imported modules or a
# helper function after import; no committed file was edited), the named -k selection was run, then the run was
# repeated clean.  (M) = tests/_minimal_order.py functions, (T) = this module's builders / constants.
#
#   plant (wrong answer)                                          -k selector                        result
#   class A order +1 (claimed) for every 1s case (T)              class_a                            FIRED
#   class A order -1 (T)                                          class_a                            FIRED
#   [text] (AA|BB) 1s order 4 at a != b                           class_a_1s and BB                  FIRED
#   [text] VB(aa) order 1 at a = b (rate 0 of 1/R not counted)    class_a_1s and VB                  FIRED
#   S + e^{-3R}E1(R)/1000 passing as class A (T)                  class_a_1s and S-                  FIRED
#   vanishing coefficient counted: coefficient_vanishes -> False  vanishing                          FIRED (reduce_logs=False path)
#   no prime-log reduction (M)                                    vanishing / dependent_constants    NOT FIRED (expected: true-
#                                                                                                      constant evaluation catches it)
#   both of the last two                                          vanishing                          FIRED
#   constants evaluated at rational stand-ins (M)                 dependent_constants                NOT FIRED (expected: prime-log
#                                                                                                      reduction catches it)
#   stand-ins AND no prime-log reduction (M)                      dependent_constants                FIRED
#   lower bound decided on the symbol polynomial (M)              dependent_constants                NOT FIRED (expected, as above)
#   same AND no reduction AND no vanishing check (M)              dependent_constants or vanishing   FIRED
#   merged rates not lowering: fn_1s(a, a) -> fn_1s(a, 2a) (T)    merging                            FIRED
#   upper bound from a point rank: generic_rank -> rank at 37/59  generic_rank_not_point_rank        FIRED
#     same plant                                                  class_b_order (9 cases)            NOT FIRED (expected: the real
#                                                                                                      hybrids are regular at 37/59)
#   [text] order_at_most: generic rank -> point rank              generic_rank_not_point_rank or     FIRED
#                                                                   class_b_order
#   [text] order_at_least tries 37/59 only                        generic_rank_not_point_rank        FIRED
#   [text] claim order 5 / rates 4 at (3, 6)                      class_b_order and p0-3-6           FIRED
#   [text] claim order 4 / rates 3 at (3, 4)                      class_b_order and p0-3-4           FIRED
#   [text] claim order 3 / rates 2 at p1 (2, 6)                   class_b_order and p1-2-6           FIRED
#   [text] claim 'gamma + ln R iff a_d = a_c'                     function_content                   FIRED (at (3, 6) and (2, 6))
#   gamma R e^{-6R} added to every hybrid (minor gamma-dependent) class_b_order                      FIRED (the gamma-free pin)
#   one E-atom dropped, e^{4R}E1(10R) at (3, 4) etc. (M)          class_b                            NOT FIRED: the planted function
#                                                                                                      has the same order (5, checked
#                                                                                                      directly); not a wrong answer
#   (3, 6) gamma coefficient x (1 + 1e-9) (T)                     limit_of_the_generic               FIRED
#     same plant                                                  class_b and 3-6 (structure only)   NOT FIRED (fast-tier structure
#                                                                                                      tests do not read values)
#     same plant, --slow                                          quadrature                         FIRED
#   (3, 6) closed form x (1 + 1e-12) (T)                          limit_of_the_generic               NOT FIRED (below its 1e-12
#                                                                                                      resolution)
#     same plant, --slow                                          quadrature                         FIRED
#   tau sum cut at tau_max = 2 for 3d orbitals (T)                class_c and 3ddel                  FIRED
#     same plant, --slow                                          class_c and 3dsig                  FIRED
#   ln R atom dropped from every Fn (M)                           class_c_order                      FIRED
#   gamma R e^{-2R} added (determinant gamma-dependent) (T)       class_c_order                      FIRED (the gamma-free pin line)
#   [text] class C order 5 accepted                               class_c_order and 1s               FIRED
#   [text] class C order 3 accepted                               class_c_order and ppsig            FIRED
#   [text] review literal numerator + 1                           reviews_number                     FIRED
#   1s tau term 0 x (1 + 1e-12) (T)                               sugiura                            FIRED
#     same plant                                                  class_c_order and 1s               FIRED (frozen determinant)
#   scaling applied to Z_A only (T)                               scaling                            FIRED
#   atom calculus: Ei / ln R derivative term sign flipped (M)     atom_calculus                      FIRED
#   [text] collector keeps E1 on its own label (M)                redundant_spelling                 FIRED
#   [text] interval half-width widened by 1 (M)                   class_b_order and p0-3-4           NOT FIRED: widening is
#                                                                                                      conservative (can only refuse);
#                                                                                                      the minor is rational x sqrt 3
# Clean runs after the plants: fast 79 passed, 3 skipped (124 s); --slow 82 passed (160 s).
