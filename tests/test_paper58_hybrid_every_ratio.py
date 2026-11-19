"""Paper 58 -- least order in R of four (AA|AB) hybrids with an L >= 1 same-centre pair, at EVERY Z_A, Z_B > 0.

Guard for the result derived and cold-reviewed 2026-11-19 (this file is the guard only, written as a separate
pass, CLAUDE.md Sec. 9).  For

    fam1 = (2p0 2p0 | 1s 1s_B),  H2 = (2p1 2p1 | 2p1 2p1_B),  H4 = (2p0 2p0 | 2p0 3d0_B),  H7 = (2p0 3d0 | 2s 1s_B)

with a_p the rate of the same-centre pair, a_c of the third A orbital, a_d of the B orbital:

    least order in R = 5  if a_d != a_c and a_d != a_c + a_p;   = 4  at a_d = a_c and at a_d = a_c + a_p.

Route.  The Proposition order(f) = d_E + k (hand proof, cold-reviewed; restated in tests/_hybrid_ratio.py), applied
to PARAMETRIC closed forms: each open exponent region (lt: a_d < a_c; mid: a_c < a_d < a_c + a_p; gt:
a_d > a_c + a_p) is swept by one symbol u > 0, the coincidences eq1 / eq2 are single forms; scaling
f(R; sZ_A, sZ_B) = s f(sR; Z_A, Z_B) makes the ratio a_d / a_c the only parameter.  The forms are FROZEN in
tests/data/paper58_hybrid_ratio_2026_11_19.json; every count is RECOMPUTED here from them (own code path:
tests/_hybrid_ratio.py; nothing from debug/), and the forms are tied back to the LIVE evaluator.

What each group pins:
  P. the Proposition itself against a DIRECT rank over Q(R) (tests/_minimal_order.py's atom calculus) on the cold
     review's eight toys, including the kernel interaction (an exponential coefficient inside L's kernel drops);
  M. the region maps: each u -> a_d/a_c is a bijection of (0, oo) onto its region; Z_A, Z_B > 0;
  A. the atom labels for every u > 0: f's rates, the E-atoms at -a_d and +a_d (the +a_d atoms are
     E1((a_d + a_c)R) and E1((a_d + a_p + a_c)R); NOT E1((a_p + a_d)R)), no two labels meeting at any u > 0;
  Q. d_E = 2, k = 3 (open) / 2 (coincidence), order 5 / 4; the E-part has the exact shape
     q(R) e^{-a_d R}[Ei((a_d-a_p-a_c)R) - Ei((a_d-a_c)R)] + eps q(-R) e^{a_d R}[E1((a_d+a_p+a_c)R) - E1((a_d+a_c)R)]
     (ln R in place of Ei(0 R) at a coincidence), eps = -1 for fam1 and H2, +1 for H4 and H7;
  Z. "for EVERY u > 0": q and every coefficient of L f are free of gamma / logs / pi, the gcd over Q[u] of their
     coefficients is a nonzero rational, and no u-factor of any denominator (f, q, L, L f) has a root u > 0;
  V. the frozen form = live hybrid_closed_form / (exact product of the four radial_norm values) at fresh u (open
     regions) and fresh scales (coincidences): a form off by any u-dependent factor has an R-constant ratio and
     must still FAIL;
  D. the order at four fresh exponent points by the direct certified route (interval minor + generic rank).

CONDITIONAL ON (not tested here): Lemma GI (atom independence; hand proof) and the Proposition's hand proof (only
its agreement with direct rank on toys is tested); that hybrid_closed_form is the integral (validated against
quadrature elsewhere: tests/test_paper58_minimal_order.py).  NOT ASSERTED: anything about the other six panel
hybrids beyond tests/test_paper58_minimal_order.py; that the frozen forms are what the symbolic build produces
today (the build takes 1-15 min per region; group V ties the frozen forms to the evaluator at points instead).

FAST / SLOW.  Fast: P, M, A (all 20 forms, labels only), Q and Z on fam1 (all five) and on the gt region plus
both coincidences of H2, H4, H7; V on fam1 mid at two fresh u.  Slow: Q and Z on the remaining six open forms
(H4 lt alone ~50 s), V everywhere else (each live build 3-32 s), and D (four live builds + certified ranks).

Each test's docstring names the wrong answer it rejects; fire-test plants are listed at the end of the file.
"""
from __future__ import annotations

import ast
from fractions import Fraction
from functools import lru_cache

import mpmath as mp
import pytest
import sympy as sp
from sympy.polys.fields import field
from sympy.polys.matrices import DomainMatrix

import _hybrid_ratio as H
import _minimal_order as M
from geovac.two_center_eri import hybrid_closed_form, radial_norm, R_s

U = H.U
DATA = H.load_data()

# ------------------------------------------------------------------ expectations (read at run time, not collected)
ORDER = {"lt": 5, "mid": 5, "gt": 5, "eq1": 4, "eq2": 4}
K_EXPECTED = {"lt": 3, "mid": 3, "gt": 3, "eq1": 2, "eq2": 2}
D_E_EXPECTED = 2
EPS = {"fam1": -1, "H2": -1, "H4": 1, "H7": 1}

TOY_EXPECTED = {   # name -> (direct order, d_E, k): the cold review's table
    "T1 generic": (5, 1, 4), "T2 r in kernel": (4, 1, 3), "T3 closure-rate exp": (3, 1, 2),
    "T4 independent at one shift": (5, 2, 3), "T5 dependent triple": (5, 2, 3), "T6 mirror pair": (5, 2, 3),
    "T7 log + Ei": (3, 1, 2), "T8 near-miss": (3, 1, 2),
}

FAST_ANALYSIS = {("fam1", r) for r in H.REGIONS} | {(c, r) for c in ("H2", "H4", "H7") for r in ("gt", "eq1", "eq2")}
CASE_REGIONS = [(c, r) for c in H.CASES for r in H.REGIONS]


def _param(c, r, fast):
    return pytest.param(c, r, id=f"{c}-{r}", marks=() if fast else pytest.mark.slow)


ANALYSIS_PARAMS = [_param(c, r, (c, r) in FAST_ANALYSIS) for c, r in CASE_REGIONS]
ALL_PARAMS = [pytest.param(c, r, id=f"{c}-{r}") for c, r in CASE_REGIONS]


def entry(case, region):
    return DATA["cases"][case]["regions"][region]


def info(case, region):
    return {k: sp.expand(sp.sympify(v, locals={"u": U})) for k, v in entry(case, region)["info"].items()}


def labels_of(case, region):
    """Atom labels (c, kind, kappa) read from the frozen keys (every frozen coefficient is nonzero)."""
    out = set()
    for ks in entry(case, region)["form"]:
        c, kind, kap, _consts = ast.literal_eval(ks)
        out.add((sp.expand(sp.sympify(c, locals={"u": U})), kind, sp.expand(sp.sympify(kap, locals={"u": U}))))
    return out


@lru_cache(maxsize=None)
def form(case, region) -> H.Form:
    return H.Form(entry(case, region))


@lru_cache(maxsize=None)
def analysis(case, region) -> dict:
    f = form(case, region)
    return H.quotient_count(f.coords, f.K, f.RG)


def e_label(c, kap):
    c, kap = sp.expand(c), sp.expand(kap)
    return (c, "log", sp.Integer(0)) if kap == 0 else (c, "Ei", kap)


# ================================================================ P. the Proposition against direct rank (toys)
def direct_order_toy(terms) -> int:
    """Rank over Q(R) of f, f', ..., f^(N) in the d/dR-closed atom basis (N = closure size), on
    tests/_minimal_order.py's own atom calculus -- independent of tests/_hybrid_ratio.py."""
    coords = {}
    for cf, (c, kind, k) in terms:
        key = ((sp.Rational(c), kind, sp.Rational(k)), ())
        coords[key] = coords.get(key, M.FR.zero) + M.FR.from_sympy(sp.sympify(cf))
    coords = {k: v for k, v in coords.items() if v != 0}
    N = len(M.closure({lb for (lb, _m) in coords})) + 1
    ds = [coords]
    for _ in range(N):
        ds.append(M.derivative(ds[-1]))
    keys = sorted(set().union(*[set(D) for D in ds]), key=str)
    Mx = DomainMatrix([[D.get(k, M.FR.zero) for D in ds] for k in keys], (len(keys), len(ds)), M.FR)
    return int(Mx.rank())


@pytest.mark.parametrize("name", list(H.toys()))
def test_proposition_equals_direct_rank_on_toys(name):
    """Rejects: the quotient count d_E + k disagreeing with the direct minimal order -- in particular k counted
    from f's own exponential rates instead of from L f (T1 would give 1 + 2 = 3, not 5; T2, whose exponential
    coefficient 5q lies in L's kernel, would give 3, not 4), d_E counted as the number of E-atoms instead of the
    dimension of their coefficient span (T5: 3, not 2; T1: 2, not 1), or a dropped closure rate."""
    terms = H.toys()[name]
    want_order, want_dE, want_k = TOY_EXPECTED[name]
    K, _u, RG = field([U, H.RS], sp.QQ)
    res = H.quotient_count(H.toy_coords(terms, K, RG), K, RG)
    direct = direct_order_toy(terms)
    assert direct == want_order, (name, direct)
    assert (res["d_E"], res["k"]) == (want_dE, want_k), (name, res["d_E"], res["k"])
    assert res["order"] == direct


def test_kernel_toy_named_in_the_task():
    """Rejects: the kernel interaction being missed.  f = q e^{-R}[Ei(R) - Ei(2R)] + 5q e^{-R} + e^{-3R} has
    order 4 = 1 + 3 (the rate -1 drops out of L f), while the generic variant with R in place of 5q has order
    5 = 1 + 4; a count of k from f's own rates gives 3 for both, so the trap is real."""
    toys = H.toys()
    K, _u, RG = field([U, H.RS], sp.QQ)
    for name, want in (("T2 r in kernel", 4), ("T1 generic", 5)):
        res = H.quotient_count(H.toy_coords(toys[name], K, RG), K, RG)
        assert res["d_E"] + res["k"] == want == direct_order_toy(toys[name])
        f_rates = {lb[0] for (_cf, lb) in toys[name] if lb[1] == "1"}
        assert res["d_E"] + len(f_rates) != want                     # counting f's rates is wrong here
    assert sp.Integer(-1) not in H.quotient_count(H.toy_coords(toys["T2 r in kernel"], K, RG), K, RG)["Lf_rates"]


# ================================================================ M. region maps
_ENDPOINTS = {"lt": ("0", "1"), "mid": ("1", "1+rho"), "gt": ("1+rho", "oo")}


@pytest.mark.parametrize("case,region", ALL_PARAMS)
def test_region_map_is_a_bijection_onto_its_region(case, region):
    """Rejects: a region parametrisation that misses part of its ratio interval (so 'every u > 0' would not be
    'every exponent choice'), a Z that is not positive for some u > 0, or rates inconsistent with the orbitals
    (a_p = Z_A (1/n_a + 1/n_b), a_c = Z_A / n_c, a_d = Z_B / n_d)."""
    a, b, c, d = (H.ORB[o] for o in H.CASES[case])
    I = info(case, region)
    ZA, ZB = I["ZA"], I["ZB"]
    assert sp.expand(I["a_p"] - ZA * (sp.Rational(1, a[0]) + sp.Rational(1, b[0]))) == 0
    assert sp.expand(I["a_c"] - ZA / c[0]) == 0 and sp.expand(I["a_d"] - ZB / d[0]) == 0
    for Z in (ZA, ZB):
        assert not H.lin_vanishes_for_some_u(Z) and Z.subs(U, 1) > 0
    rho = sp.Rational(1, a[0]) * c[0] + sp.Rational(1, b[0]) * c[0]          # a_p / a_c, a constant
    assert sp.simplify(I["a_p"] / I["a_c"] - rho) == 0
    r = sp.cancel(I["a_d"] / I["a_c"])
    if region == "eq1":
        assert r == 1
        return
    if region == "eq2":
        assert r == 1 + rho
        return
    num = sp.numer(sp.together(sp.diff(r, U)))
    assert num != 0 and not H.positive_roots(num)                      # strictly monotone on u > 0
    ends = {"0": 0, "1": 1, "1+rho": 1 + rho, "oo": sp.oo}
    lim = sorted([sp.limit(r, U, 0, "+"), sp.limit(r, U, sp.oo)], key=lambda x: (x == sp.oo, x))
    assert lim == [ends[e] for e in _ENDPOINTS[region]], (lim, region)


# ================================================================ A. atom labels, for every u > 0
def expected_labels(case, region):
    I = info(case, region)
    ap, ac, ad = I["a_p"], I["a_c"], I["a_d"]
    E_minus = {e_label(-ad, ad - ac), e_label(-ad, ad - ap - ac)}
    E_plus = {e_label(ad, -(ad + ac)), e_label(ad, -(ad + ap + ac))}
    rates = {sp.expand(x) for x in (-ad, -ac, -(ap + ac))}
    return E_minus, E_plus, rates


@pytest.mark.parametrize("case,region", ALL_PARAMS)
def test_atom_labels(case, region):
    """Rejects: the earlier memo's typo E1((a_p + a_d)R) as a +a_d atom (for H2, H4, H7 it is a different
    function at every u > 0 and absent; fam1 has a_p = a_c and cannot tell), a missing or extra E-atom at either
    shift, Ei where ln R belongs at a coincidence, or f's exponential rates other than {-a_d, -a_c, -(a_p+a_c)}
    (merged at a coincidence)."""
    labs = labels_of(case, region)
    I = info(case, region)
    E_minus, E_plus, rates = expected_labels(case, region)
    got_E = {lb for lb in labs if lb[1] != "1"}
    assert {lb for lb in got_E if lb[0] == sp.expand(-I["a_d"])} == E_minus
    assert {lb for lb in got_E if lb[0] == I["a_d"]} == E_plus
    assert got_E == E_minus | E_plus
    assert {lb[0] for lb in labs if lb[1] == "1"} == rates
    typo = e_label(I["a_d"], -(I["a_p"] + I["a_d"]))
    if case != "fam1":
        assert not H.lin_vanishes_for_some_u(I["a_p"] - I["a_c"])       # the typo label is a different function
        assert typo not in got_E
    if region == "eq1" or region == "eq2":
        assert any(lb[1] == "log" for lb in got_E)
        assert not form_has_no_gamma(case, region)                       # gamma IS in f at a coincidence


def form_has_no_gamma(case, region):
    return not any("GAM" in ks for ks in entry(case, region)["form"])


@pytest.mark.parametrize("case,region", ALL_PARAMS)
def test_no_two_labels_meet_at_any_u(case, region):
    """Rejects: an exceptional ratio inside an open region -- two rates (of f, of the d/dR closure, or of L f's
    possible rates), two E-labels at one shift, or an Ei argument meeting / vanishing at some u > 0 (Lemma GI
    needs distinct labels at the exponents in question)."""
    labs = labels_of(case, region)
    rates = {lb[0] for lb in labs if lb[1] == "1"} | {sp.expand(lb[0] + lb[2]) for lb in labs if lb[1] != "1"}
    rates = sorted(rates, key=str)
    for i, x in enumerate(rates):
        for y in rates[i + 1:]:
            assert not H.lin_vanishes_for_some_u(x - y), (x, y)
    E = sorted((lb for lb in labs if lb[1] != "1"), key=str)
    for i, x in enumerate(E):
        if x[1] == "Ei":
            assert not H.lin_vanishes_for_some_u(x[2]), x
        for y in E[i + 1:]:
            if x[0] == y[0]:
                assert not H.lin_vanishes_for_some_u(x[2] - y[2]), (x, y)
            else:
                assert not H.lin_vanishes_for_some_u(x[0] - y[0]), (x, y)


# ================================================================ Q. the quotient count on the frozen forms
@pytest.mark.parametrize("case,region", ANALYSIS_PARAMS)
def test_quotient_count_gives_five_or_four(case, region):
    """Rejects: order 4 at a generic ratio or 5 at a coincidence; d_E = 4 (the two E-atoms at one shift treated as
    independent lines: order 7); L failing to kill the E-part; k read from anything but L f's rates (the
    surviving rates must be exactly {-a_d, -a_c, -(a_p + a_c)}, merged at a coincidence)."""
    res = analysis(case, region)
    _Em, _Ep, rates = expected_labels(case, region)
    assert res["d_E"] == D_E_EXPECTED and all(v == 1 for v in res["per_shift"].values()), res["per_shift"]
    assert len(res["per_shift"]) == 2
    assert set(res["Lf_rates"]) == rates
    assert res["k"] == K_EXPECTED[region] == len(rates)
    assert res["order"] == ORDER[region]


@pytest.mark.parametrize("case,region", ANALYSIS_PARAMS)
def test_e_part_shape_and_epsilon(case, region):
    """Rejects: a wrong eps (fam1, H2 -1; H4, H7 +1), E-coefficients at one shift that are not exactly opposite
    (d_E would not be 2), the +a_d pair not the R -> -R mirror of the -a_d pair, or a q carrying gamma / a log /
    pi.  Shape: q e^{-a_d R}[Ei((a_d-a_p-a_c)R) - Ei((a_d-a_c)R)] + eps q(-R) e^{a_d R}[E1((a_d+a_p+a_c)R) -
    E1((a_d+a_c)R)], with E1(x) = -Ei(-x) and ln R in place of Ei(0 R)."""
    f = form(case, region)
    I = info(case, region)
    ap, ac, ad = I["a_p"], I["a_c"], I["a_d"]
    co = f.coords
    q = co[e_label(-ad, ad - ap - ac)]
    assert f.const_free(q)
    assert co[e_label(-ad, ad - ac)] == -q
    q_ref = f.K.from_expr(q.as_expr().subs(H.RS, -H.RS))
    eps = EPS[case]
    assert co[e_label(ad, -(ad + ap + ac))] == -eps * q_ref         # coefficient of E1((a_d+a_p+a_c)R) = eps q(-R)
    assert co[e_label(ad, -(ad + ac))] == eps * q_ref


@pytest.mark.parametrize("case,region", ANALYSIS_PARAMS)
def test_q_and_Lf_nonvanishing_and_constant_free_for_every_u(case, region):
    """Rejects: a gamma, log or pi constant surviving into q or L f (then 'nonzero' would rest on the arithmetic
    of the constants); a zero-gcd of q or of an L f coefficient that is a polynomial in u (a root u > 0 would be
    an exceptional ratio); a u-factor of any denominator (f, q, L's coefficients, L f) with a root u > 0 (the form
    would not be the integral there).  The constants ARE present in f (logs of rate ratios; gamma at a
    coincidence), so 'constant-free' is a real filter."""
    f = form(case, region)
    res = analysis(case, region)
    assert f.const_names, "f carries no constants: the constant-free check would be vacuous"
    I = info(case, region)
    q = f.coords[e_label(-I["a_d"], I["a_d"] - I["a_p"] - I["a_c"])]
    for name, v in [("q", q)] + [(f"Lf[{r}]", v) for (r, _k, _x), v in res["Lf"].items()]:
        assert f.const_free(v), name
        g = H.u_gcd(v, f.K)
        assert g.is_Rational and g != 0, (name, g)
        for fa in H.u_denominator_factors(v):
            assert not H.positive_roots(fa), (name, fa)
    for v in list(f.coords.values()) + list(res["A"]):
        for fa in H.u_denominator_factors(v):
            assert not H.positive_roots(fa), fa


# ================================================================ V. frozen forms = live evaluator / norm product
FRESH = {   # (case, region) -> fresh u (open) or fresh scale s (coincidence); none used by the memo or the review
    ("fam1", "lt"): ["19/6"], ("fam1", "mid"): ["5/13", "19/6"], ("fam1", "gt"): ["5/13"],
    ("fam1", "eq1"): ["7/5"], ("fam1", "eq2"): ["3/11"],
    ("H2", "lt"): ["5/13"], ("H2", "mid"): ["19/6"], ("H2", "gt"): ["7/17"],
    ("H2", "eq1"): ["7/5"], ("H2", "eq2"): ["3/11"],
    ("H4", "lt"): ["19/6"], ("H4", "mid"): ["5/13"], ("H4", "gt"): ["7/17"],
    ("H4", "eq1"): ["3/11"], ("H4", "eq2"): ["7/5"],
    ("H7", "lt"): ["7/17"], ("H7", "mid"): ["19/6"], ("H7", "gt"): ["5/13"],
    ("H7", "eq1"): ["7/5"], ("H7", "eq2"): ["3/11"],
}
FAST_LIVE = {("fam1", "mid")}
LIVE_PARAMS = [pytest.param(c, r, x, id=f"{c}-{r}-{x}", marks=() if (c, r) in FAST_LIVE else pytest.mark.slow)
               for (c, r), xs in FRESH.items() for x in xs]
R_LIVE = (sp.Rational(3, 7), sp.Rational(5, 2))


def _frac(x) -> Fraction:
    x = sp.Rational(x)
    return Fraction(int(x.p), int(x.q))


def norm_product(case, ZA, ZB):
    orbs = [H.ORB[o] for o in H.CASES[case]]
    out = sp.Integer(1)
    for o, Z in zip(orbs, (ZA, ZA, ZA, ZB)):
        out *= radial_norm(_frac(Z), o[0], o[1])
    return out


def test_fresh_points_are_fresh():
    """Rejects: a 'fresh' point that the memo (u = 1/3, 1/2, 1, 2, 3, 5) or the review (u = 1e-6 .. 1000 list)
    already used, and fewer than two fresh u per case."""
    used = {sp.Rational(x) for x in ("1/3", "1/2", "1", "2", "3", "5", "1/1000000", "3/7", "1000", "1/100000",
                                      "11/3", "1/10000", "9/4", "300", "2/7", "50")}
    for case in H.CASES:
        us = {sp.Rational(x) for (c, r), xs in FRESH.items() if c == case and r in H.OPEN for x in xs}
        assert len(us) >= 2 and not (us & used), (case, us)
        assert all(r in [rr for (c, rr) in FRESH if c == case] for r in H.REGIONS)


@pytest.mark.parametrize("case,region,x", LIVE_PARAMS)
def test_frozen_form_is_live_evaluator_over_norm_product(case, region, x):
    """Rejects: a frozen form that is not the integral -- including one off by a u-dependent factor, whose ratio
    to hybrid_closed_form is constant in R and would pass a 'ratio constant in R' check.  The ratio must equal the
    EXACT product of the four radial_norm values (geovac, current working copy), at two R, to 1e-40 (80 digits).
    Open region: at u = x.  Coincidence (no u): at the scaled exponents s Z, through f(R; sZ) = s f(sR; Z)."""
    I = info(case, region)
    orbs = [H.ORB[o] for o in H.CASES[case]]
    x = sp.Rational(x)
    if region in H.OPEN:
        u0, s = x, sp.Integer(1)
    else:
        assert not I["ZA"].has(U)
        u0, s = sp.Integer(1), x
    ZA, ZB = sp.Rational(I["ZA"].subs(U, u0)), sp.Rational(I["ZB"].subs(U, u0))
    live = hybrid_closed_form(_frac(s * ZA), *orbs[:3], _frac(s * ZB), orbs[3])
    N = norm_product(case, ZA, ZB)
    with mp.workdps(80):
        Nv = mp.mpf(str(sp.N(N, 90)))
        for Rv in R_LIVE:
            lv = mp.mpf(str(sp.N(live.subs(R_s, Rv), 90)))
            fv = s_mp(s) * Nv * H.form_value(entry(case, region), u0, s * Rv, dps=90)
            assert abs(lv - fv) <= mp.mpf(10) ** -40 * abs(lv), (case, region, x, Rv, lv, fv)


def s_mp(s):
    return mp.mpf(int(s.p)) / int(s.q)


# ================================================================ D. direct certified order at fresh exponents
DIRECT = [   # (case, Z_A, Z_B, where, order) -- the cold review's fresh points (not used by the memo)
    ("H7", "7", "7/2", "a_d = a_c", 4),
    ("H4", "4/3", "6", "a_d = a_c + a_p", 4),
    ("H7", "7", "13/5", "lt", 5),
    ("H2", "5", "17/3", "mid", 5),
]
DIRECT_ORDER = {(c, za, zb): o for c, za, zb, _w, o in DIRECT}


@pytest.mark.slow
@pytest.mark.parametrize("case,za,zb,where", [d[:4] for d in DIRECT], ids=[f"{d[0]}-{d[1]}-{d[2]}" for d in DIRECT])
def test_direct_certified_order_at_fresh_exponents(case, za, zb, where):
    """Rejects: order 4 at a generic ratio (H7 (7, 13/5), H2 (5, 17/3)) or 5 at a coincidence (H7 (7, 7/2):
    a_d = a_c; H4 (4/3, 6): a_d = a_c + a_p), on the LIVE closed form by the route of
    tests/test_paper58_minimal_order.py: an interval-certified r x r minor (lower bound), the generic rank over
    Q(R) (upper bound), and no (r+1)-minor nonzero.  Independent of the frozen forms and of the Proposition."""
    a, b, c, d = (H.ORB[o] for o in H.CASES[case])
    ZA, ZB = sp.Rational(za), sp.Rational(zb)
    ap, ac, ad = ZA / a[0] + ZA / b[0], ZA / c[0], ZB / d[0]
    coincidence = (ad == ac) or (ad == ac + ap)
    assert coincidence == (where.startswith("a_d ="))
    order = DIRECT_ORDER[(case, za, zb)]
    f = M.Fn(hybrid_closed_form(_frac(ZA), a, b, c, _frac(ZB), d))
    assert len(f.rates()) == order - 1                       # rates + 1
    assert M.order_at_least(f, order)
    assert M.order_at_most(f, order)
    assert not M.order_at_least(f, order + 1)
    assert (order == 4) == coincidence


# ================================================================ FIRE TESTS (2026-11-19)
# Plants applied in memory only (a pytest plugin in the session scratchpad patched this module's expectations, the
# helper tests/_hybrid_ratio.py, or this module's in-memory copy of the frozen data after collection; no file was
# edited), the named -k selection was run, then clean runs were repeated.  (D) = frozen data in memory,
# (H) = helper function, (T) = this module's expectations.
#
#   plant (wrong answer)                                         -k selector                      result
#   DIRECT_ORDER flipped 4 <-> 5 at all four fresh points (T)    direct_certified --slow          FIRED 4/4
#   ORDER open regions = 4 (T)                                   quotient_count                   FIRED 6/6 fast
#     same plant                                                 quotient_count and H4-lt --slow  FIRED
#   ORDER coincidences = 5 (T)                                   quotient_count                   FIRED 8/8
#   d_E = number of E-atoms (span_basis keeps every q) (H)       quotient_count or toys           FIRED 21 (T1..T8 bar one,
#                                                                                                   every analysed form)
#   k counted from f's own exponential rates (H)                 toys or kernel_toy               FIRED 7/9 (T1-T5, T8 + named)
#     same plant                                                 quotient_count                   NOT FIRED (expected: on the
#                                                                                                   real forms every rate of f
#                                                                                                   survives into L f; only the
#                                                                                                   toys separate the two counts)
#   +a_d atom = E1((a_p + a_d)R) (the old typo) (T)              atom_labels                      FIRED 15/20 (H2, H4, H7);
#                                                                                                   fam1 NOT FIRED (a_p = a_c:
#                                                                                                   same function; expected)
#   fam1 mid form x 1/(u + 2) (R-constant ratio) (D)             fam1-mid                         FIRED: live test 2/2 only (the
#                                                                                                   structure is unchanged: expected)
#   every open H7 form x 1/(u + 2) (D)                           live and H7 --slow               FIRED 3/3
#   coincidence forms x 3/2 (D)                                  live and fam1-eq --slow          FIRED 2/2
#   one exponential coefficient x (1 + 1e-20), fam1 mid (D)      fam1-mid                         FIRED: live 2/2
#   gamma R e^{-a_c R} added, fam1 mid (D)                       fam1-mid                         FIRED: constant-free test +
#                                                                                                   live 2/2; quotient_count NOT
#                                                                                                   FIRED (order still 5: expected)
#   gamma on one E-coefficient, fam1 mid (D)                     fam1-mid                         FIRED: e_part, constant-free,
#                                                                                                   live 2/2
#   whole fam1 mid form x (2u - 1): gcds get root u = 1/2 (D)    nonvanishing and fam1-mid        FIRED
#   whole fam1 mid form / (u - 3): denominator root u = 3 (D)    nonvanishing and fam1-mid        FIRED
#   EPS[H7] = -1 (T)                                             e_part                           FIRED 3/3 fast (gt, eq1, eq2)
#   EPS[fam1] = +1 (T)                                           e_part                           FIRED 5/5
#   fam1 lt map covering only (0, 1/2) (D)                       region_map                       FIRED (endpoint check)
#   exp rate -u-1 relabelled -3u+1 (meets -2u-1 at u = 2) (D)    no_two_labels                    FIRED
#   TOY_EXPECTED T2 = order 5 (T)                                toys or kernel                   FIRED
#   Wronskian rows built without the e^{cR} shift (H)            toys or quotient_count           FIRED 22 (L leaves E-atoms)
#
# Clean runs after the plants: fast 114 passed, 41 skipped (~51 s); --slow 155 passed (~453 s).
