"""Paper 58 -- least order in R of six more (AA|AB) hybrids of the l <= 2 panel, at EVERY Z_A, Z_B > 0.

Guard (written as a separate pass, CLAUDE.md Sec. 9) extending tests/test_paper58_hybrid_every_ratio.py from
fam1, H2, H4, H7 to

    H1  = (3d0 3d0 | 3d0 1s_B),   H3 = (3d0 3d0 | 1s 1s_B),    H5  = (3d1 3d1 | 3d1 2p1_B),
    H6  = (3d2 3d2 | 3s 2s_B),    H8 = (3p1 3d1 | 3d0 2p0_B),  H10 = (3d2 3d0 | 3d-1 3d1_B)

with a_p the rate of the same-centre pair, a_c of the third A orbital, a_d of the B orbital:

    least order in R = 5  if a_d != a_c and a_d != a_c + a_p;   = 4  at a_d = a_c and at a_d = a_c + a_p.

Route: the same Proposition order(f) = d_E + k (tests/_hybrid_ratio.py, reused UNCHANGED), on PARAMETRIC closed
forms, one symbol u > 0 per open exponent region (lt, mid, gt) plus the two coincidences (eq1, eq2).  All 30 forms
are frozen IN FULL (0.8 MB; no summaries) in tests/data/paper58_hybrid_ratio_panel_2026_11_19.json, copied verbatim
from the 2026-11-19 parallel build (fields result.info, form of the run files); every count below is recomputed
here from them.  Nothing is imported from debug/.

Groups (names as in the existing guard):
  S. data integrity: six cases, five regions each, orbitals as named; plus the kernel-toy sentinel (the only place
     a "k counted from f's own rates" count can be told apart: on the real forms every rate of f survives into L f);
  M. region maps: u -> a_d/a_c a bijection of (0, oo) onto the region; Z_A, Z_B > 0; rates consistent with orbitals;
  A. atom labels for every u > 0: E-atoms at -a_d and +a_d exactly, f's rates {-a_d, -a_c, -(a_p+a_c)}, no two
     labels (or rates, or Ei arguments) meeting at any u > 0;
  Q. d_E = 2 (one line per shift), k = 3 open / 2 coincidence, order 5 / 4, k read from L f; E-part
     q e^{-a_d R}[Ei((a_d-a_p-a_c)R) - Ei((a_d-a_c)R)] + eps q(-R) e^{a_d R}[E1((a_d+a_p+a_c)R) - E1((a_d+a_c)R)],
     eps = +1 for H5, -1 for H1, H3, H6, H8, H10 (MEASURED from the frozen forms 2026-11-19, then frozen here);
  Z. q and every L f coefficient free of gamma / logs, zero-gcd over Q[u] a nonzero rational, no u-factor of any
     denominator (f, q, L, L f) with a root u > 0;
  V. frozen form = live hybrid_closed_form / exact product of the four radial_norm values, at ONE fresh point per
     case (90 digits, tolerance 1e-40; the build memo found 40 digits too few for H5 and H10).

CONDITIONAL ON (not tested here): Lemma GI and the Proposition's hand proof (its agreement with a direct rank on the
review's toys is tested in tests/test_paper58_hybrid_every_ratio.py); that hybrid_closed_form is the integral
(tests/test_paper58_minimal_order.py).  GAPS: (i) V ties only ONE of the five forms per case to the evaluator (one
fresh point each; the other 24 forms are tied only through the shared build and through their mutual consistency
under the label checks); (ii) no direct certified-rank check at fresh exponents for these six cases (the existing
guard's group D covers fam1, H2, H4, H7 only); (iii) the frozen forms are not rebuilt from the symbolic algebra
here (the parallel build took ~97 min on 16 cores).

FAST / SLOW.  Fast: S, M, A on all 30 forms; Q and Z on both coincidences of every case and on one open region per
case for H1, H3, H5, H6, H8 (the cheapest).  Slow: Q and Z on the remaining 13 open forms (up to ~35 s each);
V (each live build 14-80 s).

Each test's docstring names the wrong answer it rejects; fire-test plants are listed at the end of the file.
"""
from __future__ import annotations

import ast
import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp
from sympy.polys.fields import field
from sympy.polys.matrices import DomainMatrix

import _hybrid_ratio as H
import _minimal_order as M
from geovac.two_center_eri import hybrid_closed_form, radial_norm, R_s

U = H.U
DATA_PATH = Path(__file__).parent / "data" / "paper58_hybrid_ratio_panel_2026_11_19.json"
DATA = json.loads(DATA_PATH.read_text())

ORB = {"1s": (1, 0, 0), "2s": (2, 0, 0), "2p0": (2, 1, 0), "2p1": (2, 1, 1), "3s": (3, 0, 0), "3p1": (3, 1, 1),
       "3d0": (3, 2, 0), "3d1": (3, 2, 1), "3d2": (3, 2, 2), "3d-1": (3, 2, -1)}
CASES = {"H1": ("3d0", "3d0", "3d0", "1s"), "H3": ("3d0", "3d0", "1s", "1s"), "H5": ("3d1", "3d1", "3d1", "2p1"),
         "H6": ("3d2", "3d2", "3s", "2s"), "H8": ("3p1", "3d1", "3d0", "2p0"),
         "H10": ("3d2", "3d0", "3d-1", "3d1")}
REGIONS = H.REGIONS
OPEN = H.OPEN

# ------------------------------------------------------------------ expectations (read at run time)
ORDER = {"lt": 5, "mid": 5, "gt": 5, "eq1": 4, "eq2": 4}
K_EXPECTED = {"lt": 3, "mid": 3, "gt": 3, "eq1": 2, "eq2": 2}
D_E_EXPECTED = 2
EPS = {"H1": -1, "H3": -1, "H5": 1, "H6": -1, "H8": -1, "H10": -1}

FAST_OPEN = {"H1": "gt", "H3": "mid", "H5": "gt", "H6": "gt", "H8": "gt"}     # H10 open forms 14-27 s each: slow
FAST_ANALYSIS = {(c, r) for c in CASES for r in ("eq1", "eq2")} | set(FAST_OPEN.items())
CASE_REGIONS = [(c, r) for c in CASES for r in REGIONS]


def _param(c, r, fast):
    return pytest.param(c, r, id=f"{c}-{r}", marks=() if fast else pytest.mark.slow)


ANALYSIS_PARAMS = [_param(c, r, (c, r) in FAST_ANALYSIS) for c, r in CASE_REGIONS]
ALL_PARAMS = [pytest.param(c, r, id=f"{c}-{r}") for c, r in CASE_REGIONS]


def entry(case, region):
    return DATA["cases"][case]["regions"][region]


def info(case, region):
    return {k: sp.expand(sp.sympify(v, locals={"u": U})) for k, v in entry(case, region)["info"].items()}


def labels_of(case, region):
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


# ================================================================ S. data integrity + kernel sentinel
def test_data_has_six_cases_five_regions_named_orbitals():
    """Rejects: a missing case or region (a silently skipped sweep), or a case whose frozen orbitals are not the
    hybrid named in the docstring (a form of a different integral passing as H5, say)."""
    assert set(DATA["cases"]) == set(CASES)
    for c, orbs in CASES.items():
        assert tuple(DATA["cases"][c]["orbitals"]) == orbs
        assert set(DATA["cases"][c]["regions"]) == set(REGIONS)
        for r in REGIONS:
            assert entry(c, r)["form"] and set(entry(c, r)["info"]) == {"ZA", "ZB", "a_p", "a_c", "a_d"}


def _direct_order_toy(terms) -> int:
    """Rank over Q(R) of f, f', ..., f^(N) on tests/_minimal_order.py's atom calculus (independent of the helper)."""
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


def test_kernel_sentinel_k_from_Lf_not_from_f():
    """Rejects: k counted from f's own exponential rates.  On the six panel forms every rate of f survives into
    L f, so that wrong count gives the right answer there; this sentinel runs the SAME quotient_count on the
    kernel toy (5q e^{-R} inside L's kernel: order 4 = 1 + 3, direct rank 4; f's own rates would give 1 + 2 = 3)
    and its generic twin (order 5; f's rates would give 3), against a direct rank."""
    toys = H.toys()
    K, _u, RG = field([U, H.RS], sp.QQ)
    for name, want in (("T2 r in kernel", 4), ("T1 generic", 5)):
        res = H.quotient_count(H.toy_coords(toys[name], K, RG), K, RG)
        assert res["order"] == want == _direct_order_toy(toys[name]), (name, res["order"])


# ================================================================ M. region maps
_ENDPOINTS = {"lt": ("0", "1"), "mid": ("1", "1+rho"), "gt": ("1+rho", "oo")}


@pytest.mark.parametrize("case,region", ALL_PARAMS)
def test_region_map_is_a_bijection_onto_its_region(case, region):
    """Rejects: a region parametrisation missing part of its ratio interval, a Z not positive for some u > 0, or
    rates inconsistent with the orbitals (a_p = Z_A (1/n_a + 1/n_b), a_c = Z_A / n_c, a_d = Z_B / n_d)."""
    a, b, c, d = (ORB[o] for o in CASES[case])
    I = info(case, region)
    ZA, ZB = I["ZA"], I["ZB"]
    assert sp.expand(I["a_p"] - ZA * (sp.Rational(1, a[0]) + sp.Rational(1, b[0]))) == 0
    assert sp.expand(I["a_c"] - ZA / c[0]) == 0 and sp.expand(I["a_d"] - ZB / d[0]) == 0
    for Z in (ZA, ZB):
        assert not H.lin_vanishes_for_some_u(Z) and Z.subs(U, 1) > 0
    rho = sp.Rational(c[0], a[0]) + sp.Rational(c[0], b[0])
    assert sp.simplify(I["a_p"] / I["a_c"] - rho) == 0
    r = sp.cancel(I["a_d"] / I["a_c"])
    if region == "eq1":
        assert r == 1
        return
    if region == "eq2":
        assert r == 1 + rho
        return
    num = sp.numer(sp.together(sp.diff(r, U)))
    assert num != 0 and not H.positive_roots(num)
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
    """Rejects: the typo E1((a_p + a_d)R) as a +a_d atom (a different function at every u > 0 for all six: a_p !=
    a_c), a missing or extra E-atom at either shift, Ei where ln R belongs at a coincidence, or f's exponential
    rates other than {-a_d, -a_c, -(a_p+a_c)} (merged at a coincidence)."""
    labs = labels_of(case, region)
    I = info(case, region)
    E_minus, E_plus, rates = expected_labels(case, region)
    got_E = {lb for lb in labs if lb[1] != "1"}
    assert {lb for lb in got_E if lb[0] == sp.expand(-I["a_d"])} == E_minus
    assert {lb for lb in got_E if lb[0] == I["a_d"]} == E_plus
    assert got_E == E_minus | E_plus
    assert {lb[0] for lb in labs if lb[1] == "1"} == rates
    assert not H.lin_vanishes_for_some_u(I["a_p"] - I["a_c"])
    assert e_label(I["a_d"], -(I["a_p"] + I["a_d"])) not in got_E
    if region in ("eq1", "eq2"):
        assert any(lb[1] == "log" for lb in got_E)
        assert any("GAM" in ks for ks in entry(case, region)["form"])     # gamma IS in f at a coincidence


@pytest.mark.parametrize("case,region", ALL_PARAMS)
def test_no_two_labels_meet_at_any_u(case, region):
    """Rejects: an exceptional ratio inside an open region -- two rates (of f, of the d/dR closure), two E-labels at
    one shift, or an Ei argument meeting / vanishing at some u > 0."""
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
    """Rejects: order 4 at a generic ratio or 5 at a coincidence; d_E != 2 (e.g. 4: the two E-atoms at one shift
    treated as independent lines, order 7); L failing to kill the E-part; L f's surviving rates other than
    {-a_d, -a_c, -(a_p + a_c)} (merged at a coincidence)."""
    res = analysis(case, region)
    _Em, _Ep, rates = expected_labels(case, region)
    assert res["d_E"] == D_E_EXPECTED and all(v == 1 for v in res["per_shift"].values()), res["per_shift"]
    assert len(res["per_shift"]) == 2
    assert set(res["Lf_rates"]) == rates
    assert res["k"] == K_EXPECTED[region] == len(rates)
    assert res["order"] == ORDER[region]


@pytest.mark.parametrize("case,region", ANALYSIS_PARAMS)
def test_e_part_shape_and_epsilon(case, region):
    """Rejects: a wrong eps, E-coefficients at one shift not exactly opposite (+-q), the +a_d pair not eps times the
    R -> -R mirror of the -a_d pair (+- eps q(-R)), or a q carrying gamma / a log.  This is the only test that
    sees a non-opposite pair: q and -2q still span one line, so d_E and the order are unchanged by it."""
    f = form(case, region)
    I = info(case, region)
    ap, ac, ad = I["a_p"], I["a_c"], I["a_d"]
    co = f.coords
    q = co[e_label(-ad, ad - ap - ac)]
    assert f.const_free(q)
    assert co[e_label(-ad, ad - ac)] == -q
    q_ref = f.K.from_expr(q.as_expr().subs(H.RS, -H.RS))
    eps = EPS[case]
    assert co[e_label(ad, -(ad + ap + ac))] == -eps * q_ref
    assert co[e_label(ad, -(ad + ac))] == eps * q_ref


@pytest.mark.parametrize("case,region", ANALYSIS_PARAMS)
def test_q_and_Lf_nonvanishing_and_constant_free_for_every_u(case, region):
    """Rejects: a gamma or log constant in q or L f (then 'nonzero' would rest on the arithmetic of the
    constants); a zero-gcd over Q[u] of q or an L f coefficient that is a polynomial in u (a root u > 0 would be an
    exceptional ratio) or zero; a u-factor of any denominator (f, q, L, L f) with a root u > 0.  The constants ARE
    present in f (logs; gamma at a coincidence), so 'constant-free' is a real filter."""
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


# ================================================================ V. frozen form = live evaluator / norm product
# One fresh point per case: u for an open region, a scale s for a coincidence.  Not used by the build memo
# (u = 1, 1/3 per open region, the single point of each coincidence) nor by the existing guard's own list of used
# points.
FRESH = {"H1": ("mid", "5/13"), "H3": ("lt", "19/6"), "H5": ("gt", "7/17"), "H6": ("lt", "5/13"),
         "H8": ("mid", "19/6"), "H10": ("gt", "7/17")}
R_LIVE = (sp.Rational(3, 7), sp.Rational(5, 2))
DPS = 90


def _frac(x) -> Fraction:
    x = sp.Rational(x)
    return Fraction(int(x.p), int(x.q))


def norm_product(case, ZA, ZB):
    out = sp.Integer(1)
    for o, Z in zip((ORB[x] for x in CASES[case]), (ZA, ZA, ZA, ZB)):
        out *= radial_norm(_frac(Z), o[0], o[1])
    return out


def test_fresh_points_are_fresh():
    """Rejects: a 'fresh' point the build memo (u = 1, 1/3) or the existing guard's used list already used, or a
    case without a live point."""
    used = {sp.Rational(x) for x in ("1/3", "1/2", "1", "2", "3", "5", "1/1000000", "3/7", "1000", "1/100000",
                                      "11/3", "1/10000", "9/4", "300", "2/7", "50")}
    assert set(FRESH) == set(CASES)
    for case, (region, x) in FRESH.items():
        assert region in REGIONS and sp.Rational(x) not in used, (case, x)


@pytest.mark.slow
@pytest.mark.parametrize("case", list(FRESH))
def test_frozen_form_is_live_evaluator_over_norm_product(case):
    """Rejects: a frozen form that is not the integral -- including one off by a u-dependent factor, whose ratio to
    hybrid_closed_form is constant in R and would pass a 'ratio constant in R' check.  The ratio must equal the
    EXACT product of the four radial_norm values (geovac, current working copy), at two R, to 1e-40, with 90
    working digits (40 were too few for H5 and H10 in the build memo)."""
    region, x = FRESH[case]
    I = info(case, region)
    orbs = [ORB[o] for o in CASES[case]]
    x = sp.Rational(x)
    if region in OPEN:
        u0, s = x, sp.Integer(1)
    else:
        assert not I["ZA"].has(U)
        u0, s = sp.Integer(1), x
    ZA, ZB = sp.Rational(I["ZA"].subs(U, u0)), sp.Rational(I["ZB"].subs(U, u0))
    live = hybrid_closed_form(_frac(s * ZA), *orbs[:3], _frac(s * ZB), orbs[3])
    N = norm_product(case, ZA, ZB)
    with mp.workdps(DPS):
        Nv = mp.mpf(str(sp.N(N, DPS + 5)))
        for Rv in R_LIVE:
            lv = mp.mpf(str(sp.N(live.subs(R_s, Rv), DPS + 5)))
            fv = mp.mpf(int(s.p)) / int(s.q) * Nv * H.form_value(entry(case, region), u0, s * Rv, dps=DPS + 5)
            assert abs(lv - fv) <= mp.mpf(10) ** -40 * abs(lv), (case, region, x, Rv, lv, fv)


# ================================================================ FIRE TESTS (2026-11-19)
# Plants applied in memory only: a pytest plugin in the session scratchpad (-p plantplug, env PLANT) patched, after
# collection, this module's expectations (T), the helper tests/_hybrid_ratio.py's functions (H), or a deep copy of
# this module's DATA (D); no file on disk was edited.  The -k selection was run per plant, then clean runs repeated.
#
#   plant (wrong answer)                                          -k selector                       result
#   ORDER open regions = 4 (T)                                    quotient_count                    FIRED 5/5 fast
#     same plant                                                  quotient_count and (H6-lt or      FIRED 2/2 --slow
#                                                                   H10-lt)
#   ORDER coincidences = 5 (T)                                    quotient_count                    FIRED 12/12
#   d_E = number of E-atoms (span_basis keeps every q) (H)        quotient_count or kernel          FIRED 18/18 (17 forms
#                                                                                                     + sentinel)
#   k counted from f's own exponential rates (H)                  kernel or quotient_count          FIRED 1: the sentinel
#                                                                                                     only; the 17 forms
#                                                                                                     NOT FIRED (expected:
#                                                                                                     every rate of f
#                                                                                                     survives into L f)
#   EPS[H5] = -1 (T)                                              e_part                            FIRED 3/3 (gt, eq1, eq2)
#   EPS[H1] = +1 (T)                                              e_part                            FIRED 3/3
#   one -a_d E-coefficient x 2, H1 gt (D)                         H1-gt                             FIRED: e_part only
#                                                                                                     (order unchanged:
#                                                                                                     expected, see docstr.)
#   +a_d pair given q(R) instead of q(-R), H3 mid (D)             H3-mid                            FIRED: e_part
#   gamma on both -a_d E-coefficients, H3 mid (D)                 H3-mid                            FIRED: e_part +
#                                                                                                     constant-free
#   gamma R e^{-a_c R} added, H1 gt (D)                           H1-gt                             FIRED: constant-free
#                                                                                                     (Lf[-1]); order NOT
#                                                                                                     FIRED (expected)
#   whole H6 gt form x (2u - 1): gcd root u = 1/2 (D)             nonvanishing and H6-gt            FIRED
#   whole H8 gt form / (u - 3): denominator root u = 3 (D)        nonvanishing and H8-gt            FIRED
#   H8 lt rate -1 relabelled -3u + 1 (meets -u - 1 at u = 1) (D)  no_two_labels and H8-lt           FIRED
#   H1 gt map a_d = u + 4: covers (4, oo), not (3, oo) (D)        region_map and H1-gt              FIRED (endpoints)
#   H1 gt map a_d = (u - 1)^2 + 3: not monotone (D)               region_map and H1-gt              FIRED (monotonicity)
#   H5 orbitals 2p1 -> 2p0 (D)                                    data_has                          FIRED
#   H10 eq2 region dropped (D)                                    data_has                          FIRED
#   each FRESH form x 1/(u + 2) (R-constant ratio) (D)            live --slow                       FIRED 6/6
#   one exponential coefficient x (1 + 1e-20), H10 gt (D)         live and H10 --slow               FIRED
#
# Clean runs: fast 144 passed, 45 skipped (~31 s); --slow 189 passed (402-492 s over two runs).
