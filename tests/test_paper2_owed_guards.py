"""Paper 2 owed guards (group5 guard pass, 2026-11-07).

Backs papers/group5_qed_gauge/Paper_2_Alpha.tex ("The Fine Structure Constant
as a Spectral Coincidence on the Three-Sphere").  The combination rule
K = pi(B + F - Delta) is an Observation; these tests check arithmetic and
counts only.  Nothing here bears on why the sum lies near 1/alpha.

Each test names its paper locus and the wrong answer it rejects.  Expected
values come from a route independent of the quantity under test: exact
arithmetic over representation dimensions, a closed form obtained separately,
or a cited literature value (CODATA).

The look-elsewhere search (sec:statistical) is re-implemented here from the
paper's description of the search space (term pool, 14 prefactors, linear and
cubic readings, unordered triples with replacement).  The pool is built in
exact rational arithmetic and the enumeration is vectorized with a window
search, a different algorithm from the preserved driver
tests/paper2_support/phase2_combinatorial.py (which is not imported).  Only
formulas near the target need evaluating, so the whole 1.92e9-formula space is
counted exactly in about a second; no reduced slice is needed.

Not repeated here (already backed):
  Dirac/Weyl cumulative sums and (7/2)*12 = 42 .... tests/test_paper2_dirac_cumulative.py
  chiral (Hopf) grouping gives B/2 = 21 ........... tests/test_paper25_hopf_grouping.py
  g_3^D = 40 ...................................... tests/test_dirac_s3.py
  Dirac heat trace (sqrt(pi)/2) t^-3/2 - (sqrt(pi)/4) t^-1/2 (row 12)
                                                    tests/test_paper18_computed_pins.py
  B_formal(m)/N(m) = d, symbolic in d ............. tests/test_paper2_corrections.py
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from math import comb
from typing import Dict, List, Tuple

import numpy as np
import pytest
import sympy as sp
from mpmath import mp, mpf, nsum, inf, polyroots

# Cited literature values of 1/alpha.
#   CODATA 2018: Tiesinga, Mohr, Newell, Taylor, Rev. Mod. Phys. 93, 025010 (2021).
#   CODATA 2022: Mohr, Newell, Taylor, Tiesinga, Rev. Mod. Phys. 97, 025002 (2025).
ALPHA_INV_2018 = "137.035999084"
ALPHA_INV_2022 = "137.035999177"


@pytest.fixture(autouse=True)
def _mp_dps_hermetic():
    saved = mp.dps
    mp.dps = 50
    yield
    mp.dps = saved


# ---------------------------------------------------------------------------
# Representation-theoretic building blocks (exact integers)
# ---------------------------------------------------------------------------

def sphere_degeneracy(k: int, D: int) -> int:
    """Multiplicity of the degree-k harmonics on S^D (D >= 1)."""
    if D == 1:
        return 1 if k == 0 else 2
    return comb(k + D, D) - (comb(k + D - 2, D) if k >= 2 else 0)


def scalar_s3_eigenvalue(n: int) -> int:
    """|lambda| of the n-th scalar S^3 shell, from the S^D formula k(k+D-1), k = n-1, D = 3."""
    k = n - 1
    return k * (k + 2)


def dirac_degeneracy_ch(k: int) -> int:
    """Camporesi-Higuchi multiplicity of level k: (j1, j2) = ((k+1)/2, k/2) plus its swap."""
    j1, j2 = Fraction(k + 1, 2), Fraction(k, 2)
    return int(2 * (2 * j1 + 1) * (2 * j2 + 1))


def B_shells(n_max: int) -> int:
    """Eq. (eq:casimir): sum over shells n <= n_max of (2l+1) l(l+1), l < n."""
    return sum((2 * l + 1) * l * (l + 1) for n in range(1, n_max + 1) for l in range(n))


def K_value() -> mpf:
    return mp.pi * (42 + mp.pi ** 2 / 6 - mpf(1) / 40)


def cubic_small_root(K: mpf) -> mpf:
    """Smallest positive root of a^3 - K a + 1 = 0, by the contraction a = (1 + a^3)/K."""
    a = 1 / K
    for _ in range(60):
        a = (1 + a ** 3) / K
    return a


# ---------------------------------------------------------------------------
# sec. "Delta: a reciprocal degeneracy"
# ---------------------------------------------------------------------------

def test_delta_two_forms_agree_only_at_cutoff_three():
    """Paper 2 eqs. (eq:boundary), (eq:Delta_dirac) and the sentence after them: as functions
    of the cutoff m, |lambda_m| N(m-1) and g_m^D agree among positive integers only at m = 3.

    Rejects: the two forms being one function of the cutoff (agreement at every m), and an
    agreement at any second cutoff (m = 1..8 checked directly, all m by the polynomial).
    Route: |lambda| and N from S^D harmonic counts, g^D from SU(2) x SU(2) dimensions; the
    'all m' clause from the factored difference polynomial.
    """
    first = {m: scalar_s3_eigenvalue(m) * sum(sphere_degeneracy(k, 3) for k in range(m - 1))
             for m in range(1, 9)}
    second = {m: dirac_degeneracy_ch(m) for m in range(1, 9)}
    assert first[3] == second[3] == 40
    assert [m for m in range(1, 9) if first[m] == second[m]] == [3]
    m = sp.Symbol("m")
    diff = sp.expand((m ** 2 - 1) * (m - 1) * m * (2 * m - 1) / 6 - 2 * (m + 1) * (m + 2))
    assert diff != 0
    for mm in range(1, 9):  # the polynomials are the same functions as the counts above
        assert diff.subs(m, mm) == first[mm] - second[mm]
    # rational root theorem: every rational root is found by sympy's factorization over Q
    rational_roots = [r for r in sp.roots(sp.Poly(diff, m), filter="Q")]
    assert [r for r in rational_roots if r.is_integer and r > 0] == [3]


# ---------------------------------------------------------------------------
# sec. "S^3 specificity of the criterion"
# ---------------------------------------------------------------------------

def sub_casimir_trace(k: int, d: int) -> Fraction:
    """Trace of the SO(d) Casimir on the degree-k harmonics of S^d, by the branching
    SO(d+1) -> SO(d): degree k restricts to degrees j = 0..k on S^(d-1), Casimir j(j+d-2)."""
    return Fraction(sum(sphere_degeneracy(j, d - 1) * j * (j + d - 2) for j in range(k + 1)))


def formal_half_trace(k: int, d: int) -> Fraction:
    """(1/2) g_k |lambda_k| on S^d, |lambda_k| = k(k+d-1): the summand of B_formal."""
    return Fraction(sphere_degeneracy(k, d) * k * (k + d - 1), 2)


def test_B_equals_B_formal_only_on_S3():
    """Paper 2 sec. 'S^3 specificity of the criterion': B_formal(2)/N(2) = d on every S^d,
    while the sub-Casimir trace B equals B_formal only for d = 3.

    Rejects: the universal reading (B = B_formal on every S^d, so B = 42 would be an S^d-generic
    number), and a B_formal(2)/N(2) that is not d.  Route: the sub-Casimir trace is summed over
    the SO(d+1) -> SO(d) branching, independent of the eigenvalue formula; the expected ratio
    2(d-1)/(d+1) is the generator count dim SO(d)/dim SO(d+1) times 2 (equal trace per generator
    on an irreducible representation), obtained separately.
    """
    equal_d = []
    for d in range(2, 11):
        for k in range(1, 7):
            ratio = sub_casimir_trace(k, d) / formal_half_trace(k, d)
            assert ratio == Fraction(2 * (d - 1), d + 1), (d, k)
        if all(sub_casimir_trace(k, d) == formal_half_trace(k, d) for k in range(7)):
            equal_d.append(d)
    assert equal_d == [3]
    for d in range(1, 11):
        Bf = sum(formal_half_trace(k, d) for k in range(3))
        N2 = sum(sphere_degeneracy(k, d) for k in range(3))
        assert Bf / N2 == d, d
    # on S^3 the k <= 2 sub-Casimir trace is the shell sum B(3) = 42
    assert sum(sub_casimir_trace(k, 3) for k in range(3)) == B_shells(3) == 42


# ---------------------------------------------------------------------------
# sec. "B: a Casimir trace over the scalar S^3 harmonics"
# ---------------------------------------------------------------------------

S2_TRACE_L_LE_2 = 36


def test_S2_trace_is_36_not_B():
    """Paper 2 sec. ingredients, 'B is an S^3 quantity, not an S^2 one': the S^2 Laplacian trace
    over l <= 2 is 0 + 6 + 30 = 36, and B = 42 exceeds it by one extra l = 1 multiplet (l = 1
    occurs in shells n = 2 and n = 3).

    Rejects: the retired 'S^2 base' attribution of B (S^2 trace = 42).  Route: the S^2 trace uses
    the S^D multiplicity formula and k(k+1); B is the (n, l) shell sum of eq:casimir.
    """
    s2 = sum(sphere_degeneracy(l, 2) * l * (l + 1) for l in range(3))
    assert s2 == S2_TRACE_L_LE_2
    assert B_shells(3) - s2 == 6 == sphere_degeneracy(1, 2) * 1 * 2
    assert sum(1 for n in range(1, 4) for l in range(n) if l == 1) == 2
    assert s2 != B_shells(3)


# ---------------------------------------------------------------------------
# sec:formula
# ---------------------------------------------------------------------------

def test_K_and_cubic_root_relative_differences():
    """Paper 2 abstract, sec:formula, conclusion: K differs from CODATA 2018 1/alpha by 4.77e-7
    relative and the cubic's smallest positive root by 8.8e-8 (8.81e-8 in sec:statistical); with
    CODATA 2022 the two are 4.76e-7 and 8.75e-8; both lie above the measured value, the root
    closer by a factor of about five.

    Rejects: swapping the two differences (K at 8.8e-8, root at 4.8e-7), swapping the CODATA
    editions, and a sign error in Delta (pi(B + F + Delta) is 1.1e-3 off).  Route: K at 50
    digits; the root by fixed-point contraction, cross-checked against polyroots; 1/alpha cited.
    """
    K = K_value()
    a = cubic_small_root(K)
    assert abs(a ** 3 - K * a + 1) < mpf(10) ** -45
    roots = [r.real for r in polyroots([1, 0, -K, 1], maxsteps=200, extraprec=200)]
    assert abs(min(r for r in roots if r > 0) - a) < mpf(10) ** -40
    assert mp.nstr(K, 9) == "137.036064"
    assert mp.nstr(1 / a, 9) == "137.036011"
    out = {}
    for tag, lit in (("2018", ALPHA_INV_2018), ("2022", ALPHA_INV_2022)):
        c = mpf(lit)
        out[tag] = ((K - c) / c, (1 / a - c) / c)
        assert out[tag][0] > 0 and out[tag][1] > 0
    assert mp.nstr(out["2018"][0], 3) == "4.77e-7"
    assert mp.nstr(out["2018"][1], 3) == "8.81e-8"
    assert mp.nstr(out["2018"][0], 2) == "4.8e-7"
    assert mp.nstr(out["2018"][1], 2) == "8.8e-8"
    assert mp.nstr(out["2022"][0], 3) == "4.76e-7"
    assert mp.nstr(out["2022"][1], 3) == "8.75e-8"
    assert 4.5 < out["2018"][0] / out["2018"][1] < 6


# ---------------------------------------------------------------------------
# sec:ledger, the sign check
# ---------------------------------------------------------------------------

def test_eight_sign_patterns():
    """Paper 2 sec:ledger (Table tab:ledger 'Signs' row and the paragraph after it): of the eight
    patterns pi(+-B +-F +-Delta), (+,+,-) gives 4.8e-7, the next best (+,+,+) 1.1e-3, the rest
    7.4e-2 or worse; flipping Delta moves K by 2 pi Delta = 1.1e-3 relative; pi(42 + pi^2/6) alone
    lies within 5.8e-4.

    Rejects: a different best pattern, a runner-up other than (+,+,+), a remaining pattern
    closer than 7.4e-2.  Route: direct evaluation of all eight against cited CODATA 2018.
    """
    c = mpf(ALPHA_INV_2018)
    B, F, D = mpf(42), mp.pi ** 2 / 6, mpf(1) / 40
    err = {s: abs(mp.pi * (s[0] * B + s[1] * F + s[2] * D) - c) / c
           for s in itertools.product((1, -1), repeat=3)}
    ranked = sorted(err, key=lambda s: err[s])
    assert ranked[0] == (1, 1, -1) and f"{float(err[ranked[0]]):.1e}" == "4.8e-07"
    assert ranked[1] == (1, 1, 1) and f"{float(err[ranked[1]]):.1e}" == "1.1e-03"
    assert all(err[s] >= mpf("7.4e-2") for s in ranked[2:])
    assert f"{float(min(err[s] for s in ranked[2:])):.1e}" == "7.4e-02"
    assert f"{float(2 * mp.pi * D / c):.1e}" == "1.1e-03"
    assert abs(mp.pi * (B + F) - c) / c < mpf("5.8e-4")


# ---------------------------------------------------------------------------
# sec:statistical, the look-elsewhere count
# ---------------------------------------------------------------------------

def build_pool() -> Tuple[np.ndarray, int]:
    """The paper's term pool: S^3 eigenvalues and degeneracies (n <= 4), S^2 Casimirs and
    degeneracies (l <= 3), cumulative counts N(n) (n <= 4), small integers 1..4, the two Casimir
    sums 42 and 10, pairwise products up to 50, signed ratios of these, and +- the six
    transcendentals.  Returns (sorted float terms, number of distinct rationals)."""
    base = set()
    for n in range(1, 5):
        base |= {n * n - 1, n * n, sum(i * i for i in range(1, n + 1))}
    for l in range(4):
        base |= {l * (l + 1), 2 * l + 1}
    base |= {1, 2, 3, 4, 42, 10}
    nz = sorted(x for x in base if x > 0)
    ints = base | {p * q for p in nz for q in nz if p * q <= 50}
    rats = {sgn * Fraction(p, q) for p in ints for q in ints if q > 0 for sgn in (1, -1)}
    tr = [mp.pi, mp.pi ** 2, mp.zeta(2), mp.zeta(3), mp.zeta(4), mp.log(2)]
    vals = [float(r) for r in rats] + [float(s * t) for t in tr for s in (1, -1)]
    return np.array(sorted(vals)), len(rats)


def prefactors() -> List[float]:
    p = np.pi
    return [p, 2 * p, p ** 2, p / 2, 4 * p, 4 * p ** 2, 1 / (2 * p), 1.0, 2.0, 0.5,
            p / 4, p / 6, 1 / p, np.sqrt(p)]


SEARCH_THRESHOLDS = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 8.8e-8]
PRINTED_HITS = [1_283_592, 125_027, 14_476, 682, 37, 11, 10]
PRINTED_N = 1_921_873_520
PRINTED_FRACTIONS = ["6.7e-04", "6.5e-05", "7.5e-06", "3.5e-07", "1.9e-08", "5.7e-09", "5.2e-09"]


def enumerate_near_hits(window: float = 0.02) -> Dict[str, np.ndarray]:
    """Every formula of the search within `window` (relative) of 1/alpha_2018, by a window
    search over sorted terms.  Returns arrays: err, mode (0 linear, 1 cubic), prefactor index,
    and the triple (i <= j <= k) of term indices."""
    terms, _ = build_pool()
    c = float(ALPHA_INV_2018)
    K_target = c + 1 / c ** 2           # K whose cubic root is exactly 1/c
    I, J = np.triu_indices(len(terms))
    S2 = terms[I] + terms[J]
    out = {key: [] for key in ("err", "mode", "pf", "i", "j", "k")}
    for t_idx, T in enumerate(prefactors()):
        for mode in (0, 1):
            centre = (c if mode == 0 else K_target) / T
            tol = window * c / T
            need = centre - S2
            lo = np.maximum(np.searchsorted(terms, need - tol, "left"), J)
            hi = np.searchsorted(terms, need + tol, "right")
            cnt = np.clip(hi - lo, 0, None)
            pair = np.repeat(np.arange(len(S2)), cnt)
            k = np.repeat(lo, cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))
            val = T * (S2[pair] + terms[k])
            if mode == 1:
                a = 1 / val
                for _ in range(8):
                    a = (1 + a ** 3) / val
                val = 1 / a
            out["err"].append(np.abs(val - c) / c)
            out["mode"].append(np.full(len(k), mode))
            out["pf"].append(np.full(len(k), t_idx))
            out["i"].append(I[pair])
            out["j"].append(J[pair])
            out["k"].append(k)
    return {key: np.concatenate(v) for key, v in out.items()}


@pytest.fixture(scope="module")
def search():
    return enumerate_near_hits()


def test_search_space_size():
    """Paper 2 sec:statistical: 743 distinct term values, 14 prefactors, linear and cubic
    readings, N = 1.92e9 formulas (1,921,873,520 = 2 x 14 x C(745, 3), unordered triples with
    replacement).

    Rejects: ordered triples (N = 2 x 14 x 743^3), a single reading (N/2), or a pool size other
    than 743.  Route: the rational part of the pool is deduplicated exactly (Fraction), and its
    float image must not merge any two values.
    """
    terms, n_rat = build_pool()
    assert n_rat + 12 == len(terms) == 743
    assert len(np.unique(np.round(terms, 12))) == 743
    n_formulas = 2 * len(prefactors()) * comb(len(terms) + 2, 3)
    assert n_formulas == PRINTED_N
    assert f"{n_formulas:.2e}" == "1.92e+09"


def test_search_hit_table(search):
    """Paper 2 Table tab:hits: the seven hit counts at relative thresholds 1e-2 .. 8.8e-8.

    Rejects: a count from one reading only (e.g. 4 linear hits at 8.8e-8), counting distinct
    values instead of formulas (3 at 8.8e-8), and a pool or prefactor list that differs from the
    paper's.  Route: an independent vectorized enumeration with a 2% window, wider than any
    threshold, so no hit is lost to the window.
    """
    err = search["err"]
    assert [int((err < t).sum()) for t in SEARCH_THRESHOLDS] == PRINTED_HITS


def test_search_rank_and_chance_count(search):
    """Paper 2 sec:statistical results: ten formulas match more closely than the cubic's root
    (8.81e-8), so the formula ranks eleventh; the ten take three distinct values, the closest
    4 pi (5/49 + 14/15 + pi^2) at 2.1e-8 (largest 7.6e-8); p = 10/N = 5.2e-9, and the
    Fraction column of Table tab:hits (hits / N at each threshold); at K's own precision
    (4.8e-7, linear reading) 27 entries match more closely, one of them this formula's cubic
    reading.

    Rejects: rank 10 or 12 (counting the formula itself or ties), a Fraction column computed
    against another N (ordered triples, or one reading only: every entry doubles), a closest
    competitor other than the printed one, and 26 or 28 entries beating K.  Route: the formula's
    own errors are located in the enumeration by its term indices and checked against an
    independent 50-digit evaluation.
    """
    terms, _ = build_pool()
    err, mode, pf = search["err"], search["mode"], search["pf"]
    tri = np.stack([search["i"], search["j"], search["k"]], axis=1)

    def locate(target_vals, m, t_idx):
        want = sorted(int(np.argmin(np.abs(terms - v))) for v in target_vals)
        hit = np.where((mode == m) & (pf == t_idx) & np.all(tri == want, axis=1))[0]
        assert len(hit) == 1
        return float(err[hit[0]])

    ingredients = [42.0, float(mp.zeta(2)), -1 / 40]
    e_cubic = locate(ingredients, 1, 0)
    e_lin = locate(ingredients, 0, 0)
    K = K_value()
    c = mpf(ALPHA_INV_2018)
    assert abs(e_cubic - float((1 / cubic_small_root(K) - c) / c)) < 1e-12
    assert abs(e_lin - float((K - c) / c)) < 1e-12

    better = err < e_cubic
    n_better = int(better.sum())
    assert n_better == 10
    assert n_better + 1 == 11                     # the rank
    vals = np.sort(err[better])
    distinct = vals[np.concatenate(([True], np.diff(vals) > 1e-15))]
    assert len(distinct) == 3
    assert f"{distinct.min():.1e}" == "2.1e-08" and f"{distinct.max():.1e}" == "7.6e-08"
    closest = mp.mpf(4) * mp.pi * (mpf(5) / 49 + mpf(14) / 15 + mp.pi ** 2)
    assert abs(float(abs(closest - c) / c) - distinct.min()) < 1e-12
    p = n_better / PRINTED_N
    assert f"{p:.1e}" == "5.2e-09"
    # Table tab:hits "Fraction" column: enumerated hits / N, as printed.
    fracs = [f"{int((err < t).sum()) / PRINTED_N:.1e}" for t in SEARCH_THRESHOLDS]
    assert fracs == PRINTED_FRACTIONS
    beat_K = err < e_lin
    assert int(beat_K.sum()) == 27
    assert e_cubic < e_lin                         # the cubic reading is among the 27


# ---------------------------------------------------------------------------
# sec:mechanisms, Table tab:mechanisms rows that are checkable identities
# ---------------------------------------------------------------------------

def test_casimir_power_traces_row2():
    """Paper 2 Table tab:mechanisms row 2: sum_{n<=3} n^2 C_4(n)^k, C_4(n) = (n^2 - 1)/2, gives
    a component of K (42 = B) only at k = 1.

    Rejects: a higher Casimir power reproducing B or 1/Delta = 40, and the normalization
    C_4 = n^2 - 1 (whose k = 1 trace is 84).  Route: B from the (n, l) shell sum, C_4 from the
    S^3 eigenvalue formula.
    """
    traces = {k: sum(Fraction(n * n) * Fraction(scalar_s3_eigenvalue(n), 2) ** k
                     for n in range(1, 4)) for k in range(0, 9)}
    assert [k for k, v in traces.items() if v == B_shells(3)] == [1]
    assert all(v not in (40, Fraction(1, 40)) for v in traces.values())
    assert sum(n * n * scalar_s3_eigenvalue(n) for n in range(1, 4)) == 2 * B_shells(3)


def test_sm_charge_trace_row7():
    """Paper 2 Table tab:mechanisms row 7: sum N_c Q_f^2 = 8 = |lambda_3| (three generations;
    8/3 in every generation, the charge-universal piece).

    Rejects: a per-generation sum of 8 (total 24), and the sum without the colour factor
    (N_c = 1: 14/9 per generation, total 14/3).  Route:
    Standard Model charges, Q in units of e: up-type 2/3, down-type -1/3 (N_c = 3), charged
    lepton -1, neutrino 0.
    """
    gen = [(3, Fraction(2, 3)), (3, Fraction(-1, 3)), (1, Fraction(-1)), (1, Fraction(0))]
    per_gen = [sum(nc * q * q for nc, q in gen) for _ in range(3)]
    assert per_gen == [Fraction(8, 3)] * 3
    assert sum(per_gen) == 8 == scalar_s3_eigenvalue(3)


def test_dirac_dirichlet_series_obstruction2():
    """Paper 2 sec:obstructions, Obstruction 2, eq. (eq:D_dirac_zeta3): sum_{m>=1} 2m(m+1)/m^4
    = 2 zeta(2) + 2 zeta(3); in the Camporesi-Higuchi index the series is
    sum_{k>=0} 2(k+1)(k+2)/(k+3/2)^4 = pi^2 - pi^4/12.

    Rejects: a Dirac series with no zeta(3) piece (2 zeta(2), which would lift F), and the CH
    form with a rational remainder (pi^2 - pi^4/12 - 8, the Hurwitz offsets left uncancelled).
    Route: direct summation at 50 digits against the closed forms; the CH closed form also via
    the Hurwitz identity zeta(s, 3/2) = (2^s - 1) zeta(s) - 2^s.
    """
    m = sp.Symbol("m", positive=True, integer=True)
    sym = sp.summation(2 * m * (m + 1) / m ** 4, (m, 1, sp.oo))
    assert sp.simplify(sym - (2 * sp.zeta(2) + 2 * sp.zeta(3))) == 0
    direct = nsum(lambda n: 2 * n * (n + 1) / n ** 4, [1, inf])
    assert abs(direct - (2 * mp.zeta(2) + 2 * mp.zeta(3))) < mpf(10) ** -30
    assert abs(direct - 2 * mp.zeta(2)) > mpf("2")
    ch = nsum(lambda k: 2 * (k + 1) * (k + 2) / (k + mpf(3) / 2) ** 4, [0, inf])
    closed = mp.pi ** 2 - mp.pi ** 4 / 12
    assert abs(ch - closed) < mpf(10) ** -30

    def hz(s):
        return (2 ** s - 1) * mp.zeta(s) - 2 ** s
    assert abs((2 * hz(2) - hz(4) / 2) - closed) < mpf(10) ** -40


def test_hopf_charge_split_of_level_three_row8():
    """Paper 2 sec:obstructions, Obstruction 3: the 40 states of Dirac level k = 3 split 20 + 20
    by chirality, equivalently into integer and half-integer charge sectors of the U(1) in one
    SU(2) factor.

    Rejects: an unequal split, and charge sectors that cut across chirality.  Route: enumerate
    the states of ((k+1)/2, k/2) and (k/2, (k+1)/2) with their left-factor weight m_L.
    """
    k = 3
    states = []
    for chir, (j1, j2) in enumerate(((Fraction(k + 1, 2), Fraction(k, 2)),
                                     (Fraction(k, 2), Fraction(k + 1, 2)))):
        for a in range(int(2 * j1) + 1):
            for b in range(int(2 * j2) + 1):
                states.append((chir, -j1 + a))
    assert len(states) == dirac_degeneracy_ch(k) == 40
    by_chir = [sum(1 for ch, _ in states if ch == c) for c in (0, 1)]
    integer_charge = {ch for ch, mL in states if mL.denominator == 1}
    half_charge = {ch for ch, mL in states if mL.denominator == 2}
    assert by_chir == [20, 20]
    assert integer_charge == {0} and half_charge == {1}
