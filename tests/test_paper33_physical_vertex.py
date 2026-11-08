"""Paper 33 Sec. 2 (sec:physical): the physical vector-photon vertex on round S^3.

Backs, in exact arithmetic:
  * the content controls (Dirac spectrum, coclosed eigenforms)  -- Sec. 2.1, 2.3 Controls;
  * Theorem thm:physical_rule (same-sign and opposite-sign channels) and the
    open-cell census of Sec. 2.3 "Results";
  * Corollaries cor:model_comparison and cor:ground_state;
  * the two-route agreement (C4) and the 9j closed form eq:ninej (all 896 cells);
  * Proposition prop:diagonal (the diagonal vertex does not vanish), by both routes;
  * Proposition prop:ward (eq:ward, Ward-Takahashi identity);
  * m conservation in each SU(2) factor (Sec. 2.3 Results);
  * the unnormalized one-loop sum of Sec. 2.6 (0.2669, 0.2702, 0.2715, 0.2719 at
    n_max = 3, 5, 8, 12);
  * geovac.qed_vertex.so4_channel_count against the theorem on the model window.

Port of VV stage A (docs/predictions/2026-11-07_physical_vector_vertex.md).
Route 1 (group theory) and route 2 (R^4 coordinates; imports neither route 1
nor geovac) live in tests/paper33_vertex_support/.  The theorem's rule is
transcribed below from the paper's text, independently of both routes.
"""
from __future__ import annotations

import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Dict, Tuple

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "paper33_vertex_support"))
import p33_route1 as R1  # noqa: E402
import p33_route2 as R2  # noqa: E402

SIGNS = (1, -1)
COMPS = ("A", "B")
N1MAX, Q1MAX = 3, 7          # route 1 range (Sec. 2.3): 896 cells

# Sec. 2.6: unnormalized one-loop sum at n_ext = 0, n_max = 3 (exact rational of
# the stage-A data; the paper prints 0.2669).  Recomputed below by two routes.
SIGMA0_NMAX3 = Fraction(2628144752, 9845600625)
SIGMA0_PRINTED = {3: 0.2669, 5: 0.2702, 8: 0.2715, 12: 0.2719}


# ---------------------------------------------------------------------------
# Theorem thm:physical_rule, transcribed from the paper (not from either route)
# ---------------------------------------------------------------------------
def theorem_rule(s1: int, s2: int, comp: str, n1: int, n2: int, q: int) -> bool:
    if q < 1:
        return False
    d = abs(n1 - n2)
    if s1 == s2:
        if (n1 + n2 + q) % 2 != 1:
            return False
        c = comp if s1 > 0 else {"A": "B", "B": "A"}[comp]      # -- is the mirror A <-> B
        if c == "A":
            return d + 1 <= q <= n1 + n2 + 1
        return d + 1 <= q <= n1 + n2 - 1
    if (n1 + n2 + q) % 2 != 0:
        return False
    if s1 < 0:                                                  # -+ is the transpose of +-
        n1, n2 = n2, n1
    if comp == "A":
        return (abs(n1 + 1 - n2) <= q + 1 <= n1 + n2 + 1) and (abs(n1 - n2 - 1) <= q - 1 <= n1 + n2 + 1)
    return (abs(n1 + 1 - n2) <= q - 1 <= n1 + n2 + 1) and (abs(n1 - n2 - 1) <= q + 1 <= n1 + n2 + 1)


def model_rule(n1: int, n2: int, q: int) -> bool:
    """Paper 28 eq:vertex_rule (the vertex model)."""
    return q >= 1 and abs(n1 - n2) <= q <= n1 + n2 and (n1 + n2 + q) % 2 == 1


def old_mixed_count(n1: int, n2: int, q: int) -> int:
    """The retired count: + content for n1, - content for n2, triangle only, no integrality."""
    if not model_rule(n1, n2, q):
        return 0
    j1 = (Fraction(n1 + 1, 2), Fraction(n1, 2))
    j2 = (Fraction(n2, 2), Fraction(n2 + 1, 2))
    tri = lambda a, b, c: abs(a - b) <= c <= a + b  # noqa: E731
    return sum(1 for g in ((Fraction(q + 1, 2), Fraction(q - 1, 2)), (Fraction(q - 1, 2), Fraction(q + 1, 2)))
               if tri(j1[0], g[0], j2[0]) and tri(j1[1], g[1], j2[1]))


# ---------------------------------------------------------------------------
# cached computations
# ---------------------------------------------------------------------------
Cell = Tuple[int, int, str, int, int, int]   # (s1, s2, comp, n1, n2, q)


@lru_cache(maxsize=None)
def route1_raw() -> Dict[Cell, Dict]:
    """Route 1's full per-cell record (open, S, and m_conserved on open cells)."""
    out = {}
    for s1 in SIGNS:
        for s2 in SIGNS:
            for comp in COMPS:
                for n1 in range(N1MAX + 1):
                    for n2 in range(N1MAX + 1):
                        for q in range(1, Q1MAX + 1):
                            out[(s1, s2, comp, n1, n2, q)] = R1.vertex_route1(s1, n1, comp, q, s2, n2)
    return out


@lru_cache(maxsize=None)
def route1_table() -> Dict[Cell, Tuple[bool, Fraction]]:
    return {c: (v["open"], v["S"]) for c, v in route1_raw().items()}


@lru_cache(maxsize=None)
def route2(nmax: int, qmax: int) -> "R2.Route2":
    return R2.Route2(nmax, qmax)


def ninej_sum(s1: int, n1: int, comp: str, q: int, s2: int, n2: int) -> Fraction:
    """eq:ninej, evaluated with sympy's Wigner 9j (independent of route 1's CG sums)."""
    from sympy import Rational, nsimplify
    from sympy.physics.wigner import wigner_9j
    tJ1, tj1 = R1.spinor_content(s1, n1)
    tJA, tjA = R1.vector_content(comp, q)
    tJ2, tj2 = R1.spinor_content(s2, n2)

    def tri(a: int, b: int, c: int) -> bool:
        return abs(a - b) <= c <= a + b and (a + b + c) % 2 == 0

    if not (tri(tjA, tj2, tj1) and tri(tJA, tJ2, tJ1)):
        return Fraction(0)
    h = lambda x: Rational(x, 2)  # noqa: E731
    w = wigner_9j(h(tjA), h(tj2), h(tj1), 1, h(1), h(1), h(tJA), h(tJ2), h(tJ1), prec=None)
    val = Rational(6 * (tjA + 1) * (tj2 + 1) * (tj1 + 1) * (tJA + 1) * (tJ2 + 1) * (tJ1 + 1)) * Rational(nsimplify(w ** 2))
    return Fraction(int(val.p), int(val.q))


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------
def test_control_dirac_spectrum_route2() -> None:
    """Sec. 2.1 / 2.3 Controls: D has eigenvalues +-(n+3/2), (n+1)(n+2) states per sign (route 2, n <= 1).

    Rejects: the constant term +3/2 in D (spectrum then contains 1/2, not of the
    form +-(n+3/2), and is not symmetric), and any multiplicity other than (n+1)(n+2)
    per sign. Independent route: exact null spaces of D - lambda on polynomial
    spinors, then D applied again to every basis field.
    """
    sp = R2.spinor_spaces(1)
    assert set(sp) == {(s, n) for s in SIGNS for n in (0, 1)}
    for (s, n), fields in sp.items():
        assert len(fields) == (n + 1) * (n + 2)
        lam = Fraction(s * (2 * n + 3), 2)
        assert all(R2.check_spinor_eigen(f, R2.C0, lam) for f in fields)
    # the other constant is not a symmetric +-(n+3/2) spectrum
    with pytest.raises(AssertionError):
        R2.spinor_spaces(1, c0=-R2.C0)


def test_control_coclosed_eigenforms_route2() -> None:
    """Sec. 2.1 eq:vector_content: coclosed 1-forms at q have Hodge eigenvalue (q+1)^2, q(q+2) states per component (q <= 2).

    Rejects: Hodge eigenvalue q(q+2) (the photon weight mu_q of eq:self_energy),
    and a per-component count other than q(q+2). Independent route: polynomial
    eigenforms of *d in R^4; tangency, coclosedness and (*d)^2 checked afterwards.
    """
    vs = R2.vector_spaces(2)
    for (sd, q), fields in vs.items():
        assert len(fields) == q * (q + 2)
        for A in fields:
            assert not R2.x_dot(A)                       # tangent
            assert not R2.div(A)                         # coclosed
            r = R2.eigen_ratio(A, R2.star_d(R2.star_d(A)))
            assert r is not None and r.im == 0
            assert r.re == (q + 1) ** 2 and r.re != q * (q + 2)


# ---------------------------------------------------------------------------
# The rule (SA1, SA2) and the census
# ---------------------------------------------------------------------------
def test_theorem_physical_rule_route1() -> None:
    """thm:physical_rule + Sec. 2.3 Results: route 1's open set on n1, n2 <= 3, q <= 7 equals the theorem exactly.

    Rejects: the model rule (q <= n1+n2, n1+n2+q odd only: it closes q = n1+n2+1
    and every opposite-sign channel), and the old mixed-sign count with no
    integrality (differs at 9 model-window triples). Independent route: route 1's
    exact CG amplitudes vs the rule transcribed from the paper.
    """
    tab = route1_table()
    assert len(tab) == 896
    mism = [c for c, (o, _) in tab.items() if o != theorem_rule(*c)]
    assert mism == []
    counts = {(s1, s2): sum(o for (a, b, *_), (o, _) in tab.items() if (a, b) == (s1, s2))
              for s1 in SIGNS for s2 in SIGNS}
    assert counts == {(1, 1): 44, (-1, -1): 44, (1, -1): 40, (-1, 1): 40}
    # every open opposite-sign cell has n1+n2+q even; at least one exists
    opp = [c for c, (o, _) in tab.items() if o and c[0] != c[1]]
    assert opp and all((c[3] + c[4] + c[5]) % 2 == 0 for c in opp)
    assert tab[(1, -1, "B", 0, 1, 1)][0]
    # same-sign: 32 of 88 open cells outside the model window, all at q = n1+n2+1, one component each
    same_open = [c for c, (o, _) in tab.items() if o and c[0] == c[1]]
    outside = [c for c in same_open if not model_rule(c[3], c[4], c[5])]
    assert len(same_open) == 88 and len(outside) == 32
    assert all(c[5] == c[3] + c[4] + 1 for c in outside)
    for s in SIGNS:
        for n1 in range(N1MAX + 1):
            for n2 in range(N1MAX + 1):
                assert sum(tab[(s, s, c, n1, n2, n1 + n2 + 1)][0] for c in COMPS) == 1
    # the retired count is not the physical one inside the model window
    diff = [(n1, n2, q) for n1 in range(4) for n2 in range(4) for q in range(1, 8)
            if model_rule(n1, n2, q)
            and old_mixed_count(n1, n2, q) != sum(tab[(1, 1, c, n1, n2, q)][0] for c in COMPS)]
    assert len(diff) == 9


def test_corollary_model_comparison_route1() -> None:
    """cor:model_comparison: on the model window both same-sign components are open (W = 2); every opposite-sign channel is outside it.

    Rejects: a same-sign count W of 0 or 1 anywhere on the model window, and any
    opposite-sign channel inside the model window.
    """
    tab = route1_table()
    for n1 in range(N1MAX + 1):
        for n2 in range(N1MAX + 1):
            for q in range(1, Q1MAX + 1):
                for s in SIGNS:
                    w = sum(tab[(s, s, c, n1, n2, q)][0] for c in COMPS)
                    if model_rule(n1, n2, q):
                        assert w == 2
                    for c in COMPS:
                        if tab[(s, -s, c, n1, n2, q)][0]:
                            assert not model_rule(n1, n2, q)


def test_corollary_ground_state_route1() -> None:
    """cor:ground_state (SA3): for n_ext = 0 the open channels are exactly same-sign q = n_int+1 and opposite-sign q = n_int (n_int >= 1), one component each.

    Rejects: the ground-state zero (the model's empty set at n_ext = 0, Paper 28
    Theorem 4), and the ground state not coupling to itself at q = 1.
    Independent route: route 1 amplitudes vs the corollary's stated set.
    """
    tab = route1_table()
    for se in SIGNS:
        for si in SIGNS:
            for comp in COMPS:
                for ni in range(N1MAX + 1):
                    for q in range(1, Q1MAX + 1):
                        if se == si:
                            want = q == ni + 1 and comp == ("A" if se > 0 else "B")
                        else:
                            want = q == ni and ni >= 1 and comp == ("B" if se > 0 else "A")
                        assert tab[(se, si, comp, 0, ni, q)][0] == want
                        assert tab[(si, se, comp, ni, 0, q)][0] == want
    assert tab[(1, 1, "A", 0, 0, 1)][0] and tab[(-1, -1, "B", 0, 0, 1)][0]


# ---------------------------------------------------------------------------
# Two routes and the closed form
# ---------------------------------------------------------------------------
def _c4(nmax: int, qmax: int) -> int:
    tab = route1_table()
    r2 = route2(nmax, qmax)
    n_open = 0
    for s1 in SIGNS:
        for s2 in SIGNS:
            for comp in COMPS:
                for n1 in range(nmax + 1):
                    for n2 in range(nmax + 1):
                        for q in range(1, qmax + 1):
                            o2, S2 = r2.cell(s1, n1, comp, q, s2, n2)
                            o1, S1 = tab[(s1, s2, comp, n1, n2, q)]
                            assert (o1, S1) == (o2, S2), (s1, s2, comp, n1, n2, q, S1, S2)
                            n_open += o1
    return n_open


def test_two_routes_agree_fast() -> None:
    """Sec. 2.3 (C4): routes 1 and 2 agree on open/closed and on exact sum_m |V|^2 (n <= 1, q <= 3: 96 cells).

    Rejects: a construction error in either route (e.g. a sign flip in the frame
    contraction or in Clifford multiplication changes sum_m |V|^2 on some cell).
    Independent route: route 2 uses polynomial fields in R^4 and monomial moments,
    no representation theory.
    """
    assert _c4(1, 3) == 20


@pytest.mark.slow
def test_two_routes_agree_full() -> None:
    """Sec. 2.3 (C4) on the paper's full common range n <= 2, q <= 5 (360 cells, 70 open).

    Rejects: as test_two_routes_agree_fast, on the paper's stated range.
    """
    assert _c4(2, 5) == 70


def test_ninej_closed_form() -> None:
    """eq:ninej + Sec. 2.3: sum_m |V|^2 equals the 9j closed form on all 896 route-1 cells (n1, n2 <= 3, q <= 7, all sign pairs and components); 168 nonzero.

    Rejects: the closed form with a different prefactor (e.g. 3 or 12 instead of 6)
    or with the frame spins swapped, and a closed form that holds only on the
    360-cell route-2 range; route 1 brute force vs sympy's 9j.  About 1.5 s.
    """
    tab = route1_table()
    assert len(tab) == 896
    nonzero = 0
    for (s1, s2, comp, n1, n2, q), (_, S) in tab.items():
        assert ninej_sum(s1, n1, comp, q, s2, n2) == S, (s1, s2, comp, n1, n2, q)
        nonzero += S != 0
    assert nonzero == 168


# ---------------------------------------------------------------------------
# Diagonal vertex and Ward-Takahashi
# ---------------------------------------------------------------------------
def test_prop_diagonal_vertex_route2() -> None:
    """prop:diagonal (SA4): sum over q = 1 photon m of |<a|c(A)|a>|^2 = 1 for every n = 0 state, open component only.

    Rejects: V(a, a) = 0 (the mechanism hard-coded in geovac/vector_qed.py) and the
    ground-state zero. Independent route: route 2 exact monomial integration.
    """
    r2 = route2(1, 3)
    for s, open_c, closed_c in ((1, "A", "B"), (-1, "B", "A")):
        d = r2.diagonal_ground(s)
        assert len(d[open_c]) == 2
        assert all(x == 1 for x in d[open_c])
        assert all(x == 0 for x in d[closed_c])


def test_prop_diagonal_vertex_route1() -> None:
    """prop:diagonal (SA4) by route 1 ("Both routes give 1 for every basis state"): for each
    n = 0 state of either sign, sum over the q = 1 photon m of |<a|c(A)|a>|^2 = 1 on the open
    component (A for +, B for -) and 0 on the other.

    Rejects: V(a, a) = 0 (the mechanism hard-coded in geovac/vector_qed.py), the ground-state
    zero, and the open/closed components swapped between the two signs.  Independent of the
    route-2 test above: exact CG amplitudes in the left-invariant frame.
    """
    for s, open_c, closed_c in ((1, "A", "B"), (-1, "B", "A")):
        d = R1.diagonal_ground(s)
        per_state = d[open_c]["sum_mA_|V(a,A,a)|^2"]
        assert d[open_c]["left_open"] and len(per_state) == 2      # (n+1)(n+2) = 2 states
        assert all(Fraction(x) == 1 for x in per_state.values())
        closed = d[closed_c]
        assert (not closed["left_open"]
                or all(Fraction(x) == 0 for x in closed["sum_mA_|V(a,A,a)|^2"].values()))


def test_m_conservation_route1() -> None:
    """Sec. 2.3 Results: "Every nonzero amplitude conserves m in each SU(2) factor separately",
    on all 168 open route-1 cells.

    Rejects: a frame-side (right SU(2)) amplitude with M1 != M_A + M2, e.g. from a Pauli
    matrix element <s1|sigma_mu|s2> with the wrong Delta m; and a left-factor Clebsch-Gordan
    routine without its m selection.  Route: route 1's computed frame-side tensor Rt is
    scanned for nonzero entries off M1 = M_A + M2 (m_conserved); the left factor is route 1's
    Clebsch-Gordan coefficient, scanned for nonzero values off m1 = mA + m2 (with a count of
    nonzero on-shell values so the scan is not vacuous).
    """
    raw = route1_raw()
    opened = [v for v in raw.values() if v["open"]]
    assert len(opened) == 168
    assert all(v["m_conserved"] for v in opened)
    n_on = 0
    for tjA in range(0, 9, 2):
        for tj2 in range(0, 5):
            for tj1 in range(abs(tjA - tj2), tjA + tj2 + 1, 2):
                for tmA in R1.mrange(tjA):
                    for tm2 in R1.mrange(tj2):
                        for tm1 in R1.mrange(tj1):
                            F = R1.cg2(tjA, tmA, tj2, tm2, tj1, tm1)[1]
                            if tm1 != tmA + tm2:
                                assert F == 0
                            else:
                                n_on += F != 0
    assert n_on > 100


def test_prop_ward_takahashi_route2() -> None:
    """prop:ward eq:ward: <u, c(df) v> = (lambda_u - lambda_v) <u, f v>, exact, f harmonic of degree 1, 2, spinor levels n <= 1.

    Rejects: the n-block-diagonality proxy, which would accept any c(df) vanishing
    between different levels (including c(df) = 0) and is false here: c(df) couples
    n_u != n_v. Also rejects the identity with the opposite sign (lambda_v - lambda_u).
    Independent route: route 2 polynomial spinors, gradients and moments.
    """
    sp = route2(1, 3).SP
    n_pairs = n_nonzero = n_cross_level = 0
    for p in (1, 2):
        for f in R2.harmonic_basis(p)[:3]:
            grad = R2.tangential_gradient(f, p)
            for (su, nu) in sp:
                for (sv, nv) in sp:
                    lu = Fraction(su * (2 * nu + 3), 2)
                    lv = Fraction(sv * (2 * nv + 3), 2)
                    for u in sp[(su, nu)][:2]:
                        for v in sp[(sv, nv)][:2]:
                            lhs = R2.inner(u, R2.clifford(grad, v))
                            rhs = R2.inner(u, [R2.p_mul(f, v[0]), R2.p_mul(f, v[1])]) * (lu - lv)
                            assert (lhs - rhs).iszero()
                            n_pairs += 1
                            if not lhs.iszero():
                                n_nonzero += 1
                                assert lu != lv
                                n_cross_level += nu != nv
    assert n_pairs == 384 and n_nonzero == 56
    assert n_cross_level > 0


# ---------------------------------------------------------------------------
# Unnormalized one-loop sum
# ---------------------------------------------------------------------------
def _sigma0(nmax: int, S_fn) -> Fraction:
    """Sec. 2.6: n_ext = 0 (sign +), both internal signs, sum over all m of |V|^2 per
    external state, unit coupling, Haar volume one, photon weight q(q+2), |lambda|^4."""
    d_ext = 2
    tot = Fraction(0)
    for s_int in SIGNS:
        for n_int in range(nmax + 1):
            lam4 = Fraction(2 * n_int + 3, 2) ** 4
            for q in range(1, n_int + 3):      # theorem: nothing open beyond q = n_int + 1 at n_ext = 0
                for comp in COMPS:
                    S = S_fn(1, 0, comp, q, s_int, n_int)
                    if S:
                        tot += S / (d_ext * lam4 * q * (q + 2))
    return tot


def test_sigma_phys_ground_nmax3() -> None:
    """Sec. 2.6: the unnormalized one-loop sum at n_ext = 0, n_max = 3 is 2628144752/9845600625 = 0.2669... (nonzero).

    Rejects: zero (the model's Theorem 4 value), and the same sum with the Hodge
    weight (q+1)^2 in place of q(q+2) (= 0.2069, which the paper says is owed, not
    used). Two routes: route 1 brute-force amplitudes and the 9j closed form.
    """
    tab = route1_table()
    v1 = _sigma0(3, lambda s1, n1, c, q, s2, n2: tab[(s1, s2, c, n1, n2, q)][1])
    v9 = _sigma0(3, ninej_sum)
    assert v1 == v9 == SIGMA0_NMAX3
    assert round(float(v1), 4) == SIGMA0_PRINTED[3]


def test_sigma_phys_ground_nmax5() -> None:
    """Sec. 2.6: one cutoff past n_max = 3 (route 1's range), n_max = 5 gives 0.2702 by the 9j closed form.

    Rejects: a printed n_max = 5 value that the closed form does not reproduce
    (guard-asymptotics check one cutoff past the n_max = 3 data).
    """
    assert round(float(_sigma0(5, ninej_sum)), 4) == SIGMA0_PRINTED[5]


def test_sigma_phys_ground_nmax8_nmax12() -> None:
    """Sec. 2.6: the printed values 0.2715 and 0.2719 at n_max = 8 and 12, by the 9j closed
    form (exact rationals), and the sequence 3, 5, 8, 12 increasing.

    Rejects: printed n_max = 8, 12 values the closed form does not reproduce (e.g. 0.2702
    carried forward, or a sum that has saturated at n_max = 5), and a decreasing sequence.
    About 0.3 s.
    """
    v = {n: _sigma0(n, ninej_sum) for n in (8, 12)}
    assert round(float(v[8]), 4) == SIGMA0_PRINTED[8]
    assert round(float(v[12]), 4) == SIGMA0_PRINTED[12]
    assert SIGMA0_NMAX3 < _sigma0(5, ninej_sum) < v[8] < v[12]


# ---------------------------------------------------------------------------
# geovac.qed_vertex.so4_channel_count
# ---------------------------------------------------------------------------
def test_so4_channel_count_matches_theorem() -> None:
    """cor:model_comparison in code: geovac.qed_vertex.so4_channel_count equals the theorem's same-sign count on the model window, 0 outside; W in {0, 2}; sum over q = 2 min(n1, n2).

    Rejects: the old count (j2 content (n2/2, (n2+1)/2), no integrality), whose q-sum
    is 2 min - 1 - delta and which disagrees at 9 triples with n <= 3.
    Independent route: the theorem's rule transcribed from the paper.
    """
    from geovac.qed_vertex import so4_channel_count
    for n1 in range(7):
        for n2 in range(7):
            tot = 0
            for q in range(1, n1 + n2 + 4):
                w = so4_channel_count(n1, n2, q)
                if model_rule(n1, n2, q):
                    want = sum(theorem_rule(1, 1, c, n1, n2, q) for c in COMPS)
                    assert want == 2
                else:
                    want = 0
                assert w == want, (n1, n2, q, w, want)
                tot += w
            assert tot == 2 * min(n1, n2)
