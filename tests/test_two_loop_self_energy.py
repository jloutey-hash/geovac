"""Tests for geovac/two_loop_self_energy.py (LS-8a).

These tests verify the structural properties of the iterated CC
spectral action two-loop SE on Dirac-S^3.  They do NOT assert
agreement with the literature value B_40 = 0.539 (Eides 2001, Eq. (59)) because the bare
spectral sum is UV-divergent — see debug/ls8a_two_loop_self_energy_memo.md
for the WEAK verdict.
"""
from __future__ import annotations

import mpmath
import pytest

from geovac.two_loop_self_energy import (
    ALPHA,
    c_2s_convergence,
    c_2s_native,
    crossed_se_spectral_sum,
    iterated_se_dimensionless,
    ls8a_prefactor_MHz,
    rainbow_se_spectral_sum,
)


def test_prefactor_matches_LS7():
    """LS-7 prefactor (alpha/pi)^2 (Z alpha)^4 / n^3 = 0.2363 MHz per unit
    coefficient for Z=1, n=2 (cited by Paper 34's two-loop row and Paper 35)."""
    pref = ls8a_prefactor_MHz(n=2, Z=1)
    assert abs(float(pref) - 0.2363) < 1e-3


def test_alpha_value():
    """CODATA 2018 fine-structure constant."""
    assert abs(float(ALPHA) - 1/137.035999084) < 1e-15


# ---------------------------------------------------------------------------
# Structural tests of the spectral sum
# ---------------------------------------------------------------------------

def test_rainbow_spectral_sum_positive():
    """Rainbow spectral sum is positive for n_ext=1."""
    result = rainbow_se_spectral_sum(n_ext=1, n_max=3)
    assert result > 0


def test_crossed_spectral_sum_positive():
    """Crossed spectral sum is positive for n_ext=1."""
    result = crossed_se_spectral_sum(n_ext=1, n_max=3)
    assert result > 0


def test_iterated_se_decomposition():
    """iterated_se_dimensionless returns rainbow + crossed = total."""
    sums = iterated_se_dimensionless(n_ext=1, n_max=3)
    R, C, T = sums["rainbow"], sums["crossed"], sums["total"]
    assert abs(R + C - T) < mpmath.mpf("1e-20")


def test_rainbow_only_excludes_crossed():
    """include_crossed=False zeros the crossed contribution."""
    sums = iterated_se_dimensionless(n_ext=1, n_max=3, include_crossed=False)
    assert sums["crossed"] == 0
    assert sums["rainbow"] > 0
    assert sums["total"] == sums["rainbow"]


def test_crossed_only_excludes_rainbow():
    sums = iterated_se_dimensionless(n_ext=1, n_max=3, include_rainbow=False)
    assert sums["rainbow"] == 0
    assert sums["crossed"] > 0


# ---------------------------------------------------------------------------
# Sign of C_2S^native
# ---------------------------------------------------------------------------

def test_c_2s_sign_correct_at_small_n_max():
    """C_2S^GeoVac has the correct positive sign matching literature."""
    r = c_2s_native(n_max=3)
    assert r["c_2s_native"] > 0
    assert r["c_2s_literature"] > 0
    assert r["ratio_native_over_lit"] > 0


# ---------------------------------------------------------------------------
# Convergence behavior (DIVERGENT — verified)
# ---------------------------------------------------------------------------

def test_spectral_sum_diverges_with_n_max():
    """Bare two-loop SE sum is UV-divergent: monotonically growing in n_max.

    This is the key structural finding of LS-8a — the framework faithfully
    reproduces the flat-space two-loop QED UV divergence.  See the memo.
    """
    sums = []
    for n_max in [3, 4, 5]:
        sums.append(float(c_2s_native(n_max=n_max)["raw_total"]))
    # Strictly increasing
    assert sums[1] > sums[0]
    assert sums[2] > sums[1]
    # Corrected W (2026-11-07): the raw totals are pinned exactly against an
    # independent Fraction route (W = 2 on every window triple;
    # debug/w_channel_count_correction.py, section E).  The old count gave
    # step ratios > 2; the corrected ones are 2.15 and 1.85, i.e. local
    # power-law slopes 2.66 and 2.76 -- growth faster than n_max^2.5.
    expected = [2994.777739248162, 6439.294240027892, 11902.54730011679]
    for got, exp in zip(sums, expected):
        assert abs(got - exp) < 1e-8 * exp
    import math
    assert math.log(sums[1] / sums[0]) / math.log(4 / 3) > 2.5
    assert math.log(sums[2] / sums[1]) / math.log(5 / 4) > 2.5


def _raw_total_independent(n_max: int, n_ext: int = 1):
    """Rainbow + crossed raw two-loop SE sum by an independent exact route
    (no geovac import): W = 2 on every window triple (|a-b| <= q <= a+b,
    a+b+q odd, q >= 1), photon factor d_T/mu = 1, electron factor
    g_n/lambda_n^2 = 2(n+1)(n+2) / ((2n+3)/2)^2."""
    from fractions import Fraction

    def win(a, b, q):
        return q >= 1 and abs(a - b) <= q <= a + b and (a + b + q) % 2 == 1

    def e(n):
        return Fraction(8 * (n + 1) * (n + 2), (2 * n + 3) ** 2)
    rng = range(n_max + 1)
    qs = range(1, 2 * n_max + 2)
    tot = Fraction(0)
    for n1 in rng:
        for n2 in rng:
            for n3 in rng:
                w = e(n1) * e(n2) * e(n3) * 16
                # rainbow: q_out at (ext,n1),(n3,ext); q_in at (n1,n2),(n2,n3)
                a = sum(1 for q in qs if win(n_ext, n1, q) and win(n3, n_ext, q))
                b = sum(1 for q in qs if win(n1, n2, q) and win(n2, n3, q))
                # crossed: q1 at (ext,n1),(n2,n3); q2 at (n1,n2),(n3,ext)
                c = sum(1 for q in qs if win(n_ext, n1, q) and win(n2, n3, q))
                d = sum(1 for q in qs if win(n1, n2, q) and win(n3, n_ext, q))
                tot += w * (a * b + c * d)
    return tot


def test_raw_totals_match_independent_route():
    """LS-8a raw totals (Paper 36 tab:ls8a_native; the pins in
    test_spectral_sum_diverges_with_n_max), backed inline.

    Rejects: the old-count totals (W in {0,1,2} with the mixed-sign
    pairing; e.g. W(1,1,1) = 0), and any frozen literal that drifts from
    the definition.  Route: exact Fraction sum over the window with W = 2
    (no geovac), against the frozen literals and against c_2s_native
    (which calls so4_channel_count)."""
    expected = {3: 2994.777739248162, 4: 6439.294240027892}
    for n_max, lit in expected.items():
        ind = float(_raw_total_independent(n_max))
        assert abs(ind - lit) < 1e-9 * lit
        got = float(c_2s_native(n_max=n_max)["raw_total"])
        assert abs(got - ind) < 1e-9 * ind


def test_c_2s_convergence_helper():
    """Sweep across n_max values returns ordered list."""
    results = c_2s_convergence(n_max_values=[2, 3])
    assert len(results) == 2
    assert results[0]["n_max"] == 2
    assert results[1]["n_max"] == 3


# ---------------------------------------------------------------------------
# Normalization conventions
# ---------------------------------------------------------------------------

def test_norm_conventions_rescale_correctly():
    """Different normalization conventions rescale C_2S by their ratio."""
    r_4pi = c_2s_native(n_max=3, norm_convention="hopf_4pi_squared")
    r_2pi = c_2s_native(n_max=3, norm_convention="hopf_2pi_squared")
    r_unit = c_2s_native(n_max=3, norm_convention="unit")

    # 1/(4pi)^2 to 1/(2pi)^2 ratio is 1/4
    ratio = r_4pi["c_2s_native"] / r_2pi["c_2s_native"]
    assert abs(float(ratio) - 0.25) < 1e-10

    # unit norm gives the raw spectral sum
    assert r_unit["c_2s_native"] == r_unit["raw_total"]


def test_invalid_norm_raises():
    with pytest.raises(ValueError):
        c_2s_native(n_max=2, norm_convention="bogus")


# ---------------------------------------------------------------------------
# Cross-validation against existing one-loop self-energy infrastructure
# ---------------------------------------------------------------------------

def test_n_ext_zero_rainbow_well_defined():
    """For n_ext=0 (1S), the rainbow sum is computable; at the test's own
    parameters it is exactly 0.0 (the assert below is >= 0 accordingly)."""
    R = rainbow_se_spectral_sum(n_ext=0, n_max=3)
    assert R >= 0  # At very small n_max R might be 0 due to vertex parity


@pytest.mark.slow
def test_c_2s_native_at_n_max_4():
    """Paper 36 tab:ls8a_native, n_max = 4: C = raw total / (4 pi)^2 = 40.78.

    Rejects: the pre-correction 4.82 (old channel count) and any literature
    ratio (Paper 36 makes no literature comparison for the model sum).
    Route: the exact Fraction total of _raw_total_independent (no geovac).
    """
    r = c_2s_native(n_max=4, norm_convention="hopf_4pi_squared")
    expected = float(_raw_total_independent(4)) / float((4 * mpmath.pi) ** 2)
    assert abs(float(r["c_2s_native"]) - expected) < 1e-9 * expected
    assert abs(expected - 40.78) < 0.01


# ---------------------------------------------------------------------------
# Vertex selection rule consistency
# ---------------------------------------------------------------------------

def test_vertex_parity_inheritance():
    """The two-loop SE inherits the vertex parity selection rule from
    one-loop. 1st-cert strengthening (2026-07-03): the old body only
    asserted R,C > 0, which a permissive vertex also satisfies. The rule
    is enforced REDUNDANTLY -- both by the _vertex_allowed guard and by
    the SO(4) channel weights, which vanish on exactly the forbidden
    triples (monkeypatching the guard alone provably changes nothing).
    Inheritance is therefore asserted at the weight level: every triple
    the one-loop vertex forbids carries ZERO weight in the two-loop sum."""
    import geovac.two_loop_self_energy as mod

    R = rainbow_se_spectral_sum(n_ext=1, n_max=2)
    C = crossed_se_spectral_sum(n_ext=1, n_max=2)
    assert R > 0
    assert C > 0

    checked = forbidden = 0
    for a in range(0, 6):
        for b in range(0, 6):
            for q in range(1, a + b + 1):
                checked += 1
                if not mod._vertex_allowed(a, b, q):
                    forbidden += 1
                    assert mod._so4_channel_count(a, b, q) == 0, (a, b, q)
    assert forbidden > 0, "rule never fired in range -- test is vacuous"
    assert checked > forbidden  # and it does not forbid everything
