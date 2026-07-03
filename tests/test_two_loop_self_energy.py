"""Tests for geovac/two_loop_self_energy.py (LS-8a).

These tests verify the structural properties of the iterated CC
spectral action two-loop SE on Dirac-S^3.  They do NOT assert
agreement with the literature value C_2S = +3.63 because the bare
spectral sum is UV-divergent — see debug/ls8a_two_loop_self_energy_memo.md
for the WEAK verdict.
"""
from __future__ import annotations

import mpmath
import pytest

from geovac.two_loop_self_energy import (
    ALPHA,
    C_2S_LITERATURE,
    LS7_REFERENCE_2S_MHZ,
    c_2s_convergence,
    c_2s_native,
    crossed_se_spectral_sum,
    iterated_se_dimensionless,
    ls8a_prefactor_MHz,
    rainbow_se_spectral_sum,
    verdict,
)


def test_prefactor_matches_LS7():
    """LS-7 prefactor 0.2363 MHz/dim for Z=1, n=2."""
    pref = ls8a_prefactor_MHz(n=2, Z=1)
    assert abs(float(pref) - 0.2363) < 1e-3


def test_C_2S_literature_target():
    """Literature C_2S = +3.6266 from Eides 2001 Tab. 7.3."""
    assert abs(float(C_2S_LITERATURE) - 3.6266) < 1e-3


def test_LS7_reference_2S_MHz():
    """Eides 2001 Tab. 7.3 two-loop SE 2S = +0.857 MHz."""
    assert abs(float(LS7_REFERENCE_2S_MHZ) - 0.857) < 1e-3


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
    # By more than a factor of 2 each step (power-law growth)
    assert sums[1] / sums[0] > 2.0
    assert sums[2] / sums[1] > 2.0


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
# Verdict logic
# ---------------------------------------------------------------------------

def test_verdict_positive_tier():
    v = verdict(C_2S_LITERATURE * mpmath.mpf("1.05"))
    assert v["sign_correct"]
    assert "POSITIVE" in v["verdict_tier"]


def test_verdict_mixed_tier():
    v = verdict(C_2S_LITERATURE * mpmath.mpf("1.30"))
    assert v["sign_correct"]
    assert "WEAK: structural form correct" in v["verdict_tier"]  # renamed from MIXED, 1st-cert 2026-07-03


def test_verdict_weak_tier():
    """C_2S off by factor of ~5 lands in WEAK tier (within order of magnitude)."""
    v = verdict(C_2S_LITERATURE * mpmath.mpf("5.0"))
    assert "WEAK" in v["verdict_tier"]


def test_verdict_negative_tier():
    """C_2S off by 100x lands in NEGATIVE tier."""
    v = verdict(C_2S_LITERATURE * mpmath.mpf("100"))
    assert "NEGATIVE" in v["verdict_tier"]


def test_verdict_sign_check():
    """Verdict catches sign errors."""
    v = verdict(-C_2S_LITERATURE)
    assert not v["sign_correct"]


# ---------------------------------------------------------------------------
# Cross-validation against existing one-loop self-energy infrastructure
# ---------------------------------------------------------------------------

def test_n_ext_zero_rainbow_well_defined():
    """For n_ext=0 (1S), the rainbow sum is computable; at the test's own
    parameters it is exactly 0.0 (the assert below is >= 0 accordingly)."""
    R = rainbow_se_spectral_sum(n_ext=0, n_max=3)
    assert R >= 0  # At very small n_max R might be 0 due to vertex parity


@pytest.mark.slow
def test_c_2s_at_n_max_4_within_factor_of_2():
    """At n_max=4 with hopf_4pi_squared, C_2S ~ 4.82.

    NOT sub-percent agreement.  Slow because n_max=4 takes ~15s with the
    5-fold spectral sum.  This test documents the n_max=4 numerical value
    for regression purposes; it is NOT a claim of physical agreement.
    """
    r = c_2s_native(n_max=4, norm_convention="hopf_4pi_squared")
    # At n_max=4: C_2S = 4.82, ratio over literature = 1.33
    assert abs(float(r["ratio_native_over_lit"]) - 1.33) < 0.10


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
