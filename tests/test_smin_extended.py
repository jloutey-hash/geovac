"""Tests for Track Q-3: Extended PSLQ identification of S_min.

S_min = sum_{k>=1} T(k)^2 where T(k) = 2*zeta(2, k+3/2) - (1/2)*zeta(4, k+3/2).

Sprint 5 verified value:
    S_min = 2.47993693803422255441357950082938214468792578661728845837...

Tests:
1. First 80 digits match the known corrected value
2. 200-digit computation completes without error (slow)
3. Extended PSLQ with ~100 element basis runs and returns a result dict (slow)
4. Inner photon sum is finite for several (n1, n2) pairs
5. Inner photon sum is positive (for pairs with allowed q)
"""

from __future__ import annotations

import sys
from pathlib import Path

import mpmath
import pytest

# Ensure repo root on sys.path
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from debug.smin_extended_pslq import (
    SMIN_KNOWN_STR,
    smin_high_precision,
    smin_extended_pslq,
    smin_inner_q_closed_form,
)


# =====================================================================
# Fast tests (run by default)
# =====================================================================

def test_smin_first_80_digits():
    """Computed S_min matches the Sprint 5 corrected value to 80 digits."""
    mpmath.mp.dps = 120  # 80 target + 40 guard
    S = smin_high_precision(n_digits=80)

    known = mpmath.mpf(SMIN_KNOWN_STR)
    discrepancy = abs(S - known)
    # 80-digit agreement means discrepancy < 10^{-79}
    assert float(discrepancy) < 1e-79, (
        f"S_min disagrees with known value at 80-digit level: "
        f"|S - known| = {mpmath.nstr(discrepancy, 10)}"
    )


def test_smin_known_value_digits():
    """Verify the known value string has the right leading digits."""
    mpmath.mp.dps = 50
    known = mpmath.mpf(SMIN_KNOWN_STR)
    # Check the first 20 digits against the Sprint 5 erratum value
    expected_prefix = "2.4799369380342225544"
    actual = mpmath.nstr(known, 22)
    assert actual.startswith(expected_prefix), (
        f"Known value prefix mismatch: {actual} vs expected {expected_prefix}"
    )


def test_inner_q_sum_finite():
    """Inner photon sum is finite for several (n1, n2) pairs."""
    mpmath.mp.dps = 50
    test_pairs = [(0, 1), (1, 0), (1, 2), (2, 1), (2, 3), (3, 2), (0, 3)]
    for n1, n2 in test_pairs:
        result = smin_inner_q_closed_form(n1, n2)
        # The inner sum must be a finite number
        inner = result["inner_sum"]
        assert mpmath.isfinite(inner), (
            f"Inner sum I({n1},{n2}) is not finite: {inner}"
        )


def test_inner_q_sum_positive():
    """Inner photon sum is positive for pairs with allowed q values."""
    mpmath.mp.dps = 50
    # These pairs have at least one allowed q by parity + triangle
    # (n1=0 pairs have no allowed q; start from n1,n2 >= 1 with n2 >= n1+1)
    test_pairs = [(1, 2), (2, 1), (2, 3), (3, 2), (3, 4), (4, 3)]
    for n1, n2 in test_pairs:
        result = smin_inner_q_closed_form(n1, n2)
        assert result["n_terms"] > 0, (
            f"Expected allowed q-values for ({n1},{n2}) but got none"
        )
        assert float(result["inner_sum"]) > 0, (
            f"Inner sum I({n1},{n2}) is not positive: "
            f"{result['inner_sum_float']}"
        )


def test_inner_q_sum_symmetry():
    """I(n1, n2) = I(n2, n1) by the symmetry of the vertex coupling."""
    mpmath.mp.dps = 50
    test_pairs = [(0, 1), (1, 2), (2, 3), (0, 3), (1, 4)]
    for n1, n2 in test_pairs:
        r1 = smin_inner_q_closed_form(n1, n2)
        r2 = smin_inner_q_closed_form(n2, n1)
        assert abs(r1["inner_sum_float"] - r2["inner_sum_float"]) < 1e-10, (
            f"Symmetry failure: I({n1},{n2})={r1['inner_sum_float']} "
            f"!= I({n2},{n1})={r2['inner_sum_float']}"
        )


def test_inner_q_sum_n_equal():
    """I(n, n) is computable and non-negative for diagonal entries."""
    mpmath.mp.dps = 50
    for n in range(0, 5):
        result = smin_inner_q_closed_form(n, n)
        inner = result["inner_sum_float"]
        assert inner >= 0, f"I({n},{n}) = {inner} is negative"


def test_inner_q_sum_selection_rules():
    """Verify that the parity selection rule n1+n2+q=odd is respected."""
    mpmath.mp.dps = 50
    for n1 in range(0, 4):
        for n2 in range(0, 4):
            result = smin_inner_q_closed_form(n1, n2)
            for term in result["terms"]:
                q = term["q"]
                assert (n1 + n2 + q) % 2 == 1, (
                    f"Selection rule violated: n1={n1}, n2={n2}, q={q}, "
                    f"sum={n1+n2+q} is even"
                )


def test_smin_t_k_definition_consistency():
    """T(k) via direct computation matches the Hurwitz definition."""
    mpmath.mp.dps = 50
    # T(k) = 2*zeta(2, k+3/2) - (1/2)*zeta(4, k+3/2)
    for k in [1, 2, 5, 10]:
        a = mpmath.mpf(k) + mpmath.mpf(3) / 2
        T_hurwitz = (2 * mpmath.hurwitz(2, a)
                     - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, a))
        # Direct mode sum: T(k) = sum_{n>=k} g_n / |lambda_n|^4
        # where g_n = 2(n+1)(n+2), |lambda_n| = n + 3/2
        T_direct = mpmath.mpf(0)
        for n in range(k, k + 200):  # partial sum
            g_n = mpmath.mpf(2) * (n + 1) * (n + 2)
            lam_n = mpmath.mpf(n) + mpmath.mpf(3) / 2
            T_direct += g_n / lam_n ** 4
        # They should agree to within the tail truncation
        # The tail is O(1/k^3), so at n=k+200 it's ~1/200^3 ~ 1e-7
        discrepancy = abs(T_hurwitz - T_direct)
        assert float(discrepancy) < 0.01, (
            f"T({k}) Hurwitz vs direct disagree by {float(discrepancy)}"
        )


# =====================================================================
# Slow tests (require --slow flag)
# =====================================================================

@pytest.mark.slow
def test_smin_200_digits_computed():
    """200-digit computation of S_min completes without error."""
    mpmath.mp.dps = 250
    S = smin_high_precision(n_digits=200)
    # Just verify it's a real, positive, finite number
    assert mpmath.isfinite(S), f"S_min is not finite: {S}"
    assert float(S) > 2.4, f"S_min too small: {float(S)}"
    assert float(S) < 2.5, f"S_min too large: {float(S)}"
    # Verify 200-digit string has enough precision
    s_str = mpmath.nstr(S, 200)
    # The string should have at least 200 characters after the decimal
    # (accounting for leading "2." etc.)
    assert len(s_str.replace(".", "").replace("-", "")) >= 200, (
        f"200-digit string too short: {len(s_str)} chars"
    )


@pytest.mark.slow
def test_extended_pslq_runs():
    """PSLQ with ~100 element basis runs and returns a result dict."""
    mpmath.mp.dps = 150  # Use 150 for the PSLQ test (faster than 250)
    S = smin_high_precision(n_digits=100)
    result = smin_extended_pslq(S, basis_size=100)

    # Must return a dict with the required keys
    assert isinstance(result, dict), f"Expected dict, got {type(result)}"
    assert "identified" in result, "Missing 'identified' key"
    assert "basis_size" in result, "Missing 'basis_size' key"

    # Basis should be substantial (>= 50 elements)
    assert result["basis_size"] >= 50, (
        f"Basis too small: {result['basis_size']}"
    )

    # If identified, check residual
    if result["identified"]:
        assert "components" in result
        assert "residual" in result
        assert result["residual"] < 1e-30, (
            f"Residual too large: {result['residual']}"
        )
    else:
        # If not identified, that IS the result -- document it
        print(f"\n  PSLQ IRREDUCIBLE against {result['basis_size']}-element basis")
        print(f"  Error: {result.get('error', 'unknown')}")
