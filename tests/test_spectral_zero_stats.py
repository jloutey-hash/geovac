"""Tests for Track RH-M spectral-side zero statistics driver.

Verifies:
(1) Dirichlet series D(s), D_even(s), D_odd(s) evaluate correctly at known
    values from Paper 28.
(2) Argument principle correctly counts zeros of a test polynomial.
(3) Synthetic Poisson spacings have CV ~ 1 and synthetic GUE spacings have CV ~ 0.42.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

# Path hack so test runs standalone
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import mpmath

from debug.compute_spectral_zero_stats import (
    D_full,
    D_even,
    D_odd,
    count_zeros_in_box,
    find_zeros_in_strip,
    synthetic_poisson_spacings,
    synthetic_gue_spacings,
    compute_spacing_stats,
    _cv_of_list,
)


@pytest.fixture(autouse=True)
def set_precision():
    """Ensure the mpmath precision is high enough for all tests."""
    mpmath.mp.dps = 40
    yield
    mpmath.mp.dps = 15  # reset


def test_D_full_at_s_equal_4():
    """D(4) = pi^2 - pi^4 / 12 (Paper 28, Table 1)."""
    val = D_full(4)
    ref = mpmath.pi ** 2 - mpmath.pi ** 4 / 12
    assert abs(val - ref) < 1e-30, \
        f"D(4) = {mpmath.nstr(val,30)}, expected {mpmath.nstr(ref,30)}"


def test_D_even_plus_D_odd_equals_D_full():
    """D_even(s) + D_odd(s) = D(s) for any Re(s) > 3 (where the series converges)."""
    for s_test in [4, 5, 6, mpmath.mpc(4, 1), mpmath.mpc(3.5, 2)]:
        val_sum = D_even(s_test) + D_odd(s_test)
        val_full = D_full(s_test)
        diff = abs(val_sum - val_full)
        assert diff < 1e-25, \
            f"D_even + D_odd != D at s={s_test}: diff={float(diff)}"


def test_D_even_matches_paper28_formula():
    """D_even(4) = pi^2/2 - pi^4/24 - 4G + 4*beta(4) (Paper 28 Eq. 13)."""
    val = D_even(4)
    G = mpmath.catalan
    beta4 = (mpmath.hurwitz(4, mpmath.mpf(1) / 4)
             - mpmath.hurwitz(4, mpmath.mpf(3) / 4)) / mpmath.power(4, 4)
    ref = (mpmath.pi ** 2 / 2 - mpmath.pi ** 4 / 24
           - 4 * G + 4 * beta4)
    assert abs(val - ref) < 1e-25, \
        f"D_even(4) != expected: diff = {abs(val - ref)}"


def test_argument_principle_on_polynomial():
    """Argument principle should correctly count zeros of (s-2)(s-4+i)."""
    def poly(s):
        s = mpmath.mpc(s)
        return (s - 2) * (s - 4 - mpmath.mpc(0, 1))

    # Both zeros are inside [1,5]x[-1,2]
    n = count_zeros_in_box(poly, 1.0, 5.0, -1.0, 2.0, n_per_edge=40)
    assert n == 2, f"Expected 2 zeros, got {n}"

    # Only zero at s=2 is inside [1,3]x[-1,1]
    n = count_zeros_in_box(poly, 1.0, 3.0, -1.0, 1.0, n_per_edge=40)
    assert n == 1, f"Expected 1 zero, got {n}"

    # No zeros in [5,7]x[-1,1]
    n = count_zeros_in_box(poly, 5.0, 7.0, -1.0, 1.0, n_per_edge=40)
    assert n == 0, f"Expected 0 zeros, got {n}"


def test_synthetic_poisson_cv():
    """Synthetic Poisson (exponential) spacings have CV close to 1."""
    poi = synthetic_poisson_spacings(1000, seed=17)
    cv = _cv_of_list(poi)
    # Theory: exponential dist has CV = 1 exactly. Finite-sample variance
    # gives roughly sqrt(2/n) for the estimator. At n=1000, expect CV within 0.1 of 1.
    assert 0.85 < cv < 1.15, f"Poisson synthetic CV = {cv}, expected ~1.00"


def test_synthetic_gue_cv():
    """Synthetic GUE spacings should have CV close to 0.42."""
    gue = synthetic_gue_spacings(400, seed=17)
    cv = _cv_of_list(gue)
    # Theory: Wigner surmise for beta=2 gives CV approx 0.4233.
    # At n=400, expect CV within 0.05 of 0.42.
    assert 0.35 < cv < 0.50, f"GUE synthetic CV = {cv}, expected ~0.42"


def test_compute_spacing_stats_insufficient_zeros():
    """compute_spacing_stats should handle empty / singleton cases gracefully."""
    stats = compute_spacing_stats([])
    assert stats["n_zeros"] == 0
    assert stats.get("cv") is None

    stats = compute_spacing_stats([mpmath.mpc(1, 5)])
    assert stats.get("cv") is None


def test_compute_spacing_stats_synthetic_poisson():
    """Feeding synthetic Poisson spacings through the stats pipeline."""
    # Build synthetic zeros: start at Im=10 and add Poisson-distributed spacings
    poi_spacings = synthetic_poisson_spacings(80, seed=2024)
    zeros = []
    y = 10.0
    for ds in poi_spacings:
        y += ds
        zeros.append(mpmath.mpc(0.5, y))

    stats = compute_spacing_stats(zeros, unfold_window=15)
    cv = stats["cv"]
    # After unfolding (which does nothing to stationary Poisson), CV should still be ~1
    assert cv is not None
    assert 0.7 < cv < 1.3, f"Unfolded Poisson CV = {cv}, expected ~1.00"


def test_D_full_has_pole_at_s_equals_1():
    """Residue at simple pole s=1 of D_full (Hurwitz zeta has residue 1)."""
    # The -1/2 * zeta(s, 3/2) term has residue -1/2 (Hurwitz has residue 1 at s=1).
    s_near = mpmath.mpc(1, 1e-12)
    val = D_full(s_near)
    residue = val * (s_near - 1)
    # Expected residue: -1/2
    assert abs(residue - mpmath.mpf(-1) / 2) < 1e-5, \
        f"Residue at s=1 = {residue}, expected -1/2"


def test_D_full_has_pole_at_s_equals_3():
    """Residue at simple pole s=3 of D_full (from 2*zeta(s-2, 3/2))."""
    s_near = mpmath.mpc(3, 1e-12)
    val = D_full(s_near)
    residue = val * (s_near - 3)
    # Expected residue: 2 (since 2 * zeta(s-2, 3/2) at s=3 has residue 2)
    assert abs(residue - mpmath.mpf(2)) < 1e-5, \
        f"Residue at s=3 = {residue}, expected 2"


def test_D_full_trivial_zero_at_s_equals_negative_2():
    """D_full(-2) = 0 (trivial zero; verified by direct computation)."""
    # We verify numerically.
    val = D_full(-2)
    assert abs(val) < 1e-25, f"D_full(-2) = {val}, expected 0"


def test_D_full_trivial_zero_at_s_equals_0():
    """D_full(0) = 0 (another trivial zero)."""
    val = D_full(0)
    assert abs(val) < 1e-25, f"D_full(0) = {val}, expected 0"


def test_find_zeros_on_polynomial():
    """Small-scale zero-finding on a cubic polynomial."""
    def poly(s):
        s = mpmath.mpc(s)
        return (s - mpmath.mpc(2, 3)) * (s - mpmath.mpc(1, 5)) * (s + 1)

    zeros = find_zeros_in_strip(
        poly,
        re_range=(-3.0, 4.0),
        im_range=(-1.0, 6.0),
        n_seeds_re=5,
        n_seeds_im=15,
        pole_near=[],  # no poles
    )
    # Expect 3 zeros: (-1, 0), (2, 3), (1, 5)
    assert len(zeros) >= 3, f"Expected at least 3 zeros, got {len(zeros)}"

    # Check each known zero was found
    expected = [(-1.0, 0.0), (2.0, 3.0), (1.0, 5.0)]
    for (er, ei) in expected:
        dist = min(abs(complex(float(mpmath.re(z)), float(mpmath.im(z)))
                       - complex(er, ei))
                   for z in zeros)
        assert dist < 1e-6, f"Zero at ({er},{ei}) not found (min dist = {dist})"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
