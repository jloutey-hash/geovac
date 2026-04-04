"""
Tests for ab initio R_ref derivation (pk_rref.py).

Validates that R_ref candidates are computed correctly from core screening
data and that the R-dependent PK weight function behaves as expected.
"""

import numpy as np
import pytest

from geovac.core_screening import CoreScreening
from geovac.pk_rref import (
    compute_rref_candidates,
    select_rref,
    pk_weight_r_dependent,
)
from geovac.ab_initio_pk import AbInitioPK


@pytest.fixture(scope="module")
def li_core():
    """Solved Li 1s^2 core."""
    core = CoreScreening(Z=3, l_max=1, n_alpha=100)
    core.solve()
    return core


@pytest.fixture(scope="module")
def li_pk(li_core):
    """Ab initio PK for Li."""
    return AbInitioPK(li_core, n_core=2)


class TestRrefCandidates:
    """Test computation of R_ref candidates."""

    def test_all_candidates_positive(self, li_core, li_pk):
        """All R_ref candidates must be positive."""
        candidates = compute_rref_candidates(li_core, li_pk)
        for name, val in candidates.items():
            assert val > 0, f"{name} = {val} is not positive"

    def test_all_candidates_reasonable_range(self, li_core, li_pk):
        """All R_ref candidates should be in the range 0.1 - 2.0 bohr for Li."""
        candidates = compute_rref_candidates(li_core, li_pk)
        for name, val in candidates.items():
            assert 0.1 < val < 2.0, (
                f"{name} = {val:.4f} outside expected range [0.1, 2.0] bohr"
            )

    def test_ordering_inv2_lt_avg(self, li_core, li_pk):
        """inv2-weighted radius should be smaller than mean radius."""
        candidates = compute_rref_candidates(li_core, li_pk)
        assert candidates['r_inv2'] < candidates['r_avg']

    def test_ordering_avg_lt_rms(self, li_core, li_pk):
        """Mean radius should be smaller than RMS radius (Jensen inequality)."""
        candidates = compute_rref_candidates(li_core, li_pk)
        assert candidates['r_avg'] < candidates['r_rms']

    def test_r_half_screen_exists(self, li_core, li_pk):
        """r_half_screen should exist for Li (Z=3, 2 core electrons)."""
        candidates = compute_rref_candidates(li_core, li_pk)
        assert 'r_half_screen' in candidates

    def test_r_half_screen_is_screening_transition(self, li_core, li_pk):
        """Z_eff at r_half_screen should equal Z - 1."""
        candidates = compute_rref_candidates(li_core, li_pk)
        r_half = candidates['r_half_screen']
        z_eff_at_half = li_core.z_eff(r_half)
        assert abs(z_eff_at_half - 2.0) < 0.05  # Z=3, Z-1=2

    def test_pk_width_matches_pk_B(self, li_core, li_pk):
        """r_pk_width should be 1/sqrt(B) from AbInitioPK."""
        candidates = compute_rref_candidates(li_core, li_pk)
        expected = 1.0 / np.sqrt(li_pk.B)
        assert abs(candidates['r_pk_width'] - expected) < 1e-10


class TestSelectRref:
    """Test R_ref selection."""

    def test_select_default(self, li_core, li_pk):
        """Default method (r_half_screen) should work."""
        r_ref = select_rref(li_core, li_pk)
        assert r_ref > 0

    def test_select_all_methods(self, li_core, li_pk):
        """All named methods should be selectable."""
        for method in ['r_half_screen', 'r_pk_width', 'r_avg',
                       'r_rms', 'r_median', 'r_inv2']:
            r_ref = select_rref(li_core, li_pk, method=method)
            assert r_ref > 0, f"method={method} gave non-positive R_ref"

    def test_select_invalid_raises(self, li_core, li_pk):
        """Invalid method should raise ValueError."""
        with pytest.raises(ValueError, match="not available"):
            select_rref(li_core, li_pk, method='nonexistent')


class TestPKWeightFunction:
    """Test R-dependent PK weight."""

    def test_weight_at_zero(self):
        """w_PK(0) = 0."""
        assert pk_weight_r_dependent(0.0, 1.0) == 0.0

    def test_weight_at_rref(self):
        """w_PK(R_ref) = 1.0 when cap=1."""
        assert pk_weight_r_dependent(1.0, 1.0, cap=1.0) == 1.0

    def test_weight_capped(self):
        """w_PK(R) should not exceed cap."""
        assert pk_weight_r_dependent(10.0, 1.0, cap=1.0) == 1.0

    def test_weight_linear_below_cap(self):
        """w_PK should be linear in R for R < R_ref * cap."""
        R_ref = 2.0
        R = 1.0
        expected = R / R_ref
        assert abs(pk_weight_r_dependent(R, R_ref) - expected) < 1e-10

    def test_weight_monotonic(self):
        """w_PK should be monotonically non-decreasing in R."""
        R_ref = 0.5
        R_vals = np.linspace(0, 5, 100)
        weights = [pk_weight_r_dependent(R, R_ref) for R in R_vals]
        assert all(w2 >= w1 for w1, w2 in zip(weights, weights[1:]))
