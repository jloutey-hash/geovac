"""
Tests for general bond angle inter-fiber coupling.

Validates the generalization from linear (θ=π) to arbitrary bond angle θ,
using Legendre polynomial rotation phases P_l(cos θ) in place of (-1)^l.

Tests:
  1. bond_angle=π reproduces existing results (regression)
  2. S_avg at bond_angle=104.5° < S_avg at bond_angle=π (weaker coupling)
  3. bond_angle=0 gives S_avg = 1.0 (aligned fibers, maximum overlap)
  4. _channel_rotation_phases consistency at known angles
  5. ComposedTriatomicSolver accepts bond_angle parameter
"""

import numpy as np
import pytest

from geovac.inter_fiber_coupling import (
    _channel_rotation_phases,
    compute_overlap_diagnostic,
    compute_overlap_from_channel_data,
    extract_channel_data,
    full_exchange_inter_fiber_energy,
)
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.composed_triatomic import ComposedTriatomicSolver
from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK


# =====================================================================
# Fixtures
# =====================================================================

@pytest.fixture(scope='module')
def beh2_bond_result():
    """Solve a single BeH bond pair at R=2.5 for testing overlap."""
    # Set up core screening for Be
    core = CoreScreening(Z=4, l_max=2, n_alpha=200)
    core.solve(verbose=False)
    pk = AbInitioPK(core, n_core=2)
    pk_potentials = [{'C_core': pk.A, 'beta_core': pk.B, 'atom': 'A'}]

    Z_eff = 2.0  # Be(4) - 2 core electrons
    Z_ligand = 1.0
    R = 2.5

    result = solve_level4_h2_multichannel(
        R=R, Z_A=Z_eff, Z_B=Z_ligand, l_max=2, n_alpha=100,
        n_Re=300, verbose=False, pk_potentials=pk_potentials,
    )
    return {
        'result': result,
        'R': R,
        'Z_A': Z_eff,
        'Z_B': Z_ligand,
        'l_max': 2,
        'n_alpha': 100,
        'pk_potentials': pk_potentials,
    }


# =====================================================================
# 1. _channel_rotation_phases unit tests
# =====================================================================

class TestChannelRotationPhases:
    """Verify rotation phase computation at known angles."""

    def test_linear_pi_matches_parity(self):
        """At θ=π, P_l(-1) = (-1)^l, so phase = (-1)^{l1+l2}."""
        channels = [(0, 0), (0, 1), (1, 0), (1, 1), (0, 2), (2, 0)]
        phases = _channel_rotation_phases(channels, bond_angle=np.pi)
        expected = np.array([
            (-1) ** (l1 + l2) for l1, l2 in channels
        ], dtype=float)
        np.testing.assert_allclose(phases, expected, atol=1e-14)

    def test_aligned_zero_gives_ones(self):
        """At θ=0, P_l(1) = 1, so phase = 1 for all channels."""
        channels = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 2)]
        phases = _channel_rotation_phases(channels, bond_angle=0.0)
        np.testing.assert_allclose(phases, 1.0, atol=1e-14)

    def test_right_angle(self):
        """At θ=π/2, P_0(0)=1, P_1(0)=0, P_2(0)=-1/2."""
        channels = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1)]
        phases = _channel_rotation_phases(channels, bond_angle=np.pi / 2)
        # P_0(0)=1, P_1(0)=0, P_2(0)=-0.5
        expected = np.array([
            1.0 * 1.0,    # (0,0): P_0*P_0 = 1
            0.0 * 1.0,    # (1,0): P_1*P_0 = 0
            -0.5 * 1.0,   # (2,0): P_2*P_0 = -0.5
            1.0 * 0.0,    # (0,1): P_0*P_1 = 0
            0.0 * 0.0,    # (1,1): P_1*P_1 = 0
        ])
        np.testing.assert_allclose(phases, expected, atol=1e-14)

    def test_h2o_angle(self):
        """At θ=104.5° (1.824 rad), phases should be between 0 and 1."""
        channels = [(0, 0), (1, 0), (0, 1)]
        theta = np.radians(104.5)
        phases = _channel_rotation_phases(channels, bond_angle=theta)
        # (0,0) channel should always be 1 (P_0 = 1)
        assert abs(phases[0] - 1.0) < 1e-14
        # (1,0) and (0,1) should have |phase| < 1
        assert abs(phases[1]) < 1.0
        assert abs(phases[2]) < 1.0


# =====================================================================
# 2. Overlap diagnostic with bond angle
# =====================================================================

class TestOverlapBondAngle:
    """Test compute_overlap_diagnostic at various bond angles."""

    def test_pi_reproduces_existing(self, beh2_bond_result):
        """bond_angle=π should give same S_avg as the original formula."""
        d = beh2_bond_result
        # Explicit π
        ovlp_pi = compute_overlap_diagnostic(
            R=d['R'], Z_A=d['Z_A'], Z_B=d['Z_B'],
            l_max=d['l_max'], n_alpha=d['n_alpha'],
            level4_result=d['result'],
            pk_potentials=d['pk_potentials'],
            n_sample_Re=15,
            bond_angle=np.pi,
        )
        # Default (should also be π)
        ovlp_default = compute_overlap_diagnostic(
            R=d['R'], Z_A=d['Z_A'], Z_B=d['Z_B'],
            l_max=d['l_max'], n_alpha=d['n_alpha'],
            level4_result=d['result'],
            pk_potentials=d['pk_potentials'],
            n_sample_Re=15,
        )
        assert abs(ovlp_pi['S_avg'] - ovlp_default['S_avg']) < 1e-12

    def test_h2o_angle_weaker_than_pi(self, beh2_bond_result):
        """S_avg at 104.5° should be less than S_avg at π."""
        d = beh2_bond_result
        ovlp_pi = compute_overlap_diagnostic(
            R=d['R'], Z_A=d['Z_A'], Z_B=d['Z_B'],
            l_max=d['l_max'], n_alpha=d['n_alpha'],
            level4_result=d['result'],
            pk_potentials=d['pk_potentials'],
            bond_angle=np.pi,
        )
        ovlp_h2o = compute_overlap_diagnostic(
            R=d['R'], Z_A=d['Z_A'], Z_B=d['Z_B'],
            l_max=d['l_max'], n_alpha=d['n_alpha'],
            level4_result=d['result'],
            pk_potentials=d['pk_potentials'],
            bond_angle=np.radians(104.5),
        )
        # Weaker coupling at smaller angle: |S_avg(104.5°)| < |S_avg(π)|
        assert abs(ovlp_h2o['S_avg']) < abs(ovlp_pi['S_avg'])

    def test_zero_angle_max_overlap(self, beh2_bond_result):
        """bond_angle=0 (aligned fibers) should give S_avg ≈ 1.0."""
        d = beh2_bond_result
        ovlp_zero = compute_overlap_diagnostic(
            R=d['R'], Z_A=d['Z_A'], Z_B=d['Z_B'],
            l_max=d['l_max'], n_alpha=d['n_alpha'],
            level4_result=d['result'],
            pk_potentials=d['pk_potentials'],
            bond_angle=0.0,
        )
        # All phases are +1, so S = sum |c|^2 = 1 (normalization)
        assert abs(ovlp_zero['S_avg'] - 1.0) < 0.01


# =====================================================================
# 3. Channel data overlap (cached path)
# =====================================================================

class TestOverlapFromChannelData:
    """Test compute_overlap_from_channel_data with bond angle."""

    def test_consistent_with_direct(self, beh2_bond_result):
        """Cached and direct overlap should agree at H₂O angle."""
        d = beh2_bond_result
        theta = np.radians(104.5)

        # Direct
        ovlp_direct = compute_overlap_diagnostic(
            R=d['R'], Z_A=d['Z_A'], Z_B=d['Z_B'],
            l_max=d['l_max'], n_alpha=d['n_alpha'],
            level4_result=d['result'],
            pk_potentials=d['pk_potentials'],
            bond_angle=theta,
        )

        # Via channel data
        ch_data = extract_channel_data(
            d['result'], d['R'], d['Z_A'], d['Z_B'],
            d['l_max'], d['n_alpha'],
            pk_potentials=d['pk_potentials'],
        )
        ovlp_cached = compute_overlap_from_channel_data(
            ch_data, bond_angle=theta)

        # Should be close (different sampling, but same physics)
        assert abs(ovlp_direct['S_avg'] - ovlp_cached['S_avg']) < 0.05


# =====================================================================
# 4. Full exchange with bond angle
# =====================================================================

class TestFullExchangeBondAngle:
    """Test full_exchange_inter_fiber_energy with bond angle."""

    def test_pi_default_backward_compatible(self, beh2_bond_result):
        """Default (π) should give same result as explicit π."""
        d = beh2_bond_result
        exch_default = full_exchange_inter_fiber_energy(
            d['result'], d['R'],
            Z_A=d['Z_A'], Z_B=d['Z_B'],
            l_max=d['l_max'], n_alpha=d['n_alpha'],
            pk_potentials=d['pk_potentials'],
            n_sample_Re=10,
        )
        exch_pi = full_exchange_inter_fiber_energy(
            d['result'], d['R'],
            Z_A=d['Z_A'], Z_B=d['Z_B'],
            l_max=d['l_max'], n_alpha=d['n_alpha'],
            pk_potentials=d['pk_potentials'],
            n_sample_Re=10,
            bond_angle=np.pi,
        )
        assert abs(exch_default['E_exchange'] - exch_pi['E_exchange']) < 1e-12

    def test_h2o_weaker_exchange(self, beh2_bond_result):
        """Exchange at 104.5° should be weaker (less negative) than at π."""
        d = beh2_bond_result
        exch_pi = full_exchange_inter_fiber_energy(
            d['result'], d['R'],
            Z_A=d['Z_A'], Z_B=d['Z_B'],
            l_max=d['l_max'], n_alpha=d['n_alpha'],
            pk_potentials=d['pk_potentials'],
            n_sample_Re=10,
            bond_angle=np.pi,
        )
        exch_h2o = full_exchange_inter_fiber_energy(
            d['result'], d['R'],
            Z_A=d['Z_A'], Z_B=d['Z_B'],
            l_max=d['l_max'], n_alpha=d['n_alpha'],
            pk_potentials=d['pk_potentials'],
            n_sample_Re=10,
            bond_angle=np.radians(104.5),
        )
        # Exchange is negative (attractive). Weaker = less negative = larger.
        assert exch_h2o['E_exchange'] > exch_pi['E_exchange']


# =====================================================================
# 5. ComposedTriatomicSolver bond_angle parameter
# =====================================================================

class TestComposedTriatomicBondAngle:
    """Test that ComposedTriatomicSolver accepts and stores bond_angle."""

    def test_default_is_pi(self):
        """Default bond_angle should be π."""
        solver = ComposedTriatomicSolver.BeH2(l_max=1, verbose=False)
        assert abs(solver.bond_angle - np.pi) < 1e-14

    def test_custom_angle_stored(self):
        """Custom bond_angle should be stored correctly."""
        theta = np.radians(104.5)
        solver = ComposedTriatomicSolver.BeH2(
            l_max=1, verbose=False, bond_angle=theta)
        assert abs(solver.bond_angle - theta) < 1e-14

    def test_explicit_pi_matches_default(self):
        """Explicit bond_angle=π should behave identically to default."""
        solver_default = ComposedTriatomicSolver.BeH2(
            l_max=1, verbose=False)
        solver_pi = ComposedTriatomicSolver.BeH2(
            l_max=1, verbose=False, bond_angle=np.pi)
        assert abs(solver_default.bond_angle - solver_pi.bond_angle) < 1e-14
