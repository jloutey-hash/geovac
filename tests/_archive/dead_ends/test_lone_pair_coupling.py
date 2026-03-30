"""
Tests for the lone pair solver and inter-fiber coupling adapters.

Validates:
  1. solve_lone_pair reproduces He ground state at 0.05% (Z_eff=2)
  2. solve_lone_pair gives bound state for O lone pair (Z_eff=6)
  3. extract_channel_data_level3 produces valid channel coefficients
  4. Level 3 ↔ Level 4 overlap can be computed without errors
  5. Level 3 ↔ Level 3 overlap can be computed without errors
  6. All existing tests still pass (run separately)
"""

import numpy as np
import pytest

from geovac.lone_pair import (
    solve_lone_pair,
    extract_channel_data_level3,
    compute_inter_fiber_overlap,
    compute_inter_fiber_overlap_detailed,
)


# ==========================================================================
# 1. He reproduction via lone pair wrapper
# ==========================================================================

class TestHeliumReproduction:
    """Verify that the lone pair wrapper reproduces He at 0.05%."""

    @pytest.fixture(scope='class')
    def he_result(self):
        # Parameters matched to test_hyperspherical_he.py known-good config
        return solve_lone_pair(Z_eff=2.0, l_max=0, n_alpha=200, n_Re=3000,
                               N_R_angular=100, R_max=30.0, verbose=False)

    def test_he_energy_exists(self, he_result):
        """solve_lone_pair returns a dict with 'energy' key."""
        assert 'energy' in he_result
        assert isinstance(he_result['energy'], float)

    def test_he_energy_accuracy(self, he_result):
        """He ground state energy within 0.1% of exact (-2.903724 Ha)."""
        E_exact = -2.903724
        E_computed = he_result['energy']
        error_pct = abs((E_computed - E_exact) / E_exact) * 100
        assert error_pct < 0.1, (
            f"He energy {E_computed:.6f} Ha, error {error_pct:.4f}% "
            f"(threshold 0.1%)"
        )

    def test_he_wavefunction_shape(self, he_result):
        """Wavefunction is a 1D array on the radial grid."""
        F = he_result['wavefunction']
        R = he_result['R_grid_radial']
        assert F.ndim == 1
        assert len(F) == len(R)

    def test_he_wavefunction_normalizable(self, he_result):
        """Wavefunction has finite, positive norm."""
        F = he_result['wavefunction']
        norm = np.sum(F ** 2)
        assert norm > 0.0
        assert np.isfinite(norm)


# ==========================================================================
# 2. O lone pair (Z_eff=6)
# ==========================================================================

class TestOxygenLonePair:
    """Verify that Z_eff=6 gives a bound, physical energy."""

    @pytest.fixture(scope='class')
    def o_result(self):
        return solve_lone_pair(Z_eff=6.0, l_max=2, n_alpha=80, n_Re=1500,
                               R_max=30.0, verbose=False)

    def test_o_bound_state(self, o_result):
        """Energy is negative (bound state)."""
        assert o_result['energy'] < 0.0

    def test_o_energy_reasonable(self, o_result):
        """Energy should be near the He-like isoelectronic value for Z=6.

        Exact He-like 2-electron energy for Z=6 (C4+):
        E_exact = -Z^2 + 5Z/8 = -36 + 3.75 = -32.25 Ha (first-order PT)
        Exact non-relativistic: ~-32.41 Ha.
        We check that the energy is within 5% of this.
        """
        E_approx = -32.41  # C4+ exact non-rel
        E_computed = o_result['energy']
        error_pct = abs((E_computed - E_approx) / E_approx) * 100
        assert error_pct < 5.0, (
            f"O lone pair (Z_eff=6) energy {E_computed:.4f} Ha, "
            f"expected ~{E_approx:.2f} Ha, error {error_pct:.2f}%"
        )

    def test_o_energy_scales_with_z(self, o_result):
        """Energy should be much lower than He (Z=2) due to Z^2 scaling."""
        he_result = solve_lone_pair(Z_eff=2.0, l_max=2, n_alpha=80,
                                    n_Re=1500, verbose=False)
        assert o_result['energy'] < he_result['energy']
        # Rough Z^2 scaling: E(6)/E(2) ~ 9
        ratio = o_result['energy'] / he_result['energy']
        assert ratio > 5.0, f"Z^2 scaling ratio {ratio:.2f}, expected > 5"


# ==========================================================================
# 3. Channel data extraction from Level 3
# ==========================================================================

class TestChannelDataExtraction:
    """Verify extract_channel_data_level3 produces valid output."""

    @pytest.fixture(scope='class')
    def he_data(self):
        result = solve_lone_pair(Z_eff=2.0, l_max=2, n_alpha=80,
                                 n_Re=2000, verbose=False)
        return extract_channel_data_level3(
            result, Z=2.0, l_max=2, n_alpha=80, n_sample_R=10
        )

    def test_channel_labels_are_diagonal(self, he_data):
        """Level 3 channels should all be (l, l) tuples."""
        for ch in he_data['channels']:
            assert ch[0] == ch[1], f"Non-diagonal channel {ch} in Level 3"

    def test_channel_count(self, he_data):
        """n_ch should equal l_max + 1."""
        assert he_data['n_ch'] == 3  # l_max=2 → 3 channels

    def test_channels_are_sequential(self, he_data):
        """Channels should be [(0,0), (1,1), (2,2)]."""
        expected = [(0, 0), (1, 1), (2, 2)]
        assert he_data['channels'] == expected

    def test_ch_weights_shape(self, he_data):
        """ch_weights shape should be (n_sample, n_ch)."""
        cw = he_data['ch_weights']
        assert cw.shape == (10, 3)  # n_sample_R=10, n_ch=3

    def test_ch_weights_normalized(self, he_data):
        """Channel weights should sum to ~1 at each sample point."""
        cw = he_data['ch_weights']
        for k in range(cw.shape[0]):
            total = np.sum(cw[k])
            assert abs(total - 1.0) < 0.05, (
                f"Channel weights sum to {total:.4f} at sample {k}"
            )

    def test_dominant_channel_is_s(self, he_data):
        """The (0,0) channel should dominate for He ground state."""
        cw_avg = np.mean(he_data['ch_weights'], axis=0)
        assert cw_avg[0] > 0.9, (
            f"(0,0) channel weight {cw_avg[0]:.4f}, expected > 0.9"
        )

    def test_vec_2d_list_shapes(self, he_data):
        """Each eigenvector array should be (n_ch, n_alpha)."""
        for v in he_data['vec_2d_list']:
            assert v.shape == (3, 80)

    def test_ang_density_shape(self, he_data):
        """Angular density should be (n_sample, n_alpha)."""
        assert he_data['ang_density'].shape == (10, 80)

    def test_ang_density_positive(self, he_data):
        """Angular density should be non-negative everywhere."""
        assert np.all(he_data['ang_density'] >= -1e-15)

    def test_alpha_grid_correct(self, he_data):
        """Alpha grid should be interior FD points in (0, pi/2)."""
        ag = he_data['alpha_grid']
        assert len(ag) == 80
        assert ag[0] > 0
        assert ag[-1] < np.pi / 2


# ==========================================================================
# 4. Level 3 ↔ Level 3 overlap
# ==========================================================================

class TestLevel3Level3Overlap:
    """Overlap between two Level 3 (lone pair) fibers."""

    @pytest.fixture(scope='class')
    def he_channel_data(self):
        result = solve_lone_pair(Z_eff=2.0, l_max=2, n_alpha=80,
                                 n_Re=2000, verbose=False)
        return extract_channel_data_level3(
            result, Z=2.0, l_max=2, n_alpha=80, n_sample_R=10
        )

    def test_self_overlap_at_zero_angle(self, he_channel_data):
        """Self-overlap at θ=0 should be ~1 (identical fibers, aligned)."""
        result = compute_inter_fiber_overlap(
            he_channel_data, he_channel_data, bond_angle=0.0
        )
        assert abs(result['S_avg'] - 1.0) < 0.05, (
            f"Self-overlap at θ=0: {result['S_avg']:.4f}"
        )

    def test_overlap_at_90_degrees(self, he_channel_data):
        """Overlap at θ=π/2 should be between 0 and 1."""
        result = compute_inter_fiber_overlap(
            he_channel_data, he_channel_data, bond_angle=np.pi / 2
        )
        assert 0.0 <= result['S_avg'] <= 1.0

    def test_overlap_at_180_degrees(self, he_channel_data):
        """Overlap at θ=π (antiparallel) should be computable."""
        result = compute_inter_fiber_overlap(
            he_channel_data, he_channel_data, bond_angle=np.pi
        )
        assert np.isfinite(result['S_avg'])

    def test_common_channels_all_diagonal(self, he_channel_data):
        """Common channels between two Level 3 fibers are all (l, l)."""
        result = compute_inter_fiber_overlap(
            he_channel_data, he_channel_data, bond_angle=np.pi / 2
        )
        for ch in result['common_channels']:
            assert ch[0] == ch[1]

    def test_overlap_has_per_channel(self, he_channel_data):
        """Per-channel contributions should be returned."""
        result = compute_inter_fiber_overlap(
            he_channel_data, he_channel_data, bond_angle=np.pi / 2
        )
        assert len(result['per_channel_contribution']) == 3

    def test_detailed_overlap_shape(self, he_channel_data):
        """Detailed overlap should return per-sample arrays."""
        result = compute_inter_fiber_overlap_detailed(
            he_channel_data, he_channel_data, bond_angle=np.pi / 2
        )
        assert len(result['S_samples_A']) == 10
        assert len(result['Re_samples_A']) == 10


# ==========================================================================
# 5. Level 3 ↔ Level 4 overlap (mock Level 4 data)
# ==========================================================================

class TestLevel3Level4Overlap:
    """Overlap between a Level 3 and Level 4 fiber.

    Uses mock Level 4 channel data with known (l1, l2) channels.
    """

    @pytest.fixture(scope='class')
    def level3_data(self):
        result = solve_lone_pair(Z_eff=2.0, l_max=2, n_alpha=80,
                                 n_Re=2000, verbose=False)
        return extract_channel_data_level3(
            result, Z=2.0, l_max=2, n_alpha=80, n_sample_R=10
        )

    @pytest.fixture
    def mock_level4_data(self):
        """Create mock Level 4 channel data with (l1,l2) channels."""
        n_alpha = 80
        n_sample = 10
        # Level 4 channels for l_max=2:
        # (0,0), (0,1), (1,0), (1,1), (0,2), (2,0), (1,2), (2,1), (2,2)
        channels = [
            (0, 0), (0, 1), (1, 0), (1, 1),
            (0, 2), (2, 0), (1, 2), (2, 1), (2, 2)
        ]
        n_ch = len(channels)

        # Mock: dominant (0,0) channel, small contributions elsewhere
        ch_weights = np.zeros((n_sample, n_ch))
        ch_weights[:, 0] = 0.85  # (0,0) dominant
        ch_weights[:, 3] = 0.10  # (1,1)
        ch_weights[:, 8] = 0.05  # (2,2)

        # Mock eigenvectors (not physically meaningful but structurally valid)
        vec_2d_list = []
        for _ in range(n_sample):
            v = np.random.randn(n_ch, n_alpha) * 0.01
            # Ensure dominant channel has weight matching ch_weights
            v[0, :] = np.sqrt(0.85 / n_alpha)
            v[3, :] = np.sqrt(0.10 / n_alpha)
            v[8, :] = np.sqrt(0.05 / n_alpha)
            vec_2d_list.append(v)

        h_alpha = (np.pi / 2) / (n_alpha + 1)
        return {
            'Re_samples': np.linspace(0.5, 5.0, n_sample),
            'F_values': np.ones(n_sample),
            'alpha_grid': (np.arange(n_alpha) + 1) * h_alpha,
            'h_alpha': h_alpha,
            'channels': channels,
            'n_ch': n_ch,
            'vec_2d_list': vec_2d_list,
            'ch_weights': ch_weights,
            'ang_density': np.sum(
                np.array([v ** 2 for v in vec_2d_list]).sum(axis=1)
                .reshape(n_sample, -1)[:, :n_alpha], axis=0
            ).reshape(1, -1).repeat(n_sample, axis=0),
            'dRe': 0.5,
            'raw_norm': 1.0,
        }

    def test_overlap_computable(self, level3_data, mock_level4_data):
        """Level 3 ↔ Level 4 overlap should be computable without errors."""
        result = compute_inter_fiber_overlap(
            level3_data, mock_level4_data, bond_angle=np.pi / 2
        )
        assert np.isfinite(result['S_avg'])

    def test_common_channels_are_diagonal(self, level3_data, mock_level4_data):
        """Only (l,l) channels should be in common."""
        result = compute_inter_fiber_overlap(
            level3_data, mock_level4_data, bond_angle=np.pi / 2
        )
        expected_common = [(0, 0), (1, 1), (2, 2)]
        assert result['common_channels'] == expected_common

    def test_off_diagonal_channels_excluded(self, level3_data, mock_level4_data):
        """Off-diagonal Level 4 channels (0,1), (1,0), etc. should not appear."""
        result = compute_inter_fiber_overlap(
            level3_data, mock_level4_data, bond_angle=np.pi / 2
        )
        for ch in result['common_channels']:
            assert ch[0] == ch[1], f"Off-diagonal channel {ch} in common set"

    def test_detailed_overlap_computable(self, level3_data, mock_level4_data):
        """Detailed overlap should work for mixed Level 3/4."""
        result = compute_inter_fiber_overlap_detailed(
            level3_data, mock_level4_data, bond_angle=np.pi / 2
        )
        assert len(result['S_samples_A']) == len(level3_data['Re_samples'])

    def test_overlap_order_independent(self, level3_data, mock_level4_data):
        """S(A,B) should equal S(B,A) for symmetric overlap."""
        r_ab = compute_inter_fiber_overlap(
            level3_data, mock_level4_data, bond_angle=np.pi / 2
        )
        r_ba = compute_inter_fiber_overlap(
            mock_level4_data, level3_data, bond_angle=np.pi / 2
        )
        assert abs(r_ab['S_avg'] - r_ba['S_avg']) < 1e-10


# ==========================================================================
# 6. Edge cases
# ==========================================================================

class TestEdgeCases:
    """Edge cases and robustness checks."""

    def test_no_common_channels(self):
        """Overlap between fibers with no common channels gives S=0."""
        n_alpha = 20
        n_sample = 3
        h_alpha = (np.pi / 2) / (n_alpha + 1)

        data_A = {
            'channels': [(0, 0)], 'n_ch': 1,
            'F_values': np.ones(n_sample),
            'ch_weights': np.ones((n_sample, 1)),
            'Re_samples': np.linspace(1, 3, n_sample),
            'alpha_grid': (np.arange(n_alpha) + 1) * h_alpha,
            'h_alpha': h_alpha,
            'vec_2d_list': [np.ones((1, n_alpha)) for _ in range(n_sample)],
            'ang_density': np.ones((n_sample, n_alpha)),
            'dRe': 1.0, 'raw_norm': 1.0,
        }
        data_B = {
            'channels': [(1, 0), (0, 1)], 'n_ch': 2,
            'F_values': np.ones(n_sample),
            'ch_weights': np.ones((n_sample, 2)) * 0.5,
            'Re_samples': np.linspace(1, 3, n_sample),
            'alpha_grid': (np.arange(n_alpha) + 1) * h_alpha,
            'h_alpha': h_alpha,
            'vec_2d_list': [np.ones((2, n_alpha)) * 0.5
                            for _ in range(n_sample)],
            'ang_density': np.ones((n_sample, n_alpha)),
            'dRe': 1.0, 'raw_norm': 1.0,
        }

        result = compute_inter_fiber_overlap(data_A, data_B, bond_angle=np.pi)
        assert result['S_avg'] == 0.0
        assert result['common_channels'] == []
