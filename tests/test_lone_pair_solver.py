"""
Tests for the lone pair solver (solve_lone_pair).

Validates that the Level 3 hyperspherical wrapper used for lone pair blocks
in composed H₂O (Paper 17) produces correct energies and channel structure.

The inter-fiber coupling tests (S·F⁰) are archived in
tests/_archive/dead_ends/test_lone_pair_coupling.py — coupling at Z_eff=6
was found to be unphysical (CLAUDE.md §3).
"""

import numpy as np
import pytest

from geovac.lone_pair import solve_lone_pair, extract_channel_data_level3


@pytest.fixture(scope='module')
def he_result():
    """He ground state via lone pair solver (Z_eff=2)."""
    return solve_lone_pair(Z_eff=2.0, l_max=0, n_alpha=200, n_Re=3000,
                           N_R_angular=100, R_max=30.0, verbose=False)


@pytest.fixture(scope='module')
def o_result():
    """O lone pair via lone pair solver (Z_eff=6)."""
    return solve_lone_pair(Z_eff=6.0, l_max=2, n_alpha=80, n_Re=1500,
                           R_max=30.0, verbose=False)


class TestLonePairEnergy:
    """Energy accuracy for the lone pair solver."""

    def test_he_energy_accuracy(self, he_result):
        """He ground state within 0.1% of exact (-2.903724 Ha)."""
        E_exact = -2.903724
        error_pct = abs((he_result['energy'] - E_exact) / E_exact) * 100
        assert error_pct < 0.1, (
            f"He energy {he_result['energy']:.6f} Ha, error {error_pct:.4f}%"
        )

    def test_o_bound_state(self, o_result):
        """O lone pair (Z_eff=6) is bound."""
        assert o_result['energy'] < 0.0

    def test_o_energy_scales_with_z(self, o_result):
        """Z^2 scaling: E(Z=6)/E(Z=2) > 5."""
        he = solve_lone_pair(Z_eff=2.0, l_max=2, n_alpha=80,
                             n_Re=1500, verbose=False)
        ratio = o_result['energy'] / he['energy']
        assert ratio > 5.0, f"Z^2 scaling ratio {ratio:.2f}, expected > 5"


class TestChannelStructure:
    """Channel decomposition from extract_channel_data_level3."""

    @pytest.fixture(scope='class')
    def he_data(self):
        result = solve_lone_pair(Z_eff=2.0, l_max=2, n_alpha=80,
                                 n_Re=2000, verbose=False)
        return extract_channel_data_level3(
            result, Z=2.0, l_max=2, n_alpha=80, n_sample_R=10
        )

    def test_channel_count(self, he_data):
        """n_ch = l_max + 1."""
        assert he_data['n_ch'] == 3

    def test_dominant_channel_is_s(self, he_data):
        """(0,0) channel weight > 0.9 for He ground state."""
        cw_avg = np.mean(he_data['ch_weights'], axis=0)
        assert cw_avg[0] > 0.9

    def test_ch_weights_normalized(self, he_data):
        """Channel weights sum to ~1 at each sample point."""
        cw = he_data['ch_weights']
        for k in range(cw.shape[0]):
            assert abs(np.sum(cw[k]) - 1.0) < 0.05
