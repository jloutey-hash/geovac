"""
Tests for Track AK: 4-electron spectral compression analysis.

Validates channel counting, SO(12) spectrum, compression ratios,
and coupling matrix classification against the Track K pattern.
"""

import pytest
import numpy as np

from geovac.n_electron_spectral import (
    level4_spectral_dimensions,
    four_electron_spectral_dimensions,
    _count_s4_singlet_channels,
    _count_s4_22_sigma_channels,
    _count_s4_22_full_channels,
    _count_even_sum_tuples,
    so12_free_spectrum_by_channel,
    jacobi_tree_basis_structure,
    compression_analysis,
    coupling_matrix_classification,
    cross_pair_vee_analysis,
    tractability_verdict,
    run_full_analysis,
)
from geovac.n_electron_scope import (
    _so_d_degeneracy,
    casimir_eigenvalue,
    four_electron_channel_count_molecular,
    level4_channel_count,
)


# ==========================================================================
# (a) Level 4 reference (Track K pattern verification)
# ==========================================================================

class TestLevel4Reference:
    """Verify Level 4 spectral dimensions match Track K results."""

    def test_level4_lmax2_homonuclear(self):
        """Track K: 5 channels, 50 spectral dim, 20x compression."""
        dims = level4_spectral_dimensions(l_max=2, n_basis=10, homonuclear=True)
        assert dims['n_channels'] == 5  # (0,0), (0,2), (1,1), (2,0), (2,2)
        assert dims['spectral_dim'] == 50
        assert dims['fd_dim'] == 5 * 200  # 1000
        assert dims['compression'] == 20.0

    def test_level4_lmax0(self):
        """l_max=0: 1 channel."""
        dims = level4_spectral_dimensions(l_max=0, n_basis=10, homonuclear=True)
        assert dims['n_channels'] == 1
        assert dims['spectral_dim'] == 10

    def test_level4_lmax1(self):
        """l_max=1: (0,0), (1,1) = 2 channels for homonuclear."""
        dims = level4_spectral_dimensions(l_max=1, n_basis=10, homonuclear=True)
        assert dims['n_channels'] == 2
        assert dims['spectral_dim'] == 20

    def test_level4_angular_dim(self):
        """Level 4: S^5 has 5 angular dimensions."""
        dims = level4_spectral_dimensions(l_max=0)
        assert dims['angular_dim'] == 5
        assert dims['n_hyperangles'] == 1
        assert dims['n_direction_angles'] == 4


# ==========================================================================
# (b) SO(12) degeneracy verification
# ==========================================================================

class TestSO12Spectrum:
    """Verify SO(12) Casimir spectrum and degeneracies."""

    def test_so12_casimir_values(self):
        """mu_free(nu) = nu*(nu+10)/2."""
        assert casimir_eigenvalue(4, 0) == 0.0
        assert casimir_eigenvalue(4, 1) == 5.5
        assert casimir_eigenvalue(4, 2) == 12.0
        assert casimir_eigenvalue(4, 3) == 19.5
        assert casimir_eigenvalue(4, 4) == 28.0

    def test_so12_degeneracies(self):
        """SO(12) degeneracy formula: g(12, nu)."""
        assert _so_d_degeneracy(12, 0) == 1
        assert _so_d_degeneracy(12, 1) == 12
        assert _so_d_degeneracy(12, 2) == 77
        assert _so_d_degeneracy(12, 3) == 352

    def test_so6_degeneracies_reference(self):
        """SO(6) degeneracy for Level 4 reference."""
        assert _so_d_degeneracy(6, 0) == 1
        assert _so_d_degeneracy(6, 1) == 6
        assert _so_d_degeneracy(6, 2) == 20
        assert _so_d_degeneracy(6, 3) == 50

    def test_so12_degeneracy_growth(self):
        """SO(12) degeneracies grow much faster than SO(6)."""
        for nu in range(1, 6):
            assert _so_d_degeneracy(12, nu) > _so_d_degeneracy(6, nu)

    def test_so12_minimum_nu_singlet(self):
        """Singlet [2,2] has nu_min = 2 (from Paper 16: N - lambda_1 = 4 - 2)."""
        # The lowest allowed SO(12) eigenvalue for singlet
        mu_min = casimir_eigenvalue(4, 2)
        assert mu_min == 12.0


# ==========================================================================
# (c) 4-electron channel counting
# ==========================================================================

class TestFourElectronChannels:
    """Verify 4-electron channel counts."""

    def test_sigma_lmax0(self):
        """l_max=0: only (0,0,0,0) -> 1 channel."""
        result = four_electron_channel_count_molecular(
            0, sigma_only=True, homonuclear=False
        )
        assert result['sigma_channels'] == 1

    def test_sigma_lmax1_heteronuclear(self):
        """l_max=1, heteronuclear: 2^4 = 16 sigma channels."""
        result = four_electron_channel_count_molecular(
            1, sigma_only=True, homonuclear=False
        )
        assert result['sigma_channels'] == 16

    def test_sigma_lmax2_heteronuclear(self):
        """l_max=2, heteronuclear: 3^4 = 81 sigma channels."""
        result = four_electron_channel_count_molecular(
            2, sigma_only=True, homonuclear=False
        )
        assert result['sigma_channels'] == 81

    def test_full_lmax0(self):
        """l_max=0, full: all m_i=0 forced, so same as sigma."""
        result = four_electron_channel_count_molecular(
            0, sigma_only=False, homonuclear=False
        )
        assert result['full_channels'] == 1


# ==========================================================================
# (d) S4 singlet symmetry reduction
# ==========================================================================

class TestS4Reduction:
    """Verify S4 [2,2] symmetry reduction for singlet state."""

    def test_sigma_lmax0(self):
        """l_max=0: 1 channel, S4 irrelevant (trivial rep)."""
        n = _count_s4_22_sigma_channels(0, homonuclear=False)
        assert n >= 1

    def test_sigma_lmax1(self):
        """l_max=1: 16 raw sigma channels, [2,2] reduces."""
        n = _count_s4_22_sigma_channels(1, homonuclear=False)
        # 16/6 ~ 2.67, so expect 2-3 channels
        assert 1 <= n <= 16

    def test_sigma_lmax2(self):
        """l_max=2: 81 raw sigma, ~14 after S4."""
        n = _count_s4_22_sigma_channels(2, homonuclear=False)
        assert 5 <= n <= 81  # ~81/6 ~ 14

    def test_reduction_factor_approaches_sixth(self):
        """For large l_max, S4 reduction should approach 1/6."""
        for l_max in [3, 4]:
            raw = (l_max + 1) ** 4
            reduced = _count_s4_22_sigma_channels(l_max, homonuclear=False)
            ratio = reduced / raw
            # Should be approximately 1/6 ~ 0.167 for large l_max
            assert 0.05 < ratio < 0.5, f"l_max={l_max}: ratio={ratio}"

    def test_even_sum_tuples(self):
        """Verify even-sum tuple counting."""
        # l_max=1: values {0,1}, 4-tuples with even sum
        # Total 16, even sum: (0,0,0,0), (0,0,1,1), ..., (1,1,1,1)
        n = _count_even_sum_tuples(1, 4)
        assert n == 8  # half of 16


# ==========================================================================
# (e) 4-electron spectral dimensions
# ==========================================================================

class TestSpectralDimensions:
    """Verify spectral dimension computations."""

    def test_angular_dim_11(self):
        """4 electrons: angular space is S^11."""
        dims = four_electron_spectral_dimensions(0)
        assert dims['angular_dim'] == 11
        assert dims['n_hyperangles'] == 3
        assert dims['n_direction_angles'] == 8

    def test_alpha_basis_cubic(self):
        """Alpha basis size is n_basis^3 per channel."""
        for nb in [3, 5, 7]:
            dims = four_electron_spectral_dimensions(0, n_basis_per_angle=nb)
            assert dims['n_alpha_basis'] == nb ** 3

    def test_spectral_dim_consistent(self):
        """Spectral dim = n_channels_s4 * n_basis^3."""
        for l_max in [0, 1, 2]:
            for nb in [3, 5]:
                dims = four_electron_spectral_dimensions(
                    l_max, n_basis_per_angle=nb, sigma_only=True
                )
                expected = dims['n_channels_s4_exact'] * nb ** 3
                assert dims['spectral_dim_s4_exact'] == expected

    def test_lmax0_minimal(self):
        """l_max=0: 1 channel, minimal spectral dim."""
        dims = four_electron_spectral_dimensions(
            0, n_basis_per_angle=5, sigma_only=True
        )
        assert dims['n_channels_raw'] == 1
        assert dims['spectral_dim_s4_exact'] == dims['n_channels_s4_exact'] * 125


# ==========================================================================
# (f) SO(12) free spectrum with channel decomposition
# ==========================================================================

class TestFreeSpectrum:
    """Verify SO(12) free spectrum construction."""

    def test_lmax0_spectrum(self):
        """l_max=0: single channel (0,0,0,0), modes indexed by (k1,k2,k3)."""
        spec = so12_free_spectrum_by_channel(0, n_basis_per_angle=3)
        assert spec['n_channels'] == 1
        assert spec['n_modes_total'] == 27  # 3^3
        assert spec['n_modes_per_channel'] == 27

    def test_ground_mode_nu(self):
        """Ground mode (all k=0, all l=0): nu=0."""
        spec = so12_free_spectrum_by_channel(0, n_basis_per_angle=1)
        assert len(spec['modes_sorted']) == 1
        assert spec['modes_sorted'][0]['nu'] == 0
        assert spec['modes_sorted'][0]['mu_free'] == 0.0

    def test_modes_sorted_by_energy(self):
        """Modes are sorted by increasing mu_free."""
        spec = so12_free_spectrum_by_channel(0, n_basis_per_angle=5)
        energies = [m['mu_free'] for m in spec['modes_sorted']]
        for i in range(len(energies) - 1):
            assert energies[i] <= energies[i + 1]

    def test_nu_formula(self):
        """Verify nu = n12 + n34 + 2*k2 with n12 = l1+l2+2k1, etc."""
        spec = so12_free_spectrum_by_channel(1, n_basis_per_angle=2)
        for m in spec['modes_sorted']:
            l1, l2, l3, l4 = m['channel']
            k1, k2, k3 = m['k']
            n12 = l1 + l2 + 2 * k1
            n34 = l3 + l4 + 2 * k3
            nu = n12 + n34 + 2 * k2
            assert m['nu'] == nu
            assert m['mu_free'] == nu * (nu + 10) / 2.0


# ==========================================================================
# (g) Compression analysis
# ==========================================================================

class TestCompressionAnalysis:
    """Verify compression ratio computations."""

    def test_compression_table_structure(self):
        """Compression analysis returns correct structure."""
        results = compression_analysis([0, 1, 2])
        assert len(results) == 3
        for r in results:
            assert 'l_max' in r
            assert '4e_spectral_dim' in r
            assert 'l4_spectral_dim' in r
            assert 'ratio_to_l4' in r

    def test_4e_larger_than_l4(self):
        """4-electron spectral dim always exceeds Level 4."""
        results = compression_analysis([0, 1, 2])
        for r in results:
            if r['l_max'] > 0:
                assert r['4e_spectral_dim'] > r['l4_spectral_dim']

    def test_compression_positive(self):
        """FD-to-spectral compression ratio is positive."""
        results = compression_analysis([0, 1, 2])
        for r in results:
            assert r['4e_compression'] > 0


# ==========================================================================
# (h) Coupling matrix classification
# ==========================================================================

class TestCouplingClassification:
    """Verify coupling matrix classification."""

    def test_five_term_types(self):
        """Should have 5 types of coupling terms."""
        terms = coupling_matrix_classification()
        assert len(terms) == 5

    def test_total_vee_count(self):
        """Total V_ee terms: 2 intra + 4 cross = C(4,2) = 6."""
        terms = coupling_matrix_classification()
        vee_count = sum(t['count'] for t in terms if 'V_ee' in t['term'])
        assert vee_count == 6

    def test_vnuc_count(self):
        """V_nuc: 4 electrons x 2 nuclei = 8."""
        terms = coupling_matrix_classification()
        vnuc = [t for t in terms if 'V_nuc' in t['term']]
        assert len(vnuc) == 1
        assert vnuc[0]['count'] == 8

    def test_rho_dependence(self):
        """Only V_nuc is rho-dependent."""
        terms = coupling_matrix_classification()
        rho_dep = [t for t in terms if t['rho_dependent']]
        assert len(rho_dep) == 1
        assert 'V_nuc' in rho_dep[0]['term']

    def test_cross_pair_partial_algebraic(self):
        """Cross-pair V_ee is only partially algebraic."""
        terms = coupling_matrix_classification()
        cross = [t for t in terms if 'cross' in t['term']]
        assert len(cross) == 1
        assert cross[0]['algebraic'] == 'PARTIAL'

    def test_intra_pair_fully_algebraic(self):
        """Intra-pair V_ee is fully algebraic."""
        terms = coupling_matrix_classification()
        intra = [t for t in terms if 'intra' in t['term']]
        assert len(intra) == 1
        assert intra[0]['algebraic'] is True


# ==========================================================================
# (i) Cross-pair V_ee analysis
# ==========================================================================

class TestCrossPairVee:
    """Verify cross-pair V_ee structural analysis."""

    def test_has_gaunt_structure(self):
        """Cross-pair V_ee angular part uses Gaunt integrals."""
        analysis = cross_pair_vee_analysis()
        assert 'Gaunt' in analysis['gaunt_structure']

    def test_not_fully_separable(self):
        """Cross-pair hyperangular part is not fully separable."""
        analysis = cross_pair_vee_analysis()
        assert 'not fully algebraic' in analysis['separability']

    def test_reducible_to_nested_1d(self):
        """Cross-pair is reducible from 3D to nested 1D integrals."""
        analysis = cross_pair_vee_analysis()
        assert 'nested 1D' in analysis['separability']


# ==========================================================================
# (j) Tractability verdict
# ==========================================================================

class TestTractability:
    """Verify tractability assessment."""

    def test_lmax0_practical(self):
        """l_max=0 should be practical (< 1 minute)."""
        v = tractability_verdict(0, n_basis_per_angle=5)
        assert v['verdict'] in ['PRACTICAL', 'TRACTABLE']

    def test_lmax2_verdict(self):
        """l_max=2 should return a verdict."""
        v = tractability_verdict(2, n_basis_per_angle=5)
        assert v['verdict'] in ['PRACTICAL', 'TRACTABLE', 'INTRACTABLE']

    def test_cost_increases_with_lmax(self):
        """Cost should increase with l_max (at least for l_max >= 1)."""
        v1 = tractability_verdict(1, n_basis_per_angle=5)
        v2 = tractability_verdict(2, n_basis_per_angle=5)
        v3 = tractability_verdict(3, n_basis_per_angle=5)
        assert v1['total_time_s'] <= v2['total_time_s'] < v3['total_time_s']

    def test_cost_increases_with_nbasis(self):
        """Cost should increase with n_basis_per_angle."""
        v3 = tractability_verdict(1, n_basis_per_angle=3)
        v5 = tractability_verdict(1, n_basis_per_angle=5)
        v7 = tractability_verdict(1, n_basis_per_angle=7)
        assert v3['total_time_s'] < v5['total_time_s'] < v7['total_time_s']

    def test_ratio_to_l4(self):
        """Ratio to Level 4 should be computed and positive."""
        v = tractability_verdict(2, n_basis_per_angle=5)
        assert v['ratio_to_l4_dim'] > 1.0


# ==========================================================================
# (k) Full analysis runner
# ==========================================================================

class TestFullAnalysis:
    """Verify the full analysis report."""

    def test_report_runs(self):
        """Full analysis should complete without error."""
        report = run_full_analysis(n_basis_per_angle=3)
        assert len(report) > 100

    def test_report_contains_sections(self):
        """Report should contain all required sections."""
        report = run_full_analysis(n_basis_per_angle=3)
        assert 'COMPRESSION TABLE' in report
        assert 'COUPLING MATRIX' in report
        assert 'CROSS-PAIR' in report
        assert 'TRACTABILITY' in report
        assert 'SPECTRAL BASIS' in report
        assert 'SO(12)' in report


# ==========================================================================
# (l) Jacobi tree basis structure
# ==========================================================================

class TestJacobiTree:
    """Verify Jacobi tree basis structure analysis."""

    def test_alpha2_depends_on_k1_k3(self):
        """The alpha2 Jacobi parameters depend on pair quantum numbers."""
        info = jacobi_tree_basis_structure(l_max=0)
        assert info['alpha2_depends_on_k1_k3'] is True

    def test_tensor_product_dim(self):
        """Tensor product dimension is n_basis^3."""
        for nb in [3, 5, 7]:
            info = jacobi_tree_basis_structure(l_max=0, n_basis_per_angle=nb)
            assert info['tensor_product_dim'] == nb ** 3


# ==========================================================================
# (m) Consistency checks
# ==========================================================================

class TestConsistency:
    """Cross-check 4-electron results against known limits."""

    def test_so12_nu0_degeneracy_is_1(self):
        """Only one state at nu=0 (trivial representation)."""
        assert _so_d_degeneracy(12, 0) == 1

    def test_so12_degeneracy_formula(self):
        """Verify g(d, nu) = C(nu+d-1,d-1) - C(nu+d-3,d-1) for nu >= 1."""
        from math import comb
        d = 12
        for nu in range(1, 6):
            expected = comb(nu + d - 1, d - 1) - comb(nu + d - 3, d - 1)
            actual = _so_d_degeneracy(d, nu)
            assert actual == expected, f"nu={nu}: {actual} != {expected}"

    def test_level4_channel_count_consistency(self):
        """Level 4 channel count from scope.py matches our computation."""
        for l_max in [0, 1, 2, 3]:
            from_scope = level4_channel_count(l_max, sigma_only=True, homonuclear=True)
            dims = level4_spectral_dimensions(l_max, homonuclear=True)
            assert dims['n_channels'] == from_scope
