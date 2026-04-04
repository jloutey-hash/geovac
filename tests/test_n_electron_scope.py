"""
Tests for N-electron scoping module (Track AG).

Validates dimension counts, symmetry analysis, and feasibility estimates
for the 4-electron mol-frame hyperspherical solver scoping.
"""

import pytest
import numpy as np
from geovac.n_electron_scope import (
    n_electron_dimensions,
    casimir_eigenvalue,
    effective_angular_momentum,
    level4_channel_count,
    four_electron_channel_count_atomic,
    four_electron_channel_count_molecular,
    _count_m_zero_states,
    _count_L0_couplings,
    s4_irrep_dimensions,
    spin_spatial_pairing,
    so12_casimir_spectrum,
    _so_d_degeneracy,
    feasibility_estimate,
    separation_analysis_4e,
    hyperspherical_tree_4e,
    antisymmetry_analysis,
    pauli_exclusion_analysis,
)


class TestCoordinates:
    """Test coordinate dimension analysis."""

    def test_1_electron(self):
        d = n_electron_dimensions(1)
        assert d['config_dim'] == 3
        assert d['angular_dim'] == 2
        assert d['sphere'] == 'S^2'
        assert d['isometry_group'] == 'SO(3)'

    def test_2_electron(self):
        d = n_electron_dimensions(2)
        assert d['config_dim'] == 6
        assert d['angular_dim'] == 5
        assert d['sphere'] == 'S^5'
        assert d['isometry_group'] == 'SO(6)'

    def test_3_electron(self):
        d = n_electron_dimensions(3)
        assert d['config_dim'] == 9
        assert d['angular_dim'] == 8
        assert d['sphere'] == 'S^8'
        assert d['isometry_group'] == 'SO(9)'

    def test_4_electron(self):
        d = n_electron_dimensions(4)
        assert d['config_dim'] == 12
        assert d['angular_dim'] == 11
        assert d['sphere'] == 'S^11'
        assert d['isometry_group'] == 'SO(12)'


class TestCasimir:
    """Test SO(3N) Casimir eigenvalues."""

    def test_so3_casimir(self):
        """SO(3) l(l+1)/2 for N=1."""
        assert casimir_eigenvalue(1, 0) == 0.0
        assert casimir_eigenvalue(1, 1) == 1.0
        assert casimir_eigenvalue(1, 2) == 3.0

    def test_so6_casimir(self):
        """SO(6) nu(nu+4)/2 for N=2 — matches Paper 13."""
        assert casimir_eigenvalue(2, 0) == 0.0
        assert casimir_eigenvalue(2, 2) == 6.0
        assert casimir_eigenvalue(2, 4) == 16.0

    def test_so12_casimir(self):
        """SO(12) nu(nu+10)/2 for N=4."""
        assert casimir_eigenvalue(4, 0) == 0.0
        assert casimir_eigenvalue(4, 1) == 5.5
        assert casimir_eigenvalue(4, 2) == 12.0
        assert casimir_eigenvalue(4, 3) == 19.5

    def test_so12_ground_state_singlet(self):
        """4-electron singlet has nu_min=2, mu_free=12 (from Paper 16)."""
        mu = casimir_eigenvalue(4, 2)
        assert mu == 12.0

    def test_effective_angular_momentum(self):
        """Check l_eff for 4 electrons."""
        # l_eff = nu + (3N-4)/2 = nu + 4
        assert effective_angular_momentum(4, 0) == 4.0
        assert effective_angular_momentum(4, 2) == 6.0


class TestChannelCounting:
    """Test channel dimension computations."""

    def test_level4_sigma_homo_known(self):
        """Verify Level 4 channel counts against Paper 15 table."""
        assert level4_channel_count(0, sigma_only=True, homonuclear=True) == 1
        assert level4_channel_count(1, sigma_only=True, homonuclear=True) == 2
        assert level4_channel_count(2, sigma_only=True, homonuclear=True) == 5
        assert level4_channel_count(3, sigma_only=True, homonuclear=True) == 8
        assert level4_channel_count(4, sigma_only=True, homonuclear=True) == 13

    def test_m_zero_count_basic(self):
        """Basic M=0 counting."""
        # All l=0: only (0,0,0,0) -> 1
        assert _count_m_zero_states(0, 0, 0, 0) == 1
        # (1,0,0,0): m1 in {-1,0,1}, rest 0, need m1=0 -> 1
        assert _count_m_zero_states(1, 0, 0, 0) == 1
        # (1,1,0,0): m1+m2=0 -> 3 states
        assert _count_m_zero_states(1, 1, 0, 0) == 3

    def test_L0_couplings_basic(self):
        """Basic L=0 coupling counting."""
        # (0,0,0,0): L12=0, L34=0 -> 1
        assert _count_L0_couplings(0, 0, 0, 0) == 1
        # (1,1,0,0): L12=0,1,2; L34=0 -> L12=L34=0 -> 1
        assert _count_L0_couplings(1, 1, 0, 0) == 1
        # (1,1,1,1): L12=0,1,2; L34=0,1,2 -> L12=L34 -> 3
        assert _count_L0_couplings(1, 1, 1, 1) == 3

    def test_4e_sigma_lmax0(self):
        """At l_max=0, only 1 sigma channel."""
        r = four_electron_channel_count_molecular(0, sigma_only=True, homonuclear=False)
        assert r['sigma_channels'] == 1

    def test_4e_sigma_lmax1(self):
        """At l_max=1, heteronuclear sigma: (l_max+1)^4 = 16."""
        r = four_electron_channel_count_molecular(1, sigma_only=True, homonuclear=False)
        assert r['sigma_channels'] == 16

    def test_4e_sigma_lmax2(self):
        """At l_max=2, heteronuclear sigma: 3^4 = 81."""
        r = four_electron_channel_count_molecular(2, sigma_only=True, homonuclear=False)
        assert r['sigma_channels'] == 81

    def test_4e_atomic_lmax0(self):
        """Atomic 4e at l_max=0."""
        r = four_electron_channel_count_atomic(0)
        assert r['coupled_L0'] == 1

    def test_4e_atomic_lmax1(self):
        """Atomic 4e at l_max=1."""
        r = four_electron_channel_count_atomic(1)
        assert r['coupled_L0'] > 1


class TestSymmetry:
    """Test symmetry reduction analysis."""

    def test_s4_irrep_dimensions_sum(self):
        """Sum of dim^2 = |S_4| = 24."""
        irreps = s4_irrep_dimensions()
        total = sum(d**2 for d in irreps.values())
        assert total == 24

    def test_singlet_spatial_irrep(self):
        """LiH singlet ground state has [2,2] spatial irrep."""
        sp = spin_spatial_pairing()
        assert sp['S=0 (singlet)']['spatial_irrep'] == '[2,2]'
        assert sp['S=0 (singlet)']['spatial_dim'] == 2

    def test_so_degeneracy_nu0(self):
        """Degeneracy at nu=0 is always 1."""
        for d in [3, 6, 9, 12]:
            assert _so_d_degeneracy(d, 0) == 1

    def test_so_degeneracy_nu1(self):
        """Degeneracy at nu=1 is d."""
        for d in [3, 6, 9, 12]:
            assert _so_d_degeneracy(d, 1) == d

    def test_so6_degeneracy_known(self):
        """SO(6) degeneracies from Paper 13."""
        # nu=0: 1, nu=1: 6, nu=2: 20
        assert _so_d_degeneracy(6, 0) == 1
        assert _so_d_degeneracy(6, 1) == 6
        # nu=2: C(7,2) - C(3,0) = 21 - 1 = 20
        assert _so_d_degeneracy(6, 2) == 20


class TestSO12Spectrum:
    """Test SO(12) Casimir spectrum."""

    def test_ground_state(self):
        """nu=0 has mu_free=0, degeneracy=1."""
        spec = so12_casimir_spectrum(5)
        assert spec[0]['mu_free'] == 0.0
        assert spec[0]['degeneracy'] == 1

    def test_nu2_singlet_ground(self):
        """nu=2 is the lowest allowed state for 4e singlet."""
        spec = so12_casimir_spectrum(5)
        assert spec[2]['mu_free'] == 12.0
        assert spec[2]['nu'] == 2

    def test_so12_degeneracy_nu1(self):
        """SO(12) nu=1 degeneracy = 12."""
        spec = so12_casimir_spectrum(5)
        assert spec[1]['degeneracy'] == 12

    def test_so12_degeneracy_nu2(self):
        """SO(12) nu=2 degeneracy.

        g(12, 2) = (2*2+10)/10 * C(11,2) = 14/10 * 55 = 77.
        """
        spec = so12_casimir_spectrum(5)
        assert spec[2]['degeneracy'] == 77


class TestSeparation:
    """Test separation structure analysis."""

    def test_rho_collapse(self):
        """rho = R/(2*R_e) collapse should still work for 4 electrons."""
        s = separation_analysis_4e()
        assert s['rho_collapse'] is True
        assert s['V_ee_rho_independent'] is True

    def test_tree_structure(self):
        """Verify 4e tree gives correct total coordinates."""
        tree = hyperspherical_tree_4e()
        # 1 + 3 + 8 = 12 = 3*4
        assert '12 = 3*4' in tree['tree_structure']


class TestFeasibility:
    """Test feasibility estimates."""

    def test_lmax0_small(self):
        """l_max=0 should be trivially small."""
        f = feasibility_estimate(0, n_basis_alpha=5)
        assert f['sigma_channels_raw'] == 1
        assert f['angular_dim_sigma'] <= 200

    def test_lmax2_dimension(self):
        """l_max=2 angular dimension check."""
        f = feasibility_estimate(2, n_basis_alpha=10)
        # 81 sigma channels / 6 ~ 13 antisym channels
        # 13 * 10^3 = 13000 angular dim
        assert f['angular_dim_sigma'] > 1000

    def test_level4_comparison(self):
        """Level 4 should be much smaller than 4e."""
        f = feasibility_estimate(2, n_basis_alpha=10)
        assert f['angular_dim_sigma'] > f['level4_dim']

    def test_feasibility_returns_times(self):
        """Check that time estimates are computed."""
        f = feasibility_estimate(2, n_basis_alpha=10)
        assert 'total_time_sigma_s' in f
        assert f['total_time_sigma_s'] > 0


class TestAntisymmetry:
    """Test antisymmetry analysis."""

    def test_antisymmetry_has_mixed_constraint(self):
        a = antisymmetry_analysis()
        assert 'MIXED constraint' in a['gerade_analog']

    def test_pauli_exclusion(self):
        p = pauli_exclusion_analysis()
        assert 'PARTIAL' in p['answer']


class TestCrossPairCusps:
    """Test that cross-pair cusp analysis is correct."""

    def test_cusp_locations(self):
        tree = hyperspherical_tree_4e()
        assert 'cross-pair cusps' in tree['cusp_locations'].lower() or \
               'Cross pairs' in tree['cusp_locations']


class TestDimensionTables:
    """Generate dimension tables for the analysis document."""

    def test_generate_tables(self):
        """Compute and verify all dimension tables."""
        # Coordinate dimensions
        for N in [1, 2, 3, 4]:
            d = n_electron_dimensions(N)
            assert d['angular_dim'] == 3 * N - 1

        # SO(12) spectrum
        spec = so12_casimir_spectrum(6)
        assert spec[0]['mu_free'] == 0.0
        assert spec[2]['mu_free'] == 12.0  # singlet ground

        # Channel counts: Level 4 (2e) vs Full (4e)
        results = {}
        for lm in [0, 1, 2, 3, 4]:
            l4h = level4_channel_count(lm, sigma_only=True, homonuclear=True)
            l4het = level4_channel_count(lm, sigma_only=True, homonuclear=False)
            e4s = four_electron_channel_count_molecular(lm, sigma_only=True, homonuclear=False)
            results[lm] = {
                'l4_homo': l4h,
                'l4_hetero': l4het,
                'e4_sigma': e4s['sigma_channels'],
            }

        # Verify specific counts
        assert results[0]['e4_sigma'] == 1    # (l_max+1)^4 = 1
        assert results[1]['e4_sigma'] == 16   # 2^4 = 16
        assert results[2]['e4_sigma'] == 81   # 3^4 = 81
        assert results[3]['e4_sigma'] == 256  # 4^4 = 256
        assert results[4]['e4_sigma'] == 625  # 5^4 = 625

        # Full M-coupled counts (expensive at l_max >= 3)
        for lm in [0, 1, 2]:
            e4f = four_electron_channel_count_molecular(lm, sigma_only=False, homonuclear=False)
            assert e4f['full_channels'] >= results[lm]['e4_sigma']

        # Full count at l_max=1
        e4f1 = four_electron_channel_count_molecular(1, sigma_only=False, homonuclear=False)
        # Expected: for each (l1,l2,l3,l4) with l_i in {0,1},
        # count M=0 states, sum over all 16 l-sets
        assert e4f1['full_channels'] > 16

        # Full count at l_max=2
        e4f2 = four_electron_channel_count_molecular(2, sigma_only=False, homonuclear=False)
        assert e4f2['full_channels'] > 81

        # Atomic L=0 coupled
        for lm in [0, 1, 2]:
            at = four_electron_channel_count_atomic(lm)
            assert at['coupled_L0'] >= 1

        # Feasibility estimates
        for lm in [0, 1, 2]:
            f = feasibility_estimate(lm, n_basis_alpha=10)
            assert f['angular_dim_sigma'] > 0
            assert f['total_time_sigma_s'] > 0

    def test_print_key_numbers(self):
        """Print key numbers for manual inspection and analysis writing."""
        import sys
        p = lambda *a: print(*a, file=sys.stderr)
        # Channel counts
        p("\n=== CHANNEL COUNT TABLE ===")
        p(f"{'l_max':>5} | {'L4 homo':>8} | {'L4 het':>7} | {'4e sigma':>9} | {'4e full':>8} | {'4e L0-coupled':>13}")
        for lm in [0, 1, 2, 3, 4]:
            l4h = level4_channel_count(lm, sigma_only=True, homonuclear=True)
            l4het = level4_channel_count(lm, sigma_only=True, homonuclear=False)
            e4s = four_electron_channel_count_molecular(lm, sigma_only=True, homonuclear=False)
            if lm <= 2:
                e4f = four_electron_channel_count_molecular(lm, sigma_only=False, homonuclear=False)
                at = four_electron_channel_count_atomic(lm)
                p(f"{lm:>5} | {l4h:>8} | {l4het:>7} | {e4s['sigma_channels']:>9} | {e4f['full_channels']:>8} | {at['coupled_L0']:>13}")
            else:
                p(f"{lm:>5} | {l4h:>8} | {l4het:>7} | {e4s['sigma_channels']:>9} | {'(slow)':>8} | {'(slow)':>13}")

        # Feasibility
        p("\n=== FEASIBILITY TABLE (n_basis=10 per alpha) ===")
        p(f"{'l_max':>5} | {'4e ang dim':>10} | {'time/rho(s)':>11} | {'total(s)':>10} | {'L4 dim':>7} | {'L4 total(s)':>11} | {'ratio':>6}")
        for lm in [0, 1, 2, 3, 4]:
            f = feasibility_estimate(lm, n_basis_alpha=10)
            ratio = f['angular_dim_sigma'] / max(1, f['level4_dim'])
            p(f"{lm:>5} | {f['angular_dim_sigma']:>10} | {f['time_per_rho_sigma_s']:>11.2e} | {f['total_time_sigma_s']:>10.2e} | {f['level4_dim']:>7} | {f['level4_total_s']:>11.2e} | {ratio:>6.0f}x")

        # SO(12) spectrum
        p("\n=== SO(12) CASIMIR SPECTRUM ===")
        p(f"{'nu':>3} | {'mu_free':>8} | {'l_eff':>6} | {'degeneracy':>11}")
        for item in so12_casimir_spectrum(6):
            p(f"{item['nu']:>3} | {item['mu_free']:>8.1f} | {item['l_eff']:>6.1f} | {item['degeneracy']:>11}")

        # S_4 irreps
        p("\n=== S_4 IRREPS ===")
        for name, dim in s4_irrep_dimensions().items():
            p(f"  {name}: dim={dim}")

        # Verify something to make this a real test
        assert True
