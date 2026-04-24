"""
Tests for the native Dirac lattice on S³.

Verifies:
  - State counts match Camporesi-Higuchi degeneracy 2(n+1)(n+2) per CH level
  - Chirality ℤ₂ pairs states correctly
  - E1 dipole adjacency is symmetric and consistent with selection rules
  - Dirac eigenvalues match |λ_n| = n + 3/2
  - Ratio to scalar lattice is exactly 2 (spinor doubling) at each n-shell
  - Graph is connected for n_max >= 2
"""

import numpy as np
import pytest
from geovac.dirac_lattice import DiracLattice
from geovac.dirac_matrix_elements import DiracLabel, kappa_to_l


class TestStateGeneration:
    """Verify that the Dirac lattice generates the correct states."""

    def test_n_max_1_state_count(self):
        lat = DiracLattice(n_max=1)
        # n=1: l=0 only, κ=-1 only (j=1/2), m_j = ±1/2 → 2 states
        assert lat.num_states == 2

    def test_n_max_2_state_count(self):
        lat = DiracLattice(n_max=2)
        # n=1: 2 states (s_{1/2})
        # n=2: l=0→κ=-1 (2 states), l=1→κ=-2,+1 (4+2=6 states) → 8 states
        # Total: 10
        assert lat.num_states == 10

    def test_n_max_3_state_count(self):
        lat = DiracLattice(n_max=3)
        # n=1: 2, n=2: 8, n=3: 18 → 28
        # Verify: Σ 2n² = 2·(1+4+9) = 28
        assert lat.num_states == 28

    def test_state_count_formula(self):
        """Total states = 2 · Σ n² = n(n+1)(2n+1)/3."""
        for n_max in range(1, 7):
            lat = DiracLattice(n_max=n_max)
            expected = n_max * (n_max + 1) * (2 * n_max + 1) // 3
            assert lat.num_states == expected, f"n_max={n_max}: {lat.num_states} != {expected}"

    def test_spinor_doubling(self):
        """Dirac lattice has exactly 2× the states of the scalar lattice."""
        for n_max in range(1, 6):
            lat = DiracLattice(n_max=n_max)
            assert lat.num_states == 2 * lat.scalar_degeneracy

    def test_per_shell_count(self):
        """Each n-shell has 2n² Dirac states."""
        lat = DiracLattice(n_max=5)
        shells = lat.shell_structure()
        for n in range(1, 6):
            assert shells[n]["num_states"] == 2 * n * n


class TestChirality:
    """Verify the chirality ℤ₂ structure."""

    def test_chirality_values(self):
        """All chiralities are ±1."""
        lat = DiracLattice(n_max=3)
        assert set(np.unique(lat.chirality)).issubset({-1, 1})

    def test_chirality_sign_convention(self):
        """κ < 0 → χ = +1, κ > 0 → χ = -1."""
        lat = DiracLattice(n_max=3)
        for i, lab in enumerate(lat.labels):
            if lab.kappa < 0:
                assert lat.chirality[i] == +1
            else:
                assert lat.chirality[i] == -1

    def test_chirality_balance_per_shell(self):
        """Within each n-shell, count χ=+1 and χ=-1 states."""
        lat = DiracLattice(n_max=4)
        for n in range(1, 5):
            indices = [i for i, lab in enumerate(lat.labels) if lab.n_fock == n]
            chi_plus = sum(1 for i in indices if lat.chirality[i] == +1)
            chi_minus = sum(1 for i in indices if lat.chirality[i] == -1)
            # n=1: only κ=-1 → all χ=+1 (2 states, 0 χ=-1)
            # n≥2: both κ signs present
            if n == 1:
                assert chi_plus == 2 and chi_minus == 0
            else:
                assert chi_plus > 0 and chi_minus > 0
                # More χ=+1 than χ=-1 because κ=-(l+1) has j=l+1/2 (larger multiplicity)
                assert chi_plus > chi_minus

    def test_chirality_partner_exists(self):
        """States with l ≥ 1 have a chirality partner in the lattice."""
        lat = DiracLattice(n_max=3)
        for lab in lat.labels:
            if lab.l >= 1:
                partner = lat.chirality_partner(lab)
                assert partner is not None, f"No partner for {lab}"

    def test_s_wave_no_partner(self):
        """l=0 states with κ=-1 have no κ>0 partner (would need κ=0)."""
        lat = DiracLattice(n_max=3)
        for lab in lat.labels:
            if lab.l == 0:
                assert lab.kappa == -1
                partner = lat.chirality_partner(lab)
                assert partner is None


class TestAdjacency:
    """Verify E1 dipole adjacency structure."""

    def test_symmetry(self):
        """Adjacency matrix is symmetric."""
        lat = DiracLattice(n_max=3)
        diff = lat.adjacency - lat.adjacency.T
        assert diff.nnz == 0

    def test_no_self_loops(self):
        """Diagonal is zero."""
        lat = DiracLattice(n_max=3)
        assert np.all(lat.adjacency.diagonal() == 0)

    def test_n_max_1_no_edges(self):
        """At n_max=1, only s_{1/2} states exist — no E1 transitions."""
        lat = DiracLattice(n_max=1)
        assert lat.num_edges == 0

    def test_n_max_2_has_edges(self):
        """At n_max=2, E1 transitions between s and p states exist."""
        lat = DiracLattice(n_max=2)
        assert lat.num_edges > 0

    def test_selection_rules(self):
        """Every edge satisfies E1 selection rules."""
        lat = DiracLattice(n_max=3)
        rows, cols = lat.adjacency.nonzero()
        for r, c in zip(rows, cols):
            if r >= c:
                continue
            a, b = lat.labels[r], lat.labels[c]
            assert abs(b.n_fock - a.n_fock) <= 1, "Δn > 1"
            assert abs(b.l - a.l) == 1, "Δl ≠ ±1"
            assert abs(b.j_times_2 - a.j_times_2) <= 2, "|Δj| > 1"
            assert abs(b.two_m_j - a.two_m_j) <= 2, "|Δm_j| > 1"

    def test_connected_n_max_2(self):
        """Graph is connected at n_max ≥ 2."""
        lat = DiracLattice(n_max=2)
        from scipy.sparse.csgraph import connected_components
        n_comp, _ = connected_components(lat.adjacency, directed=False)
        assert n_comp == 1

    def test_connected_n_max_3(self):
        lat = DiracLattice(n_max=3)
        from scipy.sparse.csgraph import connected_components
        n_comp, _ = connected_components(lat.adjacency, directed=False)
        assert n_comp == 1


class TestSpectrum:
    """Verify Dirac eigenvalues match Camporesi-Higuchi."""

    def test_eigenvalue_magnitudes(self):
        """Each state has |λ| = n_fock + 3/2."""
        lat = DiracLattice(n_max=4)
        for i, lab in enumerate(lat.labels):
            expected_mag = lab.n_fock + 1.5
            assert abs(abs(lat.dirac_eigenvalues[i]) - expected_mag) < 1e-14

    def test_eigenvalue_signs(self):
        """Sign of eigenvalue matches chirality."""
        lat = DiracLattice(n_max=4)
        for i in range(lat.num_states):
            if lat.chirality[i] == +1:
                assert lat.dirac_eigenvalues[i] > 0
            else:
                assert lat.dirac_eigenvalues[i] < 0

    def test_unique_eigenvalue_count(self):
        """Each n-shell has one positive and one negative eigenvalue magnitude."""
        lat = DiracLattice(n_max=4)
        mags = sorted(set(abs(lat.dirac_eigenvalues[i]) for i in range(lat.num_states)))
        assert len(mags) == lat.n_max
        for idx, n in enumerate(range(1, lat.n_max + 1)):
            assert abs(mags[idx] - (n + 1.5)) < 1e-14


class TestScalarComparison:
    """Compare Dirac lattice with the scalar GeometricLattice."""

    def test_fock_weights_match_scalar(self):
        """Fock weights -Z/n² are the same physics as the scalar lattice."""
        lat = DiracLattice(n_max=3, nuclear_charge=2)
        weights = lat.fock_weights()
        for i, lab in enumerate(lat.labels):
            expected = -2.0 / (lab.n_fock ** 2)
            assert abs(weights[i] - expected) < 1e-14


class TestS3Mode:
    """Tests for mode='s3' — full S³ spinor harmonics including boundary states."""

    def test_s3_state_count_formula(self):
        """Total S³ Dirac states = 2n(n+1)(n+2)/3."""
        for n_max in range(1, 6):
            lat = DiracLattice(n_max=n_max, mode='s3')
            expected = 2 * n_max * (n_max + 1) * (n_max + 2) // 3
            assert lat.num_states == expected, f"n_max={n_max}"

    def test_s3_per_shell_count(self):
        """Each S³ n-shell has 2n(n+1) Dirac states."""
        lat = DiracLattice(n_max=5, mode='s3')
        shells = lat.shell_structure()
        for n in range(1, 6):
            assert shells[n]["num_states"] == 2 * n * (n + 1)

    def test_s3_exceeds_atomic(self):
        """S³ has 2n more states per shell than atomic (the boundary states)."""
        for n_max in range(1, 5):
            s3 = DiracLattice(n_max=n_max, mode='s3')
            at = DiracLattice(n_max=n_max, mode='atomic')
            assert s3.num_states > at.num_states
            assert s3.num_states - at.num_states == sum(2 * n for n in range(1, n_max + 1))

    def test_s3_boundary_states_negative_chirality(self):
        """Boundary states (l = n_fock) all have κ > 0 → χ = -1."""
        lat = DiracLattice(n_max=3, mode='s3')
        for i, lab in enumerate(lat.labels):
            if hasattr(lab, 'is_boundary') and lab.is_boundary:
                assert lab.kappa > 0
                assert lat.chirality[i] == -1

    def test_s3_eigenvalue_magnitudes(self):
        """Each S³ state has |λ| = n_fock + 3/2."""
        lat = DiracLattice(n_max=4, mode='s3')
        for i, lab in enumerate(lat.labels):
            assert abs(abs(lat.dirac_eigenvalues[i]) - (lab.n_fock + 1.5)) < 1e-14

    def test_s3_adjacency_symmetric(self):
        lat = DiracLattice(n_max=3, mode='s3')
        diff = lat.adjacency - lat.adjacency.T
        assert diff.nnz == 0

    def test_s3_selection_rules(self):
        """Every edge in S³ mode satisfies E1 dipole selection rules."""
        lat = DiracLattice(n_max=3, mode='s3')
        rows, cols = lat.adjacency.nonzero()
        for r, c in zip(rows, cols):
            if r >= c:
                continue
            a, b = lat.labels[r], lat.labels[c]
            assert abs(b.n_fock - a.n_fock) <= 1
            assert abs(b.l - a.l) == 1
            assert abs(b.j_times_2 - a.j_times_2) <= 2
            assert abs(b.two_m_j - a.two_m_j) <= 2

    def test_s3_degeneracy_property(self):
        """s3_dirac_degeneracy property matches actual count."""
        for n_max in range(1, 5):
            lat = DiracLattice(n_max=n_max, mode='s3')
            assert lat.s3_dirac_degeneracy == lat.num_states

    def test_s3_vertex_structural_match(self):
        """State-resolved vertex correction matches shell-summed at machine precision."""
        from geovac.dirac_vertex import validate_state_resolved
        result = validate_state_resolved(n_max=3, mode='s3')
        assert result["structural_match"], f"max rel error: {result['max_rel_error']:.2e}"


class TestRepr:
    def test_repr(self):
        lat = DiracLattice(n_max=2)
        s = repr(lat)
        assert "DiracLattice" in s
        assert "states=10" in s
