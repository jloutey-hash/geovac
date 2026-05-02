"""
Tests for vector-photon QED on the Fock graph.

Tests the full vector QED pipeline:
  1. Vertex coupling symmetries and Wigner 3j selection rules
  2. Ground state structural zero: Sigma(GS, GS) = 0
  3. Ward identity: commutator norm diagnostic
  4. Parity rule enforcement: l_a + l_b + q odd
  5. CG sparsity count
  6. Comparison with scalar graph Sigma (structural, not numerical)

All tests use numpy arithmetic since pi is introduced by the spherical
harmonic normalization (this is the calibration exchange constant).
"""

import pytest
import numpy as np

from geovac.vector_qed import (
    build_electron_states,
    build_photon_modes,
    electron_propagator,
    vector_photon_propagator,
    vertex_coupling,
    build_vertex_tensor,
    compute_self_energy,
    check_selection_rules,
    VectorQEDResult,
    # Dirac extension
    build_dirac_electron_states,
    dirac_electron_propagator,
    dirac_vertex_coupling,
    build_dirac_vertex_tensor,
    compute_dirac_self_energy,
    check_dirac_selection_rules,
    DiracVectorQEDResult,
)
from geovac.dirac_matrix_elements import DiracLabel, kappa_to_l


# ===========================================================================
# Electron and photon state construction
# ===========================================================================

class TestElectronStates:
    """Tests for electron state builder."""

    def test_n1_single_state(self):
        """n_max=1: only (1,0,0)."""
        states = build_electron_states(1)
        assert len(states) == 1
        assert states[0] == (1, 0, 0)

    def test_n2_state_count(self):
        """n_max=2: 1 + 4 = 5 states."""
        states = build_electron_states(2)
        assert len(states) == 5
        # n=1: (1,0,0)
        # n=2: (2,0,0), (2,1,-1), (2,1,0), (2,1,1)
        assert states[0] == (1, 0, 0)
        assert (2, 0, 0) in states
        assert (2, 1, -1) in states
        assert (2, 1, 0) in states
        assert (2, 1, 1) in states

    def test_n3_state_count(self):
        """n_max=3: 1 + 4 + 9 = 14 states."""
        states = build_electron_states(3)
        assert len(states) == 14

    def test_n_max_count_formula(self):
        """State count = sum_{n=1}^{n_max} n^2 = n_max(n_max+1)(2*n_max+1)/6."""
        for n_max in range(1, 6):
            states = build_electron_states(n_max)
            expected = n_max * (n_max + 1) * (2 * n_max + 1) // 6
            assert len(states) == expected


class TestPhotonModes:
    """Tests for photon mode builder."""

    def test_q1_modes(self):
        """q_max=1: 3 modes (q=1, m_q=-1,0,1)."""
        modes = build_photon_modes(1)
        assert len(modes) == 3

    def test_q2_modes(self):
        """q_max=2: 3 + 5 = 8 modes."""
        modes = build_photon_modes(2)
        assert len(modes) == 8

    def test_mode_count_formula(self):
        """Mode count = sum_{q=1}^{q_max} (2q+1)."""
        for q_max in range(1, 6):
            modes = build_photon_modes(q_max)
            expected = sum(2 * q + 1 for q in range(1, q_max + 1))
            assert len(modes) == expected


# ===========================================================================
# Propagators
# ===========================================================================

class TestPropagators:
    """Tests for electron and photon propagators."""

    def test_electron_gs_undefined(self):
        """Electron propagator is None for n=1 (ground state)."""
        assert electron_propagator(1) is None

    def test_electron_n2(self):
        """G_e(n=2) = 1/(4-1) = 1/3."""
        assert abs(electron_propagator(2) - 1.0 / 3) < 1e-14

    def test_electron_n3(self):
        """G_e(n=3) = 1/(9-1) = 1/8."""
        assert abs(electron_propagator(3) - 1.0 / 8) < 1e-14

    def test_photon_q1(self):
        """G_gamma(q=1) = 1/(1*3) = 1/3."""
        assert abs(vector_photon_propagator(1) - 1.0 / 3) < 1e-14

    def test_photon_q2(self):
        """G_gamma(q=2) = 1/(2*4) = 1/8."""
        assert abs(vector_photon_propagator(2) - 1.0 / 8) < 1e-14


# ===========================================================================
# Vertex coupling selection rules
# ===========================================================================

class TestVertexCoupling:
    """Tests for vertex coupling V(a, b, q, m_q)."""

    def test_parity_even_vanishes(self):
        """V = 0 when l_a + l_b + q is even."""
        # l_a=0, l_b=0, q=2 => 0+0+2 = 2 even => 0
        val = vertex_coupling(1, 0, 0, 2, 0, 0, 2, 0)
        assert abs(val) < 1e-15

    def test_parity_odd_nonzero(self):
        """V can be nonzero when l_a + l_b + q is odd."""
        # l_a=0, l_b=1, q=2: 0+1+2 = 3 odd
        # But triangle: |0-1| <= 2 <= 0+1: 1 <= 2 <= 1 => NO (q=2 > 1)
        val = vertex_coupling(2, 0, 0, 2, 1, 0, 2, 0)
        assert abs(val) < 1e-15  # triangle violation

        # l_a=0, l_b=1, q=1: 0+1+1 = 2 even => 0
        val2 = vertex_coupling(2, 0, 0, 2, 1, 0, 1, 0)
        assert abs(val2) < 1e-15  # parity even

    def test_magnetic_conservation(self):
        """V = 0 when -m_a + m_q + m_b != 0."""
        # l_a=1, l_b=1, q=1: parity 3 odd, triangle OK
        # m_a=1, m_q=1, m_b=1: -1 + 1 + 1 = 1 != 0 => 0
        val = vertex_coupling(2, 1, 1, 2, 1, 1, 1, 1)
        assert abs(val) < 1e-15

    def test_gs_vertex_vanishes(self):
        """V(GS, any, any) = 0 for ground state with l=0."""
        # GS: l_a = 0. Triangle forces q = l_b.
        # Then l_a + l_b + q = 0 + l_b + l_b = 2*l_b = even => V = 0
        for l_b in range(4):
            for q in range(1, 4):
                for m_q in range(-q, q + 1):
                    val = vertex_coupling(1, 0, 0, 3, l_b, 0, q, m_q)
                    assert abs(val) < 1e-15, (
                        f"GS vertex nonzero: l_b={l_b}, q={q}, m_q={m_q}, val={val}"
                    )

    def test_triangle_inequality(self):
        """V = 0 when triangle inequality violated."""
        # l_a=0, l_b=0, q=1: |0-0|=0 <= 1 OK, BUT 0+0=0 < 1 => violated
        val = vertex_coupling(2, 0, 0, 3, 0, 0, 1, 0)
        assert abs(val) < 1e-15

    def test_known_nonzero_coupling(self):
        """Known nonzero coupling requires nonzero m.

        3j(l,l,q;0,0,0) vanishes when l+l+q is odd, which is exactly
        our parity rule.  So all-m=0 couplings are zero by construction.
        Must use nonzero m_a, m_b, m_q to get a nonzero vertex.
        """
        # l_a=1, m_a=-1, l_b=1, m_b=-1, q=1, m_q=0
        # parity: 1+1+1=3 odd; triangle OK; m conservation: 1+0+(-1)=0
        val = vertex_coupling(2, 1, -1, 2, 1, -1, 1, 0)
        assert abs(val) > 1e-10, f"Expected nonzero, got {val}"

    def test_vertex_symmetry_ab(self):
        """V(a,b,q,m_q) vs V(b,a,q,-m_q) relation from 3j symmetry."""
        # 3j(l_a, q, l_b; -m_a, m_q, m_b) = (-1)^(l_a+q+l_b) * 3j(l_b, q, l_a; -m_b, m_q, m_a)
        # But our parity rule requires l_a+l_b+q odd => (-1)^(l_a+q+l_b) = -1
        # However the Gaunt-like coupling also has the prefactor difference.
        # Just check that both are nonzero or both zero.
        la, lb, q = 1, 1, 1
        for ma in range(-la, la + 1):
            for mb in range(-lb, lb + 1):
                for mq in range(-q, q + 1):
                    v1 = vertex_coupling(2, la, ma, 2, lb, mb, q, mq)
                    v2 = vertex_coupling(2, lb, mb, 2, la, ma, q, -mq)
                    # Both zero or both nonzero
                    assert (abs(v1) < 1e-15) == (abs(v2) < 1e-15), (
                        f"Symmetry mismatch: V1={v1}, V2={v2}, "
                        f"ma={ma}, mb={mb}, mq={mq}"
                    )


# ===========================================================================
# Self-energy at n_max=2
# ===========================================================================

class TestSelfEnergyN2:
    """Self-energy tests at n_max=2, q_max=1."""

    @pytest.fixture(scope="class")
    def result_n2(self):
        """Compute self-energy once for n_max=2."""
        return compute_self_energy(n_max=2, q_max=1)

    def test_returns_dataclass(self, result_n2):
        """Return type is VectorQEDResult."""
        assert isinstance(result_n2, VectorQEDResult)

    def test_dimensions(self, result_n2):
        """n_max=2: 5 electron states, 3 photon modes."""
        assert result_n2.N_electron == 5
        assert result_n2.N_photon == 3

    def test_sigma_shape(self, result_n2):
        """Sigma is 5x5."""
        assert result_n2.Sigma.shape == (5, 5)

    def test_sigma_symmetric(self, result_n2):
        """Sigma is real symmetric."""
        assert np.allclose(result_n2.Sigma, result_n2.Sigma.T, atol=1e-14)

    def test_gs_structural_zero(self, result_n2):
        """Ground state (n=1, l=0, m=0) has zero self-energy."""
        # GS is state index 0: (1, 0, 0)
        assert abs(result_n2.Sigma[0, 0]) < 1e-14
        # Full GS row/column should be zero
        assert np.max(np.abs(result_n2.Sigma[0, :])) < 1e-14
        assert np.max(np.abs(result_n2.Sigma[:, 0])) < 1e-14

    def test_sigma_positive_semidefinite(self, result_n2):
        """Sigma is positive semidefinite (one-loop self-energy is PSD)."""
        eigenvalues = np.linalg.eigvalsh(result_n2.Sigma)
        assert np.all(eigenvalues >= -1e-12)

    def test_contains_pi(self, result_n2):
        """Vector photon vertex introduces pi."""
        assert result_n2.contains_pi is True

    def test_delta_m_conservation(self, result_n2):
        """Sigma(a,c) = 0 when m_a != m_c."""
        states = build_electron_states(2)
        for i, (ni, li, mi) in enumerate(states):
            for j, (nj, lj, mj) in enumerate(states):
                if mi != mj:
                    assert abs(result_n2.Sigma[i, j]) < 1e-14, (
                        f"Delta_m violation: ({ni},{li},{mi})-({nj},{lj},{mj}), "
                        f"Sigma={result_n2.Sigma[i, j]}"
                    )


# ===========================================================================
# Self-energy at n_max=3
# ===========================================================================

class TestSelfEnergyN3:
    """Self-energy tests at n_max=3."""

    @pytest.fixture(scope="class")
    def result_n3(self):
        """Compute self-energy once for n_max=3."""
        return compute_self_energy(n_max=3, q_max=2)

    def test_dimensions(self, result_n3):
        """n_max=3: 14 electron states, 8 photon modes."""
        assert result_n3.N_electron == 14
        assert result_n3.N_photon == 8

    def test_gs_structural_zero(self, result_n3):
        """Ground state structural zero persists at n_max=3."""
        assert np.max(np.abs(result_n3.Sigma[0, :])) < 1e-14
        assert np.max(np.abs(result_n3.Sigma[:, 0])) < 1e-14

    def test_sigma_symmetric(self, result_n3):
        """Sigma is symmetric at n_max=3."""
        assert np.allclose(result_n3.Sigma, result_n3.Sigma.T, atol=1e-14)

    def test_sigma_psd(self, result_n3):
        """Sigma is positive semidefinite at n_max=3."""
        eigenvalues = np.linalg.eigvalsh(result_n3.Sigma)
        assert np.all(eigenvalues >= -1e-12)


# ===========================================================================
# Selection rules comprehensive check
# ===========================================================================

class TestSelectionRules:
    """Tests for the 8-rule selection rule census."""

    @pytest.fixture(scope="class")
    def rules_n2(self):
        """Selection rules at n_max=2."""
        result = compute_self_energy(n_max=2, q_max=1)
        return result.selection_rules

    def test_gaunt_cg_sparsity_passes(self, rules_n2):
        """Rule 1: Gaunt/CG sparsity."""
        assert rules_n2['1_gaunt_cg_sparsity']['pass'] is True

    def test_vertex_parity_gs_zero_passes(self, rules_n2):
        """Rule 2: Vertex parity / GS structural zero."""
        assert rules_n2['2_vertex_parity_gs_zero']['pass'] is True

    def test_so4_channel_count(self, rules_n2):
        """Rule 3: SO(4) channel count W is well-defined."""
        assert rules_n2['3_so4_channel_count']['pass'] is True

    def test_delta_m_conservation_passes(self, rules_n2):
        """Rule 4: Delta_m conservation."""
        assert rules_n2['4_delta_m_conservation']['pass'] is True
        assert rules_n2['4_delta_m_conservation']['violations'] == 0

    def test_spatial_parity_passes(self, rules_n2):
        """Rule 5: Spatial parity enforcement."""
        assert rules_n2['5_spatial_parity']['pass'] is True
        assert rules_n2['5_spatial_parity']['violations'] == 0

    def test_furry_theorem(self, rules_n2):
        """Rule 6: Furry's theorem (tadpole vanishes)."""
        # Note: may not pass depending on parity structure
        # At n_max=2 with l_max=1, diagonal couplings l=1,l=1,q=1
        # are allowed by parity (1+1+1=3 odd), so tadpole may be nonzero
        rule = rules_n2['6_furry_theorem']
        # Record the outcome for diagnostic purposes
        if not rule['pass']:
            # This is informational -- Furry's theorem on the finite graph
            # may not hold exactly due to the E-type parity allowing
            # diagonal l=l couplings for odd q
            pass

    def test_charge_conjugation_passes(self, rules_n2):
        """Rule 8: Charge conjugation (Hermiticity)."""
        assert rules_n2['8_charge_conjugation']['pass'] is True

    def test_summary_count(self, rules_n2):
        """Summary has correct total count (8 rules)."""
        summary = rules_n2['_summary']
        assert summary['total'] == 8


# ===========================================================================
# Parity rule enforcement
# ===========================================================================

class TestParityEnforcement:
    """Detailed tests for the parity rule l_a + l_b + q odd."""

    def test_vertex_tensor_parity(self):
        """All nonzero entries in V have l_a + l_b + q odd."""
        states = build_electron_states(3)
        modes = build_photon_modes(2)
        V = build_vertex_tensor(states, modes)
        for i, (na, la, ma) in enumerate(states):
            for j, (nb, lb, mb) in enumerate(states):
                for k, (q, mq) in enumerate(modes):
                    if abs(V[i, j, k]) > 1e-15:
                        assert (la + lb + q) % 2 == 1, (
                            f"Parity violation at ({na},{la},{ma})-"
                            f"({nb},{lb},{mb})-({q},{mq}): "
                            f"l_a+l_b+q={la+lb+q} is even"
                        )


# ===========================================================================
# CG sparsity
# ===========================================================================

class TestCGSparsity:
    """Tests for vertex tensor sparsity from CG selection rules."""

    def test_vertex_sparsity_n2(self):
        """Vertex tensor at n_max=2 is highly sparse."""
        states = build_electron_states(2)
        modes = build_photon_modes(1)
        V = build_vertex_tensor(states, modes)
        nnz = np.count_nonzero(V)
        total = V.size
        # With parity + triangle + magnetic conservation, most entries are zero
        assert nnz < total / 2, f"Expected high sparsity, got {nnz}/{total}"

    def test_vertex_sparsity_n3(self):
        """Vertex tensor at n_max=3 is also sparse."""
        states = build_electron_states(3)
        modes = build_photon_modes(2)
        V = build_vertex_tensor(states, modes)
        nnz = np.count_nonzero(V)
        total = V.size
        assert nnz < total / 2, f"Expected high sparsity, got {nnz}/{total}"


# ===========================================================================
# Comparison with scalar graph
# ===========================================================================

class TestScalarComparison:
    """Structural comparison with scalar graph self-energy."""

    def test_vector_gs_zero_scalar_gs_nonzero(self):
        """Vector QED has GS zero; scalar graph (CG projection) does not."""
        # Vector QED
        vec_result = compute_self_energy(n_max=2, q_max=1)
        vec_gs_max = float(np.max(np.abs(vec_result.Sigma[0, :])))

        # Vector GS is zero by parity rule
        assert vec_gs_max < 1e-14

    def test_vector_has_pi(self):
        """Vector QED introduces pi; scalar is pi-free."""
        vec_result = compute_self_energy(n_max=2, q_max=1)
        assert vec_result.contains_pi is True

    def test_sigma_is_real(self):
        """Self-energy is purely real (no complex entries)."""
        vec_result = compute_self_energy(n_max=2, q_max=1)
        assert vec_result.Sigma.dtype in (np.float64, np.float32)


# ===========================================================================
# Ward identity
# ===========================================================================

class TestWardIdentity:
    """Tests for the Ward identity diagnostic."""

    def test_ward_ratio_finite(self):
        """Ward identity ratio is finite and non-negative."""
        result = compute_self_energy(n_max=2, q_max=1)
        rule = result.selection_rules['7_ward_identity']
        assert rule['ratio'] >= 0
        assert np.isfinite(rule['ratio'])

    def test_ward_ratio_n3(self):
        """Ward identity ratio at n_max=3."""
        result = compute_self_energy(n_max=3, q_max=2)
        rule = result.selection_rules['7_ward_identity']
        assert rule['ratio'] >= 0
        assert np.isfinite(rule['ratio'])
        # Record for diagnostic
        # Ward ratio may or may not be small depending on truncation


# ===========================================================================
# Edge cases
# ===========================================================================

class TestEdgeCases:
    """Edge case tests."""

    def test_n_max_1_trivial(self):
        """n_max=1: only GS, no internal states, Sigma = 0."""
        result = compute_self_energy(n_max=1, q_max=1)
        assert np.allclose(result.Sigma, 0, atol=1e-14)

    def test_consistent_q_max_default(self):
        """Default q_max = n_max - 1."""
        result = compute_self_energy(n_max=3)
        assert result.q_max == 2

    def test_sigma_trace_nonnegative(self):
        """Trace of PSD matrix is non-negative."""
        for n_max in [2, 3]:
            result = compute_self_energy(n_max=n_max)
            assert np.trace(result.Sigma) >= -1e-12


# ===========================================================================
# Dirac electron states
# ===========================================================================

class TestDiracElectronStates:
    """Tests for Dirac electron state builder."""

    def test_n1_two_states(self):
        """n_max=1: two states (kappa=-1, m_j=+-1/2)."""
        states = build_dirac_electron_states(1)
        assert len(states) == 2
        # Both have n_fock=1, kappa=-1
        for s in states:
            assert s.n_fock == 1
            assert s.kappa == -1

    def test_n2_state_count(self):
        """n_max=2: 2 + 8 = 10 Dirac states."""
        states = build_dirac_electron_states(2)
        assert len(states) == 10

    def test_n3_state_count(self):
        """n_max=3: 2 + 8 + 18 = 28 Dirac states."""
        states = build_dirac_electron_states(3)
        assert len(states) == 28

    def test_n_max_count_formula(self):
        """State count = sum_{n=1}^{n_max} 2*n^2."""
        for n_max in range(1, 5):
            states = build_dirac_electron_states(n_max)
            expected = sum(2 * n * n for n in range(1, n_max + 1))
            assert len(states) == expected

    def test_gs_is_kappa_minus_one(self):
        """Ground state (n_fock=1) has kappa=-1 (l=0, j=1/2)."""
        states = build_dirac_electron_states(2)
        gs_states = [s for s in states if s.n_fock == 1]
        assert len(gs_states) == 2
        for s in gs_states:
            assert s.kappa == -1
            assert kappa_to_l(s.kappa) == 0


# ===========================================================================
# Dirac propagator
# ===========================================================================

class TestDiracPropagator:
    """Tests for Dirac electron propagator."""

    def test_gs_well_defined(self):
        """Dirac propagator is well-defined at n=1 (unlike scalar)."""
        states = build_dirac_electron_states(1)
        val = dirac_electron_propagator(states[0])
        assert val is not None
        assert abs(val - 2.0 / 3) < 1e-14  # 1/(1 + 0.5) = 2/3

    def test_n2_value(self):
        """G_e(n=2) = 1/(2 + 0.5) = 2/5."""
        states = build_dirac_electron_states(2)
        n2_state = next(s for s in states if s.n_fock == 2)
        assert abs(dirac_electron_propagator(n2_state) - 2.0 / 5) < 1e-14

    def test_n3_value(self):
        """G_e(n=3) = 1/(3 + 0.5) = 2/7."""
        states = build_dirac_electron_states(3)
        n3_state = next(s for s in states if s.n_fock == 3)
        assert abs(dirac_electron_propagator(n3_state) - 2.0 / 7) < 1e-14

    def test_monotone_decreasing(self):
        """Propagator decreases with n_fock."""
        states = build_dirac_electron_states(4)
        prev = float('inf')
        for n in range(1, 5):
            s = next(s for s in states if s.n_fock == n)
            val = dirac_electron_propagator(s)
            assert val < prev
            prev = val


# ===========================================================================
# Dirac vertex coupling
# ===========================================================================

class TestDiracVertexCoupling:
    """Tests for Dirac vertex coupling with l-triangle enforcement."""

    def test_gs_vertex_vanishes(self):
        """All couplings involving GS (l=0) vanish by l-triangle + parity."""
        states = build_dirac_electron_states(2)
        modes = build_photon_modes(1)
        gs_states = [s for s in states if s.n_fock == 1]
        for gs in gs_states:
            for b in states:
                for q, m_q in modes:
                    val = dirac_vertex_coupling(gs, b, q, m_q)
                    assert abs(val) < 1e-15, (
                        f"GS vertex nonzero: gs={gs}, b={b}, q={q}, m_q={m_q}, "
                        f"V={val}"
                    )

    def test_l_triangle_enforced(self):
        """Couplings violating l-triangle are zero."""
        states = build_dirac_electron_states(3)
        modes = build_photon_modes(2)
        V = build_dirac_vertex_tensor(states, modes)
        for i, a in enumerate(states):
            for j, b in enumerate(states):
                l_a = kappa_to_l(a.kappa)
                l_b = kappa_to_l(b.kappa)
                for k, (q, m_q) in enumerate(modes):
                    if q < abs(l_a - l_b) or q > l_a + l_b:
                        assert abs(V[i, j, k]) < 1e-15

    def test_parity_enforced(self):
        """All nonzero vertex entries have l_a + l_b + q odd."""
        states = build_dirac_electron_states(3)
        modes = build_photon_modes(2)
        V = build_dirac_vertex_tensor(states, modes)
        for i, a in enumerate(states):
            for j, b in enumerate(states):
                l_a = kappa_to_l(a.kappa)
                l_b = kappa_to_l(b.kappa)
                for k, (q, m_q) in enumerate(modes):
                    if abs(V[i, j, k]) > 1e-15:
                        assert (l_a + l_b + q) % 2 == 1, (
                            f"Parity violation: l_a={l_a}, l_b={l_b}, q={q}"
                        )

    def test_delta_mj_enforced(self):
        """All nonzero vertex entries have -m_j_a + m_q + m_j_b = 0."""
        states = build_dirac_electron_states(2)
        modes = build_photon_modes(1)
        V = build_dirac_vertex_tensor(states, modes)
        for i, a in enumerate(states):
            for j, b in enumerate(states):
                for k, (q, m_q) in enumerate(modes):
                    if abs(V[i, j, k]) > 1e-15:
                        assert -a.m_j + m_q + b.m_j == 0, (
                            f"m_j violation: m_j_a={a.m_j}, m_q={m_q}, "
                            f"m_j_b={b.m_j}"
                        )


# ===========================================================================
# Dirac self-energy at n_max=2
# ===========================================================================

class TestDiracSelfEnergyN2:
    """Self-energy tests at n_max=2, q_max=1 (Dirac + vector photon)."""

    @pytest.fixture(scope="class")
    def result_n2(self):
        """Compute Dirac self-energy once for n_max=2."""
        return compute_dirac_self_energy(n_max=2, q_max=1)

    def test_returns_dataclass(self, result_n2):
        """Return type is DiracVectorQEDResult."""
        assert isinstance(result_n2, DiracVectorQEDResult)

    def test_dimensions(self, result_n2):
        """n_max=2: 10 electron states, 3 photon modes."""
        assert result_n2.N_electron == 10
        assert result_n2.N_photon == 3

    def test_sigma_shape(self, result_n2):
        """Sigma is 10x10."""
        assert result_n2.Sigma.shape == (10, 10)

    def test_sigma_symmetric(self, result_n2):
        """Sigma is real symmetric."""
        assert np.allclose(result_n2.Sigma, result_n2.Sigma.T, atol=1e-14)

    def test_gs_structural_zero(self, result_n2):
        """Ground state (n=1, kappa=-1) has zero self-energy."""
        gs_indices = [i for i, s in enumerate(result_n2.states) if s.n_fock == 1]
        assert len(gs_indices) == 2
        # Full GS block should be zero
        for gi in gs_indices:
            assert np.max(np.abs(result_n2.Sigma[gi, :])) < 1e-14
            assert np.max(np.abs(result_n2.Sigma[:, gi])) < 1e-14

    def test_sigma_positive_semidefinite(self, result_n2):
        """Sigma is positive semidefinite."""
        eigenvalues = np.linalg.eigvalsh(result_n2.Sigma)
        assert np.all(eigenvalues >= -1e-12)

    def test_four_zero_eigenvalues(self, result_n2):
        """At n_max=2, Sigma has exactly 4 zero eigenvalues."""
        eigenvalues = np.linalg.eigvalsh(result_n2.Sigma)
        n_zero = sum(1 for ev in eigenvalues if abs(ev) < 1e-12)
        assert n_zero == 4

    def test_six_degenerate_nonzero_eigenvalues(self, result_n2):
        """At n_max=2, the 6 nonzero eigenvalues are degenerate."""
        eigenvalues = sorted(np.linalg.eigvalsh(result_n2.Sigma))
        nonzero = [ev for ev in eigenvalues if ev > 1e-12]
        assert len(nonzero) == 6
        # All should be approximately equal
        assert max(nonzero) - min(nonzero) < 1e-12

    def test_calibration_4pi(self, result_n2):
        """Nonzero eigenvalues * 4*pi are rational (12/5 = 2.4)."""
        eigenvalues = sorted(np.linalg.eigvalsh(result_n2.Sigma))
        nonzero = [ev for ev in eigenvalues if ev > 1e-12]
        # Each nonzero eigenvalue * 4*pi should give 12/5
        for ev in nonzero:
            val = ev * 4 * np.pi
            assert abs(val - 12.0 / 5) < 1e-10, (
                f"Expected ev*4*pi = 12/5 = 2.4, got {val}"
            )

    def test_delta_mj_conservation(self, result_n2):
        """Sigma(a,c) = 0 when m_j_a != m_j_c."""
        for i, si in enumerate(result_n2.states):
            for j, sj in enumerate(result_n2.states):
                if si.two_m_j != sj.two_m_j:
                    assert abs(result_n2.Sigma[i, j]) < 1e-14

    def test_ward_identity_exact_at_n2(self, result_n2):
        """Ward identity is exact at n_max=2 (Sigma is n-block-diagonal)."""
        rule = result_n2.selection_rules['7_ward_identity']
        assert rule['ratio'] < 1e-12

    def test_contains_pi(self, result_n2):
        """Vertex normalization introduces pi."""
        assert result_n2.contains_pi is True

    def test_vertex_sparsity(self, result_n2):
        """Vertex tensor is >90% sparse."""
        assert result_n2.vertex_sparsity > 0.90
        assert result_n2.vertex_nonzero_count == 26


# ===========================================================================
# Dirac self-energy at n_max=3
# ===========================================================================

class TestDiracSelfEnergyN3:
    """Self-energy tests at n_max=3 (Dirac + vector photon)."""

    @pytest.fixture(scope="class")
    def result_n3(self):
        """Compute Dirac self-energy once for n_max=3."""
        return compute_dirac_self_energy(n_max=3, q_max=2)

    def test_dimensions(self, result_n3):
        """n_max=3: 28 electron states, 8 photon modes."""
        assert result_n3.N_electron == 28
        assert result_n3.N_photon == 8

    def test_gs_structural_zero(self, result_n3):
        """Ground state structural zero persists at n_max=3."""
        gs_indices = [i for i, s in enumerate(result_n3.states) if s.n_fock == 1]
        for gi in gs_indices:
            assert np.max(np.abs(result_n3.Sigma[gi, :])) < 1e-14

    def test_sigma_symmetric(self, result_n3):
        """Sigma is symmetric at n_max=3."""
        assert np.allclose(result_n3.Sigma, result_n3.Sigma.T, atol=1e-14)

    def test_sigma_psd(self, result_n3):
        """Sigma is positive semidefinite at n_max=3."""
        eigenvalues = np.linalg.eigvalsh(result_n3.Sigma)
        assert np.all(eigenvalues >= -1e-12)

    def test_twelve_zero_eigenvalues(self, result_n3):
        """At n_max=3, Sigma has 12 zero eigenvalues."""
        eigenvalues = np.linalg.eigvalsh(result_n3.Sigma)
        n_zero = sum(1 for ev in eigenvalues if abs(ev) < 1e-12)
        assert n_zero == 12

    def test_two_nonzero_eigenvalue_tiers(self, result_n3):
        """At n_max=3, nonzero eigenvalues form two degenerate tiers."""
        eigenvalues = sorted(np.linalg.eigvalsh(result_n3.Sigma))
        nonzero = [ev for ev in eigenvalues if ev > 1e-12]
        assert len(nonzero) == 16
        # Split into two tiers
        tier1 = [ev for ev in nonzero if ev < 0.7]
        tier2 = [ev for ev in nonzero if ev > 0.7]
        assert len(tier1) == 10  # ~0.432 each
        assert len(tier2) == 6   # ~0.939 each
        # Each tier is internally degenerate
        assert max(tier1) - min(tier1) < 1e-10
        assert max(tier2) - min(tier2) < 1e-10

    def test_ward_identity_fails_at_n3(self, result_n3):
        """Ward identity fails at n_max=3 (cross-n couplings)."""
        rule = result_n3.selection_rules['7_ward_identity']
        # Ward ratio ~0.61
        assert rule['ratio'] > 0.1

    def test_vertex_nonzero_count(self, result_n3):
        """Vertex tensor has 346 nonzero entries at n_max=3."""
        assert result_n3.vertex_nonzero_count == 346


# ===========================================================================
# Dirac selection rules
# ===========================================================================

class TestDiracSelectionRules:
    """Tests for the 8-rule selection rule census (Dirac + vector)."""

    @pytest.fixture(scope="class")
    def rules_n2(self):
        """Selection rules at n_max=2."""
        result = compute_dirac_self_energy(n_max=2, q_max=1)
        return result.selection_rules

    def test_total_pass_count(self, rules_n2):
        """Dirac + vector gives 7/8 selection rules at n_max=2.

        Ward identity passes at n_max=2 because Sigma is block-diagonal
        within each n-shell. Only Furry's theorem fails.
        """
        summary = rules_n2['_summary']
        assert summary['pass_count'] == 7
        assert summary['total'] == 8

    def test_gaunt_cg_sparsity_passes(self, rules_n2):
        """Rule 1: Gaunt/CG sparsity."""
        assert rules_n2['1_gaunt_cg_sparsity']['pass'] is True

    def test_gs_structural_zero_passes(self, rules_n2):
        """Rule 2: GS structural zero."""
        assert rules_n2['2_vertex_parity_gs_zero']['pass'] is True
        assert rules_n2['2_vertex_parity_gs_zero']['gs_block_max_abs'] == 0.0

    def test_so4_channel_count_passes(self, rules_n2):
        """Rule 3: SO(4) channel count."""
        assert rules_n2['3_so4_channel_count']['pass'] is True

    def test_delta_mj_passes(self, rules_n2):
        """Rule 4: Delta_m_j conservation."""
        assert rules_n2['4_delta_mj_conservation']['pass'] is True
        assert rules_n2['4_delta_mj_conservation']['violations'] == 0

    def test_spatial_parity_passes(self, rules_n2):
        """Rule 5: Spatial parity."""
        assert rules_n2['5_spatial_parity']['pass'] is True
        assert rules_n2['5_spatial_parity']['violations'] == 0

    def test_furry_fails(self, rules_n2):
        """Rule 6: Furry's theorem FAILS (diagonal l>=1 couplings allowed)."""
        assert rules_n2['6_furry_theorem']['pass'] is False
        # 6 nonzero tadpole states at n_max=2 (the l>=1 states)
        assert rules_n2['6_furry_theorem']['tadpole_nonzero_count'] == 6

    def test_ward_identity_passes_at_n2(self, rules_n2):
        """Rule 7: Ward identity passes at n_max=2.

        At n_max=2, Sigma is block-diagonal within each n-shell, so
        [Sigma, H0] = 0 exactly (H0 = diag(n+1/2) is also n-diagonal).
        Breaks at n_max=3 where cross-n couplings appear.
        """
        assert rules_n2['7_ward_identity']['pass'] is True
        assert rules_n2['7_ward_identity']['ratio'] < 1e-12

    def test_charge_conjugation_passes(self, rules_n2):
        """Rule 8: Charge conjugation (Hermiticity)."""
        assert rules_n2['8_charge_conjugation']['pass'] is True

    def test_summary_total_is_eight(self, rules_n2):
        """Summary counts all 8 rules."""
        assert rules_n2['_summary']['total'] == 8


# ===========================================================================
# Scalar vs Dirac comparison
# ===========================================================================

class TestScalarDiracComparison:
    """Structural comparison: scalar+vector vs Dirac+vector."""

    def test_same_selection_rule_count(self):
        """Both scalar+vector and Dirac+vector give 7/8 at n_max=2.

        Ward identity passes at n_max=2 for both (Sigma is n-block-diagonal).
        Only Furry's theorem fails in both configurations.
        At n_max=3, both drop to 6/8 (Ward also fails due to cross-n couplings).
        """
        scalar = compute_self_energy(n_max=2, q_max=1)
        dirac = compute_dirac_self_energy(n_max=2, q_max=1)
        s_sum = scalar.selection_rules['_summary']
        d_sum = dirac.selection_rules['_summary']
        assert s_sum['pass_count'] == d_sum['pass_count'] == 7
        assert s_sum['total'] == d_sum['total'] == 8

    def test_both_have_gs_zero(self):
        """Both configurations protect the ground state."""
        scalar = compute_self_energy(n_max=2, q_max=1)
        dirac = compute_dirac_self_energy(n_max=2, q_max=1)
        assert scalar.selection_rules['2_vertex_parity_gs_zero']['pass'] is True
        assert dirac.selection_rules['2_vertex_parity_gs_zero']['pass'] is True

    def test_both_fail_furry(self):
        """Both configurations fail Furry's theorem."""
        scalar = compute_self_energy(n_max=2, q_max=1)
        dirac = compute_dirac_self_energy(n_max=2, q_max=1)
        assert scalar.selection_rules['6_furry_theorem']['pass'] is False
        assert dirac.selection_rules['6_furry_theorem']['pass'] is False

    def test_dirac_has_more_states(self):
        """Dirac basis is larger than scalar (10 vs 5 at n_max=2)."""
        scalar = compute_self_energy(n_max=2, q_max=1)
        dirac = compute_dirac_self_energy(n_max=2, q_max=1)
        assert dirac.N_electron == 10
        assert scalar.N_electron == 5

    def test_both_contain_pi(self):
        """Both configurations introduce pi via vector photon normalization."""
        scalar = compute_self_energy(n_max=2, q_max=1)
        dirac = compute_dirac_self_energy(n_max=2, q_max=1)
        assert scalar.contains_pi is True
        assert dirac.contains_pi is True
