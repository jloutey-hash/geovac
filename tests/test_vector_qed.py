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
        """Rule 3: vertex-support channel count matches the closed form.

        De-tautologized 2026-07-03 (was a hardcoded pass=True): the count of
        (l_ext, l_int, q) channels realized in the vertex tensor must equal
        the parity + l-triangle closed form for every (n_ext, n_int) pair.
        At n_max=2, q_max=1 the only allowed channel is (1,1,1) at (2,2).
        """
        rule = rules_n2['3_so4_channel_count']
        assert rule['pass'] is True
        assert rule['mismatches'] == []
        assert rule['w_counts']['(2,2)'] == 1
        assert rule['w_observed']['(2,2)'] == 1
        # every pair other than (2,2) has zero allowed channels
        assert all(v == 0 for k, v in rule['w_counts'].items() if k != '(2,2)')

    def test_delta_m_conservation_passes(self, rules_n2):
        """Rule 4: Delta_m conservation."""
        assert rules_n2['4_delta_m_conservation']['pass'] is True
        assert rules_n2['4_delta_m_conservation']['violations'] == 0

    def test_spatial_parity_passes(self, rules_n2):
        """Rule 5: Spatial parity enforcement."""
        assert rules_n2['5_spatial_parity']['pass'] is True
        assert rules_n2['5_spatial_parity']['violations'] == 0

    def test_furry_theorem(self, rules_n2):
        """Rule 7 (Paper 33 Table 1): Furry's theorem (tadpole vanishes).

        Key renamed 6_furry_theorem -> 7_furry_theorem 2026-07-03 to align
        with the paper's Table 1 numbering (7 = charge conjugation / Furry).
        """
        # Note: may not pass depending on parity structure
        # At n_max=2 with l_max=1, diagonal couplings l=1,l=1,q=1
        # are allowed by parity (1+1+1=3 odd), so tadpole may be nonzero
        rule = rules_n2['7_furry_theorem']
        # Paper 33 Thm 1 (scalar_7_of_8): Furry is the ONE rule that FAILS
        # for scalar electrons (E-type parity allows diagonal l=l couplings
        # for odd q); the Dirac construction recovers it (test_furry_passes).
        # 8th-cert fix 2026-07-02: this body was a vacuous no-assert stub.
        assert rule['pass'] is False, (
            "scalar Furry unexpectedly PASSES -- the 7/8 census and Paper 33 "
            "Thm 1 would both be wrong; investigate before trusting")

    def test_triangle_on_n_passes(self, rules_n2):
        """Rule 8 (Paper 33 Table 1): triangle inequality on SO(4),
        |n_a - n_b| <= q <= n_a + n_b - 2 on the vertex tensor support.

        Replaced 2026-07-03: the old key 8_charge_conjugation checked
        Hermiticity of Sigma, which is symmetric BY CONSTRUCTION for any
        real vertex (tautology). The paper's actual Rule 8 had no check.
        """
        rule = rules_n2['8_triangle_on_n']
        assert rule['pass'] is True
        assert rule['violations'] == 0
        assert rule['nonzero_checked'] > 0  # non-vacuous

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
        """Ward identity (Rule 6, Paper 33 Table 1) ratio is finite and
        non-negative. Key renamed 7_ward_identity -> 6_ward_identity
        2026-07-03 to align with the paper's numbering (6 = Ward)."""
        result = compute_self_energy(n_max=2, q_max=1)
        rule = result.selection_rules['6_ward_identity']
        assert rule['ratio'] >= 0
        assert np.isfinite(rule['ratio'])

    def test_ward_ratio_n3(self):
        """Ward identity ratio at n_max=3."""
        result = compute_self_energy(n_max=3, q_max=2)
        rule = result.selection_rules['6_ward_identity']
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

    def test_six_nonzero_in_three_pairs(self, result_n2):
        """At n_max=2, the 6 nonzero eigenvalues form 3 doubly-degenerate pairs.

        The Dirac spinor phase constraint (diagonal V=0) breaks the original
        6-fold degeneracy into 3 pairs corresponding to distinct kappa channels.
        """
        eigenvalues = sorted(np.linalg.eigvalsh(result_n2.Sigma))
        nonzero = [ev for ev in eigenvalues if ev > 1e-12]
        assert len(nonzero) == 6
        # Three doubly-degenerate pairs
        assert abs(nonzero[0] - nonzero[1]) < 1e-12
        assert abs(nonzero[2] - nonzero[3]) < 1e-12
        assert abs(nonzero[4] - nonzero[5]) < 1e-12
        # Pairs are distinct
        assert nonzero[2] - nonzero[1] > 0.01

    def test_calibration_4pi_trace(self, result_n2):
        """Trace(Sigma) * 4*pi = 176/15 (rational), preserving 1/(4pi) calibration.

        The Dirac spinor phase constraint redistributes eigenvalues but the
        trace retains the 1/(4pi) calibration from the S^2 solid angle
        normalization of vector spherical harmonics.
        """
        trace_val = np.trace(result_n2.Sigma) * 4 * np.pi
        assert abs(trace_val - 176.0 / 15) < 1e-10, (
            f"Expected Tr(Sigma)*4*pi = 176/15 = {176/15:.6f}, got {trace_val}"
        )

    def test_delta_mj_conservation(self, result_n2):
        """Sigma(a,c) = 0 when m_j_a != m_j_c."""
        for i, si in enumerate(result_n2.states):
            for j, sj in enumerate(result_n2.states):
                if si.two_m_j != sj.two_m_j:
                    assert abs(result_n2.Sigma[i, j]) < 1e-14

    def test_ward_identity_exact_at_n2(self, result_n2):
        """Ward identity is exact at n_max=2 (Sigma is n-block-diagonal)."""
        rule = result_n2.selection_rules['6_ward_identity']
        assert rule['ratio'] < 1e-12

    def test_contains_pi(self, result_n2):
        """Vertex normalization introduces pi."""
        assert result_n2.contains_pi is True

    def test_vertex_sparsity(self, result_n2):
        """Vertex tensor is >90% sparse (20 nonzero after diagonal removal)."""
        assert result_n2.vertex_sparsity > 0.90
        assert result_n2.vertex_nonzero_count == 20


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

    def test_six_zero_eigenvalues(self, result_n3):
        """At n_max=3, Sigma has 6 zero eigenvalues (GS block + kernel)."""
        eigenvalues = np.linalg.eigvalsh(result_n3.Sigma)
        n_zero = sum(1 for ev in eigenvalues if abs(ev) < 1e-12)
        assert n_zero == 6

    def test_twentytwo_nonzero_eigenvalues(self, result_n3):
        """At n_max=3, 22 nonzero eigenvalues in 11 doubly-degenerate pairs."""
        eigenvalues = sorted(np.linalg.eigvalsh(result_n3.Sigma))
        nonzero = [ev for ev in eigenvalues if ev > 1e-12]
        assert len(nonzero) == 22
        # All pairs are doubly degenerate (m_j degeneracy)
        for i in range(0, len(nonzero), 2):
            assert abs(nonzero[i] - nonzero[i + 1]) < 1e-10

    def test_ward_identity_fails_at_n3(self, result_n3):
        """Ward identity fails at n_max=3 (cross-n couplings)."""
        rule = result_n3.selection_rules['6_ward_identity']
        # Ward ratio ~0.59
        assert rule['ratio'] > 0.1

    def test_vertex_nonzero_count(self, result_n3):
        """Vertex tensor has 324 nonzero entries at n_max=3 (346 - 22 diagonal)."""
        assert result_n3.vertex_nonzero_count == 324


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
        """Dirac + vector gives 8/8 selection rules at n_max=2.

        The Dirac spinor phase constraint (diagonal V=0) recovers Furry's
        theorem via a kinematic identity of the single-particle Dirac vertex.
        Combined with vector photon angular momentum, all 8 rules pass.
        """
        summary = rules_n2['_summary']
        assert summary['pass_count'] == 8
        assert summary['total'] == 8

    def test_gaunt_cg_sparsity_passes(self, rules_n2):
        """Rule 1: Gaunt/CG sparsity."""
        assert rules_n2['1_gaunt_cg_sparsity']['pass'] is True

    def test_gs_structural_zero_passes(self, rules_n2):
        """Rule 2: GS structural zero."""
        assert rules_n2['2_vertex_parity_gs_zero']['pass'] is True
        assert rules_n2['2_vertex_parity_gs_zero']['gs_block_max_abs'] == 0.0

    def test_so4_channel_count_passes(self, rules_n2):
        """Rule 3: vertex-support (kappa_ext, kappa_int, q) channel count
        matches the l-parity + l-triangle + j-triangle closed form.

        De-tautologized 2026-07-03 (was a hardcoded pass=True).
        """
        rule = rules_n2['3_so4_channel_count']
        assert rule['pass'] is True
        assert rule['mismatches'] == []
        # non-vacuous: the (2,2) pair carries realized channels at n_max=2
        assert rule['w_observed']['(2,2)'] > 0
        assert rule['w_counts']['(2,2)'] == rule['w_observed']['(2,2)']

    def test_delta_mj_passes(self, rules_n2):
        """Rule 4: Delta_m_j conservation."""
        assert rules_n2['4_delta_mj_conservation']['pass'] is True
        assert rules_n2['4_delta_mj_conservation']['violations'] == 0

    def test_spatial_parity_passes(self, rules_n2):
        """Rule 5: Spatial parity."""
        assert rules_n2['5_spatial_parity']['pass'] is True
        assert rules_n2['5_spatial_parity']['violations'] == 0

    def test_furry_passes(self, rules_n2):
        """Rule 7 (Paper 33 Table 1): Furry's theorem passes (Dirac spinor
        phase kills the tadpole). Symbolic derivation of the mechanism:
        tests/test_paper33_furry_derivation.py."""
        assert rules_n2['7_furry_theorem']['pass'] is True
        assert rules_n2['7_furry_theorem']['tadpole_nonzero_count'] == 0

    def test_ward_identity_passes_at_n2(self, rules_n2):
        """Rule 6 (Paper 33 Table 1): Ward identity passes at n_max=2.

        At n_max=2, Sigma is block-diagonal within each n-shell, so
        [Sigma, H0] = 0 exactly (H0 = diag(n+1/2) is also n-diagonal).
        Breaks at n_max=3 where cross-n couplings appear.
        """
        assert rules_n2['6_ward_identity']['pass'] is True
        assert rules_n2['6_ward_identity']['ratio'] < 1e-12

    def test_triangle_on_n_passes(self, rules_n2):
        """Rule 8 (Paper 33 Table 1): triangle inequality on SO(4),
        |n_a - n_b| <= q <= n_a + n_b - 2 on the vertex tensor support.
        Replaced the tautological Hermiticity check 2026-07-03."""
        rule = rules_n2['8_triangle_on_n']
        assert rule['pass'] is True
        assert rule['violations'] == 0
        assert rule['nonzero_checked'] == 20  # known support at n_max=2

    def test_summary_total_is_eight(self, rules_n2):
        """Summary counts all 8 rules."""
        assert rules_n2['_summary']['total'] == 8


# ===========================================================================
# Scalar vs Dirac comparison
# ===========================================================================

class TestScalarDiracComparison:
    """Structural comparison: scalar+vector vs Dirac+vector."""

    def test_dirac_recovers_furry(self):
        """Dirac+vector gives 8/8 while scalar+vector gives 7/8.

        The Dirac spinor phase constraint (diagonal V=0) is a property of
        the Dirac alpha-matrix structure that has no scalar analog. Scalar
        electrons lack the off-diagonal large/small component structure,
        so their diagonal vertex coupling is nonzero.
        """
        scalar = compute_self_energy(n_max=2, q_max=1)
        dirac = compute_dirac_self_energy(n_max=2, q_max=1)
        s_sum = scalar.selection_rules['_summary']
        d_sum = dirac.selection_rules['_summary']
        assert s_sum['pass_count'] == 7
        assert d_sum['pass_count'] == 8
        assert s_sum['total'] == d_sum['total'] == 8

    def test_both_have_gs_zero(self):
        """Both configurations protect the ground state."""
        scalar = compute_self_energy(n_max=2, q_max=1)
        dirac = compute_dirac_self_energy(n_max=2, q_max=1)
        assert scalar.selection_rules['2_vertex_parity_gs_zero']['pass'] is True
        assert dirac.selection_rules['2_vertex_parity_gs_zero']['pass'] is True

    def test_furry_scalar_fails_dirac_passes(self):
        """Scalar fails Furry, Dirac passes via spinor phase constraint."""
        scalar = compute_self_energy(n_max=2, q_max=1)
        dirac = compute_dirac_self_energy(n_max=2, q_max=1)
        assert scalar.selection_rules['7_furry_theorem']['pass'] is False
        assert dirac.selection_rules['7_furry_theorem']['pass'] is True

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


# ===========================================================================
# Paper 33 census aggregates: the four-column partition table
# (Table tab:partition: scalar Fock 1/8, Dirac graph 4/8,
#  vector+scalar 7/8, vector+Dirac 8/8)
# ===========================================================================

def _scalar_fock_graph_tally(n_max: int = 2) -> dict:
    """Honest 8-rule tally for the GRAPH-NATIVE scalar-photon construction
    (Paper 28 / Paper 33 sec. 2: photon = 1-cochain on the Fock graph,
    G_gamma = L_1^+, CG-projected vertex).  Paper 33 Table 1 numbering.

    Added 2026-07-03: the 1/8 census leg previously had no aggregate test
    anywhere (its provenance was the archived driver
    debug/archive/qed_arc/gn_selection_rule_census.py). Each rule below is
    a real computation on production geovac modules; the per-rule criteria
    follow that archived census (the source of the paper's column).
    """
    import sympy as sp
    from sympy import Rational
    from geovac.graph_qed_vertex import (
        build_projection_matrix,
        build_vertex_tensor as gq_build_vertex,
        vertex_tensor_to_matrices,
    )
    from geovac.graph_qed_photon import build_fock_graph
    from geovac.graph_qed_propagator import (
        DiracGraphOperator,
        electron_propagator as gq_electron_propagator,
    )
    from geovac.graph_qed_self_energy import (
        compute_self_energy as gq_self_energy,
        compute_vertex_correction as gq_vertex_correction,
    )
    from geovac.qed_vertex import so4_channel_count

    P, dirac_labels, _ = build_projection_matrix(n_max)
    fock_data = build_fock_graph(n_max)
    entries, N_dirac, _, E_fock = gq_build_vertex(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data)
    tally = {}

    # Rule 1 (Gaunt/CG sparsity): every Fock edge is a ladder move
    # (L-edge dm=+-1 or T-edge dn=+-1) and the vertex tensor is sparse.
    ladder_ok = True
    for (v1, v2) in fock_data.edges:
        s1, s2 = fock_data.states[v1], fock_data.states[v2]
        dn, dl, dm = (abs(s1[0] - s2[0]), abs(s1[1] - s2[1]),
                      abs(s1[2] - s2[2]))
        if not ((dn == 0 and dl == 0 and dm == 1)
                or (dn == 1 and dl == 0 and dm == 0)):
            ladder_ok = False
    tally[1] = ladder_ok and (len(entries) < N_dirac * N_dirac * E_fock)

    # Rule 2 (GS structural zero): Sigma(GS block) == 0?
    # Pendant-edge theorem (Paper 28 prop): it is NOT (= 2(n_max-1)/n_max).
    se = gq_self_energy(n_max, t=Rational(0), exact=True)
    tally[2] = bool(se.ground_state_zero)

    # Rule 3 (SO(4) channel count): every realized coupling must carry a
    # continuum SO(4) channel W > 0 (CH convention, edge photon level
    # q_CH = min endpoint level - 1, per the archived census convention).
    w0 = 0
    for a_idx, b_idx, e_idx, _val in entries:
        da, db = dirac_labels[a_idx], dirac_labels[b_idx]
        v1, v2 = fock_data.edges[e_idx]
        q_ch = min(fock_data.states[v1][0], fock_data.states[v2][0]) - 1
        if so4_channel_count(da.n_fock - 1, db.n_fock - 1, q_ch) == 0:
            w0 += 1
    tally[3] = (w0 == 0)

    # Rule 4 (Delta m_j): each coupling's electron dm_j must be carried by
    # the photon edge's magnetic transfer (+- dm_edge).
    mj_viol = 0
    for a_idx, b_idx, e_idx, _val in entries:
        da, db = dirac_labels[a_idx], dirac_labels[b_idx]
        v1, v2 = fock_data.edges[e_idx]
        dm_edge = fock_data.states[v2][2] - fock_data.states[v1][2]
        d2mj = db.two_m_j - da.two_m_j
        if d2mj != 2 * dm_edge and d2mj != -2 * dm_edge:
            mj_viol += 1
    tally[4] = (mj_viol == 0)

    # Rule 5 (spatial parity E1): every realized coupling has l_a + l_b odd.
    par_viol = sum(
        1 for a_idx, b_idx, _e, _v in entries
        if (kappa_to_l(dirac_labels[a_idx].kappa)
            + kappa_to_l(dirac_labels[b_idx].kappa)) % 2 == 0)
    tally[5] = (par_viol == 0)

    # Rule 6 (Ward): commutator Ward identity [D, Lambda] == [Sigma, D].
    vc = gq_vertex_correction(n_max, t=Rational(0), exact=True)
    op = DiracGraphOperator(n_max=n_max, t=Rational(0))
    D_np = np.array(op.matrix_sympy().tolist(), dtype=float)
    S_np, L_np = se.Sigma_numpy, vc.Lambda_total_numpy
    ward_diff = np.max(np.abs((D_np @ L_np - L_np @ D_np)
                              - (S_np @ D_np - D_np @ S_np)))
    tally[6] = bool(ward_diff < 1e-10)

    # Rule 7 (Furry): tadpole sum_e Tr(V_e G_e) == 0 (sympy-exact).
    V_mats = vertex_tensor_to_matrices(entries, N_dirac, E_fock)
    G_e, _ = gq_electron_propagator(op, exact=True)
    tadpole = sum(sp.nsimplify((V_mats[e] * G_e).trace())
                  for e in range(E_fock))
    tally[7] = (sp.simplify(tadpole) == 0)

    # Rule 8 (triangle on n): q_edge >= 1 and |dn| <= q_edge <= n_a+n_b-2
    # with the edge photon level q_CH = min endpoint level - 1.
    tri_viol = 0
    for a_idx, b_idx, e_idx, _val in entries:
        da, db = dirac_labels[a_idx], dirac_labels[b_idx]
        v1, v2 = fock_data.edges[e_idx]
        q_ch = min(fock_data.states[v1][0], fock_data.states[v2][0]) - 1
        if not (q_ch >= 1 and abs(da.n_fock - db.n_fock) <= q_ch
                <= da.n_fock + db.n_fock - 2):
            tri_viol += 1
    tally[8] = (tri_viol == 0)

    return tally


def _native_dirac_graph_tally(n_max: int = 2) -> dict:
    """Honest 8-rule tally for the NATIVE DIRAC GRAPH construction
    (Paper 33 sec. 3: nodes (n, kappa, m_j), Rule B / E1 adjacency,
    scalar 1-cochain photon).  Paper 33 Table 1 numbering.

    Added 2026-07-03: the 4/8 census leg previously had no aggregate test
    anywhere (provenance: archived driver
    debug/archive/qed_arc/dirac_graph_qed_sprint.py, whose per-rule
    criteria this tally follows, on the production Rule B graph builder).
    """
    from geovac.ihara_zeta_dirac import build_dirac_s3_graph

    A, labels, _deg, _desc = build_dirac_s3_graph(n_max, "B")
    Vn = len(labels)
    edges = sorted((i, j) for i in range(Vn) for j in range(i + 1, Vn)
                   if A[i, j] != 0)
    E = len(edges)
    tally = {}

    # Photon propagator G_gamma = L_1^+ and symmetrized edge vertices.
    B_inc = np.zeros((Vn, E))
    for k, (i, j) in enumerate(edges):
        B_inc[i, k] = 1.0
        B_inc[j, k] = -1.0
    L1 = B_inc.T @ B_inc
    G_gamma = np.linalg.pinv(L1)
    V_mats = []
    for (i, j) in edges:
        Ve = np.zeros((Vn, Vn))
        Ve[i, j] = 1.0
        Ve[j, i] = 1.0
        V_mats.append(Ve)
    Sigma = np.zeros((Vn, Vn))
    for e1 in range(E):
        for e2 in range(E):
            g = G_gamma[e1, e2]
            if abs(g) > 1e-15:
                Sigma += g * (V_mats[e1] @ V_mats[e2].T)

    # Rule 1 (Gaunt/CG): every Rule B edge is E1/CG-allowed
    # (|dl| = 1, |dm_j| <= 1) and the edge set is sparse.
    ok1 = all(
        abs(kappa_to_l(labels[i].kappa) - kappa_to_l(labels[j].kappa)) == 1
        and abs(labels[i].two_m_j - labels[j].two_m_j) <= 2
        for (i, j) in edges)
    tally[1] = ok1 and E < Vn * (Vn - 1) // 2

    # Rule 2 (GS zero): Sigma(GS block) == 0? (GS couples through its
    # radial edges, so it is NOT.)
    gs_idx = [i for i, lab in enumerate(labels) if lab.n_fock == 1]
    tally[2] = bool(np.allclose(Sigma[np.ix_(gs_idx, gs_idx)], 0,
                                atol=1e-12))

    # Rule 3 (SO(4) channel count): every edge coupling must carry a
    # continuum W > 0 (CH convention, edge level q_CH = min n_fock - 1).
    from geovac.qed_vertex import so4_channel_count
    w0 = sum(
        1 for (i, j) in edges
        if so4_channel_count(
            labels[i].n_fock - 1, labels[j].n_fock - 1,
            min(labels[i].n_fock, labels[j].n_fock) - 1) == 0)
    tally[3] = (w0 == 0)

    # Rule 4 (Delta m_j): every photon edge carries |dm_j| <= 1
    # (a spin-1 photon can balance it) -- built into the node labels.
    tally[4] = all(abs(labels[i].two_m_j - labels[j].two_m_j) <= 2
                   for (i, j) in edges)

    # Rule 5 (spatial parity E1): every edge has l_i + l_j odd.
    tally[5] = all(
        (kappa_to_l(labels[i].kappa) + kappa_to_l(labels[j].kappa)) % 2 == 1
        for (i, j) in edges)

    # Rule 6 (Ward): [D, Lambda] == [Sigma, D] with the CH Dirac operator.
    D_op = np.zeros((Vn, Vn))
    for i, lab in enumerate(labels):
        chi = 1 if lab.kappa < 0 else -1
        D_op[i, i] = chi * ((lab.n_fock - 1) + 1.5)
    G_e = np.linalg.inv(D_op)
    Lam = np.zeros((Vn, Vn))
    for e1 in range(E):
        for e2 in range(E):
            g = G_gamma[e1, e2]
            if abs(g) > 1e-15:
                Lam += g * (V_mats[e1] @ G_e @ V_mats[e2].T)
    ward_diff = np.max(np.abs((D_op @ Lam - Lam @ D_op)
                              - (Sigma @ D_op - D_op @ Sigma)))
    tally[6] = bool(ward_diff < 1e-10)

    # Rule 7 (Furry): tadpole AND 3-vertex closed loop both vanish
    # (the off-diagonal edge vertex + l-parity bipartite structure).
    tadpole = sum(G_e[i, j] + G_e[j, i] for (i, j) in edges)
    triangle = 0.0
    for e1 in range(E):
        VG1 = V_mats[e1] @ G_e
        for e2 in range(E):
            VG12 = VG1 @ (V_mats[e2] @ G_e)
            for e3 in range(E):
                triangle += np.trace(VG12 @ V_mats[e3] @ G_e)
    tally[7] = bool(abs(tadpole) < 1e-12 and abs(triangle) < 1e-12)

    # Rule 8 (triangle on n): q_edge >= 1 and |dn| <= q_edge <= n_i+n_j-2.
    tally[8] = all(
        (min(labels[i].n_fock, labels[j].n_fock) - 1) >= 1
        and abs(labels[i].n_fock - labels[j].n_fock)
        <= (min(labels[i].n_fock, labels[j].n_fock) - 1)
        <= labels[i].n_fock + labels[j].n_fock - 2
        for (i, j) in edges)

    return tally


class TestSelectionRulesN3Regression:
    """Paper 33 R3/R8 at n_max = 3: the de-tautologized census checks
    (SO(4) channel count vs closed form; triangle-on-n on the realized
    support) hold at the larger truncation the paper claims ('verified
    at n_max = 2, 3').  Added 2026-07-04 (cert-run coverage-gap
    closure: only n_max = 2 was regression-tested for these two
    rule keys)."""

    def test_scalar_rules_3_and_8_pass_at_n3(self):
        rules = compute_self_energy(3, 1).selection_rules
        assert rules['3_so4_channel_count']['pass'] is True
        assert rules['8_triangle_on_n']['pass'] is True

    def test_dirac_rules_3_and_8_pass_at_n3(self):
        rules = compute_dirac_self_energy(3, 1).selection_rules
        assert rules['3_so4_channel_count']['pass'] is True
        assert rules['8_triangle_on_n']['pass'] is True


class TestCensusAggregates:
    """Aggregate census tests for the two graph-photon legs of Paper 33's
    partition table (the vector-photon legs 7/8 and 8/8 are asserted in
    TestScalarDiracComparison / TestDiracSelectionRules)."""

    def test_scalar_fock_census_is_1_of_8(self):
        """Paper 33 sec. 2 / Table tab:partition column 1: the graph-native
        scalar-photon baseline recovers EXACTLY 1/8 rules (Gaunt/CG only)."""
        tally = _scalar_fock_graph_tally(2)
        expected = {1: True, 2: False, 3: False, 4: False,
                    5: False, 6: False, 7: False, 8: False}
        assert tally == expected, f"scalar-Fock tally moved: {tally}"
        assert sum(tally.values()) == 1

    def test_native_dirac_graph_census_is_4_of_8(self):
        """Paper 33 sec. 3 / Table tab:partition column 2: the native Dirac
        graph (Rule B) recovers EXACTLY 4/8 rules (Gaunt/CG, Delta m_j,
        spatial parity, Furry)."""
        tally = _native_dirac_graph_tally(2)
        expected = {1: True, 2: False, 3: False, 4: True,
                    5: True, 6: False, 7: True, 8: False}
        assert tally == expected, f"native-Dirac tally moved: {tally}"
        assert sum(tally.values()) == 4

    def test_partition_table_totals(self):
        """The full bottom row of Table tab:partition: 1/8, 4/8, 7/8, 8/8."""
        totals = (
            sum(_scalar_fock_graph_tally(2).values()),
            sum(_native_dirac_graph_tally(2).values()),
            compute_self_energy(2, 1).selection_rules['_summary']['pass_count'],
            compute_dirac_self_energy(2, 1).selection_rules['_summary']['pass_count'],
        )
        assert totals == (1, 4, 7, 8)


def test_scalar_l1_self_energy_is_one_over_4pi():
    """Paper 33 Prop [1/(4pi)]: each l=1 diagonal self-energy entry of the
    scalar+vector construction at n_max=2, q_max=1 equals exactly 1/(4pi)
    (the S^2 Weyl exchange constant). Added 1st-cert remediation 2026-07-02:
    the proposition previously cited a nonexistent test name."""
    import math

    result = compute_self_energy(n_max=2, q_max=1)
    diag = [float(result.Sigma[i, i].real if hasattr(result.Sigma[i, i], "real")
                  else result.Sigma[i, i]) for i in range(result.Sigma.shape[0])]
    target = 1.0 / (4.0 * math.pi)
    # indices 0,1 = 1s,2s (structural zeros); 2,3,4 = the 2p triplet
    assert diag[0] == 0.0 and diag[1] == 0.0
    for i in (2, 3, 4):
        assert abs(diag[i] - target) < 1e-14, (i, diag[i], target)
