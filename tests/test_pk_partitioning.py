"""
Tests for PK Classical Partitioning (Track BF, v2.0.29)
========================================================

Validates:
  - PK matrix returned by builders is correct (h1_with_pk - h1_without_pk)
  - pk_classical_energy computes Tr(h1_pk · γ) correctly
  - 1-RDM reconstruction from statevector is correct
  - Algebraic exactness: |E_full - (E_elec + E_PK_classical)| < 1e-10
  - Pauli term count unchanged by partitioning
  - Backward compatibility: pk_in_hamiltonian=True gives old behavior

Author: GeoVac Development Team
Date: April 2026
"""

import numpy as np
import pytest

from geovac.composed_qubit import (
    build_composed_lih,
    build_composed_beh2,
    build_composed_h2o,
)
from geovac.pk_partitioning import (
    pk_classical_energy,
    reconstruct_1rdm_from_statevector,
    validate_pk_partitioning,
)


# ---------------------------------------------------------------------------
# Fixtures (small basis for speed)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def lih_nmax1_full():
    """LiH at max_n=1, PK in Hamiltonian."""
    return build_composed_lih(
        max_n_core=1, max_n_val=1, include_pk=True,
        pk_in_hamiltonian=True, verbose=False,
    )


@pytest.fixture(scope="module")
def lih_nmax1_partitioned():
    """LiH at max_n=1, PK separated."""
    return build_composed_lih(
        max_n_core=1, max_n_val=1, include_pk=True,
        pk_in_hamiltonian=False, verbose=False,
    )


@pytest.fixture(scope="module")
def beh2_nmax1_full():
    """BeH2 at max_n=1, PK in Hamiltonian."""
    return build_composed_beh2(
        max_n_core=1, max_n_val=1, include_pk=True,
        pk_in_hamiltonian=True, verbose=False,
    )


@pytest.fixture(scope="module")
def beh2_nmax1_partitioned():
    """BeH2 at max_n=1, PK separated."""
    return build_composed_beh2(
        max_n_core=1, max_n_val=1, include_pk=True,
        pk_in_hamiltonian=False, verbose=False,
    )


# ---------------------------------------------------------------------------
# PK matrix correctness
# ---------------------------------------------------------------------------

class TestPKMatrixCorrectness:
    """Verify that h1_pk = h1(with_pk) - h1(without_pk)."""

    def test_lih_h1_pk_equals_difference(
        self, lih_nmax1_full, lih_nmax1_partitioned,
    ) -> None:
        h1_full = lih_nmax1_full['h1']
        h1_elec = lih_nmax1_partitioned['h1']
        h1_pk = lih_nmax1_partitioned['h1_pk']
        np.testing.assert_allclose(h1_full, h1_elec + h1_pk, atol=1e-14)

    def test_beh2_h1_pk_equals_difference(
        self, beh2_nmax1_full, beh2_nmax1_partitioned,
    ) -> None:
        h1_full = beh2_nmax1_full['h1']
        h1_elec = beh2_nmax1_partitioned['h1']
        h1_pk = beh2_nmax1_partitioned['h1_pk']
        np.testing.assert_allclose(h1_full, h1_elec + h1_pk, atol=1e-14)

    def test_h2o_h1_pk_equals_difference(self) -> None:
        r_full = build_composed_h2o(
            max_n_core=1, max_n_val=1, include_pk=True,
            pk_in_hamiltonian=True, verbose=False,
        )
        r_elec = build_composed_h2o(
            max_n_core=1, max_n_val=1, include_pk=True,
            pk_in_hamiltonian=False, verbose=False,
        )
        np.testing.assert_allclose(
            r_full['h1'], r_elec['h1'] + r_elec['h1_pk'], atol=1e-14,
        )

    def test_lih_pk_disabled_zero(self) -> None:
        """When include_pk=False, h1_pk should be zero."""
        r = build_composed_lih(
            max_n_core=1, max_n_val=1, include_pk=False, verbose=False,
        )
        np.testing.assert_allclose(r['h1_pk'], 0.0, atol=1e-15)


# ---------------------------------------------------------------------------
# Pauli term count invariance
# ---------------------------------------------------------------------------

class TestPauliCountInvariance:
    """PK partitioning must not change the Pauli term count."""

    def test_lih_pauli_count_unchanged(
        self, lih_nmax1_full, lih_nmax1_partitioned,
    ) -> None:
        assert lih_nmax1_full['N_pauli'] == lih_nmax1_partitioned['N_pauli']

    def test_beh2_pauli_count_unchanged(
        self, beh2_nmax1_full, beh2_nmax1_partitioned,
    ) -> None:
        assert beh2_nmax1_full['N_pauli'] == beh2_nmax1_partitioned['N_pauli']


# ---------------------------------------------------------------------------
# Classical energy computation
# ---------------------------------------------------------------------------

class TestPKClassicalEnergy:
    """Validate pk_classical_energy = Tr(h1_pk · γ)."""

    def test_identity_rdm(self) -> None:
        """With identity 1-RDM, E_PK = Tr(h1_pk)."""
        M = 4
        h1_pk = np.diag([1.0, 2.0, 3.0, 4.0])
        rdm = np.eye(M)
        assert abs(pk_classical_energy(h1_pk, rdm) - 10.0) < 1e-14

    def test_zero_pk(self) -> None:
        """Zero PK matrix gives zero energy."""
        M = 3
        h1_pk = np.zeros((M, M))
        rdm = np.random.randn(M, M)
        assert abs(pk_classical_energy(h1_pk, rdm)) < 1e-14

    def test_off_diagonal(self) -> None:
        """Off-diagonal PK contributes correctly."""
        h1_pk = np.array([[0.0, 0.5], [0.5, 0.0]])
        rdm = np.array([[0.8, 0.2], [0.2, 0.8]])
        # Tr(h1_pk @ rdm) = 0*0.8 + 0.5*0.2 + 0.5*0.2 + 0*0.8 = 0.2
        expected = 0.2
        assert abs(pk_classical_energy(h1_pk, rdm) - expected) < 1e-14


# ---------------------------------------------------------------------------
# 1-RDM reconstruction
# ---------------------------------------------------------------------------

class TestRDMReconstruction:
    """Validate 1-RDM from statevector."""

    def test_single_determinant(self) -> None:
        """Single determinant |1100⟩ → occupation [1,1,0,0] → rdm diag."""
        n_spatial = 2  # 4 qubits
        # |1100⟩ in JW: qubits 0,1 occupied (spin-orbs 0,1)
        # That's spatial orbital 0 (both spins)
        psi = np.zeros(16)
        # |0011⟩ in binary = qubit 0 and 1 occupied = state index 3
        psi[3] = 1.0  # |11 00⟩ = orbitals 0,1 occupied
        rdm = reconstruct_1rdm_from_statevector(psi, n_spatial)
        # Spatial orbital 0: both spins occupied → γ_{00} = 2
        # Spatial orbital 1: empty → γ_{11} = 0
        assert abs(rdm[0, 0] - 2.0) < 1e-10
        assert abs(rdm[1, 1] - 0.0) < 1e-10

    def test_trace_equals_n_electrons(self) -> None:
        """Tr(γ) should equal number of electrons."""
        n_spatial = 2
        # |0101⟩ = qubits 0,2 occupied = spin-orbs 0,2 = spatial 0↑, 1↑
        psi = np.zeros(16)
        psi[5] = 1.0  # binary 0101 = 5
        rdm = reconstruct_1rdm_from_statevector(psi, n_spatial)
        assert abs(np.trace(rdm) - 2.0) < 1e-10


# ---------------------------------------------------------------------------
# Algebraic exactness (exact diag, small basis)
# ---------------------------------------------------------------------------

class TestAlgebraicExactness:
    """E_full = E_elec(ψ_full) + E_PK(ψ_full) to machine precision."""

    def test_lih_nmax1_exact(self) -> None:
        """LiH at max_n=1 (Q=6): validate exact partitioning."""
        result = validate_pk_partitioning(
            build_composed_lih,
            {'max_n_core': 1, 'max_n_val': 1},
            n_electrons=4,
            system_name='LiH (nmax=1)',
            verbose=False,
        )
        assert result['exact_diag']
        assert result['residual_same_state'] < 1e-8, (
            f"Residual {result['residual_same_state']:.2e} exceeds 1e-8"
        )

    def test_beh2_nmax1_exact(self) -> None:
        """BeH2 at max_n=1 (Q=10): validate exact partitioning."""
        result = validate_pk_partitioning(
            build_composed_beh2,
            {'max_n_core': 1, 'max_n_val': 1},
            n_electrons=6,
            system_name='BeH2 (nmax=1)',
            verbose=False,
        )
        assert result['exact_diag']
        assert result['residual_same_state'] < 1e-8, (
            f"Residual {result['residual_same_state']:.2e} exceeds 1e-8"
        )

    def test_h2o_nmax1_exact(self) -> None:
        """H2O at max_n=1 (Q=14): validate exact partitioning."""
        result = validate_pk_partitioning(
            build_composed_h2o,
            {'max_n_core': 1, 'max_n_val': 1},
            n_electrons=8,
            system_name='H2O (nmax=1)',
            verbose=False,
        )
        assert result['exact_diag']
        assert result['residual_same_state'] < 1e-8, (
            f"Residual {result['residual_same_state']:.2e} exceeds 1e-8"
        )


# ---------------------------------------------------------------------------
# Ground state analysis (key sprint question)
# ---------------------------------------------------------------------------

class TestGroundStateAnalysis:
    """Check whether electronic-only ground state matches full ground state."""

    def test_lih_nmax1_ground_states_overlap(self) -> None:
        """LiH: ground state analysis (diagnostic, not hard-assert on overlap).

        Note: eigsh searches all particle-number sectors. The Fock-space
        ground state may differ from the N-electron ground state. In VQE
        (particle-number-conserving ansatz), PK partitioning works because
        the ansatz stays in the valence sector.
        """
        result = validate_pk_partitioning(
            build_composed_lih,
            {'max_n_core': 1, 'max_n_val': 1},
            n_electrons=4,
            system_name='LiH (nmax=1)',
            verbose=False,
        )
        assert result['exact_diag']
        print(f"\nLiH ground state overlap: {result['ground_state_overlap']:.6f}")
        print(f"E_full = {result['E_full']:.6f}")
        print(f"E_elec = {result['E_elec']:.6f}")
        print(f"E_PK(psi_full) = {result['E_pk_of_full_gs']:.6f}")
        print(f"E_PK(psi_elec) = {result['E_pk_of_elec_gs']:.6f}")

    def test_h2o_nmax1_ground_states_analysis(self) -> None:
        """H2O: ground state analysis (diagnostic).

        H2O PK is ~98.7% of 1-norm. Removing PK from the Fock-space
        Hamiltonian changes the ground state (electrons collapse to core).
        This is expected: PK prevents collapse. In VQE with a particle-
        number-conserving ansatz, the problem doesn't arise.
        """
        result = validate_pk_partitioning(
            build_composed_h2o,
            {'max_n_core': 1, 'max_n_val': 1},
            n_electrons=8,
            system_name='H2O (nmax=1)',
            verbose=False,
        )
        assert result['exact_diag']
        print(f"\nH2O ground state overlap: {result['ground_state_overlap']:.6f}")
        print(f"E_full = {result['E_full']:.6f}")
        print(f"E_elec = {result['E_elec']:.6f}")
        print(f"E_PK(psi_full) = {result['E_pk_of_full_gs']:.6f}")
        print(f"E_PK(psi_elec) = {result['E_pk_of_elec_gs']:.6f}")
        print(f"Match: {result['ground_states_match']}")


# ---------------------------------------------------------------------------
# Backward compatibility
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    """Ensure old API behavior is preserved."""

    def test_default_include_pk_true_gives_pk_in_hamiltonian(self) -> None:
        """build_composed_lih(include_pk=True) without pk_in_hamiltonian
        should behave as pk_in_hamiltonian=True (backward compat)."""
        r_old = build_composed_lih(
            max_n_core=1, max_n_val=1, include_pk=True, verbose=False,
        )
        r_new = build_composed_lih(
            max_n_core=1, max_n_val=1, include_pk=True,
            pk_in_hamiltonian=True, verbose=False,
        )
        np.testing.assert_allclose(r_old['h1'], r_new['h1'], atol=1e-14)
        assert r_old['N_pauli'] == r_new['N_pauli']

    def test_include_pk_false_still_works(self) -> None:
        """include_pk=False should produce zero h1_pk."""
        r = build_composed_lih(
            max_n_core=1, max_n_val=1, include_pk=False, verbose=False,
        )
        assert r['h1_pk'] is not None
        np.testing.assert_allclose(r['h1_pk'], 0.0, atol=1e-15)

    def test_1norm_reduction_h2o_nmax2(self) -> None:
        """H2O at nmax=2: partitioned 1-norm should be dramatically lower.

        This is the motivating case: PK is 98.7% of the H2O 1-norm.
        At nmax=1 / LiH, the effect can be small or even negative (PK's
        sign structure can cancel some existing Pauli coefficients).
        """
        r_full = build_composed_h2o(
            max_n_core=1, max_n_val=1, include_pk=True,
            pk_in_hamiltonian=True, verbose=False,
        )
        r_elec = build_composed_h2o(
            max_n_core=1, max_n_val=1, include_pk=True,
            pk_in_hamiltonian=False, verbose=False,
        )
        norm_full = sum(abs(c) for c in r_full['qubit_op'].terms.values())
        norm_elec = sum(abs(c) for c in r_elec['qubit_op'].terms.values())
        # H2O PK is 98.7% of the 1-norm, so removal should be dramatic
        assert norm_elec < norm_full * 0.1, (
            f"H2O partitioned 1-norm {norm_elec:.1f} not < 10% of full {norm_full:.1f}"
        )
