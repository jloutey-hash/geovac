"""
Tests for geovac.qubit_encoding — Jordan-Wigner qubit Hamiltonian encoding.

Validates:
  1. FermionOperator construction from LatticeIndex integrals
  2. Pauli term counts match known values
  3. ERI sparsity decreases with basis size (selection rules)
  4. Max Pauli weight is bounded
  5. Molecular (MolecularLatticeIndex) encoding works
  6. Energy consistency: JW Hamiltonian eigenvalue matches FCI solver

Author: GeoVac Development Team
Date: March 2026
"""

import warnings
import pytest
import numpy as np

from openfermion import jordan_wigner, get_sparse_operator

from geovac.lattice_index import LatticeIndex, MolecularLatticeIndex
from geovac.qubit_encoding import (
    JordanWignerEncoder,
    PauliAnalysis,
    build_fermion_op_from_integrals,
    fit_pauli_scaling,
)

warnings.filterwarnings('ignore', category=UserWarning)


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def he_nmax2() -> LatticeIndex:
    """He atom (Z=2, 2 electrons) at nmax=2."""
    return LatticeIndex(
        n_electrons=2, max_n=2, nuclear_charge=2,
        vee_method='slater_full', h1_method='hybrid', fci_method='matrix',
    )


@pytest.fixture(scope="module")
def he_nmax3() -> LatticeIndex:
    """He atom (Z=2, 2 electrons) at nmax=3."""
    return LatticeIndex(
        n_electrons=2, max_n=3, nuclear_charge=2,
        vee_method='slater_full', h1_method='hybrid', fci_method='matrix',
    )


@pytest.fixture(scope="module")
def h_nmax2() -> LatticeIndex:
    """H atom (Z=1, 1 electron) at nmax=2."""
    return LatticeIndex(
        n_electrons=1, max_n=2, nuclear_charge=1,
        vee_method='slater_full', h1_method='hybrid', fci_method='matrix',
    )


# --------------------------------------------------------------------------
# Basic construction tests
# --------------------------------------------------------------------------

class TestJordanWignerEncoder:
    """Tests for JordanWignerEncoder class."""

    def test_build_fermion_op_he(self, he_nmax2: LatticeIndex) -> None:
        """FermionOperator is constructed without error."""
        enc = JordanWignerEncoder(he_nmax2)
        fop = enc.build_fermion_operator()
        assert len(fop.terms) > 0

    def test_build_qubit_op_he(self, he_nmax2: LatticeIndex) -> None:
        """QubitOperator is constructed via JW transform."""
        enc = JordanWignerEncoder(he_nmax2)
        qop = enc.build_qubit_operator()
        assert len(qop.terms) > 0

    def test_qubit_count_matches_spinorbs(self, he_nmax2: LatticeIndex) -> None:
        """Number of qubits equals number of spin-orbitals."""
        enc = JordanWignerEncoder(he_nmax2)
        analysis = enc.analyze()
        assert analysis.n_qubits == he_nmax2.n_sp

    def test_pauli_terms_positive(self, he_nmax2: LatticeIndex) -> None:
        """Pauli term count is positive."""
        enc = JordanWignerEncoder(he_nmax2)
        assert enc.count_pauli_terms() > 0

    def test_max_weight_bounded(self, he_nmax2: LatticeIndex) -> None:
        """Max Pauli weight <= 2 * n_qubits (JW bound for 2-body)."""
        enc = JordanWignerEncoder(he_nmax2)
        analysis = enc.analyze()
        assert analysis.max_pauli_weight <= 2 * analysis.n_qubits

    def test_caching(self, he_nmax2: LatticeIndex) -> None:
        """Repeated calls return cached results."""
        enc = JordanWignerEncoder(he_nmax2)
        a1 = enc.analyze()
        a2 = enc.analyze()
        assert a1 is a2

    def test_single_electron_h(self, h_nmax2: LatticeIndex) -> None:
        """Single-electron H atom produces valid encoding."""
        enc = JordanWignerEncoder(h_nmax2)
        analysis = enc.analyze()
        assert analysis.n_pauli_terms > 0
        # LatticeIndex precomputes ERIs even for 1e (used if n_electrons>1),
        # so we just check that encoding completes successfully
        assert analysis.n_qubits == h_nmax2.n_sp


# --------------------------------------------------------------------------
# ERI sparsity tests — validates selection rule structure
# --------------------------------------------------------------------------

class TestERISparsity:
    """Validates that ERI density decreases with basis size."""

    def test_eri_density_decreases(self, he_nmax2: LatticeIndex, he_nmax3: LatticeIndex) -> None:
        """ERI density should decrease from nmax=2 to nmax=3."""
        enc2 = JordanWignerEncoder(he_nmax2)
        enc3 = JordanWignerEncoder(he_nmax3)
        a2 = enc2.analyze()
        a3 = enc3.analyze()
        assert a3.eri_density < a2.eri_density, (
            f"ERI density should decrease: nmax=2 {a2.eri_density:.4f} "
            f"vs nmax=3 {a3.eri_density:.4f}"
        )

    def test_eri_density_below_50pct(self, he_nmax3: LatticeIndex) -> None:
        """GeoVac ERI density at nmax=3 should be well below Gaussian ~50-100%."""
        enc = JordanWignerEncoder(he_nmax3)
        analysis = enc.analyze()
        assert analysis.eri_density < 0.10, (
            f"ERI density {analysis.eri_density:.4f} should be < 10% at nmax=3"
        )


# --------------------------------------------------------------------------
# Energy consistency — JW eigenvalue vs FCI solver
# --------------------------------------------------------------------------

class TestEnergyConsistency:
    """JW Hamiltonian ground state matches FCI solver."""

    def test_he_nmax2_energy(self, he_nmax2: LatticeIndex) -> None:
        """
        He at nmax=2: JW ground state eigenvalue matches FCI.

        This validates that the FermionOperator encodes the same
        Hamiltonian as LatticeIndex.assemble_hamiltonian().
        """
        # FCI energy
        eigvals, _ = he_nmax2.compute_ground_state(n_states=1)
        e_fci = eigvals[0]

        # JW energy via sparse matrix diagonalization
        enc = JordanWignerEncoder(he_nmax2)
        fop = enc.build_fermion_operator()
        qop = jordan_wigner(fop)
        sparse_h = get_sparse_operator(qop, n_qubits=he_nmax2.n_sp)
        from scipy.sparse.linalg import eigsh
        e_jw = eigsh(sparse_h, k=1, which='SA', return_eigenvectors=False)[0]

        assert abs(e_jw - e_fci) < 1e-4, (
            f"JW energy {e_jw:.8f} != FCI energy {e_fci:.8f}"
        )


# --------------------------------------------------------------------------
# Molecular encoding test
# --------------------------------------------------------------------------

class TestMolecularEncoding:
    """Tests for MolecularLatticeIndex encoding."""

    @pytest.fixture(scope="class")
    def h2_mol(self) -> MolecularLatticeIndex:
        """H2 molecule at R=1.4 bohr, nmax=2."""
        return MolecularLatticeIndex(
            Z_A=1, Z_B=1, nmax_A=2, nmax_B=2,
            R=1.4, n_electrons=2, n_bridges=10,
            vee_method='slater_full',
        )

    def test_h2_encoding(self, h2_mol: MolecularLatticeIndex) -> None:
        """H2 molecule produces valid qubit encoding."""
        enc = JordanWignerEncoder(h2_mol)
        analysis = enc.analyze()
        assert analysis.n_pauli_terms > 0
        assert analysis.n_qubits == h2_mol.n_sp

    def test_h2_nuclear_repulsion(self, h2_mol: MolecularLatticeIndex) -> None:
        """Nuclear repulsion is automatically extracted."""
        enc = JordanWignerEncoder(h2_mol)
        assert enc._v_nn == pytest.approx(1.0 / 1.4, abs=1e-10)


# --------------------------------------------------------------------------
# Utility function tests
# --------------------------------------------------------------------------

class TestUtilities:
    """Tests for utility functions."""

    def test_build_from_integrals_sto3g(self) -> None:
        """build_fermion_op_from_integrals works for STO-3G H2."""
        h1 = np.array([[-1.2528, 0.0], [0.0, -0.4756]])
        eri = np.zeros((2, 2, 2, 2))
        eri[0, 0, 0, 0] = 0.6746
        eri[1, 1, 1, 1] = 0.6975
        eri[0, 0, 1, 1] = 0.6632
        eri[1, 1, 0, 0] = 0.6632
        eri[0, 1, 0, 1] = 0.1813
        eri[1, 0, 1, 0] = 0.1813
        eri[0, 1, 1, 0] = 0.1813
        eri[1, 0, 0, 1] = 0.1813

        fop = build_fermion_op_from_integrals(h1, eri, 1.0 / 1.4)
        qop = jordan_wigner(fop)
        n_pauli = len(qop.terms)
        assert n_pauli == 15, f"STO-3G H2 should have 15 Pauli terms, got {n_pauli}"

    def test_fit_scaling(self) -> None:
        """fit_pauli_scaling recovers known exponent."""
        q = np.array([10, 28, 60, 110])
        # Synthetic data with exponent ~3.0
        p = 2.0 * q ** 3.0
        exp, _ = fit_pauli_scaling(q, p)
        assert abs(exp - 3.0) < 0.01, f"Expected exponent ~3.0, got {exp:.2f}"


# --------------------------------------------------------------------------
# PauliAnalysis dataclass tests
# --------------------------------------------------------------------------

class TestPauliAnalysis:
    """Tests for PauliAnalysis summary."""

    def test_summary_string(self, he_nmax2: LatticeIndex) -> None:
        """Summary string contains key metrics."""
        enc = JordanWignerEncoder(he_nmax2)
        analysis = enc.analyze()
        s = analysis.summary()
        assert "Q=" in s
        assert "Pauli=" in s
        assert "max_weight=" in s
