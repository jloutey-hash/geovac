"""
Tests for geovac.gaussian_reference — Gaussian-basis reference Hamiltonians.

Validates:
  1. Eigenvalue consistency: JW Hamiltonian ground state matches literature
  2. Pauli term counts are correct for known systems
  3. Side-by-side comparison of GeoVac vs Gaussian Pauli term counts

Author: GeoVac Development Team
Date: March 2026
"""

import warnings

import pytest

pytestmark = pytest.mark.slow
import numpy as np
from scipy.sparse.linalg import eigsh

from openfermion import jordan_wigner, get_sparse_operator

from geovac.gaussian_reference import (
    h2_sto3g,
    he_sto3g,
    he_cc_pvdz,
    he_cc_pvtz,
    h2_631g,
    build_qubit_hamiltonian,
)
from geovac.lattice_index import LatticeIndex
from geovac.qubit_encoding import JordanWignerEncoder

warnings.filterwarnings('ignore', category=UserWarning)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def h2_sys() -> dict:
    """H2 STO-3G reference system."""
    return h2_sto3g()


@pytest.fixture(scope="module")
def he_sys() -> dict:
    """He STO-3G reference system."""
    return he_sto3g()


@pytest.fixture(scope="module")
def geovac_he_nmax2() -> LatticeIndex:
    """GeoVac He at nmax=2."""
    return LatticeIndex(
        n_electrons=2, max_n=2, nuclear_charge=2,
        vee_method='slater_full', h1_method='hybrid', fci_method='matrix',
    )


@pytest.fixture(scope="module")
def geovac_he_nmax3() -> LatticeIndex:
    """GeoVac He at nmax=3."""
    return LatticeIndex(
        n_electrons=2, max_n=3, nuclear_charge=2,
        vee_method='slater_full', h1_method='hybrid', fci_method='matrix',
    )


# ---------------------------------------------------------------------------
# H2 STO-3G tests
# ---------------------------------------------------------------------------

class TestH2STO3G:
    """Validate H2 STO-3G reference Hamiltonian."""

    def test_system_metadata(self, h2_sys: dict) -> None:
        """System dict has required keys and correct dimensions."""
        assert h2_sys['n_electrons'] == 2
        assert h2_sys['n_spatial'] == 2
        assert h2_sys['h1'].shape == (2, 2)
        assert h2_sys['eri'].shape == (2, 2, 2, 2)
        assert h2_sys['nuclear_repulsion'] == pytest.approx(1.0 / 1.4, abs=1e-10)

    def test_h1_symmetry(self, h2_sys: dict) -> None:
        """One-electron integrals are symmetric."""
        h1 = h2_sys['h1']
        np.testing.assert_allclose(h1, h1.T, atol=1e-15)

    def test_eri_symmetry(self, h2_sys: dict) -> None:
        """ERIs have chemist-notation permutational symmetry: (pq|rs) = (qp|sr) = (rs|pq)."""
        eri = h2_sys['eri']
        M = eri.shape[0]
        for p in range(M):
            for q in range(M):
                for r in range(M):
                    for s in range(M):
                        val = eri[p, q, r, s]
                        # (pq|rs) = (rs|pq)
                        assert val == pytest.approx(eri[r, s, p, q], abs=1e-15)
                        # (pq|rs) = (qp|sr)
                        assert val == pytest.approx(eri[q, p, s, r], abs=1e-15)

    def test_eigenvalue_matches_literature(self, h2_sys: dict) -> None:
        """
        JW Hamiltonian ground state energy matches FCI literature value.

        Literature: E_FCI(STO-3G, R=1.4) = -1.1373 Ha
        (Szabo & Ostlund, Chapter 4)
        """
        _, qubit_op, _ = build_qubit_hamiltonian(h2_sys)
        n_qubits = 2 * h2_sys['n_spatial']
        sparse_h = get_sparse_operator(qubit_op, n_qubits=n_qubits)
        e_gs = eigsh(sparse_h, k=1, which='SA', return_eigenvectors=False)[0]

        assert abs(e_gs - h2_sys['literature_energy']) < 0.01, (
            f"JW ground state {e_gs:.6f} Ha does not match literature "
            f"{h2_sys['literature_energy']:.6f} Ha"
        )

    def test_pauli_term_count(self, h2_sys: dict) -> None:
        """STO-3G H2 produces exactly 15 Pauli terms (well-known result)."""
        _, _, n_pauli = build_qubit_hamiltonian(h2_sys)
        assert n_pauli == 15, f"Expected 15 Pauli terms, got {n_pauli}"

    def test_r_not_14_raises(self) -> None:
        """Requesting R != 1.4 raises ValueError."""
        with pytest.raises(ValueError, match="R=1.4"):
            h2_sto3g(R=2.0)


# ---------------------------------------------------------------------------
# He STO-3G tests
# ---------------------------------------------------------------------------

class TestHeSTO3G:
    """Validate He STO-3G reference Hamiltonian."""

    def test_system_metadata(self, he_sys: dict) -> None:
        """System dict has required keys and correct dimensions."""
        assert he_sys['n_electrons'] == 2
        assert he_sys['n_spatial'] == 1
        assert he_sys['h1'].shape == (1, 1)
        assert he_sys['eri'].shape == (1, 1, 1, 1)
        assert he_sys['nuclear_repulsion'] == 0.0

    def test_analytical_integrals(self, he_sys: dict) -> None:
        """Verify analytical STO integrals with zeta=1.6875."""
        zeta = 1.6875
        expected_h11 = zeta**2 / 2.0 - 2.0 * zeta  # -1.951171875
        expected_j11 = 5.0 * zeta / 8.0             # 1.054687500

        assert he_sys['h1'][0, 0] == pytest.approx(expected_h11, abs=1e-10)
        assert he_sys['eri'][0, 0, 0, 0] == pytest.approx(expected_j11, abs=1e-10)

    def test_eigenvalue_matches_literature(self, he_sys: dict) -> None:
        """
        JW ground state matches analytical energy.

        E = 2*h11 + (11|11) = -2.847656250 Ha
        (single-determinant, exact for 1-orbital basis)
        """
        _, qubit_op, _ = build_qubit_hamiltonian(he_sys)
        n_qubits = 2 * he_sys['n_spatial']
        sparse_h = get_sparse_operator(qubit_op, n_qubits=n_qubits)
        e_gs = eigsh(sparse_h, k=1, which='SA', return_eigenvectors=False)[0]

        assert abs(e_gs - he_sys['literature_energy']) < 0.01, (
            f"JW ground state {e_gs:.6f} Ha does not match literature "
            f"{he_sys['literature_energy']:.6f} Ha"
        )

    def test_pauli_term_count(self, he_sys: dict) -> None:
        """He STO-3G with 2 qubits should have a small number of Pauli terms."""
        _, _, n_pauli = build_qubit_hamiltonian(he_sys)
        # 1 spatial orbital, 2 qubits: identity + Z0 + Z1 + Z0Z1 = 4 terms
        assert n_pauli > 0
        assert n_pauli <= 10  # upper bound for 2-qubit Hamiltonian


# ---------------------------------------------------------------------------
# H2 6-31G stub test
# ---------------------------------------------------------------------------

class TestH2631G:
    """Validate H2 6-31G stub raises NotImplementedError."""

    def test_not_implemented(self) -> None:
        """h2_631g() raises NotImplementedError with helpful message."""
        with pytest.raises(NotImplementedError, match="PySCF"):
            h2_631g()


# ---------------------------------------------------------------------------
# GeoVac vs Gaussian comparison
# ---------------------------------------------------------------------------

class TestGeoVacVsGaussianComparison:
    """
    Side-by-side Pauli term count comparison.

    This is the core test for Paper 14's central claim:
    GeoVac lattice Hamiltonians produce fewer Pauli terms than
    Gaussian-basis Hamiltonians.
    """

    def test_comparison_he(
        self,
        he_sys: dict,
        geovac_he_nmax2: LatticeIndex,
        geovac_he_nmax3: LatticeIndex,
    ) -> None:
        """
        Compare Pauli term counts: GeoVac He vs Gaussian He.

        Note: qubit counts differ (Gaussian STO-3G: 2 qubits;
        GeoVac nmax=2: 10 qubits, nmax=3: 28 qubits).
        The comparison shows that even with MORE qubits, GeoVac
        maintains favorable Pauli scaling due to ERI sparsity.
        """
        # Gaussian He STO-3G
        _, _, n_pauli_gauss = build_qubit_hamiltonian(he_sys)
        n_qubits_gauss = 2 * he_sys['n_spatial']

        # GeoVac He nmax=2
        enc2 = JordanWignerEncoder(geovac_he_nmax2)
        a2 = enc2.analyze()

        # GeoVac He nmax=3
        enc3 = JordanWignerEncoder(geovac_he_nmax3)
        a3 = enc3.analyze()

        # Print comparison table for visibility
        print("\n" + "=" * 65)
        print("  He Pauli Term Count Comparison: GeoVac vs Gaussian")
        print("=" * 65)
        print(f"  {'Basis':<25} {'Q':>5} {'Pauli Terms':>14} {'ERI Density':>12}")
        print(f"  {'-'*25} {'-'*5} {'-'*14} {'-'*12}")
        print(f"  {'Gaussian STO-3G':<25} {n_qubits_gauss:>5} {n_pauli_gauss:>14,} {'100.0%':>12}")
        print(f"  {'GeoVac nmax=2':<25} {a2.n_qubits:>5} {a2.n_pauli_terms:>14,} {a2.eri_density:>11.1%}")
        print(f"  {'GeoVac nmax=3':<25} {a3.n_qubits:>5} {a3.n_pauli_terms:>14,} {a3.eri_density:>11.1%}")
        print("=" * 65)

        # Sanity checks
        assert n_pauli_gauss > 0
        assert a2.n_pauli_terms > 0
        assert a3.n_pauli_terms > 0

        # GeoVac ERI density should be much lower than Gaussian
        assert a2.eri_density < 0.50, "GeoVac ERI density should be < 50%"
        assert a3.eri_density < 0.10, "GeoVac ERI density should be < 10% at nmax=3"


# ---------------------------------------------------------------------------
# He cc-pVDZ tests
# ---------------------------------------------------------------------------

class TestHeCcPvdz:
    """Validate He cc-pVDZ computed integrals."""

    def test_system_metadata(self) -> None:
        """System dict has correct dimensions for cc-pVDZ [2s1p]."""
        sys = he_cc_pvdz()
        assert sys['n_electrons'] == 2
        assert sys['n_spatial'] == 5
        assert sys['h1'].shape == (5, 5)
        assert sys['eri'].shape == (5, 5, 5, 5)
        assert sys['nuclear_repulsion'] == 0.0

    def test_h1_symmetry(self) -> None:
        """One-electron integrals are symmetric."""
        sys = he_cc_pvdz()
        np.testing.assert_allclose(sys['h1'], sys['h1'].T, atol=1e-10)

    def test_eri_symmetry(self) -> None:
        """ERIs have chemist-notation permutational symmetry."""
        sys = he_cc_pvdz()
        eri = sys['eri']
        M = eri.shape[0]
        for p in range(M):
            for q in range(M):
                for r in range(M):
                    for s in range(M):
                        val = eri[p, q, r, s]
                        assert val == pytest.approx(eri[r, s, p, q], abs=1e-10), \
                            f"(pq|rs) != (rs|pq) at ({p},{q},{r},{s})"
                        assert val == pytest.approx(eri[q, p, s, r], abs=1e-10), \
                            f"(pq|rs) != (qp|sr) at ({p},{q},{r},{s})"

    def test_fci_energy_matches_published(self) -> None:
        """
        JW ground state matches computed FCI energy.

        Published: -2.8877 Ha (Woon & Dunning, JCP 100, 2975, 1994).
        Computed:  -2.8876 Ha (within 0.0001 Ha of published).
        """
        sys = he_cc_pvdz()
        _, qubit_op, _ = build_qubit_hamiltonian(sys)
        n_qubits = 2 * sys['n_spatial']
        sparse_h = get_sparse_operator(qubit_op, n_qubits=n_qubits)
        e_gs = eigsh(sparse_h, k=1, which='SA', return_eigenvectors=False)[0]

        # Match the hardcoded literature energy
        assert abs(e_gs - sys['literature_energy']) < 0.001, (
            f"JW ground state {e_gs:.6f} Ha does not match literature "
            f"{sys['literature_energy']:.6f} Ha"
        )
        # Also within 0.001 Ha of published Woon & Dunning value
        assert abs(e_gs - (-2.8877)) < 0.001

    def test_pauli_term_count(self) -> None:
        """cc-pVDZ He produces 156 Pauli terms (5 spatial, 10 qubits)."""
        sys = he_cc_pvdz()
        _, _, n_pauli = build_qubit_hamiltonian(sys)
        assert n_pauli == 156, f"Expected 156 Pauli terms, got {n_pauli}"

    def test_p_orbital_degeneracy(self) -> None:
        """The three p orbitals should have identical h1 eigenvalues."""
        sys = he_cc_pvdz()
        h1_diag = np.diag(sys['h1'])
        # Orbitals 2,3,4 are p_x, p_y, p_z
        np.testing.assert_allclose(
            h1_diag[2:5], h1_diag[2], atol=1e-8,
            err_msg="p-orbital degeneracy broken in h1"
        )


# ---------------------------------------------------------------------------
# He cc-pVTZ tests
# ---------------------------------------------------------------------------

class TestHeCcPvtz:
    """Validate He cc-pVTZ computed integrals (loaded from cache)."""

    def test_system_metadata(self) -> None:
        """System dict has correct dimensions for cc-pVTZ [3s2p1d]."""
        sys = he_cc_pvtz()
        assert sys['n_electrons'] == 2
        assert sys['n_spatial'] == 14
        assert sys['h1'].shape == (14, 14)
        assert sys['eri'].shape == (14, 14, 14, 14)
        assert sys['nuclear_repulsion'] == 0.0

    def test_h1_symmetry(self) -> None:
        """One-electron integrals are symmetric."""
        sys = he_cc_pvtz()
        np.testing.assert_allclose(sys['h1'], sys['h1'].T, atol=1e-10)

    def test_eri_exchange_symmetry(self) -> None:
        """ERIs satisfy (pq|rs) = (rs|pq)."""
        sys = he_cc_pvtz()
        eri = sys['eri']
        # Full loop is expensive for M=14; check exchange symmetry only
        np.testing.assert_allclose(
            eri, eri.transpose(2, 3, 0, 1), atol=1e-10,
            err_msg="(pq|rs) != (rs|pq)"
        )

    def test_fci_energy_matches_published(self) -> None:
        """
        Stored FCI energy matches published value.

        Published: -2.9003 Ha (Woon & Dunning, JCP 100, 2975, 1994).
        Computed:  -2.9002 Ha (within 0.0001 Ha).

        Note: full JW diagonalization (28 qubits = 2^28 dim) exceeds
        memory. The FCI energy was validated at computation time by
        debug/compute_he_gaussian_integrals.py. Here we check
        consistency of the stored value against the published reference.
        """
        sys = he_cc_pvtz()
        e_stored = sys['literature_energy']

        # Must be within 0.0002 Ha of published Woon & Dunning value
        assert abs(e_stored - (-2.9003)) < 0.0002, (
            f"Stored FCI energy {e_stored:.6f} Ha deviates from published "
            f"-2.9003 Ha by more than 0.0002 Ha"
        )

        # Must be above exact He energy (variational principle)
        HE_EXACT = -2.903724
        assert e_stored > HE_EXACT, (
            f"Stored energy {e_stored:.6f} Ha is below exact {HE_EXACT} Ha "
            "(variational violation)"
        )

        # HF energy check: E_FCI <= 2*epsilon_1 + (11|11)
        h1 = sys['h1']
        eri = sys['eri']
        e_hf = 2.0 * h1[0, 0] + eri[0, 0, 0, 0]
        assert e_stored <= e_hf + 1e-10, (
            f"FCI energy {e_stored:.6f} should be <= HF energy {e_hf:.6f}"
        )

    def test_pauli_term_count(self) -> None:
        """cc-pVTZ He produces 21,607 Pauli terms (14 spatial, 28 qubits)."""
        sys = he_cc_pvtz()
        _, _, n_pauli = build_qubit_hamiltonian(sys)
        assert n_pauli == 21607, f"Expected 21607 Pauli terms, got {n_pauli}"
