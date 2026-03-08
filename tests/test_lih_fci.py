"""
Tests for LiH heteronuclear FCI via MolecularLatticeIndex.

LiH is the first heteronuclear molecule in GeoVac FCI.
Uses same-atom V_ee approximation with cross-atom s-orbital V_ee.
BSSE (Basis Set Superposition Error) is quantified via Boys-Bernardi
counterpoise correction. CP-corrected binding energies are the
physically meaningful quantities.

Experimental references:
    R_eq = 3.015 Bohr (1.595 Angstrom)
    E(LiH, exact NR) = -8.0705 Ha
    E(Li, exact) = -7.4781 Ha
    E(H, exact) = -0.5000 Ha
    Binding energy = 0.0924 Ha (2.515 eV)
"""

import warnings
import numpy as np
import pytest

from geovac.lattice_index import (
    MolecularLatticeIndex, LatticeIndex, compute_bsse_correction,
)


# Suppress V_ee method warnings in tests
pytestmark = pytest.mark.filterwarnings("ignore::UserWarning")


def _build_lih(R: float, nmax: int = 2, n_bridges: int = 10,
               fci_method: str = 'matrix') -> MolecularLatticeIndex:
    """Helper to build LiH at given R with small basis for test speed."""
    return MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=4,
        n_bridges=n_bridges, vee_method='slater_full',
        fci_method=fci_method,
    )


class TestLiHBasic:
    """Basic functionality tests for LiH FCI."""

    def test_lih_runs_at_req(self):
        """LiH FCI completes at experimental R_eq without error."""
        mol = _build_lih(R=3.015, nmax=2)
        eigvals, eigvecs = mol.compute_ground_state(n_states=1)
        E = eigvals[0]
        # Must be below sum of *approximate* atomic energies
        # Li nmax=2 exact h1 gives ~ -7.10, H nmax=2 ~ -0.49
        assert E < -7.5, f"E={E:.4f} Ha, expected below -7.5 Ha"
        assert eigvecs.shape[0] == mol.n_sd

    def test_lih_nuclear_repulsion(self):
        """Nuclear repulsion Z_Li*Z_H/R = 3/R is correctly included."""
        R = 3.0
        mol = _build_lih(R=R, nmax=2)
        expected_vnn = 3.0 / R
        assert abs(mol.V_NN - expected_vnn) < 1e-12, \
            f"V_NN={mol.V_NN}, expected {expected_vnn}"

    def test_lih_state_counts(self):
        """Verify combined basis size: n_spatial = n_A + n_B."""
        mol = _build_lih(R=3.0, nmax=2)
        # nmax=2: 1+4=5 spatial states per atom
        assert mol._n_spatial_A == 5
        assert mol._n_spatial_B == 5
        assert mol._n_spatial == 10
        assert mol.n_sp == 20
        # C(20, 4) = 4845 SDs
        from math import comb
        assert mol.n_sd == comb(20, 4)

    def test_lih_heteronuclear_h1(self):
        """H1 diagonal has different scales for Li (Z=3) vs H (Z=1)."""
        mol = _build_lih(R=5.0, nmax=2)
        # First spatial state is Li 1s: h1 ~ -Z²/2 + cross ≈ -4.5 + small
        # H 1s (index nA): h1 ~ -Z²/2 + cross ≈ -0.5 + small
        nA = mol._n_spatial_A
        li_1s_diag = mol._h1_diag[0]
        h_1s_diag = mol._h1_diag[nA]
        # Li 1s should be much more negative than H 1s
        assert li_1s_diag < h_1s_diag, \
            f"Li 1s ({li_1s_diag:.3f}) should be < H 1s ({h_1s_diag:.3f})"
        assert li_1s_diag < -4.0  # ~ -4.5 for Li 1s
        assert h_1s_diag < 0  # H 1s is negative

    def test_lih_bridges_exist(self):
        """Inter-atomic bridges are created."""
        mol = _build_lih(R=3.0, nmax=2, n_bridges=10)
        assert mol._n_bridges_actual > 0


class TestLiHBinding:
    """Tests for LiH binding properties."""

    def test_lih_binding_positive(self):
        """LiH is bound: E(LiH) < E(Li,nmax) + E(H,nmax) at nmax=2."""
        mol = _build_lih(R=3.0, nmax=2)
        E_lih, _ = mol.compute_ground_state(n_states=1)

        # Compute atomic energies at same basis
        li = LatticeIndex(n_electrons=3, max_n=2, nuclear_charge=3,
                          vee_method='slater_full', h1_method='exact')
        E_li, _ = li.compute_ground_state(n_states=1)

        h = LatticeIndex(n_electrons=1, max_n=2, nuclear_charge=1,
                         vee_method='slater_full', h1_method='exact')
        E_h, _ = h.compute_ground_state(n_states=1)

        E_sep = E_li[0] + E_h[0]
        assert E_lih[0] < E_sep, \
            f"LiH ({E_lih[0]:.4f}) should be bound (< {E_sep:.4f})"

    def test_lih_pes_is_bound(self):
        """PES shows binding: E(R=3) < E(R=8) (molecule more stable than atoms)."""
        mol_eq = _build_lih(R=3.0, nmax=2)
        E_eq, _ = mol_eq.compute_ground_state(n_states=1)
        mol_far = _build_lih(R=8.0, nmax=2)
        E_far, _ = mol_far.compute_ground_state(n_states=1)

        assert E_eq[0] < E_far[0], \
            f"No binding: E(R=3)={E_eq[0]:.4f} >= E(R=8)={E_far[0]:.4f}"

    def test_lih_dissociation_limit(self):
        """E(LiH, R=10) approaches E(Li) + E(H) within 20%."""
        mol_far = _build_lih(R=10.0, nmax=2)
        E_far, _ = mol_far.compute_ground_state(n_states=1)

        # Separated atom energies at nmax=2
        li = LatticeIndex(n_electrons=3, max_n=2, nuclear_charge=3,
                          vee_method='slater_full', h1_method='exact')
        E_li, _ = li.compute_ground_state(n_states=1)
        h = LatticeIndex(n_electrons=1, max_n=2, nuclear_charge=1,
                         vee_method='slater_full', h1_method='exact')
        E_h, _ = h.compute_ground_state(n_states=1)

        E_sep = E_li[0] + E_h[0]
        # At R=10, V_NN = 0.3, cross attraction ~0, so E ≈ E_sep + 0.3
        # Allow 20% tolerance on the binding energy difference
        diff = abs(E_far[0] - E_sep)
        # The difference should be small compared to total energy
        assert diff / abs(E_sep) < 0.20, \
            f"E(R=10)={E_far[0]:.4f} too far from E_sep={E_sep:.4f}"


class TestLiHGhostAtom:
    """Tests for Z=0 ghost atom support (Boys-Bernardi counterpoise)."""

    def test_ghost_atom_no_nuclear_repulsion(self):
        """Ghost atom (Z=0) produces V_NN = 0."""
        mol = MolecularLatticeIndex(
            Z_A=1, Z_B=0, nmax_A=2, nmax_B=2,
            R=3.0, n_electrons=1,
            vee_method='slater_full', fci_method='matrix',
        )
        assert mol.V_NN == 0.0

    def test_ghost_atom_no_bridges(self):
        """Ghost atom produces zero bridges."""
        mol = MolecularLatticeIndex(
            Z_A=1, Z_B=0, nmax_A=2, nmax_B=2,
            R=3.0, n_electrons=1,
            vee_method='slater_full', fci_method='matrix',
        )
        assert mol._n_bridges_actual == 0

    def test_ghost_h_energy_lowered(self):
        """H with ghost orbitals has E <= E(H, own basis) (BSSE lowers energy)."""
        # H in own basis
        h_own = LatticeIndex(
            n_electrons=1, max_n=2, nuclear_charge=1,
            vee_method='slater_full', h1_method='exact', fci_method='matrix',
        )
        E_own = h_own.compute_ground_state(n_states=1)[0][0]

        # H with ghost orbitals
        h_ghost = MolecularLatticeIndex(
            Z_A=0, Z_B=1, nmax_A=2, nmax_B=2,
            R=3.0, n_electrons=1,
            vee_method='slater_full', fci_method='matrix',
        )
        E_ghost = h_ghost.compute_ground_state(n_states=1)[0][0]

        assert E_ghost <= E_own + 1e-10, \
            f"E_ghost ({E_ghost:.6f}) should be <= E_own ({E_own:.6f})"

    def test_ghost_atom_h1_diagonal_zero(self):
        """Ghost atom orbitals have zero h1 diagonal (no nuclear attraction)."""
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=0, nmax_A=2, nmax_B=2,
            R=3.0, n_electrons=3,
            vee_method='slater_full', fci_method='matrix',
        )
        nA = mol._n_spatial_A
        # Ghost atom B orbitals (indices nA..n-1) should have h1_diag = 0
        ghost_diag = mol._h1_diag[nA:]
        np.testing.assert_allclose(ghost_diag, 0.0, atol=1e-15,
                                   err_msg="Ghost atom h1 diagonal should be zero")


class TestLiHCounterpoise:
    """Tests for BSSE counterpoise correction."""

    def test_bsse_is_negative(self):
        """BSSE lowers energy (negative correction) at nmax=2."""
        result = compute_bsse_correction(
            Z_A=3, Z_B=1, nmax_A=2, nmax_B=2, R=3.015,
            n_electrons_A=3, n_electrons_B=1,
            vee_method='slater_full', fci_method='matrix',
        )
        assert result['BSSE'] < 0, \
            f"BSSE should be negative, got {result['BSSE']:.6f}"
        assert result['BSSE_A'] < 0, \
            f"BSSE_A should be negative, got {result['BSSE_A']:.6f}"
        # H BSSE can be very small but should be <= 0
        assert result['BSSE_B'] <= 1e-10, \
            f"BSSE_B should be non-positive, got {result['BSSE_B']:.6f}"

    def test_cp_corrected_binding_positive(self):
        """
        CP-corrected LiH binding energy is positive (molecule is bound).
        Experimental D_e = 0.0924 Ha.

        NOTE: Variational principle is violated at small nmax due to BSSE.
        The molecular basis allows electrons to borrow orbitals from the other
        center, lowering E_mol below the true FCI limit for this basis.
        CP-corrected energies are the reliable quantity.
        See CHANGELOG v0.9.9 for full diagnosis.
        """
        result = compute_bsse_correction(
            Z_A=3, Z_B=1, nmax_A=2, nmax_B=2, R=3.015,
            n_electrons_A=3, n_electrons_B=1,
            vee_method='slater_full', fci_method='matrix',
        )

        # Molecular energy
        mol = _build_lih(R=3.015, nmax=2)
        E_mol = mol.compute_ground_state(n_states=1)[0][0]

        E_ghost_sum = result['E_A_ghost'] + result['E_B_ghost']
        D_e_cp = E_ghost_sum - E_mol

        print(f"\nCP-corrected D_e = {D_e_cp:.4f} Ha (expt: 0.0924)")
        print(f"BSSE = {result['BSSE']:.4f} Ha")

        assert D_e_cp > 0, \
            f"CP-corrected D_e should be positive, got {D_e_cp:.4f} Ha"
