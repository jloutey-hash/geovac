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
    compute_cross_atom_J, compute_cross_atom_K, compute_overlap_element,
)


# Suppress V_ee method warnings in tests
pytestmark = pytest.mark.filterwarnings("ignore::UserWarning")


def _build_lih(R: float, nmax: int = 2, n_bridges: int = 10,
               fci_method: str = 'matrix',
               cross_atom_vee: object = True) -> MolecularLatticeIndex:
    """Helper to build LiH at given R with small basis for test speed."""
    return MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=4,
        n_bridges=n_bridges, vee_method='slater_full',
        fci_method=fci_method, cross_atom_vee=cross_atom_vee,
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


class TestCrossAtomJ:
    """Tests for cross-atom direct Coulomb integrals J_AB(R)."""

    # Reference values computed via compute_cross_atom_J at R=3.015, ZA=3, ZB=1
    # These are the 9 s-s pairs (×2 for ERI symmetry = 18 entries)
    REFERENCE_J = {
        (1, 1): 0.327873,
        (1, 2): 0.185578,
        (1, 3): 0.092579,
        (2, 1): 0.312119,
        (2, 2): 0.180687,
        (2, 3): 0.091118,
        (3, 1): 0.222891,
        (3, 2): 0.157987,
        (3, 3): 0.085702,
    }

    def test_j_asymptotic_zero(self):
        """J_AB(R -> inf) -> 0: orbitals don't interact at large separation."""
        j_far = compute_cross_atom_J(1, 0, 1, 0, R=200.0, ZA=3.0, ZB=1.0)
        assert abs(j_far) < 0.01, \
            f"J(1s,1s,R=200) = {j_far:.6f}, expected ~0"

    def test_j_point_charge_limit(self):
        """J_AB(R) -> 1/R for compact orbitals at large R."""
        R = 100.0
        j = compute_cross_atom_J(1, 0, 1, 0, R=R, ZA=3.0, ZB=1.0)
        # 1s orbitals are compact; at R=100 should be very close to 1/R
        np.testing.assert_allclose(j, 1.0 / R, rtol=0.01,
                                   err_msg="J(1s,1s) should approach 1/R")

    def test_j_reference_values(self):
        """J_AB at R=3.015 matches reference values within 1%."""
        R = 3.015
        for (na, nb), j_ref in self.REFERENCE_J.items():
            j = compute_cross_atom_J(na, 0, nb, 0, R=R, ZA=3.0, ZB=1.0)
            np.testing.assert_allclose(
                j, j_ref, rtol=0.01,
                err_msg=f"J(Li {na}s, H {nb}s) = {j:.6f}, ref = {j_ref:.6f}"
            )

    def test_j_monotone_in_R(self):
        """J_AB decreases with increasing R (less overlap at larger separation)."""
        R_vals = [2.0, 3.0, 5.0, 10.0]
        j_vals = [compute_cross_atom_J(1, 0, 1, 0, R=R, ZA=3.0, ZB=1.0)
                  for R in R_vals]
        for i in range(len(j_vals) - 1):
            assert j_vals[i] > j_vals[i + 1], \
                f"J not monotone: J(R={R_vals[i]})={j_vals[i]:.6f} <= " \
                f"J(R={R_vals[i+1]})={j_vals[i+1]:.6f}"

    def test_j_accepts_nonzero_l(self):
        """compute_cross_atom_J handles l>0 via monopole (spherical average)."""
        j = compute_cross_atom_J(2, 1, 1, 0, R=3.0, ZA=3.0, ZB=1.0)
        # 2p-1s cross-atom J should be positive and less than 1/R
        assert 0 < j < 1.0 / 3.0 + 0.1, \
            f"J(2p,1s,R=3) = {j:.6f}, expected in (0, 0.43)"

    def test_18_cross_atom_entries_in_eri(self):
        """MolecularLatticeIndex generates 36 cross-atom ERIs in s_only mode.

        nmax=3: 3 s-orbitals per atom (1s, 2s, 3s).
        J entries: 3*3 = 9 pairs * 2 symmetry = 18 Coulomb entries.
        K entries: 3*3 = 9 pairs * 2 = 18 exchange entries.
        Total: 36 cross-atom ERI entries.
        """
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            n_bridges=10, vee_method='slater_full',
            fci_method='matrix', cross_atom_vee='s_only',
        )
        nA = mol._n_spatial_A
        n_cross_J = 0
        n_cross_K = 0
        for (a, b, c, d), val in mol._eri.items():
            a_is_A = a < nA
            b_is_A = b < nA
            if a_is_A != b_is_A:
                if a == c and b == d:
                    n_cross_J += 1
                elif a == d and b == c:
                    n_cross_K += 1
        assert n_cross_J == 18, \
            f"Expected 18 cross-atom J entries, got {n_cross_J}"
        assert n_cross_K == 18, \
            f"Expected 18 cross-atom K entries (s-s only), got {n_cross_K}"


class TestOverlapMatrix:
    """Tests for cross-atom overlap integrals and Lowdin orthogonalization."""

    def test_same_center_overlap_is_one(self):
        """Overlap of 1s orbital with itself on the same center = 1."""
        s = compute_overlap_element(1, 0, 1, 0, ZA=1.0, ZB=1.0, R=0.0)
        np.testing.assert_allclose(s, 1.0, atol=1e-6,
                                   err_msg="Self-overlap should be 1")

    def test_overlap_h_h_1s(self):
        """H-H 1s overlap at R=1.4 matches STO analytical value."""
        R = 1.4
        s = compute_overlap_element(1, 0, 1, 0, ZA=1.0, ZB=1.0, R=R)
        # Analytical: S = e^{-R}(1 + R + R²/3)
        s_exact = np.exp(-R) * (1.0 + R + R**2 / 3.0)
        np.testing.assert_allclose(s, s_exact, rtol=0.02,
                                   err_msg=f"H-H 1s overlap: got {s:.4f}, "
                                           f"expected {s_exact:.4f}")

    def test_overlap_li_h_positive(self):
        """Li-H 1s overlap at R=3.015 is positive and < 1."""
        s = compute_overlap_element(1, 0, 1, 0, ZA=3.0, ZB=1.0, R=3.015)
        assert 0 < s < 1, f"Cross-atom overlap should be in (0,1), got {s:.4f}"

    def test_overlap_decreases_with_R(self):
        """Overlap decreases with internuclear distance."""
        s_near = compute_overlap_element(1, 0, 1, 0, ZA=3.0, ZB=1.0, R=2.0)
        s_far = compute_overlap_element(1, 0, 1, 0, ZA=3.0, ZB=1.0, R=5.0)
        assert s_near > s_far, \
            f"Overlap should decrease: S(R=2)={s_near:.4f} <= S(R=5)={s_far:.4f}"

    def test_overlap_rejects_nonzero_l(self):
        """compute_overlap_element raises NotImplementedError for l>0."""
        with pytest.raises(NotImplementedError):
            compute_overlap_element(2, 1, 1, 0, ZA=1.0, ZB=1.0, R=3.0)

    def test_overlap_matrix_identity_blocks(self):
        """Overlap matrix has identity diagonal blocks (same-atom is orthonormal)."""
        mol = _build_lih(R=3.015, nmax=2)
        S = mol._compute_overlap_matrix()
        nA = mol._n_spatial_A
        np.testing.assert_allclose(S[:nA, :nA], np.eye(nA), atol=1e-12,
                                   err_msg="Same-atom block A should be identity")
        nB = mol._n_spatial_B
        np.testing.assert_allclose(S[nA:, nA:], np.eye(nB), atol=1e-12,
                                   err_msg="Same-atom block B should be identity")

    def test_overlap_matrix_symmetric(self):
        """Overlap matrix is symmetric."""
        mol = _build_lih(R=3.015, nmax=2)
        S = mol._compute_overlap_matrix()
        np.testing.assert_allclose(S, S.T, atol=1e-12,
                                   err_msg="Overlap matrix should be symmetric")


class TestLowdinOrthogonalization:
    """Tests for Lowdin symmetric orthogonalization of molecular basis."""

    def test_lowdin_runs_without_error(self):
        """MolecularLatticeIndex with orthogonalize=True completes."""
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=1, nmax_A=2, nmax_B=2,
            R=3.015, n_electrons=4,
            n_bridges=10, vee_method='slater_full',
            fci_method='matrix', orthogonalize=True,
        )
        eigvals, _ = mol.compute_ground_state(n_states=1)
        assert eigvals[0] < -7.0, f"E={eigvals[0]:.4f}, expected below -7.0"

    def test_lowdin_preserves_basis_size(self):
        """Lowdin transform preserves spatial basis size (no functions dropped)."""
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=1, nmax_A=2, nmax_B=2,
            R=3.015, n_electrons=4,
            n_bridges=10, vee_method='slater_full',
            fci_method='matrix', orthogonalize=True,
        )
        assert mol._n_spatial == 10  # 5+5 unchanged
        assert mol.n_sp == 20

    def test_lowdin_reduces_bsse(self):
        """BSSE with orthogonalization should be smaller than without.

        Lowdin removes basis linear dependencies that allow electrons to
        borrow orbitals from the other center (the root cause of BSSE).
        """
        # BSSE without orthogonalization
        result_no = compute_bsse_correction(
            Z_A=3, Z_B=1, nmax_A=2, nmax_B=2, R=3.015,
            n_electrons_A=3, n_electrons_B=1,
            vee_method='slater_full', fci_method='matrix',
            orthogonalize=False,
        )
        # BSSE with orthogonalization
        result_orth = compute_bsse_correction(
            Z_A=3, Z_B=1, nmax_A=2, nmax_B=2, R=3.015,
            n_electrons_A=3, n_electrons_B=1,
            vee_method='slater_full', fci_method='matrix',
            orthogonalize=True,
        )

        print(f"\nBSSE without orth: {result_no['BSSE']:.6f} Ha")
        print(f"BSSE with orth:    {result_orth['BSSE']:.6f} Ha")

        # Orthogonalized BSSE should be smaller in magnitude
        assert abs(result_orth['BSSE']) <= abs(result_no['BSSE']) + 1e-6, \
            f"|BSSE_orth| ({abs(result_orth['BSSE']):.6f}) should be <= " \
            f"|BSSE_no| ({abs(result_no['BSSE']):.6f})"

    def test_ghost_energy_unchanged_with_lowdin(self):
        """Ghost atom with no cross-atom overlap gives same energy."""
        # At very large R, cross-atom overlap → 0, so Lowdin = identity
        mol_no = MolecularLatticeIndex(
            Z_A=1, Z_B=0, nmax_A=2, nmax_B=2,
            R=50.0, n_electrons=1,
            vee_method='slater_full', fci_method='matrix',
            orthogonalize=False,
        )
        E_no = mol_no.compute_ground_state(n_states=1)[0][0]

        mol_orth = MolecularLatticeIndex(
            Z_A=1, Z_B=0, nmax_A=2, nmax_B=2,
            R=50.0, n_electrons=1,
            vee_method='slater_full', fci_method='matrix',
            orthogonalize=True,
        )
        E_orth = mol_orth.compute_ground_state(n_states=1)[0][0]

        np.testing.assert_allclose(E_orth, E_no, atol=1e-8,
                                   err_msg="Ghost at large R should be unchanged")


class TestCrossAtomExchange:
    """Tests for cross-atom exchange integrals via Mulliken approximation."""

    def test_k_zero_at_large_R(self):
        """K_AB -> 0 as R -> inf (overlap vanishes)."""
        mol = _build_lih(R=100.0, nmax=2)
        nA = mol._n_spatial_A
        for (a, b, c, d), val in mol._eri.items():
            a_on_A = a < nA
            b_on_A = b < nA
            # Exchange pattern: (a,b,b,a) with a,b on different atoms
            if a_on_A != b_on_A and a == d and b == c:
                assert abs(val) < 1e-6, \
                    f"Exchange ERI {(a,b,c,d)} = {val:.8f} should be ~0 at R=100"

    def test_k_positive(self):
        """K_AB > 0 at equilibrium (exchange integrals are positive)."""
        mol = _build_lih(R=3.015, nmax=2)
        nA = mol._n_spatial_A
        found_exchange = False
        for (a, b, c, d), val in mol._eri.items():
            a_on_A = a < nA
            b_on_A = b < nA
            if a_on_A != b_on_A and a == d and b == c:
                assert val > 0, \
                    f"Exchange ERI {(a,b,c,d)} = {val:.6f} should be > 0"
                found_exchange = True
        assert found_exchange, "No cross-atom exchange ERIs found"

    def test_k_symmetric(self):
        """K(aA, bB) stored as both (a,b,b,a) and (b,a,a,b)."""
        mol = _build_lih(R=3.015, nmax=2)
        nA = mol._n_spatial_A
        for (a, b, c, d), val in mol._eri.items():
            a_on_A = a < nA
            b_on_A = b < nA
            if a_on_A != b_on_A and a == d and b == c:
                partner = mol._eri.get((b, a, a, b), None)
                assert partner is not None, \
                    f"Missing symmetric partner for exchange ERI {(a,b,c,d)}"
                np.testing.assert_allclose(val, partner, atol=1e-12)

    def test_k_smaller_than_j(self):
        """K_AB < J_AB for same (a, b) pair (exchange < Coulomb)."""
        mol = _build_lih(R=3.015, nmax=2)
        nA = mol._n_spatial_A
        for (a, b, c, d), val in mol._eri.items():
            a_on_A = a < nA
            b_on_A = b < nA
            if a_on_A != b_on_A and a == d and b == c:
                j_val = mol._eri.get((a, b, a, b), 0.0)
                assert val < j_val, \
                    f"K({a},{b}) = {val:.6f} >= J({a},{b}) = {j_val:.6f}"

    def test_cross_atom_eri_counts(self):
        """Cross-atom ERI counts with all-l default (v0.9.35).

        nmax=2: 5 spatial states per atom (1s, 2s, 2p_{-1,0,1}).
        J: 5*5 = 25 pairs * 2 symmetry = 50 Coulomb entries (all l).
        K: 2*2 = 4 s-s pairs * 2 = 8 exchange entries (s-only).
        """
        mol = _build_lih(R=3.015, nmax=2)
        nA = mol._n_spatial_A
        n_coulomb = 0
        n_exchange = 0
        for (a, b, c, d) in mol._eri:
            a_on_A = a < nA
            b_on_A = b < nA
            if a_on_A != b_on_A:
                if a == c and b == d:
                    n_coulomb += 1
                elif a == d and b == c:
                    n_exchange += 1
        assert n_coulomb == 50, \
            f"Expected 50 cross-atom J entries (all l), got {n_coulomb}"
        assert n_exchange == 8, \
            f"Expected 8 cross-atom K entries (s-s only), got {n_exchange}"

    def test_mulliken_formula_manual(self):
        """Verify Mulliken formula: K = S^2 * (F0_aa + F0_bb) / 2."""
        R = 3.015
        ZA, ZB = 3.0, 1.0
        S = compute_overlap_element(1, 0, 1, 0, ZA, ZB, R)

        mol = _build_lih(R=R, nmax=2)
        # F0(Li 1s, Li 1s) = same-atom ERI (0, 0, 0, 0) on atom A
        f0_a = mol._li_A._eri.get((0, 0, 0, 0), 0.0)
        # F0(H 1s, H 1s) = same-atom ERI (0, 0, 0, 0) on atom B
        f0_b = mol._li_B._eri.get((0, 0, 0, 0), 0.0)

        k_expected = S * S * (f0_a + f0_b) / 2.0

        nA = mol._n_spatial_A
        # Exchange ERI for Li 1s (idx=0) - H 1s (idx=nA): stored as (0, nA, nA, 0)
        k_actual = mol._eri.get((0, nA, nA, 0), 0.0)

        assert k_expected > 0, f"Expected K should be positive, got {k_expected}"
        np.testing.assert_allclose(
            k_actual, k_expected, rtol=1e-6,
            err_msg=f"K(1s,1s) = {k_actual:.8f}, expected {k_expected:.8f}"
        )

    def test_energy_lowered_by_exchange(self):
        """Exchange lowers energy: E(with K) < E(without K).

        Exchange enters Slater-Condon as -delta(sigma)*<pq|qp>.
        Since K > 0, this is a negative contribution to same-spin pairs,
        lowering the total energy.
        """
        mol = _build_lih(R=3.015, nmax=2)
        E_with_K, _ = mol.compute_ground_state(n_states=1)

        # Zero out cross-atom exchange ERIs and re-solve
        nA = mol._n_spatial_A
        exchange_keys = []
        for (a, b, c, d) in list(mol._eri.keys()):
            a_on_A = a < nA
            b_on_A = b < nA
            if a_on_A != b_on_A and a == d and b == c:
                exchange_keys.append((a, b, c, d))

        assert len(exchange_keys) > 0, "No exchange ERIs to remove"

        for key in exchange_keys:
            del mol._eri[key]

        # Re-assemble and solve without exchange
        H_no_K = mol.assemble_hamiltonian()
        from scipy.sparse.linalg import eigsh
        eigvals_no_K, _ = eigsh(H_no_K, k=1, which='SA',
                                v0=np.random.RandomState(42).randn(H_no_K.shape[0]))
        E_no_K = eigvals_no_K[0] + mol.V_NN

        assert E_with_K[0] < E_no_K, \
            f"E(with K) = {E_with_K[0]:.6f} should be < E(without K) = {E_no_K:.6f}"

    def test_energy_bound(self):
        """LiH is still bound after adding exchange."""
        mol = _build_lih(R=3.015, nmax=2)
        E_mol, _ = mol.compute_ground_state(n_states=1)

        # Separated atoms
        li = LatticeIndex(n_electrons=3, max_n=2, nuclear_charge=3,
                          vee_method='slater_full', h1_method='exact')
        E_li = li.compute_ground_state(n_states=1)[0][0]
        h = LatticeIndex(n_electrons=1, max_n=2, nuclear_charge=1,
                         vee_method='slater_full', h1_method='exact')
        E_h = h.compute_ground_state(n_states=1)[0][0]
        E_sep = E_li + E_h

        assert E_mol[0] < E_sep, \
            f"E_mol = {E_mol[0]:.6f} should be < E_sep = {E_sep:.6f}"


# ===========================================================================
# Energy decomposition tests (v0.9.24)
# ===========================================================================

class TestEnergyDecomposition:
    """Tests for Hamiltonian term decomposition instrumentation."""

    @pytest.fixture(scope="class")
    def mol_and_civec(self):
        """Build MolecularLatticeIndex at R=3.015 and solve FCI (nmax=2)."""
        mol = _build_lih(R=3.015, nmax=2)
        eigvals, eigvecs = mol.compute_ground_state(n_states=1)
        return mol, eigvecs[:, 0], eigvals[0]

    def test_decomposition_sum(self, mol_and_civec):
        """Verify decomposed components sum to E_total within 1e-5 Ha."""
        mol, civec, E_fci = mol_and_civec
        decomp = mol.decompose_energy(civec, E_fci)

        component_sum = (decomp['T'] + decomp['V_nA'] + decomp['V_nB']
                         + decomp['V_cross_A'] + decomp['V_cross_B']
                         + decomp['V_bridge'] + decomp['V_ee']
                         + decomp['V_NN'])

        assert abs(component_sum - decomp['E_total']) < 1e-5, \
            (f"Component sum {component_sum:.8f} != "
             f"E_total {decomp['E_total']:.8f}")

        # Also verify E_total matches FCI eigenvalue
        assert abs(decomp['E_total'] - E_fci) < 1e-5, \
            (f"E_total {decomp['E_total']:.8f} != "
             f"FCI eigenvalue {E_fci:.8f}")

    def test_virial_ratio(self, mol_and_civec):
        """Verify virial ratio is finite and document its value.

        NOTE: In this Hamiltonian, V_nA = -Z^2/(2n^2) contains both kinetic
        and nuclear attraction (exact atomic eigenvalues). The graph Laplacian
        hopping 'T' is NOT the physical kinetic energy, so the standard virial
        theorem (-2T/V = 1) does not apply to this decomposition. The test
        verifies the decomposition is self-consistent (finite, non-NaN).
        """
        mol, civec, E_fci = mol_and_civec
        decomp = mol.decompose_energy(civec, E_fci)

        V_total = (decomp['V_nA'] + decomp['V_nB']
                   + decomp['V_cross_A'] + decomp['V_cross_B']
                   + decomp['V_ee'] + decomp['V_NN'])

        assert abs(V_total) > 1e-10, f"V_total is near zero: {V_total}"
        virial = -2.0 * decomp['T'] / V_total

        assert np.isfinite(virial), f"Virial ratio is not finite: {virial}"
        # Document the actual value (expected ~0 because T is graph hopping,
        # not physical kinetic energy)
        print(f"\n  Virial ratio -2<T_graph>/<V> = {virial:.4f} "
              f"(T_graph={decomp['T']:.6f}, V_total={V_total:.6f})")


# ===========================================================================
# Cross-atom V_ee l>0 extension tests (v0.9.35)
# ===========================================================================

class TestCrossAtomVeeAllL:
    """Tests for cross-atom V_ee with all (n,l) orbital pairs (v0.9.35)."""

    def test_j_cross_p_orbital_physical_range(self):
        """J(2p_Li, 1s_H) at R=3.015 is positive and < 1/R."""
        j = compute_cross_atom_J(2, 1, 1, 0, R=3.015, ZA=3.0, ZB=1.0)
        assert 0 < j < 1.0 / 3.015 + 0.05, \
            f"J(2p,1s,R=3.015) = {j:.6f}, expected in (0, {1/3.015 + 0.05:.3f})"
        print(f"\n  J(Li 2p, H 1s, R=3.015) = {j:.6f} Ha")

    def test_j_cross_monotone_with_l(self):
        """J monotonically decreases with R for l>0 pairs too."""
        R_vals = [2.0, 3.015, 6.0]
        j_vals = [compute_cross_atom_J(2, 1, 1, 0, R=R, ZA=3.0, ZB=1.0)
                  for R in R_vals]
        for i in range(len(j_vals) - 1):
            assert j_vals[i] > j_vals[i + 1], \
                f"J(2p,1s) not monotone: J(R={R_vals[i]})={j_vals[i]:.6f} " \
                f"<= J(R={R_vals[i+1]})={j_vals[i+1]:.6f}"
        print(f"\n  J(2p,1s): R=2.0→{j_vals[0]:.4f}, "
              f"R=3.015→{j_vals[1]:.4f}, R=6.0→{j_vals[2]:.4f}")

    def test_cross_vee_raises_energy_at_short_R(self):
        """Cross-atom V_ee (all l) raises energy vs disabled at R=2.0.

        Adding cross-atom electron-electron repulsion makes the molecule
        less bound at short R (more positive energy).
        """
        mol_on = _build_lih(R=2.0, nmax=2, cross_atom_vee=True)
        E_on, _ = mol_on.compute_ground_state(n_states=1)

        mol_off = _build_lih(R=2.0, nmax=2, cross_atom_vee=False)
        E_off, _ = mol_off.compute_ground_state(n_states=1)

        assert E_on[0] > E_off[0], \
            f"E(cross_vee=True) = {E_on[0]:.6f} should be > " \
            f"E(cross_vee=False) = {E_off[0]:.6f} at R=2.0"
        print(f"\n  R=2.0: E(all_l)={E_on[0]:.4f}, E(off)={E_off[0]:.4f}, "
              f"delta={E_on[0] - E_off[0]:.4f} Ha")
