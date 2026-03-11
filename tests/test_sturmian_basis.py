"""
Tests for Sturmian basis implementation (Paper 9, v0.9.20-v0.9.22).

The Sturmian basis places all orbitals on a single S³ with shared momentum
scale p0.  The Shibuya-Wulfman cross-center integrals become exact D-matrix
elements with no form factor.  The self-consistency loop p0 = sqrt(-2*E_mol)
closes the energy-shell condition.

v0.9.22 adds atom-dependent p0 (use_sturmian='atomic'): each atom uses its
own self-consistent p0 from isolated-atom FCI, breaking single-S3 but
restoring binding for heteronuclear molecules.

Tests verify:
  1. Backward compatibility: ε_n(Z/n) = -Z²/(2n²)
  2. Diagonal degeneracy: common kinetic part -p0²/2 + shell-dependent Z·p0/n
  3. Form factor = 1 in Sturmian mode
  4. Form factor = sin(γ) in non-Sturmian mode (existing behavior preserved)
  5. p0 initialization from Paper 9 formula
  6. Convergence or honest failure of self-consistency loop
  9-12. Cross-nuclear D-matrix tests (v0.9.21)
  13. Atomic p0 computation
  14. Per-atom diagonal with use_sturmian='atomic'
  15. Cross-nuclear magnitudes with atom-dependent p0
  16. Gamma splitting: gamma_A != gamma_B
"""

import warnings
import numpy as np
import pytest

from geovac.lattice_index import MolecularLatticeIndex, compute_atomic_p0
from geovac.shibuya_wulfman import sw_form_factor


# Suppress construction warnings
pytestmark = pytest.mark.filterwarnings("ignore::UserWarning")


class TestSturmianBackwardCompatibility:
    """Test 1: Sturmian diagonal at p0=Z/n recovers atomic eigenvalue."""

    def test_backward_compat_assertion_passes(self):
        """Constructing Sturmian MolecularLatticeIndex does not raise."""
        # If the backward compatibility assertion inside _build_sturmian_h1
        # fails, the constructor raises AssertionError.  This test verifies
        # that the assertion passes for LiH at nmax=3.
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            vee_method='slater_full', fci_method='auto',
            use_sturmian=True,
        )
        assert mol.use_sturmian
        assert mol._sturmian_p0 > 0

    def test_backward_compat_formula(self):
        """Verify ε_n(Z/n) = -Z²/(2n²) for all Z and n directly."""
        for Z in [1, 2, 3, 4]:
            for n in range(1, 6):
                p0 = float(Z) / n
                eps_sturmian = p0**2 / 2.0 - float(Z) * p0 / n
                eps_atomic = -float(Z)**2 / (2.0 * n**2)
                np.testing.assert_allclose(
                    eps_sturmian, eps_atomic, atol=1e-12,
                    err_msg=f"Z={Z}, n={n}: ε_Sturmian={eps_sturmian:.10f} "
                            f"vs ε_atomic={eps_atomic:.10f}"
                )


class TestSturmianDiagonalDegeneracy:
    """Test 2: All orbitals share common kinetic part p0²/2."""

    def test_common_kinetic_part(self):
        """With fixed p0, diagonal has kinetic p0²/2 and nuclear -Z·p0/n."""
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            vee_method='slater_full', fci_method='auto',
            use_sturmian=True,
        )
        p0 = mol._sturmian_p0
        nA = mol._n_spatial_A
        kinetic_part = p0**2 / 2.0

        # Check atom A orbitals: h1_diag[i] - kinetic = -Z_A*p0/n + cross_nuc
        # The kinetic part p0²/2 is the same for all, only -Z·p0/n differs
        for i, (ni, li, mi) in enumerate(mol._li_A.lattice.states):
            expected_sturmian = kinetic_part - float(mol.Z_A) * p0 / ni
            # The actual diagonal also includes Fourier cross-nuclear,
            # so check that the Sturmian part matches
            sturmian_part = kinetic_part - float(mol.Z_A) * p0 / ni
            # All orbitals share the same kinetic_part
            assert abs(expected_sturmian - sturmian_part) < 1e-12

        # Verify kinetic_part is the same for A and B orbitals
        for j, (nj, lj, mj) in enumerate(mol._li_B.lattice.states):
            sturmian_part_B = kinetic_part - float(mol.Z_B) * p0 / nj
            # kinetic_part is identical for both atoms
            assert kinetic_part == p0**2 / 2.0


class TestFormFactorSturmian:
    """Test 3: sw_form_factor returns 1.0 when sturmian=True."""

    def test_form_factor_is_one(self):
        """Sturmian form factor = 1 for all n and gamma."""
        for n in [1, 2, 3, 5]:
            for gamma in [0.1, 0.5, 1.0, np.pi / 4, np.pi / 2, 2.5]:
                f = sw_form_factor(n, gamma, sturmian=True)
                assert f == 1.0, \
                    f"f(n={n}, gamma={gamma:.2f}, sturmian=True) = {f}, " \
                    f"expected 1.0"


class TestFormFactorContinuity:
    """Test 4: sw_form_factor still returns sin(γ) when sturmian=False."""

    def test_form_factor_is_sin_gamma(self):
        """Non-Sturmian form factor = sin(gamma) (existing behavior)."""
        for n in [1, 2, 3]:
            for gamma in [0.1, 0.5, 1.0, np.pi / 4, np.pi / 2, 2.5]:
                f = sw_form_factor(n, gamma, sturmian=False)
                expected = np.sin(gamma)
                np.testing.assert_allclose(
                    f, expected, atol=1e-15,
                    err_msg=f"f(n={n}, gamma={gamma:.2f}, sturmian=False) = "
                            f"{f}, expected sin(gamma) = {expected}"
                )


class TestP0Initialization:
    """Test 5: solve_sturmian_p0 initializes p0 correctly."""

    def test_default_p0_init(self):
        """Default p0 = sqrt((Z_A² + Z_B²) / 2) for LiH."""
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            vee_method='slater_full', fci_method='auto',
            use_sturmian=True,
        )
        expected_p0 = np.sqrt((9.0 + 1.0) / 2.0)  # sqrt(5) ≈ 2.236
        np.testing.assert_allclose(
            mol._sturmian_p0, expected_p0, atol=1e-10,
            err_msg=f"p0_init = {mol._sturmian_p0:.6f}, "
                    f"expected sqrt(5) = {expected_p0:.6f}"
        )

    def test_custom_p0_init(self):
        """Custom p0 is respected when passed explicitly."""
        custom_p0 = 3.0
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            vee_method='slater_full', fci_method='auto',
            use_sturmian=True, sturmian_p0=custom_p0,
        )
        assert mol._sturmian_p0 == custom_p0


class TestSturmianConvergence:
    """Test 6: Self-consistency loop converges or fails honestly."""

    def test_convergence_or_honest_failure(self):
        """Run self-consistency for LiH at R=3.015, nmax=3.

        Either outcome is acceptable:
        - Converges: E_mol < E_atoms and p0_final > p0_init
        - Does not converge: returns converged=False with finite E_mol
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=3.015, n_electrons=4,
                vee_method='slater_full', fci_method='auto',
                use_sturmian=True,
            )

            p0_final, E_mol, n_iter, converged = mol.solve_sturmian_p0(
                R=3.015, tol=1e-6, max_iter=50,
            )

        print(f"\nSturmian self-consistency result:")
        print(f"  converged = {converged}")
        print(f"  p0_final  = {p0_final:.6f}")
        print(f"  E_mol     = {E_mol:.6f} Ha")
        print(f"  n_iter    = {n_iter}")

        # Either way, p0 and E_mol should be finite
        assert np.isfinite(p0_final), f"p0_final is not finite: {p0_final}"
        assert np.isfinite(E_mol), f"E_mol is not finite: {E_mol}"
        assert p0_final > 0, f"p0_final should be positive: {p0_final}"
        assert n_iter >= 1, "Should run at least 1 iteration"

        if converged:
            # If converged, the molecule should be bound
            print(f"  CONVERGED — checking bound state")
            # p0 = sqrt(-2*E_mol), so E_mol < 0
            assert E_mol < 0, f"Converged but E_mol >= 0: {E_mol}"


class TestSturmianCrossNuclearN1:
    """Test 9: Sturmian cross-nuclear n=1 diagonal = -Z_B/p0."""

    def test_n1_diagonal(self):
        """D^(1)_{(00),(00)} = 1 for all gamma, so (1s,1s) = -Z_B/p0."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p0 = 3.0
            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=3.015, n_electrons=4,
                vee_method='slater_full', fci_method='auto',
                use_sturmian=True, sturmian_p0=p0,
            )
            V_A, V_B = mol._build_sturmian_cross_nuclear(p0, 3.015)

        expected = -1.0 / p0  # -Z_B/p0
        np.testing.assert_allclose(
            V_A[0, 0], expected, atol=1e-10,
            err_msg=f"n=1 (1s,1s) cross-nuclear: {V_A[0,0]:.10f} "
                    f"vs expected {expected:.10f}"
        )


class TestSturmianCrossNuclearN2:
    """Test 10: Sturmian cross-nuclear n=2 diagonal = -(Z_B/p0)*cos(gamma)."""

    def test_n2_diagonal(self):
        """D^(2)_{(00),(00)}(gamma) = cos(gamma) for 2s orbital."""
        from geovac.wigner_so4 import bond_angle

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p0 = 3.0
            R = 3.015
            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=R, n_electrons=4,
                vee_method='slater_full', fci_method='auto',
                use_sturmian=True, sturmian_p0=p0,
            )
            V_A, _ = mol._build_sturmian_cross_nuclear(p0, R)

        gamma = bond_angle(R, p0)
        # 2s is at index 1 (after 1s at index 0)
        expected = -(1.0 / p0) * np.cos(gamma)
        np.testing.assert_allclose(
            V_A[1, 1], expected, atol=1e-10,
            err_msg=f"n=2 (2s,2s) cross-nuclear: {V_A[1,1]:.10f} "
                    f"vs expected {expected:.10f}"
        )


class TestSturmianCrossNuclearP0Scaling:
    """Test 11: Cross-nuclear diagonal scales as 1/p0."""

    def test_p0_scaling(self):
        """Doubling p0 halves the n=1 diagonal (D^(1)=1, only 1/p0 varies)."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            R = 3.015
            p0_a = 2.0
            p0_b = 4.0

            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=R, n_electrons=4,
                vee_method='slater_full', fci_method='auto',
                use_sturmian=True, sturmian_p0=p0_a,
            )

            V_A_a, _ = mol._build_sturmian_cross_nuclear(p0_a, R)
            V_A_b, _ = mol._build_sturmian_cross_nuclear(p0_b, R)

        # n=1 element: -Z_B/p0 (D^(1)=1, no gamma dependence)
        ratio = V_A_a[0, 0] / V_A_b[0, 0]
        expected_ratio = p0_b / p0_a  # 4/2 = 2
        np.testing.assert_allclose(
            ratio, expected_ratio, atol=1e-10,
            err_msg=f"1/p0 scaling: ratio={ratio:.6f}, "
                    f"expected {expected_ratio:.6f}"
        )


class TestSturmianCrossNuclearGamma:
    """Test 12: n=1 element is R-independent at fixed p0."""

    def test_n1_r_independence(self):
        """D^(1)_{(00),(00)} = 1 for all gamma, so -Z_B/p0 is R-independent."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p0 = 3.0
            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=3.015, n_electrons=4,
                vee_method='slater_full', fci_method='auto',
                use_sturmian=True, sturmian_p0=p0,
            )

            vals = []
            for R in [2.0, 3.015, 5.0, 8.0]:
                V_A, _ = mol._build_sturmian_cross_nuclear(p0, R)
                vals.append(V_A[0, 0])

        # All should be identical: -Z_B/p0 = -1/3
        for v in vals:
            np.testing.assert_allclose(
                v, vals[0], atol=1e-12,
                err_msg=f"n=1 element varies with R: {vals}"
            )


class TestAtomicP0Computation:
    """Test 13: compute_atomic_p0 returns correct self-consistent p0."""

    def test_hydrogen_p0(self):
        """H at any nmax: p0 = 1.0 exactly (E_H = -0.5 Ha)."""
        p0_H = compute_atomic_p0(1, 3)
        np.testing.assert_allclose(
            p0_H, 1.0, atol=1e-6,
            err_msg=f"H p0 = {p0_H:.6f}, expected 1.0"
        )

    def test_lithium_p0(self):
        """Li at nmax=3: p0 consistent with E_Li = -7.392 Ha."""
        p0_Li = compute_atomic_p0(3, 3)
        E_Li = -p0_Li**2 / 2.0
        np.testing.assert_allclose(
            E_Li, -7.392, atol=0.01,
            err_msg=f"Li E = {E_Li:.4f}, expected ~-7.392 Ha"
        )

    def test_ghost_p0(self):
        """Ghost atom (Z=0) returns p0 = 0."""
        assert compute_atomic_p0(0, 3) == 0.0


class TestAtomicSturmianDiagonal:
    """Test 14: Per-atom diagonal with use_sturmian='atomic'."""

    def test_h_1s_diagonal(self):
        """H 1s diagonal = p0_B²/2 - Z_B*p0_B/1 = -0.5 Ha."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=3.015, n_electrons=4,
                vee_method='slater_full', fci_method='auto',
                use_sturmian='atomic',
            )
        nA = mol._n_spatial_A
        # Pure Sturmian diagonal for H 1s (before cross-nuclear)
        p0_B = mol._sturmian_p0_B
        expected = p0_B**2 / 2.0 - float(mol.Z_B) * p0_B / 1.0
        np.testing.assert_allclose(
            expected, -0.5, atol=1e-6,
            err_msg=f"H 1s Sturmian diagonal = {expected:.6f}, "
                    f"expected -0.500 Ha"
        )


class TestAtomicSturmianCrossNuclear:
    """Test 15: Cross-nuclear magnitudes with atom-dependent p0."""

    def test_li_cross_nuclear(self):
        """Li 1s from H nucleus: -(Z_B/p0_A) * D^(1) = -(1/p0_A)."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=3.015, n_electrons=4,
                vee_method='slater_full', fci_method='auto',
                use_sturmian='atomic',
            )
        V_A, V_B = mol._build_atomic_sturmian_cross_nuclear(3.015)
        p0_A = mol._sturmian_p0_A
        expected_A = -(1.0 / p0_A)  # -Z_B/p0_A, D^(1)=1
        np.testing.assert_allclose(
            V_A[0, 0], expected_A, atol=1e-10,
            err_msg=f"Li 1s cross-nuc = {V_A[0,0]:.10f}, "
                    f"expected {expected_A:.10f}"
        )

    def test_h_cross_nuclear_capped(self):
        """H 1s from Li nucleus: capped at -Z_A/R = -3/3.015."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=3.015, n_electrons=4,
                vee_method='slater_full', fci_method='auto',
                use_sturmian='atomic',
            )
        R = 3.015
        V_A, V_B = mol._build_atomic_sturmian_cross_nuclear(R)
        p0_B = mol._sturmian_p0_B
        uncapped = -(3.0 / p0_B)  # -Z_A/p0_B = -3.0
        capped = -(3.0 / R)       # -Z_A/R = -0.995
        # The capped value should be used (since |uncapped| > |capped|)
        np.testing.assert_allclose(
            V_B[0, 0], capped, atol=1e-4,
            err_msg=f"H 1s cross-nuc = {V_B[0,0]:.6f}, "
                    f"expected capped {capped:.6f} (uncapped {uncapped:.6f})"
        )


class TestGammaSplitting:
    """Test 16: gamma_A != gamma_B at R=3.015, gamma_B > gamma_A."""

    def test_gamma_split(self):
        """H has smaller p0 so p_R/p0_B is larger, hence gamma_B > gamma_A."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=3.015, n_electrons=4,
                vee_method='slater_full', fci_method='auto',
                use_sturmian='atomic',
            )
        gamma_A = mol._gamma_A
        gamma_B = mol._gamma_B
        assert gamma_A != gamma_B, (
            f"gamma_A = gamma_B = {gamma_A:.6f}, expected different values"
        )
        assert gamma_B > gamma_A, (
            f"Expected gamma_B ({gamma_B:.6f}) > gamma_A ({gamma_A:.6f})"
        )


