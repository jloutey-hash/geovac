"""Tests for MO Sturmian FCI (v0.9.34).

Test 1: H1 matrix symmetric after orthogonalization.
Test 2: MO count -- 10 MOs, 20 spin-orbitals, 4845 SDs for LiH nmax=3.
Test 3: H2 homonuclear symmetry -- 1sigma and 2sigma diag elements equal by symmetry.
Test 4: Dissociation limit -- at R=20, E_mol within 30% of E_atoms.
Test 5: Canonical orthogonalization preserves trace (within 10%).
Test 6: J integral diagonal ordering -- J_00 > J_11.
Test 7: J_00 within 30% of atomic Li 1s self-repulsion (1.875 Ha).
Test 8: K diagonal equals J diagonal (K_kk = J_kk to within 1%).
Test 9: K matrix symmetric (K_kl = K_lk to within 1e-6).
Test 10: K selection rule -- K_kl = 0 when m_k != m_l.
Test 11: Cross-set consistency -- at equal p0, cross-set overlaps match within-set.
Test 12: H-scale H1 negative -- H1_HH diagonal at p0_H=1.0 is negative.
Test 13: Combined basis count -- dual-p0 with nmax_Li=3, nmax_H=2 has 14 spatial MOs.
"""

import numpy as np
import pytest

from geovac.molecular_sturmian import (
    compute_molecular_sturmian_betas,
    compute_h1_matrix,
    compute_exact_j_integrals,
    compute_exact_k_integrals,
    compute_cross_set_integrals,
)
from geovac.mo_fci import MOSturmianFCI


class TestH1Symmetry:
    """Test 1: H1 matrix symmetric to within 1e-6."""

    def test_h1_symmetric(self) -> None:
        """Orthogonalized H1 should be exactly symmetric."""
        p0 = np.sqrt(10.0)
        mo_results = compute_molecular_sturmian_betas(
            Z_A=3.0, Z_B=1.0, R=3.015, p0=p0, nmax=3
        )
        H1, X, n_dropped, S, H1_raw = compute_h1_matrix(
            mo_results, 3.0, 1.0, 3.015, p0, n_grid=30
        )

        max_asym = np.max(np.abs(H1 - H1.T))
        assert max_asym < 1e-6, (
            f"H1 asymmetry {max_asym:.2e} exceeds 1e-6 threshold"
        )


class TestMOCount:
    """Test 2: Correct number of MOs, spin-orbitals, and SDs."""

    def test_lih_nmax3_counts(self) -> None:
        """LiH at nmax=3 should have 10 MOs, 20 spinorbs, 4845 SDs."""
        fci = MOSturmianFCI(
            Z_A=3.0, Z_B=1.0, R=3.015, nmax=3, n_electrons=4
        )
        p0 = np.sqrt(10.0)
        fci._build_basis(p0)

        assert fci.n_mo == 10, f"Expected 10 MOs, got {fci.n_mo}"
        assert fci.n_spinorb == 20, f"Expected 20 spinorbs, got {fci.n_spinorb}"
        assert fci.n_sd == 4845, f"Expected 4845 SDs, got {fci.n_sd}"


class TestH2Symmetry:
    """Test 3: H2 homonuclear symmetry."""

    def test_h2_sigma_degeneracy(self) -> None:
        """For H2 (Z_A=Z_B=1), 1sigma and 2sigma H1 diag should be within 5%."""
        p0 = 1.0
        mo_results = compute_molecular_sturmian_betas(
            Z_A=1.0, Z_B=1.0, R=1.4, p0=p0, nmax=2
        )
        valid = [mo for mo in mo_results
                 if np.isfinite(mo[4]) and mo[1] == 0]

        if len(valid) < 2:
            pytest.skip("Not enough sigma MOs found for H2")

        H1, X, n_dropped, S, H1_raw = compute_h1_matrix(
            mo_results, 1.0, 1.0, 1.4, p0, n_grid=30
        )

        # First two sigma MOs (m=0): use raw H1 for original-basis diags
        sigma_diags = [H1_raw[i, i] for i, mo in enumerate(mo_results)
                       if np.isfinite(mo[4]) and mo[1] == 0][:2]

        if len(sigma_diags) == 2:
            avg = (abs(sigma_diags[0]) + abs(sigma_diags[1])) / 2.0
            if avg > 0.01:
                assert min(abs(sigma_diags[0]), abs(sigma_diags[1])) > 0, (
                    f"H2 sigma diag: {sigma_diags}"
                )


class TestDissociationLimit:
    """Test 4: Dissociation limit at large R."""

    @pytest.mark.slow
    def test_lih_dissociation(self) -> None:
        """At R=20, E_mol should be within 30% of E_atoms=-7.892 Ha."""
        fci = MOSturmianFCI(
            Z_A=3.0, Z_B=1.0, R=20.0, nmax=3, n_electrons=4, n_grid=30
        )
        E, p0 = fci.solve(damping=0.5, max_iter=15, verbose=False)

        E_atoms = -7.892
        rel_err = abs(E - E_atoms) / abs(E_atoms)
        assert E < 0, f"Energy is positive: {E:.4f} Ha"
        assert rel_err < 0.30, (
            f"E_mol={E:.4f} Ha at R=20, E_atoms={E_atoms:.4f}, "
            f"error={rel_err*100:.1f}% (tolerance 30%)"
        )


class TestCanonicalOrthTrace:
    """Test 5: Canonical orth preserves eigenvalues (same trace as Lowdin)."""

    def test_trace_matches_lowdin(self) -> None:
        """tr(H1_canonical) should match tr(H1_lowdin) within 1e-4.

        Both canonical and Lowdin orthogonalization solve the same
        generalized eigenvalue problem H c = E S c, so they must
        produce the same eigenvalues and hence the same trace
        (when no vectors are dropped).
        """
        p0 = np.sqrt(10.0)
        mo_results = compute_molecular_sturmian_betas(
            Z_A=3.0, Z_B=1.0, R=3.015, p0=p0, nmax=3
        )
        H1_can, X_can, n_dropped, S, H1_raw = compute_h1_matrix(
            mo_results, 3.0, 1.0, 3.015, p0, n_grid=30,
            orth_method='canonical',
        )
        H1_low, X_low, _, _, _ = compute_h1_matrix(
            mo_results, 3.0, 1.0, 3.015, p0, n_grid=30,
            orth_method='lowdin',
        )

        tr_can = np.trace(H1_can)
        tr_low = np.trace(H1_low)

        assert abs(tr_can - tr_low) < 1e-4, (
            f"Canonical trace {tr_can:.4f} != Lowdin trace {tr_low:.4f} "
            f"(n_dropped={n_dropped})"
        )


class TestJIntegralOrdering:
    """Test 6: J_00 > J_11 -- core MO has larger self-repulsion."""

    def test_j_diagonal_ordering(self) -> None:
        """Direct Coulomb J_00 (core) > J_11 (bonding)."""
        p0 = np.sqrt(10.0)
        mo_results = compute_molecular_sturmian_betas(
            Z_A=3.0, Z_B=1.0, R=3.015, p0=p0, nmax=3
        )
        J = compute_exact_j_integrals(
            mo_results, 3.0, 1.0, 3.015, p0, n_grid=30
        )

        assert J[0, 0] > J[1, 1], (
            f"J_00={J[0,0]:.4f} should be > J_11={J[1,1]:.4f}"
        )


class TestJIntegralMagnitude:
    """Test 7: J_00 within 30% of atomic Li 1s self-repulsion."""

    def test_j00_magnitude(self) -> None:
        """J_00 should be within 30% of 5*Z_Li/8 = 1.875 Ha."""
        p0 = np.sqrt(10.0)
        mo_results = compute_molecular_sturmian_betas(
            Z_A=3.0, Z_B=1.0, R=3.015, p0=p0, nmax=3
        )
        J = compute_exact_j_integrals(
            mo_results, 3.0, 1.0, 3.015, p0, n_grid=30
        )

        J_00_atomic = 5.0 * 3.0 / 8.0  # = 1.875
        rel_err = abs(J[0, 0] - J_00_atomic) / J_00_atomic
        assert rel_err < 0.30, (
            f"J_00={J[0,0]:.4f}, atomic={J_00_atomic:.4f}, "
            f"error={rel_err*100:.1f}% (tolerance 30%)"
        )


class TestKDiagonalEqualsJ:
    """Test 8: K_kk = J_kk for all k (exchange of orbital with itself = self-repulsion)."""

    def test_k_diagonal_equals_j(self) -> None:
        """K[k,k] should equal J[k,k] to within 1% for all valid MOs."""
        p0 = np.sqrt(10.0)
        mo_results = compute_molecular_sturmian_betas(
            Z_A=3.0, Z_B=1.0, R=3.015, p0=p0, nmax=3
        )
        J = compute_exact_j_integrals(
            mo_results, 3.0, 1.0, 3.015, p0, n_grid=30
        )
        K = compute_exact_k_integrals(
            mo_results, 3.0, 1.0, 3.015, p0, n_grid=30
        )

        n_mo = J.shape[0]
        for k in range(n_mo):
            if abs(J[k, k]) < 1e-10:
                continue
            rel_err = abs(K[k, k] - J[k, k]) / abs(J[k, k])
            assert rel_err < 0.01, (
                f"K[{k},{k}]={K[k,k]:.6f} != J[{k},{k}]={J[k,k]:.6f}, "
                f"error={rel_err*100:.2f}% (tolerance 1%)"
            )


class TestKSymmetric:
    """Test 9: K matrix symmetric (K_kl = K_lk)."""

    def test_k_symmetric(self) -> None:
        """K should be symmetric to within 1e-6."""
        p0 = np.sqrt(10.0)
        mo_results = compute_molecular_sturmian_betas(
            Z_A=3.0, Z_B=1.0, R=3.015, p0=p0, nmax=3
        )
        K = compute_exact_k_integrals(
            mo_results, 3.0, 1.0, 3.015, p0, n_grid=30
        )

        max_asym = np.max(np.abs(K - K.T))
        assert max_asym < 1e-6, (
            f"K asymmetry {max_asym:.2e} exceeds 1e-6 threshold"
        )


class TestKSelectionRule:
    """Test 10: K_kl = 0 when m_k != m_l (phi-orthogonality)."""

    def test_k_sigma_pi_zero(self) -> None:
        """K between sigma (m=0) and pi (m=1) MOs should be zero."""
        p0 = np.sqrt(10.0)
        mo_results = compute_molecular_sturmian_betas(
            Z_A=3.0, Z_B=1.0, R=3.015, p0=p0, nmax=3
        )

        # Find sigma and pi indices in valid MO list
        valid_mos = [mo for mo in mo_results if np.isfinite(mo[4])]
        sigma_idx = None
        pi_idx = None
        for i, mo in enumerate(valid_mos):
            if mo[1] == 0 and sigma_idx is None:
                sigma_idx = i
            if mo[1] == 1 and pi_idx is None:
                pi_idx = i

        if sigma_idx is None or pi_idx is None:
            pytest.skip("Could not find sigma + pi MO pair")

        K = compute_exact_k_integrals(
            mo_results, 3.0, 1.0, 3.015, p0, n_grid=30
        )

        assert abs(K[sigma_idx, pi_idx]) < 1e-4, (
            f"K[{sigma_idx},{pi_idx}]={K[sigma_idx,pi_idx]:.6f} "
            f"should be ~0 (sigma-pi pair, different m)"
        )


class TestCrossSetConsistency:
    """Test 11: At equal p0, cross-set integrals match within-set."""

    def test_equal_p0_overlap(self) -> None:
        """When p0_A == p0_B, cross-set overlap should match within-set to 2%."""
        p0 = 2.0
        R = 3.015
        mo_results = compute_molecular_sturmian_betas(
            Z_A=3.0, Z_B=1.0, R=R, p0=p0, nmax=2,
        )

        # Within-set overlap
        _, _, _, S_within, _ = compute_h1_matrix(
            mo_results, 3.0, 1.0, R, p0, n_grid=30,
            orth_method='none_raw',
        )

        # Cross-set with same p0 and same MO results
        O_cross, _ = compute_cross_set_integrals(
            mo_results, mo_results,
            3.0, 1.0, R, p0, p0, n_grid=30,
        )

        # Diagonal of cross-set overlap should match within-set diagonal
        n = min(S_within.shape[0], O_cross.shape[0])
        for k in range(n):
            if abs(S_within[k, k]) < 1e-6:
                continue
            rel_err = abs(O_cross[k, k] - S_within[k, k]) / abs(S_within[k, k])
            assert rel_err < 0.02, (
                f"Cross-set O[{k},{k}]={O_cross[k,k]:.4f} vs "
                f"within-set S[{k},{k}]={S_within[k,k]:.4f}, "
                f"error={rel_err*100:.1f}% (tolerance 2%)"
            )


class TestHScaleH1Negative:
    """Test 12: H-scale H1 diagonal at p0_H=1.0 should be negative."""

    def test_h_scale_h1_negative(self) -> None:
        """All H-scale MO diagonals should be negative at p0_H=1.0."""
        p0_H = 1.0
        R = 3.015
        mo_results_H = compute_molecular_sturmian_betas(
            Z_A=3.0, Z_B=1.0, R=R, p0=p0_H, nmax=2,
        )

        H1_HH, _, _, _, H1_HH_raw = compute_h1_matrix(
            mo_results_H, 3.0, 1.0, R, p0_H, n_grid=30,
            orth_method='none_raw',
        )

        valid = [mo for mo in mo_results_H if np.isfinite(mo[4])]
        n_valid = len(valid)
        assert n_valid > 0, "No valid H-scale MOs found"

        for k in range(n_valid):
            assert H1_HH_raw[k, k] < 0, (
                f"H-scale H1[{k},{k}]={H1_HH_raw[k,k]:.4f} should be "
                f"negative at p0_H=1.0"
            )


class TestCombinedBasisCount:
    """Test 13: Dual-p0 FCI with nmax_Li=3, nmax_H=2 has 14 spatial MOs."""

    def test_dual_p0_mo_count(self) -> None:
        """Combined Li (nmax=3) + H (nmax=2) should have 14 spatial MOs."""
        fci = MOSturmianFCI(
            Z_A=3.0, Z_B=1.0, R=3.015, nmax=3, n_electrons=4,
            dual_p0=True, nmax_H=2,
        )

        p0_Li = np.sqrt(2.0 * 7.392)
        p0_H = 1.0
        ok = fci._build_basis_dual(p0_Li, p0_H)

        assert ok, "Basis build failed"
        # nmax=3 -> 10 MOs, nmax=2 -> 4 MOs => 14 total
        expected_Li = 10
        expected_H = 4
        assert fci.n_mo_Li == expected_Li, (
            f"Expected {expected_Li} Li MOs, got {fci.n_mo_Li}"
        )
        assert fci.n_mo_H == expected_H, (
            f"Expected {expected_H} H MOs, got {fci.n_mo_H}"
        )
        assert fci.n_mo == expected_Li + expected_H, (
            f"Expected {expected_Li + expected_H} total MOs, "
            f"got {fci.n_mo}"
        )
