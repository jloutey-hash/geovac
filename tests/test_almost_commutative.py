"""Tests for the almost-commutative extension (Sprint H1).

Sprint H1 deliverable: minimal electroweak almost-commutative extension
T = T_GV ⊗ T_F with A_F = C ⊕ H, H_F = matter (+) antimatter = C^8,
and an 8x8 Hermitian D_F with off-diagonal C ↔ H Yukawa block in each sector.

Architecture (PI directive 2026-05-06):
- Use TRUTHFUL CH only.
- H_F is doubled into matter/antimatter following Connes-Marcolli 2008 Ch. 13.
- J_F swaps matter <-> antimatter; J_F^2 = +I (KO-dim 6 of A_F).
- KO(combined) = 3 + 6 = 9 ≡ 1 (mod 8); J^2 = -I, JD = +DJ.

Test plan
=========

Section A. Finite electroweak triple (A_F, H_F, D_F, J_F).
Section B. Algebra action faithfulness and *-representation.
Section C. Combined Dirac structure.
Section D. Real structure J = J_GV ⊗ J_F and Connes axioms.
Section E. Inner fluctuations: gauge and Higgs sectors.
Section F. Falsifier check (memo §5).
Section G. Spectral action / Yukawa dependence sanity.
Section H. Cross-validation against memo §5 falsifier.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.almost_commutative import (
    SIGMA_0,
    SIGMA_1,
    SIGMA_2,
    SIGMA_3,
    AlmostCommutativeTriple,
    ElectroweakFiniteTriple,
    gauge_only_triple,
    heat_kernel_trace,
    minimal_electroweak_triple,
    quaternion_to_matrix,
    spectral_action_eigenvalues,
    yukawa_dependence_scan,
)


# ===========================================================================
# Section A. Finite electroweak triple
# ===========================================================================


class TestFiniteTriple:
    def test_default_yukawa_zero(self):
        T = ElectroweakFiniteTriple()
        assert T.yukawa_e == 0.0
        assert T.yukawa_nu == 0.0
        assert T.dim_H_F == 8
        assert T.dim_matter == 4

    def test_dirac_F_zero_yukawa(self):
        T = ElectroweakFiniteTriple()
        D_F = T.dirac_F()
        assert D_F.shape == (8, 8)
        assert np.allclose(D_F, 0.0)

    def test_dirac_F_block_structure(self):
        T = ElectroweakFiniteTriple(yukawa_e=0.1, yukawa_nu=0.05)
        D_F = T.dirac_F()
        # Matter sector (0:4, 0:4): off-diag L<->R Yukawa
        M = T.matter_dirac()
        assert np.allclose(D_F[0:4, 0:4], M)
        # Antimatter sector (4:8, 4:8): conj(M)
        assert np.allclose(D_F[4:8, 4:8], np.conj(M))
        # Off-diagonal matter<->antimatter blocks: zero
        assert np.allclose(D_F[0:4, 4:8], 0.0)
        assert np.allclose(D_F[4:8, 0:4], 0.0)

    def test_dirac_F_hermitian(self):
        T = ElectroweakFiniteTriple(yukawa_e=0.1 + 0.02j, yukawa_nu=0.05 - 0.01j)
        D_F = T.dirac_F()
        assert np.allclose(D_F, D_F.conj().T, atol=1e-14)

    def test_dirac_F_eigenvalues(self):
        # Eigenvalues are +/- y_nu, +/- y_e, each doubled (matter + antimatter).
        T = ElectroweakFiniteTriple(yukawa_e=0.3, yukawa_nu=0.2)
        D_F = T.dirac_F()
        eigs = sorted(np.linalg.eigvalsh(D_F))
        # Real Yukawa => antimatter conj is the same; spectrum is doubled.
        expected = sorted([-0.3, -0.3, -0.2, -0.2, 0.2, 0.2, 0.3, 0.3])
        assert np.allclose(eigs, expected, atol=1e-14)

    def test_J_F_unitary(self):
        T = ElectroweakFiniteTriple()
        U_F = T.real_structure_F()
        I = np.eye(8, dtype=np.complex128)
        assert np.allclose(U_F @ U_F.conj().T, I, atol=1e-14)

    def test_J_F_squared_plus_I(self):
        # Matter/antimatter doubling: J_F^2 = +I (KO-dim 6 of A_F)
        T = ElectroweakFiniteTriple()
        U_F = T.real_structure_F()
        J2 = U_F @ np.conj(U_F)
        I = np.eye(8, dtype=np.complex128)
        assert np.allclose(J2, +I, atol=1e-14)

    def test_J_F_swaps_matter_antimatter(self):
        T = ElectroweakFiniteTriple()
        U_F = T.real_structure_F()
        # A matter state |1, 0, ..., 0> should map to |0, 0, 0, 0, 1, 0, 0, 0>
        e0 = np.zeros(8, dtype=np.complex128)
        e0[0] = 1.0
        # J_F(e0) = U_F @ conj(e0) = U_F[:, 0]
        Je0 = U_F @ np.conj(e0)
        expected = np.zeros(8, dtype=np.complex128)
        expected[4] = 1.0
        assert np.allclose(Je0, expected)

    def test_J_F_D_F_commutation(self):
        # J_F D_F = +D_F J_F (KO-6 sign)
        T = ElectroweakFiniteTriple(yukawa_e=0.3, yukawa_nu=0.2)
        D_F = T.dirac_F()
        U_F = T.real_structure_F()
        # Antilinear: J_F D_F psi = U_F conj(D_F) conj(psi)
        # D_F J_F psi = D_F U_F conj(psi)
        # Equality: U_F conj(D_F) = D_F U_F
        LHS = U_F @ np.conj(D_F)
        RHS = D_F @ U_F
        err = float(np.max(np.abs(LHS - RHS)))
        # For real Yukawa, D_F is real so conj(D_F) = D_F.
        # Verify the identity exactly
        assert err < 1e-14, f"J_F D_F = D_F J_F failed at {err}"


# ===========================================================================
# Section A2. Chirality γ_F (Sprint G3-B, 2026-05-06)
# ===========================================================================
#
# Audit of γ_F on H_F. The H1 module is tuned to KO-dim 6 of T_F per
# Connes-Marcolli SM (CCM 2007, Connes-Marcolli 2008 Table 13.1):
#   γ_F² = I, {γ_F, D_F} = 0, {J_F, γ_F} = 0 (anticommutes), γ_F π_F γ_F = π_F.
#
# Convention deviation from the original G3-B directive ("γ_F identical on
# both copies" → KO-dim 0): module already documents and tests J_F at
# KO-dim 6, so γ_F^antimatter is sign-flipped to give the required
# {J_F, γ_F} = 0. Documented in the chirality_F() docstring and the
# audit memo debug/g3b_chirality_F_audit_memo.md.


class TestChiralityF:
    """γ_F on H_F: KO-dim 6 of T_F (Connes-Marcolli SM)."""

    def test_chirality_F_shape(self):
        T = ElectroweakFiniteTriple()
        G = T.chirality_F()
        assert G.shape == (8, 8)

    def test_chirality_F_diagonal(self):
        T = ElectroweakFiniteTriple()
        G = T.chirality_F()
        # γ_F is diagonal: off-diagonal entries vanish
        assert np.allclose(G - np.diag(np.diag(G)), 0.0)

    def test_chirality_F_matter_block(self):
        # Matter sector (0:4, 0:4) = diag(+1, +1, −1, −1):
        #   indices 0,1 (L, ℍ-doublet): +1
        #   indices 2,3 (R, ℂ-singlets): −1
        T = ElectroweakFiniteTriple()
        G = T.chirality_F()
        expected_matter = np.diag([+1.0, +1.0, -1.0, -1.0]).astype(np.complex128)
        assert np.allclose(G[0:4, 0:4], expected_matter)

    def test_chirality_F_antimatter_block(self):
        # Antimatter sector (4:8, 4:8) = − γ_F^matter (KO-dim 6 sign):
        T = ElectroweakFiniteTriple()
        G = T.chirality_F()
        expected_antimatter = np.diag([-1.0, -1.0, +1.0, +1.0]).astype(np.complex128)
        assert np.allclose(G[4:8, 4:8], expected_antimatter)

    def test_chirality_F_no_cross_blocks(self):
        T = ElectroweakFiniteTriple()
        G = T.chirality_F()
        assert np.allclose(G[0:4, 4:8], 0.0)
        assert np.allclose(G[4:8, 0:4], 0.0)

    def test_chirality_F_squared_identity(self):
        # γ_F² = I  (Z_2 grading axiom)
        T = ElectroweakFiniteTriple()
        G = T.chirality_F()
        assert np.allclose(G @ G, np.eye(8, dtype=np.complex128), atol=1e-14)

    def test_chirality_F_hermitian(self):
        T = ElectroweakFiniteTriple()
        G = T.chirality_F()
        assert np.allclose(G, G.conj().T, atol=1e-14)

    def test_anticommutes_with_DF_zero_yukawa(self):
        # {γ_F, D_F} = 0; trivially satisfied at Y=0 (D_F=0)
        T = ElectroweakFiniteTriple()
        G = T.chirality_F()
        D_F = T.dirac_F()
        anticomm = G @ D_F + D_F @ G
        assert np.allclose(anticomm, 0.0, atol=1e-14)

    def test_anticommutes_with_DF_nonzero_real_yukawa(self):
        # {γ_F, D_F} = 0 when Yukawa is real and non-zero
        T = ElectroweakFiniteTriple(yukawa_e=0.3, yukawa_nu=0.2)
        G = T.chirality_F()
        D_F = T.dirac_F()
        anticomm = G @ D_F + D_F @ G
        err = float(np.max(np.abs(anticomm)))
        assert err < 1e-14, f"{{γ_F, D_F}} residual: {err:.3e}"

    def test_anticommutes_with_DF_complex_yukawa(self):
        # Complex Yukawa: D_F is still off-diag L↔R in each sector
        T = ElectroweakFiniteTriple(
            yukawa_e=0.3 - 0.1j, yukawa_nu=0.05 + 0.02j
        )
        G = T.chirality_F()
        D_F = T.dirac_F()
        anticomm = G @ D_F + D_F @ G
        err = float(np.max(np.abs(anticomm)))
        assert err < 1e-14, f"{{γ_F, D_F}} (complex Y) residual: {err:.3e}"

    def test_algebra_elements_are_even(self):
        # γ_F π_F(a) γ_F = π_F(a) for any a ∈ A_F.
        # Algebra is "even" because π_F is block-diagonal in (L, R) and
        # γ_F is also block-diagonal in (L, R) within each sector.
        T = ElectroweakFiniteTriple()
        G = T.chirality_F()
        for lam, q in [
            (1.0 + 0.0j, (1.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j)),
            (0.5 + 0.3j, (0.4 + 0.0j, 0.2 + 0.0j, 0.1 + 0.0j, 0.05 + 0.0j)),
            (2.0 - 0.7j, (0.0 + 0.0j, 1.0 + 0.0j, 0.5 + 0.0j, 0.2 + 0.0j)),
        ]:
            a = T.algebra_action(lam, q)
            conj_a = G @ a @ G
            err = float(np.linalg.norm(conj_a - a))
            assert err < 1e-14, f"γ_F π_F γ_F − π_F: {err:.3e} for (λ, q)=({lam}, {q})"

    def test_J_F_anticommutes_with_gamma_F_KO6(self):
        # KO-dim 6 standard: {J_F, γ_F} = 0 (anticommutes).
        # Antilinear J_F = U_F K, so:
        #   J_F γ_F psi = U_F conj(γ_F psi) = U_F γ_F^* conj(psi)
        # γ_F is real diagonal, so γ_F^* = γ_F.
        # γ_F J_F psi = γ_F U_F conj(psi).
        # {J_F, γ_F} = 0 ⇔ U_F γ_F = − γ_F U_F.
        T = ElectroweakFiniteTriple()
        G = T.chirality_F()
        U_F = T.real_structure_F()
        LHS = U_F @ G
        RHS = -G @ U_F
        err = float(np.max(np.abs(LHS - RHS)))
        assert err < 1e-14, (
            f"{{J_F, γ_F}}=0 (KO-6) residual: {err:.3e}; "
            f"expected anticommuting J,γ for KO-dim 6 of T_F."
        )

    def test_J_F_commutes_with_gamma_F_KO0_explicitly_fails(self):
        # Sanity: confirm that KO-dim 0 convention ([J_F, γ_F]=0) is FALSE
        # for the chosen convention. This documents the deviation from the
        # directive's "γ_F identical on both copies" framing.
        T = ElectroweakFiniteTriple()
        G = T.chirality_F()
        U_F = T.real_structure_F()
        LHS = U_F @ G
        RHS_plus = G @ U_F
        # Should NOT match
        assert not np.allclose(LHS, RHS_plus), (
            "Expected {J_F, γ_F}=0 (KO-6), not [J_F, γ_F]=0 (KO-0). "
            "If this test passes, the convention has drifted to KO-0."
        )

    def test_AC_triple_gamma_F_passthrough(self):
        # AlmostCommutativeTriple.gamma_F() should return the same object
        # (per call) as ElectroweakFiniteTriple.chirality_F().
        T_ac = minimal_electroweak_triple(2, yukawa_e=0.1, yukawa_nu=0.05)
        G_ac = T_ac.gamma_F()
        G_F = T_ac.finite.chirality_F()
        assert np.allclose(G_ac, G_F)
        assert G_ac.shape == (8, 8)

    def test_chirality_F_independent_of_yukawa(self):
        # γ_F depends only on the chirality grading of basis labels,
        # not on the value of the Yukawa.
        T0 = ElectroweakFiniteTriple()
        T1 = ElectroweakFiniteTriple(yukawa_e=0.5, yukawa_nu=0.3)
        T2 = ElectroweakFiniteTriple(yukawa_e=1.0 - 0.2j, yukawa_nu=0.0)
        assert np.allclose(T0.chirality_F(), T1.chirality_F())
        assert np.allclose(T1.chirality_F(), T2.chirality_F())


# ===========================================================================
# Section B. Algebra action
# ===========================================================================


class TestAlgebraAction:
    def test_quaternion_to_matrix_unit(self):
        Q = quaternion_to_matrix(1, 0, 0, 0)
        assert np.allclose(Q, SIGMA_0)

    def test_quaternion_to_matrix_i(self):
        Q = quaternion_to_matrix(0, 1, 0, 0)
        assert np.allclose(Q, 1j * SIGMA_1)

    def test_quaternion_to_matrix_unitary_for_real_unit(self):
        from math import sqrt
        a = 1 / sqrt(2)
        Q = quaternion_to_matrix(a, a, 0, 0)
        I = np.eye(2, dtype=np.complex128)
        assert np.allclose(Q @ Q.conj().T, I, atol=1e-14)

    def test_matter_action_identity(self):
        T = ElectroweakFiniteTriple()
        M = T.matter_action(1.0, (1.0, 0.0, 0.0, 0.0))
        assert np.allclose(M, np.eye(4))

    def test_algebra_action_block_structure(self):
        # pi_F acts on matter only; antimatter block is zero
        T = ElectroweakFiniteTriple()
        action = T.algebra_action(2.0 + 1j, (1.0, 0.5, 0.3, 0.1))
        assert action.shape == (8, 8)
        # Antimatter block is zero
        assert np.allclose(action[4:8, 4:8], 0.0)
        # Matter block matches matter_action
        M = T.matter_action(2.0 + 1j, (1.0, 0.5, 0.3, 0.1))
        assert np.allclose(action[0:4, 0:4], M)
        # Off-diagonal blocks are zero
        assert np.allclose(action[0:4, 4:8], 0.0)
        assert np.allclose(action[4:8, 0:4], 0.0)


# ===========================================================================
# Section C. Combined Dirac
# ===========================================================================


class TestCombinedDirac:
    def test_dim_H_combined(self):
        for n_max in [1, 2, 3]:
            T = minimal_electroweak_triple(n_max)
            assert T.dim_H == T.dim_GV * 8

    def test_combined_dirac_shape(self):
        T = minimal_electroweak_triple(2)
        D = T.dirac_combined()
        assert D.shape == (T.dim_H, T.dim_H)

    def test_combined_dirac_hermitian(self):
        T = minimal_electroweak_triple(2, yukawa_e=0.1, yukawa_nu=0.05)
        D = T.dirac_combined()
        assert np.allclose(D, D.conj().T, atol=1e-12)

    def test_combined_dirac_zero_yukawa_block_structure(self):
        T = gauge_only_triple(2)
        D = T.dirac_combined()
        d = T.dim_GV
        D_GV = T.D_GV()
        # D = D_GV (x) I_8 since D_F = 0
        D_reshape = D.reshape(d, 8, d, 8)
        for i in range(d):
            for j in range(d):
                block_ij = D_reshape[i, :, j, :]
                expected = D_GV[i, j] * np.eye(8, dtype=np.complex128)
                assert np.allclose(block_ij, expected, atol=1e-12)

    def test_combined_dirac_yukawa_contribution(self):
        T0 = gauge_only_triple(2)
        T1 = minimal_electroweak_triple(2, yukawa_e=0.1, yukawa_nu=0.0)
        D0 = T0.dirac_combined()
        D1 = T1.dirac_combined()
        delta = D1 - D0
        assert np.linalg.norm(delta) > 1e-6
        # Diagonal in GV index; chirality * D_F nonzero only at i = j
        d = T0.dim_GV
        delta_4d = delta.reshape(d, 8, d, 8).transpose(0, 2, 1, 3)
        for i in range(d):
            for j in range(d):
                if i != j:
                    assert np.allclose(delta_4d[i, j], 0.0, atol=1e-14)

    def test_gamma_GV_squared_identity(self):
        T = minimal_electroweak_triple(2)
        gamma = T.gamma_GV()
        I = np.eye(T.dim_GV, dtype=np.complex128)
        assert np.allclose(gamma @ gamma, I)


# ===========================================================================
# Section D. Real structure
# ===========================================================================


class TestCombinedRealStructure:
    def test_J_combined_unitary(self):
        T = minimal_electroweak_triple(2)
        U = T.real_structure_combined()
        I = np.eye(T.dim_H, dtype=np.complex128)
        assert np.allclose(U @ U.conj().T, I, atol=1e-12)

    def test_J_combined_squared_minus_I(self):
        # J^2 = J_GV^2 (x) J_F^2 = (-I)(+I) = -I
        # Combined KO-dim = 3 + 6 = 9 ≡ 1 (mod 8) -> J^2 = -I
        T = minimal_electroweak_triple(2)
        U = T.real_structure_combined()
        J2 = U @ np.conj(U)
        I = np.eye(T.dim_H, dtype=np.complex128)
        assert np.allclose(J2, -I, atol=1e-12)

    def test_J_combined_KO1_sign(self):
        """Combined KO-dim 1 should give (epsilon, eps')=(-, +)."""
        T = minimal_electroweak_triple(2, yukawa_e=0.1, yukawa_nu=0.05)
        D = T.dirac_combined()
        U = T.real_structure_combined()
        # JD as antilinear: U conj(D) = sign * D U
        LHS = U @ np.conj(D)
        RHS_plus = D @ U
        RHS_minus = -D @ U
        err_plus = float(np.max(np.abs(LHS - RHS_plus)))
        err_minus = float(np.max(np.abs(LHS - RHS_minus)))
        # We expect (+) sign to fit (KO-1)
        assert err_plus < err_minus, (
            f"err_plus={err_plus:.3e}, err_minus={err_minus:.3e}; "
            f"expected JD = +DJ for KO-1"
        )
        # And it should hold to high precision
        assert err_plus < 1e-12, f"JD = +DJ residual: {err_plus:.3e}"


# ===========================================================================
# Section E. Inner fluctuations
# ===========================================================================


class TestInnerFluctuations:
    def _identity_quaternion(self):
        return (1.0, 0.0, 0.0, 0.0)

    def test_inner_fluctuation_zero_when_a_b_identity(self):
        T = minimal_electroweak_triple(2)
        I_GV = np.eye(T.dim_GV, dtype=np.complex128)
        unit_q = self._identity_quaternion()
        gens = [(I_GV, 1.0, unit_q, I_GV, 1.0, unit_q)]
        omega = T.inner_fluctuation_one_form(gens)
        assert np.linalg.norm(omega) < 1e-10

    def test_decompose_higgs_zero_yukawa(self):
        # With D_F = 0, inner fluctuations have NO Higgs content.
        T = gauge_only_triple(2)
        I_GV = np.eye(T.dim_GV, dtype=np.complex128)
        M = T.gv_multiplier(1)
        gens = [(I_GV, 0.7 + 0.3j, (0.4, 0.2, 0.1, 0.05),
                 M, 0.5, (0.6, 0.3, 0.2, 0.1))]
        omega = T.inner_fluctuation_one_form(gens)
        higgs_norm = T.higgs_norm(omega)
        assert higgs_norm < 1e-12, f"Higgs sector should be zero, got {higgs_norm}"

    def test_decompose_higgs_nonzero_with_yukawa(self):
        # Non-zero Yukawa => Higgs sector activates.
        T = minimal_electroweak_triple(2, yukawa_e=0.5, yukawa_nu=0.3)
        I_GV = np.eye(T.dim_GV, dtype=np.complex128)
        # Use generators that activate D_F off-diagonal:
        # b in C-summand: pi_F^matter(0.5, 0_q) = block_diag(0_2, diag(0.5, 0.5)).
        # [D_F, b_F] is nonzero off-diag block.
        gens = [(I_GV, 1.0, (1.0, 0.0, 0.0, 0.0),
                 I_GV, 0.5, (0.0, 0.0, 0.0, 0.0))]
        omega = T.inner_fluctuation_one_form(gens)
        higgs_norm = T.higgs_norm(omega)
        assert higgs_norm > 1e-6, f"Expected non-zero Higgs; got {higgs_norm}"

    def test_matter_antimatter_off_block_zero(self):
        """Inner fluctuations should NOT mix matter and antimatter sectors.

        Both pi_F (algebra) and D_F (Dirac) are block-diagonal in
        matter/antimatter; [D, b] is therefore also block-diagonal.
        The matter<->antimatter off-block of omega should be zero.
        """
        T = minimal_electroweak_triple(2, yukawa_e=0.5)
        I_GV = np.eye(T.dim_GV, dtype=np.complex128)
        gens = [(I_GV, 1.0, (1.0, 0.0, 0.0, 0.0),
                 I_GV, 0.5, (0.0, 0.0, 0.0, 0.0))]
        omega = T.inner_fluctuation_one_form(gens)
        ma_norm = T.matter_antimatter_off_norm(omega)
        assert ma_norm < 1e-12, f"matter<->antimatter off-block: {ma_norm}"

    def test_fluctuated_dirac_zero_omega(self):
        T = minimal_electroweak_triple(2, yukawa_e=0.1)
        omega = np.zeros((T.dim_H, T.dim_H), dtype=np.complex128)
        D_omega = T.fluctuated_dirac(omega, epsilon_prime=+1)
        D = T.dirac_combined()
        assert np.allclose(D_omega, D, atol=1e-12)


# ===========================================================================
# Section F. Falsifier check (memo §5)
# ===========================================================================


class TestFalsifier:
    def test_negative_holds_when_yukawa_zero(self):
        # With D_F = 0, no Higgs sector can appear from inner fluctuations.
        T = gauge_only_triple(2)
        negative_holds, reason, data = T.check_natural_negative(
            n_random_generators=20, seed=42,
        )
        assert negative_holds, f"Expected negative to hold; got {reason}"
        assert data["higgs_max"] < 1e-10

    def test_negative_fails_with_realistic_yukawa(self):
        T = minimal_electroweak_triple(2, yukawa_e=0.3, yukawa_nu=0.1)
        negative_holds, reason, data = T.check_natural_negative(
            n_random_generators=20, seed=42,
        )
        assert not negative_holds, f"Expected negative to fail; got {reason}"
        assert data["higgs_max"] > 1e-6

    def test_negative_holds_with_zero_DF_n_max_3(self):
        T = gauge_only_triple(3)
        negative, _, data = T.check_natural_negative(20, seed=0)
        assert negative


# ===========================================================================
# Section G. Spectral action / Yukawa scan
# ===========================================================================


class TestSpectralAction:
    def test_eigenvalue_count(self):
        T = minimal_electroweak_triple(2, yukawa_e=0.1)
        D = T.dirac_combined()
        eigs = spectral_action_eigenvalues(D)
        assert len(eigs) == T.dim_H

    def test_heat_kernel_trace_positivity(self):
        T = minimal_electroweak_triple(2, yukawa_e=0.1)
        D = T.dirac_combined()
        for t in [0.1, 1.0, 10.0]:
            tr = heat_kernel_trace(D, t)
            assert tr > 0

    def test_yukawa_dependence_zero_baseline(self):
        result = yukawa_dependence_scan(2, [0.0, 0.05, 0.1, 0.2])
        assert result["Frob_shift"][0] == 0.0
        assert all(s2 > s1 for s1, s2 in zip(result["Frob_shift"], result["Frob_shift"][1:]))

    def test_Tr_D_squared_polynomial_in_yukawa(self):
        """Tr(D^2) = Tr(D_0^2) + 4 * dim_GV * (|y_nu|^2 + |y_e|^2).

        Cross-term D_GV gamma (x) D_F integrates to 0 since Tr(D_F) = 0
        (D_F is off-diag in L<->R).

        D_F^2 has eigenvalues |y_nu|^2 (x4: matter±, antimatter±) and
        |y_e|^2 (x4), so Tr(D_F^2) = 4(|y_nu|^2 + |y_e|^2).
        Then Tr((gamma_GV)^2 (x) D_F^2) = dim_GV * Tr(D_F^2)
        = dim_GV * 4 * (|y_nu|^2 + |y_e|^2).
        """
        ys = [0.0, 0.1, 0.2, 0.3]
        result = yukawa_dependence_scan(2, ys)
        Tr_D2 = result["Tr_D_squared"]
        T = minimal_electroweak_triple(2)
        # Slope: dim_GV * 4 (since Tr(D_F^2) = 4 y_e^2 with y_nu = 0 in the scan)
        slope_predicted = 4 * T.dim_GV
        for y, val in zip(ys, Tr_D2):
            predicted = Tr_D2[0] + slope_predicted * y**2
            assert abs(val - predicted) < 1e-10


# ===========================================================================
# Section H. Cross-validation against memo §5 falsifier statement
# ===========================================================================


class TestMemoFalsifier:
    """Memo §5: falsifier statement is 'show that for every Hermitian D_F
    derivable from GeoVac structure, inner fluctuation produces only gauge.'

    H1 architecture explicitly switched away from offdiag CH 'Candidate A'
    (Track 2 ruled it out). We test the IMPOSED D_F. Verdict:

    - Falsifier HOLDS for D_F = 0 (no Higgs).
    - Falsifier FAILS for non-zero imposed D_F (Higgs appears).

    GeoVac does not select either; both candidates A and B failed.
    Result: Higgs CAN be defined on the AC extension, but it is IMPOSED,
    not selected from GeoVac structure.
    """

    def test_negative_with_zero_DF_holds_n_max_1(self):
        T = gauge_only_triple(1)
        negative, _, data = T.check_natural_negative(20, seed=0)
        assert negative
        assert data["higgs_max"] < 1e-12

    def test_negative_with_zero_DF_holds_n_max_3(self):
        T = gauge_only_triple(3)
        negative, _, data = T.check_natural_negative(20, seed=0)
        assert negative

    def test_negative_with_realistic_yukawa_fails(self):
        T = minimal_electroweak_triple(2, yukawa_e=0.01, yukawa_nu=0.001)
        negative, _, data = T.check_natural_negative(40, seed=0)
        assert not negative
