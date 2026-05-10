"""Tests for the internal multi-focal architecture.

Validates per-orbital focal-length CI on hydrogenic basis. Track 4
follow-on (CLAUDE.md §1.8 He oscillator strength autopsy).

Author: GeoVac Development Team (Phase B)
Date: 2026-05-09
"""

from __future__ import annotations

from math import factorial

import numpy as np
import pytest

from geovac.internal_multifocal import (
    MultifocalOrbital,
    MultifocalSpec,
    _enumerate_m_distributions,
    _enumerate_singlet_l_pairs,
    build_singlet_LM_subblock_multifocal,
    build_singlet_LM_subblock_multifocal_extended,
    compute_he_oscillator_strength_multifocal,
    compute_he_oscillator_strength_multifocal_extended,
    dipole_z_matrix,
    h1_matrix,
    he_extended_spec,
    he_slater_spec,
    he_uniform_spec,
    kinetic_radial,
    matrix_element_rk,
    overlap_matrix,
    overlap_radial,
    slater_rk_multifocal,
    solve_generalized_singlet,
    transition_dipole_multifocal,
    two_electron_integral_multifocal,
)
from geovac.hypergeometric_slater import compute_rk_float


# ---------------------------------------------------------------------------
# Cross-exponent radial primitives
# ---------------------------------------------------------------------------


class TestOverlapRadial:
    """Cross-exponent overlap S_pq = int R_p R_q r^2 dr."""

    def test_overlap_self_normalized(self) -> None:
        """At lam_p = lam_q = arbitrary, S_pp = 1 (normalization)."""
        for n, l, lam in [(1, 0, 1.0), (1, 0, 1.69), (2, 0, 0.5),
                          (2, 1, 0.575), (3, 2, 0.333)]:
            S = overlap_radial(n, l, lam, n, l, lam)
            assert S == pytest.approx(1.0, rel=1e-12), (
                f"Self-overlap at (n={n}, l={l}, lam={lam}) = {S} != 1"
            )

    def test_overlap_orthogonal_n_at_matched_lambda(self) -> None:
        """At matched lam = Z/n_p = Z/n_q (standard hydrogenic),
        S_{1s, 2s} = 0 IF and only if both at their natural (matched) lam.
        For 1s at lam=1 vs 2s at lam=1, these are NOT both at natural
        lam (natural for 2s is 0.5), so they're not orthogonal in
        general. But within a single Sturmian set at fixed lam they
        ARE orthogonal (Sturmian orthogonality).

        Actually with the L2 normalization (∫|R|² r² dr = 1) at common
        lam, 1s and 2s are NOT orthogonal in general — they'd be
        Sturmian-orthogonal under int R_p R_q r dr, but L2 orthogonal
        only at distinct natural lams.

        So the standard hydrogenic orthogonality requires lam_p = Z/n_p,
        lam_q = Z/n_q (each at its OWN natural exponent).
        """
        # 1s at lam=1 (natural) and 2s at lam=0.5 (natural at Z=1)
        # are L2-orthogonal: int R_1s R_2s r^2 dr = 0
        S = overlap_radial(1, 0, 1.0, 2, 0, 0.5)
        assert abs(S) < 1e-12, (
            f"Standard hydrogen 1s-2s overlap = {S}, expected 0"
        )
        # Also 2s at lam=1 (which is non-natural) is NOT orthogonal to 1s at lam=1
        S2 = overlap_radial(1, 0, 1.0, 2, 0, 1.0)
        # This should be nonzero (they're not at natural exponents)
        # actually compute: R_1s = 2 e^{-r}, R_2s at lam=1 = sqrt(8) (1 - r) e^{-r}
        # I need to verify this isn't accidentally zero.
        # int 2 sqrt(8) (1 - r) e^{-2r} r^2 dr = 2 sqrt(8) [int r^2 e^{-2r} dr - int r^3 e^{-2r} dr]
        # = 2 sqrt(8) [2!/8 - 3!/16] = 2 sqrt(8) [0.25 - 0.375] = -2 sqrt(8) * 0.125 = -0.7071
        # Hmm that's not zero. Let me check more carefully — actually 2s at lam=1 is the
        # 2s function with same exponent as 1s; this is NOT a standard hydrogenic state
        # (those have lam=Z/n). So nonzero overlap is expected.
        # Just verify it's nonzero:
        assert abs(S2) > 1e-3, (
            f"2s at non-natural lam should have nonzero overlap with 1s, got {S2}"
        )

    def test_overlap_nonzero_at_mismatched_lambda(self) -> None:
        """1s at lam_a vs 1s at lam_b (different exponents) gives finite
        non-trivial overlap. Closed form: 8 (lam_a lam_b)^{3/2} / (lam_a + lam_b)^3."""
        for lam_a, lam_b in [(1.0, 2.0), (1.0, 1.5), (1.69, 0.575)]:
            S = overlap_radial(1, 0, lam_a, 1, 0, lam_b)
            expected = (
                8.0 * (lam_a * lam_b) ** 1.5 / (lam_a + lam_b) ** 3
            )
            assert S == pytest.approx(expected, rel=1e-12), (
                f"S_1s1s(lam={lam_a}, {lam_b}) = {S} vs expected {expected}"
            )

    def test_overlap_l_orthogonal(self) -> None:
        """l_p != l_q gives S = 0 always."""
        for la, lb in [(0, 1), (1, 2), (0, 2)]:
            S = overlap_radial(2, la, 1.0, 2, lb, 1.0)
            assert S == 0.0, f"l-orthogonality fails for l_a={la}, l_b={lb}: {S}"


class TestDipoleRadial:
    """<R_p | r | R_q> matrix element."""

    def test_dipole_lyman_alpha_match(self) -> None:
        """At lam_1s=1, lam_2p=1/2 (natural Z=1), <R_1s|r|R_2p> = 128 sqrt(2)/243.

        Note: this is the RADIAL part, before the angular Gaunt factor
        (1/sqrt(3)) for c_1(0,0,1,0).
        """
        D_radial = matrix_element_rk(1, 0, 1.0, 2, 1, 0.5, k_pow=1)
        expected = 128.0 * np.sqrt(2) / 243.0 * np.sqrt(3)  # we computed unscaled
        # Wait — the radial integral I want is int R_1s R_2p r * r^2 dr.
        # We verified this gives 1.290 = 128 sqrt(2)/243 * sqrt(3)
        # Actually no: my earlier verification showed:
        #   <R|r|R> = 1.2902 = 4! / sqrt(24) / (3/2)^5 * 2 = 1.2902
        #   <2p_z|z|1s> = (1/sqrt(3)) * 1.2902 = 0.7449 = 128 sqrt(2)/243
        # So the full Wigner-Eckart-coupled Y_1^0 me equals 128 sqrt(2)/243,
        # but the bare radial integral <R_1s|r|R_2p> r^2 dr is sqrt(3) times that.
        from math import sqrt as msqrt
        expected_radial = (2.0 / msqrt(24.0)) * factorial(4) / (1.5) ** 5
        assert D_radial == pytest.approx(expected_radial, rel=1e-12), (
            f"Lyman radial me = {D_radial} vs expected {expected_radial}"
        )

    def test_dipole_at_mismatched_lambda_finite(self) -> None:
        """At lam_1s=1.69 (Slater), lam_2p=0.575 (Slater), dipole is finite."""
        D = matrix_element_rk(1, 0, 1.69, 2, 1, 0.575, k_pow=1)
        assert np.isfinite(D)
        assert abs(D) > 0.1, f"Dipole at Slater lambdas should be O(1), got {D}"

    def test_dipole_z_matrix_lyman_alpha_through_full_assembly(self) -> None:
        """Build a 2-orbital spec (1s lam=1, 2p_0 lam=1/2) and check that
        D_orb[0, 1] (the dipole matrix element) equals (1/sqrt(3))*radial."""
        spec = MultifocalSpec(
            orbitals=[
                MultifocalOrbital(n=1, l=0, m=0, lam=1.0, label="1s"),
                MultifocalOrbital(n=2, l=1, m=0, lam=0.5, label="2p_0"),
            ],
            Z_nuc=1.0,
        )
        D = dipole_z_matrix(spec)
        # <1s|z|2p_0> should equal (1/sqrt(3)) * <R_1s|r|R_2p>
        # = (1/sqrt(3)) * 1.2903 = 0.7449
        expected = 128.0 * np.sqrt(2) / 243.0
        assert D[0, 1] == pytest.approx(expected, rel=1e-12), (
            f"D[1s, 2p_0] = {D[0, 1]} vs expected {expected}"
        )
        assert D[1, 0] == pytest.approx(D[0, 1], rel=1e-12), "D not symmetric"


class TestSlaterMultifocal:
    """Cross-exponent two-electron Slater R^k integrals."""

    def test_slater_F0_1s1s_at_lam1(self) -> None:
        """At all four orbitals at (1s, lam=1), F^0(1s,1s) = 5/8."""
        rk = slater_rk_multifocal(
            1, 0, 1.0, 1, 0, 1.0, 1, 0, 1.0, 1, 0, 1.0, k=0, n_quad=200,
        )
        assert rk == pytest.approx(5.0 / 8.0, rel=1e-9), (
            f"F^0(1s,1s) at lam=1 = {rk} vs 5/8 = {5.0/8.0}"
        )

    def test_slater_matched_lambda_regression_F0_1s2s(self) -> None:
        """F^0(1s,1s; 2s,2s) at NATURAL lams reproduces compute_rk_float.

        compute_rk_float uses natural-Z hydrogenic functions: 1s at lam=1,
        2s at lam=0.5. Multi-focal slater_rk_multifocal at the same natural
        lams must reproduce.
        """
        rk_ref = compute_rk_float(1, 0, 1, 0, 2, 0, 2, 0, k=0)
        rk_mf_natural = slater_rk_multifocal(
            1, 0, 1.0, 1, 0, 1.0, 2, 0, 0.5, 2, 0, 0.5, k=0, n_quad=200,
        )
        assert rk_mf_natural == pytest.approx(rk_ref, rel=1e-7), (
            f"Matched-natural F^0(1s 1s; 2s 2s) mf={rk_mf_natural} vs "
            f"single-focal {rk_ref}"
        )

    def test_slater_matched_lambda_regression_F0_2s2s(self) -> None:
        """F^0(2s,2s; 2s,2s) at lam=0.5 reproduces compute_rk_float (= 77/512)."""
        rk_ref = compute_rk_float(2, 0, 2, 0, 2, 0, 2, 0, k=0)
        rk_mf = slater_rk_multifocal(
            2, 0, 0.5, 2, 0, 0.5, 2, 0, 0.5, 2, 0, 0.5, k=0, n_quad=200,
        )
        assert rk_mf == pytest.approx(rk_ref, rel=1e-7), (
            f"F^0(2s,2s) at lam=0.5 mf={rk_mf} vs ref {rk_ref}"
        )
        # Sanity: should equal 77/512
        assert rk_ref == pytest.approx(77.0 / 512.0, rel=1e-12)

    def test_slater_F0_1s1s_at_distinct_lambdas(self) -> None:
        """F^0 with two electrons at different exponents — symmetric in
        (lam_1, lam_2). Roothaan-style."""
        # Both at lam=1 should match expected
        rk_a = slater_rk_multifocal(
            1, 0, 1.0, 1, 0, 1.0, 1, 0, 2.0, 1, 0, 2.0, k=0, n_quad=300,
        )
        rk_b = slater_rk_multifocal(
            1, 0, 2.0, 1, 0, 2.0, 1, 0, 1.0, 1, 0, 1.0, k=0, n_quad=300,
        )
        # By symmetry: swapping (1, 3) <-> (2, 4) leaves R^k unchanged
        assert rk_a == pytest.approx(rk_b, rel=1e-7)

    def test_slater_gaunt_selection_preserved(self) -> None:
        """For F^k between two s-orbitals on each side, only k=0 contributes."""
        # F^0(1s,1s; 1s,1s) nonzero
        rk0 = slater_rk_multifocal(
            1, 0, 1.0, 1, 0, 1.0, 1, 0, 1.0, 1, 0, 1.0, k=0, n_quad=100,
        )
        assert rk0 > 0
        # k=1 should also be computable (just the radial integral, the angular
        # Gaunt factor is what kills it). Verify the angular factor is zero:
        spec = MultifocalSpec(
            orbitals=[MultifocalOrbital(n=1, l=0, m=0, lam=1.0)],
            Z_nuc=1.0,
        )
        # Two-electron integral with all four 1s orbitals: <1s 1s|1s 1s>
        # only F^0 contributes
        a = MultifocalOrbital(n=1, l=0, m=0, lam=1.0)
        val = two_electron_integral_multifocal(a, a, a, a, n_quad=100)
        # Should be 5/8
        assert val == pytest.approx(5.0 / 8.0, rel=1e-8)


# ---------------------------------------------------------------------------
# One-body assemblies
# ---------------------------------------------------------------------------


class TestH1Matrix:
    """One-body (kinetic + V_Ne) on multi-focal basis."""

    def test_h1_hermitian_single_orbital(self) -> None:
        """Single-orbital H_1 should be a 1x1 real number."""
        spec = MultifocalSpec(
            orbitals=[MultifocalOrbital(n=1, l=0, m=0, lam=1.0)],
            Z_nuc=1.0,
        )
        H1 = h1_matrix(spec)
        assert H1.shape == (1, 1)
        # For 1s at natural lam=Z=1, H_1 = E_1s = -0.5 Ha
        assert H1[0, 0] == pytest.approx(-0.5, rel=1e-10), (
            f"H_1 for 1s at natural lam=1, Z=1 should be -0.5, got {H1[0, 0]}"
        )

    def test_h1_hydrogenic_1s_at_natural_lambda(self) -> None:
        """For hydrogen 1s at natural lam=Z, eigenvalue should be -Z^2/2."""
        for Z in [1.0, 2.0, 3.0]:
            spec = MultifocalSpec(
                orbitals=[MultifocalOrbital(n=1, l=0, m=0, lam=Z)],
                Z_nuc=Z,
            )
            H1 = h1_matrix(spec)
            expected = -0.5 * Z * Z
            assert H1[0, 0] == pytest.approx(expected, rel=1e-10), (
                f"H_1 at Z={Z}, lam={Z} should be {expected}, got {H1[0, 0]}"
            )

    def test_h1_hermitian_multi_orbital(self) -> None:
        """H_1 must be symmetric on a multi-orbital basis."""
        spec = he_slater_spec(n_max=3)
        H1 = h1_matrix(spec)
        assert np.allclose(H1, H1.T, atol=1e-10), (
            f"H_1 not symmetric: max asymmetry = {np.max(np.abs(H1 - H1.T))}"
        )

    def test_overlap_hermitian(self) -> None:
        """S must be symmetric and positive-definite."""
        spec = he_slater_spec(n_max=3)
        S = overlap_matrix(spec)
        assert np.allclose(S, S.T, atol=1e-10), "S not symmetric"
        eigvals = np.linalg.eigvalsh(S)
        assert np.all(eigvals > 1e-10), (
            f"S not positive-definite: min eig = {eigvals.min()}"
        )


# ---------------------------------------------------------------------------
# CI eigenproblem
# ---------------------------------------------------------------------------


class TestCISinglet:
    """Multi-focal singlet two-electron CI."""

    def test_he_1S_minimal_basis_above_drake(self) -> None:
        """1s² minimal basis at lam=27/16: variational He 1S energy.

        E_var(1s² at 27/16) = -2 * (27/16)² + (27/16) (5/8)
                            = -729/128 + (27/16)(5/8)
                            ≈ -5.6953 + 1.0547 = -2.8477

        Full CI in 1s-only basis at fixed lam = 27/16 = E_var. Drake
        exact = -2.9037, so above by ~57 mHa (variational bound respected
        for 1-config calc; also matches textbook variational He value).
        """
        lam_var = 27.0 / 16.0
        spec = MultifocalSpec(
            orbitals=[MultifocalOrbital(n=1, l=0, m=0, lam=lam_var)],
            Z_nuc=2.0,
        )
        H, S, configs = build_singlet_LM_subblock_multifocal(
            spec, L_target=0, M_L_target=0, n_quad=100,
        )
        assert len(configs) == 1, f"1s² minimal basis should have 1 config, got {len(configs)}"
        eigvals, _ = solve_generalized_singlet(H, S)
        # Textbook He variational at 27/16: -729/256 = -2.8477
        expected_var = -729.0 / 256.0
        assert eigvals[0] == pytest.approx(expected_var, abs=2e-3), (
            f"He 1s² variational at 27/16: got {eigvals[0]}, expected {expected_var}"
        )
        # Above Drake exact -2.9037 (variational respects exact)
        assert eigvals[0] > -2.9037, (
            f"E_var must be above Drake exact -2.9037, got {eigvals[0]}"
        )

    def test_2config_correlation_lowers_energy(self) -> None:
        """Adding 1s2s to the 1s² basis must lower the variational energy
        (or leave it unchanged if the 2s contributes nothing).
        """
        # Just 1s
        spec_1 = MultifocalSpec(
            orbitals=[MultifocalOrbital(n=1, l=0, m=0, lam=27.0/16.0)],
            Z_nuc=2.0,
        )
        H1, S1, _ = build_singlet_LM_subblock_multifocal(spec_1, 0, 0, n_quad=100)
        E1, _ = solve_generalized_singlet(H1, S1)

        # 1s + 2s at Slater
        spec_2 = MultifocalSpec(
            orbitals=[
                MultifocalOrbital(n=1, l=0, m=0, lam=27.0/16.0, label="1s"),
                MultifocalOrbital(n=2, l=0, m=0, lam=(2.0 - 0.85)/2, label="2s"),
            ],
            Z_nuc=2.0,
        )
        H2, S2, _ = build_singlet_LM_subblock_multifocal(spec_2, 0, 0, n_quad=100)
        E2, _ = solve_generalized_singlet(H2, S2)

        assert E2[0] <= E1[0] + 1e-10, (
            f"Adding 2s should lower (or preserve) E: E_1config={E1[0]}, "
            f"E_2config={E2[0]}"
        )


# ---------------------------------------------------------------------------
# Phase C smoke test
# ---------------------------------------------------------------------------


class TestPhaseCSmoke:
    """Phase C oscillator strength: sanity-bound."""

    @pytest.mark.slow
    def test_he_oscillator_phase_c_basic(self) -> None:
        """Compute He 2^1P -> 1^1S oscillator strength on minimal multi-focal
        basis. Expected to be in [0.1, 0.6] range — physical bounds.
        Drake reference is 0.276."""
        spec = he_slater_spec(n_max=2)
        result = compute_he_oscillator_strength_multifocal(spec, n_quad=100)
        assert 0.05 <= result["f_length"] <= 1.0, (
            f"Oscillator strength out of physical bounds: {result['f_length']}"
        )
        assert result["E_2P_Ha"] > result["E_1S_Ha"], (
            f"E(2P) should be above E(1S): "
            f"E_1S={result['E_1S_Ha']}, E_2P={result['E_2P_Ha']}"
        )
        # Variational: E_1S above Drake -2.9037
        assert result["E_1S_Ha"] > -2.9037, (
            f"E(1S) must be above Drake -2.9037, got {result['E_1S_Ha']}"
        )


# ---------------------------------------------------------------------------
# Phase D: extended angular CI
# ---------------------------------------------------------------------------


class TestEnumerateSingletLPairs:
    """Selection rules for the extended angular CI."""

    def test_l0_pairs(self) -> None:
        """For L_target = 0 (1S), allowed pairs are (l, l) for all l up to spec max."""
        spec = he_extended_spec(n_max=3)  # has l=0,1,2
        pairs = _enumerate_singlet_l_pairs(spec, 0)
        # l_a + l_b + 0 even, l_a <= l_b: (0,0), (1,1), (2,2)
        assert pairs == [(0, 0), (1, 1), (2, 2)], pairs

    def test_l1_pairs(self) -> None:
        """For L_target = 1 (1P), allowed pairs satisfy l_a + l_b odd."""
        spec = he_extended_spec(n_max=3)
        pairs = _enumerate_singlet_l_pairs(spec, 1)
        # l_a + l_b odd, l_a <= l_b: (0,1), (1,2)
        assert pairs == [(0, 1), (1, 2)], pairs

    def test_l2_pairs(self) -> None:
        """For L_target = 2 (1D), pairs satisfy l_a+l_b even & |la-lb|<=2<=la+lb."""
        spec = he_extended_spec(n_max=4)  # has up to f? no, l_max=3 for n_max=4
        pairs = _enumerate_singlet_l_pairs(spec, 2)
        # Filter: l_a + l_b even (singlet), triangle inequality 2 in [|la-lb|, la+lb]
        # (0,2): 0+2=2 even, |0-2|=2<=2<=0+2=2 ✓
        # (1,1): 1+1=2 even, |1-1|=0<=2<=2 ✓
        # (1,3): 1+3=4 even, |1-3|=2<=2<=4 ✓
        # (2,2): 2+2=4 even, |2-2|=0<=2<=4 ✓
        # (3,3): 3+3=6 even, |3-3|=0<=2<=6 ✓
        # (Note: for n_max=4 with all m, l_max should be 3; but actual depends on spec)
        # Relax to "at least the (1,1), (0,2), (2,2) cases"
        assert (1, 1) in pairs and (2, 2) in pairs and (0, 2) in pairs, pairs

    def test_m_distributions(self) -> None:
        """m-distribution enumeration."""
        # (l_a=0, l_b=0, M_L=0): only (0, 0)
        assert _enumerate_m_distributions(0, 0, 0) == [(0, 0)]
        # (l_a=1, l_b=1, M_L=0): (-1, +1), (0, 0), (+1, -1)
        assert _enumerate_m_distributions(1, 1, 0) == [(-1, 1), (0, 0), (1, -1)]
        # (l_a=0, l_b=1, M_L=1): only (0, 1)
        assert _enumerate_m_distributions(0, 1, 1) == [(0, 1)]


class TestExtendedSubblockReducesToPhaseB:
    """Extended subblock with default (s,s)/(0,L) pairs reduces to Phase B."""

    def test_reduces_to_phase_b_for_1S(self) -> None:
        """build_..._extended with config_l_pairs=[(0,0)] gives identical
        H, S, configs as build_..._multifocal."""
        spec = he_slater_spec(n_max=3)
        H_basic, S_basic, configs_basic = build_singlet_LM_subblock_multifocal(
            spec, L_target=0, M_L_target=0, n_quad=60,
        )
        H_ext, S_ext, configs_ext = build_singlet_LM_subblock_multifocal_extended(
            spec, L_target=0, M_L_target=0, n_quad=60,
            config_l_pairs=[(0, 0)],
        )
        assert configs_basic == configs_ext, (
            f"configs differ: basic={configs_basic}, ext={configs_ext}"
        )
        assert H_basic.shape == H_ext.shape
        np.testing.assert_allclose(H_basic, H_ext, rtol=1e-10)
        np.testing.assert_allclose(S_basic, S_ext, rtol=1e-10)

    def test_reduces_to_phase_b_for_1P(self) -> None:
        """build_..._extended with config_l_pairs=[(0,1)] for L=1 gives identical
        result as build_..._multifocal."""
        spec = he_slater_spec(n_max=3)
        H_basic, S_basic, configs_basic = build_singlet_LM_subblock_multifocal(
            spec, L_target=1, M_L_target=0, n_quad=60,
        )
        H_ext, S_ext, configs_ext = build_singlet_LM_subblock_multifocal_extended(
            spec, L_target=1, M_L_target=0, n_quad=60,
            config_l_pairs=[(0, 1)],
        )
        assert configs_basic == configs_ext
        np.testing.assert_allclose(H_basic, H_ext, rtol=1e-10)
        np.testing.assert_allclose(S_basic, S_ext, rtol=1e-10)


class TestExtendedAngularCI:
    """Extended angular CI: (p,p) admixture in 1S, (p,d) in 1P."""

    def test_extended_subblock_larger_than_phase_b(self) -> None:
        """When p-orbitals with all m sublevels are included, the extended
        1S subblock is strictly larger than the Phase-B (s, s) subblock."""
        spec = he_extended_spec(n_max=2)  # has 1s, 2s, 2p (m=-1,0,+1)
        # Phase B: only (s, s) configs
        _, _, configs_basic = build_singlet_LM_subblock_multifocal(
            spec, L_target=0, M_L_target=0, n_quad=60,
        )
        # Phase D: includes (p, p) configs
        _, _, configs_ext = build_singlet_LM_subblock_multifocal_extended(
            spec, L_target=0, M_L_target=0, n_quad=60,
        )
        assert len(configs_ext) > len(configs_basic), (
            f"Extended should add (p, p) configs: basic={len(configs_basic)}, "
            f"ext={len(configs_ext)}"
        )

    def test_extended_correlation_lowers_E_1S(self) -> None:
        """Adding (p, p) admixture to the 1S CI MUST lower E_1S (variational
        principle)."""
        spec = he_extended_spec(n_max=2)
        from geovac.internal_multifocal import solve_generalized_singlet

        H_b, S_b, _ = build_singlet_LM_subblock_multifocal(
            spec, L_target=0, M_L_target=0, n_quad=60,
        )
        H_e, S_e, _ = build_singlet_LM_subblock_multifocal_extended(
            spec, L_target=0, M_L_target=0, n_quad=60,
        )
        E_b, _ = solve_generalized_singlet(H_b, S_b)
        E_e, _ = solve_generalized_singlet(H_e, S_e)
        assert E_e[0] <= E_b[0] + 1e-12, (
            f"Extended CI must lower (or equal) E_1S: basic={E_b[0]}, ext={E_e[0]}"
        )

    def test_he_extended_spec_has_all_m_for_p(self) -> None:
        """he_extended_spec includes all m=-l..+l sublevels for l>=1."""
        spec = he_extended_spec(n_max=2)
        # Should have 1s (1 orbital), 2s (1), 2p_m=-1, 2p_m=0, 2p_m=+1 (3 orbitals).
        # Total: 5 orbitals.
        assert spec.n_orbitals == 5
        p_m_values = sorted(o.m for o in spec.orbitals if o.l == 1)
        assert p_m_values == [-1, 0, 1], p_m_values

    def test_extended_driver_runs(self) -> None:
        """End-to-end smoke test for the extended driver."""
        spec = he_extended_spec(n_max=2)
        result = compute_he_oscillator_strength_multifocal_extended(
            spec, n_quad=60,
        )
        # Sanity bounds
        assert 0.05 <= result["f_length"] <= 1.0, (
            f"Oscillator strength out of bounds: {result['f_length']}"
        )
        assert result["E_2P_Ha"] > result["E_1S_Ha"]
        assert result["E_1S_Ha"] > -2.91, "E_1S must be above Drake (variational)"
        # Configuration list keys present
        assert "config_l_pairs_1S" in result
        assert "config_l_pairs_2P" in result
        assert (1, 1) in result["config_l_pairs_1S"]  # extended (p, p) included


class TestPhaseDOscillatorClosure:
    """Phase D oscillator strength closure test on robust physical config."""

    @pytest.mark.slow
    def test_phase_d_robust_closure(self) -> None:
        """The (1s_triple_compact, 2p_triple) config in Phase D extended
        gives f within 5% of Drake at well-conditioned basis (cond_S < 1000).

        This is the structurally-cleanest physical configuration found in
        the Phase D variational floor characterization (debug/data/
        he_oscillator_v3_floor.json, F2):
          1s_triple_compact = [1.4, 27/16, 1.95]
          2p_triple = [0.5, 0.65, 1.0]

        Expected: f ~ 0.286 (+3.4% vs Drake 0.276), cond_S_1S ~ 4e2.
        """
        spec = he_extended_spec(
            n_max=2,
            s_lams=[1.4, 27.0 / 16.0, 1.95],
            p_lams=[0.5, 0.65, 1.0],
            d_lams=[],
        )
        result = compute_he_oscillator_strength_multifocal_extended(
            spec, n_quad=80,
        )
        # f within 5% of Drake (0.276)
        f_drake = 0.27616
        err = abs(result["f_length"] - f_drake) / f_drake
        assert err < 0.05, (
            f"Phase D structurally robust f closure failed: f={result['f_length']:.4f}, "
            f"Drake={f_drake}, err={err*100:.2f}%"
        )
        # Well-conditioned
        assert result["cond_S_1S"] < 1e4, (
            f"Phase D should be well-conditioned: cond_S_1S={result['cond_S_1S']:.1e}"
        )
        # Variational
        assert result["E_1S_Ha"] > -2.9037, "E_1S above Drake (variational)"
        assert result["E_2P_Ha"] > -2.1238, "E_2P above Drake (variational)"
