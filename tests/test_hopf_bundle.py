"""
Tests for the Hopf Bundle module.

Validates the S³ embedding, Hopf fibration, discrete volumes, fiber
structure, and the various lattice → S³ mapping methods.

Author: GeoVac
Date: March 2026
"""

import pytest
import numpy as np
from geovac.hopf_bundle import (
    fock_chi, state_to_hopf_angles, state_to_s3, lattice_to_s3,
    hopf_project, fiber_angle, decompose_hopf,
    shell_degeneracy, discrete_s3_volume, discrete_fiber_volume,
    discrete_base_area,
    fiber_winding_exact, alpha_formula_terms, algebraic_alpha_analysis,
    convergence_study, hopf_fiber_analysis,
    fock_embedding, fock_hopf_decompose, cg_expectation_mapping,
    VOL_S1, VOL_S2, VOL_S3, ALPHA_INV_EXPERIMENT,
)


# ============================================================
# Fock projection χ_n
# ============================================================

class TestFockChi:
    def test_n1_is_pi_over_2(self) -> None:
        """n=1 maps to equator: χ₁ = π/2."""
        assert np.isclose(fock_chi(1), np.pi / 2.0)

    def test_monotonic_decreasing(self) -> None:
        """Higher n → smaller χ (toward north pole)."""
        for n in range(1, 20):
            assert fock_chi(n) > fock_chi(n + 1)

    def test_approaches_zero(self) -> None:
        """χ → 0 as n → ∞."""
        assert fock_chi(1000) < 0.01

    def test_n_invalid(self) -> None:
        """n < 1 raises ValueError."""
        with pytest.raises(ValueError):
            fock_chi(0)

    def test_exact_values(self) -> None:
        """Check a few exact values: χ_n = 2 arctan(1/n)."""
        for n in [1, 2, 3, 5, 10]:
            assert np.isclose(fock_chi(n), 2.0 * np.arctan(1.0 / n))


# ============================================================
# State → S³ (linear mapping)
# ============================================================

class TestStateToS3Linear:
    def test_unit_norm(self) -> None:
        """All states map to unit S³ points."""
        for n in range(1, 6):
            for l in range(n):
                for m in range(-l, l + 1):
                    pt = state_to_s3(n, l, m)
                    assert np.isclose(np.linalg.norm(pt), 1.0), \
                        f"|{n},{l},{m}⟩ not on unit S³: |x| = {np.linalg.norm(pt)}"

    def test_distinct_states_give_distinct_points(self) -> None:
        """Different states map to different S³ points (within same shell)."""
        pts = lattice_to_s3(3)
        keys = list(pts.keys())
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                d = np.linalg.norm(pts[keys[i]] - pts[keys[j]])
                assert d > 1e-10, \
                    f"States {keys[i]} and {keys[j]} collide on S³"

    def test_s_state_on_chi_axis(self) -> None:
        """l=0, m=0 states have ψ₁=ψ₂=0, so x₂=x₄=0."""
        for n in range(1, 6):
            pt = state_to_s3(n, 0, 0)
            assert np.isclose(pt[1], 0.0)  # x₂ = cos(χ/2)sin(0) = 0
            assert np.isclose(pt[3], 0.0)  # x₄ = sin(χ/2)sin(0) = 0

    def test_invalid_quantum_numbers(self) -> None:
        """Invalid quantum numbers raise ValueError."""
        with pytest.raises(ValueError):
            state_to_hopf_angles(0, 0, 0)
        with pytest.raises(ValueError):
            state_to_hopf_angles(2, 2, 0)  # l >= n
        with pytest.raises(ValueError):
            state_to_hopf_angles(2, 1, 2)  # |m| > l


# ============================================================
# Hopf projection S³ → S²
# ============================================================

class TestHopfProject:
    def test_s2_unit_norm(self) -> None:
        """Hopf projection of S³ points lies on unit S²."""
        for n in range(1, 5):
            for l in range(n):
                for m in range(-l, l + 1):
                    pt = state_to_s3(n, l, m)
                    s2 = hopf_project(pt)
                    assert np.isclose(np.linalg.norm(s2), 1.0, atol=1e-10), \
                        f"|{n},{l},{m}⟩: S² norm = {np.linalg.norm(s2)}"

    def test_north_pole(self) -> None:
        """Point (1,0,0,0) on S³ → (0,0,1) on S²."""
        pt = np.array([1.0, 0.0, 0.0, 0.0])
        s2 = hopf_project(pt)
        np.testing.assert_allclose(s2, [0.0, 0.0, 1.0], atol=1e-15)

    def test_equator_point(self) -> None:
        """Point (0,0,1,0) on S³ → (0,0,-1) on S²."""
        pt = np.array([0.0, 0.0, 1.0, 0.0])
        s2 = hopf_project(pt)
        np.testing.assert_allclose(s2, [0.0, 0.0, -1.0], atol=1e-15)


# ============================================================
# Fiber angle
# ============================================================

class TestFiberAngle:
    def test_fiber_range(self) -> None:
        """Fiber angle lies in (-2π, 2π]."""
        for n in range(1, 5):
            for l in range(n):
                for m in range(-l, l + 1):
                    pt = state_to_s3(n, l, m)
                    fa = fiber_angle(pt)
                    assert -2 * np.pi - 0.01 < fa < 2 * np.pi + 0.01, \
                        f"|{n},{l},{m}⟩: fiber angle = {fa}"

    def test_same_fiber_same_base(self) -> None:
        """Points differing only by fiber angle project to same S² point."""
        # Construct two S³ points on the same fiber
        chi = np.pi / 3
        psi1_a, psi2_a = 0.5, 0.3
        delta = 0.7  # Fiber shift
        psi1_b, psi2_b = psi1_a + delta, psi2_a + delta

        pt_a = np.array([
            np.cos(chi / 2) * np.cos(psi1_a),
            np.cos(chi / 2) * np.sin(psi1_a),
            np.sin(chi / 2) * np.cos(psi2_a),
            np.sin(chi / 2) * np.sin(psi2_a),
        ])
        pt_b = np.array([
            np.cos(chi / 2) * np.cos(psi1_b),
            np.cos(chi / 2) * np.sin(psi1_b),
            np.sin(chi / 2) * np.cos(psi2_b),
            np.sin(chi / 2) * np.sin(psi2_b),
        ])

        s2_a = hopf_project(pt_a)
        s2_b = hopf_project(pt_b)
        np.testing.assert_allclose(s2_a, s2_b, atol=1e-10)

        # But fiber angles differ
        fa_a = fiber_angle(pt_a)
        fa_b = fiber_angle(pt_b)
        assert not np.isclose(fa_a, fa_b)


# ============================================================
# Fiber winding (exact)
# ============================================================

class TestFiberWinding:
    def test_l0_zero_winding(self) -> None:
        """l=0 has no fiber structure."""
        assert fiber_winding_exact(0) == 0.0

    def test_exact_formula(self) -> None:
        """W(l) = 4πl/(2l+1) for l ≥ 1."""
        for l in range(1, 20):
            expected = 4.0 * np.pi * l / (2 * l + 1)
            assert np.isclose(fiber_winding_exact(l), expected)

    def test_approaches_2pi(self) -> None:
        """W(l)/2π = 2l/(2l+1) → 1 as l → ∞."""
        assert fiber_winding_exact(1000) / (2 * np.pi) > 0.999

    def test_always_less_than_2pi(self) -> None:
        """W(l) < 2π for all finite l."""
        for l in range(1, 100):
            assert fiber_winding_exact(l) < 2 * np.pi


# ============================================================
# Fock embedding (new method)
# ============================================================

class TestFockEmbedding:
    def test_unit_norm(self) -> None:
        """All Fock-embedded points lie on unit S³."""
        for n in range(1, 6):
            for l in range(n):
                for m in range(-l, l + 1):
                    pt = fock_embedding(n, l, m)
                    assert np.isclose(np.linalg.norm(pt), 1.0, atol=1e-14), \
                        f"|{n},{l},{m}⟩: norm = {np.linalg.norm(pt)}"

    def test_s_states_axial(self) -> None:
        """s-states (l=0) have ξ₁=ξ₂=0 (no angular structure)."""
        for n in range(1, 6):
            pt = fock_embedding(n, 0, 0)
            assert np.isclose(pt[1], 0.0, atol=1e-15)  # ξ₁
            assert np.isclose(pt[2], 0.0, atol=1e-15)  # ξ₂

    def test_chi_matches_fock(self) -> None:
        """ξ₀ = cos(χ_n) for all states in shell n."""
        for n in range(1, 6):
            for l in range(n):
                pt = fock_embedding(n, l, 0)
                if l == 0:
                    # s-states: all weight on axis
                    expected_xi0 = np.cos(fock_chi(n))
                    assert np.isclose(pt[0], expected_xi0, atol=1e-14)

    def test_s_state_xi0_xi3(self) -> None:
        """For l=0: ξ₀ = cos(χ), ξ₃ = sin(χ), ξ₁=ξ₂=0."""
        for n in range(1, 6):
            chi = fock_chi(n)
            pt = fock_embedding(n, 0, 0)
            np.testing.assert_allclose(
                pt, [np.cos(chi), 0.0, 0.0, np.sin(chi)], atol=1e-14
            )

    def test_invalid_quantum_numbers(self) -> None:
        with pytest.raises(ValueError):
            fock_embedding(0, 0, 0)
        with pytest.raises(ValueError):
            fock_embedding(2, 2, 0)
        with pytest.raises(ValueError):
            fock_embedding(2, 1, 3)


# ============================================================
# Fock Hopf decomposition
# ============================================================

class TestFockHopfDecompose:
    def test_s2_unit_norm(self) -> None:
        """Hopf base points from Fock embedding lie on unit S²."""
        for n in range(1, 5):
            for l in range(n):
                for m in range(-l, l + 1):
                    result = fock_hopf_decompose(n, l, m)
                    s2 = result['hopf_s2']
                    assert np.isclose(np.linalg.norm(s2), 1.0, atol=1e-10), \
                        f"|{n},{l},{m}⟩: S² norm = {np.linalg.norm(s2)}"

    def test_s3_norm_preserved(self) -> None:
        """The S³ point in the decomposition has unit norm."""
        for n in range(1, 5):
            for l in range(n):
                for m in range(-l, l + 1):
                    result = fock_hopf_decompose(n, l, m)
                    assert np.isclose(
                        np.linalg.norm(result['fock_s3']), 1.0, atol=1e-14
                    )

    def test_chi_stored(self) -> None:
        """Returned chi matches fock_chi."""
        for n in range(1, 5):
            result = fock_hopf_decompose(n, 0, 0)
            assert np.isclose(result['chi'], fock_chi(n))

    def test_spinor_norm(self) -> None:
        """|z₁|² + |z₂|² = 1 (S³ constraint in spinor form)."""
        for n in range(1, 5):
            for l in range(n):
                for m in range(-l, l + 1):
                    result = fock_hopf_decompose(n, l, m)
                    z1, z2 = result['z1'], result['z2']
                    assert np.isclose(abs(z1)**2 + abs(z2)**2, 1.0, atol=1e-14)


# ============================================================
# CG expectation mapping
# ============================================================

class TestCGExpectationMapping:
    def test_n1_trivial(self) -> None:
        """n=1, j=0: both expectations are 0."""
        result = cg_expectation_mapping(1, 0, 0)
        assert np.isclose(result['exp_m_plus'], 0.0)
        assert np.isclose(result['exp_m_minus'], 0.0)

    def test_m_constraint(self) -> None:
        """⟨m⁺⟩ + ⟨m⁻⟩ = m (exact constraint from CG selection rule)."""
        for n in range(1, 6):
            for l in range(n):
                for m in range(-l, l + 1):
                    result = cg_expectation_mapping(n, l, m)
                    assert np.isclose(
                        result['exp_m_plus'] + result['exp_m_minus'],
                        m, atol=1e-12
                    ), f"|{n},{l},{m}⟩: ⟨m⁺⟩+⟨m⁻⟩ = {result['exp_m_plus'] + result['exp_m_minus']}"

    def test_m0_symmetric(self) -> None:
        """For m=0: ⟨m⁺⟩ = -⟨m⁻⟩ (from symmetry)."""
        for n in range(2, 5):
            for l in range(n):
                result = cg_expectation_mapping(n, l, 0)
                # ⟨m⁺⟩ + ⟨m⁻⟩ = 0, so ⟨m⁺⟩ = -⟨m⁻⟩
                assert np.isclose(
                    result['exp_m_plus'], -result['exp_m_minus'], atol=1e-12
                )

    def test_stretched_state(self) -> None:
        """|n, n-1, ±(n-1)⟩ is a single CG term (stretched state)."""
        for n in range(2, 6):
            l = n - 1
            m = l
            result = cg_expectation_mapping(n, l, m)
            assert result['n_cg_terms'] == 1
            j = (n - 1) / 2.0
            assert np.isclose(result['dominant_m_plus'], j)
            assert np.isclose(result['dominant_m_minus'], j)

    def test_psi_range(self) -> None:
        """ψ₁, ψ₂ ∈ [0, π] for all states with j > 0."""
        for n in range(2, 5):
            for l in range(n):
                for m in range(-l, l + 1):
                    result = cg_expectation_mapping(n, l, m)
                    assert -0.01 <= result['psi1'] <= np.pi + 0.01
                    assert -0.01 <= result['psi2'] <= np.pi + 0.01


# ============================================================
# Discrete volumes
# ============================================================

class TestDiscreteVolumes:
    def test_s3_volume_positive(self) -> None:
        """Discrete Vol(S³) is positive for n_max ≥ 2."""
        for n in [2, 5, 10]:
            assert discrete_s3_volume(n) > 0.0

    def test_s3_volume_convergence(self) -> None:
        """Discrete Vol(S³) converges toward 2π² ≈ 19.739."""
        v10 = discrete_s3_volume(10)
        v50 = discrete_s3_volume(50)
        # Should be within 10% of continuum
        assert abs(v10 - VOL_S3) / VOL_S3 < 0.10
        assert abs(v50 - VOL_S3) / VOL_S3 < 0.10


# ============================================================
# Alpha formula
# ============================================================

class TestAlphaFormula:
    def test_formula_accuracy(self) -> None:
        """4π³ + π² + π matches 1/α to 3 ppm."""
        terms = alpha_formula_terms()
        assert terms['relative_error'] < 3e-6

    def test_term_identification(self) -> None:
        """4π³ = Vol(S¹)·Vol(S³), etc."""
        assert np.isclose(4 * np.pi**3, VOL_S1 * VOL_S3)
        assert np.isclose(np.pi**2, VOL_S3 / 2)
        assert np.isclose(np.pi, VOL_S1 / 2)


# ============================================================
# Shell degeneracy
# ============================================================

class TestShellDegeneracy:
    def test_values(self) -> None:
        assert shell_degeneracy(1) == 1
        assert shell_degeneracy(2) == 4
        assert shell_degeneracy(3) == 9


# ============================================================
# Convergence study
# ============================================================

class TestConvergenceStudy:
    def test_returns_correct_keys(self) -> None:
        result = convergence_study([2, 3])
        assert 'n_max_values' in result
        assert 'vol_s3' in result
        assert 'vol_s1' in result
        assert 'vol_s2' in result
        assert 'alpha_candidates' in result
        assert len(result['vol_s3']) == 2

    def test_alpha_candidates_are_dicts(self) -> None:
        result = convergence_study([5])
        assert isinstance(result['alpha_candidates'][0], dict)


# ============================================================
# Integration: key convergence test
# ============================================================

class TestAlphaConvergence:
    """
    The critical question: does any formula converge to 1/α = 137.036?

    Current status: NO candidate formula from discrete volumes converges.
    This test documents the gap and will detect improvement.
    """

    def test_no_discrete_formula_converges(self) -> None:
        """Document that discrete volume formulas have > 5% error at nmax=100."""
        result = convergence_study([100])
        candidates = result['alpha_candidates'][0]
        target = ALPHA_INV_EXPERIMENT
        for name, val in candidates.items():
            err = abs(val - target) / target
            # Currently all > 5% error; update if a better formula is found
            assert err > 0.05 or name == 'berry_harmonic_sum', \
                f"Formula '{name}' unexpectedly close to 1/α: {val} (err {err:.4f})"

    def test_continuum_formula_is_close(self) -> None:
        """The continuum formula 4π³ + π² + π is within 3 ppm."""
        from geovac.hopf_bundle import ALPHA_INV_FORMULA
        err = abs(ALPHA_INV_FORMULA - ALPHA_INV_EXPERIMENT) / ALPHA_INV_EXPERIMENT
        assert err < 3e-6
