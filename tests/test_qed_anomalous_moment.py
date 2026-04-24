"""Tests for the three-point vertex correction and anomalous magnetic moment on S^3.

Validates the m_j-resolved CG vertex correction in geovac.qed_anomalous_moment.

Key physical results:
  1. B = 0 for the two-point function (no probe, q_probe=0): the vertex
     correction is m_j-independent by the Wigner-Eckart theorem.
  2. B != 0 for the three-point function with q_probe=1: the magnetic
     form factor F_2 is extracted from the m_j-dependent piece.
  3. q_probe=2 gives B = 0 by parity: n_int + n_int + q_probe must be
     odd, and q_probe=2 forces this to be even.
  4. F_2/Schwinger converges to ~1.08 (curvature-shifted from flat-space 1).
  5. Tree-level V_magnetic != 0 while V_charge = 0 for q_probe=1 at n=1.

Test structure:
1. Vertex selection rule and CG infrastructure
2. Tree-level probe coupling
3. Two-point structural zero (Wigner-Eckart)
4. Three-point nonzero magnetic signal
5. Even-parity zero (q_probe=2)
6. Convergence of F_2/Schwinger
7. Cross-validation with debug data
"""

import pytest
from sympy import Rational, S

from geovac.qed_anomalous_moment import (
    half_integer_range,
    vertex_allowed,
    vertex_amp_polarized,
    get_photon_channels,
    tree_level_probe_coupling,
    tree_level_probe_magnetic,
    compute_vertex_3pt_single_level,
    compute_anomalous_magnetic_moment,
    anomalous_moment_convergence,
    mode_dependent_analysis,
)


# ---------------------------------------------------------------------------
# Vertex selection rule
# ---------------------------------------------------------------------------

class TestVertexAllowed:
    """Parity and triangle constraints on the SO(4) vertex."""

    def test_odd_parity_allowed(self) -> None:
        assert vertex_allowed(1, 1, 1)  # 1+1+1=3 odd
        assert vertex_allowed(1, 2, 2)  # 1+2+2=5 odd
        assert vertex_allowed(2, 1, 2)  # 2+1+2=5 odd
        assert not vertex_allowed(0, 1, 1)  # 0+1+1=2 even

    def test_even_parity_forbidden(self) -> None:
        assert not vertex_allowed(1, 1, 2)  # 1+1+2=4 even
        assert not vertex_allowed(2, 2, 2)  # 2+2+2=6 even
        assert not vertex_allowed(0, 0, 2)  # 0+0+2=2 even

    def test_triangle_inequality(self) -> None:
        assert not vertex_allowed(0, 0, 1)  # q > n1 + n2
        assert not vertex_allowed(0, 3, 1)  # q < |n1 - n2|

    def test_q_must_be_positive(self) -> None:
        assert not vertex_allowed(1, 1, 0)
        assert not vertex_allowed(1, 1, -1)

    def test_n_ext_0_has_no_self_coupling(self) -> None:
        """n=0 cannot couple to itself at any q (structural zero)."""
        for q in range(1, 20):
            assert not vertex_allowed(0, 0, q)


# ---------------------------------------------------------------------------
# Half-integer range helper
# ---------------------------------------------------------------------------

class TestHalfIntegerRange:

    def test_j_half(self) -> None:
        result = half_integer_range(Rational(1, 2))
        assert result == [Rational(-1, 2), Rational(1, 2)]

    def test_j_one(self) -> None:
        result = half_integer_range(Rational(1))
        assert result == [Rational(-1), Rational(0), Rational(1)]

    def test_j_zero(self) -> None:
        result = half_integer_range(Rational(0))
        assert result == [Rational(0)]


# ---------------------------------------------------------------------------
# Photon channels
# ---------------------------------------------------------------------------

class TestPhotonChannels:

    def test_q1_n1_channels(self) -> None:
        """n=1, q=1 should have allowed channels."""
        jL = Rational(2, 2)  # (n+1)/2 = 1
        jR = Rational(1, 2)  # n/2 = 1/2
        channels = get_photon_channels(1, 1, 1, jL, jR, jL, jR)
        assert len(channels) > 0

    def test_forbidden_returns_empty(self) -> None:
        """Even parity q_probe=2 at n=1 gives no channels."""
        jL = Rational(2, 2)
        jR = Rational(1, 2)
        channels = get_photon_channels(1, 1, 2, jL, jR, jL, jR)
        assert channels == []


# ---------------------------------------------------------------------------
# Tree-level probe coupling
# ---------------------------------------------------------------------------

class TestTreeLevelProbe:
    """Tree-level probe coupling properties."""

    def test_magnetic_nonzero_q1(self) -> None:
        """V_tree_magnetic at n=1, j=1/2 should be nonzero for q_probe=1."""
        v_mag = tree_level_probe_magnetic(1, Rational(1, 2), q_probe=1)
        assert v_mag != 0, "Tree-level magnetic coupling should be nonzero"
        assert float(v_mag) > 0, "Should be positive"

    def test_magnetic_matches_debug_value(self) -> None:
        """V_magnetic should match the debug data value."""
        v_mag = tree_level_probe_magnetic(1, Rational(1, 2), q_probe=1)
        assert abs(float(v_mag) - 0.5579088621) < 1e-6

    def test_charge_is_zero_q1(self) -> None:
        """V_charge = [V(+1/2) + V(-1/2)]/2 should be zero at n=1, q_probe=1.

        The q=1 probe couples purely magnetically at the ground Dirac level.
        """
        v_up = tree_level_probe_coupling(1, Rational(1, 2), Rational(1, 2))
        v_dn = tree_level_probe_coupling(1, Rational(1, 2), Rational(-1, 2))
        v_charge = (v_up + v_dn) / 2
        assert v_charge == 0, (
            f"V_charge should be zero, got {float(v_charge)}")

    def test_q_probe_2_magnetic_zero(self) -> None:
        """Even q_probe gives zero magnetic coupling by parity."""
        v_mag = tree_level_probe_magnetic(1, Rational(1, 2), q_probe=2)
        assert v_mag == 0, (
            f"q_probe=2 should give zero magnetic coupling, got {float(v_mag)}")

    def test_n0_probe_forbidden(self) -> None:
        """n=0 has no self-coupling at q_probe=1 (parity)."""
        v = tree_level_probe_coupling(0, Rational(1, 2), Rational(1, 2), q_probe=1)
        assert v == 0

    def test_different_n_states(self) -> None:
        """Probe coupling at different n levels are distinct when allowed."""
        allowed_levels = []
        for n in range(5):
            if vertex_allowed(n, n, 1):
                v = tree_level_probe_magnetic(n, Rational(1, 2), q_probe=1)
                if v != 0:
                    allowed_levels.append((n, float(v)))
        assert len(allowed_levels) >= 2, "Need at least 2 allowed levels"
        assert allowed_levels[0][1] != allowed_levels[1][1], (
            "Different n levels should give different magnetic couplings")


# ---------------------------------------------------------------------------
# Two-point structural zero (Wigner-Eckart)
# ---------------------------------------------------------------------------

class TestTwoPointZero:
    """B = 0 for the two-point function (no probe insertion).

    The two-point vertex correction L(m_j) is m_j-independent by the
    Wigner-Eckart theorem (it's a scalar operator on S^3), so the
    magnetic difference B = L(+1/2) - L(-1/2) = 0.

    In the code this manifests through q_probe: a probe that selects
    out the magnetic piece must have odd q_probe. With no probe
    insertion, the internal state self-coupling at q_probe=0 trivially
    gives zero because vertex_allowed(n, n, 0) is always False.
    """

    def test_q_probe_0_gives_zero_tree(self) -> None:
        """Tree-level at q_probe=0 gives zero."""
        v = tree_level_probe_magnetic(1, Rational(1, 2), q_probe=0)
        assert v == 0

    def test_q_probe_0_gives_zero_vertex_correction(self) -> None:
        """Three-point vertex at q_probe=0 gives zero for each level."""
        for n_int in range(4):
            c_up = compute_vertex_3pt_single_level(
                1, Rational(1, 2), Rational(1, 2), n_int, q_probe=0)
            c_dn = compute_vertex_3pt_single_level(
                1, Rational(1, 2), Rational(-1, 2), n_int, q_probe=0)
            assert c_up - c_dn == 0, (
                f"q_probe=0 should give B=0 at n_int={n_int}")


# ---------------------------------------------------------------------------
# Three-point nonzero magnetic signal
# ---------------------------------------------------------------------------

class TestThreePointNonzero:
    """B != 0 for the three-point function with q_probe=1."""

    def test_n_int_1_gives_nonzero_B(self) -> None:
        """The n_int=1 internal level contributes nonzero B."""
        c_up = compute_vertex_3pt_single_level(
            1, Rational(1, 2), Rational(1, 2), 1, q_probe=1)
        c_dn = compute_vertex_3pt_single_level(
            1, Rational(1, 2), Rational(-1, 2), 1, q_probe=1)
        b = c_up - c_dn
        assert b != 0, "n_int=1 should give nonzero B"
        assert float(b) > 0, "B should be positive at n_int=1"

    def test_n_int_0_gives_zero_B(self) -> None:
        """n_int=0 gives zero B (vertex parity forbids self-coupling)."""
        c_up = compute_vertex_3pt_single_level(
            1, Rational(1, 2), Rational(1, 2), 0, q_probe=1)
        c_dn = compute_vertex_3pt_single_level(
            1, Rational(1, 2), Rational(-1, 2), 0, q_probe=1)
        assert c_up - c_dn == 0

    def test_anomalous_moment_nonzero(self) -> None:
        """Full anomalous moment computation gives nonzero F2."""
        result = compute_anomalous_magnetic_moment(n_ext=1, n_max=1, q_probe=1)
        assert result["F2"] != 0
        assert result["F2"] > 0

    def test_b_level_1_matches_debug(self) -> None:
        """B contribution from n_int=1 matches debug data."""
        c_up = compute_vertex_3pt_single_level(
            1, Rational(1, 2), Rational(1, 2), 1, q_probe=1)
        c_dn = compute_vertex_3pt_single_level(
            1, Rational(1, 2), Rational(-1, 2), 1, q_probe=1)
        b = float(c_up - c_dn)
        assert abs(b - 6.9578e-4) < 1e-7, (
            f"B_level at n_int=1 should be ~6.958e-4, got {b}")


# ---------------------------------------------------------------------------
# Even-parity zero (q_probe=2)
# ---------------------------------------------------------------------------

class TestEvenParityZero:
    """q_probe=2 gives B = 0 because even parity is forbidden."""

    def test_q_probe_2_no_allowed_self_coupling(self) -> None:
        """No Dirac level can self-couple at q_probe=2."""
        for n in range(10):
            assert not vertex_allowed(n, n, 2), (
                f"n={n} should NOT self-couple at q_probe=2")

    def test_q_probe_2_tree_magnetic_zero(self) -> None:
        v_mag = tree_level_probe_magnetic(1, Rational(1, 2), q_probe=2)
        assert v_mag == 0

    def test_q_probe_2_vertex_correction_zero(self) -> None:
        """Full anomalous moment at q_probe=2 gives zero B."""
        result = compute_anomalous_magnetic_moment(n_ext=1, n_max=3, q_probe=2)
        assert result["B_magnetic"] == 0.0
        assert result["V_tree_magnetic"] == 0.0


# ---------------------------------------------------------------------------
# Convergence of F2/Schwinger
# ---------------------------------------------------------------------------

class TestConvergence:
    """F2/Schwinger should converge as n_max increases."""

    def test_f2_increases_with_n_max(self) -> None:
        """F2 should increase monotonically (all contributions are positive)."""
        r1 = compute_anomalous_magnetic_moment(n_ext=1, n_max=1, q_probe=1)
        r3 = compute_anomalous_magnetic_moment(n_ext=1, n_max=3, q_probe=1)
        assert r3["F2"] > r1["F2"], (
            f"F2(n_max=3)={r3['F2']} should exceed F2(n_max=1)={r1['F2']}")

    @pytest.mark.slow
    def test_f2_schwinger_ratio_converges(self) -> None:
        """F2/Schwinger at n_max=3 and n_max=5 should be within 1%."""
        r3 = compute_anomalous_magnetic_moment(n_ext=1, n_max=3, q_probe=1)
        r5 = compute_anomalous_magnetic_moment(n_ext=1, n_max=5, q_probe=1)
        ratio_3 = r3["F2_over_schwinger"]
        ratio_5 = r5["F2_over_schwinger"]
        rel_change = abs(ratio_5 - ratio_3) / abs(ratio_3)
        assert rel_change < 0.01, (
            f"F2/Schwinger not converging: n_max=3 gives {ratio_3:.6f}, "
            f"n_max=5 gives {ratio_5:.6f}, rel change {rel_change:.4f}")

    def test_convergence_function_runs(self) -> None:
        """anomalous_moment_convergence produces valid results."""
        results = anomalous_moment_convergence([1, 2, 3], n_ext=1)
        assert len(results) == 3
        assert results[0]["delta"] is None
        for r in results[1:]:
            assert r["delta"] is not None
            assert r["delta"] >= 0


# ---------------------------------------------------------------------------
# Cross-validation with debug data
# ---------------------------------------------------------------------------

class TestCrossValidation:
    """Check that production code reproduces the debug script values."""

    def test_n_max_3_f2_schwinger_ratio(self) -> None:
        """F2/Schwinger at n_max=3 should be ~1.084."""
        result = compute_anomalous_magnetic_moment(n_ext=1, n_max=3, q_probe=1)
        ratio = result["F2_over_schwinger"]
        assert abs(ratio - 1.084) < 0.002, (
            f"F2/Schwinger at n_max=3 should be ~1.084, got {ratio:.6f}")

    def test_n_max_3_b_magnetic(self) -> None:
        """B_magnetic at n_max=3 should match debug data."""
        result = compute_anomalous_magnetic_moment(n_ext=1, n_max=3, q_probe=1)
        b = result["B_magnetic"]
        assert abs(b - 7.024e-4) < 1e-6, (
            f"B_magnetic at n_max=3 should be ~7.024e-4, got {b}")

    def test_per_level_structure(self) -> None:
        """Per-level data should have correct structure."""
        result = compute_anomalous_magnetic_moment(n_ext=1, n_max=2, q_probe=1)
        levels = result["per_level"]
        assert len(levels) == 3  # n_int = 0, 1, 2
        assert levels[0]["B_level"] == 0.0  # n_int=0 is zero
        assert levels[1]["B_level"] > 0  # n_int=1 is dominant
        assert levels[2]["B_level"] > 0  # n_int=2 contributes


# ---------------------------------------------------------------------------
# CG amplitude properties
# ---------------------------------------------------------------------------

class TestCGAmplitude:
    """Properties of the polarized CG vertex amplitude."""

    def test_diagonal_element_at_j_half(self) -> None:
        """Diagonal amplitude (source=target) is computable and exact."""
        jL = Rational(1)
        jR = Rational(1, 2)
        j = Rational(1, 2)
        channels = get_photon_channels(1, 1, 1, jL, jR, jL, jR)
        assert len(channels) > 0
        jpL, jpR = channels[0]
        amp = vertex_amp_polarized(
            jL, jR, j, Rational(1, 2),
            jL, jR, j, Rational(1, 2),
            jpL, jpR, Rational(0), Rational(0))
        # Result is a sympy expression (exact algebraic)
        assert float(amp) != float('nan')
        assert float(amp) != float('inf')

    def test_symmetry_mj_flip(self) -> None:
        """CG amplitudes for m_j and -m_j are generally different."""
        jL = Rational(1)
        jR = Rational(1, 2)
        j = Rational(1, 2)
        channels = get_photon_channels(1, 1, 1, jL, jR, jL, jR)
        total_up = S.Zero
        total_dn = S.Zero
        for jpL, jpR in channels:
            for mpL in half_integer_range(jpL):
                for mpR in half_integer_range(jpR):
                    total_up += vertex_amp_polarized(
                        jL, jR, j, Rational(1, 2),
                        jL, jR, j, Rational(1, 2),
                        jpL, jpR, mpL, mpR)
                    total_dn += vertex_amp_polarized(
                        jL, jR, j, Rational(-1, 2),
                        jL, jR, j, Rational(-1, 2),
                        jpL, jpR, mpL, mpR)
        # The difference is the magnetic coupling
        diff = total_up - total_dn
        assert diff != 0, "CG difference should detect magnetic coupling"


# ---------------------------------------------------------------------------
# Mode dependence (n_ext = 2, 3)
# ---------------------------------------------------------------------------

class TestModeDependence:
    """F2/Schwinger at higher Dirac modes (n_ext > 1).

    Higher external modes have larger Dirac eigenvalue lambda = n_ext + 3/2,
    which suppresses the curvature correction. The ratio F2/Schwinger
    decreases monotonically with n_ext.

    Reference values from debug/data/alpha_g_minus_2_curvature_expansion.json
    (computed at n_max=5 with exact CG coefficients).
    """

    def test_n_ext_2_v_magnetic(self) -> None:
        """V_magnetic at n_ext=2 matches debug data within 0.1%."""
        v_mag = tree_level_probe_magnetic(2, Rational(1, 2), q_probe=1)
        v_mag_f = float(v_mag)
        expected = 0.38925844503283874
        rel_err = abs(v_mag_f - expected) / expected
        assert rel_err < 0.001, (
            f"V_magnetic at n_ext=2: got {v_mag_f:.10f}, "
            f"expected {expected:.10f}, rel err {rel_err:.2e}")

    def test_n_ext_3_v_magnetic(self) -> None:
        """V_magnetic at n_ext=3 matches debug data within 0.1%."""
        v_mag = tree_level_probe_magnetic(3, Rational(1, 2), q_probe=1)
        v_mag_f = float(v_mag)
        expected = 0.3000988014334038
        rel_err = abs(v_mag_f - expected) / expected
        assert rel_err < 0.001, (
            f"V_magnetic at n_ext=3: got {v_mag_f:.10f}, "
            f"expected {expected:.10f}, rel err {rel_err:.2e}")

    def test_n_ext_2_f2_ratio(self) -> None:
        """F2/Schwinger at n_ext=2 is below 1.0 (higher modes suppressed)."""
        result = compute_anomalous_magnetic_moment(n_ext=2, n_max=3, q_probe=1)
        ratio = result["F2_over_schwinger"]
        assert ratio < 1.0, (
            f"F2/Schwinger at n_ext=2 should be < 1.0, got {ratio:.6f}")
        assert ratio > 0.0, (
            f"F2/Schwinger at n_ext=2 should be positive, got {ratio:.6f}")

    def test_n_ext_2_b_nonzero(self) -> None:
        """B(n_ext=2, n_int=1) is positive and nonzero."""
        c_up = compute_vertex_3pt_single_level(
            2, Rational(1, 2), Rational(1, 2), 1, q_probe=1)
        c_dn = compute_vertex_3pt_single_level(
            2, Rational(1, 2), Rational(-1, 2), 1, q_probe=1)
        b = float(c_up - c_dn)
        assert b > 0, f"B(n_ext=2, n_int=1) should be positive, got {b}"

    def test_mode_dependence_decreasing(self) -> None:
        """F2/Schwinger decreases monotonically with n_ext at n_max=3."""
        ratios = []
        for n_ext in [1, 2, 3]:
            result = compute_anomalous_magnetic_moment(
                n_ext=n_ext, n_max=3, q_probe=1)
            ratios.append(result["F2_over_schwinger"])
        for i in range(len(ratios) - 1):
            assert ratios[i] > ratios[i + 1], (
                f"F2/Schwinger should decrease: "
                f"n_ext={i+1} gives {ratios[i]:.6f}, "
                f"n_ext={i+2} gives {ratios[i+1]:.6f}")

    def test_mode_dependent_analysis_runs(self) -> None:
        """mode_dependent_analysis runs without error and returns valid data."""
        result = mode_dependent_analysis(
            n_ext_values=[1, 2], n_max=2, q_probe=1)
        assert "per_mode" in result
        assert "ratios" in result
        assert "monotone_decreasing" in result
        assert "fit_coeffs" in result
        assert len(result["per_mode"]) == 2
        assert len(result["ratios"]) == 2
        # Check that fit coefficients are present
        assert "c1" in result["fit_coeffs"]
        assert "c2" in result["fit_coeffs"]
        # Check per_mode structure
        for entry in result["per_mode"]:
            assert "n_ext" in entry
            assert "lambda_ext" in entry
            assert "V_magnetic" in entry
            assert "F2_over_schwinger" in entry
            assert "correction" in entry
