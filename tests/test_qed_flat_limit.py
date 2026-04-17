"""Tests for the flat-space limit of S^3 two-loop QED spectral sums.

Validates:
1. Weyl density law (mode count vs volume integral)
2. R-scaling of the sunset sum (correct power law)
3. R-independence of D(s,R)/R^s for s=4,5
4. Structural match: zeta(3) content stable across truncation/radius
5. Even-s remains pi^{even} at all R (one-loop sanity)
"""

import pytest
import mpmath

mpmath.mp.dps = 50


# ---------------------------------------------------------------------------
# Weyl density
# ---------------------------------------------------------------------------

class TestWeylDensity:
    """Verify that the Dirac mode count matches Weyl's asymptotic law."""

    def test_mode_count_closed_form(self):
        """The exact mode count sum_{n=0}^{N} 2(n+1)(n+2) has closed form
        (N+1)(N+2)(2N+3)/3."""
        from geovac.qed_flat_limit import weyl_density_check
        res = weyl_density_check(n_max=50, R=1.0)
        assert res["closed_form_match"], "Closed form does not match direct sum"

    def test_weyl_density_improves_with_nmax(self):
        """Weyl relative error should decrease as n_max increases."""
        from geovac.qed_flat_limit import weyl_density_check
        err_20 = weyl_density_check(n_max=20, R=1.0)["relative_error"]
        err_100 = weyl_density_check(n_max=100, R=1.0)["relative_error"]
        err_500 = weyl_density_check(n_max=500, R=1.0)["relative_error"]
        assert err_100 < err_20, "Error should decrease with n_max"
        assert err_500 < err_100, "Error should decrease with n_max"

    def test_weyl_density_order_of_magnitude(self):
        """At n_max=100, Weyl should match to within a few percent."""
        from geovac.qed_flat_limit import weyl_density_check
        res = weyl_density_check(n_max=100, R=1.0)
        assert res["relative_error"] < 0.05, (
            f"Weyl error {res['relative_error']:.4f} too large at n_max=100"
        )

    def test_weyl_density_at_large_R(self):
        """Weyl law holds at any R (mode count is R-independent)."""
        from geovac.qed_flat_limit import weyl_density_check
        res_R1 = weyl_density_check(n_max=100, R=1.0)
        res_R5 = weyl_density_check(n_max=100, R=5.0)
        # Mode count is the same (g_n doesn't depend on R)
        assert abs(res_R1["mode_count_exact"] - res_R5["mode_count_exact"]) < 1e-10


# ---------------------------------------------------------------------------
# R-scaling of the sunset sum
# ---------------------------------------------------------------------------

class TestRScaling:
    """Verify that S(R) = R^power * S(1) exactly."""

    def test_sunset_R_scaling_exact(self):
        """S(R) / R^power should equal S(1) exactly (algebraic identity)."""
        from geovac.qed_flat_limit import two_loop_sunset_R_dependent

        res_R1 = two_loop_sunset_R_dependent(n_max=10, R=1.0)
        res_R2 = two_loop_sunset_R_dependent(n_max=10, R=2.0)
        res_R5 = two_loop_sunset_R_dependent(n_max=10, R=5.0)

        # All ratios S(R)/R^power should be equal
        ratio1 = res_R1["S_R_over_R_power"]
        ratio2 = res_R2["S_R_over_R_power"]
        ratio5 = res_R5["S_R_over_R_power"]

        assert abs(ratio1 - ratio2) / abs(ratio1) < 1e-10, (
            f"R-scaling broken: ratio(R=1)={ratio1}, ratio(R=2)={ratio2}"
        )
        assert abs(ratio1 - ratio5) / abs(ratio1) < 1e-10, (
            f"R-scaling broken: ratio(R=1)={ratio1}, ratio(R=5)={ratio5}"
        )

    def test_expected_power(self):
        """At s1=s2=2, photon_exponent=1, the power should be 10."""
        from geovac.qed_flat_limit import two_loop_sunset_R_dependent
        res = two_loop_sunset_R_dependent(n_max=5, R=1.0)
        assert res["expected_R_power"] == 10

    def test_expected_power_custom(self):
        """Custom exponents: s1=1, s2=1, photon_exponent=2 -> power = 2+2+4 = 8."""
        from geovac.qed_flat_limit import two_loop_sunset_R_dependent
        res = two_loop_sunset_R_dependent(n_max=5, R=1.0, s1=1, s2=1,
                                          photon_exponent=2)
        assert res["expected_R_power"] == 8

    def test_flat_space_limit_scaling_variation(self):
        """flat_space_limit_scaling should show zero variation at fixed n_max."""
        from geovac.qed_flat_limit import flat_space_limit_scaling
        res = flat_space_limit_scaling(
            n_max=8, R_values=[1.0, 2.0, 3.0, 5.0],
        )
        assert res["variation"] < 1e-10, (
            f"Non-zero variation {res['variation']} in R-scaling"
        )


# ---------------------------------------------------------------------------
# D(s, R) / R^s is R-independent
# ---------------------------------------------------------------------------

class TestDirichletRIndependence:
    """Verify D(s,R)/R^s = D(s,1) at s=4 (even) and s=5 (odd)."""

    def test_D4_R_independence(self):
        """D(4, R) / R^4 should equal D(4, 1) for any R."""
        from geovac.qed_flat_limit import two_loop_vacuum_energy_R_dependent

        for R in [1.0, 2.0, 0.5, 10.0]:
            res = two_loop_vacuum_energy_R_dependent(s=4, n_max=200, R=R)
            assert res["R_independence_check"] < 1e-10, (
                f"D(4, R={R}) / R^4 deviates from D(4, 1): "
                f"check={res['R_independence_check']}"
            )

    def test_D5_R_independence(self):
        """D(5, R) / R^5 should equal D(5, 1) for any R."""
        from geovac.qed_flat_limit import two_loop_vacuum_energy_R_dependent

        for R in [1.0, 3.0, 0.25]:
            res = two_loop_vacuum_energy_R_dependent(s=5, n_max=200, R=R)
            assert res["R_independence_check"] < 1e-10, (
                f"D(5, R={R}) / R^5 deviates from D(5, 1): "
                f"check={res['R_independence_check']}"
            )

    def test_D4_is_pi_even(self):
        """D(4) = pi^2 - pi^4/12 (even parity => pi^{even}).

        D(4) converges slowly (~1/n_max) because the degeneracy g_n ~ n^2
        makes the sum D(4) ~ sum n^2/n^4 = sum 1/n^2, which converges as
        O(1/n_max). We check structural properties, not tight convergence.
        """
        from geovac.qed_flat_limit import two_loop_vacuum_energy_R_dependent
        res = two_loop_vacuum_energy_R_dependent(s=4, n_max=500, R=1.0)
        assert res["parity"] == "even"
        assert res["transcendental_class"] == "pi^{even}"
        # D(4) converges as O(1/n_max); at n_max=500 error should be < 1%
        assert res["truncation_error"] < 0.01

    def test_D5_is_odd_zeta(self):
        """D(5) = 14*zeta(3) - 31/2*zeta(5) (odd parity => odd-zeta).

        D(5) tail ~ sum n^2/n^5 = sum 1/n^3, converges as O(1/n_max^2).
        At n_max=500: error ~ 4e-6.
        """
        from geovac.qed_flat_limit import two_loop_vacuum_energy_R_dependent
        res = two_loop_vacuum_energy_R_dependent(s=5, n_max=500, R=1.0)
        assert res["parity"] == "odd"
        assert res["transcendental_class"] == "odd-zeta"
        assert res["truncation_error"] < 1e-4


# ---------------------------------------------------------------------------
# zeta(3) structural match
# ---------------------------------------------------------------------------

class TestZeta3Structural:
    """Verify the zeta(3) content of D(5) is structurally consistent."""

    def test_D5_exact_identity(self):
        """D(5) = 14*zeta(3) - 31/2*zeta(5) to high precision."""
        from geovac.qed_flat_limit import extract_zeta3_coefficient
        res = extract_zeta3_coefficient(
            n_max_values=[100, 200], R_values=[1.0],
        )
        assert res["match_error"] < 1e-40, (
            f"D(5) identity mismatch: error = {res['match_error']}"
        )

    def test_D5_convergence(self):
        """Truncated D(5) converges monotonically to exact."""
        from geovac.qed_flat_limit import extract_zeta3_coefficient
        res = extract_zeta3_coefficient(
            n_max_values=[50, 100, 200, 500], R_values=[1.0],
        )
        errors = [c["rel_error"] for c in res["convergence"]]
        # Errors should be monotonically decreasing
        for i in range(len(errors) - 1):
            assert errors[i + 1] < errors[i], (
                f"Convergence not monotonic: err[{res['convergence'][i]['n_max']}]"
                f"={errors[i]}, err[{res['convergence'][i+1]['n_max']}]"
                f"={errors[i+1]}"
            )

    def test_D5_R_independence_structural(self):
        """D(5,R)/R^5 is the same at R=1 and R=3."""
        from geovac.qed_flat_limit import extract_zeta3_coefficient
        res = extract_zeta3_coefficient(
            n_max_values=[200], R_values=[1.0, 3.0],
        )
        r_data = res["R_independence"]
        ratio1 = r_data[0]["D5_R_over_R5"]
        ratio3 = r_data[1]["D5_R_over_R5"]
        assert abs(ratio1 - ratio3) / abs(ratio1) < 1e-10

    def test_D5_contains_zeta3(self):
        """D(5) = 14*zeta(3) - 31/2*zeta(5), confirmed by direct identity.

        The generic PSLQ decomposer may fail on this value because the
        coefficient 31/2 is large and the basis includes many competing
        transcendentals. We verify the exact identity directly instead.
        """
        from geovac.qed_flat_limit import extract_zeta3_coefficient
        res = extract_zeta3_coefficient(
            n_max_values=[200], R_values=[1.0],
        )
        # The direct identity check (14*z3 - 31/2*z5) is the reliable test
        assert res["match_error"] < 1e-40, (
            f"D(5) = 14*zeta(3) - 31/2*zeta(5) identity failed: "
            f"error = {res['match_error']}"
        )


# ---------------------------------------------------------------------------
# Weyl zeta(3) verification
# ---------------------------------------------------------------------------

class TestWeylZeta3:
    """Comprehensive structural verification."""

    def test_D4_is_pi_even_exact(self):
        """D(4) = pi^2 - pi^4/12 exactly via Hurwitz."""
        from geovac.qed_flat_limit import verify_weyl_zeta3
        res = verify_weyl_zeta3(n_max_list=[100])
        assert res["D4_match"] is True

    def test_D5_is_odd_zeta_exact(self):
        """D(5) = 14*zeta(3) - 31/2*zeta(5) exactly via Hurwitz."""
        from geovac.qed_flat_limit import verify_weyl_zeta3
        res = verify_weyl_zeta3(n_max_list=[100])
        assert res["D5_match"] is True

    def test_D6_is_pi_even_exact(self):
        """D(6) = pi^4/3 - pi^6/30 exactly via Hurwitz."""
        from geovac.qed_flat_limit import verify_weyl_zeta3
        res = verify_weyl_zeta3(n_max_list=[100])
        assert res["D6_match"] is True

    def test_parity_pattern(self):
        """Even s => pi^{even}; odd s => odd-zeta for s=4..8."""
        from geovac.qed_flat_limit import verify_weyl_zeta3
        res = verify_weyl_zeta3(n_max_list=[100])
        for entry in res["parity_check"]:
            s = entry["s"]
            if s % 2 == 0:
                assert entry["class"] == "pi^{even}", (
                    f"s={s} should be pi^{{even}}"
                )
            else:
                assert entry["class"] == "odd-zeta", (
                    f"s={s} should be odd-zeta"
                )

    def test_convergence_D4_rate(self):
        """D(4) truncation error should decrease with n_max."""
        from geovac.qed_flat_limit import verify_weyl_zeta3
        res = verify_weyl_zeta3(n_max_list=[50, 100, 200, 500])
        errors = [c["rel_error"] for c in res["convergence_D4"]]
        for i in range(len(errors) - 1):
            assert errors[i + 1] < errors[i]

    def test_convergence_D5_rate(self):
        """D(5) truncation error should decrease with n_max."""
        from geovac.qed_flat_limit import verify_weyl_zeta3
        res = verify_weyl_zeta3(n_max_list=[50, 100, 200, 500])
        errors = [c["rel_error"] for c in res["convergence_D5"]]
        for i in range(len(errors) - 1):
            assert errors[i + 1] < errors[i]


# ---------------------------------------------------------------------------
# One-loop sanity: even-s remains pi^{even}
# ---------------------------------------------------------------------------

class TestOneLoopSanity:
    """Verify that the T9 theorem (even-s => pi^{even}) holds at all R."""

    def test_even_s_pi_even_at_R1(self):
        """D(4, R=1) is pi^{even}."""
        from geovac.qed_flat_limit import two_loop_vacuum_energy_R_dependent
        res = two_loop_vacuum_energy_R_dependent(s=4, n_max=200, R=1.0)
        assert res["transcendental_class"] == "pi^{even}"

    def test_even_s_pi_even_at_R10(self):
        """D(4, R=10) / R^4 is still pi^{even} (same value as R=1)."""
        from geovac.qed_flat_limit import two_loop_vacuum_energy_R_dependent
        res = two_loop_vacuum_energy_R_dependent(s=4, n_max=200, R=10.0)
        assert res["transcendental_class"] == "pi^{even}"
        # And the value matches D(4, R=1)
        assert res["R_independence_check"] < 1e-10

    def test_D6_even(self):
        """D(6) at any R should be pi^{even}."""
        from geovac.qed_flat_limit import two_loop_vacuum_energy_R_dependent
        res = two_loop_vacuum_energy_R_dependent(s=6, n_max=200, R=1.0)
        assert res["parity"] == "even"
        assert res["transcendental_class"] == "pi^{even}"


# ---------------------------------------------------------------------------
# Slow tests: larger computations
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestSlowFlatLimit:
    """Larger computations marked as slow."""

    def test_weyl_density_high_nmax(self):
        """Weyl error at n_max=1000 should be < 1%."""
        from geovac.qed_flat_limit import weyl_density_check
        res = weyl_density_check(n_max=1000, R=1.0)
        assert res["relative_error"] < 0.01

    def test_sunset_scaling_multiple_R(self):
        """Sunset R-scaling holds across many R values."""
        from geovac.qed_flat_limit import flat_space_limit_scaling
        res = flat_space_limit_scaling(
            n_max=15,
            R_values=[0.5, 1.0, 2.0, 5.0, 10.0, 20.0],
        )
        assert res["variation"] < 1e-10

    def test_D5_convergence_fine(self):
        """D(5) convergence with fine n_max grid.

        D(5) tail ~ sum n^2/n^5 = sum 1/n^3, so truncation error is
        O(1/n_max^2). At n_max=1000: error ~ 1e-6.
        """
        from geovac.qed_flat_limit import extract_zeta3_coefficient
        res = extract_zeta3_coefficient(
            n_max_values=[50, 100, 200, 500, 1000],
            R_values=[1.0],
        )
        last_err = res["convergence"][-1]["rel_error"]
        assert last_err < 1e-4, f"D(5) error at n_max=1000: {last_err}"
        # Verify monotonic convergence
        errors = [c["rel_error"] for c in res["convergence"]]
        for i in range(len(errors) - 1):
            assert errors[i + 1] < errors[i], "Convergence not monotonic"
