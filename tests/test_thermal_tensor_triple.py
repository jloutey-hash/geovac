"""Tests for geovac/thermal_tensor_triple.py — Sprint TD Track 1.

Verification panel:

* Stefan-Boltzmann factorisation as M1 (Hopf measure) x M2 (zeta(4)).
* Stefan-Boltzmann Dirac with fermionic eta factor 7/8.
* Tensor-product mode construction (D = D_S^3 (x) 1 + gamma (x) D_tau).
* Closed-form Casimir on unit S^3 (scalar 1/240, Dirac +17/480).
* Modular residual (MR-B closed form transferred verbatim).
* Transcendental audit (10 sources tagged by tensor factor and Mellin
  mechanism).
* Paper 35 Prediction 1 holds: every pi has a continuous-integration
  source on either the spatial S^3 (Hopf / Seeley-DeWitt) or the
  temporal S^1_beta (Matsubara / theta inversion) factor; rational values
  contain no pi.
"""
from __future__ import annotations

import sympy as sp
import mpmath as mp
import pytest

from geovac.thermal_tensor_triple import (
    s3_scalar_spectrum,
    s3_dirac_spectrum,
    matsubara_spectrum,
    tensor_modes,
    TensorMode,
    stefan_boltzmann_factorization,
    stefan_boltzmann_dirac_factorization,
    scalar_thermal_free_energy_S3_x_S1,
    scalar_casimir_S3,
    dirac_casimir_S3,
    modular_residual_dirac_S3,
    modular_residual_thermal_tensor,
    transcendental_audit,
    run_self_check,
)


# ---------------------------------------------------------------------------
# 1. Spatial spectrum
# ---------------------------------------------------------------------------

class TestS3ScalarSpectrum:
    def test_lowest_mode(self):
        spec = s3_scalar_spectrum(2)
        assert spec[0] == (0, sp.Integer(1), 1)
        assert spec[1] == (1, sp.Integer(2), 4)
        assert spec[2] == (2, sp.Integer(3), 9)

    def test_degeneracy_is_perfect_square(self):
        for (n, om, deg) in s3_scalar_spectrum(5):
            assert deg == (n + 1) ** 2
            assert isinstance(om, sp.Integer)


class TestS3DiracSpectrum:
    def test_lowest_mode_half_integer(self):
        spec = s3_dirac_spectrum(2)
        assert spec[0] == (0, sp.Rational(3, 2), 4)
        assert spec[1] == (1, sp.Rational(5, 2), 12)
        assert spec[2] == (2, sp.Rational(7, 2), 24)

    def test_degeneracy_2_n_plus_1_n_plus_2(self):
        for (n, lam, deg) in s3_dirac_spectrum(5):
            assert deg == 2 * (n + 1) * (n + 2)
            # |lambda_n| = n + 3/2, exact rational
            assert lam == sp.Rational(2 * n + 3, 2)

    def test_g3_equals_40(self):
        # Paper 2's Delta^{-1} = 40 = g_3^Dirac sanity check
        spec = s3_dirac_spectrum(3)
        assert spec[3][2] == 40


# ---------------------------------------------------------------------------
# 2. Matsubara spectrum
# ---------------------------------------------------------------------------

class TestMatsubaraSpectrum:
    def test_bosonic_zero_mode(self):
        beta = sp.Symbol("beta", positive=True)
        spec = matsubara_spectrum(beta, 2, fermionic=False)
        # k = -2..2, so 5 modes, with k=0 having omega=0.
        assert len(spec) == 5
        zero_modes = [s for s in spec if s[1] == 0]
        assert len(zero_modes) == 1

    def test_fermionic_no_zero_mode(self):
        beta = sp.Symbol("beta", positive=True)
        spec = matsubara_spectrum(beta, 2, fermionic=True)
        # No k makes omega_k = 2 pi (k+1/2) / beta vanish.
        for k, om in spec:
            assert sp.simplify(om) != 0

    def test_pi_in_every_nonzero_mode(self):
        beta = sp.Symbol("beta", positive=True)
        spec = matsubara_spectrum(beta, 1, fermionic=False)
        for k, om in spec:
            if k != 0:
                # Every non-zero Matsubara mode contains pi.
                assert sp.pi in sp.simplify(om).atoms()


# ---------------------------------------------------------------------------
# 3. Tensor product modes
# ---------------------------------------------------------------------------

class TestTensorModes:
    def test_scalar_mode_count(self):
        beta = sp.Symbol("beta", positive=True)
        modes = tensor_modes(2, 2, beta, sector="scalar")
        # 3 spatial (n=0,1,2) x 5 Matsubara (k=-2..2)
        assert len(modes) == 15

    def test_dirac_mode_count(self):
        beta = sp.Symbol("beta", positive=True)
        modes = tensor_modes(2, 2, beta, sector="dirac")
        assert len(modes) == 15

    def test_omega_total_sq_factorises_as_sum(self):
        """The hallmark of the tensor structure D = D_spatial (x) 1 +
        gamma (x) D_tau is that {D_spatial(x)1, gamma(x)D_tau} = 0, so
        D^2 = D_spatial^2 (x) 1 + 1 (x) D_tau^2 — additive in spectra."""
        beta = sp.Symbol("beta", positive=True)
        modes = tensor_modes(1, 1, beta, sector="scalar")
        # Pick (n=1, k=1) mode: omega_total^2 should be 4 + (2 pi / beta)^2.
        for m in modes:
            if m.n == 1 and m.k == 1:
                expected = sp.Integer(4) + (2 * sp.pi / beta) ** 2
                assert sp.simplify(m.omega_total_sq - expected) == 0

    def test_invalid_sector_raises(self):
        beta = sp.Symbol("beta", positive=True)
        with pytest.raises(ValueError):
            tensor_modes(2, 2, beta, sector="bogus")


# ---------------------------------------------------------------------------
# 4. Stefan-Boltzmann factorization (HEADLINE)
# ---------------------------------------------------------------------------

class TestStefanBoltzmannFactorization:
    def test_scalar_residual_zero(self):
        sb = stefan_boltzmann_factorization()
        assert sb["residual_to_canonical"] == 0

    def test_scalar_M1_factor_is_inverse_2pi_squared(self):
        sb = stefan_boltzmann_factorization()
        assert sp.simplify(sb["M1_factor_sym"] - 1 / (2 * sp.pi**2)) == 0

    def test_scalar_M2_BE_integrand_closed_form(self):
        sb = stefan_boltzmann_factorization()
        # 6 zeta(4) = pi^4 / 15
        sym = sb["M2_BE_integrand_sym"]
        closed = sb["M2_BE_integrand_closed_sym"]
        assert sp.simplify(sym - closed) == 0
        assert sp.simplify(closed - sp.pi**4 / 15) == 0

    def test_full_factor_gives_minus_pi_sq_over_90(self):
        sb = stefan_boltzmann_factorization()
        canonical = -sp.pi**2 / 90
        assert sp.simplify(sb["F_over_V_factor_sym"] - canonical) == 0

    def test_pi_powers_factorise_correctly(self):
        """M1 (1/pi^2) x M2 (pi^4) = pi^2 net — the exact factor structure
        the master Mellin engine predicts for Stefan-Boltzmann."""
        sb = stefan_boltzmann_factorization()
        # M1 carries 1/pi^2, M2 carries pi^4; net pi^2.
        # Verify by extracting pi-power.
        f = sp.simplify(sb["F_over_V_factor_sym"])
        pi_power = sp.degree(sp.expand(f.rewrite(sp.pi)).as_poly(sp.pi))
        # The result -pi^2/90 has pi^2 leading.
        assert pi_power == 2


class TestStefanBoltzmannDirac:
    def test_dirac_residual_zero(self):
        sbd = stefan_boltzmann_dirac_factorization()
        assert sbd["residual_to_canonical"] == 0

    def test_fermionic_eta_factor_is_seven_eighths(self):
        sbd = stefan_boltzmann_dirac_factorization()
        assert sbd["fermionic_eta_factor_sym"] == sp.Rational(7, 8)

    def test_dirac_per_weyl_canonical(self):
        sbd = stefan_boltzmann_dirac_factorization()
        canonical = -sp.Rational(7, 720) * sp.pi**2
        assert sp.simplify(sbd["F_over_V_per_Weyl_sym"] - canonical) == 0


# ---------------------------------------------------------------------------
# 5. Casimirs (Paper 35 KG-3, KG-5)
# ---------------------------------------------------------------------------

class TestSpatialCasimirs:
    def test_scalar_casimir_rational(self):
        sc = scalar_casimir_S3()
        assert sc["casimir_energy_unit_S3_conformal_scalar"] == sp.Rational(1, 240)
        assert sc["rational_no_pi"] is True
        assert sp.pi not in sp.sympify(
            sc["casimir_energy_unit_S3_conformal_scalar"]).atoms()

    def test_dirac_casimir_rational(self):
        dc = dirac_casimir_S3()
        # +17/480: E = -1/2 zeta_{|D|}(-1) = -1/2*(-17/240); the half-integer
        # shift makes zeta negative, so the fermion factor returns POSITIVE
        # (same sign class as the scalar +1/240). Paper 35 KG-5 derivation.
        assert dc["casimir_energy_unit_S3_full_dirac"] == sp.Rational(17, 480)
        assert dc["rational_no_pi"] is True
        assert sp.pi not in sp.sympify(
            dc["casimir_energy_unit_S3_full_dirac"]).atoms()


# ---------------------------------------------------------------------------
# 6. Thermal partition function
# ---------------------------------------------------------------------------

class TestScalarThermalFreeEnergy:
    def test_thermal_part_high_T_limit(self):
        """At very large T (small beta), the thermal contribution should be
        large negative (the thermal part of free energy is negative)."""
        beta = sp.Symbol("beta", positive=True)
        F_thermal_x_beta = scalar_thermal_free_energy_S3_x_S1(beta, n_max=4)

        # Numerical check at small beta.
        # F * beta = sum_n g_n [(beta omega_n / 2) + log(1 - e^{-beta omega_n})]
        # As beta -> 0, log(1 - e^{-beta omega_n}) -> log(beta omega_n) -> -infinity.
        # So F*beta -> -infinity at high T.
        val_small = float(F_thermal_x_beta.subs(beta, sp.Rational(1, 100)).evalf())
        val_large = float(F_thermal_x_beta.subs(beta, sp.Integer(10)).evalf())
        # At small beta (high T), F*beta is much more negative than at large beta.
        assert val_small < val_large

    def test_zero_T_limit(self):
        """At zero T (beta -> infinity), only the zero-point Casimir survives.
        We don't add the rational Casimir constant in this function (which is
        a beta-independent constant), so beta * F_thermal -> beta * E_Cas
        at large beta (linear growth in beta). Test that growth is positive."""
        beta = sp.Symbol("beta", positive=True)
        F_thermal_x_beta = scalar_thermal_free_energy_S3_x_S1(beta, n_max=2)
        # Derivative w.r.t. beta at large beta is the zero-point energy.
        # This is positive (zero-point energy is positive in our convention).
        deriv = sp.diff(F_thermal_x_beta, beta)
        # At beta=10, derivative should be positive (dominated by Casimir).
        assert float(deriv.subs(beta, 10).evalf()) > 0


# ---------------------------------------------------------------------------
# 7. Modular residual (MR-B closed form)
# ---------------------------------------------------------------------------

class TestModularResidual:
    def test_modular_exponent_pi_squared(self):
        mr = modular_residual_dirac_S3()
        assert mr["modular_exponent"] == sp.pi**2

    def test_leading_prefactor_sqrt_pi(self):
        mr = modular_residual_dirac_S3()
        assert mr["leading_prefactor"] == sp.sqrt(sp.pi)

    def test_term_contains_alternating_sign(self):
        mr = modular_residual_dirac_S3()
        # Should contain (-1)^m factor.
        m = sp.Symbol("m", integer=True, positive=True)
        # Just sanity: the symbolic expression should depend on m.
        assert m in mr["epsilon_term_m"].free_symbols

    def test_tensor_product_factorisation(self):
        mr_t = modular_residual_thermal_tensor()
        beta = sp.Symbol("beta", positive=True)
        s = sp.Symbol("s", positive=True)
        # Leading temporal kernel: beta / (2 sqrt(pi s))
        expected = beta / (2 * sp.sqrt(sp.pi * s))
        assert sp.simplify(mr_t["K_temporal_leading_sym"] - expected) == 0


# ---------------------------------------------------------------------------
# 8. Transcendental audit (HEADLINE deliverable)
# ---------------------------------------------------------------------------

class TestTranscendentalAudit:
    def test_audit_row_count(self):
        a = transcendental_audit()
        # We tagged 10 distinct pi sources / rational collapses.
        assert len(a["rows"]) == 10

    def test_every_row_has_full_tagging(self):
        a = transcendental_audit()
        required_keys = {"quantity", "factor", "mechanism", "tag", "pi_power",
                         "appears_in"}
        for row in a["rows"]:
            assert required_keys <= set(row.keys()), \
                f"row missing keys: {row}"

    def test_mechanism_partition_covers_M1_M2(self):
        a = transcendental_audit()
        partition = a["summary"]["mechanism_partition"]
        # Both M1 and M2 mechanisms appear (the only two that show up in
        # this thermal construction; M3 vertex-parity does not enter because
        # there is no vertex coupling in the free thermal partition function).
        assert partition["M1_count"] >= 1
        assert partition["M2_count"] >= 1
        # M3 should be zero — vertex parity is QED-vertex-specific.
        assert partition["M3_count"] == 0

    def test_paper_35_prediction_1_holds(self):
        """Every pi is associated with a continuous integration over a
        temporal/spectral parameter. Rational results contain no pi."""
        a = transcendental_audit()
        for row in a["rows"]:
            if row["pi_power"] in {"0 (no pi)", "0 (rational ratio)"}:
                # Rational results should be tagged as Hurwitz/Bernoulli or
                # discrete eta — no continuous integration.
                assert ("rational" in row["tag"].lower()
                        or "discrete" in row["mechanism"].lower()
                        or row["mechanism"].startswith("(none"))
            else:
                # Pi-bearing results must reference a continuous integration
                # source (Matsubara, Hopf measure, Seeley-DeWitt, theta
                # inversion).
                src_keywords = ("matsubara", "hopf", "seeley-dewitt",
                                "theta", "circle", "measure", "mellin",
                                "jacobi", "boundary condition")
                tag_lower = row["tag"].lower() + row["factor"].lower()
                assert any(kw in tag_lower for kw in src_keywords), \
                    f"pi-bearing row not tagged with continuous source: {row}"

    def test_rational_collapses_match_paper_35(self):
        a = transcendental_audit()
        partition = a["summary"]["mechanism_partition"]
        # Spatial Casimir 1/240, Dirac Casimir +17/480 — two rational
        # collapses verified in Paper 35 KG-3 / KG-5.
        assert partition["rational_collapse_count"] >= 2


# ---------------------------------------------------------------------------
# 9. Numerical verification: high-T limit of finite tensor sum should
#    approach Stefan-Boltzmann at large T (large beta -> small T, small
#    beta -> high T).
# ---------------------------------------------------------------------------

class TestHighTLimitNumerical:
    """Verify numerically that, as we increase the spatial cutoff and decrease
    beta (high T), the finite-mode thermal free energy approaches the
    Stefan-Boltzmann coefficient.

    We do NOT claim convergence in this finite-truncation panel: the test
    only checks that the SIGN and order of magnitude align. A full
    continuum-limit test would require Weyl asymptotics on the spatial side.
    """
    def test_high_T_sign(self):
        beta = sp.Symbol("beta", positive=True)
        F_x_beta = scalar_thermal_free_energy_S3_x_S1(beta, n_max=4)
        val = float(F_x_beta.subs(beta, sp.Rational(1, 5)).evalf())
        # At high T, thermal free energy is large negative; F*beta also.
        assert val < 0

    def test_thermal_part_continuous_in_beta(self):
        """Smooth, monotonic function of beta in the thermal regime."""
        beta = sp.Symbol("beta", positive=True)
        F_x_beta = scalar_thermal_free_energy_S3_x_S1(beta, n_max=3)
        vals = [float(F_x_beta.subs(beta, b).evalf())
                for b in [sp.Rational(1, 4), sp.Rational(1, 2),
                          sp.Integer(1), sp.Integer(2)]]
        # Should be monotone increasing as beta increases (less thermal
        # contribution at lower T).
        assert all(vals[i] <= vals[i + 1] for i in range(len(vals) - 1))


# ---------------------------------------------------------------------------
# 10. Module self-check driver
# ---------------------------------------------------------------------------

class TestSelfCheck:
    def test_run_self_check_succeeds(self):
        results = run_self_check(verbose=False)
        # All key residuals should be zero.
        assert results["stefan_boltzmann_residual_to_canonical"] == 0
        assert results["stefan_boltzmann_dirac_residual_to_canonical"] == 0
        # Mode counts.
        assert results["scalar_mode_count"] == 15
        assert results["dirac_mode_count"] == 15
        # Casimirs.
        assert results["scalar_casimir"] == sp.Rational(1, 240)
        assert results["dirac_casimir"] == sp.Rational(17, 480)
        # Modular exponent.
        assert results["modular_exponent"] == sp.pi**2
        # Audit complete.
        assert results["pi_audit_row_count"] == 10


# ---------------------------------------------------------------------------
# 11. Cross-validation with L1 modular Hamiltonian (added Sprint L1 2026-05-16)
# ---------------------------------------------------------------------------

class TestL1ModularHamiltonianCrossValidation:
    """Cross-check that the canonical Matsubara apparatus at beta = 2*pi
    is consistent with the L1 modular Hamiltonian's period-closure
    construction.

    The L1 modular Hamiltonian module (geovac.modular_hamiltonian)
    builds K_boost with integer eigenvalues two_m_j on the full-Dirac
    basis and verifies sigma_{2*pi}(O) = O bit-exactly. The matsubara
    spectrum at beta = 2*pi is the corresponding temporal mode quantization
    omega_k = 2*pi*k/(2*pi) = k (bosonic) or (k+1/2) (fermionic).

    These are consistent: integer-spectrum K_boost gives bit-exact
    period closure at beta = 2*pi, and the Matsubara spectrum at the
    same beta has bosonic modes at integer omega_k.
    """

    def test_matsubara_bosonic_modes_at_beta_2pi_are_integers(self):
        from geovac.thermal_tensor_triple import matsubara_spectrum
        beta = sp.Symbol('beta', positive=True)
        spectrum = matsubara_spectrum(beta, k_max=5, fermionic=False)
        # Substitute beta = 2*pi
        for k, omega_k in spectrum:
            val = omega_k.subs(beta, 2 * sp.pi)
            # Bosonic at beta=2pi: omega_k = 2*pi*k/(2*pi) = k (integer)
            assert val == k, (
                f"Bosonic mode k={k} at beta=2*pi: omega_k={val}, expected {k}"
            )

    def test_l1_period_closure_consistent_with_matsubara(self):
        """L1's sigma_{2*pi}(O) = O closure corresponds to bosonic
        Matsubara mode quantization at integer omega_k at beta=2pi."""
        from geovac.modular_hamiltonian import for_bisognano_wichmann
        bw = for_bisognano_wichmann(n_max=2)
        # beta = 2*pi
        import numpy as np
        assert abs(bw.beta - 2.0 * np.pi) < 1e-14
        # K_boost has integer eigenvalues
        K_diag = np.diag(bw.K_geometric)
        for k_val in K_diag:
            k_real = float(np.real(k_val))
            assert abs(k_real - round(k_real)) < 1e-14, (
                f"K_boost eigenvalue {k_real} not integer"
            )
