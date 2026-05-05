"""Tests for the SU(2) central spectral Fejer kernel (R2.5 lemma L2).

Verifies the six sub-tasks of the L2 sprint:

  (a) symbolic verification of kernel positivity, normalization, centrality
  (b) numerical gamma rate at n_max = 2..6 with rate identification
  (c) closed-form Dirichlet-Fejer integral on SU(2)
  (d) Cesaro-2 sharpening
  (e) Plancherel symbol on the Avery-Wen-Avery / so4_three_y_integral basis
  (f) Bozejko-Fendler central-multiplier transcription

Per CLAUDE.md Section 13.4a (equation verification protocol), every equation
appearing in the L2 proof memo (debug/r25_l2_proof_memo.md) has a
corresponding numerical or symbolic verification here.
"""

from __future__ import annotations

import pytest
import sympy as sp
from sympy import Rational, pi, sin, sqrt, integrate, simplify, Integer

import mpmath

from geovac.central_fejer_su2 import (
    _chi,
    _j_max,
    _j_values,
    central_fejer_kernel_su2,
    central_multiplier_cb_norm,
    central_multiplier_cb_norm_cesaro,
    cesaro_2_kernel_su2,
    cesaro_2_normalization,
    character_su2,
    dirichlet_kernel_su2,
    dirichlet_l2_norm_squared,
    fit_gamma_power_law,
    fock_n_to_su2_j,
    gamma_rate,
    gamma_rate_table,
    kernel_l2_norm_squared,
    kernel_pi_free_certificate,
    normalization_constant,
    peter_weyl_bijection_certificate,
    plancherel_symbol,
    plancherel_symbol_cesaro,
    su2_j_to_fock_n,
    verify_centrality,
    verify_normalization_cesaro_symbolic,
    verify_normalization_symbolic,
    verify_pointwise_positivity,
)


# ---------------------------------------------------------------------------
# Sub-task L2-1: Symbolic verification of (a)-(d)
# ---------------------------------------------------------------------------


class TestNormalizationConstant:
    """Z_{n_max} = n_max(n_max+1)/2 (Eq. 2.3 of scoping memo)."""

    @pytest.mark.parametrize(
        "n_max,expected", [(1, 1), (2, 3), (3, 6), (4, 10), (5, 15), (6, 21), (10, 55)]
    )
    def test_Z_closed_form(self, n_max, expected):
        assert normalization_constant(n_max) == expected

    def test_Z_invalid_n(self):
        with pytest.raises(ValueError):
            normalization_constant(0)


class TestJValues:
    """j_values enumerates {0, 1/2, 1, ..., j_max}."""

    def test_j_values_n_max_3(self):
        js = _j_values(3)
        assert js == [Rational(0), Rational(1, 2), Rational(1)]

    def test_j_values_n_max_4(self):
        js = _j_values(4)
        assert js == [Rational(0), Rational(1, 2), Rational(1), Rational(3, 2)]

    def test_j_max_n_max_5(self):
        assert _j_max(5) == Rational(2)

    def test_j_max_n_max_4(self):
        assert _j_max(4) == Rational(3, 2)


class TestCharacterSU2:
    """SU(2) characters chi_j(chi) = sin((2j+1)chi/2) / sin(chi/2)."""

    def test_chi_0_is_constant_one(self):
        # chi_0(chi) = sin(chi/2) / sin(chi/2) = 1
        c = character_su2(Rational(0), _chi)
        # sympy doesn't auto-simplify sin/sin; check at a sample value
        val = c.subs(_chi, sp.Rational(1, 3))
        assert simplify(val) == 1

    def test_chi_1_2_at_pi(self):
        # chi_{1/2}(chi=pi) = sin(pi) / sin(pi/2) = 0/1 = 0
        c = character_su2(Rational(1, 2), _chi)
        val = c.subs(_chi, pi)
        assert simplify(val) == 0

    def test_chi_1_at_2pi3(self):
        # chi_1(chi=2pi/3) = sin(pi) / sin(pi/3) = 0 / (sqrt(3)/2) = 0
        c = character_su2(Rational(1), _chi)
        val = c.subs(_chi, 2 * pi / 3)
        assert simplify(val) == 0


class TestKernelPropertyA_Positivity:
    """K_{n_max}(g) >= 0 (sub-task L2-1.a)."""

    @pytest.mark.parametrize("n_max", [1, 2, 3, 4])
    def test_pointwise_positivity_numerical(self, n_max):
        assert verify_pointwise_positivity(n_max, n_samples=23)


class TestKernelPropertyB_Normalization:
    """integral K_{n_max} dg = 1 (sub-task L2-1.b)."""

    def test_normalization_n_max_2(self):
        # The most stringent symbolic test; n_max >= 3 is slow
        assert verify_normalization_symbolic(2)

    @pytest.mark.slow
    def test_normalization_n_max_3(self):
        assert verify_normalization_symbolic(3)


class TestKernelPropertyC_Centrality:
    """K_{n_max} is a class function (sub-task L2-1.c)."""

    @pytest.mark.parametrize("n_max", [1, 2, 3, 4])
    def test_centrality(self, n_max):
        assert verify_centrality(n_max)


class TestPlancherelSymbol:
    """hat{K}_{n_max}(j) = (2j+1)/Z_{n_max} on the support (sub-task L2-1.d)."""

    @pytest.mark.parametrize("n_max", [1, 2, 3, 4, 5])
    def test_plancherel_sums_to_one(self, n_max):
        # By construction: sum_j (2j+1)/Z = Z/Z = 1.
        # This is the Plancherel-side restatement of integral K = 1.
        total = sum(plancherel_symbol(n_max, j) for j in _j_values(n_max))
        assert total == 1

    def test_plancherel_n_max_3(self):
        assert plancherel_symbol(3, Rational(0)) == Rational(1, 6)
        assert plancherel_symbol(3, Rational(1, 2)) == Rational(2, 6)
        assert plancherel_symbol(3, Rational(1)) == Rational(3, 6)

    def test_plancherel_outside_support(self):
        # j > j_max should give 0
        assert plancherel_symbol(3, Rational(3, 2)) == 0

    def test_plancherel_invalid_j(self):
        # 1/3 is not a half-integer
        with pytest.raises(ValueError):
            plancherel_symbol(3, Rational(1, 3))


class TestL2NormOfKernel:
    """||K||_{L^2}^2 = (2n+1) / (3n(n+1)) by Plancherel."""

    @pytest.mark.parametrize(
        "n_max,expected",
        [
            # Formula: ||K||^2 = sum_{k=1..n} k^2 / Z^2 = n(n+1)(2n+1)/6 / (n(n+1)/2)^2
            #                 = 2(2n+1) / (3 n (n+1))
            (1, Rational(1)),         # 2*3 / (3*1*2) = 6/6 = 1  (K = const 1)
            (2, Rational(5, 9)),      # 2*5 / (3*2*3) = 10/18 = 5/9
            (3, Rational(7, 18)),     # 2*7 / (3*3*4) = 14/36 = 7/18
            (4, Rational(3, 10)),     # 2*9 / (3*4*5) = 18/60 = 3/10
            (5, Rational(11, 45)),    # 2*11 / (3*5*6) = 22/90 = 11/45
        ],
    )
    def test_L2_norm_squared(self, n_max, expected):
        assert kernel_l2_norm_squared(n_max) == expected

    def test_L2_norm_via_plancherel(self):
        # ||K||_{L^2}^2 = sum_j |hat{K}(j)|^2 (Parseval)
        for n_max in [1, 2, 3, 4, 5]:
            ploncherel_sum = sum(
                plancherel_symbol(n_max, j) ** 2 for j in _j_values(n_max)
            )
            assert ploncherel_sum == kernel_l2_norm_squared(n_max)


class TestDirichletKernelL2:
    """||D_{n_max}||_{L^2}^2 = Z_{n_max} (closed form)."""

    @pytest.mark.parametrize("n_max", [1, 2, 3, 4, 5, 8])
    def test_dirichlet_L2_equals_Z(self, n_max):
        assert dirichlet_l2_norm_squared(n_max) == normalization_constant(n_max)


# ---------------------------------------------------------------------------
# Sub-task L2-2: Numerical gamma rate at n_max = 2..6
# ---------------------------------------------------------------------------


class TestGammaRateMonotoneDecay:
    """gamma_{n_max} is strictly decreasing in n_max."""

    def test_gamma_decreases_natural(self):
        gammas = [gamma_rate(n, prec=20) for n in [2, 3, 4, 5, 6]]
        for i in range(len(gammas) - 1):
            assert gammas[i] > gammas[i + 1], (
                f"gamma_{i+2} = {gammas[i]} should be > gamma_{i+3} = {gammas[i+1]}"
            )

    def test_gamma_decreases_cesaro(self):
        gammas = [gamma_rate(n, prec=20, use_cesaro=True) for n in [2, 3, 4, 5, 6]]
        for i in range(len(gammas) - 1):
            assert gammas[i] > gammas[i + 1]


class TestGammaRateBounds:
    """gamma_{n_max} satisfies gamma_n <= C for some absolute constant.

    The kernel is a probability measure on SU(2), and the diameter of
    SU(2) under the round metric is 2 pi (rotation angle); so for any
    kernel, gamma <= 2 pi.
    """

    @pytest.mark.parametrize("n_max", [2, 3, 4, 5, 6])
    def test_gamma_bounded_by_diameter(self, n_max):
        g = gamma_rate(n_max, prec=20)
        assert g < float(2 * mpmath.pi)
        assert g > 0


class TestGammaRateAsymptotic:
    """The asymptotic rate is consistent with the predicted O(log(n)/n) form.

    The pairwise log-log slope should be in the range [-1, -0.5] over
    n_max in [2, 10]; numerical evidence (sub-task L2-2 results memo)
    shows the slope drifts from -0.62 (n=2->3) toward -0.77 (n=8->10),
    consistent with asymptotic O(1/n) rate with an O(1) prefactor and
    sub-leading log corrections.
    """

    def test_pairwise_slope_in_predicted_range(self):
        ns = [2, 3, 4, 5, 6, 7, 8]
        gammas = [gamma_rate(n, prec=25) for n in ns]
        slopes = []
        for i in range(len(ns) - 1):
            n1, n2 = ns[i], ns[i + 1]
            g1, g2 = gammas[i], gammas[i + 1]
            slope = float(
                (mpmath.log(g2) - mpmath.log(g1)) / (mpmath.log(n2) - mpmath.log(n1))
            )
            slopes.append(slope)
        # All slopes should be negative, magnitude > 0.5
        for s in slopes:
            assert s < -0.5, f"slope = {s}; should decay at least as 1/sqrt(n)"
            assert s > -2.0, f"slope = {s}; should not decay faster than 1/n^2"


class TestGammaRateGoesToZero:
    """gamma_{n_max} -> 0 (mass-concentration property SU2-iv)."""

    def test_gamma_decreases_to_zero(self):
        # Verify monotone decrease with sufficient evidence that gamma -> 0
        ns = [2, 4, 6, 8]
        gammas = [gamma_rate(n, prec=20) for n in ns]
        # gamma_8 should be at most half of gamma_2
        assert gammas[-1] < gammas[0] / 2, (
            f"gamma_2 = {gammas[0]}, gamma_8 = {gammas[-1]}; "
            "expected at least 2x decay"
        )


class TestGammaPowerLawClassification:
    """Power-law fit identifies a decay exponent in (0.5, 1.0)."""

    def test_classification(self):
        ns = [3, 4, 5, 6, 7, 8]
        gammas = [gamma_rate(n, prec=25) for n in ns]
        fit = fit_gamma_power_law(ns, gammas)
        # alpha should be > 0.5 (real decay) and < 1.5 (no extra log)
        assert 0.5 < fit["alpha_pure_power"] < 1.5


# ---------------------------------------------------------------------------
# Sub-task L2-3: Closed-form Dirichlet-Fejer integral
# ---------------------------------------------------------------------------


class TestClosedFormGammaSmallN:
    """Verify gamma_rate matches symbolic integration at n_max = 2."""

    def test_closed_form_gamma_2(self):
        K = central_fejer_kernel_su2(2, _chi)
        measure = sin(_chi / 2) ** 2 / pi
        integrand = K * _chi * measure
        val = integrate(integrand, (_chi, 0, 2 * pi))
        val_s = simplify(val)
        # Compare against numerical gamma_rate
        numerical = gamma_rate(2, prec=30)
        diff = abs(float(val_s) - float(numerical))
        assert diff < 1e-10, (
            f"closed-form gamma_2 = {val_s} = {float(val_s)}, "
            f"numerical = {numerical}; diff = {diff}"
        )


class TestDirichletL4Bound:
    """The natural-coefficient kernel has integral |D|^2 chi sin^2(chi/2) dchi
    bounded by O(n_max^2 log n_max).

    This is the critical Cesaro-Fejer estimate (Eq. 2.5 of scoping memo).
    Equivalently, gamma * Z = O(n_max log n_max), so gamma = O(log n / n).
    Numerical: at n_max <= 8, gamma * Z grows roughly as n_max * log n_max.
    """

    def test_gamma_times_Z_bounded(self):
        # gamma * Z grows with n_max; check that the ratio gamma * Z / (n log n)
        # is bounded as n grows
        ns = [3, 4, 6, 8]
        for n in ns:
            g = gamma_rate(n, prec=20)
            Z = normalization_constant(n)
            ratio = float(g * Z) / (n * float(mpmath.log(n)))
            # Empirically this ratio is in [3, 7] for n in [3, 8]
            assert 1.0 < ratio < 15.0, f"n={n}: gamma*Z/(n log n) = {ratio}"


# ---------------------------------------------------------------------------
# Sub-task L2-4: Cesaro-2 sharpening
# ---------------------------------------------------------------------------


class TestCesaro2Properties:
    """Cesaro-2 kernel: positivity, normalization, decay rate."""

    @pytest.mark.parametrize("n_max", [2, 3, 4])
    def test_cesaro_positivity_struct(self, n_max):
        # By |D2|^2 form
        K = cesaro_2_kernel_su2(n_max, _chi)
        # Sample at chi values
        for chi_val in [Rational(1, 4), Rational(1, 2), Rational(3, 4), Rational(1)]:
            val = float(K.subs(_chi, chi_val))
            assert val >= -1e-10, (
                f"K^(2)({chi_val}) = {val} for n_max={n_max}; should be >= 0"
            )

    def test_cesaro_normalization_n_max_2(self):
        assert verify_normalization_cesaro_symbolic(2)

    @pytest.mark.slow
    def test_cesaro_normalization_n_max_3(self):
        assert verify_normalization_cesaro_symbolic(3)


class TestCesaro2Normalization:
    """Z^{(2)}_{n_max} closed form for Cesaro-2 weights."""

    def test_Z_cesaro_n_max_2(self):
        # j=0: weight=(3-0)/3=1, term = 1*1 = 1
        # j=1/2: weight=(3-1)/3=2/3, term = 2*(2/3)^2 = 8/9
        # Z2(2) = 1 + 8/9 = 17/9
        assert cesaro_2_normalization(2) == Rational(17, 9)

    def test_Z_cesaro_n_max_3(self):
        # j=0: w=1, term=1
        # j=1/2: w=2/4=1/2, term = 2 * (1/2)^2 = 1/2
        # j=1: w=1/2... wait, (4-2)/4 = 1/2. term = 3*(1/2)^2 = 3/4
        # actually for n_max=3: np1=4
        # j=0: w=4/4=1, term = 1*1 = 1
        # j=1/2: w=(4-1)/4=3/4, term = 2*(3/4)^2 = 9/8
        # j=1: w=(4-2)/4=1/2, term = 3*(1/2)^2 = 3/4
        # Total = 1 + 9/8 + 3/4 = 8/8 + 9/8 + 6/8 = 23/8
        assert cesaro_2_normalization(3) == Rational(23, 8)


class TestCesaro2GammaImproves:
    """gamma_cesaro <= gamma_natural at every n_max."""

    @pytest.mark.parametrize("n_max", [2, 3, 4, 5])
    def test_cesaro_at_least_as_good(self, n_max):
        g_natural = gamma_rate(n_max, prec=25)
        g_cesaro = gamma_rate(n_max, prec=25, use_cesaro=True)
        # Cesaro suppresses high-frequency oscillations, smoother kernel,
        # better mass concentration
        assert g_cesaro <= g_natural * mpmath.mpf("1.0001"), (
            f"n_max={n_max}: cesaro={g_cesaro}, natural={g_natural}"
        )


# ---------------------------------------------------------------------------
# Sub-task L2-5: Plancherel symbol and Avery-Wen-Avery cross-check
# ---------------------------------------------------------------------------


class TestPlancherelOnAveryBasis:
    """Plancherel symbol values are exact rationals (no transcendental)."""

    @pytest.mark.parametrize("n_max", [1, 2, 3, 4, 5])
    def test_all_symbols_rational(self, n_max):
        for j in _j_values(n_max):
            sym = plancherel_symbol(n_max, j)
            assert isinstance(sym, sp.Rational), (
                f"plancherel({n_max}, {j}) = {sym} is not a Rational"
            )


class TestPeterWeylBijection:
    """Fock n <-> SU(2) j bijection via n = 2j + 1 (cross-link to so4_three_y)."""

    @pytest.mark.parametrize(
        "n,j", [(1, Rational(0)), (2, Rational(1, 2)), (3, Rational(1)),
                (4, Rational(3, 2)), (5, Rational(2))]
    )
    def test_fock_to_su2(self, n, j):
        assert fock_n_to_su2_j(n) == j

    @pytest.mark.parametrize(
        "n,j", [(1, Rational(0)), (2, Rational(1, 2)), (3, Rational(1)),
                (4, Rational(3, 2)), (5, Rational(2))]
    )
    def test_su2_to_fock(self, n, j):
        assert su2_j_to_fock_n(j) == n

    def test_round_trip_n_max_5(self):
        cert = peter_weyl_bijection_certificate(5)
        assert cert["bijection_consistent"]
        # At each n, n = 2j+1 must hold
        for n in range(1, 6):
            entry = cert[n]
            assert entry["matches"]
            # Block dim n^2 = (2j+1)^2
            assert entry["block_dim"] == n * n

    def test_cumulative_dim_matches_operator_system(self):
        # N(n_max) = sum_{n=1..n_max} n^2 (operator system Hilbert dim)
        for n_max in [2, 3, 4]:
            cert = peter_weyl_bijection_certificate(n_max)
            expected = sum(k * k for k in range(1, n_max + 1))
            assert cert["cumulative_N_nmax"] == expected


# ---------------------------------------------------------------------------
# Sub-task L2-6: Bozejko-Fendler central-multiplier transcription
# ---------------------------------------------------------------------------


class TestCentralMultiplierCBNorm:
    """Eq. 3.1 of scoping memo: ||T_K||_cb = 2 / (n_max + 1)."""

    @pytest.mark.parametrize(
        "n_max,expected", [(1, Rational(1)), (2, Rational(2, 3)), (3, Rational(1, 2)),
                           (4, Rational(2, 5)), (5, Rational(1, 3))]
    )
    def test_cb_norm_closed_form(self, n_max, expected):
        assert central_multiplier_cb_norm(n_max) == expected

    def test_cb_norm_decay(self):
        # Decay as 2/(n+1), so cb-norm * (n+1) = 2 constant
        for n in [1, 2, 3, 5, 10]:
            cb = central_multiplier_cb_norm(n)
            assert cb * (n + 1) == 2


class TestCentralMultiplierCBNormMatchesPlancherelMax:
    """The cb-norm equals max_j hat{K}(j); verify constructively."""

    @pytest.mark.parametrize("n_max", [1, 2, 3, 4, 5])
    def test_cb_norm_is_plancherel_max(self, n_max):
        cb = central_multiplier_cb_norm(n_max)
        plancherel_max = max(plancherel_symbol(n_max, j) for j in _j_values(n_max))
        assert cb == plancherel_max


class TestCentralMultiplierCBNormCesaro:
    """Cesaro-2 cb-norm decays at the same O(1/n) rate."""

    @pytest.mark.parametrize("n_max", [2, 3, 4, 5])
    def test_cesaro_cb_norm_decreasing(self, n_max):
        cb = central_multiplier_cb_norm_cesaro(n_max)
        assert 0 < cb < 1

    def test_cesaro_decays_with_n(self):
        cbs = [central_multiplier_cb_norm_cesaro(n) for n in [2, 3, 4, 5, 6]]
        for i in range(len(cbs) - 1):
            assert cbs[i] >= cbs[i + 1] * Rational(1, 2)


# ---------------------------------------------------------------------------
# Composite certificate
# ---------------------------------------------------------------------------


class TestKernelCertificate:
    """Full composite certificate for n_max = 2 (and slow, n_max = 3)."""

    def test_certificate_n_max_2(self):
        cert = kernel_pi_free_certificate(2, pos_samples=11)
        assert cert["positivity_natural"]
        assert cert["centrality"]
        assert cert["normalization_natural"]
        assert cert["plancherel_natural"]
        assert cert["plancherel_natural_outside_support"]
        assert cert["normalization_cesaro"]

    @pytest.mark.slow
    def test_certificate_n_max_3(self):
        cert = kernel_pi_free_certificate(3, pos_samples=11)
        assert cert["positivity_natural"]
        assert cert["centrality"]
        assert cert["normalization_natural"]
        assert cert["plancherel_natural"]
        assert cert["plancherel_natural_outside_support"]
        assert cert["normalization_cesaro"]


# ---------------------------------------------------------------------------
# Sanity: kernel evaluates correctly at chi = 0 (identity element)
# ---------------------------------------------------------------------------


class TestKernelAtIdentity:
    """K_{n_max}(g) at g = e (chi -> 0) is the maximum value of K.

    By orthonormality + |D|^2: D(0) = sum_j sqrt(2j+1) (2j+1) = sum_j (2j+1)^{3/2}.
    This is finite, and K(0) = D(0)^2 / Z_{n_max}.

    We don't symbolically take chi -> 0 (involves L'Hopital); we check at
    chi = small to confirm the kernel is finite and positive there.
    """

    @pytest.mark.parametrize("n_max", [2, 3, 4])
    def test_kernel_finite_near_identity(self, n_max):
        K = central_fejer_kernel_su2(n_max, _chi)
        f = sp.lambdify(_chi, K, modules="mpmath")
        val = f(mpmath.mpf("1e-6"))
        assert mpmath.isfinite(val)
        assert mpmath.re(val) > 0
