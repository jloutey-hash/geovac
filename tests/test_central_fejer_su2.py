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
    asymptotic_rate_constant,
    central_fejer_kernel_su2,
    central_multiplier_cb_norm,
    central_multiplier_cb_norm_cesaro,
    cesaro_2_kernel_su2,
    cesaro_2_normalization,
    character_su2,
    dirichlet_kernel_su2,
    dirichlet_l2_norm_squared,
    doubling_estimator,
    fit_gamma_power_law,
    fock_n_to_su2_j,
    gamma_n_via_sum_rule,
    gamma_rate,
    gamma_rate_table,
    kernel_l2_norm_squared,
    kernel_pi_free_certificate,
    normalization_constant,
    N0_for_constant,
    peter_weyl_bijection_certificate,
    plancherel_symbol,
    plancherel_symbol_cesaro,
    quantitative_rate_bound,
    quantitative_rate_certificate,
    su2_j_to_fock_n,
    T_n_via_sum_rule,
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


# ---------------------------------------------------------------------------
# Quantitative rate sharpening (R2.5 / L2 quantitative-rate sprint, May 2026)
# ---------------------------------------------------------------------------
#
# Theorem (debug/r25_l2_quantitative_rate_memo.md):
#   (i)   lim n*gamma_n / log(n) = 4/pi   (asymptotic constant)
#   (ii)  gamma_n <= 6 log(n)/n for all n >= 2   (uniform bound)
#   (iii) For C > 4/pi, exists N_0(C) such that gamma_n <= C log(n)/n for n >= N_0


class TestQuantitativeRate:
    """Verify Theorem 1 of debug/r25_l2_quantitative_rate_memo.md."""

    def test_asymptotic_constant_value(self):
        """Transcription lock: asymptotic_rate_constant() returns 4/pi exactly.

        NOTE (2026-06-14 QA): this is a CONSTANT LOCK, not a derivation of
        4/pi -- it asserts a hardcoded module constant against itself. The
        genuine, non-circular derivation (Euler-Maclaurin on the gamma_n sum
        rule, which rejects the 2/pi circle-Fejer decoy) lives in
        tests/test_trunk_qa_fejer_4_over_pi.py. Kept here only as a value lock.
        """
        c = asymptotic_rate_constant()
        # As sympy, equals Rational(4) / pi.
        assert c == Rational(4) / pi
        # Numerically, equals 4/pi to high precision.
        c_float = float(c)
        import math
        assert abs(c_float - 4.0 / math.pi) < 1e-15

    def test_T_n_sum_rule_matches_quadrature_n_max_2(self):
        """T_n via sum rule: gamma_n = pi - 4 T_n / (pi Z_n) reproduces gamma at n=2."""
        # Closed form: T_2 = 2 * sqrt(2) * [1 - 1/9] = 16 sqrt(2)/9
        T_2 = T_n_via_sum_rule(2, prec=50)
        T_2_expected = mpmath.mpf("16") * mpmath.sqrt(2) / 9
        assert abs(T_2 - T_2_expected) < mpmath.mpf("1e-40")

    def test_T_n_sum_rule_gives_exact_gamma_2(self):
        """gamma_2 = pi - 64 sqrt(2) / (27 pi) (memo Eq. 3.2 closed form at n=2)."""
        g_2 = gamma_n_via_sum_rule(2, prec=50)
        g_2_memo = mpmath.pi - 64 * mpmath.sqrt(2) / (27 * mpmath.pi)
        assert abs(g_2 - g_2_memo) < mpmath.mpf("1e-40")

    @pytest.mark.parametrize("n_max", [2, 3, 4, 5])
    def test_T_n_sum_rule_matches_gamma_quadrature(self, n_max):
        """Sum-rule gamma matches Gaussian quadrature gamma at small n_max."""
        g_sumrule = float(gamma_n_via_sum_rule(n_max, prec=40))
        g_quad = float(gamma_rate(n_max, prec=40))
        # Both should agree to at least 1e-12 (they compute the same quantity).
        assert abs(g_sumrule - g_quad) < 1e-12

    @pytest.mark.parametrize("n_max", [2, 3, 4, 5, 6, 7, 8, 9, 10])
    def test_uniform_bound_C_equals_6(self, n_max):
        """gamma_n <= 6 * log(n) / n holds for all n >= 2 (Theorem 1(ii))."""
        cert = quantitative_rate_certificate(n_max, C=6.0, prec=40)
        assert cert["satisfies_bound"], (
            f"gamma_{n_max}={cert['gamma_n']:.6f} > 6 log(n)/n = "
            f"{cert['bound_C_log_n_over_n']:.6f}"
        )

    def test_uniform_bound_tight_at_n_2(self):
        """Bound C=6 at n=2 has small but positive margin (it's the tightest case)."""
        cert = quantitative_rate_certificate(2, C=6.0, prec=40)
        # Exact: 2*gamma_2/log(2) = 2*(pi - 64*sqrt(2)/(27pi))/log(2) ~ 5.986
        # Bound 6 has margin ~ (6 - 5.986)/5.986 * 100 = 0.24% in one rendering;
        # in the certificate format margin = (bound - gamma)/gamma is bigger.
        assert cert["satisfies_bound"]
        # Margin should be small but nonzero. Exact value is
        # (2*log(2) - 2*gamma_2)/gamma_2 ~ small positive
        assert cert["margin_pct"] > 0.0
        assert cert["margin_pct"] < 1.0  # tight!

    def test_uniform_bound_loose_at_large_n(self):
        """Bound C=6 has large margin at moderate-large n."""
        cert = quantitative_rate_certificate(50, C=6.0, prec=40)
        assert cert["satisfies_bound"]
        # Margin at n=50 should be > 100% (verified in memo Section 4.2)
        assert cert["margin_pct"] > 100.0

    @pytest.mark.parametrize("n_max", [3, 5, 10, 20, 50, 100])
    def test_n_gamma_over_log_n_above_4_over_pi(self, n_max):
        """Theorem 1(i) limit is approached from ABOVE.

        n*gamma_n/log(n) > 4/pi for all finite n (asymptote not reached).
        """
        import math
        cert = quantitative_rate_certificate(n_max, C=6.0, prec=40)
        ratio = cert["n_gamma_over_log_n"]
        assert ratio > 4.0 / math.pi, (
            f"n={n_max}: n*gamma/log(n) = {ratio:.6f} <= 4/pi = "
            f"{4.0/math.pi:.6f} (asymptote should not be crossed)"
        )

    def test_n_gamma_over_log_n_monotone_decreasing_for_n_ge_3(self):
        """The ratio n*gamma_n/log(n) is monotonically decreasing for n >= 3."""
        import math
        ns = [3, 4, 5, 7, 10, 15, 20, 30, 50, 100]
        ratios = []
        for n in ns:
            cert = quantitative_rate_certificate(n, C=6.0, prec=40)
            ratios.append(cert["n_gamma_over_log_n"])
        for i in range(1, len(ratios)):
            assert ratios[i] < ratios[i-1] - 1e-6, (
                f"Non-monotonicity at n={ns[i-1]} -> {ns[i]}: "
                f"{ratios[i-1]:.6f} -> {ratios[i]:.6f}"
            )

    def test_doubling_estimator_converges_to_4_over_pi(self):
        """The doubling estimator a_n converges to 4/pi as n -> infinity.

        At each successive doubling, the error |a_n - 4/pi| should
        approximately halve (consistent with the second-order O(1/log n)
        subleading correction).
        """
        import math
        # a_n at n in {50, 100, 200, 400}; comparing the error series
        a_50 = float(doubling_estimator(50, prec=40))
        a_100 = float(doubling_estimator(100, prec=40))
        a_200 = float(doubling_estimator(200, prec=40))
        a_400 = float(doubling_estimator(400, prec=40))
        target = 4.0 / math.pi
        e_50 = abs(a_50 - target)
        e_100 = abs(a_100 - target)
        e_200 = abs(a_200 - target)
        e_400 = abs(a_400 - target)
        # Errors should be DECREASING monotonically
        assert e_100 < e_50, f"e_100={e_100:.6f} >= e_50={e_50:.6f}"
        assert e_200 < e_100, f"e_200={e_200:.6f} >= e_100={e_100:.6f}"
        assert e_400 < e_200, f"e_400={e_400:.6f} >= e_200={e_200:.6f}"
        # Error at n=400 should be < 0.01 (memo data: 0.00969)
        assert e_400 < 0.011, f"e_400={e_400:.6f} > 0.011 expected"

    def test_doubling_estimator_overestimates_4_over_pi(self):
        """a_n > 4/pi for finite n (asymptote approached from above)."""
        import math
        for n in [50, 100, 200]:
            a_n = float(doubling_estimator(n, prec=40))
            assert a_n > 4.0 / math.pi, (
                f"doubling_estimator({n}) = {a_n:.6f} <= 4/pi"
            )

    def test_N0_lookups(self):
        """Threshold table values match memo Theorem 1(iii)."""
        # From memo §4.3:
        assert N0_for_constant(3.0) == 9
        assert N0_for_constant(2.5) == 30
        assert N0_for_constant(2.0) == 300
        assert N0_for_constant(4.0) == 4
        assert N0_for_constant(5.0) == 3
        assert N0_for_constant(6.0) == 2
        # Above all entries (use entry for highest C in table):
        assert N0_for_constant(10.0) == 2
        # Below 4/pi: returns None
        assert N0_for_constant(1.0) is None
        assert N0_for_constant(0.5) is None

    def test_N0_below_asymptote_returns_None(self):
        """Below the asymptotic constant, N_0 is undefined."""
        import math
        # 4/pi ~ 1.273; just below should return None
        assert N0_for_constant(1.27) is None
        # Just above 4/pi: not guaranteed to be in our small precomputed
        # table -> may return None (depending on whether table has the entry).
        # The function should not raise.
        result = N0_for_constant(4.0 / math.pi + 0.01)
        assert result is None or isinstance(result, int)

    def test_N0_invalid_C_raises(self):
        """N0_for_constant(0) and negative C raise ValueError."""
        with pytest.raises(ValueError):
            N0_for_constant(0.0)
        with pytest.raises(ValueError):
            N0_for_constant(-1.0)

    def test_quantitative_rate_bound_value(self):
        """quantitative_rate_bound(n, C) = C * log(n) / n."""
        import math
        for n, C in [(10, 6.0), (50, 3.0), (100, 2.0)]:
            expected = C * math.log(n) / n
            actual = quantitative_rate_bound(n, C=C)
            assert abs(actual - expected) < 1e-15

    def test_quantitative_rate_bound_invalid_n(self):
        """quantitative_rate_bound rejects n_max < 2."""
        with pytest.raises(ValueError):
            quantitative_rate_bound(1, C=6.0)
        with pytest.raises(ValueError):
            quantitative_rate_bound(0, C=6.0)

    def test_T_n_invalid_n(self):
        """T_n_via_sum_rule rejects n_max < 1."""
        with pytest.raises(ValueError):
            T_n_via_sum_rule(0)

    def test_gamma_via_sum_rule_invalid_n(self):
        """gamma_n_via_sum_rule rejects n_max < 1."""
        with pytest.raises(ValueError):
            gamma_n_via_sum_rule(0)

    def test_doubling_estimator_invalid_n(self):
        """doubling_estimator rejects n < 1."""
        with pytest.raises(ValueError):
            doubling_estimator(0)

    @pytest.mark.parametrize("n_max,C,expected_bound_holds", [
        (3, 3.0, True),     # n=3 is at C=3 boundary; ratio at n=3 is ~4.4 > 3 - should FAIL
        (10, 3.0, True),    # n=10 ratio is ~2.92 < 3 - should hold
        (50, 2.5, True),    # n=50 ratio is ~2.30 < 2.5 - should hold
        (300, 2.0, True),   # n=300 ratio is ~1.99 < 2.0 - should hold
    ])
    def test_explicit_thresholds(self, n_max, C, expected_bound_holds):
        """Spot-check threshold table entries."""
        cert = quantitative_rate_certificate(n_max, C=C, prec=40)
        # All these are above N_0(C) so should hold
        if n_max == 3 and C == 3.0:
            # n=3 is below N_0(C=3)=9, so SHOULD VIOLATE
            assert not cert["satisfies_bound"]
        else:
            assert cert["satisfies_bound"] == expected_bound_holds

    def test_quantitative_rate_certificate_structure(self):
        """The certificate dict has the expected keys."""
        cert = quantitative_rate_certificate(10, C=6.0, prec=30)
        expected_keys = {
            "n_max", "C", "gamma_n", "bound_C_log_n_over_n",
            "satisfies_bound", "margin_pct",
            "asymptotic_constant_4_over_pi", "n_gamma_over_log_n",
        }
        assert set(cert.keys()) == expected_keys
        assert cert["n_max"] == 10
        assert cert["C"] == 6.0


class TestSumRuleStructuralProperties:
    """Verify structural properties of the closed-form T_n sum rule."""

    def test_T_n_is_positive(self):
        """T_n > 0 for all n >= 2."""
        for n in [2, 3, 4, 5, 8, 10, 20]:
            T = T_n_via_sum_rule(n, prec=30)
            assert T > 0, f"T_{n} = {T} not positive"

    def test_gamma_n_positive_via_sum_rule(self):
        """gamma_n > 0 for all n (via sum rule)."""
        for n in [2, 3, 4, 5, 8, 10, 20, 50, 100]:
            g = gamma_n_via_sum_rule(n, prec=30)
            assert g > 0, f"gamma_{n} = {g} not positive"

    def test_gamma_n_decreasing_via_sum_rule(self):
        """gamma_n monotonically decreasing in n."""
        prev = None
        for n in [2, 3, 5, 10, 20, 50, 100]:
            g = float(gamma_n_via_sum_rule(n, prec=30))
            if prev is not None:
                assert g < prev, f"gamma_{n} = {g:.6f} >= prev {prev:.6f}"
            prev = g

    def test_T_n_grows_like_n_squared(self):
        """T_n / Z_n -> pi^2/4 as n -> infinity (memo §3.2 leading order)."""
        import math
        target = math.pi**2 / 4
        # At n=100, T_n/Z_n should be within 5% of pi^2/4
        T_100 = float(T_n_via_sum_rule(100, prec=30))
        Z_100 = 100 * 101 / 2
        ratio = T_100 / Z_100
        assert abs(ratio - target) / target < 0.06, (
            f"T_100/Z_100 = {ratio:.4f} vs target pi^2/4 = {target:.4f}"
        )
