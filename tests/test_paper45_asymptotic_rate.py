r"""Asymptotic-rate verification for Paper 45's surviving spatial rate (Rem rem:asymp_rate).

Paper: ``papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex''

Surviving spatial state-space-GH rate (Paper 45's Lorentzian Lambda^L is
SU(2) x U(1)_T Krein spectral triples, Paper 45 Theorem~\\ref{thm:main},
Eq. eq:main_bound):

   Lambda^L(T^L_{n_max,N_t,T}, T^L_M)
       <= C_3^joint(n_max, N_t) * gamma^joint_{n_max,N_t,T}
       ->_{(n_max, N_t) -> (oo, oo)} 0,

with the asymptotic-rate sharpening (Paper 45 Remark~\\ref{rem:asymp_rate};
Paper 38 Section L2; Paper 40 Section 3.2):

   lim_{n_max -> oo} n_max * gamma^SU_{n_max} / log(n_max) = 4 / pi.

The existing test suite ``tests/test_lorentzian_propinquity.py'' covers only
panel cells (n_max, N_t) in {(2,3), (3,5), (4,7)}. This test extends to an
asymptotic n_max sweep at N_t = 2*n_max - 1.

Computational scope and rationale
---------------------------------

Building the joint LorentzianTunnelingPair at (n_max, N_t) costs
~O((dim_K)^4) memory (vec-stack for span/rank). Empirical timings (Windows,
default precision):

   (2, 3):   0.2 s
   (3, 5):   3.2 s
   (4, 7):  68.7 s
   (5, 9):  ~60 GB allocation, infeasible
   (6, 11): infeasible

The full propinquity panel is therefore CAPPED at n_max <= 4 (matching the
existing slow-test panel).

For the asymptotic-rate verification at n_max in {2, 3, 4, 5, 6} requested
by the task, we use the fact that the propinquity bound at every panel cell
that we CAN build equals gamma_joint_su2 exactly:

   Lambda_bound = max(C_3 * gamma_su2, gamma_su2, gamma_su2, 0) = gamma_su2.

(C_3 = 1, the Paper 38 Lemma L3 comparison constant (gradient seminorm);
the height term is dominated by the reach terms. The earlier
sqrt(1 - 1/n_max) envelope form is withdrawn (op-norm-false; see
test_paper46_c3_operator_system.py). See compute_lorentzian_propinquity_bound
in geovac/lorentzian_propinquity_compact_temporal.py.)

The SU(2) factor gamma^SU_n_max is the Paper 38 L2 quantitative rate, which
is independent of N_t and which can be computed without building the
tunneling pair. We use this rate (via gamma_rate from
geovac/central_fejer_su2.py) for the n_max in {2, 3, 4, 5, 6} sweep and for
the (4 / pi) log(n) / n leading-constant fit (extended to larger n_max
where the asymptote is approached more closely).

Per CLAUDE.md Section 13.4a, each numerical claim verified here corresponds
to a specific paper equation:

  - Theorem~\\ref{thm:main} (monotone decrease) -> test_lambda_monotone_*
  - Remark~\\ref{rem:asymp_rate} (4/pi leading constant) -> test_4_over_pi_*
  - Honest scope (asymptote approached from above) -> test_above_asymptote
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geovac.central_fejer_su2 import gamma_rate
from geovac.lorentzian_propinquity_compact_temporal import (
    LorentzianPropinquityBound,
    compute_lorentzian_propinquity_bound,
)


# Sweep specifications from the task
#   n_max in {2, 3, 4, 5, 6}, N_t = 2 * n_max - 1
SWEEP_CELLS = [(n, 2 * n - 1) for n in [2, 3, 4, 5, 6]]
SWEEP_N_MAX = [n for (n, _) in SWEEP_CELLS]

# Cells at which we CAN build the full LorentzianTunnelingPair without
# memory blowup (~60 GB at n_max = 5). The existing test_lorentzian_propinquity
# slow-test panel uses (2,3), (3,5), (4,7).
FULL_PAIR_CELLS = [(2, 3), (3, 5), (4, 7)]

# Asymptote: Paper 38 L2 / Paper 40 Section 3.2 universal constant
FOUR_OVER_PI = 4.0 / math.pi


# ---------------------------------------------------------------------------
# Sweep at restricted (n_max, N_t) cells where the full pair is buildable
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestPanelMonotoneDecrease:
    """Verify the spatial-rate monotone-decrease prediction at the buildable
    panel cells.

    The task requests monotone Lambda decrease at n_max in {2,3,4,5,6} with
    N_t = 2*n_max - 1. The (5, 9) and (6, 11) cells are infeasible to BUILD
    (memory >60 GB), so we verify the equivalent claim on the buildable
    cells (2,3), (3,5), (4,7) at the full-pair level, then verify the
    extension to n_max in {2,3,4,5,6} via the gamma^SU asymptotic
    rate (TestAsymptoticRate) which is the dominant component.
    """

    def test_lambda_monotone_full_pair(self):
        """Spatial-rate prediction: gamma^SU decreases strictly with n_max
        at the buildable panel cells (2, 3), (3, 5), (4, 7)."""
        bounds = []
        for (n_max, N_t) in FULL_PAIR_CELLS:
            b = compute_lorentzian_propinquity_bound(n_max=n_max, N_t=N_t)
            bounds.append((n_max, N_t, b.propinquity_bound))
        for (n_a, _, l_a), (n_b, _, l_b) in zip(bounds[:-1], bounds[1:]):
            assert l_b < l_a, (
                f"Non-monotone Lambda: n_max={n_a} -> Lambda={l_a:.6f}, "
                f"n_max={n_b} -> Lambda={l_b:.6f}"
            )

    def test_lambda_equals_gamma_su2_at_panel(self):
        """At every buildable cell, Lambda_bound = gamma_joint_su2 (because
        C_3 <= 1 makes height term sub-dominant; height_P = 0).

        This is the load-bearing identity that justifies using gamma_rate
        directly for the n_max in {2,3,4,5,6} sweep.
        """
        for (n_max, N_t) in FULL_PAIR_CELLS:
            b = compute_lorentzian_propinquity_bound(n_max=n_max, N_t=N_t)
            assert b.propinquity_bound == pytest.approx(
                b.gamma_joint_su2, rel=1e-12
            ), (
                f"cell ({n_max}, {N_t}): Lambda={b.propinquity_bound} "
                f"!= gamma_su2={b.gamma_joint_su2}"
            )


# ---------------------------------------------------------------------------
# Asymptotic n_max sweep via the SU(2) rate (cheap, no tunneling pair)
# ---------------------------------------------------------------------------


class TestSweepMonotoneDecrease:
    """Verify monotone decrease over the requested sweep n_max in {2,3,4,5,6}
    using the gamma^SU rate (which equals the propinquity bound exactly at
    every buildable panel cell -- see TestPanelMonotoneDecrease)."""

    def test_gamma_su2_monotone_sweep(self):
        """Lambda^L = gamma^SU_{n_max} decreases strictly across the sweep."""
        gammas = [float(gamma_rate(n, prec=30)) for n in SWEEP_N_MAX]
        for (n_a, g_a), (n_b, g_b) in zip(
            zip(SWEEP_N_MAX[:-1], gammas[:-1]),
            zip(SWEEP_N_MAX[1:], gammas[1:]),
        ):
            assert g_b < g_a, (
                f"Non-monotone gamma^SU: n_max={n_a} -> gamma={g_a:.6f}, "
                f"n_max={n_b} -> gamma={g_b:.6f}"
            )

    def test_gamma_su2_positive_and_finite(self):
        for n in SWEEP_N_MAX:
            g = float(gamma_rate(n, prec=30))
            assert g > 0.0
            assert math.isfinite(g)


# ---------------------------------------------------------------------------
# Leading-constant verification:\ n * gamma_n / log(n) -> 4 / pi
# ---------------------------------------------------------------------------


class TestAsymptoticConstant:
    """Verify Paper 45 Remark~\\ref{rem:asymp_rate} leading constant 4/pi.

    The convergence n * gamma_n / log(n) -> 4 / pi is from ABOVE and is
    logarithmically slow: at the sweep cells n_max in {2..6} the ratio is
    in the range [3.3, 6.0]. At n_max = 1000 it is still ~1.85.

    Two checks below:

    1. test_above_asymptote_at_sweep: At every n_max in {2..6}, the ratio
       n * gamma / log(n) exceeds 4 / pi (Theorem 1(i) of
       debug/r25_l2_quantitative_rate_memo.md, approach-from-above).

    2. test_above_asymptote_extended: Extended to larger n_max where the
       ratio more closely approaches 4/pi (verifying the leading constant
       is the correct asymptotic value and not e.g. 4/(3 pi) or 1).
    """

    def test_above_asymptote_at_sweep(self):
        """At every n_max in the sweep, n * gamma / log(n) > 4 / pi."""
        for n in SWEEP_N_MAX:
            g = float(gamma_rate(n, prec=30))
            ratio = n * g / math.log(n)
            assert ratio > FOUR_OVER_PI, (
                f"n_max={n}: ratio = {ratio:.6f} <= 4/pi = "
                f"{FOUR_OVER_PI:.6f} (asymptote should not be crossed)"
            )

    def test_ratio_monotone_decreasing_in_sweep(self):
        """The ratio n * gamma / log(n) decreases monotonically over the
        sweep (Theorem 1(i) approached from above)."""
        ratios = []
        for n in SWEEP_N_MAX:
            g = float(gamma_rate(n, prec=30))
            ratios.append(n * g / math.log(n))
        for n_a, r_a, n_b, r_b in zip(
            SWEEP_N_MAX[:-1], ratios[:-1],
            SWEEP_N_MAX[1:], ratios[1:],
        ):
            assert r_b < r_a, (
                f"Non-monotone ratio: n={n_a} -> {r_a:.6f}, "
                f"n={n_b} -> {r_b:.6f}"
            )

    @pytest.mark.slow
    def test_above_asymptote_extended(self):
        """Extended verification at larger n_max: ratio remains above 4/pi
        and approaches it. Verifies that 4/pi is the correct leading
        constant in Paper 45 Remark~\\ref{rem:asymp_rate}."""
        # Larger n_max where the asymptotic regime is more visible.
        extended_ns = [20, 50, 100, 200, 500]
        ratios = []
        for n in extended_ns:
            g = float(gamma_rate(n, prec=30))
            r = n * g / math.log(n)
            ratios.append(r)
            assert r > FOUR_OVER_PI, (
                f"n_max={n}: ratio = {r:.6f} <= 4/pi = "
                f"{FOUR_OVER_PI:.6f}"
            )
        # The ratio should decrease monotonically toward 4/pi.
        for r_a, r_b in zip(ratios[:-1], ratios[1:]):
            assert r_b < r_a, "extended ratios non-monotone"

    @pytest.mark.slow
    def test_leading_constant_fit(self):
        r"""Fit gamma_n at the trailing points to (4/pi) * log(n) / n + C / n.

        We fit gamma_n = a * log(n) / n + b / n (two-parameter least-squares
        in 1/n with weights log(n) and 1) at the trailing extended points
        and verify that the fitted leading constant 'a' agrees with 4/pi
        within a tolerance scaled to the finite-n truncation error.

        The convergence is logarithmically slow: at n = 1000 the apparent
        leading constant (n * gamma / log(n)) is still ~1.85, vs the true
        asymptote 4/pi ~= 1.273. We therefore use a paired-window estimator
        a_est(n, 2n) that exploits the gamma_n = a * log(n)/n + O(1/n)
        form to extract 'a':

            n * gamma_n     = a log(n)     + b + O(log n / n)
           2n * gamma_{2n}  = a log(2n)    + b + O(log n / n)
           ----------------------------------
           (subtract:) 2n*g_{2n} - n*g_n   = a log(2)      + O(log n / n)

           => a_est = (2n*gamma_{2n} - n*gamma_n) / log(2)

        which converges to 'a' as n -> infinity at rate O(log n / n).
        """
        # Doubling-pair estimator at trailing windows
        pairs = [(50, 100), (100, 200), (200, 400)]
        a_ests = []
        for (n, two_n) in pairs:
            g_n = float(gamma_rate(n, prec=40))
            g_2n = float(gamma_rate(two_n, prec=40))
            a_est = (two_n * g_2n - n * g_n) / math.log(2.0)
            a_ests.append((n, two_n, a_est))

        # The trailing-window estimator should be CLOSE to 4/pi.
        # We accept a generous tolerance because the finite-n correction
        # is O(log n / n) ~ O(0.02-0.04) at n in the range tested.
        final_n, final_two_n, final_a = a_ests[-1]
        residual = abs(final_a - FOUR_OVER_PI)
        # At n = 200, log(200)/200 ~ 0.026, so we accept residuals
        # up to ~10 * log(n)/n to cover the O(log n / n) correction term
        # and the finite-precision floor.
        tol = 10.0 * math.log(final_n) / final_n
        assert residual < tol, (
            f"Fitted leading constant {final_a:.6f} at doubling pair "
            f"({final_n}, {final_two_n}) deviates from 4/pi = "
            f"{FOUR_OVER_PI:.6f} by {residual:.6f} > tol {tol:.6f}. "
            f"a_ests: {a_ests}"
        )
        # Each successive doubling should bring a_est closer to 4/pi
        # (NOT strictly monotonically since the correction is O(log n / n),
        # but each pair should be within the cumulative tolerance window).
        for (n, _, a) in a_ests:
            tol_n = 20.0 * math.log(n) / n
            assert abs(a - FOUR_OVER_PI) < tol_n, (
                f"At doubling start n = {n}, a_est = {a:.6f} deviates "
                f"from 4/pi by {abs(a - FOUR_OVER_PI):.6f} > "
                f"tol {tol_n:.6f}"
            )


# ---------------------------------------------------------------------------
# Documented infeasibility:\ the (5, 9) and (6, 11) cells cannot be built
# at full-pair level. This is recorded as a test so future hardware upgrades
# or algorithmic improvements can lift the cap.
# ---------------------------------------------------------------------------


class TestPanelInfeasibilityCap:
    """Document the empirical cap on full-pair panel computation."""

    def test_panel_cap_at_n_max_4(self):
        """At n_max = 4, N_t = 7 the full pair builds in ~ 1 min and is
        the largest cell currently feasible. n_max = 5, N_t = 9 requires
        ~60 GB; n_max = 6, N_t = 11 is unbounded. The asymptotic rate
        verification (TestAsymptoticConstant) therefore uses the
        cheap gamma_rate function, which is the dominant component of
        Lambda by TestPanelMonotoneDecrease.test_lambda_equals_gamma_su2_at_panel.
        """
        # No computation; this documents the test architecture choice.
        max_buildable = max(n for (n, _) in FULL_PAIR_CELLS)
        assert max_buildable == 4
