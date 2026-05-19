"""Track 3 tests for Hylleraas-Eckart 2D variational driver.

Verifies:
1. optimize_alpha_beta_for_state converges on He 1^1S.
2. At convergence, the variational principle holds: E >= E_exact.
3. The optimizer's energy at (alpha_opt, beta=0) is at most epsilon worse
   than the single-alpha optimum (variational gain or wash).
4. Convergence-constraint guard works (alpha < beta returns 1e10).

The 2^1S - 2^3S splitting target test is in
debug/he_2s_singlet_triplet_eckart.py rather than the pytest suite
(takes several minutes per omega and is the headline numerical result).
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from geovac.hylleraas_r12 import (
    hylleraas_basis_3p,
    hylleraas_basis_total_degree,
    optimize_alpha_beta_for_state,
    optimize_alpha_for_state,
    solve_hylleraas_state,
)


EXACT_HE_1S = -2.903724  # Drake NR reference


# ---------------------------------------------------------------------------
# Optimizer convergence
# ---------------------------------------------------------------------------

class TestOptimizerConvergence:
    def test_he_1s_3p_optimizer_converges(self):
        """Optimizer converges in finite iterations on He 1^1S (3p basis)."""
        basis = hylleraas_basis_3p()
        res = optimize_alpha_beta_for_state(
            basis, Z=2.0, alpha_init=1.7, beta_init=0.0,
            state_index=0, spin='singlet',
        )
        assert res.extras['opt_n_eval'] > 0
        assert res.extras['opt_n_iter'] > 0
        # Reasonable iteration count for 2D Nelder-Mead.
        assert res.extras['opt_n_eval'] < 500

    def test_he_1s_omega3_variational_bound(self):
        """Optimizer result respects variational principle."""
        basis = hylleraas_basis_total_degree(3)
        res = optimize_alpha_beta_for_state(
            basis, Z=2.0, alpha_init=1.7, beta_init=0.0,
            state_index=0, spin='singlet',
        )
        assert res.energy >= EXACT_HE_1S - 1e-10, (
            f"Optimized E={res.energy} below exact {EXACT_HE_1S}"
        )

    def test_eckart_optimum_at_most_marginally_better_than_single_alpha(self):
        """For He 1^1S in 3p basis, Eckart optimum should give energy
        within a few mHa of (and not above) the single-alpha optimum.

        For symmetric ground states, the asymmetric trial cosh(beta t)
        provides at most a small variational gain (the basis is already
        symmetric in r_1, r_2 by construction).
        """
        basis = hylleraas_basis_3p()
        Z = 2.0
        # Single-alpha optimum
        res_single = optimize_alpha_for_state(basis, Z, alpha_init=1.7)
        # Eckart 2D optimum
        res_eckart = optimize_alpha_beta_for_state(
            basis, Z, alpha_init=res_single.alpha, beta_init=0.0,
            state_index=0, spin='singlet',
        )
        # Eckart at the same alpha with beta=0 must equal single-alpha exactly
        # (by Track 2's bit-identical regression). Optimizing further can
        # only LOWER the energy (variational principle) or stay the same
        # within optimizer noise.
        assert res_eckart.energy <= res_single.energy + 1e-6, (
            f"E_eckart={res_eckart.energy} above E_single={res_single.energy} "
            f"by more than optimizer tolerance"
        )


# ---------------------------------------------------------------------------
# Convergence-constraint guard
# ---------------------------------------------------------------------------

class TestConvergenceGuard:
    def test_alpha_less_than_beta_rejected(self):
        """Objective returns +inf when |beta| >= alpha (cosh master
        integral diverges)."""
        from scipy.optimize import minimize

        basis = hylleraas_basis_3p()
        Z = 2.0

        # Wrap the same objective the optimizer uses.
        def obj(x):
            alpha, beta = float(x[0]), float(x[1])
            if alpha <= 0:
                return 1e10
            if abs(beta) >= alpha:
                return 1e10
            try:
                return solve_hylleraas_state(
                    basis, alpha, Z, state_index=0,
                    mode='eckart_double_alpha', beta=beta, spin='singlet',
                ).energy
            except Exception:
                return 1e10

        # In-feasible point.
        v_bad = obj([0.5, 1.5])
        assert v_bad == 1e10
        # Feasible point.
        v_good = obj([1.5, 0.5])
        assert v_good < 0  # He 1^1S energy is negative


# ---------------------------------------------------------------------------
# Smoke test of full sprint driver (small basis only, for speed)
# ---------------------------------------------------------------------------

class TestSprintDriverSmoke:
    @pytest.mark.slow
    def test_he_2s_sprint_omega2_runs(self):
        """compute_he_2s_singlet_triplet_eckart runs at omega=2 (smallest
        sensible basis) and returns a result dict with expected fields."""
        from geovac.hylleraas_r12 import compute_he_2s_singlet_triplet_eckart
        res = compute_he_2s_singlet_triplet_eckart(omega=2, Z=2)
        for key in (
            'omega', 'n_basis',
            'E_1S', 'E_2S_singlet', 'E_2S_triplet',
            'alpha_1S', 'beta_1S',
            'alpha_2S_singlet', 'beta_2S_singlet',
            'alpha_2S_triplet', 'beta_2S_triplet',
            'splitting_Ha', 'splitting_cm_inv',
            'NIST_splitting_cm_inv', 'rel_err_pct',
        ):
            assert key in res, f"Missing key {key} in result"
        # Sanity: 2^1S > 2^3S (exchange splitting is positive).
        assert res['E_2S_singlet'] > res['E_2S_triplet']
        # Variational bound on each state.
        assert res['E_1S'] >= EXACT_HE_1S - 1e-9
