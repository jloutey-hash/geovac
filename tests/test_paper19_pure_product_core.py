r"""Paper 19, well-shape subsection, "Refinement: the atom-centered core is a
convenience, not a necessity" (Table tab:pure_product_core).

Claim backed: a PURE prolate product basis -- xi-Laguerre x eta-Legendre, one
shared exponent alpha, no atom-centered function (the Paper 11/12 basis) --
carries the Li-like 1s core on the TWO-BODY side.  For Z = 3 with the H nucleus
absent, the Galerkin 1s is exactly e^{-Z r_A} in the complete-basis limit and its
self-repulsion (cc|cc) is exactly 5Z/8 = 15/8.  At alpha = R/2 the 15x15 basis
puts (cc|cc) within 0.01 mHa of 15/8 (measured +6.2e-3 mHa at R = 3.015), and
its R-spread over R = 2.5 / 3.015 / 3.5 below 1e-3 mHa (measured 8.5e-4).

Wrong answers this rejects (fire-tested with debug/qa/fire_test.py):
  * the 10x10 value (+0.362 mHa): the basis is NOT trivially exact, so a guard
    that passes at 10x10 would be checking nothing;
  * too low a working precision.  At dps = 62 (instead of 40 + 8N = 160) the
    same computation returns +5.75e-3 mHa -- inside the 0.01 mHa gate, finite and
    plausible, but wrong.  This is Paper 19's precision trap in miniature, and
    only the pinned high-precision value catches it.

The two-body evaluation is the Paper 12 Neumann route built from geovac
primitives only (closed-form A moments and B seeds of
geovac.neumann_vee_general_m, the Laguerre/Legendre coefficients and eta moments
of geovac.prolate_recondition).  The l-sum terminates at l = 2N because the pair
densities are polynomial in eta.  The routines (core_orbital, conv2, eri_rho) are
the production ones in geovac/prolate_allelectron/pure_product_core_eri.py (promoted
from the sprint driver, LiH productionization step 2a, 2026-10-04; the test previously
carried a private copy so as not to depend on the prunable debug/ tree).

Transcendentals (tagged): the truncated value is a finite Q-combination of the
Paper 12 Neumann seeds e^{-c}, e^{c}E_1(2c), Euler gamma, ln c, ln 2 at the density
rate c = 2 alpha; they cancel to the rational 15/8 only in the complete-basis limit.
"""
from __future__ import annotations

from typing import Tuple

import mpmath as mp
import pytest

from geovac.prolate_allelectron import pure_product_core_eri as ppc

N_GATE = 15
R_GATE = "3.015"
# Converged (cc|cc) - 15/8 at N = 15, alpha = R/2, R = 3.015, in mHa.  Measured
# 2026-09-26 at dps 80, 100, 160 (identical to 17 significant digits).
ERR_15_MHA = 6.2219770788e-3
PIN_TOL_MHA = 1e-6


def _dps(N: int) -> int:
    """Working precision: ~8N digits (Paper 19's precision rule)."""
    return 40 + 8 * N


def cc_error_mha(N: int, R: str, dps: int) -> Tuple[float, float]:
    """((cc|cc) - 15/8 in mHa, 1s energy error in Ha) for Z = 3, alpha = R/2."""
    old = mp.mp.dps
    try:
        mp.mp.dps = dps
        ppc._YCACHE.clear()
        Rm = mp.mpf(R)
        alpha = Rm / 2
        Z = mp.mpf(3)
        K, E, _cond = ppc.core_orbital(N, N, alpha, Rm, Z)
        rho = ppc.conv2(K, K)
        J = ppc.eri_rho(rho, 2 * alpha, rho, 2 * alpha, Rm)
        return float((J - mp.mpf(15) / 8) * 1000), E + 4.5
    finally:
        mp.mp.dps = old
        ppc._YCACHE.clear()


# --------------------------------------------------------------- tests
@pytest.fixture(scope="module")
def err15() -> Tuple[float, float]:
    return cc_error_mha(N_GATE, R_GATE, _dps(N_GATE))


def test_one_body_core_is_resolved(err15):
    """The Galerkin 1s energy is -Z^2/2 = -4.5 to well below the two-body error."""
    _, e1 = err15
    assert 0.0 <= e1 < 1e-6


def test_cc_within_gate_at_15x15(err15):
    """(cc|cc) - 15/8 < 0.01 mHa at 15x15, alpha = R/2, R = 3.015 (Paper 19)."""
    err, _ = err15
    assert abs(err) < 0.01


def test_cc_matches_converged_value(err15):
    """Pinned to the high-precision value; rejects a too-low working precision."""
    err, _ = err15
    assert abs(err - ERR_15_MHA) < PIN_TOL_MHA


def test_smaller_basis_is_not_converged():
    """The gate discriminates: 10x10 is off by +0.36 mHa (Table tab:pure_product_core)."""
    err, _ = cc_error_mha(10, R_GATE, _dps(10))
    assert 0.3 < err < 0.4


@pytest.mark.slow
def test_r_spread_at_15x15():
    """alpha proportional to R removes the R-dependence: spread < 1e-3 mHa at 15x15."""
    errs = [cc_error_mha(N_GATE, R, _dps(N_GATE))[0] for R in ("2.5", "3.015", "3.5")]
    assert all(abs(e) < 0.01 for e in errs)
    assert max(errs) - min(errs) < 1e-3
