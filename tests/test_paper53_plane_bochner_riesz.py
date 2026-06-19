"""Paper 53 backing test: the boundaryless plane R^2_alpha Bochner-Riesz Berezin.

Paper 53's genuine convergent carrier is NOT the finite-radius Dirichlet disk
(whose Markov-Cesaro reconstruction is obstructed by the boundary -- positivity
fails, min g < 0, and the interior approximate-identity error does NOT decay:
the disk driver gives Lambda^{+0.07}, not a convergent rate).  The CLEAN route,
already invoked in the paper as the cure, is the de-compactified cone = the
boundaryless plane R^2 via the 2D Bochner-Riesz / Cesaro band-limit of the
Hankel transform:

    B_Lambda(f)(rho) = int_0^{sqrt(Lambda)} (1 - xi^2/Lambda)^s fhat(xi)
                           J_0(xi rho) xi dxi.

There is no Dirichlet boundary.  This test verifies the two numerically
load-bearing L4 properties on the plane for radial Schwartz observables with
closed-form Hankel transforms:

  (a) positivity   -- the 2D Bochner-Riesz kernel is positive for s >= 1/2
                      (classical), so f >= 0 => B_Lambda(f) >= 0;
  (d) approx id    -- ||B_Lambda(f) - f|| -> 0 at an algebraic rate
                      Lambda^{-p} with p > 0 (genuine decay), proportional to
                      ||grad f|| (band-limiting error ~ tail of fhat beyond
                      sqrt(Lambda)).

This is a GENUINE carrier test: it runs the actual plane construction and
asserts a *converging* rate (p > 0) and positivity, the opposite of what the
finite disk produces.  It replaces the previously-cited S^3 tests
(test_berezin_reconstruction / test_joint_berezin_compact_temporal), which are
wrong-carrier (they test the round-S^3 Berezin map, not the disk/plane).

Provenance: debug/archive/gravity_arc/b3_plane_bochner_riesz.py (the sprint
driver); the construction is self-contained here so the test survives debug/
pruning.  See CHANGELOG and the Paper-53 abstract (boundaryless-plane cure).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import j0


def _bochner_riesz(f, fhat, rho: np.ndarray, lam: float, s: float) -> np.ndarray:
    """Plane 2D Cesaro/Bochner-Riesz band-limit of a radial function f."""
    xc = np.sqrt(lam)
    out = np.zeros_like(rho)
    for i, r in enumerate(rho):
        integ = lambda xi: (1.0 - xi**2 / lam) ** s * fhat(xi) * j0(xi * r) * xi
        val, _ = quad(integ, 0.0, xc, limit=250)
        out[i] = val
    return out


# Radial Schwartz observables with closed-form 2D Hankel transforms.
# Gaussian is Hankel self-dual in 2D: fhat(xi) = exp(-xi^2/2).
_TESTS = {
    "gaussian": dict(
        f=lambda r: np.exp(-(r**2) / 2.0),
        fhat=lambda xi: np.exp(-(xi**2) / 2.0),
    ),
    # sigma = 1/2 : f = exp(-2 r^2), fhat = (1/4) exp(-xi^2/8)
    "gaussian_narrow": dict(
        f=lambda r: np.exp(-2.0 * r**2),
        fhat=lambda xi: 0.25 * np.exp(-(xi**2) / 8.0),
    ),
}

_RHO = np.linspace(0.0, 12.0, 400)
_LAMBDAS = [5.0, 10.0, 20.0, 40.0, 80.0, 160.0]
_S = 2.0  # Bochner-Riesz order (>= 1/2 guarantees classical positivity)


@pytest.mark.slow
@pytest.mark.parametrize("name", list(_TESTS))
def test_plane_bochner_riesz_converges(name: str) -> None:
    """The plane reconstruction error decays at a genuine algebraic rate p > 0.

    This is the load-bearing contrast with the finite disk (Lambda^{+0.07},
    NON-convergent).  We require the fitted exponent p in ||g-f|| ~ Lambda^{-p}
    to be clearly positive (p > 0.3), confirming de-compactified convergence.
    """
    d = _TESTS[name]
    f, fhat = d["f"], d["fhat"]
    fr = f(_RHO)
    gradnum = float(np.max(np.abs(np.gradient(fr, _RHO))))

    errs = []
    for lam in _LAMBDAS:
        g = _bochner_riesz(f, fhat, _RHO, lam, _S)
        errs.append(float(np.max(np.abs(g - fr)) / gradnum))

    # log-log slope; convergence <=> slope < 0 i.e. p = -slope > 0
    slope = float(np.polyfit(np.log(_LAMBDAS), np.log(errs), 1)[0])
    p = -slope
    assert p > 0.3, (
        f"{name}: plane reconstruction must CONVERGE (p > 0.3); got "
        f"Lambda^{{-{p:.3f}}} (errs={[round(e, 4) for e in errs]})"
    )
    # sanity: the error must actually shrink across the sweep
    assert errs[-1] < errs[0], f"{name}: error did not decrease over the Lambda sweep"


@pytest.mark.slow
def test_plane_bochner_riesz_positivity() -> None:
    """f >= 0 => B_Lambda(f) >= 0 (classical Bochner-Riesz positivity, s >= 1/2).

    Checked on the (everywhere-positive) Gaussian; the reconstructed g must stay
    non-negative up to small numerical/quadrature tolerance -- the opposite of
    the finite-disk Markov-Cesaro construction, where min g reaches ~ -0.13.
    """
    d = _TESTS["gaussian"]
    f, fhat = d["f"], d["fhat"]
    worst = 0.0
    for lam in _LAMBDAS:
        g = _bochner_riesz(f, fhat, _RHO, lam, _S)
        worst = min(worst, float(g.min()))
    assert worst >= -1e-2, (
        f"plane Bochner-Riesz (s={_S}) must be positivity-preserving for f>=0; "
        f"min g = {worst:.4f}"
    )
