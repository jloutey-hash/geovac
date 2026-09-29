"""Paper 53, Sec. 'Why the Bochner-Riesz map is not the lift' (subsec:not_br).

The 2D Bochner-Riesz (Cesaro) mean on the flat plane,

    B_Lambda(f) = F^{-1}[(1 - |xi|^2/Lambda)_+^s  fhat],   kernel K_s(r) ~ J_{s+1}(r)/r^{s+1},

was the reconstruction map of an earlier version of Paper 53, which claimed it is
positivity-preserving ("classical, s >= 1/2").  That is FALSE: the kernel
oscillates, so it is negative somewhere at every finite order, and a
non-negative observable concentrated where K < 0 is mapped to a function with
negative values.  The earlier "positivity test" checked one everywhere-positive
Gaussian, which cannot detect this.  (Registered as the C16 entry
p53-bochner-riesz-positivity-and-apex-backbone, 2026-09-29.)

What this file asserts:

  (1) test_bochner_riesz_kernel_is_negative_somewhere -- min K / max K is
      strictly negative at s in {0.75, 1, 2, 4}, and matches the paper's printed
      values (-0.071, -0.059, -0.029, -0.009).  Wrong answers rejected: "the BR
      kernel is non-negative (positivity-preserving) for s >= 1/2", and a drift
      of the printed ratios.

  (2) test_plane_bochner_riesz_converges -- rem:br_rate: on individual smooth
      (Gaussian) observables the BR mean does reconstruct, at an
      observable-dependent algebraic rate (~Lambda^-0.9 and ~Lambda^-0.6).  This
      is NOT the uniform Lipschitz-ball rate the Gromov-Hausdorff bound needs
      (that is the squared-window first moment, see
      test_paper53_flat_plane_squared_window.py).  Wrong answer rejected: "the
      BR error on a smooth observable does not decay".

The construction is self-contained (the sprint driver lived under debug/).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import j0, jv


def _bochner_riesz(f, fhat, rho: np.ndarray, lam: float, s: float) -> np.ndarray:
    """Plane 2D Cesaro/Bochner-Riesz band-limit of a radial function f."""
    xc = np.sqrt(lam)
    out = np.zeros_like(rho)
    for i, r in enumerate(rho):
        integ = lambda xi: (1.0 - xi**2 / lam) ** s * fhat(xi) * j0(xi * r) * xi
        val, _ = quad(integ, 0.0, xc, limit=250)
        out[i] = val
    return out


def bochner_riesz_kernel_ratio(s: float, rmax: float = 60.0, n: int = 200001) -> float:
    """min K_s / max K_s for the 2D BR kernel K_s(r) = J_{s+1}(r) / r^{s+1}.

    K_s is maximal at r -> 0 (value 1/(2^{s+1} Gamma(s+2))); its first negative
    lobe lies inside r < 60 for every s tested, so a dense grid resolves it.
    """
    r = np.linspace(1e-6, rmax, n)
    K = jv(s + 1.0, r) / r ** (s + 1.0)
    return float(K.min() / K.max())


_PRINTED_RATIOS = {0.75: -0.071, 1.0: -0.059, 2.0: -0.029, 4.0: -0.009}


@pytest.mark.parametrize("s", sorted(_PRINTED_RATIOS))
def test_bochner_riesz_kernel_is_negative_somewhere(s: float) -> None:
    """The BR kernel takes negative values, so the BR map is not positive."""
    ratio = bochner_riesz_kernel_ratio(s)
    assert ratio < -1e-3, (
        f"s={s}: BR kernel min/max = {ratio:.5f}; a non-negative kernel would make "
        "the map positivity-preserving, which the paper says it is NOT"
    )
    printed = _PRINTED_RATIOS[s]
    assert abs(ratio - printed) < 0.0015, (
        f"s={s}: min/max {ratio:.5f} drifted from the printed {printed}"
    )


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
_S = 2.0  # Bochner-Riesz order used in rem:br_rate


@pytest.mark.parametrize("name", list(_TESTS))
def test_plane_bochner_riesz_converges(name: str) -> None:
    """rem:br_rate -- per-observable BR reconstruction error decays, p > 0.3.

    Observable-specific rate only; says nothing about positivity (false, see
    above) or about a uniform Lipschitz-ball rate.
    """
    d = _TESTS[name]
    f, fhat = d["f"], d["fhat"]
    fr = f(_RHO)
    gradnum = float(np.max(np.abs(np.gradient(fr, _RHO))))

    errs = []
    for lam in _LAMBDAS:
        g = _bochner_riesz(f, fhat, _RHO, lam, _S)
        errs.append(float(np.max(np.abs(g - fr)) / gradnum))

    slope = float(np.polyfit(np.log(_LAMBDAS), np.log(errs), 1)[0])
    p = -slope
    assert p > 0.3, (
        f"{name}: BR reconstruction of a smooth observable must decay (p > 0.3); got "
        f"Lambda^{{-{p:.3f}}} (errs={[round(e, 4) for e in errs]})"
    )
    assert errs[-1] < errs[0], f"{name}: error did not decrease over the Lambda sweep"
