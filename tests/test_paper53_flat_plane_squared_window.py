"""Paper 53, Sec. 'The flat plane' -- the squared-window lift (thm:plane_gh).

Window: psi_1 = F^{-1} phi with phi ~ (1 - |xi|^2)_+^m supported in the unit
ball, psi_Lambda(x) = Lambda^{1/2} psi_1(Lambda^{1/2} x).  The lift averages
against the kernel k_Lambda = |psi_Lambda|^2, and the Gromov-Hausdorff bound is
its first moment gamma_Lambda = int |x| k_Lambda(x) dx = gamma_1 Lambda^{-1/2}.

Each test names the wrong answer it rejects:

  test_kernel_is_nonnegative_and_normalised
      k = c_m J_{m+1}(r)^2 / r^{2m+2} (prop:windows) -- checked against a direct
      numerical Hankel transform of (1 - rho^2)^m, so the Sonine closed form is
      not trusted; k >= 0; int k = 1.  Rejects: a wrong c_m / wrong Bessel order.
  test_first_moment_closed_form
      gamma_1 = 32/(5 pi) for m = 1; 2.2995 (m = 2); 2.5934 (m = 3).
      Rejects: a drift of the printed constants.
  test_first_moment_scales_as_inverse_sqrt_lambda
      the fitted exponent of gamma_Lambda over Lambda in [1, 256] is -1/2.
      Rejects: "gamma ~ Lambda^{-1}" (the band-limit scale is sqrt(Lambda)).
  test_fact_iii_on_lipschitz_observables
      ||f * k - f||_inf <= gamma_Lambda L(f) for f = |x - x0| (where it is an
      EQUALITY at x0, so the bound is sharp) and for f = sin(x_1).
      Rejects: a bound with a smaller constant (e.g. gamma/2), and a
      mis-scaled gamma.
  test_sharp_window_first_moment_diverges_logarithmically
      m = 0: partial moment 2 int_0^R J_1^2 grows by (2/pi) log 2 per doubling
      of R.  Rejects: "the sharp window has a finite first moment".
"""

from __future__ import annotations

from math import factorial, pi

import numpy as np
import pytest
from scipy import integrate, special


def c_m(m: int) -> float:
    return 4 ** m * factorial(m) ** 2 * (2 * m + 1) / pi


def k1(r, m: int):
    """Squared window kernel at Lambda = 1 (closed form, prop:windows)."""
    r = np.asarray(r, dtype=float)
    return c_m(m) * special.jv(m + 1, r) ** 2 / r ** (2 * m + 2)


def _radial_integral(fun, rmax: float, step: float = 1.0) -> float:
    edges = np.arange(0.0, rmax + step, step)
    edges[0] = 1e-9
    return float(sum(integrate.quad(fun, a, b, limit=100)[0]
                     for a, b in zip(edges[:-1], edges[1:])))


def gamma_lambda(m: int, lam: float, rmax_scaled: float = 2000.0) -> float:
    """int |x| Lambda k1(sqrt(Lambda) x) d^2x, integrated in x (not substituted)."""
    sl = np.sqrt(lam)
    fun = lambda r: r * lam * k1(sl * r, m) * 2 * pi * r
    return _radial_integral(fun, rmax_scaled / sl, step=1.0 / sl)


@pytest.mark.parametrize("m", [1, 2, 3])
def test_kernel_is_nonnegative_and_normalised(m: int) -> None:
    # independent route: unitary 2D Fourier transform of the normalised window
    raw = lambda rho: (1.0 - rho ** 2) ** m
    nrm2, _ = integrate.quad(lambda rho: raw(rho) ** 2 * 2 * pi * rho, 0.0, 1.0)
    phi = lambda rho: raw(rho) / np.sqrt(nrm2)
    for r in (0.3, 1.7, 4.2, 9.9, 23.5):
        psi, _ = integrate.quad(lambda rho: phi(rho) * special.j0(r * rho) * rho, 0.0, 1.0,
                                limit=200)
        assert abs(psi ** 2 - k1(r, m)) < 1e-10 + 1e-8 * k1(r, m), (
            f"m={m} r={r}: |psi|^2={psi**2:.3e} vs closed form {k1(r, m):.3e}"
        )
    r = np.linspace(1e-6, 80.0, 40001)
    assert np.all(k1(r, m) >= 0.0)
    total = _radial_integral(lambda t: k1(t, m) * 2 * pi * t, 400.0)
    assert abs(total - 1.0) < 1e-6, f"m={m}: int k = {total}"


@pytest.mark.parametrize("m,expected", [(1, 32.0 / (5.0 * pi)), (2, 2.29947), (3, 2.59341)])
def test_first_moment_closed_form(m: int, expected: float) -> None:
    g = gamma_lambda(m, 1.0)
    assert abs(g - expected) < 2e-5, f"m={m}: gamma_1 = {g:.6f}, expected {expected:.6f}"


def test_first_moment_scales_as_inverse_sqrt_lambda() -> None:
    lams = np.array([1.0, 4.0, 16.0, 64.0, 256.0])
    gs = np.array([gamma_lambda(2, lam) for lam in lams])
    slope = float(np.polyfit(np.log(lams), np.log(gs), 1)[0])
    assert abs(slope + 0.5) < 1e-3, f"gamma_Lambda exponent {slope:.5f}, expected -0.5"


def _conv_minus_f_at(y: float, m: int, lam: float) -> float:
    """(f * k)(x) - f(x) for f = |x - x0|, at |x - x0| = y (2D polar quadrature)."""
    sl = np.sqrt(lam)
    kk = lambda r: lam * k1(sl * r, m)

    def inner(r):
        ang = lambda th: np.sqrt(y * y + r * r + 2 * y * r * np.cos(th)) - y
        val, _ = integrate.quad(ang, 0.0, 2 * pi, limit=200)
        return val * kk(r) * r

    return _radial_integral(inner, 2000.0 / sl, step=1.0 / sl)


def test_fact_iii_on_lipschitz_observables() -> None:
    m, lam = 2, 4.0
    gam = gamma_lambda(m, lam)
    # f = |x - x0|, L(f) = 1: the defect is maximal at x0, where it EQUALS gamma
    defects = [_conv_minus_f_at(y, m, lam) for y in (0.0, 0.25, 1.0, 3.0)]
    assert abs(defects[0] - gam) < 1e-5 * gam, f"defect at x0 {defects[0]} vs gamma {gam}"
    assert max(defects) <= gam * (1 + 1e-6), f"defects {defects} exceed gamma {gam}"
    assert defects[0] > 0.9 * gam  # sharp: a bound gamma/2 would fail
    # f = sin(x_1), L(f) = 1: sup |f*k - f| = 1 - int cos(z_1) k = 1 - int J_0(r) k 2 pi r
    sl = np.sqrt(lam)
    ft = _radial_integral(lambda r: special.j0(r) * lam * k1(sl * r, m) * 2 * pi * r,
                          2000.0 / sl, step=1.0 / sl)
    defect_sin = 1.0 - ft
    assert 0.0 < defect_sin <= gam, f"sin defect {defect_sin} vs gamma {gam}"


def test_sharp_window_first_moment_diverges_logarithmically() -> None:
    # m = 0: k = J_1^2/(pi r^2); partial moment int_{|x|<R} |x| k = 2 int_0^R J_1^2
    Rs = [50.0, 100.0, 200.0, 400.0, 800.0]
    M = [_radial_integral(lambda t: 2.0 * special.jv(1, t) ** 2, R) for R in Rs]
    increments = np.diff(M) / np.log(2.0)
    assert np.all(np.abs(increments - 2.0 / pi) < 0.01), (
        f"per-log slopes {increments}, expected 2/pi = {2 / pi:.4f} (log divergence)"
    )
    assert M[-1] - M[0] > 1.7  # grows without bound, not converging
