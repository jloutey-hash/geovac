"""Paper 63, Proposition prop:tree_rung_l: the rung bound for l > 0 densities (degree <= 2 in the centre derivative).

The proof changes three steps of prop:tree_rung; this file checks the identities those steps use and the
end-to-end reduction on one chain.  Each test names the wrong answer it rejects.

(a) channel-L Poisson: Phi_L W = 4 pi/(2L+1) int r<^L r>^(-L-1) W(s) s^2 ds obeys
    (theta - L)(theta + L + 1) Phi_L W = -4 pi r^2 W, theta = r d/dr, for L = 1, 2
    (rejects a wrong constant 4 pi/(2L+1), a wrong L(L+1), and an r^(L+1) / r^(-L) kernel);
(b) Hobson for r^2 P2 on K0 at an off-axis point, and the d twist b(rho) = K0 + 2 K1/rho = rho^-2 (theta-2) theta K0
    (rejects (theta-1) theta, a wrong P2, and K0 + K1/rho);
(c) d-root transfer: int b(t r) g(r) dr = t^-2 (theta_t - 2) theta_t T(r^-2 g), T(h) = int K0(t r) h(r) dr
    (rejects (theta_t - 1), t^-1, and T(r^-1 g));
(d) the p_z leaf - p_z middle - d root chain: the radial recursion rooted at either end agrees, through different
    Gaunt coefficients (rejects cos P2 -> (3/5) P1 and a dropped 4 pi/5 P2 normalisation); the common value is pinned
    to the one a 4e7-sample Monte Carlo of the full 9-D integral agreed with (rejects common-mode errors, which the
    root swap cannot see); the Monte Carlo itself, small sample, is @slow;
(e) the twist bookkeeping: the l-type radial part r^l (r^-1 d/dr)^l K0 = r^-l prod_{j<l} (theta - 2j) K0
    = (-1)^l K_l, so r^2 x radial part is O(1 + |log r|) at 0 for l = 0, 1, 2 and the Laplacian, and
    grows like r^-1 for l = 3 -- the reason the hypothesis is degree <= 2.

Backing scripts (not imported): debug/p63rb_identities.py, debug/p63rb_chain_check.py.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest
from scipy.integrate import cumulative_simpson, quad, simpson
from scipy.special import k0, k0e, k1, k1e, kv

# ----------------------------------------------------------------------------------------------- (a) Poisson
POISSON_CONST_NUM = 4          # 4 pi / (2L + 1)
POISSON_CONST_DEN_SHIFT = 1
KERNEL_IN_POW = 0              # r<^(L + KERNEL_IN_POW) r>^(-L - 1)


def _quad(fn, a, b):
    return quad(fn, a, b, epsabs=0.0, epsrel=1e-13, limit=200)[0]


def _phi_L(f, L: int, r: float) -> float:
    Lp = L + KERNEL_IN_POW
    inner = _quad(lambda s: s ** Lp / r ** (L + 1) * f(s) * s ** 2, 0.0, r)
    outer = _quad(lambda s: r ** Lp / s ** (L + 1) * f(s) * s ** 2, r, 2 * r) +         _quad(lambda s: r ** Lp / s ** (L + 1) * f(s) * s ** 2, 2 * r, np.inf)
    return POISSON_CONST_NUM * np.pi / (2 * L + POISSON_CONST_DEN_SHIFT) * (inner + outer)


def _d12(F, x, h):
    v = {k: F(x + k * h) for k in (-2, -1, 0, 1, 2)}
    d1 = (v[-2] - 8 * v[-1] + 8 * v[1] - v[2]) / (12 * h)
    d2 = (-v[-2] + 16 * v[-1] - 30 * v[0] + 16 * v[1] - v[2]) / (12 * h * h)
    return v[0], d1, d2


_SOURCES = {"K0K1": lambda s: k0(s) * k1(s), "K1sq_over_r": lambda s: k1(s) ** 2 / s}


@pytest.mark.parametrize("L", [1, 2])
@pytest.mark.parametrize("src", ["K0K1", "K1sq_over_r"])
def test_channel_L_poisson(L: int, src: str) -> None:
    """float64 quadrature + 5-point finite differences; residual ~1e-8, wrong answers are O(1) off."""
    f = _SOURCES[src]
    for r in (0.7, 2.5):
        P0, d1, d2 = _d12(lambda u: _phi_L(f, L, u), r, 2e-3 * r)
        th = r * d1
        th2 = th + r ** 2 * d2
        lhs = th2 + th - L * (L + 1) * P0          # (theta - L)(theta + L + 1)
        rhs = -4 * np.pi * r ** 2 * f(r)
        assert abs(lhs / rhs - 1) < 1e-6, (L, src, r, lhs / rhs)


# ----------------------------------------------------------------------------------------------- (b) Hobson
def _theta_ops(F, r):
    """theta F and theta^2 F at r by mpmath differentiation."""
    d1 = mp.diff(F, r, 1)
    d2 = mp.diff(F, r, 2)
    return r * d1, r * d1 + r ** 2 * d2


D_TWIST_SHIFT = 2              # (theta - 2) theta for l = 2


def test_hobson_r2P2_on_K0() -> None:
    with mp.workdps(30):
        x, y, z = mp.mpf("0.4"), mp.mpf("-0.9"), mp.mpf("1.1")
        F = lambda X, Y, Z: mp.besselk(0, mp.sqrt(X * X + Y * Y + Z * Z))
        lhs = (mp.diff(F, (x, y, z), (0, 0, 2))
               - (mp.diff(F, (x, y, z), (2, 0, 0)) + mp.diff(F, (x, y, z), (0, 2, 0))) / 2)
        r = mp.sqrt(x * x + y * y + z * z)
        c = z / r
        P2 = (3 * c * c - 1) / 2
        th, th2 = _theta_ops(lambda u: mp.besselk(0, u), r)
        rad = r ** -4 * (th2 - D_TWIST_SHIFT * th)     # (r^-1 d/dr)^2 K0 = r^-4 (theta - 2) theta K0
        rhs = r ** 2 * P2 * rad
        assert abs(lhs / rhs - 1) < 1e-20


def _b_closed(rho):
    return mp.besselk(0, rho) + 2 * mp.besselk(1, rho) / rho


def test_d_twist_closed_form() -> None:
    with mp.workdps(30):
        for rho in (mp.mpf("0.3"), mp.mpf("1.7"), mp.mpf("4.0")):
            th, th2 = _theta_ops(lambda u: mp.besselk(0, u), rho)
            op = rho ** -2 * (th2 - D_TWIST_SHIFT * th)
            assert abs(_b_closed(rho) / op - 1) < 1e-20
            assert abs(_b_closed(rho) / mp.besselk(2, rho) - 1) < 1e-25   # = K_2


# ----------------------------------------------------------------------------------------------- (c) transfer
TRANSFER_T_POW = -2
TRANSFER_G_POW = -2


def test_d_root_transfer() -> None:
    """float64 quadrature + 5-point finite differences in t; residual ~1e-8, wrong answers are O(1) off."""
    g = lambda r: r ** 3 * k1(r) * k0(r)
    b = lambda rho: k0(rho) + 2 * k1(rho) / rho
    for t in (0.6, 2.2):
        lhs = _quad(lambda r: b(t * r) * g(r), 0.0, 1.0) + _quad(lambda r: b(t * r) * g(r), 1.0, np.inf)

        def T(u):
            fn = lambda r: k0(u * r) * g(r) * r ** TRANSFER_G_POW
            return _quad(fn, 0.0, 1.0) + _quad(fn, 1.0, np.inf)
        _, d1, d2 = _d12(T, t, 2e-3 * t)
        th = t * d1
        th2 = th + t ** 2 * d2
        rhs = t ** TRANSFER_T_POW * (th2 - D_TWIST_SHIFT * th)
        assert abs(lhs / rhs - 1) < 1e-6, (t, lhs / rhs)


# ----------------------------------------------------------------------------------------------- (d) chain
_u = np.linspace(np.log(1e-9), np.log(80.0), 200001)
_r = np.exp(_u)
_f = -k1(_r) / (2 * np.pi ** 2)            # p_z density = f(r) cos(theta)

GAUNT_COS_COS = 2 / 3                      # cos^2 = 1/3 + (2/3) P2
GAUNT_COS_P2 = 2 / 5                       # cos P2 = (2/5) P1 + (3/5) P3
NORM_P1 = 4 * np.pi / 3                    # int P1^2 dOmega
NORM_P2 = 4 * np.pi / 5                    # int P2^2 dOmega
ROOT_POW = 5                               # root = t^5/(2 pi^2) K_2(t r) P2


def _cum0(y):
    return cumulative_simpson(y * _r, x=_u, initial=0.0)


def _tail(y):
    zz = (y * _r)[::-1]
    return cumulative_simpson(zz, x=-_u[::-1], initial=0.0)[::-1]


def _phi(L, w):
    return 4 * np.pi / (2 * L + 1) * (_r ** (-(L + 1)) * _cum0(_r ** (L + 2) * w) + _r ** L * _tail(_r ** (1 - L) * w))


def _routes(t: float):
    R3 = t ** ROOT_POW / (2 * np.pi ** 2) * (k0(t * _r) + 2 * k1(t * _r) / (t * _r))
    g1 = _phi(1, _f)
    I3 = NORM_P2 * simpson(_r ** 2 * R3 * _phi(2, GAUNT_COS_COS * _f * g1) * _r, x=_u)
    I1 = NORM_P1 * simpson(_r ** 2 * _f * _phi(1, GAUNT_COS_P2 * _f * _phi(2, R3)) * _r, x=_u)
    return I3, I1


# values of the radial recursion that a 4e7-sample Monte Carlo of the full 9-D integral agreed with
# (z = 0.48, 0.51; debug/data/p63rb_chain_check.json)
_MC_VALIDATED = {0.5: 3.3131903415067546e-4, 1.3: 3.365418149870855e-3}


@pytest.mark.parametrize("t", [0.5, 1.3])
def test_chain_root_swap(t: float) -> None:
    I3, I1 = _routes(t)
    assert abs(I3 / I1 - 1) < 1e-12, (I3, I1)


@pytest.mark.parametrize("t", [0.5, 1.3])
def test_chain_value_pinned_to_mc_validated(t: float) -> None:
    I3, _ = _routes(t)
    assert abs(I3 / _MC_VALIDATED[t] - 1) < 1e-8


def _mc(t: float, n: int, seed: int = 7):
    rng = np.random.default_rng(seed)

    def dirs(k):
        v = rng.normal(size=(k, 3))
        return v / np.linalg.norm(v, axis=1)[:, None]
    a1, a2, a3 = rng.gamma(2.0, size=n), rng.gamma(2.0, size=n), rng.exponential(1.0 / t, size=n)
    d1, d2, d3 = dirs(n), dirs(n), dirs(n)
    x1, x2, x3 = a1[:, None] * d1, a2[:, None] * d2, a3[:, None] * d3
    p1 = np.exp(-a1) / (4 * np.pi * a1)
    p2 = np.exp(-a2) / (4 * np.pi * a2)
    p3 = t * np.exp(-t * a3) / (4 * np.pi * a3 ** 2)
    rho1 = -k1e(a1) * np.exp(-a1) * d1[:, 2] / (2 * np.pi ** 2)
    rho2 = -k1e(a2) * np.exp(-a2) * d2[:, 2] / (2 * np.pi ** 2)
    tr = t * a3
    P2 = 1.5 * d3[:, 2] ** 2 - 0.5
    rho3 = t ** 5 / (2 * np.pi ** 2) * (k0e(tr) + 2 * k1e(tr) / tr) * np.exp(-tr) * P2
    w = rho1 * rho2 * rho3 / (p1 * p2 * p3)
    k12 = 1 / np.linalg.norm(x1 - x2, axis=1)
    val = w * k12 * (1 / np.linalg.norm(x2 - x3, axis=1) + 1 / np.linalg.norm(x2 + x3, axis=1)) / 2
    return float(val.mean()), float(val.std(ddof=1) / np.sqrt(n))


@pytest.mark.slow
def test_chain_monte_carlo_small_sample() -> None:
    """Independent 9-D Monte Carlo, no multipoles; ~10% relative SE at 4e6, so a 4-SE window (weak by design)."""
    t = 1.3
    I3, _ = _routes(t)
    m, se = _mc(t, 4_000_000)
    assert abs(m - I3) < 4 * se, (m, se, I3)


# ----------------------------------------------------------------------------------------------- (e) twists
MAX_DEGREE = 2                     # the proposition's hypothesis


def _radial_by_twist(ell: int, r: mp.mpf) -> mp.mpf:
    """r^-l prod_{j<l} (theta - 2j) K0, theta applied by mpmath differentiation of polynomials in theta."""
    # expand prod_{j<l}(theta - 2j) = sum c_k theta^k and theta^k F = sum_i S(k,i) r^i F^(i)
    poly = [mp.mpf(1)]
    for j in range(ell):
        poly = [(poly[k - 1] if k >= 1 else 0) - 2 * j * (poly[k] if k < len(poly) else 0) for k in range(len(poly) + 1)]
    derivs = [mp.diff(lambda u: mp.besselk(0, u), r, i) for i in range(ell + 1)]
    total = mp.mpf(0)
    for k, ck in enumerate(poly):
        for i in range(k + 1):
            total += ck * mp.stirling2(k, i) * r ** i * derivs[i]
    return r ** -ell * total


@pytest.mark.parametrize("ell", [1, 2, 3])
def test_twist_equals_K_ell(ell: int) -> None:
    with mp.workdps(30):
        for r in (mp.mpf("0.4"), mp.mpf("2.1")):
            assert abs(_radial_by_twist(ell, r) / ((-1) ** ell * mp.besselk(ell, r)) - 1) < 1e-20


def _r2_radial(kind, r):
    if kind == "lap":                      # r^-2 theta(theta + 1) K0 = K0 - K1/r
        return r ** 2 * (k0(r) - k1(r) / r)
    return r ** 2 * kv(kind, r)            # |l-type radial part| = K_l


def _bounded_by_log(kind) -> bool:
    rs = np.array([1e-3, 1e-5, 1e-7, 1e-9])
    ratio = np.abs(_r2_radial(kind, rs)) / (1 + np.abs(np.log(rs)))
    return bool(np.all(ratio < 3.0))


def test_twist_bound_holds_through_degree_two() -> None:
    for kind in [*range(MAX_DEGREE + 1), "lap"]:
        assert _bounded_by_log(kind), kind


def test_twist_bound_fails_at_degree_three() -> None:
    assert not _bounded_by_log(MAX_DEGREE + 1)
    # and the failure is the r^-1 growth of r^2 K_3 ~ 8/r
    rs = np.array([1e-4, 1e-6])
    assert np.allclose(_r2_radial(3, rs) * rs, 8.0, rtol=1e-3)
