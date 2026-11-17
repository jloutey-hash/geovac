"""Paper 63, tree fibre at l > 0: engine for the chain of three with p_z leaves (record P63L1, PM to pre-register).

OBJECT (all normalisations stated; conventions of prop:tree_reduction).
  Fingerprint s-type density of scale a:   rho_a(r) = a^3 K0(a r) / (2 pi^2)        (int rho_a d^3r = 1)
  p_z leaf = centre derivative (Lemma lem:orbital, l = 1):
      rho_p(r_vec) = -d/dA_z rho_1(|r - A|) |_{A=0} = d/dz rho_1(r) = f_p(r) cos(theta),
      f_p(r) = -K1(r) / (2 pi^2)                       (dipole moment int z rho_p = -1)
  Chain of three on the slice (t; 1, 1): root (middle) density rho_t, s-type; both leaves rho_p at scale 1;
  Coulomb kernels; coincident centres (the fibre).
      FIBRE(t) = int rho_t(r0) rho_p(r1) rho_p(r2) / (|r0 - r1| |r0 - r2|) d^3r0 d^3r1 d^3r2.

  L = 1 leaf potential (Laplace multipole; replaces the shell theorem):
      V_p(r_vec) = (4 pi / 3) [ r^-2 int_0^r s^3 f_p(s) ds + r int_r^oo f_p(s) ds ] cos(theta) = g(r) cos(theta)
  Closed form (= d/dz of the s-leaf potential (2/pi) chi(r)/r, chi = int_0^r s K1):
      g(r) = (2/pi) [ K1(r) - chi(r) / r^2 ]  = -1/r^2 + E(r),     E(r) = (2/pi) [ K1(r) + R(r)/r^2 ],
      R(r) = int_r^oo s K1 = pi/2 - chi(r).
  Gaunt at the root: cos^2 averages to 1/3 (the L = 0 part of Y10 Y10; the L = 2 part is killed by the s root).
      FIBRE(t) = (4 pi/3) int rho_t g^2 r^2 dr = t^3 int_0^oo K0(t r) H(r) dr,   H(r) = (2/(3 pi)) r^2 g(r)^2.

  The s-type control (prop:tree_reduction, n = 3): CHAIN_s(t) = (2/pi)^3 t^3 int K0(t x) chi(x)^2 dx.

BRANCH AT t = -2 (as record SA): at t = -(2 - eps) + i0, K0(-y + i0) = K0(y) - i pi I0(y); only the piece of H
decaying as e^(-2r) diverges there.  For l = 1 it is H_sing = (2/(3 pi)) r^2 E(r)^2; for s it is (2/pi)^3 R^2.
The scored quantity is the branch integral
      S(eps) = int_1^oo I0((2 - eps) x) H_sing(x) dx      (Im FIBRE = pi (2-eps)^3 S + holomorphic)
"""
from __future__ import annotations

import math
from typing import Callable, Dict, Tuple

import mpmath as mp
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import i0e, k0, k1, k1e, roots_laguerre

# ------------------------------------------------------------------ float64 log grid (as the tree guard)
def make_grid(npts: int = 200001, rmin: float = 1e-10, rmax: float = 90.0) -> Tuple[np.ndarray, np.ndarray]:
    U = np.linspace(np.log(rmin), np.log(rmax), npts)
    return U, np.exp(U)


class Grid:
    def __init__(self, npts: int = 200001, rmin: float = 1e-10, rmax: float = 90.0) -> None:
        self.U, self.R = make_grid(npts, rmin, rmax)

    def cum0(self, f: np.ndarray) -> np.ndarray:
        """int_0^r f(s) ds."""
        return cumulative_trapezoid(f * self.R, self.U, initial=0.0)

    def tail(self, f: np.ndarray) -> np.ndarray:
        """int_r^oo f(s) ds."""
        c = self.cum0(f)
        return c[-1] - c

    def integral(self, f: np.ndarray) -> float:
        return float(np.trapezoid(f * self.R, self.U))

    # densities
    def rho_s(self, a: float) -> np.ndarray:
        return a**3 * k0(a * self.R) / (2 * np.pi**2)

    def f_p(self, a: float = 1.0) -> np.ndarray:
        """radial factor of d/dz rho_a: -a^4 K1(a r)/(2 pi^2)."""
        return -(a**4) * k1(a * self.R) / (2 * np.pi**2)

    # Coulomb multipole potentials: density f(r) Y-type with angular factor P_L(cos) -> potential radial factor
    def phi0(self, W: np.ndarray) -> np.ndarray:
        R = self.R
        return 4 * np.pi * (self.cum0(W * R**2) / R + self.tail(W * R))

    def phi1(self, W: np.ndarray) -> np.ndarray:
        """density W(r) cos(theta) -> potential radial factor (times cos(theta)), Laplace L = 1."""
        R = self.R
        return 4 * np.pi / 3 * (self.cum0(W * R**3) / R**2 + R * self.tail(W))

    def chi(self) -> np.ndarray:
        return self.cum0(self.R * k1(self.R))


def g_closed_grid(G: Grid) -> np.ndarray:
    chi = G.chi()
    return (2 / np.pi) * (k1(G.R) - chi / G.R**2)


def fibre_l1_root(t: float, G: Grid) -> float:
    """rooted at the s root: (4 pi/3) int rho_t V_p^2 r^2 dr with V_p from the multipole recursion."""
    Vp = G.phi1(G.f_p(1.0))
    return 4 * np.pi / 3 * G.integral(G.rho_s(t) * Vp**2 * G.R**2)


def fibre_l1_leaf(t: float, G: Grid) -> float:
    """rooted at leaf 1: leaf 2 -> V_p (L=1); times rho_t -> an L=1 density; its L=1 potential against f_p."""
    Vp = G.phi1(G.f_p(1.0))
    W = G.rho_s(t) * Vp                      # radial factor of an L = 1 density (cos theta)
    U1 = G.phi1(W)
    return 4 * np.pi / 3 * G.integral(G.f_p(1.0) * U1 * G.R**2)


def fibre_l1_closed(t: float, G: Grid) -> float:
    g = g_closed_grid(G)
    return t**3 * G.integral(k0(t * G.R) * (2 / (3 * np.pi)) * G.R**2 * g**2)


def chain_s_recursion(t: float, G: Grid, root: int = 0) -> float:
    """s-type chain (t;1,1) by the tree recursion; root 0 = middle, root 1 = a leaf."""
    if root == 0:
        V = G.phi0(G.rho_s(1.0))
        return 4 * np.pi * G.integral(G.rho_s(t) * V**2 * G.R**2)
    V = G.phi0(G.rho_s(1.0))
    W = G.phi0(G.rho_s(t) * V)
    return 4 * np.pi * G.integral(G.rho_s(1.0) * W * G.R**2)


def chain_s_closed(t: float, G: Grid) -> float:
    chi = G.chi()
    return 8 / np.pi**3 * t**3 * G.integral(chi**2 * k0(t * G.R))


# ------------------------------------------------------------------ mpmath (high precision) building blocks
def chi_mp(r) -> mp.mpf:
    return mp.quad(lambda s: s * mp.besselk(1, s), [0, r])


def Rscaled_quad_mp(x) -> mp.mpf:
    """e^x R(x) = int_0^oo (x+u) K1(x+u) e^(x+u) e^-u du (direct quadrature)."""
    x = mp.mpf(x)
    return mp.quad(lambda u: (x + u) * mp.besselk(1, x + u) * mp.exp(x + u) * mp.exp(-u), [0, 1, 10, mp.inf])


def _R_series_coeffs(M: int = 400):
    """e^x R(x) ~ sqrt(pi/2) x^(1/2) sum_m d_m x^-m, d_m = sum_{k+j=m} a_k(K1) ff(1/2 - k, j) (exact)."""
    from fractions import Fraction
    a = [Fraction(1)]
    for k in range(1, M):
        a.append(a[-1] * Fraction(4 - (2 * k - 1) ** 2, 8 * k))
    d = [Fraction(0)] * M
    for k in range(M):
        b = Fraction(1, 2) - k
        ff = Fraction(1)
        for j in range(M - k):
            d[k + j] += a[k] * ff
            ff *= (b - j)
    return d


_DSER: list = []
_RCACHE: Dict = {}
X_SERIES = 80


def Rscaled_mp(x) -> mp.mpf:
    """e^x R(x): asymptotic series for x >= 80 (optimal-truncation error ~ e^-x < 1e-34 relative),
    direct quadrature below (cached)."""
    x = mp.mpf(x)
    if x >= X_SERIES:
        if not _DSER:
            _DSER.extend(_R_series_coeffs())
        tol = mp.mpf(10) ** (-mp.mp.dps - 4)
        s, p, prev = mp.mpf(0), mp.mpf(1), mp.inf
        for dm in _DSER:
            term = mp.mpf(dm.numerator) / dm.denominator * p
            if abs(term) > prev:
                raise RuntimeError("asymptotic series diverged before reaching tolerance")
            s += term
            if abs(term) < tol and dm != 0:
                break
            prev = abs(term) if dm != 0 else prev
            p /= x
        return mp.sqrt(mp.pi / 2) * mp.sqrt(x) * s
    return Rscaled_struve_mp(x)


def Rscaled_struve_mp(x) -> mp.mpf:
    """e^x R(x), R = x K0 + pi/2 - int_0^x K0, int_0^x K0 = (pi x/2)[K0 L_{-1} + K1 L_0] (modified Struve);
    the cancellation loses ~x/ln 10 digits, restored by raising the working precision.  Checked against
    Rscaled_quad_mp to ~1e-31 at x = 1, 3.3, 50.3, 79.9 (dps 30)."""
    dps = mp.mp.dps
    with mp.workdps(dps + int(mp.mpf(x) / 2.3) + 15):
        x = mp.mpf(x)
        v = x * mp.besselk(0, x) + mp.pi / 2 - mp.pi * x / 2 * (mp.besselk(0, x) * mp.struvel(-1, x)
                                                              + mp.besselk(1, x) * mp.struvel(0, x))
        v = v * mp.exp(x)
    return +v


def g_closed_mp(r) -> mp.mpf:
    r = mp.mpf(r)
    return 2 / mp.pi * (mp.besselk(1, r) - chi_mp(r) / r**2)


def Escaled_mp(x) -> mp.mpf:
    """e^x E(x) = (2/pi) [e^x K1(x) + e^x R(x)/x^2]."""
    x = mp.mpf(x)
    return 2 / mp.pi * (mp.besselk(1, x) * mp.exp(x) + Rscaled_mp(x) / x**2)


# ------------------------------------------------------------------ branch integrands (scaled: e^{-2x} removed)
# h(x) = e^{2x} H_sing(x); branch integrand = I0((2-eps)x) e^{-(2-eps)x} * h(x) * e^{-eps x}
_UL, _WL = roots_laguerre(80)


def Rscaled_f(x: np.ndarray) -> np.ndarray:
    X = x[:, None] + _UL[None, :]
    return (X * k1e(X)) @ _WL


def h_l1_f(x: np.ndarray) -> np.ndarray:
    Es = 2 / np.pi * (k1e(x) + Rscaled_f(x) / x**2)
    return 2 / (3 * np.pi) * x**2 * Es**2


def h_s_f(x: np.ndarray) -> np.ndarray:
    return (2 / np.pi) ** 3 * Rscaled_f(x) ** 2


def h_l1_mp(x) -> mp.mpf:
    return 2 / (3 * mp.pi) * x**2 * Escaled_mp(x) ** 2


def h_s_mp(x) -> mp.mpf:
    return (2 / mp.pi) ** 3 * Rscaled_mp(x) ** 2


def branch_float(eps: float, h: Callable[[np.ndarray], np.ndarray], npts: int = 400001) -> float:
    """S(eps) = int_1^oo I0((2-eps)x) H_sing dx, float64, log grid (SA's scheme)."""
    y = np.linspace(0.0, np.log(1.0 + 80.0 / eps), npts)
    x = np.exp(y)
    f = i0e((2 - eps) * x) * h(x) * np.exp(-eps * x)
    return float(np.trapezoid(f * x, y))


def _gl_nodes(n: int) -> Tuple[list, list]:
    xs, ws = [], []
    # mpmath Gauss-Legendre on [-1, 1] via legendre roots
    def pdp(z):
        p0, p1 = mp.mpf(1), z
        for j in range(2, n + 1):
            p0, p1 = p1, ((2 * j - 1) * z * p1 - (j - 1) * p0) / j
        return p1, n * (z * p1 - p0) / (z * z - 1)

    for k in range(1, n + 1):
        x = mp.cos(mp.pi * (k - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(100):
            p, dp = pdp(x)
            dx = p / dp
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        p, dp = pdp(x)
        xs.append(x)
        ws.append(2 / ((1 - x**2) * dp**2))
    return xs, ws


_GL: Dict = {}


def branch_mp(eps, h: Callable, ngl: int = 40, cache: Dict | None = None) -> mp.mpf:
    """S(eps) at mp.dps, composite Gauss-Legendre on [1, 100/eps]: geometric panels (ratio 2) until the panel
    length reaches 2/eps, then uniform panels of length 2/eps.  h values cached by node."""
    eps = mp.mpf(eps)
    key_gl = (ngl, mp.mp.dps)
    if key_gl not in _GL:
        _GL[key_gl] = _gl_nodes(ngl)
    xs, ws = _GL[key_gl]
    edges = [mp.mpf(1)]
    X = 100 / eps
    while edges[-1] < X:
        a = edges[-1]
        edges.append(min(a + min(a, 2 / eps), X))
    tot = mp.mpf(0)
    a2 = 2 - eps
    for a, b in zip(edges[:-1], edges[1:]):
        hm, hw = (a + b) / 2, (b - a) / 2
        for xi, wi in zip(xs, ws):
            x = hm + hw * xi
            key = None
            if cache is not None:
                key = mp.nstr(x, mp.mp.dps)
                hv = cache.get(key)
                if hv is None:
                    hv = h(x)
                    cache[key] = hv
            else:
                hv = h(x)
            tot += wi * hw * mp.besseli(0, a2 * x) * mp.exp(-a2 * x) * hv * mp.exp(-eps * x)
    return tot
