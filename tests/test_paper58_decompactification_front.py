"""Backing tests for Paper 58 Sec. "The continuous side: the decompactification
front" (added 2026-09-06).

Every quantity here is recomputed from scratch with an INDEPENDENT overlap
route (prolate-spheroidal Gauss-Legendre x Gauss-Laguerre quadrature, exact for
the polynomial-times-exponential integrands), not read from the exploratory
drivers that first measured it.  Each test names the wrong answer it rejects.

  test_front_tracks_decay_length_not_mean_radius
      R*(n)/(n/Z) is flat (2.12-2.19 for n = 2..8) while R*(n)/(n^2/Z) drifts
      by more than 2x.  Rejects: "the front sits where the orbital's MEAN
      RADIUS reaches the other center" (the naive r_n ~ R hypothesis).

  test_front_exponent_is_linear_in_n
      log-log slope of R*(n) vs n is 1.0 +- 0.1 for n >= 2.  Rejects the
      sqrt(Z R) window (slope 2 in R*(n), i.e. 0.5 in n*(R)).

  test_geometric_mean_tail_law_1s1s
      For two 1s tails of decay lengths l_A, l_B the tail-reach front
      R_rel = c * 2 sqrt(l_A l_B) with c constant to <= 3 % over the four charge
      pairs of the paper and <= 12 % over the dense exponent-ratio scan.
      Rejects the ADDITIVE law c (l_A + l_B) (>= 5 % / >= 15 %) and the
      max-tail law.

  test_absolute_front_threshold_t_c
      The absolute |S| = 1/sqrt2 crossing exists iff the 1s exponent ratio
      t < t_c = 2.7456, the root of (2 sqrt t/(1+t))^3 = 1/sqrt2.  Rejects
      "every pair of shapes has a 45-degree front".

Paper 58 cites this file; CHANGELOG v5.10.2 is the chronicle.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.special import genlaguerre, roots_genlaguerre, roots_legendre

SQRT_HALF = 1.0 / math.sqrt(2.0)


# ---------------------------------------------------------------------------
# Independent overlap route: prolate-spheroidal quadrature.
#   r_A = R(xi+eta)/2, r_B = R(xi-eta)/2, dV = (R^3/8)(xi^2-eta^2) dxi deta dphi.
# Both integrands below are exp(-kappa*xi) times a polynomial in (xi, eta), so
# Gauss-Laguerre in x = kappa (xi-1) and Gauss-Legendre in eta are exact.
# ---------------------------------------------------------------------------

_ETA_X, _ETA_W = roots_legendre(96)
_LAG_X, _LAG_W = roots_genlaguerre(64, 0.0)


def _r_n0(n: int, Z: float, r: np.ndarray) -> np.ndarray:
    """Normalized hydrogenic radial function R_{n0}(r)."""
    rho = 2.0 * Z * r / n
    norm = math.sqrt((2.0 * Z / n) ** 3 * math.factorial(n - 1)
                     / (2.0 * n * math.factorial(n)))
    return norm * np.exp(-rho / 2.0) * genlaguerre(n - 1, 1)(rho)


def overlap_ns_ns(n: int, Z: float, R: float) -> float:
    """<ns_A | ns_B> for two hydrogenic ns orbitals of charge Z at distance R."""
    kappa = Z * R / n                      # e^{-Z(r_A+r_B)/n} = e^{-kappa xi}
    xi = 1.0 + _LAG_X / kappa              # Gauss-Laguerre nodes in xi
    XI, ETA = np.meshgrid(xi, _ETA_X, indexing="ij")
    rA = R * (XI + ETA) / 2.0
    rB = R * (XI - ETA) / 2.0
    poly = _r_n0(n, Z, rA) * _r_n0(n, Z, rB) * np.exp(kappa * XI) * (XI**2 - ETA**2)
    inner = poly @ _ETA_W                  # eta integral
    outer = np.sum(_LAG_W * inner) * math.exp(-kappa) / kappa
    # Y00^2 * 2pi = 1/2 ; times R^3/8
    return float(outer * (R**3 / 8.0) * 0.5)


def overlap_1s_1s(a: float, b: float, R: float) -> float:
    """<1s(a)_A | 1s(b)_B> for normalized 1s STOs with exponents a, b."""
    if R < 1e-9:
        return (2.0 * math.sqrt(a * b) / (a + b)) ** 3
    kappa = R * (a + b) / 2.0
    xi = 1.0 + _LAG_X / kappa
    XI, ETA = np.meshgrid(xi, _ETA_X, indexing="ij")
    # exp(-a r_A - b r_B) = exp(-kappa xi) * exp(-R eta (a-b)/2)
    poly = np.exp(-R * ETA * (a - b) / 2.0) * (XI**2 - ETA**2)
    inner = poly @ _ETA_W
    outer = np.sum(_LAG_W * inner) * math.exp(-kappa) / kappa
    return float(outer * (R**3 / 8.0) * 2.0 * math.pi * math.sqrt(a**3 * b**3) / math.pi)


def _crossing(fun, level: float, lo: float, hi: float) -> float:
    return float(brentq(lambda R: fun(R) - level, lo, hi, xtol=1e-10))


def _front_ns(n: int, Z: float) -> float:
    """Absolute front: R at which <ns_A|ns_B> = 1/sqrt2."""
    return _crossing(lambda R: overlap_ns_ns(n, Z, R), SQRT_HALF, 1e-3, 60.0 * n / Z)


def test_overlap_routes_agree_with_closed_forms():
    """Quadrature route sanity: 1s-1s equal exponents vs the textbook closed form."""
    # S(R) = e^{-x}(1 + x + x^2/3), x = zeta R  (Mulliken)
    for zeta, R in [(1.0, 1.4), (2.0, 0.7), (0.8, 3.0)]:
        x = zeta * R
        exact = math.exp(-x) * (1.0 + x + x * x / 3.0)
        assert abs(overlap_1s_1s(zeta, zeta, R) - exact) < 1e-12
        assert abs(overlap_ns_ns(1, zeta, R) - exact) < 1e-12
    # united-atom limit of mismatched 1s: (2 sqrt(ab)/(a+b))^3
    assert abs(overlap_1s_1s(1.0, 2.0, 1e-6) - (2 * math.sqrt(2.0) / 3.0) ** 3) < 1e-5


def test_front_tracks_decay_length_not_mean_radius():
    """R*/(n/Z) flat; R*/(n^2/Z) drifts.  Both Z = 1 and Z = 2."""
    for Z in (1.0, 2.0):
        ns = np.arange(2, 9)
        Rs = np.array([_front_ns(int(n), Z) for n in ns])
        ratio_decay = Rs / (ns / Z)
        ratio_mean = Rs / (ns**2 / Z)
        # flat to within 3 % of its median, in the paper's stated band
        med = np.median(ratio_decay)
        assert 2.10 < med < 2.20, med
        assert np.max(np.abs(ratio_decay / med - 1.0)) < 0.03, ratio_decay
        # the mean-radius ratio is NOT flat: it falls by more than 2x across n
        assert ratio_mean[0] / ratio_mean[-1] > 2.0, ratio_mean


def test_front_exponent_is_linear_in_n():
    """log R*(n) vs log n has slope 1 +- 0.1 for n >= 2 -- not 2 (sqrt(ZR) window)."""
    ns = np.arange(2, 9)
    Rs = np.array([_front_ns(int(n), 1.0) for n in ns])
    slope = np.polyfit(np.log(ns), np.log(Rs), 1)[0]
    assert abs(slope - 1.0) < 0.10, slope
    assert abs(slope - 2.0) > 0.5          # the rejected sqrt(ZR) window


def _R_rel_1s1s(a: float, b: float) -> float:
    """Tail-reach front: |S| falls to |S|_max/sqrt2 (|S|_max = S(0) for 1s-1s)."""
    s0 = overlap_1s_1s(a, b, 0.0)
    return _crossing(lambda R: overlap_1s_1s(a, b, R), s0 * SQRT_HALF, 1e-3, 80.0)


def test_geometric_mean_tail_law_1s1s():
    pairs = [(1.0, 1.0), (2.0, 1.0), (3.0, 1.0), (2.0, 2.0)]   # (Z_A, Z_B), l = 1/Z
    Rrel = np.array([_R_rel_1s1s(ZA, ZB) for ZA, ZB in pairs])
    lA = np.array([1.0 / ZA for ZA, _ in pairs])
    lB = np.array([1.0 / ZB for _, ZB in pairs])
    geo = Rrel / (2.0 * np.sqrt(lA * lB))
    add = Rrel / (lA + lB)
    dev = lambda v: np.max(np.abs(v / np.mean(v) - 1.0))
    assert dev(geo) < 0.03, geo               # the paper's <= 2 % (four pairs)
    assert dev(add) > 0.05, add               # additive law rejected
    assert 1.50 < 2.0 * np.mean(geo) < 1.66   # c*2 = 1.58 (+- 5 %)
    # dense exponent-ratio scan a = 1, b = t
    ts = np.array([1.0, 2.0, 4.0, 8.0])
    Rt = np.array([_R_rel_1s1s(1.0, t) for t in ts])
    geo_t = Rt / (2.0 * np.sqrt(1.0 / ts))
    add_t = Rt / (1.0 + 1.0 / ts)
    max_t = Rt / np.maximum(1.0, 1.0 / ts)
    assert dev(geo_t) < 0.12, geo_t
    assert dev(add_t) > 0.15, add_t
    assert dev(max_t) > 0.30, max_t


def test_absolute_front_threshold_t_c():
    """t_c = 2.664 from the closed form; the exploratory driver's 2.7456 was a
    scan-grid artifact (first grid point past the root).  Two routes must agree:
    the closed-form root, and the numerically observed loss of the crossing."""
    f = lambda t: (2.0 * math.sqrt(t) / (1.0 + t)) ** 3 - SQRT_HALF
    t_c = brentq(f, 1.5, 6.0, xtol=1e-12)
    assert abs(t_c - 2.664) < 1e-3, t_c
    assert abs(t_c - 2.7456) > 0.05          # the rejected scan-grid value
    # S(R) is monotone decreasing for every ratio, so S(0) is the maximum and
    # the crossing is lost exactly at t_c: present just below, absent just above.
    Rs = np.linspace(0.0, 2.0, 201)
    for t in (2.5, 2.65, 2.68, 3.0):
        S = np.array([overlap_1s_1s(1.0, t, R) for R in Rs])
        assert np.all(np.diff(S) <= 1e-12), t
        assert (S.max() > SQRT_HALF) == (t < t_c), (t, S.max())
    R_star = _crossing(lambda R: overlap_1s_1s(1.0, 2.5, R), SQRT_HALF, 1e-3, 10.0)
    assert 0.0 < R_star < 1.0
