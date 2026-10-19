"""Paper 23 sec:finite_nucleus, the Dirac case near and past Z alpha = 1.

Backs the paragraph "Near and past Z alpha = 1": (i) the 1s finite-size shift scales as R^(2 gamma);
(iii) the critical charge Z_cr = 169.71 for R = 1.2 (2.5 Z)^(1/3) fm; (iv) successive diving radii at
Z alpha = 1.2 in ratio e^(pi/sigma). Record: docs/predictions/2026-10-18_relativistic_nucleus.md.

Independent of the mpmath script that produced the paper's numbers: float64, scipy DOP853 for the
exterior, and the interior series written in the scaled variable t = r/R.

Uniformly charged sphere, kappa = -1, atomic units, G = r g, F = r f, total energy W:
  dG/dr = -(kappa/r) G + (W - V + c^2)/c F,   dF/dr = (kappa/r) F - (W - V - c^2)/c G.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

C = 137.035999084
KAPPA = -1
BOHR_FM = 52917.72109


def interior_ratio(W: float, Z: float, R: float, nterms: int = 120) -> float:
    """F/G at r = R for the regular interior solution; series in t = r/R."""
    PR = (W + C**2) * R + 1.5 * Z      # (W + V0 + c^2) R
    MR = (W - C**2) * R + 1.5 * Z      # (W + V0 - c^2) R
    BR = 0.5 * Z                       # B R^3
    a = np.zeros(nterms + 3)
    b = np.zeros(nterms + 3)
    a[1] = 1.0                         # kappa = -1: G ~ r
    for k in range(2, nterms):
        b3 = b[k - 3] if k >= 3 else 0.0
        a3 = a[k - 3] if k >= 3 else 0.0
        a[k] = (PR * b[k - 1] - BR * b3) / C / (k + KAPPA)
        b[k] = -(MR * a[k - 1] - BR * a3) / C / (k - KAPPA)
    return float(np.sum(b[:nterms]) / np.sum(a[:nterms]))


def exterior_at(W: float, Z: float, s_points: list[float]) -> list[np.ndarray]:
    """Decaying exterior (G, F), normalized, at each s = ln r in s_points (descending)."""
    if W > -C**2:
        lam = math.sqrt(C**4 - W**2) / C
        r_max = 45.0 / lam
        y0 = [1.0, -math.sqrt((C**2 - W) / (C**2 + W))]
    else:
        r_max = 3600.0 / (8.0 * Z)
        dlog = -math.sqrt(2.0 * Z / r_max) - 1.0 / (4.0 * r_max)
        y0 = [1.0, (C * r_max / Z) * (dlog + KAPPA / r_max)]

    def rhs(s, y):
        r = math.exp(s)
        G, F = y
        return [-KAPPA * G + ((W + C**2) * r + Z) / C * F,
                KAPPA * F - ((W - C**2) * r + Z) / C * G]

    sol = solve_ivp(rhs, (math.log(r_max), min(s_points)), y0, method="DOP853",
                    rtol=1e-13, atol=1e-300, dense_output=True)
    assert sol.success
    out = []
    for s in s_points:
        y = sol.sol(s)
        out.append(y / np.hypot(y[0], y[1]))
    return out


def match(W: float, Z: float, R: float) -> float:
    G, F = exterior_at(W, Z, [math.log(R)])[0]
    return interior_ratio(W, Z, R) * G - F


def sommerfeld_1s(Z: float) -> float:
    return C**2 * math.sqrt(1.0 - (Z / C) ** 2)


def test_shift_exponent_is_two_gamma():
    """(i) d ln(Delta W)/d ln R = 2 gamma at Z alpha = 0.9, R = 1e-8 .. 1e-7."""
    za = 0.9
    Z = za * C
    W0 = sommerfeld_1s(Z)
    shifts = []
    for R in (1e-8, 1e-7):
        W = brentq(lambda w: match(w, Z, R), W0 + 1e-6, W0 + 5.0, xtol=1e-12, rtol=1e-15)
        shifts.append(W - W0)
    slope = math.log(shifts[1] / shifts[0]) / math.log(10.0)
    expected = 2.0 * math.sqrt(1.0 - za**2)
    assert abs(slope - expected) < 1e-3, (slope, expected)


def test_critical_charge_for_the_radius_law():
    """(iii) the 1s level reaches W = -c^2 between Z = 169.70 and 169.72."""
    Rz = lambda Z: 1.2 * (2.5 * Z) ** (1.0 / 3.0) / BOHR_FM
    above = match(-C**2, 169.70, Rz(169.70))
    below = match(-C**2, 169.72, Rz(169.72))
    assert above > 0.0 > below, (above, below)


def test_diving_radii_ratio_is_e_to_pi_over_sigma():
    """(iv) at Z alpha = 1.2 the second and third diving radii are in ratio e^(pi/sigma)."""
    za = 1.2
    Z = za * C
    sigma = math.sqrt(za**2 - 1.0)
    W = -C**2
    roots = []
    for guess in (9.13329e-7, 8.01197e-9):
        lo, hi = math.log(guess) - 0.05, math.log(guess) + 0.05

        def D(s):
            G, F = exterior_at(W, Z, [s])[0]
            return interior_ratio(W, Z, math.exp(s)) * G - F

        assert D(lo) * D(hi) < 0.0, "no diving radius in the bracket"
        roots.append(math.exp(brentq(D, lo, hi, xtol=1e-13)))
    ratio = roots[0] / roots[1]
    expected = math.exp(math.pi / sigma)
    assert abs(ratio / expected - 1.0) < 2e-4, (ratio, expected)
