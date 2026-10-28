"""Guard for Paper 59's water passage, the bond-angle singular point of the shared-centre three-centre integral
(record docs/predictions/2026-10-27_water_angle_singularities.md; WA0, WA3 held, WA1/WA2 scored as measured).

(XY|XZ) over 1s Slater orbitals, X at the origin, |XY| = D1, |XZ| = D2, bond angle gamma, x = cos gamma, is the
Legendre series sum_l F_l P_l(x) with F_l = (2(2l+1)/pi) int dk A_l^{XY}(k) A_l^{XZ}(k) (the channel formula of
tests/test_pairwise_three_center.py, re-implemented here exponent-general).  Claims guarded, as the paper states them:
  * the channel sum reproduces the Slater referee at the test geometry (control);
  * F_l keeps one sign in every case tested, so the linear geometry x = -1 is not the nearest singular point;
  * at D1 = D2 the terms fall like a power, about l^-8 (no geometric decay);
  * at D1 != D2 they fall geometrically at a rate within 3% of D</D> (above it; the registered 1% was missed).
"""
import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss, legval
from scipy.special import spherical_jn

REFEREE_NG16 = 0.204941652260   # Slater referee of tests/test_pairwise_three_center.py
WATER = 1.809
CASES = [(WATER, WATER), (WATER, 2.5), (WATER, 1.0)]
SETS = {'E1': (1.0, 1.0), 'E2': (2.25, 1.24)}     # (zeta on X, zeta on Y and Z)
LHI = {(WATER, WATER): 30, (WATER, 2.5): 30, (WATER, 1.0): 23}


def _tanh_sinh(n, h):
    j = np.arange(-n, n + 1) * h
    x = np.tanh(0.5 * np.pi * np.sinh(j))
    w = 0.5 * np.pi * np.cosh(j) / np.cosh(0.5 * np.pi * np.sinh(j)) ** 2 * h
    t, wt = 0.5 * (x + 1), 0.5 * w
    keep = (t > 0) & (t < 1)
    return t[keep], wt[keep]


def _G(k, t, D, zA, zB):
    N = np.sqrt(zA ** 3 / np.pi) * np.sqrt(zB ** 3 / np.pi)
    Dl = np.sqrt(k ** 2 * t * (1 - t) + t * zA ** 2 + (1 - t) * zB ** 2)
    return N * 2 * np.pi * t * (1 - t) * zA * zB * np.exp(-D * Dl) * (D ** 2 * Dl ** 2 + 3 * D * Dl + 3) / Dl ** 5


def F_coeffs(D1, D2, zX, zH, lmax, level=0, k_scale=3.0):
    t, wt = _tanh_sinh(60 * 2 ** level, 0.07 / 2 ** level)
    u, wu = leggauss(400 * 2 ** level)
    u, wu = 0.5 * (u + 1), 0.5 * wu
    k = k_scale * u / (1 - u)
    wk = wu * k_scale / (1 - u) ** 2
    G1 = _G(k[:, None], t[None, :], D1, zX, zH)
    G2 = _G(k[:, None], t[None, :], D2, zX, zH)
    out = []
    for l in range(lmax + 1):
        A1 = (spherical_jn(l, np.outer(k, t) * D1) * G1) @ wt
        A2 = (spherical_jn(l, np.outer(k, t) * D2) * G2) @ wt
        out.append(2 * (2 * l + 1) / np.pi * np.sum(wk * A1 * A2))
    return np.array(out)


def fit_rate(F, lmin, lmax):
    """The registered fit: log|F_l| = a - p log l + l log r."""
    ls = np.arange(lmin, lmax + 1)
    A = np.vstack([np.ones(len(ls)), -np.log(ls), ls.astype(float)]).T
    coef = np.linalg.lstsq(A, np.log(np.abs(F[ls])), rcond=None)[0]
    return coef[1], coef[2]          # p, log r


def test_channel_sum_reproduces_the_slater_referee():
    Y, Z = np.array([0, 0, 2.0]), np.array([1.5, 0, 0.5])
    D1, D2 = np.linalg.norm(Y), np.linalg.norm(Z)
    F = F_coeffs(D1, D2, 1.0, 1.0, 20)
    assert abs(legval(Y @ Z / (D1 * D2), F) - REFEREE_NG16) < 1e-9


def test_symmetric_water_terms_keep_one_sign():
    F = F_coeffs(WATER, WATER, 1.0, 1.0, 30)
    assert np.all(F[15:31] > 0)


@pytest.mark.slow
def test_bond_angle_singular_point_all_cases():
    for name, (zX, zH) in SETS.items():
        for (D1, D2) in CASES:
            lhi = LHI[(D1, D2)]
            F0 = F_coeffs(D1, D2, zX, zH, lhi, level=0)
            F1 = F_coeffs(D1, D2, zX, zH, lhi, level=1)
            ls = np.arange(15, lhi + 1)
            assert np.all(np.abs(F1 - F0)[ls] < 1e-2 * np.abs(F1[ls])), (name, D1, D2, "above the noise floor")
            assert np.all(F1[ls] > 0), (name, D1, D2, "one sign: x = -1 is not the nearest singular point")
            p, logr = fit_rate(F1, 15, lhi)
            if D1 == D2:
                assert abs(logr) < 0.02 and 7.0 < p < 8.5, (name, p, logr)
            else:
                ratio = np.exp(logr) / (min(D1, D2) / max(D1, D2))
                assert 1.0 < ratio < 1.03, (name, D1, D2, ratio)
