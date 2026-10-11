"""Pairwise test: the three-centre two-electron integral (XY|XZ) as a sum over
angular channels of PRODUCTS OF TWO-CENTRE FACTORS.

With the shared atom X at the origin, each two-centre density has a momentum
transform rho~(k) = int_0^1 dt e^{i t k.B} G(|k|, t) (Feynman reduction of the
1s-1s Slater product, B = the other centre).  Expanding the plane wave about X
and using the addition theorem on the Coulomb kernel 1/k^2 gives, exactly,

    (XY|XZ) = sum_l P_l(cos gamma) * (2(2l+1)/pi) int_0^inf dk A_l^{XY}(k) A_l^{XZ}(k),
    A_l^{XB}(k) = int_0^1 dt j_l(t k |B|) G(k, t),

gamma = the angle Y-X-Z.  Each A_l involves ONE pair of centres only; the
whole three-centre content is the bond angle, entering through P_l.  The
truncation in l is the only approximation (the elliptic content of Paper 59 is
in the infinite l-sum).  G in closed form (zeta derivatives done analytically):

    G = N_A N_B 2 pi t(1-t) zA zB e^{-D Delta} (D^2 Delta^2 + 3 D Delta + 3) / Delta^5,
    Delta^2 = k^2 t(1-t) + t zA^2 + (1-t) zB^2.

Reference: tests/test_routeC_momentum.py (X, Y, Z, 1s zeta = 1;
GROUND_TRUTH 0.2049417218 from 8-Gaussian fits, true-Slater routes agree to
~4e-7).

Driver: debug/pairwise_three_center.py.
"""
import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.special import spherical_jn

X = np.array([0.0, 0.0, 0.0])
Y = np.array([0.0, 0.0, 2.0])
Z = np.array([1.5, 0.0, 0.5])
GROUND_TRUTH = 0.2049417218262148


def tanh_sinh(n, h):
    """Double-exponential nodes/weights on [0, 1] (endpoint layers in t)."""
    j = np.arange(-n, n + 1) * h
    x = np.tanh(0.5 * np.pi * np.sinh(j))
    w = 0.5 * np.pi * np.cosh(j) / np.cosh(0.5 * np.pi * np.sinh(j)) ** 2 * h
    return 0.5 * (x + 1), 0.5 * w


def G(k, t, D, zA=1.0, zB=1.0):
    N = np.sqrt(zA ** 3 / np.pi) * np.sqrt(zB ** 3 / np.pi)
    Dl = np.sqrt(k ** 2 * t * (1 - t) + t * zA ** 2 + (1 - t) * zB ** 2)
    return N * 2 * np.pi * t * (1 - t) * zA * zB * np.exp(-D * Dl) * (D ** 2 * Dl ** 2 + 3 * D * Dl + 3) / Dl ** 5


def channels(lmax, nk=400, k_scale=3.0, nt=60, ht=0.07):
    t, wt = tanh_sinh(nt, ht)
    keep = (t > 0) & (t < 1)
    t, wt = t[keep], wt[keep]
    u, wu = leggauss(nk)
    u = 0.5 * (u + 1); wu = 0.5 * wu
    k = k_scale * u / (1 - u)
    wk = wu * k_scale / (1 - u) ** 2
    D1, D2 = np.linalg.norm(Y), np.linalg.norm(Z)
    G1 = G(k[:, None], t[None, :], D1)
    G2 = G(k[:, None], t[None, :], D2)
    cosg = float(Y @ Z / (D1 * D2))
    rows = []
    total = 0.0
    for l in range(lmax + 1):
        A1 = (spherical_jn(l, np.outer(k, t) * D1) * G1) @ wt
        A2 = (spherical_jn(l, np.outer(k, t) * D2) * G2) @ wt
        Il = 2 * (2 * l + 1) / np.pi * np.sum(wk * A1 * A2)
        Pl = np.polynomial.legendre.legval(cosg, np.eye(l + 1)[l])
        total += Pl * Il
        rows.append({'l': l, 'I_l': float(Il), 'term': float(Pl * Il), 'partial_sum': float(total)})
    return rows, cosg




# ================================================================ tests
# Referee: McMurchie-Davidson on n_gauss = 16 Gaussian fits of the 1s Slater
# orbitals (geovac.noci_engine.eri_md), measured 2026-10-10; the fit sequence
# n_gauss = 6..16 converges onto the pairwise limit (gap 4.7e-7 -> 2.7e-10).
REFEREE_NG16 = 0.204941652260


def test_pairwise_sum_matches_converged_gaussian_referee():
    rows, _ = channels(20)
    assert abs(rows[-1]['partial_sum'] - REFEREE_NG16) < 1e-9


def test_gaussian_fit_sequence_converges_onto_pairwise_value():
    """Independent route, live: more Gaussians move the MD integral TOWARD the
    pairwise value (excludes a pairwise sum that merely matches one fit)."""
    from geovac import noci_engine as GE
    target = channels(20)[0][-1]['partial_sum']
    gaps = []
    for ng in (6, 10):
        a_, d_, _ = GE.fit_sto_shape(0, 1, n_gauss=ng)
        sh = {"1s": (a_, d_)}
        f = [GE.sto_shape_basis(c, "1s", 1.0, sh, (0, 0, 0)) for c in (X, Y, X, Z)]
        gaps.append(abs(GE.eri_md(*f) - target))
    assert gaps[1] < gaps[0] / 50 and gaps[1] < 1e-8


def test_angle_enters_only_through_legendre():
    """The l-channels I_l are pair-only: rotating Z about the X-Y axis (same bond
    angle) leaves every I_l and the sum unchanged; the l = 0 term alone is NOT
    the answer (three-centre content is real)."""
    global Z
    base, _ = channels(12)
    keep = Z.copy()
    try:
        r = np.hypot(Z[0], Z[1])
        Z = np.array([r * np.cos(1.1), r * np.sin(1.1), Z[2]])
        rot, _ = channels(12)
    finally:
        Z = keep
    for a, b in zip(base, rot):
        assert abs(a['I_l'] - b['I_l']) < 1e-13
    assert abs(base[0]['partial_sum'] - base[-1]['partial_sum']) > 1e-3


def test_fast_angular_convergence():
    rows, _ = channels(20)
    lim = rows[-1]['partial_sum']
    assert abs(rows[3]['partial_sum'] - lim) < 2e-5
    assert abs(rows[10]['partial_sum'] - lim) < 3e-9
    assert abs(rows[14]['partial_sum'] - lim) < 1.5e-9
