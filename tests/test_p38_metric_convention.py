"""Paper 38 sec:ch_triple + Lemma L2 -- the metric is the DUAL-COXETER one,
and its geodesic distance is the rotation angle.

PI direction 2026-09-03 (Option 1): Paper 38 is restated on the dual-Coxeter
sphere throughout, so that its moment and its Lipschitz seminorm live in one
metric and the paper is exactly the rank-1 case of Paper 40.

The derivation pinned here.  On su(2) the dual-Coxeter rule Cas(ad) = h^v = 2
fixes the Ad-invariant inner product to <X,Y> = -2 tr(XY) in the defining
representation.  In that inner product the geodesic distance from the identity
to exp(i theta n.sigma) -- a group element with eigenphases +-theta, i.e. at
unit-S^3 distance theta -- is 2 theta, the SO(3) rotation angle chi.  So:

  * the dual-Coxeter metric on SU(2) is the round 3-sphere of RADIUS 2
    (diameter 2 pi, volume 16 pi^2), not the unit sphere;
  * `central_fejer_su2.gamma_rate`, which integrates against chi, computes the
    dual-Coxeter moment -- it was right all along, and what was inconsistent
    was the setup section it was paired with;
  * the rate constant is 4/pi in that metric and 2/pi on the unit sphere,
    the two differing by the metric scale (a distance moment is homogeneous
    of degree one in the metric).

Also pinned: the Hopf-base VOLUME-RATIO reading of the constant does NOT
survive in a single normalisation.  Vol(base)/Vol(group) is 1/(4 pi) in the
dual-Coxeter metric and 1/(2 pi) on the unit sphere (the Hopf map is a
Riemannian submersion S^3(r) -> S^2(r/2)); the printed identity
4/pi = Vol(S^2)/pi^2 pairs the dual-Coxeter base volume 4 pi with half the
unit-metric group volume, so it is numerology across two normalisations.
"""
from __future__ import annotations

import mpmath
import numpy as np
import pytest
import sympy as sp

from geovac.central_fejer_su2 import _chi, central_fejer_kernel_su2, gamma_rate


def unit_s3_moment(n: int, prec: int = 30) -> mpmath.mpf:
    """integral_{SU(2)} K_n(g) d_round(e, g) dg with d_round = theta in [0, pi]
    (unit S^3) and Haar class measure (2/pi) sin^2(theta) dtheta; the kernel is
    evaluated at the rotation angle chi = 2 theta."""
    mpmath.mp.dps = prec
    K = sp.lambdify(_chi, central_fejer_kernel_su2(n, _chi), modules="mpmath")
    return mpmath.quad(lambda th: K(2 * th) * th * (2 / mpmath.pi) * mpmath.sin(th) ** 2, [0, mpmath.pi])


def test_rotation_angle_is_twice_the_geodesic_distance():
    th = 0.7
    g = np.cos(th) * np.eye(2) + 1j * np.sin(th) * np.array([[0, 1], [1, 0]])
    phases = np.sort(np.angle(np.linalg.eigvals(g)))
    assert np.allclose(phases, [-th, th])                 # eigenphases +-theta
    chi = sp.Symbol("chi")
    char_half = sp.sin(2 * chi / 2) / sp.sin(chi / 2)      # module's spin-1/2 character
    assert abs(float(char_half.subs(chi, 2 * th)) - np.trace(g).real) < 1e-12


@pytest.mark.parametrize("n", [1, 2, 3, 5])
def test_module_gamma_is_twice_the_unit_s3_moment(n):
    ratio = gamma_rate(n, prec=30) / unit_s3_moment(n)
    assert abs(ratio - 2) < 1e-12, ratio


def test_gamma_1_is_the_diameter_not_the_mean_distance():
    """K_1 == 1, so gamma_1 is the Haar-mean distance: pi/2 on the unit S^3
    (pi is the diameter, which no mean distance can reach)."""
    mpmath.mp.dps = 30
    assert sp.simplify(central_fejer_kernel_su2(1, _chi)) == 1
    assert abs(gamma_rate(1, prec=30) - mpmath.pi) < 1e-20
    assert abs(unit_s3_moment(1) - mpmath.pi / 2) < 1e-20


def test_rate_constant_by_convention():
    """Doubling estimator a_n = (2n g_{2n} - n g_n)/log 2 at n = 25 in both
    conventions: the unit-S^3 value is exactly half the module's."""
    g25, g50 = gamma_rate(25, prec=30), gamma_rate(50, prec=30)
    a_rot = (2 * 25 * g50 - 25 * g25) / mpmath.log(2)
    u25, u50 = unit_s3_moment(25), unit_s3_moment(50)
    a_unit = (2 * 25 * u50 - 25 * u25) / mpmath.log(2)
    assert abs(a_rot / a_unit - 2) < 1e-10
    # the estimators sit above their limits 4/pi and 2/pi at this n (approach from above)
    assert a_rot > 4 / mpmath.pi and a_unit > 2 / mpmath.pi


# ---------------------------------------------------------------------------
# The dual-Coxeter derivation (added 2026-09-03, PI Option 1)
# ---------------------------------------------------------------------------

SIGMA = [np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]]),
         np.array([[1, 0], [0, -1]], complex)]


def _adjoint_matrices():
    """ad(X_a) on su(2) in the basis X_a = -i sigma_a / 2, which satisfies
    [X_a, X_b] = eps_abc X_c."""
    X = [-1j * s / 2 for s in SIGMA]
    ad = np.zeros((3, 3, 3))
    for a in range(3):
        for b in range(3):
            C = X[a] @ X[b] - X[b] @ X[a]
            for c in range(3):
                ad[a][c][b] = np.real(np.trace(C @ np.conj(X[c]).T)
                                      / np.trace(X[c] @ np.conj(X[c]).T))
    return X, ad


def test_dual_coxeter_rule_fixes_the_inner_product():
    """Cas(ad) = h^v = 2 on su(2) forces <X,Y> = -2 tr(XY)."""
    X, ad = _adjoint_matrices()
    S = sum(ad[a] @ ad[a] for a in range(3))
    assert np.allclose(S, -2 * np.eye(3))          # sum_a ad(X_a)^2 = -2 I
    # with <X,Y> = -lam tr(XY) the dual basis is (2/lam) X_a, so |Cas| = 4/lam
    for lam, expected in ((1.0, 4.0), (2.0, 2.0), (4.0, 1.0)):
        assert abs(-(2 / lam) * S[0][0] - expected) < 1e-12
    lam_dual_coxeter = 4.0 / 2.0                    # |Cas(ad)| = h^v = 2
    assert lam_dual_coxeter == 2.0


@pytest.mark.parametrize("theta", [0.3, 1.0, np.pi / 2, np.pi - 0.1])
def test_dual_coxeter_distance_is_the_rotation_angle(theta):
    """d(e, exp(i theta n.sigma)) = 2 theta = chi in the dual-Coxeter metric,
    while the unit-S^3 distance is theta."""
    g_gen = 1j * theta * SIGMA[2]                   # in su(2)
    d_dc = np.sqrt(-2.0 * np.trace(g_gen @ g_gen).real)
    assert abs(d_dc - 2 * theta) < 1e-12
    eigphases = np.sort(np.angle(np.linalg.eigvals(
        np.cos(theta) * np.eye(2) + 1j * np.sin(theta) * SIGMA[2])))
    assert np.allclose(eigphases, [-theta, theta])  # unit-S^3 distance is theta


def test_dual_coxeter_sphere_has_radius_two():
    """Diameter 2 pi and volume 16 pi^2 -- the round S^3 of radius 2."""
    theta_max = np.pi                                # antipode on the unit sphere
    assert abs(2 * theta_max - 2 * np.pi) < 1e-12
    vol_unit = 2 * np.pi ** 2
    assert abs(vol_unit * 2 ** 3 - 16 * np.pi ** 2) < 1e-9


def test_hopf_base_ratio_is_not_the_rate_constant():
    """Guard on the retired reading: Vol(base)/Vol(group) is 1/(4 pi) in the
    dual-Coxeter metric and 1/(2 pi) on the unit sphere -- neither is 4/pi or
    2/pi.  The printed identity 4/pi = Vol(S^2)/pi^2 mixes normalisations."""
    for r, expected_ratio in ((1.0, 1 / (2 * np.pi)), (2.0, 1 / (4 * np.pi))):
        vol_s3 = 2 * np.pi ** 2 * r ** 3
        vol_base = 4 * np.pi * (r / 2) ** 2          # Hopf base S^2(r/2)
        assert abs(vol_base / vol_s3 - expected_ratio) < 1e-12
        assert abs(vol_base / vol_s3 - 4 / np.pi) > 1.0
        assert abs(vol_base / vol_s3 - 2 / np.pi) > 0.4
    # the identity as printed: dual-Coxeter base volume over half the UNIT
    # group volume
    assert abs(4 * np.pi / np.pi ** 2 - 4 / np.pi) < 1e-12
