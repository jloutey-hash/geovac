"""Paper 38 Lemma L2 -- the metric convention of the mass-concentration moment.

WHY THIS FILE EXISTS (trunk FULL run #3, 2026-09-03, carryforward I.0.2).
`geovac.central_fejer_su2.gamma_rate` integrates K_n(chi) * chi over the
class angle chi in [0, 2 pi] with weight sin^2(chi/2)/pi and calls chi the
"round-S^3 geodesic distance".  The characters sin((2j+1)chi/2)/sin(chi/2)
make chi the SO(3) ROTATION angle, chi = 2 theta, where theta in [0, pi] is
the geodesic distance from e on the UNIT S^3 (g = cos theta + i sin theta
n.sigma has eigenphases +-theta).  Hence

    gamma_n(module) = 2 * integral K_n d_round^{unit S^3}   for every n,

and the rate constant is 4/pi in the rotation-angle normalisation (round
S^3 of radius 2) but 2/pi on the unit S^3 that Paper 38 Sec. 2 states.  The
unconditional theorem is unaffected (a distance moment scales linearly with
the metric); the printed constant is convention-labelled from 2026-09-03.
This test pins the convention so it can never again be silent.
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
