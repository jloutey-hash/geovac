"""Row 8 part 1 (record R8, docs/predictions/2026-11-14_row8_two_centre_l2.md):
two-centre ONE-electron integrals over hydrogenic orbitals with l <= 2, exact.

Orbital = (centre, Z, n, l, m), centre A at the origin, B at R zhat, complex Y_lm,
radial rate a = Z/n (geovac.two_center_eri convention).  In prolate spheroidal
coordinates (foci on the nuclei)

    r_A = R(xi+eta)/2,  r_B = R(xi-eta)/2,  z_A = R(1+xi eta)/2,  z_B = z_A - R,
    rho^2 = (R/2)^2 (xi^2-1)(1-eta^2),  d3r = (R/2)^3 (xi^2-eta^2) dxi deta dphi,

each orbital is  poly(xi, eta) * rho^{|m|} e^{i m phi} e^{-a r_C}; a product with
m_a = m_b carries rho^{2|m|}, a polynomial, and the phi integral is 2 pi.  A single
1/r_A or 1/r_B cancels against one factor of xi^2 - eta^2 = (xi+eta)(xi-eta).
What is left is a polynomial times e^{-p xi - q eta}, integrated exactly by the
Mulliken auxiliaries A_i(p), B_j(q) (geovac.qfd_core).  No quadrature.

Integrals:  S = <a|b>,  VA = <a|1/r_A|b>,  VB = <a|1/r_B|b>,
            T = <a|-1/2 nabla^2|b>  by two routes:
              T_lap   the radial Laplacian (f'' + 2f'/r - l(l+1)f/r^2) Y_lm on b;
              T_eig   T|b> = E_b|b> + (Z_b/r_{C_b})|b>  (b hydrogenic, E_b = -Z_b^2/2n_b^2).
Independent numeric check (`numeric`): scipy quadrature in cylindrical coordinates
from scipy.special.sph_harm + genlaguerre, sharing no code with the exact route.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb, factorial

import sympy as sp

from geovac.qfd_core import A_aux, B_aux
from geovac.two_center_eri import (R_s, _legendre_rational_part, _sph_norm, radial_norm,
                                   radial_poly)

xi, eta = sp.symbols("xi eta")
R = R_s


def _solid(l, m, z, rho2):
    """r^l P_l^m(z/r) with rho^{|m|} stripped (signed m, Condon-Shortley in P_l^m)."""
    w = sp.Symbol("w")
    out = sp.Integer(0)
    for (j,), c in sp.Poly(_legendre_rational_part(l, m, w), w).terms():
        e = l - abs(m) - j
        assert e >= 0 and e % 2 == 0
        out += c * z ** j * (rho2 + z ** 2) ** (e // 2)
    return sp.expand(out)


def orbital(orb, radial_override=None):
    """(poly in xi, eta [no rho^{|m|}], rate, centre).  radial_override: {power: coeff}
    replacing N * sum c_k r^k (powers may go down to l - 1)."""
    cen, Z, n, l, m = orb
    Z = Fraction(Z)
    c, a = radial_poly(Z, n, l)
    N = radial_norm(Z, n, l)
    rad = radial_override if radial_override is not None else {k: N * v for k, v in c.items()}
    rC = R * (xi + eta) / 2 if cen == "A" else R * (xi - eta) / 2
    zA = R * (1 + xi * eta) / 2
    z = zA if cen == "A" else zA - R
    rho2 = R ** 2 / 4 * (xi ** 2 - 1) * (1 - eta ** 2)
    ang = _sph_norm(l, m) * _solid(l, m, z, rho2)
    poly = sp.Integer(0)
    for k, v in rad.items():
        e = k - l
        poly += v * (rC ** e)
    return sp.expand(poly * ang), a, cen


def _integrate(P, p, q):
    P = sp.Poly(sp.expand(P), xi, eta)
    tot = sp.Integer(0)
    for (i, j), c in P.terms():
        tot += c * A_aux(i, p) * B_aux(j, q)
    return sp.expand(tot)


def _pair(oa, ob, kernel=None, radial_b=None):
    """int conj(a) [kernel] b d3r."""
    if oa[4] != ob[4]:
        return sp.Integer(0)
    Pa, aa, ca = orbital(oa)
    Pb, ab, cb = orbital(ob, radial_b)
    mm = abs(oa[4])
    W = (xi ** 2 - 1) * (1 - eta ** 2)
    vol = (xi + eta) * (xi - eta)
    integrand = Pa * Pb * (R / 2) ** (2 * mm) * W ** mm * vol
    if kernel == "A":
        integrand = integrand * 2 / (R * (xi + eta))
    elif kernel == "B":
        integrand = integrand * 2 / (R * (xi - eta))
    integrand = sp.cancel(sp.together(integrand))
    num, den = sp.fraction(integrand)
    # any 1/r_C^k with k > 1 left over would show here
    assert not den.has(xi) and not den.has(eta), den
    integrand = sp.expand(num / den)
    # exponent: e^{-aa r_ca - ab r_cb}
    sa = 1 if ca == "A" else -1
    sb = 1 if cb == "A" else -1
    p = (aa + ab) * R / 2
    q = (sa * aa + sb * ab) * R / 2
    return sp.expand(2 * sp.pi * (R / 2) ** 3 * _integrate(integrand, p, q))


def S(oa, ob):
    return _pair(oa, ob)


def VA(oa, ob):
    return _pair(oa, ob, "A")


def VB(oa, ob):
    return _pair(oa, ob, "B")


def T_lap(oa, ob):
    cen, Z, n, l, m = ob
    c, a = radial_poly(Fraction(Z), n, l)
    N = radial_norm(Fraction(Z), n, l)
    lap = {}
    for k, v in c.items():
        for dk, dc in ((k - 2, sp.Integer(k * (k + 1) - l * (l + 1))),
                       (k - 1, -2 * a * (k + 1)), (k, a ** 2)):
            if dc != 0:
                lap[dk] = lap.get(dk, 0) + N * v * dc
    lap = {k: -v / 2 for k, v in lap.items() if v != 0}
    assert min(lap) >= l - 1
    return _pair(oa, ob, None, lap)


def T_eig(oa, ob):
    cen, Z, n, l, m = ob
    E = -sp.Rational(Fraction(Z).numerator ** 2, Fraction(Z).denominator ** 2) / (2 * n ** 2)
    return sp.expand(E * S(oa, ob) + sp.nsimplify(Fraction(Z)) * _pair(oa, ob, cen))


# ------------------------------------------------------------------ independent numeric


def numeric(oa, ob, Rv, which, dps=12):
    """scipy cylindrical quadrature, independent code path (sph_harm + genlaguerre)."""
    import numpy as np
    from scipy import integrate
    from scipy.special import genlaguerre, sph_harm

    def psi(orb, x_rho, x_z):
        cen, Z, n, l, m = orb
        Z = float(Fraction(Z))
        zc = x_z if cen == "A" else x_z - Rv
        r = np.hypot(x_rho, zc)
        th = np.arctan2(x_rho, zc)
        rho_ = 2 * Z * r / n
        Nn = np.sqrt((2 * Z / n) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
        Rnl = Nn * np.exp(-rho_ / 2) * rho_ ** l * genlaguerre(n - l - 1, 2 * l + 1)(rho_)
        # phi-independent modulus; sph_harm(m, l, azimuth, polar) at azimuth 0
        return Rnl * np.real(sph_harm(m, l, 0.0, th)), r

    def f(x_z, x_rho):
        a, rA = psi(oa, x_rho, x_z)
        b, rB = psi(ob, x_rho, x_z)
        rA_ = np.hypot(x_rho, x_z)
        rB_ = np.hypot(x_rho, x_z - Rv)
        k = {"S": 1.0, "VA": 1.0 / rA_, "VB": 1.0 / rB_}[which]
        return 2 * np.pi * x_rho * a * b * k

    val, err = integrate.dblquad(f, 0, 40, lambda _: -40, lambda _: 40 + Rv,
                                 epsabs=1e-13, epsrel=1e-11)
    return val
