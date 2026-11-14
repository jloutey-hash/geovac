"""Exact two-centre ONE-electron integrals over hydrogenic orbitals of any l.

Promoted 2026-11-14 from the record-R8 driver (docs/predictions/2026-11-14_row8_two_centre_l2.md;
Paper 58 sec:census).  `geovac.qfd_core` covers the s-type case; this module covers l > 0 and
agrees with qfd_core exactly at l = 0.

Orbital = (centre, Z, n, l, m): centre A at the origin, B at R zhat, complex Y_lm, radial
rate a = Z/n (the `geovac.two_center_eri` convention).  In prolate spheroidal coordinates
(foci on the nuclei)

    r_A = R(xi+eta)/2,  r_B = R(xi-eta)/2,  z_A = R(1+xi eta)/2,  z_B = z_A - R,
    rho^2 = (R/2)^2 (xi^2-1)(1-eta^2),  d3r = (R/2)^3 (xi^2-eta^2) dxi deta dphi,

each orbital is poly(xi, eta) * rho^{|m|} e^{i m phi} e^{-a r_C}.  A product with m_a = m_b
carries rho^{2|m|}, which is a polynomial, and the phi integral is 2 pi.  A single 1/r_A or
1/r_B cancels against one factor of xi^2 - eta^2 = (xi+eta)(xi-eta).  What remains is a
polynomial times e^{-p xi - q eta}, integrated exactly by the Mulliken auxiliaries
A_i(p), B_j(q).  No quadrature.

Integrals:  S = <a|b>,  VA = <a|1/r_A|b>,  VB = <a|1/r_B|b>,
            T = <a|-1/2 nabla^2|b>, by two routes:
              T_lap   the radial Laplacian (f'' + 2f'/r - l(l+1)f/r^2) Y_lm on b;
              T_eig   T|b> = E_b|b> + (Z_b/r_{C_b})|b>  (b hydrogenic, E_b = -Z_b^2/2n_b^2).
"""
from __future__ import annotations

from fractions import Fraction
from typing import Dict, Optional, Tuple

import sympy as sp

from geovac.qfd_core import A_aux, B_aux
from geovac.two_center_eri import (R_s, _legendre_rational_part, _sph_norm, radial_norm,
                                   radial_poly)

Orbital = Tuple[str, Fraction, int, int, int]

xi, eta = sp.symbols("xi eta")
R = R_s

__all__ = ["S", "VA", "VB", "T_lap", "T_eig", "orbital", "R"]


def _solid(l: int, m: int, z: sp.Expr, rho2: sp.Expr) -> sp.Expr:
    """r^l P_l^m(z/r) with rho^{|m|} stripped (signed m, Condon-Shortley in P_l^m)."""
    w = sp.Symbol("w")
    out = sp.Integer(0)
    for (j,), c in sp.Poly(_legendre_rational_part(l, m, w), w).terms():
        e = l - abs(m) - j
        assert e >= 0 and e % 2 == 0
        out += c * z ** j * (rho2 + z ** 2) ** (e // 2)
    return sp.expand(out)


def orbital(orb: Orbital, radial_override: Optional[Dict[int, sp.Expr]] = None
            ) -> Tuple[sp.Expr, sp.Expr, str]:
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
        poly += v * (rC ** (k - l))
    return sp.expand(poly * ang), a, cen


def _integrate(P: sp.Expr, p: sp.Expr, q: sp.Expr) -> sp.Expr:
    P = sp.Poly(sp.expand(P), xi, eta)
    tot = sp.Integer(0)
    for (i, j), c in P.terms():
        tot += c * A_aux(i, p) * B_aux(j, q)
    return sp.expand(tot)


def _pair(oa: Orbital, ob: Orbital, kernel: Optional[str] = None,
          radial_b: Optional[Dict[int, sp.Expr]] = None) -> sp.Expr:
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
    sa = 1 if ca == "A" else -1
    sb = 1 if cb == "A" else -1
    p = (aa + ab) * R / 2
    q = (sa * aa + sb * ab) * R / 2
    return sp.expand(2 * sp.pi * (R / 2) ** 3 * _integrate(integrand, p, q))


def S(oa: Orbital, ob: Orbital) -> sp.Expr:
    """Overlap <a|b>."""
    return _pair(oa, ob)


def VA(oa: Orbital, ob: Orbital) -> sp.Expr:
    """<a|1/r_A|b>."""
    return _pair(oa, ob, "A")


def VB(oa: Orbital, ob: Orbital) -> sp.Expr:
    """<a|1/r_B|b>."""
    return _pair(oa, ob, "B")


def T_lap(oa: Orbital, ob: Orbital) -> sp.Expr:
    """Kinetic <a|-1/2 nabla^2|b> from the radial Laplacian acting on b."""
    cen, Z, n, l, m = ob
    c, a = radial_poly(Fraction(Z), n, l)
    N = radial_norm(Fraction(Z), n, l)
    lap: Dict[int, sp.Expr] = {}
    for k, v in c.items():
        for dk, dc in ((k - 2, sp.Integer(k * (k + 1) - l * (l + 1))),
                       (k - 1, -2 * a * (k + 1)), (k, a ** 2)):
            if dc != 0:
                lap[dk] = lap.get(dk, 0) + N * v * dc
    lap = {k: -v / 2 for k, v in lap.items() if v != 0}
    assert min(lap) >= l - 1
    return _pair(oa, ob, None, lap)


def T_eig(oa: Orbital, ob: Orbital) -> sp.Expr:
    """Kinetic <a|-1/2 nabla^2|b> from b's own Schroedinger equation."""
    cen, Z, n, l, m = ob
    Zf = Fraction(Z)
    E = -sp.Rational(Zf.numerator ** 2, Zf.denominator ** 2) / (2 * n ** 2)
    return sp.expand(E * S(oa, ob) + sp.nsimplify(Zf) * _pair(oa, ob, cen))
