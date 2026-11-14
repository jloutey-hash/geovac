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
from math import factorial

# The exact engine was promoted to geovac/two_center_one_electron.py (2026-11-14); this driver
# keeps only the independent numeric referee.
from geovac.two_center_one_electron import S, T_eig, T_lap, VA, VB, orbital, R  # noqa: F401


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
