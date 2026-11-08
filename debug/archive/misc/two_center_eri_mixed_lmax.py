"""Brick 1: mixed-exponent, l>0 two-center ERI for ONE-CENTER densities (rho_12 on A,
rho_34 on B) -- generalizes geovac.two_center_eri.aabb_quadrature to per-orbital exponents.
Orbital = (a, n, l, m) with decay a placed on (n,l) via Z = a*n.
Validated against the exact closed form aabb_closed_form at single exponent."""
import os, sys
import numpy as np
from fractions import Fraction
import sympy as sp
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
import geovac.two_center_eri as TC
from geovac.two_center_eri import (multipole_decomposition, V_L_radial, r_s,
                                    sph_norm_numeric, plm_signed, aabb_closed_form)
from scipy import integrate


def _Z(a, n):
    """decay a on principal n  ->  Z with a = Z/n (as exact Fraction)."""
    return Fraction(a).limit_denominator(10**6) * n


def eri_aabb_mixed(o1, o2, o3, o4, R):
    """(o1 o2 | o3 o4), o1,o2 on A, o3,o4 on B. o_i = (a_i, n_i, l_i, m_i)."""
    Z1, Z2 = _Z(o1[0], o1[1]), _Z(o2[0], o2[1])
    Z3, Z4 = _Z(o3[0], o3[1]), _Z(o4[0], o4[1])
    termsA = multipole_decomposition(Z1, o1[1], o1[2], o1[3], Z2, o2[1], o2[2], o2[3])
    termsB = multipole_decomposition(Z3, o3[1], o3[2], o3[3], Z4, o4[1], o4[2], o4[3])
    total = 0.0
    for LA, MA, gA, radA, bA in termsA:
        VA = sp.lambdify(r_s, V_L_radial(radA, bA, LA), "numpy")
        nA = sph_norm_numeric(LA, MA)
        for LB, MB, gB, radB, bB in termsB:
            if MA + MB != 0:
                continue
            radB_f = sp.lambdify(r_s, sum(c * r_s ** k for k, c in radB.items()) * sp.exp(-bB * r_s), "numpy")
            nB = sph_norm_numeric(LB, MB)

            def integrand(u, rb, LA=LA, MA=MA, LB=LB, MB=MB, VA=VA, radB_f=radB_f, nA=nA, nB=nB):
                rA = np.sqrt(rb * rb + R * R + 2 * rb * R * u)
                if rA < 1e-12:
                    return 0.0
                cA = np.clip((R + rb * u) / rA, -1.0, 1.0)
                return float(np.real(radB_f(rb) * rb * rb * nB * plm_signed(LB, MB, np.array([u]))[0]
                                     * VA(rA) * nA * plm_signed(LA, MA, np.array([cA]))[0]))
            val, _ = integrate.dblquad(integrand, 0.0, 60.0, lambda _r: -1.0, lambda _r: 1.0,
                                       epsabs=1e-11, epsrel=1e-11)
            total += float(sp.re(gA * gB)) * 2 * np.pi * val
    return total


if __name__ == "__main__":
    R = 1.5
    print("VALIDATION vs exact aabb_closed_form (single exponent per center):")
    # closed form takes (ZA, orbA1=(n,l,m), orbA2, ZB, orbB1, orbB2, R)
    cases = [
        ("s s | s s", 1.2, (1, 0, 0), (1, 0, 0), 1.2, (1, 0, 0), (1, 0, 0)),
        ("pz s | s s", 1.2, (2, 1, 0), (1, 0, 0), 1.2, (1, 0, 0), (1, 0, 0)),
        ("pz pz| pz pz",1.0,(2, 1, 0), (2, 1, 0), 1.0, (2, 1, 0), (2, 1, 0)),
        ("ppi s | s s", 1.1, (2, 1, 1), (1, 0, 0), 1.1, (1, 0, 0), (1, 0, 0)),
    ]
    for label, ZA, a1, a2, ZB, b1, b2 in cases:
        try:
            ref = float(sp.N(aabb_closed_form(Fraction(ZA), a1, a2, Fraction(ZB), b1, b2, R), 20))
            src = "closed"
        except Exception:
            ref = TC.aabb_quadrature(Fraction(ZA), a1, a2, Fraction(ZB), b1, b2, R)
            src = "quad  "
        o1 = (ZA / a1[0], *a1); o2 = (ZA / a2[0], *a2)
        o3 = (ZB / b1[0], *b1); o4 = (ZB / b2[0], *b2)
        mine = eri_aabb_mixed(o1, o2, o3, o4, R)
        print(f"  {label:14s} {src}={ref:+.9f}  mixed-engine={mine:+.9f}  diff={mine-ref:+.1e}")

    print("\nMIXED-EXPONENT l>0 (no closed form to compare -- the new capability):")
    v = eri_aabb_mixed((1.0, 2, 1, 0), (1.7, 1, 0, 0), (1.3, 2, 1, 0), (0.8, 1, 0, 0), R)
    print(f"  (pz(1.0) s(1.7) | pz(1.3) s(0.8)) = {v:+.9f}")
