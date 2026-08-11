"""Full-build increment 1b: close the (AA|BB) class end-to-end.

Increment 1 built and validated the one-center pieces: the exact multipole
decomposition of an orbital product, and the potential V_A it generates. This
closes the chain:

    (ab|cd) = integral rho_B(r2) V_A(r2) d3r2

evaluated with A at the origin and B at R zhat. Because both distributions are
one-center, the phi integral is trivial (both share the z axis, so it enforces
M_A + M_B = 0) and what remains is a 2D integral over (r_B, theta_B) with

    r_A     = sqrt(r_B^2 + R^2 + 2 r_B R cos theta_B)
    cos th_A = (R + r_B cos theta_B) / r_A

DELIBERATE SEQUENCING. This increment integrates NUMERICALLY on purpose. The
closed form (increment 1c) replaces the quadrature only after the decomposition +
normalization chain is known good -- otherwise a disagreement with the reference
cannot be localized between "decomposition wrong" and "closed form wrong". Get
the physics chain right first, then make it algebraic.

THE REFERENCE, and why it is better than eri_md alone. For two 1s densities with
a common exponent zeta there is an exact closed form (the classic VB "J" integral):

    J(R) = (1/R) [ 1 - e^{-2 rho} (1 + (11/8) rho + (3/4) rho^2 + (1/6) rho^3) ],
    rho = zeta R

so this increment is checked against an EXACT analytic result, not only against
the Gaussian-fitted MD engine. eri_md is used as a second, independent check
(agreement expected at ~1e-6, limited by the STO->Gaussian fits, not by us).

Run from repo root:  python debug/eri_aabb_twocenter.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from scipy import integrate

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "debug"))

from eri_aabb_multipole import (  # noqa: E402
    V_L_radial, multipole_decomposition, r_s,
)
from geovac import noci_engine as E  # noqa: E402


def J_exact(zeta: float, R: float) -> float:
    """Exact two-center Coulomb integral between two 1s densities."""
    rho = zeta * R
    return (1.0 / R) * (1.0 - np.exp(-2 * rho)
                        * (1 + 1.375 * rho + 0.75 * rho ** 2 + rho ** 3 / 6.0))


def aabb_numeric(Z, na, la, ma, nb, lb, mb, Zb, nc, lc, mc, nd, ld, md, R):
    """(ab|cd) with a,b on A (origin) and c,d on B (R zhat), by 2D quadrature."""
    termsA = multipole_decomposition(Z, na, la, ma, Z, nb, lb, mb)
    termsB = multipole_decomposition(Zb, nc, lc, mc, Zb, nd, ld, md)

    total = 0.0
    for LA, MA, gA, radA, bA in termsA:
        VA = sp.lambdify(r_s, sp.re(V_L_radial(radA, bA, LA)), "numpy")
        for LB, MB, gB, radB, bB in termsB:
            # phi integral: both harmonics share the z axis -> needs MA + MB = 0
            if MA + MB != 0:
                continue
            radB_f = sp.lambdify(
                r_s, sum(c * r_s ** k for k, c in radB.items()) * sp.exp(-bB * r_s),
                "numpy")
            # theta parts, as functions of cos(theta_B)
            PA = sp.lambdify(sp.Symbol("x"),
                             sp.simplify(sp.assoc_legendre(LA, abs(MA),
                                                           sp.Symbol("x"))), "numpy")
            PB = sp.lambdify(sp.Symbol("x"),
                             sp.simplify(sp.assoc_legendre(LB, abs(MB),
                                                           sp.Symbol("x"))), "numpy")
            normA = float(sp.sqrt((2 * LA + 1) / (4 * sp.pi)
                                 * sp.factorial(LA - abs(MA))
                                 / sp.factorial(LA + abs(MA))))
            normB = float(sp.sqrt((2 * LB + 1) / (4 * sp.pi)
                                  * sp.factorial(LB - abs(MB))
                                  / sp.factorial(LB + abs(MB))))

            def integrand(u, rb, PA=PA, PB=PB, VA=VA, radB_f=radB_f):
                rA = np.sqrt(rb * rb + R * R + 2 * rb * R * u)
                if rA < 1e-12:
                    return 0.0
                cA = (R + rb * u) / rA
                cA = min(1.0, max(-1.0, cA))
                return (radB_f(rb) * rb * rb * PB(u) * normB
                        * VA(rA) * PA(cA) * normA)

            val, _err = integrate.dblquad(
                integrand, 0.0, 60.0, lambda _rb: -1.0, lambda _rb: 1.0,
                epsabs=1e-11, epsrel=1e-11)
            total += float(gA) * float(gB) * 2 * np.pi * val
    return total


def aabb_md(zeta_a, na, la_lmn, zeta_b, nb, lb_lmn, R, shapes):
    """Same integral through the independent McMurchie-Davidson engine."""
    pa, pb = np.array([0., 0., 0.]), np.array([0., 0., R])
    a = E.sto_shape_basis(pa, na, zeta_a, shapes, la_lmn)
    b = E.sto_shape_basis(pb, nb, zeta_b, shapes, lb_lmn)
    return E.eri_md(a, a, b, b)


def main() -> None:
    print("Increment 1b -- (AA|BB) end to end, against an EXACT reference\n")
    Z = Fraction(1)
    shapes = {}
    for kind, (l, n_r) in (("1s", (0, 1)), ("2s", (0, 2)), ("2p", (1, 2))):
        arr, dco, q = E.fit_sto_shape(l, n_r)
        shapes[kind] = (arr, dco)

    print("Case: (1s_A 1s_A | 1s_B 1s_B), zeta = 1 both centers\n")
    print(f"{'R':>6}{'this build':>18}{'exact J(R)':>18}{'|diff|':>12}"
          f"{'eri_md':>16}{'md diff':>12}")
    print("-" * 84)
    worst_exact = 0.0
    worst_md = 0.0
    for R in (1.5, 2.5, 4.0):
        got = aabb_numeric(Z, 1, 0, 0, 1, 0, 0, Z, 1, 0, 0, 1, 0, 0, R)
        ref = J_exact(1.0, R)
        md = aabb_md("1s", "1s", (0, 0, 0), 1.0, "1s", (0, 0, 0), R, shapes) \
            if False else aabb_md(1.0, "1s", (0, 0, 0), 1.0, "1s", (0, 0, 0), R, shapes)
        d1, d2 = abs(got - ref), abs(md - ref)
        worst_exact = max(worst_exact, d1)
        worst_md = max(worst_md, d2)
        print(f"{R:>6.2f}{got:>18.12f}{ref:>18.12f}{d1:>12.2e}"
              f"{md:>16.12f}{d2:>12.2e}")
    print("-" * 84)
    print(f"\nworst |this build - exact| = {worst_exact:.2e}")
    print(f"worst |eri_md      - exact| = {worst_md:.2e}   "
          f"(Gaussian-fit limited, not ours)")
    ok = worst_exact < 1e-8
    print(f"\nGATE: {'PASS' if ok else 'FAIL'} -- the decomposition + potential + "
          f"two-center chain\nreproduces the exact analytic J(R).")
    if ok:
        print("\n=> (AA|BB) is closed end to end numerically. Increment 1c can now")
        print("   replace the quadrature with a closed form, with a trustworthy")
        print("   reference to check it against at every step.")


if __name__ == "__main__":
    main()
