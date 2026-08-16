"""Bet 1: the pinned-source spheroidal-Neumann assembly of the 3-centre 1e integral.

    I = <chi_Y | -Z_X/|r-X| | chi_Z>  =  -Z_X * integral rho(r)/|r-X| d^3r,
    rho = conj(chi_Y) chi_Z  (two-centre density on the Y-Z axis, foci Y and Z).

This is the two-electron EXCHANGE assembly (exchange_value) with electron 2's
density replaced by a point charge at X:  rho_2 -> Z_X delta^3(r-X).  Its effect:
  - eta integral of electron 2  ->  P_tau(eta_X)   (evaluation, no integral)
  - xi integration of electron 2 ->  the ordered split is at a FIXED xi_X
  - prefactor loses one electron's  a^3 * (2 pi)

Neumann kernel and per-tau weight are taken verbatim from the validated
`exchange_value`, so no constant is re-derived -- only the second electron is
pinned.  sigma = m_Z - m_Y = 0 build here (cases A, B); the numerical version
validates structure + prefactor against debug/bet1_three_center_1e_reference.py.

Run:  python debug/bet1_three_center_1e_assembly.py
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

from geovac.two_center_eri import (  # noqa: E402
    _legendre_Q, _z, eta_s, integrate_poly_exp, two_center_spheroidal_product,
    xi_s,
)


def spheroidal_of_X(Yc, Zc, Xc):
    Yc, Zc, Xc = map(np.asarray, (Yc, Zc, Xc))
    R = np.linalg.norm(Zc - Yc)
    rY, rZ = np.linalg.norm(Xc - Yc), np.linalg.norm(Xc - Zc)
    return R, (rY + rZ) / R, (rY - rZ) / R, float(np.arctan2(Xc[1], Xc[0]))


def coeffs_of(P):
    d = {}
    for (j, k), c in sp.Poly(sp.expand((xi_s ** 2 - eta_s ** 2) * P),
                             xi_s, eta_s).terms():
        d[(j, k)] = c
    return d


def three_center_1e_num(ZY, orbY, Yc, ZZ, orbZ, Zc, ZX, Xc, tau_max=16):
    R_YZ, xiX, etaX, phiX = spheroidal_of_X(Yc, Zc, Xc)
    sig = orbZ[2] - orbY[2]
    assert sig == 0, "sigma=0 numerical build (cases A,B)"
    R = sp.Rational(str(R_YZ))
    P1, h1, p1, q1 = two_center_spheroidal_product(
        Fraction(ZY), orbY, Fraction(ZZ), orbZ, R)
    H1 = int(h1)
    c1 = coeffs_of(P1)
    p1f, q1f = float(p1), float(q1)

    total = 0.0
    for tau in range(tau_max + 1):
        Ptau, Qtau = sp.legendre(tau, _z), _legendre_Q(tau)
        Pf = sp.lambdify(_z, Ptau, "numpy")
        Qf = sp.lambdify(_z, Qtau, "numpy")
        PtauEtaX = float(Pf(etaX))
        QxiX, PxiX = float(Qf(xiX)), float(Pf(xiX))

        # eta integral for electron 1 (closed form), per eta-power k
        B1 = {}
        for k in {k for _j, k in c1}:
            integrand = sp.expand(eta_s ** k * (1 - eta_s ** 2) ** H1
                                  * Ptau.subs(_z, eta_s) * sp.exp(-q1 * eta_s))
            B1[k] = float(integrate_poly_exp(integrand, eta_s,
                                             sp.Integer(-1), sp.Integer(1)))

        def xi_pinned(j1):
            def inP(x):
                return x ** j1 * (x * x - 1) ** H1 * np.exp(-p1f * x) * float(Pf(x))

            def outQ(x):
                return x ** j1 * (x * x - 1) ** H1 * np.exp(-p1f * x) * float(Qf(x))
            a1, _ = integrate.quad(inP, 1.0, xiX, epsabs=1e-13, epsrel=1e-12,
                                   limit=200)
            b1, _ = integrate.quad(outQ, xiX, np.inf, epsabs=1e-13, epsrel=1e-12,
                                   limit=200)
            return QxiX * a1 + PxiX * b1

        acc = 0.0
        for (j1, k1), cc1 in c1.items():
            acc += float(cc1) * B1[k1] * PtauEtaX * xi_pinned(j1)
        total += (2 * tau + 1) * acc

    C_pinned = (R_YZ ** 3 / 8) * (2 * np.pi) * (2.0 / R_YZ)
    return -ZX * C_pinned * total


def main():
    Z1 = Fraction(1)
    Yc, Zc, Xc, ZX = (0., 0., 0.), (0., 0., 2.), (1.3, 0., 0.7), 1.0
    R, xiX, etaX, phiX = spheroidal_of_X(Yc, Zc, Xc)
    print(f"X in Y-Z spheroidal: xi_X={xiX:.4f}, eta_X={etaX:.4f}, phi_X={phiX:.4f}\n")

    refs = {"A 1s_Y x 1s_Z": -0.341962099253,
            "B 2p0_Y x 1s_Z": -0.186874906193}
    cases = {"A 1s_Y x 1s_Z": ((1, 0, 0), (1, 0, 0)),
             "B 2p0_Y x 1s_Z": ((2, 1, 0), (1, 0, 0))}
    for name, (oY, oZ) in cases.items():
        for tmax in (10, 16, 22):
            val = three_center_1e_num(Z1, oY, Yc, Z1, oZ, Zc, ZX, Xc, tau_max=tmax)
            ref = refs[name]
            print(f"  {name}  tau_max={tmax:2d}: assembly={val:+.9f}  ref={ref:+.9f}"
                  f"  |diff|={abs(val-ref):.2e}")
        print()


if __name__ == "__main__":
    main()
