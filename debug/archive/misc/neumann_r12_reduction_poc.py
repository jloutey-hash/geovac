"""PoC: odd powers of r12 reduce to the EXISTING general-m Neumann machinery.

Canonical driver for the reduction identity behind the Neumann-r12 build
(explicit correlation in the accurate prolate engine). See
debug/sprint_neumann_r12_build_memo.md.

Key identity (sigma / Dphi-averaged, prolate spheroidal, half_R=1 i.e. R=2):
    r12^2 = A - B cos(Dphi),   A = (z1-z2)^2 + rho1^2 + rho2^2,  B = 2 rho1 rho2
    rho_i = sqrt((xi_i^2-1)(1-eta_i^2)),  z_i = xi_i eta_i

    => r12^(2k+1) = (A - B cosDphi)^(k+1) * (1/r12)

Binomial-expanding (A - B cosDphi)^(k+1) and using cos^j Dphi -> sum of cos(m Dphi)
turns EVERY odd power's Dphi-average into a finite combination of
    <cos(m Dphi) / r12>_phi   (m = 0, 1, 2, ...),
i.e. the m-th Neumann term. m=0 is neumann_vee.py; general m is
neumann_vee_general_m.py (v5.12.7). No new special functions. B^j pairs with the
m=j term's P_l^j (which carry sqrt(xi^2-1) sqrt(1-eta^2)) to give a polynomial.

Proves the reduction three ways, constant-free where possible:
  (R1) r12^1 avg  =  A<1/r12> - B<cosDphi/r12>
  (R2) r12^3 avg  =  A^2<1/r12> - 2AB<cosDphi/r12> + B^2<cos^2Dphi/r12>
  (R3) <1/r12>_phi  =  (2/R) sum_l (2l+1) P_l(xi_<) Q_l(xi_>) P_l(eta1) P_l(eta2)
       (m=0 Neumann term, corpus convention) vs direct quadrature.
  (E)  even powers polynomial: <r12^2>=A, <r12^4>=A^2+B^2/2  (in (R/2)^2 units).
R1/R2 use DIRECT Dphi quadrature on the RHS averages, so they test the algebraic
reduction itself (no expansion, no constants). R3 tests the corpus m=0 machinery.

Result (2026-09-19): R1 1.8e-16, R2 1.8e-16, R3 5.3e-7 (l-trunc), E2/E4 2.6e-14.
"""
import numpy as np
from numpy.polynomial import legendre as L
from scipy.special import lqmn

R = 2.0
half_R = R / 2.0  # = 1

NPHI = 400
x, w = np.polynomial.legendre.leggauss(NPHI)
dphi = np.pi * (x + 1.0)          # [0, 2pi]
wdphi = np.pi * w                 # sum = 2pi
avg = lambda f: np.sum(f * wdphi) / (2 * np.pi)


def geom(xi1, eta1, xi2, eta2):
    rho1 = np.sqrt((xi1**2 - 1) * (1 - eta1**2))
    rho2 = np.sqrt((xi2**2 - 1) * (1 - eta2**2))
    z1, z2 = xi1 * eta1, xi2 * eta2
    A = (z1 - z2)**2 + rho1**2 + rho2**2
    B = 2 * rho1 * rho2
    return A, B


def r12_of(A, B):
    return half_R * np.sqrt(np.maximum(A - B * np.cos(dphi), 1e-300))


def Pl(l, x):
    c = np.zeros(l + 1); c[l] = 1.0
    return L.legval(x, c)


def neumann_m0_inv(xi1, eta1, xi2, eta2, lmax=40):
    xlt, xgt = (xi1, xi2) if xi1 < xi2 else (xi2, xi1)
    Q, _ = lqmn(0, lmax, xgt)
    s = 0.0
    for l in range(lmax + 1):
        s += (2 * l + 1) * Pl(l, xlt) * Q[0, l] * Pl(l, eta1) * Pl(l, eta2)
    return (2.0 / R) * s


def main():
    pts = [
        (1.6, 0.3, 2.4, -0.2), (2.1, -0.5, 1.3, 0.7), (3.0, 0.1, 1.8, 0.4),
        (1.4, 0.6, 2.9, -0.6), (2.5, -0.3, 2.0, 0.2),
    ]
    print(f"{'point':>22} | {'R1 rel':>10} | {'R2 rel':>10} | {'R3(m0) rel':>11} | "
          f"{'E2 rel':>9} | {'E4 rel':>9}")
    print("-" * 92)
    worst = {"R1": 0, "R2": 0, "R3": 0, "E2": 0, "E4": 0}
    for (xi1, eta1, xi2, eta2) in pts:
        A, B = geom(xi1, eta1, xi2, eta2)
        r12 = r12_of(A, B)
        cd = np.cos(dphi)
        inv, cinv, c2inv = avg(1 / r12), avg(cd / r12), avg(cd**2 / r12)
        r1, r3, r2, r4 = avg(r12), avg(r12**3), avg(r12**2), avg(r12**4)
        R1_rhs = half_R**2 * (A * inv - B * cinv)
        R2_rhs = half_R**4 * (A**2 * inv - 2 * A * B * cinv + B**2 * c2inv)
        m0 = neumann_m0_inv(xi1, eta1, xi2, eta2)
        E2_rhs, E4_rhs = half_R**2 * A, half_R**4 * (A**2 + B**2 / 2)
        rel = lambda a, b: abs(a - b) / max(abs(b), 1e-30)
        es = [rel(R1_rhs, r1), rel(R2_rhs, r3), rel(m0, inv),
              rel(E2_rhs, r2), rel(E4_rhs, r4)]
        for k, e in zip(worst, es):
            worst[k] = max(worst[k], e)
        print(f"({xi1},{eta1},{xi2},{eta2})".rjust(22) +
              f" | {es[0]:10.2e} | {es[1]:10.2e} | {es[2]:11.2e} | "
              f"{es[3]:9.2e} | {es[4]:9.2e}")
    print("-" * 92)
    print("worst:", {k: f"{v:.2e}" for k, v in worst.items()})


if __name__ == "__main__":
    main()
