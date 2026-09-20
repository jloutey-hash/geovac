"""H2 Neumann-r12 build, increment 5: the p0xp1 KINETIC block reduces to the odd
(m=0,m=1) machinery -- the last unvalidated block.

Prolate Green form, bra g_i (p=0, phi-indep so d_phi phi_i = 0 -> azimuthal term
vanishes), ket phi_j = g_j r12 (p=1):
  T_ij = (R/4)(2pi)^2 sum_{xi,eta} w <T1 + T2>_phi
  <T1>_phi = [(xi1^2-1) d_xi1 g_i <d_xi1 phi_j>_phi
            + (1-eta1^2) d_eta1 g_i <d_eta1 phi_j>_phi] J2   (J = (R/2)^3(xi^2-eta^2))
  <d_xi1 phi_j>_phi = d_xi1 g_j <r12>_phi
                    + g_j (R/2)^2 [d_xi1 A * K0 - d_xi1 B * K1]/2
  <r12>_phi = (R/2)^2 [A K0 - B K1],  K0=<1/r12>_phi, K1=<cosDphi/r12>_phi (Neumann sums)

Validate the p0xp1 kinetic block (pointwise Neumann kernels + xi,eta quadrature)
against the crude analytical kinetic for a mixed p={0,1} basis (its p0xp1 subblock).
"""
import numpy as np
from math import factorial
from numpy.polynomial import legendre as LEG
from scipy.special import lqmn
from geovac.hylleraas import (
    HylleraasBasisFunction, build_quadrature_grids, evaluate_basis_and_derivs,
    _full_derivs_python, _r12_and_derivs_python,
)

R = 1.4011
ALPHA = 1.0
C = 2.0 * ALPHA
half_R = R / 2.0
LMAX = 20


def _leg(l):
    c = np.zeros(l + 1); c[l] = 1.0
    return c


def _dleg(l):
    return LEG.legder(_leg(l))


def kernels(x1, e1, x2, e2):
    """K0 = <1/r12>_phi, K1 = <cosDphi/r12>_phi (pointwise Neumann sums)."""
    xlt, xgt = (x1, x2) if x1 < x2 else (x2, x1)
    Q0, Q0p = lqmn(0, LMAX, xgt)       # Q0[0,l], Q0p[0,l]=dQ_l/dxi
    K0 = 0.0; K1 = 0.0
    sx = np.sqrt(xgt**2 - 1.0)
    for l in range(LMAX + 1):
        Pl_lt = LEG.legval(xlt, _leg(l))
        Pl_e1 = LEG.legval(e1, _leg(l)); Pl_e2 = LEG.legval(e2, _leg(l))
        K0 += (2 * l + 1) * Pl_lt * Q0[0, l] * Pl_e1 * Pl_e2
        if l >= 1:
            dP_lt = LEG.legval(xlt, _dleg(l))
            dP_e1 = LEG.legval(e1, _dleg(l)); dP_e2 = LEG.legval(e2, _dleg(l))
            # P_l^1(xi_<) = sqrt(xi_<^2-1) dP; Q_l^1(xi_>) = sqrt(xi_>^2-1) dQ
            Pl1_lt = np.sqrt(xlt**2 - 1.0) * dP_lt
            Ql1_gt = sx * Q0p[0, l]
            Pl1_e1 = np.sqrt(1 - e1**2) * dP_e1
            Pl1_e2 = np.sqrt(1 - e2**2) * dP_e2
            c = (2 * l + 1) * (factorial(l - 1) / factorial(l + 1))**2
            K1 += -c * Pl1_lt * Ql1_gt * Pl1_e1 * Pl1_e2
    return (2.0 / R) * K0, (2.0 / R) * K1


def dA_dB(x1, e1, x2, e2):
    diff = x1 * e1 - x2 * e2
    rho1 = np.sqrt(max((x1**2 - 1) * (1 - e1**2), 1e-300))
    rho2 = np.sqrt(max((x2**2 - 1) * (1 - e2**2), 0.0))
    dA_x1 = 2 * diff * e1 + 2 * x1 * (1 - e1**2)
    dA_e1 = 2 * diff * x1 - 2 * e1 * (x1**2 - 1)
    dB_x1 = 2 * rho2 * x1 * (1 - e1**2) / rho1
    dB_e1 = -2 * rho2 * e1 * (x1**2 - 1) / rho1
    return dA_x1, dA_e1, dB_x1, dB_e1


def reduction_p0p1(p0_basis, p1_basis, g):
    xi, wxi, eta, weta = g['xi'], g['w_xi'], g['eta'], g['w_eta']
    n0, n1 = len(p0_basis), len(p1_basis)
    T = np.zeros((n0, n1))
    pref = (R / 4.0) * (2 * np.pi)**2
    for a in range(len(xi)):
        x1 = xi[a]
        for c in range(len(xi)):
            x2 = xi[c]
            for b in range(len(eta)):
                e1 = eta[b]; J1 = half_R**3 * (x1**2 - e1**2)
                for d in range(len(eta)):
                    e2 = eta[d]; J2 = half_R**3 * (x2**2 - e2**2)
                    K0, K1 = kernels(x1, e1, x2, e2)
                    rbar = half_R**2 * (
                        ((x1 * e1 - x2 * e2)**2 + (x1**2 - 1) * (1 - e1**2)
                         + (x2**2 - 1) * (1 - e2**2)) * K0
                        - 2 * np.sqrt(max((x1**2 - 1) * (1 - e1**2), 0.0)
                                      * max((x2**2 - 1) * (1 - e2**2), 0.0)) * K1)
                    dA_x1, dA_e1, dB_x1, dB_e1 = dA_dB(x1, e1, x2, e2)
                    dA_x2, dA_e2, dB_x2, dB_e2 = dA_dB(x2, e2, x1, e1)
                    # <d r12>_phi for each coordinate
                    dr_x1 = half_R**2 * (dA_x1 * K0 - dB_x1 * K1) / 2
                    dr_e1 = half_R**2 * (dA_e1 * K0 - dB_e1 * K1) / 2
                    dr_x2 = half_R**2 * (dA_x2 * K0 - dB_x2 * K1) / 2
                    dr_e2 = half_R**2 * (dA_e2 * K0 - dB_e2 * K1) / 2
                    w = wxi[a] * wxi[c] * weta[b] * weta[d]
                    gi = [evaluate_basis_and_derivs(bf, x1, e1, x2, e2) for bf in p0_basis]
                    gj = [evaluate_basis_and_derivs(bf, x1, e1, x2, e2) for bf in p1_basis]
                    for i in range(n0):
                        _, gix1, gie1, gix2, gie2 = gi[i]
                        for j in range(n1):
                            gjv, gjx1, gje1, gjx2, gje2 = gj[j]
                            # <d_x1 phi_j>_phi etc.
                            dphij_x1 = gjx1 * rbar + gjv * dr_x1
                            dphij_e1 = gje1 * rbar + gjv * dr_e1
                            dphij_x2 = gjx2 * rbar + gjv * dr_x2
                            dphij_e2 = gje2 * rbar + gjv * dr_e2
                            T1 = ((x1**2 - 1) * gix1 * dphij_x1
                                  + (1 - e1**2) * gie1 * dphij_e1) * J2
                            T2 = ((x2**2 - 1) * gix2 * dphij_x2
                                  + (1 - e2**2) * gie2 * dphij_e2) * J1
                            T[i, j] += w * (T1 + T2)
    return pref * T


def crude_p0p1(p0_basis, p1_basis, g):
    """Crude analytical kinetic, p0(bra) x p1(ket) block (ground truth)."""
    xi, wxi, eta, weta = g['xi'], g['w_xi'], g['eta'], g['w_eta']
    dphi, wphi = g['dphi'], g['w_phi']
    n0, n1 = len(p0_basis), len(p1_basis)
    T = np.zeros((n0, n1)); Tpref = R / 4.0
    for a in range(len(xi)):
        x1 = xi[a]; x1m1 = x1**2 - 1
        for c in range(len(xi)):
            x2 = xi[c]; x2m1 = x2**2 - 1; wac = wxi[a] * wxi[c]
            for b in range(len(eta)):
                e1 = eta[b]; om1 = 1 - e1**2; J1 = half_R**3 * (x1**2 - e1**2)
                for d in range(len(eta)):
                    e2 = eta[d]; om2 = 1 - e2**2; J2 = half_R**3 * (x2**2 - e2**2)
                    Tacc = np.zeros((n0, n1))
                    for ip in range(len(dphi)):
                        r12, drx1, dre1, drx2, dre2, drdp = _r12_and_derivs_python(
                            x1, e1, x2, e2, dphi[ip], R)
                        r12 = max(r12, 1e-15)
                        Di = [_full_derivs_python(bf, x1, e1, x2, e2, r12,
                              drx1, dre1, drx2, dre2, drdp) for bf in p0_basis]
                        Dj = [_full_derivs_python(bf, x1, e1, x2, e2, r12,
                              drx1, dre1, drx2, dre2, drdp) for bf in p1_basis]
                        for i in range(n0):
                            _, ix1, ie1, ix2, ie2, _ = Di[i]   # p0 bra: d_phi=0
                            for j in range(n1):
                                _, jx1, je1, jx2, je2, _ = Dj[j]
                                T1 = (x1m1 * ix1 * jx1 + om1 * ie1 * je1) * J2
                                T2 = (x2m1 * ix2 * jx2 + om2 * ie2 * je2) * J1
                                Tacc[i, j] += (T1 + T2) * wphi[ip]
                    T += wac * weta[b] * weta[d] * Tpref * 2 * np.pi * Tacc
    return T


def main():
    p0 = [HylleraasBasisFunction(j, k, l, m, 0, ALPHA)
          for (j, k, l, m) in [(0, 0, 0, 0), (1, 0, 0, 0)]]
    p1 = [HylleraasBasisFunction(j, k, l, m, 1, ALPHA)
          for (j, k, l, m) in [(0, 0, 0, 0), (1, 0, 0, 0)]]
    g = build_quadrature_grids(N_xi=18, N_eta=14, N_phi=24, xi_max=14.0)
    Tred = reduction_p0p1(p0, p1, g)
    Tcru = crude_p0p1(p0, p1, g)
    rel = np.abs(Tred - Tcru) / np.maximum(np.abs(Tcru), 1e-12)
    print("T p0xp1 reduction:\n", np.round(Tred, 6))
    print("T p0xp1 crude:\n", np.round(Tcru, 6))
    print(f"max rel diff = {np.max(rel):.3e}")
    print("p0xp1 KINETIC", "VALIDATED" if np.max(rel) < 5e-3 else "MISMATCH")


if __name__ == "__main__":
    main()
