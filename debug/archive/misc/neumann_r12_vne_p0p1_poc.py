"""H2 Neumann-r12 build, increment 4: the p0xp1 V_ne cross term is NOT a
two-Coulomb "hybrid" integral for HOMONUCLEAR H2. See sprint_neumann_r12_build_memo.md.

Homonuclear key: 1/r1A + 1/r1B = (2/R)[1/(xi1+eta1)+1/(xi1-eta1)]
                                = (4/R) xi1/(xi1^2-eta1^2),
so (1/r1A+1/r1B)*J1 = (4/R)xi1*(R/2)^3 = (R^2/2) xi1 -- the (xi1^2-eta1^2)
CANCELS. Hence V_ne * J1 J2 = -(R^2/2)(R/2)^3 [xi1(xi2^2-eta2^2)+xi2(xi1^2-eta1^2)]
= -(R^2/2)(R/2)^3 P_Vne, a polynomial. So

  <g_i g_j r12 V_ne>  =  odd r12^1 moment  with  P_Vne  in place of Jp1 Jp2.

Validate the p0xp1 V_ne block (bra g_i p=0, ket g_j r12 p=1) via the odd
machinery (reused from neumann_r12_oddpower_vee_poc) vs 5D quadrature.
Result (2026-09-19): max rel diff 8.3e-6. The "remaining hard integral" dissolves.
"""
import numpy as np
from math import factorial
from geovac.hylleraas import HylleraasBasisFunction, build_quadrature_grids
from geovac.neumann_vee import _get_unsym_terms
from neumann_r12_oddpower_vee_poc import (   # sibling import (run from debug/)
    mul, sym_product, Cl, Dl, build_X, _xi_grid, A_TERMS, LMAX, R, ALPHA, C,
)

# P_Vne = xi1(xi2^2-eta2^2) + xi2(xi1^2-eta1^2)  (coef, dP1,dQ1,dP2,dQ2)
P_VNE_TERMS = [
    (1.0, 1, 0, 2, 0), (-1.0, 1, 0, 0, 2),
    (1.0, 2, 0, 1, 0), (-1.0, 0, 2, 1, 0),
]


def vne_reduction(p0_basis, p1_basis):
    XI, WXI = _xi_grid()
    X0, X1 = build_X(0, XI, WXI), build_X(1, XI, WXI)
    cl1 = np.array([0.0] + [(2 * l + 1) * (factorial(l - 1) / factorial(l + 1))**2
                            for l in range(1, LMAX + 1)])
    pref = -R * (R / 2.0)**5 * (2 * np.pi)**2   # = -(R^2/2)(R/2)^3 (R/2)^2 (2pi)^2 (2/R)
    V = np.zeros((len(p0_basis), len(p1_basis)))
    for i, bi in enumerate(p0_basis):
        for j, bj in enumerate(p1_basis):
            base = mul(sym_product(bi, bj), P_VNE_TERMS)
            mono0, mono1 = mul(base, A_TERMS), base
            tot = 0.0
            for (P1, Q1, P2, Q2), c in mono0.items():
                for l in range(LMAX + 1):
                    cc = Cl(l, Q1) * Cl(l, Q2)
                    if cc:
                        tot += (2 * l + 1) * c * X0[l, P1, P2] * cc
            for (P1, Q1, P2, Q2), c in mono1.items():
                for l in range(1, LMAX + 1):
                    dd = Dl(l, Q1) * Dl(l, Q2)
                    if dd:
                        tot += 2.0 * cl1[l] * c * X1[l, P1, P2] * dd
            V[i, j] = pref * tot
    return V


def vne_quadrature(p0_basis, p1_basis):
    g = build_quadrature_grids(N_xi=28, N_eta=20, N_phi=32, xi_max=16.0)
    xi, wxi, eta, weta = g['xi'], g['w_xi'], g['eta'], g['w_eta']
    dphi, wphi = g['dphi'], g['w_phi']
    hR = R / 2.0
    V = np.zeros((len(p0_basis), len(p1_basis)))

    def sv(bf, x1, e1, x2, e2):
        return sum(x1**j * x2**k * e1**l * e2**m
                   for (j, k, l, m) in _get_unsym_terms(bf))

    for a in range(len(xi)):
        x1 = xi[a]
        for c2 in range(len(xi)):
            x2 = xi[c2]
            ef = np.exp(-C * (x1 + x2))
            for b in range(len(eta)):
                e1 = eta[b]; Jp1 = x1**2 - e1**2
                r1A = hR * (x1 + e1); r1B = hR * abs(x1 - e1)
                for d in range(len(eta)):
                    e2 = eta[d]; Jp2 = x2**2 - e2**2
                    r2A = hR * (x2 + e2); r2B = hR * abs(x2 - e2)
                    vne = -(1 / r1A + 1 / r1B + 1 / r2A + 1 / r2B)
                    rho1 = np.sqrt(max((x1**2 - 1) * (1 - e1**2), 0.0))
                    rho2 = np.sqrt(max((x2**2 - 1) * (1 - e2**2), 0.0))
                    A = (x1 * e1 - x2 * e2)**2 + rho1**2 + rho2**2
                    Bc = 2 * rho1 * rho2
                    r12 = hR * np.sqrt(np.maximum(A - Bc * np.cos(dphi), 0.0))
                    rbar = np.sum(r12 * wphi) / (2 * np.pi)
                    wgt = wxi[a] * wxi[c2] * weta[b] * weta[d] * ef * Jp1 * Jp2
                    kern = hR**6 * (2 * np.pi)**2 * rbar * vne
                    for i, bi in enumerate(p0_basis):
                        si = sv(bi, x1, e1, x2, e2)
                        for j, bj in enumerate(p1_basis):
                            V[i, j] += wgt * si * sv(bj, x1, e1, x2, e2) * kern
    return V


def main():
    p0 = [HylleraasBasisFunction(j, k, l, m, 0, ALPHA)
          for (j, k, l, m) in [(0, 0, 0, 0), (1, 0, 0, 0)]]
    p1 = [HylleraasBasisFunction(j, k, l, m, 1, ALPHA)
          for (j, k, l, m) in [(0, 0, 0, 0), (1, 0, 0, 0)]]
    Vr, Vn = vne_reduction(p0, p1), vne_quadrature(p0, p1)
    rel = np.abs(Vr - Vn) / np.maximum(np.abs(Vn), 1e-12)
    print("V_ne p0xp1 reduction:\n", np.round(Vr, 6))
    print("V_ne p0xp1 quadrature:\n", np.round(Vn, 6))
    print(f"max rel diff = {np.max(rel):.3e}")
    print("HOMONUCLEAR V_ne p0xp1", "VALIDATED" if np.max(rel) < 2e-3 else "MISMATCH")


if __name__ == "__main__":
    main()
