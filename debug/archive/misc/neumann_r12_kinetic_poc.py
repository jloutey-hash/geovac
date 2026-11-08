"""H2 Neumann-r12 build, increment 3: the KINETIC energy for p=1 collapses to the
EVEN path. See debug/sprint_neumann_r12_build_memo.md.

Claim (Cartesian Green form + IBP): for phi = g r12 (p=1),
  <phi_i|T|phi_j> = (1/2)∫(grad1 g_i·grad1 g_j + grad2 g_i·grad2 g_j) r12^2 dV
                    - 2 <g_i|g_j>_{p=0}
Both terms r12-EVEN. Derivation: grad(g r12) = (grad g) r12 + g (grad r12);
|grad_k r12|^2 = 1 exactly; r12 grad_k r12 = +-(r1-r2); the vector cross term
= (1/2)∫(r1-r2)·(grad1 - grad2)(g_i g_j) dV, and IBP with div1(r1-r2)=3,
div2(r1-r2)=-3 collapses it to -3∫g_i g_j; the |grad r12|^2 term adds +∫g_i g_j;
total -2∫g_i g_j. So the "James-Coolidge vector terms" cancel -- the kinetic is
NOT a hard piece.

Validates T_formula = KE_g_r2 - 2 OV_g against ground_T (crude analytical
prolate-Green kinetic of phi = g r12). Result (2026-09-19): max rel diff 9.4e-5.
"""
import numpy as np
from geovac.hylleraas import (
    HylleraasBasisFunction, build_quadrature_grids,
    evaluate_basis_and_derivs, _full_derivs_python, _r12_and_derivs_python,
)

R = 1.4011
ALPHA = 1.0
half_R = R / 2.0


def ground_T(basis, g):
    xi, wxi, eta, weta = g['xi'], g['w_xi'], g['eta'], g['w_eta']
    dphi, wphi = g['dphi'], g['w_phi']
    n = len(basis); T = np.zeros((n, n)); Tpref = R / 4.0
    for a in range(len(xi)):
        x1 = xi[a]; x1m1 = x1**2 - 1
        for c in range(len(xi)):
            x2 = xi[c]; x2m1 = x2**2 - 1; wac = wxi[a] * wxi[c]
            for b in range(len(eta)):
                e1 = eta[b]; om1 = 1 - e1**2; J1 = half_R**3 * (x1**2 - e1**2)
                for d in range(len(eta)):
                    e2 = eta[d]; om2 = 1 - e2**2; J2 = half_R**3 * (x2**2 - e2**2)
                    den1 = x1m1 * om1; den2 = x2m1 * om2
                    pf1 = (x1**2 - e1**2) / den1 if den1 > 1e-30 else 0.0
                    pf2 = (x2**2 - e2**2) / den2 if den2 > 1e-30 else 0.0
                    Tacc = np.zeros((n, n))
                    for ip in range(len(dphi)):
                        r12, drx1, dre1, drx2, dre2, drdp = _r12_and_derivs_python(
                            x1, e1, x2, e2, dphi[ip], R)
                        r12 = max(r12, 1e-15)
                        D = [_full_derivs_python(bf, x1, e1, x2, e2, r12,
                                                 drx1, dre1, drx2, dre2, drdp)
                             for bf in basis]
                        for i in range(n):
                            _, dix1, die1, dix2, die2, didp = D[i]
                            for j in range(n):
                                _, djx1, dje1, djx2, dje2, djdp = D[j]
                                T1 = (x1m1 * dix1 * djx1 + om1 * die1 * dje1
                                      + pf1 * didp * djdp) * J2
                                T2 = (x2m1 * dix2 * djx2 + om2 * die2 * dje2
                                      + pf2 * didp * djdp) * J1
                                Tacc[i, j] += (T1 + T2) * wphi[ip]
                    T += wac * weta[b] * weta[d] * Tpref * 2 * np.pi * Tacc
    return T


def formula_T(basis, g):
    xi, wxi, eta, weta = g['xi'], g['w_xi'], g['eta'], g['w_eta']
    n = len(basis); KE = np.zeros((n, n)); OV = np.zeros((n, n))
    Tpref = R / 4.0; phi2 = (2 * np.pi)**2
    for a in range(len(xi)):
        x1 = xi[a]; x1m1 = x1**2 - 1
        for c in range(len(xi)):
            x2 = xi[c]; x2m1 = x2**2 - 1; wac = wxi[a] * wxi[c]
            for b in range(len(eta)):
                e1 = eta[b]; om1 = 1 - e1**2; J1 = half_R**3 * (x1**2 - e1**2)
                for d in range(len(eta)):
                    e2 = eta[d]; om2 = 1 - e2**2; J2 = half_R**3 * (x2**2 - e2**2)
                    A = (x1 * e1 - x2 * e2)**2 + max(x1m1 * om1, 0.0) + max(x2m1 * om2, 0.0)
                    w = wac * weta[b] * weta[d]
                    gd = [evaluate_basis_and_derivs(bf, x1, e1, x2, e2) for bf in basis]
                    for i in range(n):
                        gi, gix1, gie1, gix2, gie2 = gd[i]
                        for j in range(n):
                            gj, gjx1, gje1, gjx2, gje2 = gd[j]
                            T1 = (x1m1 * gix1 * gjx1 + om1 * gie1 * gje1) * J2
                            T2 = (x2m1 * gix2 * gjx2 + om2 * gie2 * gje2) * J1
                            KE[i, j] += w * (T1 + T2) * A
                            OV[i, j] += w * gi * gj * J1 * J2
    return KE * Tpref * phi2 * half_R**2 - 2 * OV * phi2


def main():
    specs = [(0, 0, 0, 0), (1, 0, 0, 0)]
    basis = [HylleraasBasisFunction(j, k, l, m, 1, ALPHA) for (j, k, l, m) in specs]
    g = build_quadrature_grids(N_xi=24, N_eta=18, N_phi=28, xi_max=15.0)
    Tg = ground_T(basis, g); Tf = formula_T(basis, g)
    rel = np.abs(Tf - Tg) / np.maximum(np.abs(Tg), 1e-12)
    print("ground_T:\n", np.round(Tg, 6))
    print("formula_T:\n", np.round(Tf, 6))
    print(f"max rel diff = {np.max(rel):.3e}")
    print("KINETIC EVEN-PATH", "VALIDATED" if np.max(rel) < 5e-3 else "MISMATCH")


if __name__ == "__main__":
    main()
