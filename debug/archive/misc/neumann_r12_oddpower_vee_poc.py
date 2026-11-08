"""H2 Neumann-r12 build, increment 2: the ODD-power V_ee matrix, assembled + validated.

See debug/sprint_neumann_r12_build_memo.md. V_ee for a pure p=1 basis:
<phi_i|1/r12|phi_j> = <g_i g_j r12> (r12*(1/r12)*r12 = r12). Reduction:
  r12bar = (R/2)^2 [A K0 - B K1],
  K0 = <1/r12>_phi       = (2/R) sum_l (2l+1) P_l(xi_<)Q_l(xi_>) P_l(eta1)P_l(eta2)
  K1 = <cosDphi/r12>_phi = -(2/R) sum_l (2l+1)[(l-1)!/(l+1)!]^2 P_l^1(xi_<)Q_l^1(xi_>)P_l^1(eta1)P_l^1(eta2)

After the sqrt cancellations (B=2 rho1 rho2 pairs with the m=1 element's sqrt) the
matrix element separates into ordered-xi integrals x eta-moments:
  V_ee_ij = (R/2)^8 (2pi)^2 (2/R) [
      sum_l (2l+1) sum_mono0 c0 Xl0(P1,P2) Cl(Q1) Cl(Q2)
    + 2 sum_l cl1 sum_mono1 c1 X1l(P1,P2) Dl(Q1) Dl(Q2) ]
  mono0 = expand(Jp1 Jp2 sym_i sym_j * A),  mono1 = expand(Jp1 Jp2 sym_i sym_j)
  Xl0 = int int xi1^P1 xi2^P2 e^{-2a(xi1+xi2)} P_l(xi_<) Q_l(xi_>) dxi
  X1l = int int xi1^P1 xi2^P2 e^{-2a(xi1+xi2)} (xi1^2-1)(xi2^2-1) P_l'(xi_<)Q_l'(xi_>) dxi
  Cl(Q) = int eta^Q P_l deta,  Dl(Q) = int eta^Q (1-eta^2) P_l' deta,  cl1=(2l+1)[(l-1)!/(l+1)!]^2

Ground truth: 5D quadrature <g_i g_j r12>. The ordered-xi integrals here are 2D
Gauss (the algebraic recurrences compute_Xl / neumann_vee_general_m.build_Xtab are
already corpus-validated drop-ins); the NEW content under test is the reduction +
kernels + mono/eta assembly. Result (2026-09-19): max rel diff 9.2e-6.
"""
import numpy as np
from math import factorial
from numpy.polynomial import legendre as LEG
from scipy.special import lqmn
from geovac.hylleraas import HylleraasBasisFunction, build_quadrature_grids
from geovac.neumann_vee import _get_unsym_terms

R = 1.4011
ALPHA = 1.0
C = 2.0 * ALPHA
LMAX = 12
PMAX = 8

JP1 = [(1.0, 2, 0, 0, 0), (-1.0, 0, 2, 0, 0)]
JP2 = [(1.0, 0, 0, 2, 0), (-1.0, 0, 0, 0, 2)]
A_TERMS = [
    (1.0, 2, 2, 0, 0), (-2.0, 1, 1, 1, 1), (1.0, 0, 0, 2, 2),
    (1.0, 2, 0, 0, 0), (-1.0, 2, 2, 0, 0), (-1.0, 0, 0, 0, 0), (1.0, 0, 2, 0, 0),
    (1.0, 0, 0, 2, 0), (-1.0, 0, 0, 2, 2), (-1.0, 0, 0, 0, 0), (1.0, 0, 0, 0, 2),
]


def mul(poly, terms):
    out = {}
    for (P1, Q1, P2, Q2), c in poly.items():
        for (tc, d1, e1, d2, e2) in terms:
            k = (P1 + d1, Q1 + e1, P2 + d2, Q2 + e2)
            out[k] = out.get(k, 0.0) + c * tc
    return out


def sym_product(bf_i, bf_j):
    poly = {}
    for (ji, ki, li, mi) in _get_unsym_terms(bf_i):
        for (jj, kj, lj, mj) in _get_unsym_terms(bf_j):
            k = (ji + jj, li + lj, ki + kj, mi + mj)
            poly[k] = poly.get(k, 0.0) + 1.0
    return poly


def _leg(l):
    c = np.zeros(l + 1); c[l] = 1.0
    return c


def _integ(coef):
    return sum(2.0 * a / (k + 1) for k, a in enumerate(coef) if a != 0 and k % 2 == 0)


def Cl(l, Q):
    p = LEG.leg2poly(_leg(l))
    return _integ(np.concatenate([np.zeros(Q), p]))


def Dl(l, Q):
    if l == 0:
        return 0.0
    dp = LEG.leg2poly(LEG.legder(_leg(l)))          # P_l'
    q2 = np.concatenate([np.zeros(2), dp])          # eta^2 P_l'
    n = max(len(dp), len(q2))
    poly = np.concatenate([dp, np.zeros(n - len(dp))]) - np.concatenate([q2, np.zeros(n - len(q2))])
    return _integ(np.concatenate([np.zeros(Q), poly]))


def _xi_grid(n=90, xi_max=16.0):
    u, w = np.polynomial.legendre.leggauss(n)
    t = (u + 1) / 2
    xi = 1.0 + (xi_max - 1.0) * t**2
    return xi, w * (xi_max - 1.0) * t


def build_X(weight_pow, XI, WXI):
    ng = len(XI)
    exp = np.exp(-C * XI)
    lt = XI[:, None] <= XI[None, :]
    w2 = np.outer(WXI * exp, WXI * exp)
    Xl = np.zeros((LMAX + 1, PMAX + 1, PMAX + 1))
    for l in range(LMAX + 1):
        if weight_pow == 0:
            Pf = LEG.legval(XI, _leg(l))
            Qf = np.array([lqmn(0, l, x)[0][0, l] for x in XI])
            wa = wb = np.ones(ng)
        else:
            if l == 0:
                continue
            Pf = LEG.legval(XI, LEG.legder(_leg(l)))
            Qf = np.array([lqmn(0, l, x)[1][0, l] for x in XI])   # dQ_l/dxi
            wa = wb = (XI**2 - 1.0)
        Pcol = np.where(lt, Pf[:, None], Pf[None, :])
        Qcol = np.where(lt, Qf[None, :], Qf[:, None])
        base = w2 * (wa[:, None] * wb[None, :]) * Pcol * Qcol
        for P1 in range(PMAX + 1):
            xa = XI**P1
            for P2 in range(PMAX + 1):
                Xl[l, P1, P2] = np.sum(base * xa[:, None] * (XI**P2)[None, :])
    return Xl


def vee_reduction(basis):
    XI, WXI = _xi_grid()
    X0, X1 = build_X(0, XI, WXI), build_X(1, XI, WXI)
    cl1 = np.array([0.0] + [(2 * l + 1) * (factorial(l - 1) / factorial(l + 1))**2
                            for l in range(1, LMAX + 1)])
    pref = (R / 2.0)**8 * (2 * np.pi)**2 * (2.0 / R)
    n = len(basis)
    V = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            base = mul(mul(sym_product(basis[i], basis[j]), JP1), JP2)
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


def vee_quadrature(basis):
    g = build_quadrature_grids(N_xi=28, N_eta=20, N_phi=32, xi_max=16.0)
    xi, wxi, eta, weta = g['xi'], g['w_xi'], g['eta'], g['w_eta']
    dphi, wphi = g['dphi'], g['w_phi']
    hR = R / 2.0
    n = len(basis)
    V = np.zeros((n, n))

    def sv(bf, x1, e1, x2, e2):
        return sum(x1**j * x2**k * e1**l * e2**m
                   for (j, k, l, m) in _get_unsym_terms(bf))

    for a in range(len(xi)):
        x1 = xi[a]
        for c2 in range(len(xi)):
            x2 = xi[c2]
            ef = np.exp(-C * (x1 + x2))
            for b in range(len(eta)):
                e1 = eta[b]
                Jp1 = x1**2 - e1**2
                for d in range(len(eta)):
                    e2 = eta[d]
                    Jp2 = x2**2 - e2**2
                    rho1 = np.sqrt(max((x1**2 - 1) * (1 - e1**2), 0.0))
                    rho2 = np.sqrt(max((x2**2 - 1) * (1 - e2**2), 0.0))
                    A = (x1 * e1 - x2 * e2)**2 + rho1**2 + rho2**2
                    Bc = 2 * rho1 * rho2
                    r12 = hR * np.sqrt(np.maximum(A - Bc * np.cos(dphi), 0.0))
                    rbar = np.sum(r12 * wphi) / (2 * np.pi)
                    wgt = wxi[a] * wxi[c2] * weta[b] * weta[d] * ef * Jp1 * Jp2
                    kern = hR**6 * (2 * np.pi)**2 * rbar
                    for i in range(n):
                        si = sv(basis[i], x1, e1, x2, e2)
                        for j in range(n):
                            V[i, j] += wgt * si * sv(basis[j], x1, e1, x2, e2) * kern
    return V


def main():
    specs = [(0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 0, 0)]
    basis = [HylleraasBasisFunction(j, k, l, m, 1, ALPHA) for (j, k, l, m) in specs]
    Vr, Vn = vee_reduction(basis), vee_quadrature(basis)
    rel = np.abs(Vr - Vn) / np.maximum(np.abs(Vn), 1e-12)
    print("basis (sigma, p=1):", specs)
    print("V_red diag:", np.round(np.diag(Vr), 6))
    print("V_num diag:", np.round(np.diag(Vn), 6))
    print(f"max rel diff = {np.max(rel):.3e}")
    print("ODD-POWER V_ee", "VALIDATED" if np.max(rel) < 2e-3 else "MISMATCH")


if __name__ == "__main__":
    main()
