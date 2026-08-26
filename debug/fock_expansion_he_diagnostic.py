"""Fock-expansion diagnostic on helium (2026-08-19).

Question (per PI directive 'let's try the Fock expansion'): in the Avery/continuous
pivot (sparsity no longer the constraint that killed the TC attempts), is the FOCK
EXPANSION a genuine accuracy lever for the electron-electron cusp -- and specifically,
is its DISTINCTIVE content (the logarithmic terms R^2 ln R, r12 ln r12, which carry the
non-analyticity that orbital products cannot represent) a meaningful energy contribution,
or does 'the Fock expansion' collapse to plain explicitly-correlated linear-r12 (Kato)
with the log part negligible?

Method: a minimal Hylleraas variational He calc in (r1, r2, r12) coordinates, comparing
incremental basis families and isolating the linear-r12 lever from the log lever.
  phi_i = r1^l r2^m r12^n [x logfactor] e^{-zeta (r1+r2)}   (S-state; symmetrized l<->m)
Kinetic energy via the symmetric GRADIENT form (robust; handles log functions via sympy):
  <T> = (1/2) int [ d1_i d1_j + d2_i d2_j + 2 d12_i d12_j
                    + (d1_i d12_j + d12_i d1_j) M1 + (d2_i d12_j + d12_i d2_j) M2 ] dV
  M1 = (r1^2 - r2^2 + r12^2)/(2 r1 r12),  M2 = (r2^2 - r1^2 + r12^2)/(2 r2 r12)
  dV = r1 r2 r12 dr1 dr2 dr12 (angular 8pi^2 dropped; cancels in the generalized eig).
Domain r12 in [|r1-r2|, r1+r2] via r12 = |r1-r2| + 2 min(r1,r2) u, u in [0,1].
Radial [0,inf) via Gauss-Laguerre with the e^{-2 zeta r} weight extracted.

Validation gate: the single 1s^2 function (l=m=n=0) at zeta=Z-5/16=1.6875 must give the
textbook E = -(Z-5/16)^2 = -2.84766 Ha exactly.
"""
from __future__ import annotations
import numpy as np
import sympy as sp
from numpy.polynomial.laguerre import laggauss

r1, r2, r12, zeta = sp.symbols('r1 r2 r12 zeta', positive=True)
Z = 2


def base_expr(l, m, n, logkind):
    """Symbolic base function (no exponential): r1^l r2^m r12^n x logfactor."""
    e = r1**l * r2**m * r12**n
    if logkind == 'ln_s':          # r12-cusp Fock log analog: ln(r1+r2)
        e = e * sp.log(r1 + r2)
    elif logkind == 'ln_R':        # hyperradial Fock log: (1/2) ln(r1^2+r2^2) = ln R
        e = e * sp.log(sp.sqrt(r1**2 + r2**2))
    elif logkind == 'ln_u':        # interelectronic log: ln(r12)  (the r12 ln r12 Fock term)
        e = e * sp.log(r12)
    return e


def build_integrands(terms):
    """terms: list of (l,m,n,logkind). Returns lambdified g_S,g_T,g_Vne,g_Vee for every
    pair (i,j), symmetrized in l<->m."""
    bases = []
    for (l, m, n, lk) in terms:
        b = base_expr(l, m, n, lk) + base_expr(m, l, n, lk)     # S-state symmetrization
        bases.append(sp.expand(b))
    # reduced derivatives (exponential factored out): d1 = dB/dr1 - zeta B, etc.
    d1 = [sp.diff(b, r1) - zeta * b for b in bases]
    d2 = [sp.diff(b, r2) - zeta * b for b in bases]
    d12 = [sp.diff(b, r12) for b in bases]
    M1 = (r1**2 - r2**2 + r12**2) / (2 * r1 * r12)
    M2 = (r2**2 - r1**2 + r12**2) / (2 * r2 * r12)
    vol = r1 * r2 * r12
    N = len(bases)
    gS = [[None]*N for _ in range(N)]
    gT = [[None]*N for _ in range(N)]
    gVne = [[None]*N for _ in range(N)]
    gVee = [[None]*N for _ in range(N)]
    args = (r1, r2, r12, zeta)
    for i in range(N):
        for j in range(i, N):
            S = bases[i]*bases[j]*vol
            T = sp.Rational(1, 2)*(d1[i]*d1[j] + d2[i]*d2[j] + 2*d12[i]*d12[j]
                                   + (d1[i]*d12[j] + d12[i]*d1[j])*M1
                                   + (d2[i]*d12[j] + d12[i]*d2[j])*M2)*vol
            Vne = bases[i]*bases[j]*(-Z*(1/r1 + 1/r2))*vol
            Vee = bases[i]*bases[j]*(1/r12)*vol
            gS[i][j] = sp.lambdify(args, S, 'numpy')
            gT[i][j] = sp.lambdify(args, T, 'numpy')
            gVne[i][j] = sp.lambdify(args, Vne, 'numpy')
            gVee[i][j] = sp.lambdify(args, Vee, 'numpy')
    return N, gS, gT, gVne, gVee


def grid(zval, nlag=48):
    """PERIMETRIC coordinates u,v,w in [0,inf) (no |r1-r2| kink): r1=(v+w)/2, r2=(u+w)/2,
    r12=(u+v)/2, Jacobian 1/4, and e^{-2zeta(r1+r2)} = e^{-zeta(u+v+2w)}.  Tensor
    Gauss-Laguerre extracts the exponential in each of u,v (rate zeta) and w (rate 2zeta)."""
    xa, wa = laggauss(nlag)                 # int_0^inf e^{-x} f dx ~ sum wa f(xa)
    uu = xa/zval; vv = xa/zval; ww = xa/(2*zval)
    A, B, C = np.meshgrid(np.arange(nlag), np.arange(nlag), np.arange(nlag), indexing='ij')
    ug = uu[A]; vg = vv[B]; wg = ww[C]
    r1g = (vg + wg)/2; r2g = (ug + wg)/2; r12g = (ug + vg)/2
    # weight: (1/zeta)(1/zeta)(1/2zeta) wa wa wa * Jacobian(1/4); exponential is the GL weight
    w = (wa[A]*wa[B]*wa[C]) / (zval*zval*2*zval) * 0.25
    return r1g.ravel(), r2g.ravel(), r12g.ravel(), w.ravel()


def solve(terms, zval, nlag=48):
    N, gS, gT, gVne, gVee = build_integrands(terms)
    r1g, r2g, r12g, w = grid(zval, nlag)
    zc = np.full_like(r1g, zval)
    S = np.zeros((N, N)); H = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            s = np.sum(w*gS[i][j](r1g, r2g, r12g, zc))
            t = np.sum(w*gT[i][j](r1g, r2g, r12g, zc))
            vne = np.sum(w*gVne[i][j](r1g, r2g, r12g, zc))
            vee = np.sum(w*gVee[i][j](r1g, r2g, r12g, zc))
            S[i, j] = S[j, i] = s
            H[i, j] = H[j, i] = t + vne + vee
    from scipy.linalg import eigh
    evals = eigh(H, S, eigvals_only=True)
    return float(evals[0]), N


if __name__ == '__main__':
    EXACT = -2.9037243770
    print("VALIDATION: single 1s^2, zeta=1.6875 -> should be -2.847656")
    e0, _ = solve([(0, 0, 0, None)], 1.6875)
    print(f"   E = {e0:.6f}   (target -2.847656, err {abs(e0+2.847656):.2e})")

    zval = 1.8155                       # near-optimal single exponent for correlated He
    print(f"\nFAMILY COMPARISON (single zeta={zval}); exact = {EXACT:.6f}")

    # A: analytic orbital-product-like (NO r12)
    A = [(0,0,0,None),(1,0,0,None),(2,0,0,None),(1,1,0,None),(2,1,0,None),(2,2,0,None)]
    eA, nA = solve(A, zval)
    print(f"  A  analytic (no r12),        {nA} terms: E = {eA:.6f}  err {abs(eA-EXACT)*1e3:8.3f} mHa")

    # B: + linear r12 (the Kato cusp)
    B = A + [(0,0,1,None),(1,0,1,None),(1,1,1,None)]
    eB, nB = solve(B, zval)
    print(f"  B  A + linear r12 (Kato),    {nB} terms: E = {eB:.6f}  err {abs(eB-EXACT)*1e3:8.3f} mHa")

    # C: + Fock log terms (the DISTINCTIVE non-analytic content)
    C_s = B + [(0,0,1,'ln_s')]
    eC1, nC1 = solve(C_s, zval)
    print(f"  C1 B + r12*ln(r1+r2),        {nC1} terms: E = {eC1:.6f}  err {abs(eC1-EXACT)*1e3:8.3f} mHa")
    C_u = B + [(0,0,1,'ln_u')]
    eC2, nC2 = solve(C_u, zval)
    print(f"  C2 B + r12*ln(r12),          {nC2} terms: E = {eC2:.6f}  err {abs(eC2-EXACT)*1e3:8.3f} mHa")
    C_R = B + [(0,0,0,'ln_R'),(2,0,0,'ln_R')]
    eC3, nC3 = solve(C_R, zval)
    print(f"  C3 B + {{1,r1^2}}*lnR (Fock), {nC3} terms: E = {eC3:.6f}  err {abs(eC3-EXACT)*1e3:8.3f} mHa")

    print("\nLEVERS:")
    print(f"  linear r12 (Kato)  A->B : {(eA-eB)*1e3:8.3f} mHa   <-- the cusp lever")
    print(f"  Fock log r12 ln s  B->C1: {(eB-eC1)*1e6:8.3f} uHa")
    print(f"  Fock log r12 ln u  B->C2: {(eB-eC2)*1e6:8.3f} uHa")
    print(f"  Fock log R^2 ln R  B->C3: {(eB-eC3)*1e6:8.3f} uHa")
