"""H2 Neumann-r12 build, increment 7: fully-ALGEBRAIC S + V_ee + V_ne core
(mono-moments + ordered-xi X_l + analytic eta moments; grid-free, scalable).

Consolidates the validated blocks:
  S:   r12^0 even (poly mom), r12^1 odd (A K0 - B K1), r12^2 even (A poly mom)
  Vee: r12^-1 odd k=0 (K0), r12^0 even (=<g g>), r12^1 odd k=1
  Vne: homonuclear poly P_Vne in place of Jacobian; even (r12^0,r12^2) / odd (r12^1)
Kinetic is NOT included here (separate; p0p0/p1p1 even-poly + p0p1 odd, owed as the
last algebraic block).  Purpose of this module: a fast grid-free core to (a)
cross-check S,Vee,Vne against the grid-converged crude engine and (b) measure
cond(S) vs basis size -- the conditioning wall that gates physical accuracy.

Run from debug/ (sibling import of the oddpower primitives).
"""
import time
import numpy as np
from math import factorial
from geovac.hylleraas import HylleraasBasisFunction, generate_basis
from geovac.neumann_vee import _get_unsym_terms, compute_An_table
from neumann_r12_oddpower_vee_poc import (
    mul, sym_product, Cl, Dl, build_X, _xi_grid, A_TERMS, JP1, JP2, R, ALPHA, C, LMAX,
)
from neumann_r12_vne_p0p1_poc import P_VNE_TERMS

hR = R / 2.0
_AN = None
_X0 = None
_X1 = None
_CL1 = None


def _init(nmom=60):
    global _AN, _X0, _X1, _CL1
    _AN = compute_An_table(nmom, C)            # A_n(2 alpha)
    XI, WXI = _xi_grid()
    _X0 = build_X(0, XI, WXI)
    _X1 = build_X(1, XI, WXI)
    _CL1 = np.array([0.0] + [(2 * l + 1) * (factorial(l - 1) / factorial(l + 1))**2
                             for l in range(1, LMAX + 1)])


def _M(q):
    return 2.0 / (q + 1) if q % 2 == 0 else 0.0


def even_mom(mono):
    """sum coef A_{P1}(2a) A_{P2}(2a) M(Q1) M(Q2)."""
    s = 0.0
    for (P1, Q1, P2, Q2), c in mono.items():
        if c == 0.0:
            continue
        m = _M(Q1) * _M(Q2)
        if m != 0.0 and P1 < len(_AN) and P2 < len(_AN):
            s += c * _AN[P1] * _AN[P2] * m
    return s


def odd_mom(base, k):
    """k=0: <1/r12>-type (K0 only).  k=1: <r12>-type (A K0 - B K1)."""
    tot = 0.0
    mono0 = mul(base, A_TERMS) if k == 1 else base
    for (P1, Q1, P2, Q2), c in mono0.items():
        if c == 0.0:
            continue
        for l in range(LMAX + 1):
            cc = Cl(l, Q1) * Cl(l, Q2)
            if cc != 0.0 and P1 < _X0.shape[1] and P2 < _X0.shape[2]:
                tot += (2 * l + 1) * c * _X0[l, P1, P2] * cc
    if k == 1:
        for (P1, Q1, P2, Q2), c in base.items():
            if c == 0.0:
                continue
            for l in range(1, LMAX + 1):
                dd = Dl(l, Q1) * Dl(l, Q2)
                if dd != 0.0 and P1 < _X1.shape[1] and P2 < _X1.shape[2]:
                    tot += 2.0 * _CL1[l] * c * _X1[l, P1, P2] * dd
    return tot


def build_S_Vee_Vne(basis):
    n = len(basis)
    S = np.zeros((n, n)); Vee = np.zeros((n, n)); Vne = np.zeros((n, n))
    pref_even = hR**6 * (2 * np.pi)**2               # r12^0 even
    pref_even2 = hR**8 * (2 * np.pi)**2              # r12^2 even (extra (R/2)^2)
    pref_odd = hR**8 * (2 * np.pi)**2 * (2.0 / R)    # r12^1 odd
    pref_vee0 = hR**6 * (2 * np.pi)**2 * (2.0 / R)   # r12^-1 (K0) = pi^2 R^5/8
    pv_even = -(R**2 / 2) * hR**3 * (2 * np.pi)**2       # Vne r12^0
    pv_even2 = -(R**2 / 2) * hR**5 * (2 * np.pi)**2      # Vne r12^2
    pv_odd = -R * hR**5 * (2 * np.pi)**2                 # Vne r12^1
    for i in range(n):
        for j in range(i, n):
            sym = sym_product(basis[i], basis[j])
            jac = mul(mul(sym, JP1), JP2)
            pvn = mul(sym, P_VNE_TERMS)
            Pij = basis[i].p + basis[j].p
            if Pij == 0:
                S[i, j] = pref_even * even_mom(jac)
                Vee[i, j] = pref_vee0 * odd_mom(jac, 0)
                Vne[i, j] = pv_even * even_mom(pvn)
            elif Pij == 1:
                S[i, j] = pref_odd * odd_mom(jac, 1)
                Vee[i, j] = pref_even * even_mom(jac)          # r12^0 = <g g>
                Vne[i, j] = pv_odd * odd_mom(pvn, 1)
            else:  # Pij == 2
                S[i, j] = pref_even2 * even_mom(mul(jac, A_TERMS))
                Vee[i, j] = pref_odd * odd_mom(jac, 1)          # r12^1
                Vne[i, j] = pv_even2 * even_mom(mul(pvn, A_TERMS))
            S[j, i] = S[i, j]; Vee[j, i] = Vee[i, j]; Vne[j, i] = Vne[i, j]
    return S, Vee, Vne


def main():
    _init()
    # --- validate S,Vee,Vne vs grid-converged crude (small basis) ---
    from geovac.hylleraas import build_quadrature_grids, compute_overlap_matrix, compute_hamiltonian_matrix
    basis = generate_basis(j_max=1, l_max=0, p_max=1, alpha=ALPHA)
    S, Vee, Vne = build_S_Vee_Vne(basis)
    gc = build_quadrature_grids(N_xi=90, N_eta=60, N_phi=64, xi_max=16.0)
    S_cr = compute_overlap_matrix(basis, R, gc)
    print(f"S vs crude(nx=90): max rel = {np.max(np.abs(S - S_cr)/np.maximum(np.abs(S_cr),1e-10)):.2e}")

    # --- cond(S) vs basis size (the conditioning wall) ---
    print("\ncond(S) vs basis (float64):")
    for (jm, lm, pm) in [(1, 0, 1), (2, 1, 1), (3, 2, 1), (4, 2, 1), (4, 3, 1)]:
        b = generate_basis(j_max=jm, l_max=lm, p_max=pm, alpha=ALPHA)
        t0 = time.perf_counter()
        Sb, _, _ = build_S_Vee_Vne(b)
        dt = time.perf_counter() - t0
        cond = np.linalg.cond(Sb)
        print(f"  (j={jm},l={lm},p<=1) N={len(b):3d}  cond(S)={cond:.2e}  build {dt:.1f}s", flush=True)


if __name__ == "__main__":
    main()
