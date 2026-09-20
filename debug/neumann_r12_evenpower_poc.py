"""H2 Neumann-r12 build, increment 1: the EVEN-power path, validated.

See debug/sprint_neumann_r12_build_memo.md. Overlap of two p=1 Hylleraas
functions needs <g_i g_j r12^2>. Since
    r12^2 = (R/2)^2 (A - B cos(phi1-phi2)),
    A = (xi1 eta1 - xi2 eta2)^2 + (xi1^2-1)(1-eta1^2) + (xi2^2-1)(1-eta2^2),
the azimuthal integral kills the B term, leaving a PURE POLYNOMIAL moment
(R/2)^2 (2pi)^2 <g_i g_j A> -- no Neumann sum. This computes that overlap
ALGEBRAICALLY (products of xi-moments A_n(2a) and eta-moments 2/(q+1)) and checks
it against the crude 5D-quadrature compute_overlap_matrix for a PURE p=1 basis
(every pair is r12^2 = even).

Result (2026-09-19): max rel diff 5.4e-6 vs a grid-limited ground truth
(algebraic side exact). EVEN-POWER PATH VALIDATED.
"""
import numpy as np
from geovac.hylleraas import (
    HylleraasBasisFunction, build_quadrature_grids, compute_overlap_matrix,
)
from geovac.neumann_vee import compute_An_table, _get_unsym_terms

R = 1.4011
ALPHA = 1.0

# A as separable terms: (coeff, e_xi1, e_eta1, e_xi2, e_eta2)
A_TERMS = [
    (1.0, 2, 2, 0, 0), (-2.0, 1, 1, 1, 1), (1.0, 0, 0, 2, 2),          # (xi1 eta1 - xi2 eta2)^2
    (1.0, 2, 0, 0, 0), (-1.0, 2, 2, 0, 0), (-1.0, 0, 0, 0, 0), (1.0, 0, 2, 0, 0),  # (xi1^2-1)(1-eta1^2)
    (1.0, 0, 0, 2, 0), (-1.0, 0, 0, 2, 2), (-1.0, 0, 0, 0, 0), (1.0, 0, 0, 0, 2),  # (xi2^2-1)(1-eta2^2)
]
J1_TERMS = [(1.0, 2, 0), (-1.0, 0, 2)]   # xi1^2 - eta1^2
J2_TERMS = [(1.0, 2, 0), (-1.0, 0, 2)]   # xi2^2 - eta2^2


def eta_moment(q):
    return 2.0 / (q + 1) if q % 2 == 0 else 0.0


def algebraic_overlap_p1(bf_i, bf_j, R, alpha):
    c = 2.0 * alpha
    An = compute_An_table(40, c)   # A_n(2 alpha) = int_1^inf xi^n e^{-2a xi} dxi
    prod_monos = {}
    for (ji, ki, li, mi) in _get_unsym_terms(bf_i):
        for (jj, kj, lj, mj) in _get_unsym_terms(bf_j):
            key = (ji + jj, li + lj, ki + kj, mi + mj)
            prod_monos[key] = prod_monos.get(key, 0.0) + 1.0
    total = 0.0
    for (sx1, se1, sx2, se2), scoef in prod_monos.items():
        for (ac, ax1, ae1, ax2, ae2) in A_TERMS:
            for (j1c, jx1, je1) in J1_TERMS:
                for (j2c, jx2, je2) in J2_TERMS:
                    total += (scoef * ac * j1c * j2c
                              * An[sx1 + ax1 + jx1] * eta_moment(se1 + ae1 + je1)
                              * An[sx2 + ax2 + jx2] * eta_moment(se2 + ae2 + je2))
    half_R = R / 2.0
    return half_R**8 * (2 * np.pi)**2 * total


def main():
    specs = [(0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 0, 0), (2, 0, 0, 0), (0, 0, 1, 1)]
    basis = [HylleraasBasisFunction(j, k, l, m, 1, ALPHA) for (j, k, l, m) in specs]
    grids = build_quadrature_grids(N_xi=26, N_eta=18, N_phi=24, xi_max=14.0)
    S_num = compute_overlap_matrix(basis, R, grids)
    n = len(basis)
    S_alg = np.array([[algebraic_overlap_p1(basis[i], basis[j], R, ALPHA)
                       for j in range(n)] for i in range(n)])
    rel = np.abs(S_alg - S_num) / np.maximum(np.abs(S_num), 1e-12)
    print("basis (all p=1):", specs)
    print(f"max abs diff = {np.max(np.abs(S_alg - S_num)):.3e}")
    print(f"max rel diff = {np.max(rel):.3e}   (grid-limited ground truth)")
    print("EVEN-POWER PATH", "VALIDATED" if np.max(rel) < 1e-4 else "MISMATCH")


if __name__ == "__main__":
    main()
