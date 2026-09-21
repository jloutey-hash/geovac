"""Increment 2b: the Huzinaga core-orthogonality projector for frozen-core LiH.

The first look (lih_frozen_core_first.py) showed the valence over-binds 0.6 Ha
(E_val -1.73 vs expected -1.12) by collapsing into the -3/r_A core well without a
keep-out constraint. Fix: add a level-shift H -> H + lambda(P_1s(1)+P_1s(2)),
P_1s = |1s_Li><1s_Li| (normalized), which pushes any valence component overlapping
the core 1s up by ~lambda, forcing orthogonality as lambda -> large.

|1s_Li>(r_A) = exp(-zc*r_A), r_A=(R/2)(xi+eta), zc = He-like Li^2+ 1s exponent.
c_i = <g_i|1s_Li>/sqrt(<1s|1s>) in the prolate sigma basis (quadrature).
2e projector (factorized, p=0): lambda[ c_a1 c_a2 S1_{b1 b2} + S1_{a1 a2} c_b1 c_b2 ].

TEST: E_val(lambda) at R=3.015 should rise from -1.73 (collapsed) and PLATEAU near
the expected ~-1.12 as lambda grows (collapse removed, orthogonality enforced).

Run from root:  python debug/lih_projector_test.py
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from geovac import prolate_recondition as pr                 # noqa: E402
from geovac.hylleraas import build_quadrature_grids          # noqa: E402
import prolate_r12_mpf as m                                  # noqa: E402
from heh_probe import build_basis_full                       # noqa: E402
from r12ci_first_energy import solve_canonical               # noqa: E402
from prolate_core_hartree import V_H_closed, ZC_LI           # noqa: E402
import lih_frozen_core_first as L                            # noqa: E402

R = 3.015
ALPHA = 1.0
JMAX, LMAX = 2, 2


def _one_elec_pieces(jl_pairs, R, alpha, grid):
    """1e V_H, overlap S, and c_i=<g_i|1s_Li> (unnormalized) + <1s|1s>."""
    xi, wxi = grid['xi'], grid['w_xi']
    eta, weta = grid['eta'], grid['w_eta']
    hR = R / 2.0
    n = len(jl_pairs)
    VH = np.zeros((n, n)); S = np.zeros((n, n)); c = np.zeros(n); norm1s = 0.0
    for a in range(len(xi)):
        x = xi[a]; efa = np.exp(-2.0 * alpha * x)
        s1s_x = np.exp(-ZC_LI * hR * x)                     # 1s exp in xi
        for b in range(len(eta)):
            e = eta[b]; r_A = hR * (x + e)
            if r_A <= 0:
                continue
            jac = (x**2 - e**2)
            wgeom = wxi[a] * weta[b] * jac * (hR**3) * 2.0 * np.pi
            vh = V_H_closed(r_A, ZC_LI)
            phi1s = np.exp(-ZC_LI * r_A)                     # = s1s_x * exp(-zc*hR*e)
            norm1s += wgeom * phi1s * phi1s
            for i, (ji, li) in enumerate(jl_pairs):
                gi = x**ji * e**li
                # <g_i|1s>: g_i carries e^{-alpha xi}; 1s carries e^{-zc r_A}
                c[i] += wgeom * (gi * np.exp(-alpha * x)) * phi1s
                for jj, (jx, lx) in enumerate(jl_pairs):
                    gj = x**jx * e**lx
                    wg = wgeom * efa * gi * gj
                    S[i, jj] += wg
                    VH[i, jj] += wg * vh
    return VH, S, c / np.sqrt(norm1s)


def build(basis_p, R, alpha, lam):
    idx = {}
    for (b, p) in basis_p:
        idx.setdefault((b.j, b.l), None); idx.setdefault((b.k, b.m), None)
    jl = list(idx.keys()); pos = {p: i for i, p in enumerate(jl)}
    grid = build_quadrature_grids(N_xi=40, N_eta=30, N_phi=4, xi_max=20.0)
    VH1, S1, cn = _one_elec_pieces(jl, R, alpha, grid)
    n = len(basis_p)
    W = np.zeros((n, n)); P = np.zeros((n, n))
    for i, (bi, _) in enumerate(basis_p):
        a1, b1 = pos[(bi.j, bi.l)], pos[(bi.k, bi.m)]
        for jj, (bj, _) in enumerate(basis_p):
            a2, b2 = pos[(bj.j, bj.l)], pos[(bj.k, bj.m)]
            W[i, jj] = VH1[a1, a2] * S1[b1, b2] + S1[a1, a2] * VH1[b1, b2]
            P[i, jj] = cn[a1] * cn[a2] * S1[b1, b2] + S1[a1, a2] * cn[b1] * cn[b2]
    return W, lam * P


def main():
    t0 = time.time()
    basis = build_basis_full(JMAX, LMAX, ALPHA, p_set=(0,))
    S, H = m.assemble_hetero(basis, R, ALPHA, L.Z_LI, L.Z_H, l_neumann=16, dps=30)
    print(f"Huzinaga projector test, LiH R={R}, (j,l)=({JMAX},{LMAX}), alpha={ALPHA}")
    print(f"  expected valence E_val ~ -1.122 (gives E_tot = -8.070 exact)")
    print(f"  {'lambda':>8}  {'E_val':>10}  {'E_tot':>10}")
    for lam in (0.0, 1.0, 10.0, 100.0, 1000.0, 1e4, 1e5):
        VH2, Plam = build(basis, R, ALPHA, lam)
        e_val = solve_canonical(S, H + VH2 + Plam)[0]
        E_tot = e_val + L.Z_LI * L.Z_H / R - V_H_closed(R, ZC_LI) + L.E_CORE
        print(f"  {lam:8.0f}  {e_val:10.5f}  {E_tot:10.5f}   [{time.time()-t0:.0f}s]",
              flush=True)


if __name__ == '__main__':
    main()
