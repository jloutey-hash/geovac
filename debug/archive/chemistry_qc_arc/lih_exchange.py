"""Increment 4: core-valence EXCHANGE for frozen-core prolate LiH.

For a closed-shell 1s^2 core, each valence electron sees the core Fock potential
2J - K:  2J = V_H (the Hartree screening, increment 1, already in);  -K is the
non-local exchange with the SAME-spin core electron. K is the two-electron integral

    K_ij = <g_i(1) 1s_A(2) | 1/r12 | 1s_A(1) g_j(2)>
         = INT INT [g_i(r1) 1s_A(r1)] (1/r12) [1s_A(r2) g_j(r2)] dr1 dr2,

the Coulomb interaction of the overlap densities g_i*1s and 1s*g_j. In prolate
coords with sigma (mu=0) basis, the phi-integral of 1/r12 is done ANALYTICALLY:
  INT_0^2pi INT_0^2pi dphi1 dphi2 / r12 = (2pi/hR) * 4/sqrt(A+B) * K(sqrt(2B/(A+B)))
(complete elliptic K), which captures the integrable 1/r12 singularity (log as the
electrons coincide). The remaining (xi,eta)x(xi,eta) quadrature is a matrix triple
product K1e = F Gphi F^T with F[i,a] = g_i(a) 1s(a) w(a).

Adds -K (per valence electron) to the valence one-body Hamiltonian:
  H_val -> H_val + V_H(1)+V_H(2) - K(1)-K(2),  K2e = K1e (x) S1 + S1 (x) K1e.
Expect: E lowers from the Hartree-only -8.012 TOWARD exact -8.070 (exchange binds),
by a sane amount (tens of mHa); if it overshoots below exact, sign/coeff is wrong.

FIRST LOOK: moderate grid; the 1/r12 near-singularity makes this less accurate than
the closed-form pieces -- magnitude + sign + trend are the deliverable.

Run from root:  python debug/lih_exchange.py
"""
import os
import sys
import time
import numpy as np
from scipy.special import ellipk

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from geovac import prolate_recondition as pr                 # noqa: E402
from geovac.hylleraas import build_quadrature_grids          # noqa: E402
import prolate_r12_mpf as m                                  # noqa: E402
from heh_probe import build_basis_full                       # noqa: E402
from r12ci_first_energy import solve_canonical               # noqa: E402
from prolate_core_hartree import V_H_closed, ZC_LI           # noqa: E402
import lih_frozen_core_first as L                            # noqa: E402
from lih_projector_test import build, _one_elec_pieces       # noqa: E402

R = 3.015
ALPHA = 1.0
JMAX, LMAX = 2, 2


def compute_K1e(jl_pairs, R, alpha, grid):
    """1-electron exchange matrix K1e[i,j] over the (j,l) index set (sigma)."""
    xi, wxi = grid['xi'], grid['w_xi']
    eta, weta = grid['eta'], grid['w_eta']
    hR = R / 2.0
    # flatten the (xi,eta) grid
    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    WXI, WETA = np.meshgrid(wxi, weta, indexing='ij')
    xg = XI.ravel(); eg = ETA.ravel()
    wg = (WXI * WETA).ravel() * (xg**2 - eg**2) * (hR**3)     # per-electron weight
    rho = np.sqrt(np.clip((xg**2 - 1.0) * (1.0 - eg**2), 0.0, None))
    ng = xg.size
    # Gphi[a,b] = (2pi/hR) * 4/sqrt(A+B) * K(2B/(A+B)),  A=(xi_a eta_a - xi_b eta_b)^2
    #             + rho_a^2 + rho_b^2 ; B = 2 rho_a rho_b
    xe = xg * eg
    A = (xe[:, None] - xe[None, :])**2 + rho[:, None]**2 + rho[None, :]**2
    B = 2.0 * rho[:, None] * rho[None, :]
    mparam = np.clip(2.0 * B / (A + B + 1e-300), 0.0, 1.0 - 1e-9)
    Gphi = (2.0 * np.pi / hR) * (4.0 / np.sqrt(A + B + 1e-300)) * ellipk(mparam)
    # F[i,a] = g_i(a) * 1s(a) * w(a)
    n = len(jl_pairs)
    F = np.empty((n, ng))
    s1s = np.exp(-ZC_LI * hR * (xg + eg))                     # 1s_A on the grid
    base = np.exp(-alpha * xg) * s1s * wg
    for i, (ji, li) in enumerate(jl_pairs):
        F[i] = (xg**ji) * (eg**li) * base
    return F @ Gphi @ F.T


def _jl_and_S(basis_p, R, alpha):
    idx = {}
    for (b, p) in basis_p:
        idx.setdefault((b.j, b.l), None); idx.setdefault((b.k, b.m), None)
    jl = list(idx.keys()); pos = {p: i for i, p in enumerate(jl)}
    grid = build_quadrature_grids(N_xi=40, N_eta=30, N_phi=4, xi_max=20.0)
    _, S1, _ = _one_elec_pieces(jl, R, alpha, grid)
    return jl, pos, S1, grid


def K2e(basis_p, jl, pos, K1, S1):
    n = len(basis_p)
    W = np.zeros((n, n))
    for i, (bi, _) in enumerate(basis_p):
        a1, b1 = pos[(bi.j, bi.l)], pos[(bi.k, bi.m)]
        for jj, (bj, _) in enumerate(basis_p):
            a2, b2 = pos[(bj.j, bj.l)], pos[(bj.k, bj.m)]
            W[i, jj] = K1[a1, a2] * S1[b1, b2] + S1[a1, a2] * K1[b1, b2]
    return W


def main():
    t0 = time.time()
    basis = build_basis_full(JMAX, LMAX, ALPHA, p_set=(0,))
    S, H = m.assemble_hetero(basis, R, ALPHA, L.Z_LI, L.Z_H, l_neumann=16, dps=30)
    jl, pos, S1, grid = _jl_and_S(basis, R, ALPHA)
    VH2, Plam = build(basis, R, ALPHA, 1000.0)
    K1 = compute_K1e(jl, R, ALPHA, grid)
    Kx = K2e(basis, jl, pos, K1, S1)
    print(f"LiH exchange increment, R={R}, (j,l)=({JMAX},{LMAX})   [{time.time()-t0:.0f}s]")
    print(f"  K1e symmetric: {np.allclose(K1, K1.T)}; K1e diag (should be >0): "
          f"{np.round(np.diag(K1)[:4], 4)}")

    def E_of(Hmat):
        e_val = solve_canonical(S, Hmat)[0]
        return e_val + L.Z_LI * L.Z_H / R - V_H_closed(R, ZC_LI) + L.E_CORE, e_val

    E_hartree, ev_h = E_of(H + VH2 + Plam)
    E_exch, ev_x = E_of(H + VH2 + Plam - Kx)
    print(f"  Hartree-only:   E_val={ev_h:.5f}  E_tot={E_hartree:.5f}")
    print(f"  + exchange(-K): E_val={ev_x:.5f}  E_tot={E_exch:.5f}  "
          f"(dE = {(E_exch-E_hartree)*1e3:+.2f} mHa)")
    print(f"  exact LiH ~ -8.070;  variational check: E_tot > exact? "
          f"{'YES' if E_exch > -8.070 else 'NO (overshoot -> sign/coeff bug)'}")


if __name__ == '__main__':
    main()
