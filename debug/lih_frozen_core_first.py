"""Frozen-core prolate LiH — FIRST energy + R-scan (increments 2-3, Hartree-level).

The core question: does the variational-prolate (H2-style) recipe put LiH's bond
near R_e = 3.015 bohr, instead of the composed/balanced recipe's +8.8% drift?

Model (frozen Li 1s^2 core, 2 valence electrons, prolate two-center basis):
  E_tot(R) = E_val_elec(R) + E_core + V_NN(R) + E_coreH(R)
    E_val_elec = <val| T + V_ne(-3/r_A - 1/r_B) + V_H(core screening) + V_ee |val>
                 [assemble_hetero(Z_A=3, Z_B=1) gives T+V_ne+V_ee; add the 2e V_H]
    E_core     = frozen Li^2+ 1s^2 energy (He-like), R-INDEPENDENT (~ -7.28 Ha)
    V_NN(R)    = Z_Li Z_H / R = 3/R
    E_coreH(R) = 2 core electrons attracted to the H nucleus = -V_Hdens(R)
                 (potential of the 1s^2 density at distance R)
  Only E_val_elec + 3/R - V_Hdens(R) is R-dependent -> sets R_eq.

FIRST LOOK: Hartree only. NO core-valence exchange (increment 4) and NO Huzinaga
orthogonality projector yet (increment 2b) -- so watch for variational COLLAPSE of
the valence into the -3/r_A core well (E_val absurdly low). If it collapses, the
projector is required; if a modest basis stays stable, we get a first R_eq.

Run from root:  python debug/lih_frozen_core_first.py [jmax] [lmax]
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

Z_LI, Z_H = 3.0, 1.0
E_CORE = -7.2799            # He-like Li^2+ 1s^2 (R-independent), _FIRST_ROW_CORE_ENERGY[3]
R_E_EXP = 3.015            # LiH experimental R_e (bohr)
E_LIH_REF = -8.070         # LiH BO total energy near R_e (Ha, approx)

JMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2
LMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2
ALPHA = 1.0               # valence is diffuse (net Li^2+ ~ +1); small exponent
R_GRID = [2.60, 2.85, 3.015, 3.20, 3.45, 3.75]


def _one_elec_VH_and_S(jl_pairs, R, alpha, grid):
    """1-electron V_H(r_A) matrix and overlap over the (j,l) index set (sigma)."""
    xi, wxi = grid['xi'], grid['w_xi']
    eta, weta = grid['eta'], grid['w_eta']
    hR = R / 2.0
    n = len(jl_pairs)
    VH = np.zeros((n, n)); S = np.zeros((n, n))
    for a in range(len(xi)):
        x = xi[a]; ef = np.exp(-2.0 * alpha * x)
        for b in range(len(eta)):
            e = eta[b]; r_A = hR * (x + e)
            if r_A <= 0:
                continue
            jac = (x**2 - e**2)
            wbase = wxi[a] * weta[b] * jac * (hR**3) * ef * 2.0 * np.pi
            vh = V_H_closed(r_A, ZC_LI)
            for i, (ji, li) in enumerate(jl_pairs):
                gi = x**ji * e**li
                for jj, (jx, lx) in enumerate(jl_pairs):
                    gj = x**jx * e**lx
                    S[i, jj] += wbase * gi * gj
                    VH[i, jj] += wbase * vh * gi * gj
    return VH, S


def _twoelec_VH(basis_p, R, alpha):
    """2-electron core-Hartree matrix <g_i|V_H(r1A)+V_H(r2A)|g_j> at p=0 (factorizes:
    V_H acts one-body, no r12 coupling). Same ordering as build_basis_full."""
    # unique (j,l) for electron 1 and (k,m) for electron 2 -> a shared (a,b) index set
    idx = {}
    for (b, p) in basis_p:
        idx.setdefault((b.j, b.l), None); idx.setdefault((b.k, b.m), None)
    jl = list(idx.keys())
    pos = {p: i for i, p in enumerate(jl)}
    grid = build_quadrature_grids(N_xi=40, N_eta=30, N_phi=4, xi_max=20.0)
    VH1, S1 = _one_elec_VH_and_S(jl, R, alpha, grid)
    n = len(basis_p)
    W = np.zeros((n, n))
    for i, (bi, _) in enumerate(basis_p):
        a1, b1 = pos[(bi.j, bi.l)], pos[(bi.k, bi.m)]
        for jj, (bj, _) in enumerate(basis_p):
            a2, b2 = pos[(bj.j, bj.l)], pos[(bj.k, bj.m)]
            W[i, jj] = VH1[a1, a2] * S1[b1, b2] + S1[a1, a2] * VH1[b1, b2]
    return W


def E_tot(R):
    basis = build_basis_full(JMAX, LMAX, ALPHA, p_set=(0,))   # p=0 first look
    S, H = m.assemble_hetero(basis, R, ALPHA, Z_LI, Z_H, l_neumann=16, dps=30)
    VH2 = _twoelec_VH(basis, R, ALPHA)
    e_val = solve_canonical(S, H + VH2)[0]
    V_Hdens_at_H = V_H_closed(R, ZC_LI)          # potential of 1s^2 core at the H nucleus
    E_Rdep = e_val + Z_LI * Z_H / R - V_Hdens_at_H
    return E_Rdep, e_val


def main():
    t0 = time.time()
    print(f"Frozen-core prolate LiH FIRST look  (j,l)=({JMAX},{LMAX}) alpha={ALPHA} "
          f"(Hartree only, no projector/exchange)")
    E, Ev = [], []
    for R in R_GRID:
        er, ev = E_tot(R)
        E.append(er); Ev.append(ev)
        print(f"  R={R:.3f}  E_val={ev:.5f}  E_Rdep={er:.5f}  "
              f"E_tot={er+E_CORE:.5f}  [{time.time()-t0:.0f}s]", flush=True)
    E = np.array(E)
    # collapse check: E_val should be a sane valence energy (~ -0.8..-1.5 Ha), not -inf
    print(f"\n  min E_val = {min(Ev):.4f}  (collapse if << -2 Ha)")
    p = np.poly1d(np.polyfit(np.array(R_GRID) - R_E_EXP, E, min(4, len(R_GRID) - 1)))
    dp, ddp = p.deriv(1), p.deriv(2)
    roots = dp.r[np.isreal(dp.r)].real
    cand = [r + R_E_EXP for r in roots if ddp(r) > 0 and min(R_GRID) < r + R_E_EXP < max(R_GRID)]
    Req = float(min(cand, key=lambda rr: p(rr - R_E_EXP))) if cand else float('nan')
    err = (Req - R_E_EXP) / R_E_EXP * 100 if np.isfinite(Req) else float('nan')
    print("\n" + "=" * 60)
    print(f"R_eq(computed) = {Req:.4f} bohr   exp = {R_E_EXP}   drift = {err:+.1f}%")
    print(f"  (composed/balanced recipe drifts +8.8%; prolate HeH+ PoC was -0.5% at l=3)")
    print("=" * 60)


if __name__ == '__main__':
    main()
