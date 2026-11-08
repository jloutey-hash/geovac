"""Phase 0 (relaxable-analytic-core route): the CORE-BREATHING test for prolate LiH.

Increment 6 pinned the frozen core as the R_eq-drift culprit: WITH r12 the drift
STILL grows outward with angular basis (+3.4% l=2 -> +5.5% l=3), the OPPOSITE of
HeH+ (which converges), and the ONE variable differing is the FROZEN (rigid,
zc=2.6875-fixed) Li 1s^2 core.  The all-electron cure hit a NUMERICAL grid wall
(tight-core R-dependent resolution).  This driver tries the untried middle path:
keep the clean ANALYTIC/atom-centered core (no grid wall) but let it RELAX.

Phase 0 = the cheapest relaxation: let the core BREATHE (one variational knob,
the 1s exponent zc), staying spherical.  If the total-energy-optimal zc*(R) varies
with R and pulls R_eq inward from +5.5% toward the true 3.015 bohr, radial
relaxation is a real lever and we continue (Phase 1 = core polarization).  If zc*
is pinned at 2.6875 and R_eq is unmoved, breathing is inert and the lever is
polarization -- and we know that in ONE run.

zc enters the frozen-core model in exactly three consistent places:
  1. core-Hartree screening felt by the valence     V_H(r_A; zc)      (grid + Vdens)
  2. the 1s core in the Huzinaga projector            1s(r_A) = e^{-zc r_A}
  3. the isolated core energy                          E_core(zc)
The bare valence Hamiltonian (T + V_ne[-3/r_A,-1/r_B] + V_ee) is zc-INDEPENDENT, so
assemble_hetero (the mpf bottleneck) is built ONCE per (R,alpha) and reused across
all zc -- the breathing sweep is nearly free.

E_core(zc): a single-exponent 1s^2 He-like core at Z=3 (Slater):
  E_core_var(zc) = 2*(1/2 zc^2 - 3 zc) + (5/8) zc = zc^2 - 5.375 zc      (min at 2.6875)
anchored so E_core(2.6875) == the validated L.E_CORE (-7.2799), i.e. the frozen
baseline is reproduced EXACTLY at zc=zc0 and breathing enters as the zc-variation.
(The (5/8)zc self-repulsion is the SAME Hartree energy V_H represents -> consistent.)

E_tot(R) = e_val(R,alpha; zc, lambda) + Z_LI*Z_H/R - V_H_closed(R, zc) + E_core(zc)

Run from root:  python debug/lih_r12_breathing.py [jmax] [lmax] [lambda] [l_neumann]
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from geovac import prolate_recondition as pr                 # noqa: E402,F401
import prolate_r12_mpf as m                                  # noqa: E402
from heh_probe import build_basis_full                       # noqa: E402
from r12ci_first_energy import solve_canonical               # noqa: E402
from prolate_core_hartree import V_H_closed, ZC_LI           # noqa: E402
import lih_frozen_core_first as L                            # noqa: E402
import lih_r12_coupled as C                                  # noqa: E402

R_E_EXP = 3.015
JMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2
LMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2
LAM = float(sys.argv[3]) if len(sys.argv) > 3 else 1000.0
LNEU = int(sys.argv[4]) if len(sys.argv) > 4 else 20
ALPHAS = [1.0, 1.4]
R_GRID = [2.70, 2.85, 3.015, 3.20, 3.45]
# core-breathing grid: contraction (>2.6875) and expansion (<2.6875) of the 1s^2 core
ZC_GRID = np.array([2.30, 2.45, 2.55, 2.60, 2.6875, 2.75, 2.85, 3.00, 3.20])
ZC0 = ZC_LI  # 2.6875, the isolated (frozen) optimum


def E_core_var(zc):
    """Single-exponent 1s^2 He-like core energy at Z=3 (Slater), min at zc=2.6875."""
    return zc**2 - 5.375 * zc


# anchor so E_core(2.6875) == the validated frozen L.E_CORE (-7.2799)
E_CORE_ANCHOR = L.E_CORE - E_core_var(ZC0)


def E_core(zc):
    return E_core_var(zc) + E_CORE_ANCHOR


def assemble_mpf(R, alpha):
    """The zc-INDEPENDENT bottleneck: bare valence S, H (T+V_ne+V_ee).  Built once."""
    basis = build_basis_full(JMAX, LMAX, alpha, p_set=(0, 1))
    S, H = m.assemble_hetero(basis, R, alpha, L.Z_LI, L.Z_H, l_neumann=LNEU, dps=30)
    return basis, S, H


def e_val_at_zc(basis, S, H, R, alpha, zc):
    """Valence eigenvalue for a given core exponent zc (cheap float64 grid ops)."""
    G = C.grid_arrays(R, alpha, N_xi=32, N_eta=24, xi_max=18.0, zc=zc)
    VH2 = C.vh_coupled(basis, G, alpha)
    Plam = LAM * C.projector_coupled(basis, G, alpha)
    return solve_canonical(S, H + VH2 + Plam)[0]


def _parabolic_vertex(zc3, e3):
    """Vertex (zc*, E*) of the parabola through 3 points; falls back to the min."""
    (x0, x1, x2), (y0, y1, y2) = zc3, e3
    denom = (x0 - x1) * (x0 - x2) * (x1 - x2)
    if abs(denom) < 1e-18:
        i = int(np.argmin(e3)); return zc3[i], e3[i]
    a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
    b = (x2**2 * (y0 - y1) + x1**2 * (y2 - y0) + x0**2 * (y1 - y2)) / denom
    if a <= 0:
        i = int(np.argmin(e3)); return zc3[i], e3[i]
    xv = -b / (2 * a)
    lo, hi = min(zc3), max(zc3)
    if not (lo <= xv <= hi):
        i = int(np.argmin(e3)); return zc3[i], e3[i]
    c = y1 - a * x1**2 - b * x1
    return xv, a * xv**2 + b * xv + c


def E_tot_breathing(basis, S, H, R, alpha):
    """min over the zc grid (parabolic-refined) -> (E_tot*, zc*).  Also returns the
    fixed-zc0 value for the apples-to-apples baseline."""
    evs = np.array([e_val_at_zc(basis, S, H, R, alpha, zc) for zc in ZC_GRID])
    etot = evs + L.Z_LI * L.Z_H / R - np.array([V_H_closed(R, zc) for zc in ZC_GRID]) \
        + np.array([E_core(zc) for zc in ZC_GRID])
    i = int(np.argmin(etot))
    if 0 < i < len(ZC_GRID) - 1:
        zc_star, e_star = _parabolic_vertex(ZC_GRID[i-1:i+2], etot[i-1:i+2])
    else:
        zc_star, e_star = ZC_GRID[i], etot[i]
    # fixed-zc0 baseline on the SAME assembly
    ev0 = e_val_at_zc(basis, S, H, R, alpha, ZC0)
    e0 = ev0 + L.Z_LI * L.Z_H / R - V_H_closed(R, ZC0) + E_core(ZC0)
    return e_star, zc_star, e0


def _req_from_curve(R_grid, E):
    p = np.poly1d(np.polyfit(np.array(R_grid) - R_E_EXP, E, min(4, len(R_grid) - 1)))
    dp, ddp = p.deriv(1), p.deriv(2)
    roots = dp.r[np.isreal(dp.r)].real
    cand = [r + R_E_EXP for r in roots
            if ddp(r) > 0 and min(R_grid) < r + R_E_EXP < max(R_grid)]
    if not cand:
        return float('nan')
    return float(min(cand, key=lambda rr: p(rr - R_E_EXP)))


def main():
    t0 = time.time()
    print(f"CORE-BREATHING prolate LiH R_eq  (j,l)=({JMAX},{LMAX})  lambda={LAM:.0f} "
          f" l_neu={LNEU}   exp R_e={R_E_EXP}", flush=True)
    print(f"  frozen baseline (zc={ZC0}): +3.4% (l=2) / +5.5% (l=3).  "
          f"GO if breathing pulls R_eq inward.", flush=True)
    print(f"  zc grid: {ZC_GRID}", flush=True)
    E_free, E_fix, ZCS = [], [], []
    for R in R_GRID:
        best_free = None  # (E*, zc*, alpha)
        best_fix = None   # (E0, alpha)
        for a in ALPHAS:
            basis, S, H = assemble_mpf(R, a)
            e_star, zc_star, e0 = E_tot_breathing(basis, S, H, R, a)
            if best_free is None or e_star < best_free[0]:
                best_free = (e_star, zc_star, a)
            if best_fix is None or e0 < best_fix[0]:
                best_fix = (e0, a)
        E_free.append(best_free[0]); ZCS.append(best_free[1]); E_fix.append(best_fix[0])
        print(f"  R={R:.3f}  E_fix={best_fix[0]:.5f}  E_free={best_free[0]:.5f}  "
              f"zc*={best_free[1]:.4f}  dE(breath)={ (best_free[0]-best_fix[0])*1e3:+.2f} mHa "
              f" (a_free={best_free[2]})  [{time.time()-t0:.0f}s]", flush=True)
    E_free = np.array(E_free); E_fix = np.array(E_fix)
    Req_fix = _req_from_curve(R_GRID, E_fix)
    Req_free = _req_from_curve(R_GRID, E_free)
    err_fix = (Req_fix - R_E_EXP) / R_E_EXP * 100 if np.isfinite(Req_fix) else float('nan')
    err_free = (Req_free - R_E_EXP) / R_E_EXP * 100 if np.isfinite(Req_free) else float('nan')
    print("\n" + "=" * 68, flush=True)
    print(f"  zc*(R):  " + "  ".join(f"{r:.2f}:{z:.3f}" for r, z in zip(R_GRID, ZCS)),
          flush=True)
    print(f"  R_eq FROZEN (zc={ZC0}) = {Req_fix:.4f} bohr   drift = {err_fix:+.2f}%",
          flush=True)
    print(f"  R_eq BREATHING (zc free) = {Req_free:.4f} bohr   drift = {err_free:+.2f}%",
          flush=True)
    print(f"  --> breathing moved R_eq by {(Req_free-Req_fix):+.4f} bohr "
          f"({err_free-err_fix:+.2f} pp)", flush=True)
    if np.isfinite(err_free) and np.isfinite(err_fix):
        verdict = ("GO: breathing pulls inward toward truth"
                   if (err_free < err_fix - 0.3) else
                   "NO-GO: breathing inert/wrong-way -> lever is polarization (Phase 1)")
        print(f"  VERDICT: {verdict}", flush=True)
    print("=" * 68, flush=True)


if __name__ == '__main__':
    main()
