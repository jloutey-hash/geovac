"""Increment 6 PAYOFF: frozen-core prolate LiH R_eq scan WITH r12 coupling
(p_set={0,1}), the apples-to-apples test against the HeH+ PoC that converged.

The load-bearing question (PI's bet: the FROZEN CORE is the culprit, not r12):
  LiH WITHOUT r12 (p=0):  R_eq +2.2% (l=2) -> +7.1% (l=3)   [DIVERGES outward]
  HeH+  WITH   r12:       R_eq -4.5% (l=2) -> -0.5% (l=3)   [CONVERGES]
  LiH WITH r12 (this):    R_eq  ??? (l=2) ->  ??? (l=3)
    - shrinks toward 0  -> r12 was the lever (PI's bet WRONG)
    - grows toward +8.8% -> frozen core is the culprit (PI's bet RIGHT)

Same model as lih_req_scan (Hartree V_H + Huzinaga projector, NO exchange/r12)
EXCEPT V_H and the projector are now r12-COUPLED (lih_r12_coupled) so the basis
can carry p={0,1}.  Isolates the r12 effect: the ONLY change from the p=0 baseline
is turning r12 on.

E_tot(R) = e_val(R;lambda) + E_core + 3/R - V_Hdens(R),  lambda in the projector plateau.

Run from root:  python debug/lih_r12_req_scan.py [jmax] [lmax] [lambda] [l_neumann]
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from geovac import prolate_recondition as pr                 # noqa: E402
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


def E_tot(R, alpha):
    basis = build_basis_full(JMAX, LMAX, alpha, p_set=(0, 1))
    S, H = m.assemble_hetero(basis, R, alpha, L.Z_LI, L.Z_H, l_neumann=LNEU, dps=30)
    G = C.grid_arrays(R, alpha, N_xi=32, N_eta=24, xi_max=18.0)
    VH2 = C.vh_coupled(basis, G, alpha)
    Plam = LAM * C.projector_coupled(basis, G, alpha)
    e_val = solve_canonical(S, H + VH2 + Plam)[0]
    return e_val + L.Z_LI * L.Z_H / R - V_H_closed(R, ZC_LI) + L.E_CORE, e_val, len(basis)


def main():
    t0 = time.time()
    print(f"Frozen-core prolate LiH R_eq scan WITH r12  (j,l)=({JMAX},{LMAX})  "
          f"lambda={LAM:.0f}  l_neu={LNEU}   exp R_e={R_E_EXP}", flush=True)
    print(f"  p=0 baseline: +2.2% (l=2) / +7.1% (l=3);  HeH+ r12: -4.5%->-0.5%", flush=True)
    E = []
    for R in R_GRID:
        best = None; na = None
        for a in ALPHAS:
            er, ev, na = E_tot(R, a)
            if best is None or er < best[0]:
                best = (er, ev, a)
        E.append(best[0])
        print(f"  R={R:.3f}  E_tot={best[0]:.5f}  E_val={best[1]:.5f}  "
              f"(a*={best[2]}, n={na})  [{time.time()-t0:.0f}s]", flush=True)
    E = np.array(E)
    p = np.poly1d(np.polyfit(np.array(R_GRID) - R_E_EXP, E, min(4, len(R_GRID) - 1)))
    dp, ddp = p.deriv(1), p.deriv(2)
    roots = dp.r[np.isreal(dp.r)].real
    cand = [r + R_E_EXP for r in roots if ddp(r) > 0 and min(R_GRID) < r + R_E_EXP < max(R_GRID)]
    Req = float(min(cand, key=lambda rr: p(rr - R_E_EXP))) if cand else float('nan')
    err = (Req - R_E_EXP) / R_E_EXP * 100 if np.isfinite(Req) else float('nan')
    print("\n" + "=" * 64, flush=True)
    print(f"R_eq(WITH r12) = {Req:.4f} bohr   exp = {R_E_EXP}   drift = {err:+.1f}%", flush=True)
    print(f"E_min(grid)    = {E.min():.5f} Ha   LiH ref ~ -8.070", flush=True)
    print(f"  variational? E_min > -8.070 : {'YES' if E.min() > -8.070 else 'NO (overshoot)'}", flush=True)
    print("=" * 64, flush=True)


if __name__ == '__main__':
    main()
