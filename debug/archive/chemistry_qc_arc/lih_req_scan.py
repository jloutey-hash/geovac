"""Frozen-core prolate LiH R_eq scan WITH the Huzinaga orthogonality projector
(increments 2-3 complete at Hartree level; no core-valence exchange, no r12 yet).

The payoff question: does the variational-prolate recipe put LiH's bond near the
experimental R_e = 3.015 bohr, instead of the composed/balanced recipe's +8.8%
(3.28 bohr) OUTWARD drift that worsens with basis?

E_tot(R) = E_val(R; lambda) + E_core + 3/R - V_Hdens(R), lambda in the projector
plateau (>=100, from lih_projector_test.py). Angular basis is the geometry lever
(HeH+ PoC: -4.5% at l=2 -> -0.5% at l=3), so we run l=2 and l=3.

Run from root:  python debug/lih_req_scan.py [jmax] [lmax] [lambda]
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
from lih_projector_test import build                         # noqa: E402

R_E_EXP = 3.015
JMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2
LMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2
LAM = float(sys.argv[3]) if len(sys.argv) > 3 else 1000.0
ALPHAS = [1.0, 1.3, 1.6, 2.0]   # widened (a* pegged at 1.3 in the first run)
R_GRID = [2.60, 2.85, 3.015, 3.20, 3.45, 3.75]


def E_tot(R, alpha):
    basis = build_basis_full(JMAX, LMAX, alpha, p_set=(0,))
    S, H = m.assemble_hetero(basis, R, alpha, L.Z_LI, L.Z_H, l_neumann=16, dps=30)
    VH2, Plam = build(basis, R, alpha, LAM)
    e_val = solve_canonical(S, H + VH2 + Plam)[0]
    return e_val + L.Z_LI * L.Z_H / R - V_H_closed(R, ZC_LI) + L.E_CORE


def main():
    t0 = time.time()
    print(f"Frozen-core prolate LiH R_eq scan  (j,l)=({JMAX},{LMAX})  lambda={LAM:.0f}  "
          f"(Hartree + projector; no exchange/r12)   exp R_e={R_E_EXP}")
    E = []
    for R in R_GRID:
        vals = [E_tot(R, a) for a in ALPHAS]
        E.append(min(vals))
        print(f"  R={R:.3f}  E_tot={min(vals):.5f}  (a*={ALPHAS[int(np.argmin(vals))]})  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    E = np.array(E)
    p = np.poly1d(np.polyfit(np.array(R_GRID) - R_E_EXP, E, min(4, len(R_GRID) - 1)))
    dp, ddp = p.deriv(1), p.deriv(2)
    roots = dp.r[np.isreal(dp.r)].real
    cand = [r + R_E_EXP for r in roots if ddp(r) > 0 and min(R_GRID) < r + R_E_EXP < max(R_GRID)]
    Req = float(min(cand, key=lambda rr: p(rr - R_E_EXP))) if cand else float('nan')
    err = (Req - R_E_EXP) / R_E_EXP * 100 if np.isfinite(Req) else float('nan')
    print("\n" + "=" * 62)
    print(f"R_eq(computed) = {Req:.4f} bohr   exp = {R_E_EXP} bohr   drift = {err:+.1f}%")
    print(f"E_min(grid)    = {E.min():.5f} Ha   LiH ref ~ -8.070")
    print(f"  composed/balanced = +8.8% (3.28) worsening;  prolate HeH+ PoC = -0.5% @ l=3")
    print("=" * 62)


if __name__ == '__main__':
    main()
