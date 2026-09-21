"""VERDICT scan: all-electron (variational-core) prolate LiH R_eq.

The frozen-core prolate LiH drifts +3.4%->+5.5% (l=2->3) because the core is
frozen (v5.15.6). This scans the ALL-ELECTRON multi-exponent FCI (all 4 electrons
active, R-adaptive core, no PK, no frozen core) at two basis sizes:
  - does R_eq converge toward experiment (3.015 bohr) -> variational core CURES the drift
  - or still drift outward -> deeper than the frozen core

Absolute E is basis-limited (sigma-only, ~0.3 Ha high) but the core energy is
~R-independent, so R_eq (valence + V_NN balance) is the meaningful read-out;
stability of R_eq across basis size is the trust check.

Run from root:  python debug/lih_allelectron_req_scan.py
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import prolate_allelectron_fci as A                          # noqa: E402

R_E_EXP = 3.015
R_GRID = [2.60, 2.85, 3.015, 3.20, 3.45]
SCALE_SETS = {
    'Nsc3': [(2.9, 1.0), (1.5, 0.9), (0.8, 0.7)],
    'Nsc4': [(2.9, 1.0), (1.8, 1.0), (1.0, 0.8), (0.65, 0.65)],
}


def E_tot(R, scales):
    h1, eri, M, cond = A.build_mo_integrals_multiexp(
        R, scales, 3.0, 1.0, n_ang_max=1, N_xi_solve=6000, N_grid=44, xi_max=15.0)
    E, nd = A.fci_energy(h1, eri, M, nelec=4)
    return E + 3.0 / R


def req_from_grid(E):
    E = np.array(E)
    p = np.poly1d(np.polyfit(np.array(R_GRID) - R_E_EXP, E, min(4, len(R_GRID) - 1)))
    dp, ddp = p.deriv(1), p.deriv(2)
    roots = dp.r[np.isreal(dp.r)].real
    cand = [r + R_E_EXP for r in roots
            if ddp(r) > 0 and min(R_GRID) < r + R_E_EXP < max(R_GRID)]
    if not cand:
        return float('nan'), float('nan')
    Req = float(min(cand, key=lambda rr: p(rr - R_E_EXP)))
    return Req, (Req - R_E_EXP) / R_E_EXP * 100


def main():
    t0 = time.time()
    print(f"All-electron prolate LiH R_eq scan  exp R_e={R_E_EXP}", flush=True)
    print(f"  frozen-core (v5.15.6): +3.4% (l=2) -> +5.5% (l=3), drifts outward", flush=True)
    for tag, scales in SCALE_SETS.items():
        E = []
        for R in R_GRID:
            e = E_tot(R, scales)
            E.append(e)
            print(f"  [{tag}] R={R:.3f}  E_tot={e:.5f}  [{time.time()-t0:.0f}s]", flush=True)
        Req, err = req_from_grid(E)
        print(f"  ==> [{tag}] R_eq = {Req:.4f} bohr   drift = {err:+.1f}%   "
              f"E_min={min(E):.5f}", flush=True)
        print("", flush=True)


if __name__ == '__main__':
    main()
