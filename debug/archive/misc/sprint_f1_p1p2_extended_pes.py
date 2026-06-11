"""Sprint F1 Phase 2 follow-on: extended PES scan to verify internal minimum.

The main F1 P2 driver runs only 3 R-points for the mini-PES. To check
whether the W1c+mz combined architecture truly has an internal minimum
(or just looks like binding at R_eq vs R_diss without one), run all four
architectures on a denser R-grid.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(line_buffering=True)

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import nah_spec


def main():
    out_path = Path('debug/data/sprint_f1_p1p2_extended_pes.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    R_grid = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0]
    archs = [
        ('bare', False, False),
        ('bare+mz', False, True),
        ('W1c', True, False),
        ('W1c+mz', True, True),
    ]

    results = {label: {'R': [], 'E': []} for label, _, _ in archs}

    spec = nah_spec(max_n=2)
    n_e = sum(b.n_electrons for b in spec.blocks)

    for label, screened, mz in archs:
        print(f"--- {label} ---", flush=True)
        for R in R_grid:
            t0 = time.perf_counter()
            ham = build_balanced_hamiltonian(
                spec, R=R, screened_cross_center=screened,
                multi_zeta_basis=mz, verbose=False,
                n_grid_vne=4000, L_max=4,
            )
            fci = coupled_fci_energy(ham, n_electrons=n_e)
            E = fci['E_coupled']
            results[label]['R'].append(R)
            results[label]['E'].append(E)
            print(f"  R={R:6.2f}  E={E:12.6f}  [{time.perf_counter()-t0:.1f}s]",
                  flush=True)
        print()

    # Summary table
    print(f"{'R':>6}  {'bare':>12}  {'bare+mz':>12}  {'W1c':>12}  "
          f"{'W1c+mz':>12}  {'mz-vs-W1c':>14}")
    for i, R in enumerate(R_grid):
        E_b = results['bare']['E'][i]
        E_bm = results['bare+mz']['E'][i]
        E_w = results['W1c']['E'][i]
        E_wm = results['W1c+mz']['E'][i]
        print(f"{R:6.2f}  {E_b:12.6f}  {E_bm:12.6f}  {E_w:12.6f}  "
              f"{E_wm:12.6f}  {E_wm-E_w:+14.4e}")
    print()

    print("--- Architecture descent depth diagnostic ---")
    print(f"{'arch':>10}  {'R_min':>6}  {'E_min':>12}  {'descent depth':>14}  "
          f"{'internal min':>12}")
    for label in results:
        Es = np.array(results[label]['E'])
        Rs = np.array(results[label]['R'])
        i_min = int(np.argmin(Es))
        R_min = float(Rs[i_min])
        E_min = float(Es[i_min])
        descent = float(Es[-1] - E_min)
        internal = bool(0 < i_min < len(Es) - 1)
        results[label]['R_min'] = R_min
        results[label]['E_min'] = E_min
        results[label]['descent_depth'] = descent
        results[label]['internal_min'] = internal
        print(f"  {label:>10}  {R_min:6.2f}  {E_min:12.6f}  "
              f"{descent:+14.4f}  {str(internal):>12}")
    print()

    # Save
    out = {
        'metadata': {
            'sprint': 'F1 Phase 2 extended PES',
            'date': '2026-05-23',
            'R_grid': R_grid,
            'comment': 'Extended PES scan to verify internal minimum (or lack thereof) for all four architectures.',
        },
        'results': results,
    }
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"[saved] {out_path}")


if __name__ == '__main__':
    main()
