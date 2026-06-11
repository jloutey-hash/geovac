"""Extended NaH PES scan for Phase C-W1c, to find equilibrium and asymptote.

The default (R=2..6) scan still showed the screened PES descending toward
small R, so we extend the scan in both directions:
  - down to R=1.5 to look for overshoot / turn
  - up to R=12 to find the asymptote
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from openfermion.linalg import get_sparse_operator
from scipy.sparse.linalg import eigsh

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import nah_spec


def main():
    spec = nah_spec()
    R_grid = (
        list(np.arange(1.5, 4.5, 0.25))
        + list(np.arange(4.5, 12.5, 0.5))
    )
    rows = []
    print(f"NaH extended PES: {len(R_grid)} R-points")
    print(f"{'R':>6} {'E_bare':>13} {'E_scr':>13} "
          f"{'vne_b':>9} {'vne_s':>9}")
    for R in R_grid:
        t0 = time.perf_counter()
        r_bare = build_balanced_hamiltonian(spec, R=R, screened_cross_center=False)
        r_scr = build_balanced_hamiltonian(spec, R=R, screened_cross_center=True)
        H_b = get_sparse_operator(r_bare['qubit_op'])
        H_s = get_sparse_operator(r_scr['qubit_op'])
        e_b = float(np.min(eigsh(H_b, k=2, which='SA', return_eigenvectors=False)))
        e_s = float(np.min(eigsh(H_s, k=2, which='SA', return_eigenvectors=False)))
        elapsed = time.perf_counter() - t0
        vneb = float(np.trace(r_bare['h1_cross_vne']))
        vnes = float(np.trace(r_scr['h1_cross_vne']))
        rows.append({
            'R_bohr': float(R),
            'E_bare_Ha': e_b, 'E_screened_Ha': e_s,
            'vne_trace_bare': vneb, 'vne_trace_screened': vnes,
            'wall_time_s': elapsed,
        })
        print(f"{R:>6.2f} {e_b:>13.5f} {e_s:>13.5f} {vneb:>9.3f} {vnes:>9.3f}")

    # Analyze
    R_arr = np.array([r['R_bohr'] for r in rows])
    E_s = np.array([r['E_screened_Ha'] for r in rows])
    E_b = np.array([r['E_bare_Ha'] for r in rows])
    idx_s = int(np.argmin(E_s))
    idx_b = int(np.argmin(E_b))
    print()
    print(f"BARE PES min: R={R_arr[idx_b]:.2f} bohr, E={E_b[idx_b]:.5f} Ha (idx {idx_b}/{len(R_arr)})")
    print(f"SCREENED min: R={R_arr[idx_s]:.2f} bohr, E={E_s[idx_s]:.5f} Ha (idx {idx_s}/{len(R_arr)})")
    print(f"Asymptotic E_screened at R={R_arr[-1]:.2f}: {E_s[-1]:.5f} Ha")
    print(f"Asymptotic E_bare at R={R_arr[-1]:.2f}: {E_b[-1]:.5f} Ha")
    is_bound_s = (idx_s > 0) and (idx_s < len(R_arr) - 1)
    print(f"SCREENED has interior minimum: {is_bound_s}")

    out = {
        'molecule': 'NaH',
        'R_grid_bohr': [float(r['R_bohr']) for r in rows],
        'rows': rows,
        'R_eq_screened_bohr': float(R_arr[idx_s]),
        'E_eq_screened_Ha': float(E_s[idx_s]),
        'E_inf_screened_Ha': float(E_s[-1]),
        'D_e_screened_Ha': float(E_s[-1] - E_s[idx_s]),
        'is_bound_screened': bool(is_bound_s),
    }
    out_path = Path('debug/data/multifocal_c_w1c_pes_nah_extended.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"Wrote {out_path}")


if __name__ == '__main__':
    main()
