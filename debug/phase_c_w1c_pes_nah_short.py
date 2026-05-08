"""Phase C-W1c short NaH PES: 8 strategic R-points to characterize the curve.

Purpose: in 8 R-points, determine whether the screened NaH PES has an
interior minimum (BOUND), descends monotonically toward small R (still
overattracted, residual W1c), or descends toward large R (artifact / not
yet asymptotic).

Sampling: dense at small/medium R (where equilibrium would be), sparse at
large R (asymptote).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from openfermion.linalg import get_sparse_operator
from scipy.sparse.linalg import eigsh

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import nah_spec


def main():
    spec = nah_spec()
    R_grid = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0, 15.0]
    print(f"NaH short PES: {len(R_grid)} R-points")
    print(f"{'R':>6} {'E_bare':>13} {'E_scr':>13} {'vne_b':>9} {'vne_s':>9}", flush=True)
    rows = []
    for R in R_grid:
        rb = build_balanced_hamiltonian(spec, R=R, screened_cross_center=False)
        rs = build_balanced_hamiltonian(spec, R=R, screened_cross_center=True)
        Eb = float(np.min(eigsh(get_sparse_operator(rb['qubit_op']), k=2, which='SA', return_eigenvectors=False)))
        Es = float(np.min(eigsh(get_sparse_operator(rs['qubit_op']), k=2, which='SA', return_eigenvectors=False)))
        vneb = float(np.trace(rb['h1_cross_vne']))
        vnes = float(np.trace(rs['h1_cross_vne']))
        rows.append({'R_bohr': float(R), 'E_bare_Ha': Eb, 'E_screened_Ha': Es,
                     'vne_trace_bare': vneb, 'vne_trace_screened': vnes})
        print(f"{R:>6.2f} {Eb:>13.5f} {Es:>13.5f} {vneb:>9.3f} {vnes:>9.3f}", flush=True)

    R_arr = np.array([r['R_bohr'] for r in rows])
    E_s = np.array([r['E_screened_Ha'] for r in rows])
    E_b = np.array([r['E_bare_Ha'] for r in rows])
    idx_s = int(np.argmin(E_s))
    idx_b = int(np.argmin(E_b))
    print()
    print(f"BARE min:     R={R_arr[idx_b]:.2f}, E={E_b[idx_b]:.5f} Ha (idx {idx_b}/{len(R_arr)})")
    print(f"SCREENED min: R={R_arr[idx_s]:.2f}, E={E_s[idx_s]:.5f} Ha (idx {idx_s}/{len(R_arr)})")
    print(f"E_screened range: [{E_s.min():.4f}, {E_s.max():.4f}] = {E_s.max()-E_s.min():.4f} Ha")
    print(f"E_screened at largest R={R_arr[-1]:.1f}: {E_s[-1]:.5f}")
    is_bound = (idx_s > 0) and (idx_s < len(R_arr) - 1)
    print(f"SCREENED interior min: {is_bound}")
    print(f"D_e_proxy_screened = E(R_max) - E(R_min) = {E_s[-1] - E_s[idx_s]:+.5f} Ha")

    out = {
        'molecule': 'NaH',
        'R_grid_bohr': [float(r) for r in R_grid],
        'rows': rows,
        'idx_min_screened': idx_s,
        'R_min_screened_bohr': float(R_arr[idx_s]),
        'E_min_screened_Ha': float(E_s[idx_s]),
        'E_max_R_screened_Ha': float(E_s[-1]),
        'D_e_proxy_screened': float(E_s[-1] - E_s[idx_s]),
        'is_bound_screened_interior': bool(is_bound),
    }
    out_path = Path('debug/data/multifocal_c_w1c_pes_nah_short.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"Wrote {out_path}", flush=True)


if __name__ == '__main__':
    main()
