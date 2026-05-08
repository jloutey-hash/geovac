"""Phase C-W1c HCl PES regression: bare vs screened cross-center V_ne.

HCl is 8-electron with [Ne] frozen core on Cl side and 3 lone pair blocks
+ 1 bond block. Q=50, so the FCI matrix is 2^50 sparse — only the lowest
eigenvalues are extracted via eigsh. Faster to scan a small number of R
to test whether the screened path produces equilibrium.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from openfermion.linalg import get_sparse_operator
from scipy.sparse.linalg import eigsh

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import hcl_spec


def run_pes(R_grid):
    spec = hcl_spec()
    rows = []
    print(f"HCl PES: {len(R_grid)} R-points (each can take a few minutes)")
    print(f"{'R':>6} {'N_pauli':>9} {'E_bare':>12} {'E_scr':>12} "
          f"{'dE':>10} {'vne_bare':>10} {'vne_scr':>10}")
    for R in R_grid:
        t0 = time.perf_counter()
        r_bare = build_balanced_hamiltonian(
            spec, R=R, screened_cross_center=False, verbose=False,
        )
        r_scr = build_balanced_hamiltonian(
            spec, R=R, screened_cross_center=True, verbose=False,
        )
        try:
            H_b = get_sparse_operator(r_bare['qubit_op'])
            H_s = get_sparse_operator(r_scr['qubit_op'])
            e_b = eigsh(H_b, k=2, which='SA', return_eigenvectors=False, maxiter=5000)
            e_s = eigsh(H_s, k=2, which='SA', return_eigenvectors=False, maxiter=5000)
            E_b = float(np.min(e_b))
            E_s = float(np.min(e_s))
        except Exception as ex:
            print(f"  Eigenvalue solve failed at R={R}: {ex}")
            E_b, E_s = float('nan'), float('nan')
        elapsed = time.perf_counter() - t0
        vneb = float(np.trace(r_bare['h1_cross_vne']))
        vnes = float(np.trace(r_scr['h1_cross_vne']))
        rows.append({
            'R_bohr': float(R),
            'E_bare_Ha': E_b,
            'E_screened_Ha': E_s,
            'shift_bare_minus_scr': E_b - E_s if not (np.isnan(E_b) or np.isnan(E_s)) else float('nan'),
            'vne_trace_bare': vneb,
            'vne_trace_screened': vnes,
            'N_pauli_bare': r_bare['N_pauli'],
            'N_pauli_screened': r_scr['N_pauli'],
            'one_norm_bare': float(r_bare['one_norm']),
            'one_norm_screened': float(r_scr['one_norm']),
            'wall_time_s': elapsed,
        })
        print(f"{R:>6.2f} {r_bare['N_pauli']:>9} {E_b:>12.5f} {E_s:>12.5f} "
              f"{E_b - E_s:>+10.5f} {vneb:>10.3f} {vnes:>10.3f}")

    return spec, rows


def main():
    R_grid = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0]
    spec, rows = run_pes(R_grid)

    R_arr = np.array([r['R_bohr'] for r in rows])
    E_b = np.array([r['E_bare_Ha'] for r in rows])
    E_s = np.array([r['E_screened_Ha'] for r in rows])
    summary = {}
    for tag, E in [('bare', E_b), ('screened', E_s)]:
        if np.all(np.isnan(E)):
            summary[tag] = {'note': 'all NaN'}
            continue
        idx = int(np.nanargmin(E))
        E_min = float(E[idx])
        R_min = float(R_arr[idx])
        E_inf = float(E[-1]) if not np.isnan(E[-1]) else float('nan')
        summary[tag] = {
            'R_min': R_min, 'E_min': E_min, 'E_inf_proxy': E_inf,
            'D_e_proxy': E_inf - E_min, 'idx_min': idx,
            'is_bound': bool((idx > 0) and (idx < len(R_arr) - 1)
                             and (E_inf - E_min > 0)),
        }

    print()
    print("=" * 72)
    for tag, d in summary.items():
        if 'note' in d:
            print(f"  {tag}: {d['note']}")
            continue
        bound = "BOUND" if d.get('is_bound') else "NOT BOUND"
        print(f"  {tag}:  R_min={d['R_min']:.2f} bohr,  E_min={d['E_min']:.4f} Ha,  "
              f"D_e_proxy={d['D_e_proxy']:+.4f} Ha,  {bound}")

    out_dir = Path('debug/data')
    out_dir.mkdir(parents=True, exist_ok=True)
    out = {
        'molecule': 'HCl',
        'spec_name': spec.name,
        'description': 'Phase C-W1c balanced FCI PES, bare vs screened V_ne.',
        'R_grid_bohr': R_grid,
        'rows': rows,
        'summary': summary,
    }
    out_path = out_dir / 'multifocal_c_w1c_pes_hcl.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote: {out_path}")


if __name__ == '__main__':
    main()
