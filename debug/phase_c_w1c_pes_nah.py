"""Phase C-W1c NaH PES regression: bare vs screened cross-center V_ne.

Recomputes the Sprint 7 NaH balanced-FCI PES with the new screened
cross-center V_ne path. Predicted outcome (from Phase B-W1c-diag):
the screened path should reduce the H-side overattraction by ~10x
and produce equilibrium binding (D_e > 0) instead of monotonic
overattraction.

Output: ``debug/data/multifocal_c_w1c_pes_nah.json``.
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


def run_pes(R_grid):
    spec = nah_spec()
    rows = []
    print(f"NaH PES: {len(R_grid)} R-points")
    print(f"{'R':>6} {'N_pauli':>9} {'E_bare':>12} {'E_scr':>12} "
          f"{'dE_b-s':>10} {'vne_bare':>10} {'vne_scr':>10}")
    for R in R_grid:
        t0 = time.perf_counter()
        r_bare = build_balanced_hamiltonian(
            spec, R=R, screened_cross_center=False, verbose=False,
        )
        r_scr = build_balanced_hamiltonian(
            spec, R=R, screened_cross_center=True, verbose=False,
        )
        # 2-electron NaH: ground state via lowest eigenvalue
        H_b = get_sparse_operator(r_bare['qubit_op'])
        H_s = get_sparse_operator(r_scr['qubit_op'])
        e_b = eigsh(H_b, k=3, which='SA', return_eigenvectors=False)
        e_s = eigsh(H_s, k=3, which='SA', return_eigenvectors=False)
        E_b = float(np.min(e_b))
        E_s = float(np.min(e_s))
        elapsed = time.perf_counter() - t0
        vneb = float(np.trace(r_bare['h1_cross_vne']))
        vnes = float(np.trace(r_scr['h1_cross_vne']))
        rows.append({
            'R_bohr': float(R),
            'E_bare_Ha': E_b,
            'E_screened_Ha': E_s,
            'shift_bare_minus_scr': E_b - E_s,
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


def analyze_pes(rows):
    """Find equilibrium R, D_e from the PES rows."""
    R_arr = np.array([r['R_bohr'] for r in rows])
    E_b = np.array([r['E_bare_Ha'] for r in rows])
    E_s = np.array([r['E_screened_Ha'] for r in rows])

    out = {}
    for tag, E in [('bare', E_b), ('screened', E_s)]:
        idx = int(np.argmin(E))
        E_min = float(E[idx])
        R_min = float(R_arr[idx])
        # Estimate dissociation energy: D_e = E_atomic_sum - E_min
        # NaH dissociates to Na + H. With balanced framework, asymptotic
        # energy at R -> infinity should ~ E(R = R_max).
        E_inf = float(E[-1])  # last R-point as proxy for R -> infinity
        D_e = E_inf - E_min
        # Bound check: equilibrium if min is interior and D_e > 0
        is_bound = (idx > 0) and (idx < len(R_arr) - 1) and (D_e > 0)
        out[tag] = {
            'R_min': R_min,
            'E_min': E_min,
            'E_inf_proxy': E_inf,
            'D_e_proxy': D_e,
            'is_bound': bool(is_bound),
            'idx_min': idx,
        }
    return out


def main():
    R_grid = np.arange(2.0, 6.5, 0.5).tolist()
    spec, rows = run_pes(R_grid)
    summary = analyze_pes(rows)

    print()
    print("=" * 72)
    print("Summary:")
    for tag, d in summary.items():
        bound = "BOUND" if d['is_bound'] else "NOT BOUND"
        print(f"  {tag}:  R_min={d['R_min']:.2f} bohr,  E_min={d['E_min']:.4f} Ha,  "
              f"D_e_proxy={d['D_e_proxy']:+.4f} Ha,  {bound}")

    out_dir = Path('debug/data')
    out_dir.mkdir(parents=True, exist_ok=True)
    out = {
        'molecule': 'NaH',
        'spec_name': spec.name,
        'description': (
            'Phase C-W1c: balanced FCI PES of NaH with bare vs screened '
            'cross-center V_ne. Bare path uses Z_Na = 11 regardless of [Ne] '
            'frozen core; screened path uses FrozenCore Z_eff^Na(rho) tail.'
        ),
        'R_grid_bohr': R_grid,
        'rows': rows,
        'summary': summary,
    }
    out_path = out_dir / 'multifocal_c_w1c_pes_nah.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote: {out_path}")


if __name__ == '__main__':
    main()
