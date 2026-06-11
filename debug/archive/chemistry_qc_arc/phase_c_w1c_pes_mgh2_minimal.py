"""Minimal Phase C-W1c MgH2 PES: 4 strategic R-points only.

Goal: characterize whether the screened PES has an interior minimum or
remains monotonic, with minimal compute time. Q=40 (4-electron), so each
R-point is several minutes via eigsh.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from openfermion.linalg import get_sparse_operator
from scipy.sparse.linalg import eigsh

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import mgh2_spec


def main():
    spec = mgh2_spec()
    # Symmetric R = R(Mg-H1) = R(Mg-H2). Equilibrium ~3.3 bohr (1.7 Å).
    R_grid = [2.5, 3.5, 5.0, 8.0]
    print(f"MgH2 minimal PES: {len(R_grid)} R-points", flush=True)
    print(f"{'R':>6} {'E_bare':>13} {'E_scr':>13} {'vne_b':>10} {'vne_s':>10} {'wall':>7}", flush=True)
    rows = []
    for R in R_grid:
        t0 = time.perf_counter()
        rb = build_balanced_hamiltonian(spec, R=R, screened_cross_center=False, verbose=False)
        rs = build_balanced_hamiltonian(spec, R=R, screened_cross_center=True, verbose=False)
        try:
            Eb = float(np.min(eigsh(get_sparse_operator(rb['qubit_op']),
                                    k=2, which='SA', return_eigenvectors=False, maxiter=5000)))
            Es = float(np.min(eigsh(get_sparse_operator(rs['qubit_op']),
                                    k=2, which='SA', return_eigenvectors=False, maxiter=5000)))
        except Exception as ex:
            print(f"  Eigsh failure at R={R}: {ex}", flush=True)
            Eb, Es = float('nan'), float('nan')
        vneb = float(np.trace(rb['h1_cross_vne']))
        vnes = float(np.trace(rs['h1_cross_vne']))
        elapsed = time.perf_counter() - t0
        rows.append({
            'R_bohr': float(R),
            'E_bare_Ha': Eb,
            'E_screened_Ha': Es,
            'vne_trace_bare': vneb,
            'vne_trace_screened': vnes,
            'one_norm_bare': float(rb['one_norm']),
            'one_norm_screened': float(rs['one_norm']),
            'N_pauli_bare': rb['N_pauli'],
            'N_pauli_screened': rs['N_pauli'],
            'wall_time_s': elapsed,
        })
        print(f"{R:>6.2f} {Eb:>13.5f} {Es:>13.5f} {vneb:>10.3f} {vnes:>10.3f} {elapsed:>7.1f}", flush=True)

    R_arr = np.array([r['R_bohr'] for r in rows])
    E_s = np.array([r['E_screened_Ha'] for r in rows])
    E_b = np.array([r['E_bare_Ha'] for r in rows])
    finite_s = np.isfinite(E_s)
    finite_b = np.isfinite(E_b)
    if not np.any(finite_s):
        verdict = 'all_nan_screened'
        idx_s = -1
    else:
        idx_s = int(np.argmin(np.where(finite_s, E_s, np.inf)))
        if idx_s == 0:
            verdict = 'monotonic_descending_to_small_R'
        elif idx_s == len(R_arr) - 1:
            verdict = 'monotonic_descending_to_large_R_artifact'
        else:
            verdict = 'interior_minimum'

    if np.any(finite_b):
        idx_b = int(np.argmin(np.where(finite_b, E_b, np.inf)))
        print(f"\nBARE   min: R={R_arr[idx_b]:.2f}, E={E_b[idx_b]:.5f} Ha (idx {idx_b})", flush=True)
    if idx_s >= 0:
        print(f"SCREEN min: R={R_arr[idx_s]:.2f}, E={E_s[idx_s]:.5f} Ha (idx {idx_s})", flush=True)
    print(f"Verdict: {verdict}", flush=True)
    if np.any(finite_s):
        print(f"E_scr range across R: [{E_s[finite_s].min():.4f}, {E_s[finite_s].max():.4f}] = "
              f"{E_s[finite_s].max() - E_s[finite_s].min():.4f} Ha", flush=True)

    # vne ratio across R
    vne_b = np.array([r['vne_trace_bare'] for r in rows])
    vne_s = np.array([r['vne_trace_screened'] for r in rows])
    safe = vne_s != 0
    ratios = np.where(safe, vne_b / np.where(safe, vne_s, 1), np.nan)
    print(f"vne ratios (bare/screened): {ratios}", flush=True)

    out = {
        'molecule': 'MgH2',
        'R_grid_bohr': R_grid,
        'rows': rows,
        'verdict_screened': verdict,
        'idx_min_screened': idx_s,
        'R_min_screened_bohr': float(R_arr[idx_s]) if idx_s >= 0 else None,
        'E_min_screened_Ha': float(E_s[idx_s]) if idx_s >= 0 else None,
        'E_at_R_max_screened_Ha': float(E_s[-1]) if np.isfinite(E_s[-1]) else None,
        'D_e_proxy_screened': float(E_s[-1] - E_s[idx_s]) if (idx_s >= 0 and np.isfinite(E_s[-1])) else None,
        'vne_ratio_mean_bare_over_scr': float(np.nanmean(ratios)),
    }
    out_path = Path('debug/data/multifocal_c_w1c_pes_mgh2_minimal.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"Wrote {out_path}", flush=True)


if __name__ == '__main__':
    main()
