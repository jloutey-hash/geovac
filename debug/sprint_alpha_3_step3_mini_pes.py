"""Sprint alpha-PES Step 3: mini-PES at R in {3.0, 3.5, 4.0, 5.0} bohr.

Since Step 2 revealed that multi-zeta has BIT-ZERO effect on the FCI
energy (because the FCI GS lives on the H side, not the Na 3s),
Step 3 documents this for the record across R values and computes
the comparison to the W1c-only Track 3 baseline.

This confirms the diagnostic-level finding: multi-zeta is the wrong
W1c-residual axis on the FCI variational state at NaH max_n=2.
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


def run_one(R: float, screened: bool, multi_zeta: bool) -> dict:
    spec = nah_spec(R=R, max_n=2)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=4000, L_max=4,
        screened_cross_center=screened,
        multi_zeta_basis=multi_zeta,
        verbose=False,
    )
    fci = coupled_fci_energy(ham, n_electrons=n_e, verbose=False)
    return {
        'R': float(R),
        'screened': bool(screened),
        'multi_zeta': bool(multi_zeta),
        'E_total': float(fci['E_coupled']),
        'V_NN_plus_E_core': float(ham['nuclear_repulsion']),
        'h1_trace_cross_vne': float(np.trace(ham['h1_cross_vne'])),
    }


def has_internal_min(pts):
    if len(pts) < 3:
        return None
    E_arr = [p['E_total'] for p in pts]
    i_min = int(np.argmin(E_arr))
    is_internal = 0 < i_min < len(E_arr) - 1
    return {
        'i_min': i_min,
        'R_min': float(pts[i_min]['R']),
        'E_min': float(E_arr[i_min]),
        'has_internal_min': is_internal,
    }


def main():
    out_path = Path('debug/data/sprint_alpha_3_step3_mini_pes.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    R_grid = [2.5, 3.0, 3.5, 4.0, 5.0]
    R_eq_exp = 3.566

    print("=== Step 3 — Mini-PES with three architectures ===\n", flush=True)
    print(f"R_eq_exp = {R_eq_exp:.3f} bohr,  R_grid = {R_grid}\n", flush=True)

    # --- A: bare baseline ---
    print(f"--- A: bare (Track 3 baseline, no W1c, no multi-zeta) ---", flush=True)
    pts_a = []
    for R in R_grid:
        t0 = time.perf_counter()
        pt = run_one(R, screened=False, multi_zeta=False)
        print(f"  R={R:.3f}: E={pt['E_total']:.6f} Ha  "
              f"tr(h1_cross_vne)={pt['h1_trace_cross_vne']:.4f}  "
              f"[{time.perf_counter()-t0:.1f}s]", flush=True)
        pts_a.append(pt)

    # --- B: bare + multi-zeta (M-Y Path A naive) ---
    print(f"\n--- B: bare + multi_zeta (M-Y Path A naive) ---", flush=True)
    pts_b = []
    for R in R_grid:
        t0 = time.perf_counter()
        pt = run_one(R, screened=False, multi_zeta=True)
        print(f"  R={R:.3f}: E={pt['E_total']:.6f} Ha  "
              f"tr(h1_cross_vne)={pt['h1_trace_cross_vne']:.4f}  "
              f"[{time.perf_counter()-t0:.1f}s]", flush=True)
        pts_b.append(pt)

    # --- C: W1c screened only (Track 3) ---
    print(f"\n--- C: W1c screened (Track 3, no multi-zeta) ---", flush=True)
    pts_c = []
    for R in R_grid:
        t0 = time.perf_counter()
        pt = run_one(R, screened=True, multi_zeta=False)
        print(f"  R={R:.3f}: E={pt['E_total']:.6f} Ha  "
              f"tr(h1_cross_vne)={pt['h1_trace_cross_vne']:.4f}  "
              f"[{time.perf_counter()-t0:.1f}s]", flush=True)
        pts_c.append(pt)

    # --- Compare ---
    print(f"\n=== Comparison ===", flush=True)
    print(f"{'R':>6}  {'E_bare':>12}  {'E_bare+mz':>12}  {'E_W1c':>12}  "
          f"{'diff_mz_vs_bare':>16}", flush=True)
    for i, R in enumerate(R_grid):
        print(f"  {R:5.3f}  {pts_a[i]['E_total']:12.6f}  {pts_b[i]['E_total']:12.6f}  "
              f"{pts_c[i]['E_total']:12.6f}  "
              f"{pts_b[i]['E_total']-pts_a[i]['E_total']:+16.3e}", flush=True)

    # Diagnostics
    diag_a = has_internal_min(pts_a)
    diag_b = has_internal_min(pts_b)
    diag_c = has_internal_min(pts_c)
    print(f"\n--- Internal-minimum diagnostic ---", flush=True)
    if diag_a:
        print(f"  A (bare):       R_min={diag_a['R_min']:.3f}, E_min={diag_a['E_min']:.4f}, "
              f"internal_min={diag_a['has_internal_min']}", flush=True)
    if diag_b:
        print(f"  B (bare+mz):    R_min={diag_b['R_min']:.3f}, E_min={diag_b['E_min']:.4f}, "
              f"internal_min={diag_b['has_internal_min']}", flush=True)
    if diag_c:
        print(f"  C (W1c):        R_min={diag_c['R_min']:.3f}, E_min={diag_c['E_min']:.4f}, "
              f"internal_min={diag_c['has_internal_min']}", flush=True)

    # Headline finding
    delta_b_vs_a = max(
        abs(pts_b[i]['E_total'] - pts_a[i]['E_total']) for i in range(len(R_grid))
    )
    print(f"\n--- HEADLINE STRUCTURAL FINDING ---", flush=True)
    print(f"  max |E_bare+mz - E_bare| across all R = {delta_b_vs_a:.3e} Ha", flush=True)
    if delta_b_vs_a < 1e-10:
        print(f"  Multi-zeta substitution has BIT-ZERO effect on the FCI energy.", flush=True)
        print(f"  The 2-electron FCI ground state at NaH max_n=2 does NOT live on", flush=True)
        print(f"  the Na 3s orbital; it lives on the H side (overattracted by Na Z=11).", flush=True)
        print(f"  M-Y Path A is the wrong axis on this basis.", flush=True)

    out = {
        'metadata': {
            'sprint': 'alpha-PES Step 3',
            'date': '2026-05-23',
            'R_grid': R_grid,
            'R_eq_exp': R_eq_exp,
            'max_n': 2,
        },
        'pts_a_bare_baseline': pts_a,
        'pts_b_bare_plus_multizeta': pts_b,
        'pts_c_w1c_screened': pts_c,
        'diag_a': diag_a,
        'diag_b': diag_b,
        'diag_c': diag_c,
        'max_abs_delta_b_vs_a_Ha': delta_b_vs_a,
        'multi_zeta_dispatch_effect': 'bit-zero on FCI energy' if delta_b_vs_a < 1e-10 else f'~{delta_b_vs_a:.3e} Ha',
    }
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] {out_path}", flush=True)


if __name__ == '__main__':
    main()
