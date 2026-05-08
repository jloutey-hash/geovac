"""
Second-row PES test: NaH balanced coupled with screened W1c.

Phase C D-PES regression confirmed cross-V_ne reduction by 5.4-6.0× at
NaH equilibrium under screened_cross_center=True, but reportedly could
not determine binding because a particle-number-projected FCI driver
was not yet wired in. coupled_fci_energy IS that driver — this script
runs it.

Question: with screened W1c, does NaH balanced FCI find a bound
equilibrium R_eq, or does the PES still descend monotonically?

Q=20 for NaH at max_n=2 (frozen [Ne] core), 2-electron FCI sector.
This is a minutes-scale compute, not hours.
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


def run_single(R: float, screened: bool):
    spec = nah_spec(R=R, max_n=2)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=4000, L_max=4,
        screened_cross_center=screened, verbose=False,
    )
    fci = coupled_fci_energy(ham, n_electrons=n_e, verbose=False)
    return {
        'R': float(R),
        'screened': screened,
        'M': int(ham['M']),
        'Q': int(ham['Q']),
        'n_electrons': int(n_e),
        'E_total': float(fci['E_coupled']),
        'V_NN_plus_E_core': float(ham['nuclear_repulsion']),
    }


def main():
    out_path = Path('debug/data/chemistry_solver_retest_nah.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    R_grid = [2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0]

    R_eq_exp = 3.566   # bohr (experimental NaH ground state)

    print("=== NaH balanced FCI, max_n=2 ===", flush=True)
    print(f"  Reference R_eq = {R_eq_exp:.3f} bohr (experimental)", flush=True)

    print(f"\n--- screened=False (bare cross-V_ne, Sprint 7 baseline) ---",
          flush=True)
    pts_off = []
    for R in R_grid:
        t0 = time.perf_counter()
        try:
            pt = run_single(R, False)
            print(f"  R={R:6.3f}: E={pt['E_total']:.6f} Ha "
                  f"(V_NN+E_core={pt['V_NN_plus_E_core']:.4f}) "
                  f"[{time.perf_counter()-t0:.1f}s]", flush=True)
            pts_off.append(pt)
        except Exception as e:
            print(f"  R={R:6.3f}: FAIL {e}", flush=True)

    print(f"\n--- screened=True (W1c FrozenCore-aware) ---", flush=True)
    pts_on = []
    for R in R_grid:
        t0 = time.perf_counter()
        try:
            pt = run_single(R, True)
            print(f"  R={R:6.3f}: E={pt['E_total']:.6f} Ha "
                  f"(V_NN+E_core={pt['V_NN_plus_E_core']:.4f}) "
                  f"[{time.perf_counter()-t0:.1f}s]", flush=True)
            pts_on.append(pt)
        except Exception as e:
            print(f"  R={R:6.3f}: FAIL {e}", flush=True)

    def has_minimum(pts):
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

    diag_off = has_minimum(pts_off)
    diag_on = has_minimum(pts_on)

    print(f"\n=== Verdict ===", flush=True)
    if diag_off:
        print(f"  screened OFF: R_min={diag_off['R_min']:.3f}, "
              f"internal_min={diag_off['has_internal_min']}", flush=True)
    if diag_on:
        print(f"  screened ON:  R_min={diag_on['R_min']:.3f}, "
              f"internal_min={diag_on['has_internal_min']}", flush=True)

    # Compute cross-V_ne magnitude at R=3.5 for verification
    if pts_off and pts_on:
        E_diff_at_3 = (
            next((p['E_total'] for p in pts_off if abs(p['R'] - 3.5) < 0.01),
                 None)
            or pts_off[1]['E_total']
        )
        E_diff_screen = (
            next((p['E_total'] for p in pts_on if abs(p['R'] - 3.5) < 0.01),
                 None)
            or pts_on[1]['E_total']
        )
        print(f"  W1c effect at R=3.5: ΔE = {E_diff_at_3 - E_diff_screen:.4f} Ha "
              f"(off-on)", flush=True)

    out = {
        'R_grid': R_grid,
        'pts_off': pts_off,
        'pts_on': pts_on,
        'diag_off': diag_off,
        'diag_on': diag_on,
        'R_eq_exp': R_eq_exp,
        'sprint_7_baseline': 'monotonically descending PES, no equilibrium',
    }
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] {out_path}", flush=True)


if __name__ == '__main__':
    main()
