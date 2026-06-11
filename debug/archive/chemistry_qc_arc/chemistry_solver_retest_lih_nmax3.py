"""
LiH max_n=3 single-point + small grid scan.

Builds on debug/chemistry_solver_retest_lih.py. Runs a focused max_n=3 study
to determine drift slope between Track CD's reported (n_max=2: R_eq=3.226,
n_max=3: R_eq=3.280) and the current production architecture.
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
from geovac.molecular_spec import lih_spec, _FIRST_ROW_CORE_ENERGY


def run_single(R: float, max_n: int, screened: bool):
    t0 = time.perf_counter()
    spec = lih_spec(R=R, max_n=max_n)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=4000, L_max=4,
        screened_cross_center=screened, verbose=False,
    )
    t_build = time.perf_counter() - t0
    t0 = time.perf_counter()
    fci = coupled_fci_energy(ham, n_electrons=n_e, verbose=False)
    t_fci = time.perf_counter() - t0
    return {
        'R': R, 'max_n': max_n, 'screened': screened,
        'M': int(ham['M']), 'Q': int(ham['Q']),
        'n_electrons': int(n_e),
        'E_total': float(fci['E_coupled']),
        'V_NN_plus_E_core': float(ham['nuclear_repulsion']),
        'cross_block_eri_count': int(ham.get('cross_block_eri_count', 0)),
        'cross_vne_count': int(ham.get('cross_vne_count', 0)),
        't_build_s': float(t_build),
        't_fci_s': float(t_fci),
    }


def main():
    out_path = Path('debug/data/chemistry_solver_retest_lih_nmax3.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    E_core_off = _FIRST_ROW_CORE_ENERGY.get(3, -(3 - 5/16) ** 2)
    print(f"E_core offset (Li): {E_core_off:.4f} Ha")

    # Single point timing test at R = 3.015
    print("\n=== max_n=3 timing test at R=3.015 ===", flush=True)
    t0 = time.perf_counter()
    pt0 = run_single(3.015, 3, screened=True)
    print(f"  M={pt0['M']}, Q={pt0['Q']}, n_e={pt0['n_electrons']}", flush=True)
    print(f"  E_total = {pt0['E_total']:.6f} Ha "
          f"(TC-conv: {pt0['E_total'] - E_core_off:.6f})", flush=True)
    print(f"  t_build={pt0['t_build_s']:.1f}s, t_fci={pt0['t_fci_s']:.1f}s",
          flush=True)

    elapsed = time.perf_counter() - t0
    points = [pt0]

    # Decide how many additional points based on per-point time
    if elapsed < 90:
        R_extras = [2.7, 3.3, 3.6, 3.0]
    elif elapsed < 240:
        R_extras = [2.7, 3.3]
    else:
        R_extras = []

    for R in R_extras:
        print(f"\n  R={R:.3f} ...", flush=True)
        t0 = time.perf_counter()
        try:
            pt = run_single(R, 3, screened=True)
            print(f"    E={pt['E_total']:.6f} Ha "
                  f"(TC-conv: {pt['E_total'] - E_core_off:.6f}) "
                  f"[t_build={pt['t_build_s']:.1f}s, "
                  f"t_fci={pt['t_fci_s']:.1f}s]", flush=True)
            points.append(pt)
        except Exception as e:
            print(f"    FAIL: {type(e).__name__}: {e}", flush=True)
            points.append({'R': R, 'error': str(e)})

    # Parabolic fit
    R_arr = np.array([p['R'] for p in points if 'E_total' in p])
    E_arr = np.array([p['E_total'] for p in points if 'E_total' in p])

    R_eq = float('nan')
    E_eq = float('nan')
    if len(R_arr) >= 3:
        # Sort by R, find minimum
        idx = np.argsort(R_arr)
        R_arr = R_arr[idx]
        E_arr = E_arr[idx]
        i_min = int(np.argmin(E_arr))
        if 0 < i_min < len(R_arr) - 1:
            R3 = R_arr[i_min - 1: i_min + 2]
            E3 = E_arr[i_min - 1: i_min + 2]
            p = np.polyfit(R3, E3, 2)
            R_eq = -p[1] / (2 * p[0])
            E_eq = p[0] * R_eq ** 2 + p[1] * R_eq + p[2]
            print(f"\n  Parabolic fit: R_eq = {R_eq:.4f} bohr, "
                  f"E_min = {E_eq:.6f} Ha "
                  f"(TC-conv: {E_eq - E_core_off:.6f})", flush=True)
        else:
            print(f"\n  Minimum at edge: R = {R_arr[i_min]:.3f}, "
                  f"E = {E_arr[i_min]:.6f} Ha", flush=True)

    result = {
        'points': points,
        'R_eq_fit': R_eq,
        'E_min_fit': E_eq,
        'E_min_track_cd_conv': E_eq - E_core_off,
        'R_arr': [float(x) for x in R_arr],
        'E_arr': [float(x) for x in E_arr],
        'baseline_track_cd_nmax3': {'R_eq': 3.280, 'E_min': -8.055},
    }

    if not np.isnan(R_eq):
        result['R_err_pct'] = abs(R_eq - 3.015) / 3.015 * 100
        result['E_err_pct'] = abs((E_eq - E_core_off) - (-8.071)) / 8.071 * 100

    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2, default=str)
    print(f"\n[saved] {out_path}", flush=True)


if __name__ == '__main__':
    main()
