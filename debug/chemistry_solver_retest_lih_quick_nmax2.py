"""
Recompute LiH max_n=2 PES (5 points) with full result-saving and parabolic
fit, then save consolidated JSON. Reproduces from the killed-run results.
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


def run_single(R, max_n, screened):
    spec = lih_spec(R=R, max_n=max_n)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=4000, L_max=4,
        screened_cross_center=screened, verbose=False,
    )
    fci = coupled_fci_energy(ham, n_electrons=n_e, verbose=False)
    return float(fci['E_coupled'])


def parab_fit(R_arr, E_arr):
    idx = np.argsort(R_arr)
    R_arr = np.array(R_arr)[idx]
    E_arr = np.array(E_arr)[idx]
    i = int(np.argmin(E_arr))
    if 0 < i < len(R_arr) - 1:
        p = np.polyfit(R_arr[i - 1: i + 2], E_arr[i - 1: i + 2], 2)
        Re = -p[1] / (2 * p[0])
        Ee = p[0] * Re * Re + p[1] * Re + p[2]
        return float(Re), float(Ee)
    return float(R_arr[i]), float(E_arr[i])


def main():
    E_core_off = _FIRST_ROW_CORE_ENERGY.get(3, -(3 - 5/16) ** 2)
    print(f"E_core offset (Li): {E_core_off:.4f} Ha", flush=True)

    R_grid = [2.4, 2.7, 3.0, 3.3, 3.6]
    print(f"\n=== LiH max_n=2, screened=False ===", flush=True)
    E_off = []
    for R in R_grid:
        t0 = time.perf_counter()
        E = run_single(R, 2, False)
        print(f"  R={R:.3f}: E={E:.6f} Ha ({time.perf_counter()-t0:.1f}s)",
              flush=True)
        E_off.append(E)

    print(f"\n=== LiH max_n=2, screened=True ===", flush=True)
    E_on = []
    for R in R_grid:
        t0 = time.perf_counter()
        E = run_single(R, 2, True)
        print(f"  R={R:.3f}: E={E:.6f} Ha ({time.perf_counter()-t0:.1f}s)",
              flush=True)
        E_on.append(E)

    # W1c bit-identical?
    diff = max(abs(a - b) for a, b in zip(E_off, E_on))
    print(f"\nW1c max abs diff (LiH first-row): {diff:.2e} Ha", flush=True)

    # Parabolic R_eq fit
    Re_off, Ee_off = parab_fit(R_grid, E_off)
    Re_on, Ee_on = parab_fit(R_grid, E_on)
    print(f"\n  screened OFF: R_eq = {Re_off:.4f} bohr, "
          f"E_min = {Ee_off:.6f} Ha (TC-conv: {Ee_off-E_core_off:.6f})",
          flush=True)
    print(f"  screened ON:  R_eq = {Re_on:.4f} bohr, "
          f"E_min = {Ee_on:.6f} Ha (TC-conv: {Ee_on-E_core_off:.6f})",
          flush=True)

    R_err = abs(Re_off - 3.015) / 3.015 * 100
    E_err = abs((Ee_off - E_core_off) - (-8.071)) / 8.071 * 100
    print(f"\n  R_eq error vs experiment (3.015): {R_err:.3f}%", flush=True)
    print(f"  E_min error vs ref (-8.071): {E_err:.3f}%", flush=True)

    out = {
        'R_grid': R_grid,
        'E_screened_off': E_off,
        'E_screened_on': E_on,
        'w1c_max_abs_diff': float(diff),
        'w1c_bit_identical': diff < 1e-10,
        'fit_off': {'R_eq': Re_off, 'E_min': Ee_off,
                    'E_min_tc_conv': Ee_off - E_core_off},
        'fit_on': {'R_eq': Re_on, 'E_min': Ee_on,
                   'E_min_tc_conv': Ee_on - E_core_off},
        'R_err_pct': R_err,
        'E_err_pct': E_err,
        'experimental_R_eq': 3.015,
        'reference_E_min': -8.071,
        'track_cd_baseline_nmax2': {'R_eq': 3.226, 'E_min': -7.924,
                                     'R_err_pct': 7.0, 'E_err_pct': 1.8},
    }

    out_path = Path('debug/data/chemistry_solver_retest_lih_nmax2.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] {out_path}", flush=True)


if __name__ == '__main__':
    main()
