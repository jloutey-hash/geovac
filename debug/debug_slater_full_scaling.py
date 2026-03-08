#!/usr/bin/env python3
"""
debug_slater_full_scaling.py
============================
Scaling study for slater_full FCI: He and Li at max_n = 2..5.
Generates data for publication figures.

Author: GeoVac Development Team, March 2026
"""

import os
import sys
import time
import json
import warnings

import numpy as np
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

HE_EXACT_HA: float = -2.9037
LI_EXACT_HA: float = -7.4781

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'debug', 'data')
os.makedirs(DATA_DIR, exist_ok=True)


def run_single(
    label: str,
    n_electrons: int,
    nuclear_charge: int,
    max_n: int,
    exact_energy: float,
    vee_method: str = 'slater_full',
    h1_method: str = 'hybrid',
) -> dict:
    """Run a single FCI audit."""
    from geovac.lattice_index import LatticeIndex

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t0 = time.perf_counter()
        idx = LatticeIndex(
            n_electrons=n_electrons,
            max_n=max_n,
            nuclear_charge=nuclear_charge,
            vee_method=vee_method,
            h1_method=h1_method,
        )
        H = idx.assemble_hamiltonian()
        t_asm = time.perf_counter() - t0

    k = min(1, idx.n_sd - 2)
    t_eig = time.perf_counter()
    eigvals_raw, _ = eigsh(H, k=k, which="SA")
    t_eig = time.perf_counter() - t_eig
    E = float(eigvals_raw[np.argsort(eigvals_raw)][0])

    abs_err = abs(E - exact_energy)
    pct_err = abs_err / abs(exact_energy) * 100.0

    result = {
        "label": label,
        "Z": nuclear_charge,
        "n_el": n_electrons,
        "max_n": max_n,
        "n_sd": idx.n_sd,
        "nnz": H.nnz,
        "E_geovac": E,
        "E_exact": exact_energy,
        "abs_err": abs_err,
        "pct_err": pct_err,
        "t_assembly": t_asm,
        "t_eigsolve": t_eig,
        "t_total": t_asm + t_eig,
        "vee_method": vee_method,
        "h1_method": h1_method,
        "n_eri": len(idx._eri) if hasattr(idx, '_eri') else 0,
    }
    print(f"  {label:>25}  max_n={max_n}  n_sd={idx.n_sd:>8,}  "
          f"E={E:>12.6f}  err={pct_err:>6.2f}%  "
          f"t={t_asm:.1f}+{t_eig:.1f}s  nnz={H.nnz:,}")
    return result


def main() -> None:
    print("=" * 80)
    print("SLATER_FULL FCI SCALING STUDY")
    print("=" * 80)

    all_results = []

    # --- He convergence ---
    print("\n--- Helium (Z=2, 2 electrons) ---")
    for nmax in [2, 3, 4, 5]:
        try:
            r = run_single(f"He slater_full+hybrid", 2, 2, nmax, HE_EXACT_HA,
                           'slater_full', 'hybrid')
            all_results.append(r)
        except Exception as e:
            print(f"  He max_n={nmax} FAILED: {e}")
            break

    # He with exact h1 for comparison
    print("\n--- Helium (exact h1 comparison) ---")
    for nmax in [2, 3, 4]:
        try:
            r = run_single(f"He slater_full+exact", 2, 2, nmax, HE_EXACT_HA,
                           'slater_full', 'exact')
            all_results.append(r)
        except Exception as e:
            print(f"  He exact max_n={nmax} FAILED: {e}")
            break

    # He with slater F0 (baseline)
    print("\n--- Helium (slater F0 baseline) ---")
    for nmax in [2, 3, 4]:
        try:
            r = run_single(f"He slater+hybrid", 2, 2, nmax, HE_EXACT_HA,
                           'slater', 'hybrid')
            all_results.append(r)
        except Exception as e:
            print(f"  He slater max_n={nmax} FAILED: {e}")
            break

    # --- Li convergence ---
    print("\n--- Lithium (Z=3, 3 electrons) ---")
    for nmax in [2, 3, 4, 5]:
        try:
            r = run_single(f"Li slater_full+hybrid", 3, 3, nmax, LI_EXACT_HA,
                           'slater_full', 'hybrid')
            all_results.append(r)
        except Exception as e:
            print(f"  Li max_n={nmax} FAILED: {e}")
            break

    # Li with exact h1 for comparison
    print("\n--- Lithium (exact h1 comparison) ---")
    for nmax in [2, 3]:
        try:
            r = run_single(f"Li slater_full+exact", 3, 3, nmax, LI_EXACT_HA,
                           'slater_full', 'exact')
            all_results.append(r)
        except Exception as e:
            print(f"  Li exact max_n={nmax} FAILED: {e}")
            break

    # --- Summary table ---
    print(f"\n{'=' * 80}")
    print("CONVERGENCE SUMMARY")
    print(f"{'=' * 80}")
    print(f"{'Method':>28}  {'n':>3}  {'n_sd':>8}  {'E(Ha)':>12}  "
          f"{'err%':>7}  {'nnz':>10}  {'t(s)':>7}")
    print("-" * 80)
    for r in all_results:
        print(f"{r['label']:>28}  {r['max_n']:>3}  {r['n_sd']:>8,}  "
              f"{r['E_geovac']:>12.6f}  {r['pct_err']:>6.2f}%  "
              f"{r['nnz']:>10,}  {r['t_total']:>6.1f}")

    # Save results to JSON
    out_path = os.path.join(DATA_DIR, 'slater_full_scaling.json')
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {out_path}")
    print("[DONE]")


if __name__ == "__main__":
    main()
