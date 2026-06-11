"""
H₂O coupling decomposition diagnostic.

Decomposes the inter-fiber coupling at R=1.81 bohr (experimental R_eq)
into individual pair contributions:
  - Bond 1 ↔ Bond 2  (1 pair, θ_HOH = 104.5°)
  - Bond ↔ Lone pair  (4 pairs, θ_bond_lone = 120°)
  - Lone pair 1 ↔ Lone pair 2  (1 pair, θ_lone_lone = 114°)

Then runs bond-bond-only PES scan to isolate validated physics.

Output: debug/data/h2o_coupling_decomposition.json
"""

import json
import time
import numpy as np
from pathlib import Path

from geovac.composed_water import ComposedWaterSolver


def run_decomposition():
    results = {}

    # ---------------------------------------------------------------
    # Part 1: Coupling decomposition at R = 1.81 bohr
    # ---------------------------------------------------------------
    print("=" * 64)
    print("Part 1: Coupling decomposition at R = 1.81 bohr")
    print("=" * 64)

    solver = ComposedWaterSolver(
        l_max=2, n_alpha=100,
        include_coupling=True,
        coupling_pairs='all',
        verbose=True,
    )
    solver.solve_core()
    solver.solve_lone_pair()

    R_eq = 1.81  # experimental

    # Solve bond pair at R_eq (need full result for coupling)
    bond_result = solver._solve_bond_at_R(R_eq, n_Re=300, return_full=True)
    E_bond = bond_result['E_elec']
    print(f"\n  Bond pair energy at R={R_eq}: {E_bond:.6f} Ha")

    # Compute all 6 couplings individually
    coupling = solver._compute_all_coupling(R_eq, bond_result)

    print(f"\n  Coupling decomposition at R = {R_eq} bohr:")
    print(f"    Bond-Bond (1 pair):        {coupling['E_bond_bond']:.6f} Ha")
    print(f"    Bond-Lone (1 pair):        {coupling['E_bond_lone_one']:.6f} Ha")
    print(f"    Bond-Lone (4 pairs total): {coupling['E_bond_lone_total']:.6f} Ha")
    print(f"    Lone-Lone (1 pair):        {coupling['E_lone_lone']:.6f} Ha")
    print(f"    TOTAL (6 pairs):           {coupling['E_coupling_total']:.6f} Ha")

    # For comparison: compute at a few more R values
    R_decompose = [1.4, 1.6, 1.81, 2.0, 2.4, 2.8, 3.2, 4.0]
    decomposition_data = []

    print(f"\n  Coupling decomposition across R:")
    print(f"  {'R':>6s}  {'Bond-Bond':>10s}  {'4×Bond-Lone':>12s}"
          f"  {'Lone-Lone':>10s}  {'Total':>10s}")
    print(f"  {'-'*6}  {'-'*10}  {'-'*12}  {'-'*10}  {'-'*10}")

    for R in R_decompose:
        try:
            br = solver._solve_bond_at_R(R, n_Re=300, return_full=True)
            c = solver._compute_all_coupling(R, br)
            decomposition_data.append({
                'R': R,
                'E_bond_bond': c['E_bond_bond'],
                'E_bond_lone_one': c['E_bond_lone_one'],
                'E_bond_lone_total': c['E_bond_lone_total'],
                'E_lone_lone': c['E_lone_lone'],
                'E_coupling_total': c['E_coupling_total'],
            })
            print(f"  {R:6.3f}  {c['E_bond_bond']:10.4f}"
                  f"  {c['E_bond_lone_total']:12.4f}"
                  f"  {c['E_lone_lone']:10.4f}"
                  f"  {c['E_coupling_total']:10.4f}")
        except Exception as e:
            print(f"  {R:6.3f}  FAILED: {e}")

    results['decomposition'] = decomposition_data
    results['R_eq_expt'] = R_eq

    # ---------------------------------------------------------------
    # Part 2: Bond-bond-only PES scan
    # ---------------------------------------------------------------
    print("\n" + "=" * 64)
    print("Part 2: Bond-bond-only PES scan")
    print("=" * 64)

    solver_bb = ComposedWaterSolver(
        l_max=2, n_alpha=100,
        include_coupling=True,
        coupling_pairs='bond_bond',
        verbose=True,
    )
    solver_bb.solve_core()
    solver_bb.solve_lone_pair()

    R_grid = np.arange(1.2, 4.5, 0.2)
    pes_bb = solver_bb.scan_pes(R_grid=R_grid, n_Re=300)
    spectro_bb = solver_bb.fit_spectroscopic_constants()

    results['bond_bond_only'] = {
        'pes': pes_bb,
        'spectro': spectro_bb,
    }

    # ---------------------------------------------------------------
    # Part 3: Uncoupled PES scan for comparison
    # ---------------------------------------------------------------
    print("\n" + "=" * 64)
    print("Part 3: Uncoupled PES scan (reference)")
    print("=" * 64)

    solver_uc = ComposedWaterSolver(
        l_max=2, n_alpha=100,
        include_coupling=False,
        verbose=True,
    )
    solver_uc.solve_core()
    solver_uc.solve_lone_pair()

    pes_uc = solver_uc.scan_pes(R_grid=R_grid, n_Re=300)
    spectro_uc = solver_uc.fit_spectroscopic_constants()

    results['uncoupled'] = {
        'pes': pes_uc,
        'spectro': spectro_uc,
    }

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    print("\n" + "=" * 64)
    print("SUMMARY")
    print("=" * 64)

    ref_R = 1.809

    R_eq_uc = spectro_uc['R_eq']
    R_eq_bb = spectro_bb['R_eq']
    err_uc = abs(R_eq_uc - ref_R) / ref_R * 100
    err_bb = abs(R_eq_bb - ref_R) / ref_R * 100

    print(f"\n  Reference R_eq:      {ref_R:.3f} bohr")
    print(f"  Uncoupled R_eq:      {R_eq_uc:.3f} bohr  ({err_uc:.1f}% error)")
    print(f"  Bond-bond only R_eq: {R_eq_bb:.3f} bohr  ({err_bb:.1f}% error)")

    # Identify dominant pair
    if decomposition_data:
        r181 = [d for d in decomposition_data if abs(d['R'] - 1.81) < 0.01]
        if r181:
            d = r181[0]
            total = abs(d['E_coupling_total'])
            bb_frac = abs(d['E_bond_bond']) / total * 100
            bl_frac = abs(d['E_bond_lone_total']) / total * 100
            ll_frac = abs(d['E_lone_lone']) / total * 100
            print(f"\n  At R=1.81 bohr:")
            print(f"    Bond-Bond:  {d['E_bond_bond']:.4f} Ha  ({bb_frac:.1f}%)")
            print(f"    Bond-Lone:  {d['E_bond_lone_total']:.4f} Ha  ({bl_frac:.1f}%)")
            print(f"    Lone-Lone:  {d['E_lone_lone']:.4f} Ha  ({ll_frac:.1f}%)")
            print(f"    Total:      {d['E_coupling_total']:.4f} Ha")

    results['summary'] = {
        'ref_R_eq': ref_R,
        'uncoupled_R_eq': float(R_eq_uc),
        'uncoupled_R_eq_error_pct': float(err_uc),
        'bond_bond_only_R_eq': float(R_eq_bb),
        'bond_bond_only_R_eq_error_pct': float(err_bb),
    }

    # Save
    out_path = Path(__file__).parent / 'data' / 'h2o_coupling_decomposition.json'
    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved to {out_path}")


if __name__ == '__main__':
    run_decomposition()
