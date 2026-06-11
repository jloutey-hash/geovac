"""
H₂O PES scan with charge_center origin.

Compares to previous midpoint results:
  - Uncoupled: R_eq = 2.72 bohr
  - Coupled:   R_eq = 2.56 bohr
  - Reference:  R_eq = 1.809 bohr

Uses same R-grid as original: 1.4 to 4.4 bohr, 16 points.
"""

import json
import time
import numpy as np
from geovac.composed_water import ComposedWaterSolver


def run_scan(include_coupling: bool, label: str) -> dict:
    """Run H₂O PES scan and return results."""
    solver = ComposedWaterSolver(
        l_max=2,
        n_alpha=100,
        include_coupling=include_coupling,
        verbose=True,
    )

    R_grid = np.linspace(1.4, 4.4, 16)

    result = solver.run_all(R_grid=R_grid, n_Re=300)

    return {
        'label': label,
        'origin': 'charge_center',
        'include_coupling': include_coupling,
        'R_eq_grid': result['pes']['R_eq'],
        'R_eq_fit': result['spectro']['R_eq'],
        'E_min': result['pes']['E_min'],
        'D_e': result['pes']['D_e'],
        'R': result['pes']['R'],
        'E_total': result['pes']['E_total'],
    }


if __name__ == '__main__':
    results = {}

    print("=" * 64)
    print("H₂O PES SCAN: charge_center origin")
    print("=" * 64)

    # Uncoupled scan
    t0 = time.time()
    results['uncoupled'] = run_scan(False, 'uncoupled_charge_center')
    t_uncoupled = time.time() - t0

    # Coupled scan
    t0 = time.time()
    results['coupled'] = run_scan(True, 'coupled_charge_center')
    t_coupled = time.time() - t0

    # Summary
    print("\n" + "=" * 64)
    print("COMPARISON SUMMARY")
    print("=" * 64)
    print(f"\n{'':20s} {'Previous':>12s} {'New (CC)':>12s} {'Reference':>12s}")
    print(f"{'':20s} {'(midpoint)':>12s} {'(charge ctr)':>12s}")
    print("-" * 58)

    R_eq_unc = results['uncoupled']['R_eq_fit']
    R_eq_coup = results['coupled']['R_eq_fit']
    ref_R = 1.809

    print(f"{'R_eq uncoupled':20s} {'2.72':>12s} {R_eq_unc:12.3f} {ref_R:12.3f}")
    print(f"{'R_eq coupled':20s} {'2.56':>12s} {R_eq_coup:12.3f} {ref_R:12.3f}")

    err_old_unc = abs(2.72 - ref_R) / ref_R * 100
    err_new_unc = abs(R_eq_unc - ref_R) / ref_R * 100
    err_old_coup = abs(2.56 - ref_R) / ref_R * 100
    err_new_coup = abs(R_eq_coup - ref_R) / ref_R * 100

    print(f"\n{'Err uncoupled':20s} {err_old_unc:11.1f}% {err_new_unc:11.1f}%")
    print(f"{'Err coupled':20s} {err_old_coup:11.1f}% {err_new_coup:11.1f}%")

    print(f"\nTimings: uncoupled={t_uncoupled:.0f}s, coupled={t_coupled:.0f}s")

    # Save
    out = {
        'origin': 'charge_center',
        'reference_R_eq': ref_R,
        'uncoupled': results['uncoupled'],
        'coupled': results['coupled'],
        'previous_midpoint': {
            'R_eq_uncoupled': 2.72,
            'R_eq_coupled': 2.56,
        },
    }
    with open('debug/data/h2o_charge_center_scan.json', 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print("\nSaved to debug/data/h2o_charge_center_scan.json")
