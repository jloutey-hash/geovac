"""
B formula selection: compute all B candidates and run LiH with each.

Determines which core radius definition gives LiH R_eq closest to experiment.
"""

import numpy as np
import time

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.composed_diatomic import ComposedDiatomicSolver


def main() -> None:
    # ---- Step 1: Compute B candidates for Li and Be ----
    print("=" * 70)
    print("Step 1: B candidates for Li (Z=3) and Be (Z=4)")
    print("=" * 70)

    for Z, label in [(3, 'Li'), (4, 'Be')]:
        cs = CoreScreening(Z=Z, l_max=2, n_alpha=200)
        cs.solve(verbose=False)
        pk = AbInitioPK(cs, n_core=2)
        print(f"\n{label} (Z={Z}):")
        print(pk.summary())
        print()

    # ---- Step 2: Run LiH with each B candidate ----
    print("=" * 70)
    print("Step 2: LiH R_eq for each B formula")
    print("=" * 70)

    # First solve core once to get A and candidates
    li_core = CoreScreening(Z=3, l_max=2, n_alpha=200)
    li_core.solve(verbose=False)
    pk_ref = AbInitioPK(li_core, n_core=2)
    A_fixed = pk_ref.A
    candidates = pk_ref.B_candidates

    R_grid = np.concatenate([
        np.linspace(2.0, 2.5, 3),
        np.linspace(2.7, 4.0, 10),
        np.linspace(4.5, 7.0, 5),
    ])

    results = {}
    for name, info in candidates.items():
        B_val = info['B']
        r_val = info['r']
        t0 = time.time()
        print(f"\n  Running LiH with B_method='{name}': B={B_val:.4f}"
              f" (r={r_val:.4f})...")

        s = ComposedDiatomicSolver(
            Z_A=3, Z_B=1, n_core=2,
            M_A=7.016003, M_B=1.00782503,
            label='LiH', pk_mode='manual',
            pk_A=A_fixed, pk_B=B_val,
            l_max=2, verbose=False,
        )
        s.solve_core()
        s.scan_pes(R_grid=R_grid, n_Re=300)
        s.fit_spectroscopic_constants(fit_window=1.5)

        R_eq = s.spectro['R_eq']
        D_e = s.pes_result['D_e']
        omega_e = s.spectro['omega_e']
        dt = time.time() - t0

        results[name] = {
            'r': r_val, 'B': B_val, 'R_eq': R_eq,
            'D_e': D_e, 'omega_e': omega_e,
        }
        print(f"    R_eq={R_eq:.3f}, D_e={D_e:.4f}, omega_e={omega_e:.0f}"
              f"  ({dt:.1f}s)")

    # ---- Summary table ----
    print("\n" + "=" * 70)
    print("B formula comparison (LiH, A={:.4f} fixed)".format(A_fixed))
    print("=" * 70)
    print(f"  {'Method':8s} {'r':>8s} {'B':>8s} {'R_eq':>8s}"
          f" {'D_e':>8s} {'omega_e':>8s}")
    for name, res in results.items():
        marker = ""
        err = abs(res['R_eq'] - 3.015) / 3.015 * 100
        print(f"  {name:8s} {res['r']:8.4f} {res['B']:8.4f}"
              f" {res['R_eq']:8.3f} {res['D_e']:8.4f}"
              f" {res['omega_e']:8.0f}  ({err:.1f}% from expt)")
    print(f"  {'Expt':8s} {'---':>8s} {'---':>8s} {'3.015':>8s}"
          f" {'0.092':>8s} {'1406':>8s}")
    print("=" * 70)


if __name__ == '__main__':
    main()
