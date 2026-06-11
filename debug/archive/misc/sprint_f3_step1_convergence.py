"""Sprint F3 Step 1 convergence check on quadrature resolution.

Tests how the cross-block h1 matrix elements converge under increasing
(n_rho, n_z) and (rho_max, z_max) values for the hydrogenic Na 3s + H 1s
diagnostic at NaH R = 3.566 bohr.
"""

from __future__ import annotations

import json
import sys
import time

sys.path.insert(0, '.')

from debug.sprint_f3_step1_diagnostic import (
    step1_diagnostic, sanity_check_h1s_normalization,
)


def main():
    settings = [
        (80, 100, 20.0, 20.0),
        (160, 200, 30.0, 30.0),
        (250, 300, 40.0, 40.0),
        (400, 400, 50.0, 50.0),
    ]

    print("Sprint F3 Step 1 convergence check")
    print("=" * 78)
    print(f"{'n_rho':>6} {'n_z':>5} {'rho_max':>8} {'z_max':>7} "
          f"{'<H1s|H1s>':>11} {'h12 (Ha)':>11} {'splitting (Ha)':>14}")
    print("-" * 78)
    results = []
    for n_rho, n_z, rho_max, z_max in settings:
        t0 = time.time()
        S_HH = sanity_check_h1s_normalization(
            n_rho=n_rho, n_z=n_z, rho_max=rho_max, z_max=z_max,
        )
        r = step1_diagnostic(
            R_AB=3.566, n_rho=n_rho, n_z=n_z, rho_max=rho_max, z_max=z_max,
            use_multi_zeta_for_na3s=False,
        )
        elapsed = time.time() - t0
        h12 = r['h1_2x2_matrix']['h12_offdiag_Ha']
        split = r['eigenvalues_orthogonal_basis']['splitting_Ha']
        print(f"{n_rho:>6} {n_z:>5} {rho_max:>8.1f} {z_max:>7.1f} "
              f"{S_HH:>11.7f} {h12:>+11.5f} {split:>+14.5f}  (t={elapsed:.1f}s)")
        results.append({
            'n_rho': n_rho, 'n_z': n_z, 'rho_max': rho_max, 'z_max': z_max,
            'S_HH': S_HH, 'h12': h12, 'splitting': split,
            'h11': r['h1_2x2_matrix']['h11_Na3s_diag_Ha'],
            'h22': r['h1_2x2_matrix']['h22_H1s_diag_Ha'],
            'overlap_S_AB': r['matrix_elements']['overlap_S_AB'],
            'elapsed_s': elapsed,
        })

    with open("debug/data/sprint_f3_step1_convergence.json", "w") as f:
        json.dump({'settings_tested': results}, f, indent=2)


if __name__ == "__main__":
    main()
