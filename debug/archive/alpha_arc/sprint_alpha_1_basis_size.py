"""Basis-size consistency check for the bimodule distance diagnostic.

Tests whether d_R/d_L for NaH (i)-vs-(iii) is stable under variation of the
FrozenCore radial-solver grid (n_grid in _solve_screened_radial) and the
physical-valence interpolation onto the diagnostic 1D bond grid.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def main():
    import io
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass

    print("=" * 70)
    print(" Sprint alpha-1 basis-size consistency: NaH (i) vs (iii)")
    print("=" * 70)

    from geovac.neon_core import _solve_screened_radial

    # Vary FrozenCore solver grid size
    results = []
    for n_grid in [4000, 8000, 12000, 16000]:
        energy, u_vec, r_grid = _solve_screened_radial(
            Z=11, l=0, n_target=3,
            core_type='Ne',
            n_grid=n_grid, r_max=80.0,
            allow_l0=True,
        )
        # Compute mean radius
        R_vec = np.where(r_grid > 0, u_vec / r_grid, 0.0)
        norm_sq = np.trapezoid(R_vec**2 * r_grid**2, r_grid)
        if norm_sq > 0:
            R_vec /= np.sqrt(norm_sq)
        r_mean = float(np.trapezoid(R_vec**2 * r_grid * r_grid**2, r_grid))
        # Compute d_R using the consistency-script function
        from sprint_alpha_1_consistency import bimodule_distance_LR_w
        from sprint_alpha_1_bimodule_distance import (
            build_pin_states, hydrogenic_1s
        )

        # Build pin states using current grid
        from scipy.special import genlaguerre

        def H_hyd_Z1(r):
            return hydrogenic_1s(1.0, r)

        def M_hyd_Z1(r):
            n = 3
            Z_eff = 1.0
            rho = 2.0 * Z_eff * r / n
            L_poly = genlaguerre(n - 1, 1)(rho)
            wf = np.exp(-rho / 2.0) * L_poly
            norm_sq = np.trapezoid(wf**2 * r**2, r)
            if norm_sq > 0:
                wf = wf / np.sqrt(norm_sq)
            return wf

        def M_physical(r):
            return np.interp(r, r_grid, R_vec, left=0.0, right=0.0)

        # Build PinState-like objects
        from dataclasses import dataclass
        from typing import Callable, Tuple
        @dataclass
        class P:
            H_side: Callable
            M_side: Callable
            coefficients: Tuple[float, float]

        coef = (1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0))
        psi_i = P(H_hyd_Z1, M_hyd_Z1, coef)
        psi_iii = P(H_hyd_Z1, M_physical, coef)

        d_L, d_R = bimodule_distance_LR_w(psi_i, psi_iii, 3.566, width=1.0)
        ratio = d_R / d_L if d_L > 1e-12 else float('inf')
        results.append({
            'n_grid': n_grid,
            'energy_Na_3s': float(energy),
            'r_mean_Na_3s': r_mean,
            'd_L': d_L,
            'd_R': d_R,
            'ratio': ratio,
        })
        print(f"  n_grid={n_grid:5d}: E_3s={energy:.5f} Ha, <r>={r_mean:.4f}, "
              f"d_L={d_L:.4f}, d_R={d_R:.4f}, ratio={ratio:.2f}")

    out_path = os.path.join(os.path.dirname(__file__), 'data',
                              'sprint_alpha_1_basis_size.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved to {out_path}")


if __name__ == '__main__':
    main()
