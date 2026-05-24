"""Consistency checks for sprint_alpha_1_bimodule_distance.

Tests:
  - Sensitivity to Gaussian-probe width (1.0 vs 1.5 vs 2.0 bohr)
  - Diagnostic of d_R as a function of |r_phys - r_hyd| (shape-difference
    magnitude) rather than r_phys alone — this is the actual scaling axis.
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sprint_alpha_1_bimodule_distance import (
    build_pin_states, mean_radius
)


def bimodule_distance_LR_w(psi_a, psi_b, R_bond, width=1.0):
    """Like bimodule_distance_LR but with adjustable Gaussian-probe width."""
    n_z = 400
    z_grid = np.linspace(-1.0, R_bond + 1.0, n_z)
    r_M_at_z = np.abs(z_grid)
    r_H_at_z = np.abs(z_grid - R_bond)

    def xi_at_z(psi):
        return psi.coefficients[0] * psi.H_side(r_H_at_z) + \
               psi.coefficients[1] * psi.M_side(r_M_at_z)

    xi_a = xi_at_z(psi_a)
    xi_b = xi_at_z(psi_b)

    weight_L = np.exp(-(z_grid - R_bond)**2 / (2 * width**2))
    weight_R = np.exp(-(z_grid - 0.0)**2     / (2 * width**2))
    weight_L /= np.trapezoid(weight_L, z_grid)
    weight_R /= np.trapezoid(weight_R, z_grid)

    d_L = float(np.sqrt(np.trapezoid(weight_L * (xi_a - xi_b)**2, z_grid)))
    d_R = float(np.sqrt(np.trapezoid(weight_R * (xi_a - xi_b)**2, z_grid)))
    return d_L, d_R


def main():
    import io
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass

    print("=" * 70)
    print(" Sprint alpha-1 Consistency: probe-width and shape-difference scaling")
    print("=" * 70)

    print("\n--- Sensitivity to Gaussian-probe width ---")
    print(f"{'Mol':>6} {'width':>6}  {'d_L':>8} {'d_R':>8} {'ratio R/L':>10}")
    width_results = {}
    for M_symbol, Z_M, n_val_M, R_bond in [('Na', 11, 3, 3.566)]:
        states = build_pin_states(M_symbol, Z_M, n_val_M)
        for width in [0.5, 1.0, 1.5, 2.0, 3.0]:
            d_L, d_R = bimodule_distance_LR_w(states['i'], states['iii'],
                                              R_bond, width=width)
            ratio = d_R / d_L if d_L > 1e-12 else float('inf')
            print(f"{M_symbol+'H':>6} {width:>6.1f}  {d_L:>8.4f} {d_R:>8.4f} {ratio:>10.2f}")
            width_results[f"{M_symbol}H_w{width}"] = {
                'd_L': d_L, 'd_R': d_R, 'ratio': ratio,
            }

    # Scaling: d_R vs |r_phys - r_hyd|
    print("\n--- Scaling vs shape-difference magnitude ---")
    print(f"{'Mol':>6} {'<r>_M_phys':>12} {'<r>_M_hyd':>12} {'|r_phys-r_hyd|':>16} {'d_R(i,iii)':>12}")
    scaling = []
    for M_symbol, Z_M, n_val_M, R_bond in [
        ('Li', 3, 2, 3.015),
        ('Na', 11, 3, 3.566),
        ('K', 19, 4, 4.243),
    ]:
        states = build_pin_states(M_symbol, Z_M, n_val_M)
        r_phys = mean_radius(states['iii'].M_side)
        r_hyd = mean_radius(states['i'].M_side)
        d_L, d_R = bimodule_distance_LR_w(states['i'], states['iii'],
                                          R_bond, width=1.0)
        shape_diff = abs(r_phys - r_hyd)
        print(f"{M_symbol+'H':>6} {r_phys:>12.3f} {r_hyd:>12.3f} {shape_diff:>16.3f} {d_R:>12.4f}")
        scaling.append({
            'molecule': f"{M_symbol}H",
            'r_phys': r_phys,
            'r_hyd': r_hyd,
            'shape_diff': shape_diff,
            'd_R': d_R,
            'd_L': d_L,
        })

    # Fit log d_R vs log |shape_diff|
    print("\n--- Linear regression: d_R vs shape-difference ---")
    sd = np.array([s['shape_diff'] for s in scaling])
    dr = np.array([s['d_R']        for s in scaling])
    # Linear fit d_R = a * shape_diff + b
    A_lin = np.vstack([sd, np.ones_like(sd)]).T
    a_lin, b_lin = np.linalg.lstsq(A_lin, dr, rcond=None)[0]
    print(f"  d_R = {a_lin:.4f} * shape_diff + {b_lin:.4f}")
    # Try sqrt scaling
    A_sqrt = np.vstack([np.sqrt(sd), np.ones_like(sd)]).T
    a_sqrt, b_sqrt = np.linalg.lstsq(A_sqrt, dr, rcond=None)[0]
    print(f"  d_R = {a_sqrt:.4f} * sqrt(shape_diff) + {b_sqrt:.4f}")

    out = {
        'width_sensitivity': width_results,
        'scaling_vs_shape_diff': scaling,
        'linear_fit_a_b': [float(a_lin), float(b_lin)],
        'sqrt_fit_a_b': [float(a_sqrt), float(b_sqrt)],
    }
    out_path = os.path.join(os.path.dirname(__file__), 'data',
                              'sprint_alpha_1_consistency.json')
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n  Saved to {out_path}")


if __name__ == '__main__':
    main()
