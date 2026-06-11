"""
Track EP-2c: Paper 27 Prediction 2 pilot — entropy vs V_ee off-diagonal mass.

Paper 27 §VII.B predicts S_B ~ w_B^alpha where
  S_B = spatial 1-RDM von Neumann entropy of the block ground state,
  w_B = ||off-diag_{H_1 eigenbasis}(V_ee)||_F (the 6-8% cross-shell
        residue that EP-1 identified as the cusp correlation).

Pilot: single-block He-like isoelectronic sequence at n_max=3, varying
Z in {2, 3, 4, 6, 8, 10}. This is the cleanest one-block baseline
(uniform geometry, uniform basis, zero PK) and reuses the EP-1
pipeline with float-arithmetic Slater integrals (hypergeometric_slater).

Confirmation: tight monotone S_B(w_B), recoverable power-law fit.
Refutation: scatter incompatible with a single-variable relation.

Output: debug/data/ep2c_entropy_vs_vee_offdiag.json
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.energy_entanglement_decoupling import (
    build_decomposed_hamiltonians,
    solve_and_entangle,
)


def _off_diag_mass_in_h1_eigenbasis(H_vee: np.ndarray,
                                    H_kin: np.ndarray) -> dict:
    """Return w_B = ||off-diag_{H_1 eigenbasis}(V_ee)||_F and companions."""
    eH, U = np.linalg.eigh(H_kin)
    V_in = U.T @ H_vee @ U
    V_frob = float(np.linalg.norm(V_in))
    V_diag_frob = float(np.sqrt(np.sum(np.diag(V_in) ** 2)))
    V_off_frob = float(np.sqrt(max(0.0, V_frob ** 2 - V_diag_frob ** 2)))
    return {
        'w_B': V_off_frob,
        'V_frobenius_total': V_frob,
        'V_diag_frobenius': V_diag_frob,
        'V_diag_fraction_unsquared': V_diag_frob / V_frob if V_frob > 0 else 0.0,
    }


def run_point(Z: float, n_max: int = 3) -> dict:
    t0 = time.time()
    data = build_decomposed_hamiltonians(Z, n_max)
    configs = data['configs']
    n_spatial = data['n_spatial']

    H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
    H_vee = data['H_vee_full']
    H_full = H_kin + H_vee

    E_full, _, _, ent = solve_and_entangle(H_full, configs, n_spatial)
    S_B = float(ent['von_neumann_entropy'])
    occ = ent['occupation_numbers']

    off_info = _off_diag_mass_in_h1_eigenbasis(H_vee, H_kin)
    H_kin_frob = float(np.linalg.norm(H_kin))
    w_B_dimless = off_info['w_B'] / H_kin_frob if H_kin_frob > 0 else 0.0

    return {
        'Z': float(Z),
        'n_max': n_max,
        'n_configs': data['n_configs'],
        'E_full': float(E_full),
        'S_B': S_B,
        'w_B': off_info['w_B'],
        'w_B_dimless': w_B_dimless,  # w_B / ||H_kin||
        'H_kin_frobenius': H_kin_frob,
        'V_frobenius': off_info['V_frobenius_total'],
        'V_diag_fraction_unsquared': off_info['V_diag_fraction_unsquared'],
        'occ_top4': [float(x) for x in occ[:4]],
        'wall_seconds': time.time() - t0,
    }


def _fit_power_law(xs, ys):
    """Log-log least squares: log y = alpha * log x + const."""
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    mask = (xs > 0) & (ys > 0)
    if mask.sum() < 2:
        return None
    lx = np.log(xs[mask]); ly = np.log(ys[mask])
    alpha, intercept = np.polyfit(lx, ly, 1)
    pred = alpha * lx + intercept
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return {
        'alpha': float(alpha),
        'intercept': float(intercept),
        'r2': float(r2),
        'A': float(np.exp(intercept)),
    }


def main():
    Z_list = [2.0, 3.0, 4.0, 6.0, 8.0, 10.0]
    print("=" * 70)
    print("Track EP-2c: entropy S_B vs V_ee off-diagonal mass w_B (pilot)")
    print("He-like isoelectronic sequence, n_max=3, float Slater integrals")
    print("=" * 70)

    results = []
    for Z in Z_list:
        r = run_point(Z)
        results.append(r)
        print(f"  Z = {Z:4.1f}  S_B = {r['S_B']:.5e}  "
              f"w_B = {r['w_B']:.5e}  "
              f"V_diag_frac = {r['V_diag_fraction_unsquared']:.4f}  "
              f"({r['wall_seconds']:.1f}s)")

    # Fit S_B = A * w_B^alpha  (raw mass)
    S_vals = [r['S_B'] for r in results]
    w_vals = [r['w_B'] for r in results]
    w_dim_vals = [r['w_B_dimless'] for r in results]
    fit = _fit_power_law(w_vals, S_vals)

    # Also Z-scaling: S_B ~ Z^(-alpha_Z) (Paper 26 reported -2.56)
    Z_vals = [r['Z'] for r in results]
    fit_Z = _fit_power_law(Z_vals, S_vals)

    # Dimensionless: S_B vs w_B / ||H_kin||
    fit_dim = _fit_power_law(w_dim_vals, S_vals)

    print("\n" + "-" * 70)
    if fit:
        print(f"Power-law fit S_B = A * w_B^alpha:")
        print(f"  alpha = {fit['alpha']:.4f}")
        print(f"  A     = {fit['A']:.4e}")
        print(f"  R^2   = {fit['r2']:.6f}")
    if fit_Z:
        print(f"\nZ-scaling fit S_B = A_Z * Z^alpha_Z:")
        print(f"  alpha_Z = {fit_Z['alpha']:.4f}  (Paper 26 reported -2.56)")
        print(f"  R^2     = {fit_Z['r2']:.6f}")
    if fit_dim:
        print(f"\nDimensionless fit S_B = A * (w_B/||H_kin||)^alpha:")
        print(f"  alpha = {fit_dim['alpha']:.4f}")
        print(f"  A     = {fit_dim['A']:.4e}")
        print(f"  R^2   = {fit_dim['r2']:.6f}")

    # Verdict
    print("\nVerdict:")
    if fit and fit['r2'] > 0.98:
        print(f"  CONFIRMED: tight power law S_B ~ w_B^{fit['alpha']:.3f} "
              f"(R^2 = {fit['r2']:.4f})")
    elif fit and fit['r2'] > 0.9:
        print(f"  WEAK: loose power law S_B ~ w_B^{fit['alpha']:.3f} "
              f"(R^2 = {fit['r2']:.4f})")
    else:
        r2 = fit['r2'] if fit else float('nan')
        print(f"  REFUTED: no clean power law (R^2 = {r2:.4f})")

    summary = {
        'sweep': 'isoelectronic_He_like_n_max3',
        'Z_list': Z_list,
        'results': results,
        'power_law_S_vs_w': fit,
        'power_law_S_vs_Z': fit_Z,
        'power_law_S_vs_w_dimless': fit_dim,
    }

    out = os.path.join(os.path.dirname(__file__), 'data',
                       'ep2c_entropy_vs_vee_offdiag.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved -> {out}")


if __name__ == '__main__':
    main()
