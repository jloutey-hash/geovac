"""Final multi-zeta fit for Na 3s and 3p with mixed-n primitives.

Produces the production parameters for Z=11 physical-fit entry in
geovac/multi_zeta_orbitals.py.

Uses mixed Slater n primitives (n in {1, 2, 3}) so the basis can
represent the inner-shell tail of Na 3s (which has 2 radial nodes
and a large amplitude at the origin from core penetration).

Selected fit:
  3s: K=5, n_list=[1,1,2,3,3]
  3p: K=4, n_list=[2,2,3,3]

Both achieve overlap with target > 0.999999, L2 error < 1e-5, and
the correct number of radial nodes (3s=2, 3p=1).
"""

import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

import sys
sys.path.insert(0, str(Path(__file__).parent))
from sprint_alpha_2_multizeta_compute import (
    compute_physical_na_orbitals,
    sto_radial,
    multi_zeta_radial,
)
from sprint_alpha_2_node_diagnose import fit_multi_zeta_mixed_n


def main():
    print("Final multi-zeta fit for Na 3s and 3p")
    print("=" * 60)

    phys = compute_physical_na_orbitals(n_grid=16000, r_max=60.0)

    # ---- Na 3s: K=5, n=[1,1,2,3,3] ----
    r_3s = np.array(phys['3s']['r_grid'])
    R_3s = np.array(phys['3s']['R_r'])

    fit_3s = fit_multi_zeta_mixed_n(
        r_3s, R_3s,
        n_list=[1, 1, 2, 3, 3],
        zeta_init=[10.0, 5.0, 3.0, 1.0, 0.5],
        l_label=0,
    )

    print("\nNa 3s fit (K=5, n=[1,1,2,3,3]):")
    for i, (n, z, c) in enumerate(zip([1, 1, 2, 3, 3], fit_3s['zetas'], fit_3s['coefficients'])):
        print(f"  primitive {i}: n={n}, zeta={z:.6f}, c={c:+.6f}")
    print(f"  overlap = {fit_3s['overlap']:.8f}")
    print(f"  L2 error = {fit_3s['L2_error']:.4e}")
    print(f"  nodes = {fit_3s['sign_changes']}")

    # Verify physics: compute mean radius from the fit
    R_3s_fit = multi_zeta_radial(
        list(zip(fit_3s['n_list'], fit_3s['zetas'])),
        np.array(fit_3s['coefficients']),
        r_3s,
    )
    mean_r_3s_fit = float(np.trapezoid(r_3s * R_3s_fit ** 2 * r_3s ** 2, r_3s))
    print(f"  mean radius (fit) = {mean_r_3s_fit:.4f} bohr (physical {phys['3s']['mean_radius_bohr']:.4f})")
    print(f"  R(0) approx = {R_3s_fit[0]:.4e}")

    # ---- Na 3p: K=4, n=[2,2,3,3] ----
    r_3p = np.array(phys['3p']['r_grid'])
    R_3p = np.array(phys['3p']['R_r'])

    fit_3p = fit_multi_zeta_mixed_n(
        r_3p, R_3p,
        n_list=[2, 2, 3, 3],
        zeta_init=[6.0, 3.0, 1.0, 0.4],
        l_label=1,
    )

    print("\nNa 3p fit (K=4, n=[2,2,3,3]):")
    for i, (n, z, c) in enumerate(zip([2, 2, 3, 3], fit_3p['zetas'], fit_3p['coefficients'])):
        print(f"  primitive {i}: n={n}, zeta={z:.6f}, c={c:+.6f}")
    print(f"  overlap = {fit_3p['overlap']:.8f}")
    print(f"  L2 error = {fit_3p['L2_error']:.4e}")
    print(f"  nodes = {fit_3p['sign_changes']}")

    R_3p_fit = multi_zeta_radial(
        list(zip(fit_3p['n_list'], fit_3p['zetas'])),
        np.array(fit_3p['coefficients']),
        r_3p,
    )
    mean_r_3p_fit = float(np.trapezoid(r_3p * R_3p_fit ** 2 * r_3p ** 2, r_3p))
    print(f"  mean radius (fit) = {mean_r_3p_fit:.4f} bohr (physical {phys['3p']['mean_radius_bohr']:.4f})")

    # ---- Orthogonality (informational) ----
    # Use a common dense grid; both wavefunctions have to be evaluated on it
    r_common = np.geomspace(1e-4, 60.0, 8000)
    R_3s_common = multi_zeta_radial(
        list(zip(fit_3s['n_list'], fit_3s['zetas'])),
        np.array(fit_3s['coefficients']),
        r_common,
    )
    R_3p_common = multi_zeta_radial(
        list(zip(fit_3p['n_list'], fit_3p['zetas'])),
        np.array(fit_3p['coefficients']),
        r_common,
    )
    norm_3s = float(np.trapezoid(R_3s_common ** 2 * r_common ** 2, r_common))
    norm_3p = float(np.trapezoid(R_3p_common ** 2 * r_common ** 2, r_common))
    radial_overlap = float(np.trapezoid(R_3s_common * R_3p_common * r_common ** 2, r_common))

    print(f"\nOrthogonality:")
    print(f"  <3s|3s>_radial = {norm_3s:.6f}")
    print(f"  <3p|3p>_radial = {norm_3p:.6f}")
    print(f"  <3s|3p>_radial = {radial_overlap:+.6e}  (informational; full <3s|3p>=0 by l-orthog)")

    # ---- Save final production parameters ----
    out_dir = Path(__file__).parent / 'data'
    out_dir.mkdir(parents=True, exist_ok=True)

    save_data = {
        'sprint': 'alpha-2 final mixed-n multi-zeta fit for Na 3s and 3p',
        'date': '2026-05-23',
        'description': (
            'Production parameters for Z=11 physical-fit entry in '
            'geovac/multi_zeta_orbitals.py. Mixed Slater n primitives with '
            'optimized zetas, fitted to physical Na 3s and 3p screened '
            'wavefunctions from FrozenCore([Ne]) + radial Schrodinger solver.'
        ),
        'physical_observables': {
            'Na_3s_eigenvalue_Ha': phys['3s']['energy_Ha'],
            'Na_3p_eigenvalue_Ha': phys['3p']['energy_Ha'],
            'Na_3s_mean_radius_bohr_physical': phys['3s']['mean_radius_bohr'],
            'Na_3p_mean_radius_bohr_physical': phys['3p']['mean_radius_bohr'],
            'Na_3s_mean_radius_bohr_fit': mean_r_3s_fit,
            'Na_3p_mean_radius_bohr_fit': mean_r_3p_fit,
            'Na_3s_radial_nodes': 2,
            'Na_3p_radial_nodes': 1,
        },
        'production_fits': {
            'Na_3s': {
                'n_orbital': 3,
                'l_orbital': 0,
                'K': 5,
                'n_slater_list': [1, 1, 2, 3, 3],
                'zetas': fit_3s['zetas'],
                'coefficients': fit_3s['coefficients'],
                'overlap_with_physical': fit_3s['overlap'],
                'L2_error': fit_3s['L2_error'],
                'max_pointwise_error': fit_3s['max_pointwise_error'],
                'radial_nodes': fit_3s['sign_changes'],
            },
            'Na_3p': {
                'n_orbital': 3,
                'l_orbital': 1,
                'K': 4,
                'n_slater_list': [2, 2, 3, 3],
                'zetas': fit_3p['zetas'],
                'coefficients': fit_3p['coefficients'],
                'overlap_with_physical': fit_3p['overlap'],
                'L2_error': fit_3p['L2_error'],
                'max_pointwise_error': fit_3p['max_pointwise_error'],
                'radial_nodes': fit_3p['sign_changes'],
            },
        },
        'orthonormality': {
            'norm_3s': norm_3s,
            'norm_3p': norm_3p,
            'radial_overlap_3s_3p': radial_overlap,
        },
    }

    out_path = out_dir / 'sprint_alpha_2_multizeta_fits.json'
    with open(out_path, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\nSaved final production fits to {out_path}")


if __name__ == '__main__':
    main()
