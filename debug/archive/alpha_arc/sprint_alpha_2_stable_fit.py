"""Stable multi-zeta fit for Na 3s and 3p — bounded coefficients.

Adds bounds to avoid the K=5 near-degenerate-zeta instability we saw.
Constraint: |c_k| < 5 for all primitives (rules out the |c|~24 degenerate
solution). Forces optimizer to find a well-conditioned basis.
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


def fit_multi_zeta_bounded(
    r, R_target,
    n_list,
    zeta_init,
    zeta_bounds=(0.05, 30.0),  # broad zeta bounds
    coef_bound=5.0,  # |c_k| < coef_bound to avoid degenerate solutions
    max_nfev=20000,
):
    """Fit with bounded coefficients to enforce well-conditioned basis."""
    K = len(n_list)
    assert len(zeta_init) == K

    # Initial c from lstsq projection
    A0 = np.column_stack([sto_radial(n, z, r) for n, z in zip(n_list, zeta_init)])
    weights = r
    Aw = A0 * weights[:, None]
    bw = R_target * weights
    c0, *_ = np.linalg.lstsq(Aw, bw, rcond=None)
    # Clip initial c to bounds
    c0 = np.clip(c0, -coef_bound, coef_bound)

    x0 = np.concatenate([c0, np.array(zeta_init)])

    # Bounds: c in [-coef_bound, coef_bound], zeta in zeta_bounds
    lb = np.concatenate([
        np.full(K, -coef_bound),
        np.full(K, zeta_bounds[0]),
    ])
    ub = np.concatenate([
        np.full(K, coef_bound),
        np.full(K, zeta_bounds[1]),
    ])

    def residual(x):
        c = x[:K]
        z = x[K:]
        R_fit = np.zeros_like(r)
        for k in range(K):
            R_fit = R_fit + c[k] * sto_radial(n_list[k], z[k], r)
        return (R_fit - R_target) * r

    result = least_squares(
        residual, x0, bounds=(lb, ub),
        method='trf', max_nfev=max_nfev,
    )

    c_fit = result.x[:K]
    z_fit = result.x[K:]

    primitives = list(zip(n_list, z_fit.tolist()))
    R_fit = multi_zeta_radial(primitives, c_fit, r)

    # Renormalize
    norm_sq = np.trapezoid(R_fit ** 2 * r ** 2, r)
    if norm_sq > 0:
        scale = 1.0 / math.sqrt(norm_sq)
        c_fit = c_fit * scale
        R_fit = R_fit * scale

    overlap = np.trapezoid(R_fit * R_target * r ** 2, r)
    L2_err = np.trapezoid((R_fit - R_target) ** 2 * r ** 2, r)
    max_err = float(np.max(np.abs(R_fit - R_target)))
    sign_changes = int(np.sum(np.diff(np.sign(R_fit)) != 0))

    return {
        'n_list': n_list,
        'zetas': z_fit.tolist(),
        'coefficients': c_fit.tolist(),
        'overlap': float(overlap),
        'L2_error': float(L2_err),
        'max_pointwise_error': max_err,
        'sign_changes': sign_changes,
        'optimizer_cost': float(result.cost),
    }


def main():
    print("Stable multi-zeta fit for Na 3s and 3p (bounded |c| < 5)")
    print("=" * 60)

    phys = compute_physical_na_orbitals(n_grid=16000, r_max=60.0)

    # Na 3s: K=5 with well-spaced initial zetas
    r_3s = np.array(phys['3s']['r_grid'])
    R_3s = np.array(phys['3s']['R_r'])

    print("\n[Na 3s] K=5, n=[1,1,2,3,3], bounded coefs")
    fit_3s = fit_multi_zeta_bounded(
        r_3s, R_3s,
        n_list=[1, 1, 2, 3, 3],
        zeta_init=[15.0, 6.0, 3.0, 1.0, 0.5],
        coef_bound=5.0,
    )
    print(f"  primitives (n, zeta, c):")
    for n, z, c in zip(fit_3s['n_list'], fit_3s['zetas'], fit_3s['coefficients']):
        print(f"    n={n}, zeta={z:.6f}, c={c:+.6f}")
    print(f"  overlap = {fit_3s['overlap']:.6f}")
    print(f"  L2 err = {fit_3s['L2_error']:.4e}")
    print(f"  max err = {fit_3s['max_pointwise_error']:.4e}")
    print(f"  nodes = {fit_3s['sign_changes']}")

    # Try K=4 also (simpler, cleaner if it works)
    print("\n[Na 3s] K=4, n=[1,2,3,3], bounded coefs")
    fit_3s_K4 = fit_multi_zeta_bounded(
        r_3s, R_3s,
        n_list=[1, 2, 3, 3],
        zeta_init=[10.0, 3.0, 1.0, 0.5],
        coef_bound=5.0,
    )
    print(f"  primitives:")
    for n, z, c in zip(fit_3s_K4['n_list'], fit_3s_K4['zetas'], fit_3s_K4['coefficients']):
        print(f"    n={n}, zeta={z:.6f}, c={c:+.6f}")
    print(f"  overlap = {fit_3s_K4['overlap']:.6f}")
    print(f"  L2 err = {fit_3s_K4['L2_error']:.4e}")
    print(f"  max err = {fit_3s_K4['max_pointwise_error']:.4e}")
    print(f"  nodes = {fit_3s_K4['sign_changes']}")

    # Na 3p
    r_3p = np.array(phys['3p']['r_grid'])
    R_3p = np.array(phys['3p']['R_r'])

    print("\n[Na 3p] K=4, n=[2,2,3,3], bounded coefs")
    fit_3p = fit_multi_zeta_bounded(
        r_3p, R_3p,
        n_list=[2, 2, 3, 3],
        zeta_init=[6.0, 3.0, 1.0, 0.4],
        coef_bound=5.0,
    )
    print(f"  primitives:")
    for n, z, c in zip(fit_3p['n_list'], fit_3p['zetas'], fit_3p['coefficients']):
        print(f"    n={n}, zeta={z:.6f}, c={c:+.6f}")
    print(f"  overlap = {fit_3p['overlap']:.6f}")
    print(f"  L2 err = {fit_3p['L2_error']:.4e}")
    print(f"  max err = {fit_3p['max_pointwise_error']:.4e}")
    print(f"  nodes = {fit_3p['sign_changes']}")

    # Pick the best 3s fit: K=4 if it gives nodes=2 and overlap > 0.99, otherwise K=5
    print("\n[Selection] choose best 3s fit:")
    if fit_3s_K4['sign_changes'] == 2 and fit_3s_K4['overlap'] > 0.99:
        chosen_3s = fit_3s_K4
        choice_3s = 'K=4'
    elif fit_3s['sign_changes'] == 2 and fit_3s['overlap'] > 0.99:
        chosen_3s = fit_3s
        choice_3s = 'K=5'
    else:
        # Pick the higher-overlap one
        if fit_3s['overlap'] > fit_3s_K4['overlap']:
            chosen_3s = fit_3s
            choice_3s = 'K=5 (fallback)'
        else:
            chosen_3s = fit_3s_K4
            choice_3s = 'K=4 (fallback)'
    print(f"  Chosen: {choice_3s}, overlap={chosen_3s['overlap']:.6f}, nodes={chosen_3s['sign_changes']}")

    # Verify physics
    R_3s_fit = multi_zeta_radial(
        list(zip(chosen_3s['n_list'], chosen_3s['zetas'])),
        np.array(chosen_3s['coefficients']),
        r_3s,
    )
    mean_r_3s_fit = float(np.trapezoid(r_3s * R_3s_fit ** 2 * r_3s ** 2, r_3s))
    print(f"  mean radius (fit) = {mean_r_3s_fit:.4f} bohr (physical {phys['3s']['mean_radius_bohr']:.4f})")

    # Orthogonality on common grid
    r_common = np.geomspace(1e-4, 60.0, 8000)
    R_3s_common = multi_zeta_radial(
        list(zip(chosen_3s['n_list'], chosen_3s['zetas'])),
        np.array(chosen_3s['coefficients']),
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

    print(f"\nOrthogonality on common grid:")
    print(f"  <3s|3s>_radial = {norm_3s:.6f}")
    print(f"  <3p|3p>_radial = {norm_3p:.6f}")
    print(f"  <3s|3p>_radial = {radial_overlap:+.6e} (informational)")

    # Final JSON
    out_dir = Path(__file__).parent / 'data'
    out_dir.mkdir(parents=True, exist_ok=True)

    save_data = {
        'sprint': 'alpha-2 final stable multi-zeta fit for Na 3s and 3p',
        'date': '2026-05-23',
        'description': (
            'Production parameters with bounded coefficients (|c| <= 5) to '
            'prevent near-degenerate basis instability. Selected based on '
            'overlap > 0.99 and correct node count.'
        ),
        'physical_observables': {
            'Na_3s_eigenvalue_Ha': phys['3s']['energy_Ha'],
            'Na_3p_eigenvalue_Ha': phys['3p']['energy_Ha'],
            'Na_3s_mean_radius_bohr_physical': phys['3s']['mean_radius_bohr'],
            'Na_3p_mean_radius_bohr_physical': phys['3p']['mean_radius_bohr'],
            'Na_3s_radial_nodes_physical': 2,
            'Na_3p_radial_nodes_physical': 1,
        },
        'production_fits': {
            'Na_3s': {
                'n_orbital': 3,
                'l_orbital': 0,
                'K': len(chosen_3s['n_list']),
                'n_slater_list': chosen_3s['n_list'],
                'zetas': chosen_3s['zetas'],
                'coefficients': chosen_3s['coefficients'],
                'overlap_with_physical': chosen_3s['overlap'],
                'L2_error': chosen_3s['L2_error'],
                'max_pointwise_error': chosen_3s['max_pointwise_error'],
                'radial_nodes_fit': chosen_3s['sign_changes'],
                'mean_radius_bohr_fit': mean_r_3s_fit,
                'selection_note': choice_3s,
            },
            'Na_3p': {
                'n_orbital': 3,
                'l_orbital': 1,
                'K': len(fit_3p['n_list']),
                'n_slater_list': fit_3p['n_list'],
                'zetas': fit_3p['zetas'],
                'coefficients': fit_3p['coefficients'],
                'overlap_with_physical': fit_3p['overlap'],
                'L2_error': fit_3p['L2_error'],
                'max_pointwise_error': fit_3p['max_pointwise_error'],
                'radial_nodes_fit': fit_3p['sign_changes'],
            },
        },
        'orthonormality': {
            'norm_3s': norm_3s,
            'norm_3p': norm_3p,
            'radial_overlap_3s_3p': radial_overlap,
        },
        'alternative_fits': {
            'Na_3s_K5': {
                'n_slater_list': fit_3s['n_list'],
                'zetas': fit_3s['zetas'],
                'coefficients': fit_3s['coefficients'],
                'overlap': fit_3s['overlap'],
                'L2_error': fit_3s['L2_error'],
                'radial_nodes_fit': fit_3s['sign_changes'],
            },
            'Na_3s_K4': {
                'n_slater_list': fit_3s_K4['n_list'],
                'zetas': fit_3s_K4['zetas'],
                'coefficients': fit_3s_K4['coefficients'],
                'overlap': fit_3s_K4['overlap'],
                'L2_error': fit_3s_K4['L2_error'],
                'radial_nodes_fit': fit_3s_K4['sign_changes'],
            },
        },
    }

    out_path = out_dir / 'sprint_alpha_2_multizeta_fits.json'
    with open(out_path, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\nSaved final production fits to {out_path}")


if __name__ == '__main__':
    main()
