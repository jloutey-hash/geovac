"""Diagnose the 3s node count and near-origin behavior of the multi-zeta fit.

The 3s physical wavefunction has 2 radial nodes. With the Slater basis at
fixed Slater n=3, we have r^{n-1} = r^2 prefactor, so the basis functions
look like r^2 * exp(-zeta r). A linear combination with two sign-change
zeros requires enough flexibility — K=3 may be insufficient.

We test:
 - K=4 with better initial zetas (and increase the importance of inner region)
 - K=5
 - Whether using mixed Slater n (n=1, 2, 3) helps capture node structure
 - The near-origin behavior: physical Na 3s has very large value at r->0 from
   the inner electron tail; fit may miss this without an inner-shell-class
   primitive.
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


def fit_multi_zeta_mixed_n(
    r, R_target,
    n_list,  # list of Slater n's, length K
    zeta_init,  # list of initial zetas, length K
    l_label=0,
):
    """Fit with mixed Slater n values."""
    K = len(n_list)
    assert len(zeta_init) == K

    primitives_init = list(zip(n_list, zeta_init))

    # Initial coefficients via least squares against initial primitives
    A0 = np.column_stack([sto_radial(n, z, r) for n, z in primitives_init])
    weights = r
    Aw = A0 * weights[:, None]
    bw = R_target * weights
    c0, *_ = np.linalg.lstsq(Aw, bw, rcond=None)

    x0 = np.concatenate([c0, np.array(zeta_init)])

    def residual(x):
        c = x[:K]
        z = x[K:]
        if np.any(z <= 0.0):
            return np.full_like(R_target, 1e6)
        R_fit = np.zeros_like(r)
        for k in range(K):
            R_fit = R_fit + c[k] * sto_radial(n_list[k], z[k], r)
        return (R_fit - R_target) * r

    result = least_squares(residual, x0, method='lm', max_nfev=10000)

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
    }


def main():
    print("Diagnose 3s node count and near-origin behavior")
    print("=" * 60)

    phys = compute_physical_na_orbitals(n_grid=16000, r_max=60.0)
    r = np.array(phys['3s']['r_grid'])
    R_target = np.array(phys['3s']['R_r'])

    # Look at the wavefunction structure
    print(f"\nPhysical Na 3s on grid r=[{r[0]:.4e}, {r[-1]:.2f}]:")
    print(f"  R(r=0+) = {R_target[0]:.6e}")
    print(f"  R range  = [{R_target.min():.4e}, {R_target.max():.4e}]")
    # Find sign changes and zero crossings
    signs = np.sign(R_target)
    nonzero_mask = signs != 0
    sign_pts = []
    for i in range(1, len(R_target)):
        if signs[i] != signs[i-1] and signs[i] != 0 and signs[i-1] != 0:
            sign_pts.append((r[i-1], r[i]))
    print(f"  Sign changes at intervals:")
    for r0, r1 in sign_pts:
        print(f"    r in [{r0:.4f}, {r1:.4f}]")
    print(f"  Number of radial nodes  = {len(sign_pts)}")

    # Show profile at various radii
    print(f"\nProfile R(r):")
    for r_target in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0, 15.0]:
        idx = np.argmin(np.abs(r - r_target))
        print(f"  r = {r[idx]:6.3f}  R = {R_target[idx]:+.4e}")

    # ----- K=4 mixed n: include n=1 (1s-class) for inner shell, n=2 (2s/2p-class),
    # and n=3 (3s/3p-class) primitives -----
    print("\nFit 3s with mixed (n_Slater) primitives:")

    # Test 1: K=4 mixed, n in (1, 2, 3, 3) with broad zetas
    print("\n[Test 1] K=4, n_list=[1, 2, 3, 3], inner-to-outer zetas")
    fit1 = fit_multi_zeta_mixed_n(
        r, R_target,
        n_list=[1, 2, 3, 3],
        zeta_init=[8.0, 4.0, 1.5, 0.5],
    )
    print(f"  zetas: {fit1['zetas']}")
    print(f"  coefs: {fit1['coefficients']}")
    print(f"  overlap = {fit1['overlap']:.6f}")
    print(f"  L2 err  = {fit1['L2_error']:.4e}")
    print(f"  max err = {fit1['max_pointwise_error']:.4e}")
    print(f"  nodes   = {fit1['sign_changes']}")

    # Test 2: K=5 mixed
    print("\n[Test 2] K=5, n_list=[1, 1, 2, 3, 3]")
    fit2 = fit_multi_zeta_mixed_n(
        r, R_target,
        n_list=[1, 1, 2, 3, 3],
        zeta_init=[10.0, 5.0, 3.0, 1.5, 0.5],
    )
    print(f"  zetas: {fit2['zetas']}")
    print(f"  coefs: {fit2['coefficients']}")
    print(f"  overlap = {fit2['overlap']:.6f}")
    print(f"  L2 err  = {fit2['L2_error']:.4e}")
    print(f"  max err = {fit2['max_pointwise_error']:.4e}")
    print(f"  nodes   = {fit2['sign_changes']}")

    # Test 3: K=6 mixed for fullness
    print("\n[Test 3] K=6, n_list=[1, 1, 2, 2, 3, 3]")
    fit3 = fit_multi_zeta_mixed_n(
        r, R_target,
        n_list=[1, 1, 2, 2, 3, 3],
        zeta_init=[12.0, 6.0, 4.0, 2.0, 1.0, 0.4],
    )
    print(f"  zetas: {fit3['zetas']}")
    print(f"  coefs: {fit3['coefficients']}")
    print(f"  overlap = {fit3['overlap']:.6f}")
    print(f"  L2 err  = {fit3['L2_error']:.4e}")
    print(f"  max err = {fit3['max_pointwise_error']:.4e}")
    print(f"  nodes   = {fit3['sign_changes']}")

    # Now Na 3p with mixed n
    print("\n" + "=" * 60)
    print("Na 3p with mixed-n fit")
    r_3p = np.array(phys['3p']['r_grid'])
    R_3p = np.array(phys['3p']['R_r'])

    print("\n[Test 4] K=4 3p, n_list=[2, 2, 3, 3]")
    fit4 = fit_multi_zeta_mixed_n(
        r_3p, R_3p,
        n_list=[2, 2, 3, 3],
        zeta_init=[6.0, 3.0, 1.0, 0.4],
    )
    print(f"  zetas: {fit4['zetas']}")
    print(f"  coefs: {fit4['coefficients']}")
    print(f"  overlap = {fit4['overlap']:.6f}")
    print(f"  L2 err  = {fit4['L2_error']:.4e}")
    print(f"  max err = {fit4['max_pointwise_error']:.4e}")
    print(f"  nodes   = {fit4['sign_changes']}")

    # Save best fits as the "Z=11 physical-fit" production parameters
    out_dir = Path(__file__).parent / 'data'
    save_data = {
        'sprint': 'alpha-2 mixed-n diagnose for 3s nodes',
        'phys_3s_nodes': len(sign_pts),
        'fits_3s': {
            'K4_n1233': fit1,
            'K5_n11233': fit2,
            'K6_n112233': fit3,
        },
        'fits_3p': {
            'K4_n2233': fit4,
        },
    }
    out_path = out_dir / 'sprint_alpha_2_mixed_n_diagnose.json'
    with open(out_path, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\nSaved {out_path}")


if __name__ == '__main__':
    main()
