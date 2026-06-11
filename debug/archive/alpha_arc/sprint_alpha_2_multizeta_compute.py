"""
Sprint α-Multi-zeta (Track α-2) — Fit physical Na 3s and 3p screened
wavefunctions to a multi-zeta Slater expansion.

Architecture preparation step for the α-PES test track. This script:
  1. Computes physical Na 3s and 3p screened wavefunctions via
     FrozenCore(Z=11) + _solve_screened_radial / _solve_screened_radial_log.
  2. Fits each to a K-zeta Slater STO expansion via least squares (allowing
     negative coefficients so radial nodes can form).
  3. Verifies normalization, orthogonality, mean radius, density at origin,
     and node count.
  4. Emits the fitted (c_k, ζ_k) parameters as JSON for the production-code
     entry to consume.

The fitted parameters will be added to ``geovac.multi_zeta_orbitals`` under a
Z=11 physical-fit registry, distinguishable from the heuristic two-zeta
infrastructure of the Xe-core branch.

Author: GeoVac Development Team (Track α-Multi-zeta, post-M-Y bimodule sprint)
Date: 2026-05-23
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import least_squares

from geovac.neon_core import (
    FrozenCore,
    _solve_screened_radial,
    _solve_screened_radial_log,
)


# ---------------------------------------------------------------------------
# Slater primitive (matches geovac.multi_zeta_orbitals.STO)
# ---------------------------------------------------------------------------

def sto_normalization(n: int, zeta: float) -> float:
    """N_{n, zeta} = (2 zeta)^{n + 1/2} / sqrt((2n)!)."""
    return (2.0 * zeta) ** (n + 0.5) / math.sqrt(float(math.factorial(2 * n)))


def sto_radial(n: int, zeta: float, r: np.ndarray) -> np.ndarray:
    """chi_n^{zeta}(r) = N * r^{n-1} * exp(-zeta r), STO convention.
    Integral chi^2 r^2 dr = 1.
    """
    N = sto_normalization(n, zeta)
    return N * r ** (n - 1) * np.exp(-zeta * r)


def multi_zeta_radial(
    primitives: List[Tuple[int, float]],
    coefficients: np.ndarray,
    r: np.ndarray,
) -> np.ndarray:
    """R(r) = sum_k c_k * chi_{n_k, zeta_k}(r)."""
    result = np.zeros_like(r, dtype=float)
    for (n, zeta), c in zip(primitives, coefficients):
        result = result + c * sto_radial(n, zeta, r)
    return result


# ---------------------------------------------------------------------------
# Physical Na 3s and 3p screened wavefunctions
# ---------------------------------------------------------------------------

def compute_physical_na_orbitals(
    n_grid: int = 16000, r_max: float = 60.0,
) -> Dict[str, Dict]:
    """Compute physical Na 3s and 3p screened wavefunctions on uniform grid.

    Returns dict with keys '3s', '3p' holding eigenvalue, R(r), r_grid,
    and key diagnostics (mean radius, density at origin, node count).
    """
    results: Dict[str, Dict] = {}

    # ----- Na 3s (l=0) via log-grid solver (handles l=0 correctly) -----
    energy_3s, u_3s, r_3s, R0_3s = _solve_screened_radial_log(
        Z=11, l=0, n_target=3,
        n_grid=n_grid, r_max=r_max,
    )
    # Convert reduced u(r) = r * R(r) to R(r)
    R_3s = u_3s / r_3s
    # Fix sign: choose R(0) > 0 convention for s-states
    if R_3s[0] < 0:
        R_3s = -R_3s
        u_3s = -u_3s

    # Renormalize so int |R|^2 r^2 dr = 1 (already from u, but check)
    norm_sq = np.trapezoid(R_3s ** 2 * r_3s ** 2, r_3s)
    R_3s = R_3s / math.sqrt(norm_sq)

    # Diagnostics
    rho_3s = R_3s ** 2 * r_3s ** 2
    mean_r_3s = float(np.trapezoid(r_3s * rho_3s, r_3s))
    # density at origin: extrapolate to r=0 via R(0)
    rho_origin_3s = float(R_3s[0] ** 2)  # R(0)^2 at smallest grid point
    # Count radial nodes by sign changes in R
    sign_changes_3s = int(np.sum(np.diff(np.sign(R_3s)) != 0))

    results['3s'] = {
        'l': 0,
        'n_phys': 3,
        'energy_Ha': float(energy_3s),
        'r_grid': r_3s.tolist(),
        'R_r': R_3s.tolist(),
        'mean_radius_bohr': mean_r_3s,
        'density_at_origin_au': rho_origin_3s,
        'sign_changes': sign_changes_3s,
        'normalization': float(np.trapezoid(R_3s ** 2 * r_3s ** 2, r_3s)),
    }

    # ----- Na 3p (l=1) via uniform-FD solver -----
    energy_3p, u_3p, r_3p = _solve_screened_radial(
        Z=11, l=1, n_target=3,
        n_grid=n_grid, r_max=r_max,
    )
    R_3p = u_3p / r_3p

    # Sign convention: choose first non-zero peak > 0
    # Find first index where |R| > 0.05 * max(|R|) and use that as sign reference
    idx_pk = int(np.argmax(np.abs(R_3p)))
    if R_3p[idx_pk] < 0:
        R_3p = -R_3p
        u_3p = -u_3p

    norm_sq_3p = np.trapezoid(R_3p ** 2 * r_3p ** 2, r_3p)
    R_3p = R_3p / math.sqrt(norm_sq_3p)

    rho_3p = R_3p ** 2 * r_3p ** 2
    mean_r_3p = float(np.trapezoid(r_3p * rho_3p, r_3p))
    # density at origin: 3p is l=1, R(r) ~ r^l at origin = r near 0
    rho_origin_3p = float(R_3p[0] ** 2)
    sign_changes_3p = int(np.sum(np.diff(np.sign(R_3p)) != 0))

    results['3p'] = {
        'l': 1,
        'n_phys': 3,
        'energy_Ha': float(energy_3p),
        'r_grid': r_3p.tolist(),
        'R_r': R_3p.tolist(),
        'mean_radius_bohr': mean_r_3p,
        'density_at_origin_au': rho_origin_3p,
        'sign_changes': sign_changes_3p,
        'normalization': float(np.trapezoid(R_3p ** 2 * r_3p ** 2, r_3p)),
    }

    return results


# ---------------------------------------------------------------------------
# Multi-zeta Slater fitting
# ---------------------------------------------------------------------------

def fit_multi_zeta(
    r: np.ndarray,
    R_target: np.ndarray,
    n_orbital: int,  # principal n of the orbital being fit (3 for Na 3s/3p)
    l_orbital: int,  # 0 for s, 1 for p
    K: int = 3,
    zeta_init: List[float] = None,
) -> Dict:
    """Fit R_target(r) on grid r to sum_k c_k * chi_{n_orbital, zeta_k}(r).

    The Slater n is fixed at n_orbital (uniform expansion); coefficients c_k
    and zetas zeta_k are optimized jointly via least squares with weights
    r^2 (so the fit minimizes integral |R_fit - R_target|^2 r^2 dr).

    Allows negative coefficients (essential for radial node formation).

    Parameters
    ----------
    r : ndarray
        Radial grid.
    R_target : ndarray
        Target wavefunction R(r) on the grid.
    n_orbital : int
        Slater n (3 for Na 3s/3p in this sprint).
    l_orbital : int
        Just used as a label; STO does not enforce l explicitly.
    K : int
        Number of zeta primitives.
    zeta_init : list of floats
        Initial zeta values (if None, distribute logarithmically).

    Returns
    -------
    dict with fitted params and quality-of-fit metrics.
    """
    if zeta_init is None:
        # Distribute logarithmically across a range appropriate for screened
        # valence: Na 3s mean r ~ 4 bohr => peak zeta ~ 1/4 = 0.25; allow
        # a broad inner zeta to capture core penetration.
        zeta_init = list(np.geomspace(0.4, 2.0, K))

    if len(zeta_init) != K:
        raise ValueError(f"zeta_init must have length {K}")

    # Parameter vector: [c_0, ..., c_{K-1}, zeta_0, ..., zeta_{K-1}]
    primitives_init = [(n_orbital, z) for z in zeta_init]
    # Initial coefficients: estimate as least-squares projection onto initial STOs
    A0 = np.column_stack([
        sto_radial(n_orbital, z, r) for z in zeta_init
    ])
    # Solve weighted least squares for initial c: minimize sum w_i (R_target - A c)^2
    weights = r  # use r weight (since R is amplitude and quadrature has r^2 dr)
    Aw = A0 * weights[:, None]
    bw = R_target * weights
    c0, _resid, _rank, _sv = np.linalg.lstsq(Aw, bw, rcond=None)

    x0 = np.concatenate([c0, np.array(zeta_init)])

    def residual(x):
        c = x[:K]
        z = x[K:]
        # Stabilize: enforce positive zetas
        if np.any(z <= 0.0):
            return np.full_like(R_target, 1e6)
        # Build R_fit
        R_fit = np.zeros_like(r)
        for k in range(K):
            R_fit = R_fit + c[k] * sto_radial(n_orbital, z[k], r)
        # Weight by r so the fit reflects integral |R - R_target|^2 r^2 dr
        return (R_fit - R_target) * r

    result = least_squares(
        residual, x0, method='lm', max_nfev=5000,
    )

    c_fit = result.x[:K]
    z_fit = result.x[K:]

    # Sort by zeta ascending (canonical form)
    order = np.argsort(z_fit)
    c_fit = c_fit[order]
    z_fit = z_fit[order]

    primitives = [(n_orbital, float(z)) for z in z_fit]
    coefficients = c_fit.astype(float)

    # Quality of fit
    R_fit = multi_zeta_radial(primitives, coefficients, r)
    norm_fit = float(np.trapezoid(R_fit ** 2 * r ** 2, r))
    # Renormalize fit to 1
    if norm_fit > 0:
        scale = 1.0 / math.sqrt(norm_fit)
        coefficients = coefficients * scale
        R_fit = R_fit * scale

    # Radial overlap with target (after renormalization both are unit-normalized)
    overlap = float(np.trapezoid(R_fit * R_target * r ** 2, r))
    # L2 error
    L2_err = float(np.trapezoid((R_fit - R_target) ** 2 * r ** 2, r))
    # Pointwise max abs deviation
    max_err = float(np.max(np.abs(R_fit - R_target)))

    return {
        'n_orbital': n_orbital,
        'l_orbital': l_orbital,
        'K': K,
        'zetas': z_fit.tolist(),
        'coefficients': coefficients.tolist(),
        'norm_pre_renorm': norm_fit,
        'overlap_with_target': overlap,
        'L2_error': L2_err,
        'max_pointwise_error': max_err,
        'optimizer_cost': float(result.cost),
        'optimizer_status': int(result.status),
    }


# ---------------------------------------------------------------------------
# Orthogonality verification
# ---------------------------------------------------------------------------

def verify_orthogonality(
    fit_3s: Dict,
    fit_3p: Dict,
    r: np.ndarray,
) -> Dict[str, float]:
    """Compute ⟨3s|3s⟩, ⟨3p|3p⟩, ⟨3s|3p⟩ on the radial integral with r^2 dr.

    Note: ⟨3s|3p⟩ should vanish by angular orthogonality regardless of radial
    overlap (different l). The "radial overlap" here is just a sanity check
    that the radial parts integrate properly; it does not need to be zero.
    """
    R_3s = multi_zeta_radial(
        [(fit_3s['n_orbital'], z) for z in fit_3s['zetas']],
        np.array(fit_3s['coefficients']),
        r,
    )
    R_3p = multi_zeta_radial(
        [(fit_3p['n_orbital'], z) for z in fit_3p['zetas']],
        np.array(fit_3p['coefficients']),
        r,
    )

    norm_3s = float(np.trapezoid(R_3s ** 2 * r ** 2, r))
    norm_3p = float(np.trapezoid(R_3p ** 2 * r ** 2, r))
    radial_overlap_3s_3p = float(np.trapezoid(R_3s * R_3p * r ** 2, r))

    # Mean radius
    mean_r_3s = float(np.trapezoid(r * R_3s ** 2 * r ** 2, r))
    mean_r_3p = float(np.trapezoid(r * R_3p ** 2 * r ** 2, r))

    # Sign changes for node count
    sc_3s = int(np.sum(np.diff(np.sign(R_3s)) != 0))
    sc_3p = int(np.sum(np.diff(np.sign(R_3p)) != 0))

    return {
        'norm_3s': norm_3s,
        'norm_3p': norm_3p,
        'radial_overlap_3s_3p': radial_overlap_3s_3p,
        'mean_radius_3s_bohr': mean_r_3s,
        'mean_radius_3p_bohr': mean_r_3p,
        'sign_changes_3s_fit': sc_3s,
        'sign_changes_3p_fit': sc_3p,
    }


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("Sprint alpha-Multi-zeta -- Fit physical Na 3s/3p to multi-zeta Slater")
    print("=" * 78)

    # --- §1. Compute physical Na 3s and 3p ---
    print("\n§1. Computing physical Na 3s and 3p screened wavefunctions...")
    phys = compute_physical_na_orbitals(n_grid=16000, r_max=60.0)

    print(f"\n  Na 3s:")
    print(f"    Eigenvalue         = {phys['3s']['energy_Ha']:.6f} Ha")
    print(f"    Mean radius <r>    = {phys['3s']['mean_radius_bohr']:.4f} bohr")
    print(f"    R(r_min)^2 (~|psi_origin|^2 norm) = {phys['3s']['density_at_origin_au']:.4e}")
    print(f"    Sign changes       = {phys['3s']['sign_changes']}")
    print(f"    Normalization      = {phys['3s']['normalization']:.6f}")

    print(f"\n  Na 3p:")
    print(f"    Eigenvalue         = {phys['3p']['energy_Ha']:.6f} Ha")
    print(f"    Mean radius <r>    = {phys['3p']['mean_radius_bohr']:.4f} bohr")
    print(f"    R(r_min)^2         = {phys['3p']['density_at_origin_au']:.4e}")
    print(f"    Sign changes       = {phys['3p']['sign_changes']}")
    print(f"    Normalization      = {phys['3p']['normalization']:.6f}")

    # --- §2. Multi-zeta Slater fits ---
    print("\n§2. Multi-zeta Slater fits (K=3 primitives per orbital)...")

    r_3s = np.array(phys['3s']['r_grid'])
    R_3s_target = np.array(phys['3s']['R_r'])

    fit_3s = fit_multi_zeta(
        r_3s, R_3s_target,
        n_orbital=3, l_orbital=0, K=3,
        # Initial zetas spanning core penetration (~5) to diffuse tail (~0.4)
        zeta_init=[5.0, 1.5, 0.5],
    )

    print(f"\n  Na 3s fit (K={fit_3s['K']}):")
    for i, (z, c) in enumerate(zip(fit_3s['zetas'], fit_3s['coefficients'])):
        print(f"    primitive {i}: zeta = {z:.6f}, c = {c:+.6f}")
    print(f"    overlap with target = {fit_3s['overlap_with_target']:.6f}")
    print(f"    L2 error            = {fit_3s['L2_error']:.4e}")
    print(f"    max pointwise error = {fit_3s['max_pointwise_error']:.4e}")

    r_3p = np.array(phys['3p']['r_grid'])
    R_3p_target = np.array(phys['3p']['R_r'])

    fit_3p = fit_multi_zeta(
        r_3p, R_3p_target,
        n_orbital=3, l_orbital=1, K=3,
        # Initial zetas for 3p: a bit smaller than 3s (3p more diffuse)
        zeta_init=[3.0, 1.0, 0.4],
    )

    print(f"\n  Na 3p fit (K={fit_3p['K']}):")
    for i, (z, c) in enumerate(zip(fit_3p['zetas'], fit_3p['coefficients'])):
        print(f"    primitive {i}: zeta = {z:.6f}, c = {c:+.6f}")
    print(f"    overlap with target = {fit_3p['overlap_with_target']:.6f}")
    print(f"    L2 error            = {fit_3p['L2_error']:.4e}")
    print(f"    max pointwise error = {fit_3p['max_pointwise_error']:.4e}")

    # Try K=4 for 3s for cross-check (3s has 2 nodes -> more flexibility helps)
    fit_3s_K4 = fit_multi_zeta(
        r_3s, R_3s_target,
        n_orbital=3, l_orbital=0, K=4,
        zeta_init=[8.0, 3.0, 1.0, 0.4],
    )
    print(f"\n  Na 3s fit (K={fit_3s_K4['K']} cross-check):")
    print(f"    overlap with target = {fit_3s_K4['overlap_with_target']:.6f}")
    print(f"    L2 error            = {fit_3s_K4['L2_error']:.4e}")
    print(f"    max pointwise error = {fit_3s_K4['max_pointwise_error']:.4e}")
    for i, (z, c) in enumerate(zip(fit_3s_K4['zetas'], fit_3s_K4['coefficients'])):
        print(f"    primitive {i}: zeta = {z:.6f}, c = {c:+.6f}")

    # --- §3. Orthogonality verification ---
    print("\n§3. Orthogonality verification...")
    # Use r_3s grid for both (they were solved on the same uniform grid here)
    # Actually both 3s and 3p were on different grids (log vs uniform); use 3s grid
    ortho = verify_orthogonality(fit_3s, fit_3p, r_3s)
    print(f"  <3s|3s>_radial             = {ortho['norm_3s']:.6f}")
    print(f"  <3p|3p>_radial             = {ortho['norm_3p']:.6f}")
    print(f"  <3s|3p>_radial (informational, l-orthog kills full) = {ortho['radial_overlap_3s_3p']:+.6e}")
    print(f"  Mean radius 3s (fit)       = {ortho['mean_radius_3s_bohr']:.4f} bohr")
    print(f"  Mean radius 3p (fit)       = {ortho['mean_radius_3p_bohr']:.4f} bohr")
    print(f"  Sign changes 3s (fit)      = {ortho['sign_changes_3s_fit']}")
    print(f"  Sign changes 3p (fit)      = {ortho['sign_changes_3p_fit']}")

    # --- Save results ---
    out_dir = Path(__file__).parent / 'data'
    out_dir.mkdir(parents=True, exist_ok=True)

    # Strip out the bulky r_grid / R_r arrays before saving the fit summary
    # (keep them only in a separate physical-orbital section so the JSON is reasonable size)
    save_data = {
        'sprint': 'alpha-2 Multi-zeta physical Na 3s/3p fit',
        'date': '2026-05-23',
        'physical_orbitals': {
            '3s': {
                'l': phys['3s']['l'],
                'n_phys': phys['3s']['n_phys'],
                'energy_Ha': phys['3s']['energy_Ha'],
                'mean_radius_bohr': phys['3s']['mean_radius_bohr'],
                'sign_changes': phys['3s']['sign_changes'],
                'normalization': phys['3s']['normalization'],
                # Save first 50 points of r_grid and R_r for record
                'r_grid_first50': phys['3s']['r_grid'][:50],
                'R_r_first50': phys['3s']['R_r'][:50],
            },
            '3p': {
                'l': phys['3p']['l'],
                'n_phys': phys['3p']['n_phys'],
                'energy_Ha': phys['3p']['energy_Ha'],
                'mean_radius_bohr': phys['3p']['mean_radius_bohr'],
                'sign_changes': phys['3p']['sign_changes'],
                'normalization': phys['3p']['normalization'],
                'r_grid_first50': phys['3p']['r_grid'][:50],
                'R_r_first50': phys['3p']['R_r'][:50],
            },
        },
        'fits': {
            '3s_K3': fit_3s,
            '3s_K4': fit_3s_K4,
            '3p_K3': fit_3p,
        },
        'orthogonality': ortho,
    }

    out_path = out_dir / 'sprint_alpha_2_multizeta_fits.json'
    with open(out_path, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\nSaved {out_path}")

    # --- §6. Verdict ---
    print("\n§6. Verdict:")
    ok_phys = (
        phys['3s']['normalization'] > 0.99 and phys['3s']['normalization'] < 1.01
        and phys['3p']['normalization'] > 0.99 and phys['3p']['normalization'] < 1.01
        and phys['3s']['energy_Ha'] < 0.0
        and phys['3p']['energy_Ha'] < 0.0
    )
    print(f"  Physical Na 3s/3p computed and bound          : {'Y' if ok_phys else 'N'}")
    ok_fit = (
        fit_3s['overlap_with_target'] > 0.99
        and fit_3p['overlap_with_target'] > 0.99
        and fit_3s['L2_error'] < 0.01
        and fit_3p['L2_error'] < 0.01
    )
    print(f"  Multi-zeta basis reproduces input to <1%       : {'Y' if ok_fit else 'N'}")
    nodes_ok = (
        ortho['sign_changes_3s_fit'] == 2 and ortho['sign_changes_3p_fit'] == 1
    )
    print(f"  Radial node count correct (3s=2, 3p=1)         : {'Y' if nodes_ok else 'N'}")

    if ok_phys and ok_fit and nodes_ok:
        print("\n  STATUS: All physical fits clean. Ready to wire into multi_zeta_orbitals.py.")
    else:
        print("\n  STATUS: At least one quality check failed; investigate before wiring.")

    return save_data


if __name__ == '__main__':
    main()
