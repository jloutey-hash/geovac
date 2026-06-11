"""
Task 10: Self-Consistent Projected PK Iteration.

Tests whether iteratively updating PK weights based on l=0 content of the
eigenstate stabilizes R_eq for LiH at higher l_max.

Also tests R-dependent PK scaling as an alternative (Part 4).

Method: Wrapper around the existing angular eigenvalue solver. At each
iteration:
  1. Build angular Hamiltonian with NO PK
  2. Add PK manually with per-channel weights w[l1,l2] = (w1, w2)
  3. Solve -> eigenvector psi_0
  4. Extract per-electron l=0 content (p1, p2) from psi_0
  5. Update weights: w_new[l1,l2] = (delta_{l1,0}*p1, delta_{l2,0}*p2)
  6. Damp: w = lambda*w_new + (1-lambda)*w_old
  7. Repeat until converged

References:
  Paper 17 Sec III.B (PK pseudopotential), Sec VI.A (l_max divergence)
  Paper 15 Sec IV (angular eigenvalue problem)
  pk_saturation_real.json (Task 9 results: +0.23 bohr/l_max drift)
"""

import sys
import json
import time
import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from geovac.level4_multichannel import (
    _channel_list,
    build_angular_hamiltonian,
    compute_pk_pseudopotential,
)

# === Physical parameters (match pk_saturation_real.py / Task 9) ===
Z_A_EFF = 1.0         # Screened charge for valence (Z_A - n_core = 3-2)
Z_B = 1.0             # H
Z_A_BARE = 3.0        # Li bare nuclear charge
PK_A = 6.93           # Ha*bohr^2 (ab initio PK amplitude)
PK_B_PARAM = 6.80     # bohr^-2 (ab initio PK Gaussian exponent)
R_EXP = 3.015         # Experimental R_eq (bohr)

# Iteration parameters
MAX_ITER = 20
DAMPING = 0.5
TOL = 1e-4

# Grid parameters
N_ALPHA = 100
N_RE = 25
R_E_GRID = np.logspace(np.log10(0.3), np.log10(5.0), N_RE)

# PES scan grid (same as Task 9)
R_VALUES = np.array([2.0, 2.4, 2.7, 3.0, 3.2, 3.5, 3.8, 4.2, 4.7, 5.5, 6.5])


def get_channels(l_max: int) -> list:
    """Get channel list. Z_A_eff=Z_B=1 → homonuclear symmetry."""
    return _channel_list(l_max, homonuclear=(Z_A_EFF == Z_B))


def solve_angular_selfconsistent(
    R: float,
    R_e: float,
    l_max: int,
    n_alpha: int = N_ALPHA,
    max_iter: int = MAX_ITER,
    damping: float = DAMPING,
    tol: float = TOL,
    verbose: bool = False,
) -> dict:
    """
    Solve angular eigenvalue problem with self-consistent PK weights.

    Returns dict with converged mu0, p1, p2, iteration history.
    """
    rho = R / (2.0 * R_e)
    z0 = R * (Z_A_EFF - Z_B) / (2.0 * (Z_A_EFF + Z_B))

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    channels = get_channels(l_max)
    n_ch = len(channels)

    # Compute PK potentials (V_pk_e1, V_pk_e2) once
    rho_A = (R / 2.0 - z0) / R_e
    rho_B = (R / 2.0 + z0) / R_e
    pk_pots_raw = [{
        'C_core': PK_A,
        'beta_core': PK_B_PARAM,
        'atom': 'A',
    }]
    V_pk_e1, V_pk_e2 = compute_pk_pseudopotential(
        alpha, rho_A, R_e, pk_pots_raw, rho_B=rho_B,
    )

    # Build base Hamiltonian WITHOUT PK
    H_base = build_angular_hamiltonian(
        alpha, rho, R_e, l_max=l_max, Z=1.0,
        m_max=0, Z_A=Z_A_EFF, Z_B=Z_B, z0=z0,
        pk_potentials=None,
    )

    # Initialize PK weights: l-dependent (delta_{l_i,0})
    # w_pk[ic] = (w1, w2) for channel ic
    w_pk = np.zeros((n_ch, 2))
    for ic, (l1, l2) in enumerate(channels):
        w_pk[ic, 0] = 1.0 if l1 == 0 else 0.0
        w_pk[ic, 1] = 1.0 if l2 == 0 else 0.0

    history = []

    for iteration in range(max_iter):
        # Build H = H_base + PK with current weights
        H = H_base.copy()
        for ic in range(n_ch):
            w1, w2 = w_pk[ic]
            if w1 == 0.0 and w2 == 0.0:
                continue
            for i in range(n_alpha):
                ii = ic * n_alpha + i
                H[ii, ii] += R_e * w1 * V_pk_e1[i]
                H[ii, ii] += R_e * w2 * V_pk_e2[i]

        # Solve
        evals, evecs = eigh(H)
        mu0 = evals[0]
        psi = evecs[:, 0]

        # Extract per-channel norms
        channel_norms = np.zeros(n_ch)
        for ic in range(n_ch):
            start = ic * n_alpha
            end = start + n_alpha
            channel_norms[ic] = np.sum(psi[start:end] ** 2) * h

        total_norm = np.sum(channel_norms)
        if total_norm < 1e-30:
            total_norm = 1.0

        # Per-electron l=0 content
        p1 = 0.0  # fraction with electron 1 in l=0
        p2 = 0.0  # fraction with electron 2 in l=0
        for ic, (l1, l2) in enumerate(channels):
            if l1 == 0:
                p1 += channel_norms[ic]
            if l2 == 0:
                p2 += channel_norms[ic]
        p1 /= total_norm
        p2 /= total_norm

        # Compute new weights
        w_new = np.zeros((n_ch, 2))
        for ic, (l1, l2) in enumerate(channels):
            w_new[ic, 0] = p1 if l1 == 0 else 0.0
            w_new[ic, 1] = p2 if l2 == 0 else 0.0

        # Check convergence
        delta = np.max(np.abs(w_new - w_pk))

        history.append({
            'iteration': iteration,
            'mu0': float(mu0),
            'p1': float(p1),
            'p2': float(p2),
            'delta': float(delta),
        })

        if verbose:
            print(f"    iter {iteration:2d}: mu0={mu0:.6f}, "
                  f"p1={p1:.4f}, p2={p2:.4f}, delta={delta:.6f}")

        # Damp and update
        w_pk = damping * w_new + (1.0 - damping) * w_pk

        if delta < tol:
            if verbose:
                print(f"    Converged at iteration {iteration}")
            break

    converged = (delta < tol) if history else False

    return {
        'mu0': float(mu0),
        'p1': float(p1),
        'p2': float(p2),
        'converged': converged,
        'n_iter': len(history),
        'history': history,
        'R': R,
        'R_e': R_e,
        'l_max': l_max,
    }


def solve_angular_ldependent(
    R: float, R_e: float, l_max: int, n_alpha: int = N_ALPHA,
) -> float:
    """Standard l-dependent PK (w=1.0 for l=0). Returns mu0."""
    rho = R / (2.0 * R_e)
    z0 = R * (Z_A_EFF - Z_B) / (2.0 * (Z_A_EFF + Z_B))

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    pk_pots = [{
        'C_core': PK_A,
        'beta_core': PK_B_PARAM,
        'atom': 'A',
        'channel_mode': 'l_dependent',
    }]

    H = build_angular_hamiltonian(
        alpha, rho, R_e, l_max=l_max, Z=1.0,
        m_max=0, Z_A=Z_A_EFF, Z_B=Z_B, z0=z0,
        pk_potentials=pk_pots,
    )
    evals = eigh(H, eigvals_only=True)
    return float(evals[0])


def solve_angular_r_dependent_pk(
    R: float, R_e: float, l_max: int, R_ref: float = R_EXP,
    pk_cap: float = 2.0, n_alpha: int = N_ALPHA,
) -> float:
    """
    R-dependent PK scaling: strengthen PK at large R.

    w_PK(R) = delta_{l_i,0} * min(pk_cap, R/R_ref) for R > R_ref
    w_PK(R) = delta_{l_i,0} for R <= R_ref

    Returns mu0.
    """
    rho = R / (2.0 * R_e)
    z0 = R * (Z_A_EFF - Z_B) / (2.0 * (Z_A_EFF + Z_B))

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    channels = get_channels(l_max)
    n_ch = len(channels)

    # PK potentials
    rho_A = (R / 2.0 - z0) / R_e
    rho_B = (R / 2.0 + z0) / R_e
    pk_pots_raw = [{
        'C_core': PK_A,
        'beta_core': PK_B_PARAM,
        'atom': 'A',
    }]
    V_pk_e1, V_pk_e2 = compute_pk_pseudopotential(
        alpha, rho_A, R_e, pk_pots_raw, rho_B=rho_B,
    )

    # Build base H without PK
    H = build_angular_hamiltonian(
        alpha, rho, R_e, l_max=l_max, Z=1.0,
        m_max=0, Z_A=Z_A_EFF, Z_B=Z_B, z0=z0,
        pk_potentials=None,
    )

    # R-dependent scaling factor
    scale = min(pk_cap, max(1.0, R / R_ref))

    # Add PK with scaled l-dependent weights
    for ic, (l1, l2) in enumerate(channels):
        w1 = scale if l1 == 0 else 0.0
        w2 = scale if l2 == 0 else 0.0
        if w1 == 0.0 and w2 == 0.0:
            continue
        for i in range(n_alpha):
            ii = ic * n_alpha + i
            H[ii, ii] += R_e * w1 * V_pk_e1[i]
            H[ii, ii] += R_e * w2 * V_pk_e2[i]

    evals = eigh(H, eigvals_only=True)
    return float(evals[0])


def compute_pes(
    l_max: int,
    R_grid: np.ndarray,
    Re_grid: np.ndarray,
    mode: str = 'l_dependent',
    R_ref: float = R_EXP,
    pk_cap: float = 2.0,
    max_iter: int = MAX_ITER,
    verbose: bool = False,
) -> dict:
    """
    Compute PES curve for a given l_max and PK mode.

    mode: 'l_dependent', 'selfconsistent', 'r_dependent'
    """
    n_R = len(R_grid)
    E_elec = np.zeros(n_R)

    for ir, R in enumerate(R_grid):
        U_vals = np.zeros(len(Re_grid))

        for ire, R_e in enumerate(Re_grid):
            if mode == 'l_dependent':
                mu0 = solve_angular_ldependent(R, R_e, l_max)
            elif mode == 'selfconsistent':
                res = solve_angular_selfconsistent(
                    R, R_e, l_max, max_iter=max_iter, verbose=False,
                )
                mu0 = res['mu0']
            elif mode == 'r_dependent':
                mu0 = solve_angular_r_dependent_pk(
                    R, R_e, l_max, R_ref=R_ref, pk_cap=pk_cap,
                )
            else:
                raise ValueError(f"Unknown mode: {mode}")

            U_vals[ire] = (mu0 + 15.0 / 8.0) / (R_e ** 2)

        E_elec[ir] = np.min(U_vals)

    V_nn = Z_A_EFF * Z_B / R_grid
    E_total = E_elec + V_nn

    return {
        'R': R_grid.tolist(),
        'E_total': E_total.tolist(),
        'E_elec': E_elec.tolist(),
    }


def find_r_eq(R_vals: np.ndarray, E_vals: np.ndarray) -> tuple:
    """Find R_eq by cubic interpolation near the minimum."""
    i_min = np.argmin(E_vals)

    if i_min == 0 or i_min == len(R_vals) - 1:
        return float(R_vals[i_min]), float(E_vals[i_min]), False

    i_lo = max(0, i_min - 2)
    i_hi = min(len(R_vals), i_min + 3)
    R_sub = R_vals[i_lo:i_hi]
    E_sub = E_vals[i_lo:i_hi]

    if len(R_sub) < 3:
        return float(R_vals[i_min]), float(E_vals[i_min]), True

    cs = CubicSpline(R_sub, E_sub)
    R_fine = np.linspace(R_sub[0], R_sub[-1], 500)
    E_fine = cs(R_fine)
    i_fine_min = np.argmin(E_fine)

    return float(R_fine[i_fine_min]), float(E_fine[i_fine_min]), True


def main():
    print("=" * 76)
    print("Task 10: Self-Consistent Projected PK Iteration")
    print("=" * 76)
    print(f"Z_A_eff = {Z_A_EFF}, Z_B = {Z_B}")
    print(f"PK: A = {PK_A}, B = {PK_B_PARAM}")
    print(f"Damping = {DAMPING}, max_iter = {MAX_ITER}, tol = {TOL}")
    print(f"Experimental R_eq = {R_EXP} bohr")
    print()

    results = {}

    # ==================================================================
    # PART 1 & 2: Fixed-geometry convergence test
    # ==================================================================
    print("=" * 76)
    print("PART 2: Fixed-Geometry Convergence Test")
    print("=" * 76)

    test_points = [
        (R_EXP, 1.5),   # R_eq, R_e ~ R/2
        (5.0, 2.5),      # Extended geometry
    ]

    for l_max in [2, 4]:
        n_ch = len(get_channels(l_max))
        print(f"\n--- l_max = {l_max} ({n_ch} channels) ---")

        for R, R_e in test_points:
            print(f"\n  R = {R:.3f} bohr, R_e = {R_e:.3f} bohr:")

            # Reference: l-dependent (w=1.0)
            mu0_ldep = solve_angular_ldependent(R, R_e, l_max)
            print(f"    l-dependent (w=1.0):  mu0 = {mu0_ldep:.6f}")

            # Self-consistent iteration
            res = solve_angular_selfconsistent(
                R, R_e, l_max, verbose=True,
            )

            print(f"    Self-consistent:      mu0 = {res['mu0']:.6f}")
            print(f"    Converged: {res['converged']} "
                  f"({res['n_iter']} iterations)")
            print(f"    p1 (e1 l=0 content) = {res['p1']:.4f}")
            print(f"    p2 (e2 l=0 content) = {res['p2']:.4f}")
            print(f"    delta_mu0 = {res['mu0'] - mu0_ldep:.6f} "
                  f"(SC - l_dep)")

            key = f"fixed_lmax{l_max}_R{R:.1f}_Re{R_e:.1f}"
            results[key] = {
                'R': R,
                'R_e': R_e,
                'l_max': l_max,
                'mu0_ldependent': mu0_ldep,
                'mu0_selfconsistent': res['mu0'],
                'p1': res['p1'],
                'p2': res['p2'],
                'converged': res['converged'],
                'n_iter': res['n_iter'],
                'history': res['history'],
            }

    # ==================================================================
    # PART 2b: p1, p2 vs R (at optimal R_e for each R)
    # ==================================================================
    print("\n" + "=" * 76)
    print("PART 2b: Converged PK Weights (p1, p2) vs R")
    print("=" * 76)

    for l_max in [2, 4]:
        print(f"\n--- l_max = {l_max} ---")
        pk_vs_R = []

        for R in [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.5]:
            # Find optimal R_e by scanning
            best_U = np.inf
            best_Re = R_E_GRID[0]
            for R_e in R_E_GRID:
                res = solve_angular_selfconsistent(
                    R, R_e, l_max, verbose=False,
                )
                U = (res['mu0'] + 15.0 / 8.0) / (R_e ** 2)
                if U < best_U:
                    best_U = U
                    best_Re = R_e
                    best_p1 = res['p1']
                    best_p2 = res['p2']
                    best_mu0 = res['mu0']

            pk_vs_R.append({
                'R': R, 'R_e_opt': float(best_Re),
                'p1': float(best_p1), 'p2': float(best_p2),
                'mu0': float(best_mu0),
            })
            print(f"  R={R:.1f}: R_e*={best_Re:.2f}, "
                  f"p1={best_p1:.4f}, p2={best_p2:.4f}")

        results[f'pk_vs_R_lmax{l_max}'] = pk_vs_R

    # ==================================================================
    # PART 3: PES Sweep — Self-Consistent vs L-Dependent
    # ==================================================================
    print("\n" + "=" * 76)
    print("PART 3: PES Sweep — Self-Consistent vs L-Dependent")
    print("=" * 76)

    l_max_values = [2, 3, 4]
    pes_results = {}

    for l_max in l_max_values:
        n_ch = len(get_channels(l_max))
        print(f"\n--- l_max = {l_max} ({n_ch} channels) ---")

        # L-dependent PES
        t0 = time.time()
        print(f"  l-dependent PES...", end=" ", flush=True)
        pes_ldep = compute_pes(l_max, R_VALUES, R_E_GRID, mode='l_dependent')
        t_ldep = time.time() - t0
        R_eq_ldep, E_min_ldep, _ = find_r_eq(
            np.array(pes_ldep['R']), np.array(pes_ldep['E_total']))
        print(f"R_eq={R_eq_ldep:.3f}, E_min={E_min_ldep:.6f}, {t_ldep:.1f}s")

        # Self-consistent PES (limit iterations for speed)
        t0 = time.time()
        print(f"  self-consistent PES (max_iter=10)...", end=" ", flush=True)
        pes_sc = compute_pes(
            l_max, R_VALUES, R_E_GRID, mode='selfconsistent',
            max_iter=10,
        )
        t_sc = time.time() - t0
        R_eq_sc, E_min_sc, _ = find_r_eq(
            np.array(pes_sc['R']), np.array(pes_sc['E_total']))
        print(f"R_eq={R_eq_sc:.3f}, E_min={E_min_sc:.6f}, {t_sc:.1f}s")

        pes_results[l_max] = {
            'l_dependent': {
                'R_eq': R_eq_ldep, 'E_min': E_min_ldep,
                'time': t_ldep, **pes_ldep,
            },
            'selfconsistent': {
                'R_eq': R_eq_sc, 'E_min': E_min_sc,
                'time': t_sc, **pes_sc,
            },
        }

    results['pes_sweep'] = pes_results

    # ==================================================================
    # PART 4: R-Dependent PK Scaling
    # ==================================================================
    print("\n" + "=" * 76)
    print("PART 4: R-Dependent PK Scaling (w = min(cap, R/R_ref))")
    print("=" * 76)

    for pk_cap in [1.5, 2.0]:
        print(f"\n--- pk_cap = {pk_cap} ---")
        rdep_results = {}

        for l_max in l_max_values:
            n_ch = len(get_channels(l_max))
            t0 = time.time()
            print(f"  l_max={l_max} ({n_ch} ch)...", end=" ", flush=True)

            pes_rdep = compute_pes(
                l_max, R_VALUES, R_E_GRID, mode='r_dependent',
                R_ref=R_EXP, pk_cap=pk_cap,
            )
            t_rdep = time.time() - t0
            R_eq_rdep, E_min_rdep, _ = find_r_eq(
                np.array(pes_rdep['R']), np.array(pes_rdep['E_total']))
            print(f"R_eq={R_eq_rdep:.3f}, E_min={E_min_rdep:.6f}, "
                  f"{t_rdep:.1f}s")

            rdep_results[l_max] = {
                'R_eq': R_eq_rdep, 'E_min': E_min_rdep,
                'time': t_rdep, **pes_rdep,
            }

        results[f'r_dependent_cap{pk_cap}'] = rdep_results

    # ==================================================================
    # Summary Table
    # ==================================================================
    print("\n" + "=" * 76)
    print("SUMMARY TABLE")
    print("=" * 76)
    print(f"{'Mode':<25s} {'l_max':>5s} {'R_eq':>8s} {'dR_eq':>8s} "
          f"{'% error':>8s} {'E_min':>12s}")
    print("-" * 76)

    for l_max in l_max_values:
        if l_max in pes_results:
            for mode_name, mode_key in [('l-dependent', 'l_dependent'),
                                         ('self-consistent', 'selfconsistent')]:
                data = pes_results[l_max][mode_key]
                R_eq = data['R_eq']
                dR = R_eq - R_EXP
                pct = abs(dR) / R_EXP * 100
                print(f"{mode_name:<25s} {l_max:>5d} {R_eq:>8.3f} "
                      f"{dR:>+8.3f} {pct:>8.1f} {data['E_min']:>12.6f}")

    for pk_cap in [1.5, 2.0]:
        cap_key = f'r_dependent_cap{pk_cap}'
        if cap_key in results:
            for l_max in l_max_values:
                if l_max in results[cap_key]:
                    data = results[cap_key][l_max]
                    R_eq = data['R_eq']
                    dR = R_eq - R_EXP
                    pct = abs(dR) / R_EXP * 100
                    label = f"R-dep (cap={pk_cap})"
                    print(f"{label:<25s} {l_max:>5d} {R_eq:>8.3f} "
                          f"{dR:>+8.3f} {pct:>8.1f} {data['E_min']:>12.6f}")

    print(f"{'Experiment':<25s} {'':>5s} {R_EXP:>8.3f} "
          f"{0.0:>+8.3f} {0.0:>8.1f}")

    # ==================================================================
    # Convergence analysis: R_eq vs l_max slopes
    # ==================================================================
    print("\n" + "=" * 76)
    print("R_eq DRIFT ANALYSIS")
    print("=" * 76)

    for mode_name, mode_key in [('l-dependent', 'l_dependent'),
                                 ('self-consistent', 'selfconsistent')]:
        lm_arr = []
        req_arr = []
        for l_max in l_max_values:
            if l_max in pes_results:
                lm_arr.append(l_max)
                req_arr.append(pes_results[l_max][mode_key]['R_eq'])
        if len(lm_arr) >= 2:
            lm_arr = np.array(lm_arr)
            req_arr = np.array(req_arr)
            slope = np.polyfit(lm_arr, req_arr, 1)[0]
            print(f"  {mode_name}: slope = {slope:+.3f} bohr/l_max")

    for pk_cap in [1.5, 2.0]:
        cap_key = f'r_dependent_cap{pk_cap}'
        if cap_key in results:
            lm_arr = []
            req_arr = []
            for l_max in l_max_values:
                if l_max in results[cap_key]:
                    lm_arr.append(l_max)
                    req_arr.append(results[cap_key][l_max]['R_eq'])
            if len(lm_arr) >= 2:
                slope = np.polyfit(lm_arr, req_arr, 1)[0]
                print(f"  R-dep (cap={pk_cap}): slope = {slope:+.3f} "
                      f"bohr/l_max")

    # ==================================================================
    # Save results
    # ==================================================================
    outfile = Path(__file__).parent / 'data' / 'selfconsistent_pk.json'
    outfile.parent.mkdir(exist_ok=True)

    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {outfile}")

    # ==================================================================
    # Generate plots
    # ==================================================================
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        plot_dir = Path(__file__).parent / 'plots'
        plot_dir.mkdir(exist_ok=True)

        # --- Plot A: Iteration convergence at fixed geometry ---
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Self-Consistent PK Iteration Convergence', fontsize=14)

        for idx, l_max in enumerate([2, 4]):
            for jdx, (R, R_e) in enumerate(test_points):
                ax = axes[idx, jdx]
                key = f"fixed_lmax{l_max}_R{R:.1f}_Re{R_e:.1f}"
                if key in results:
                    hist = results[key]['history']
                    iters = [h['iteration'] for h in hist]
                    mu0s = [h['mu0'] for h in hist]
                    p1s = [h['p1'] for h in hist]
                    p2s = [h['p2'] for h in hist]

                    ax2 = ax.twinx()
                    ax.plot(iters, mu0s, 'b-o', markersize=4, label='μ₀')
                    ax2.plot(iters, p1s, 'r--s', markersize=3, label='p₁')
                    ax2.plot(iters, p2s, 'g--^', markersize=3, label='p₂')

                    ax.set_xlabel('Iteration')
                    ax.set_ylabel('μ₀', color='b')
                    ax2.set_ylabel('p₁, p₂', color='r')
                    ax.set_title(f'l_max={l_max}, R={R}, R_e={R_e}')

                    # Add l-dependent reference
                    mu0_ref = results[key]['mu0_ldependent']
                    ax.axhline(mu0_ref, color='b', linestyle=':', alpha=0.5,
                               label=f'l-dep μ₀={mu0_ref:.3f}')
                    ax.legend(loc='upper left', fontsize=8)
                    ax2.legend(loc='upper right', fontsize=8)

        plt.tight_layout()
        plt.savefig(plot_dir / 'selfconsistent_pk_convergence.png', dpi=150)
        plt.close()
        print(f"Plot A saved: {plot_dir}/selfconsistent_pk_convergence.png")

        # --- Plot B: Converged PK weights vs R ---
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle('Converged PK Weights (p₁, p₂) vs R', fontsize=14)

        for idx, l_max in enumerate([2, 4]):
            ax = axes[idx]
            key = f'pk_vs_R_lmax{l_max}'
            if key in results:
                data = results[key]
                Rs = [d['R'] for d in data]
                p1s = [d['p1'] for d in data]
                p2s = [d['p2'] for d in data]

                ax.plot(Rs, p1s, 'ro-', label='p₁ (electron 1)')
                ax.plot(Rs, p2s, 'bs-', label='p₂ (electron 2)')
                ax.axvline(R_EXP, color='k', linestyle='--', alpha=0.5,
                           label=f'R_exp={R_EXP}')
                ax.set_xlabel('R (bohr)')
                ax.set_ylabel('l=0 content')
                ax.set_title(f'l_max = {l_max}')
                ax.legend()
                ax.set_ylim(0, 1.1)
                ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(plot_dir / 'selfconsistent_pk_weights_vs_R.png', dpi=150)
        plt.close()
        print(f"Plot B saved: {plot_dir}/selfconsistent_pk_weights_vs_R.png")

        # --- Plot C: R_eq vs l_max comparison ---
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.set_title('R_eq vs l_max: PK Mode Comparison', fontsize=14)

        markers = {'l_dependent': 'o', 'selfconsistent': 's'}
        colors = {'l_dependent': 'blue', 'selfconsistent': 'red'}
        labels = {'l_dependent': 'l-dependent (w=1)',
                  'selfconsistent': 'self-consistent'}

        for mode_key in ['l_dependent', 'selfconsistent']:
            lm_arr = []
            req_arr = []
            for l_max in l_max_values:
                if l_max in pes_results:
                    lm_arr.append(l_max)
                    req_arr.append(pes_results[l_max][mode_key]['R_eq'])
            ax.plot(lm_arr, req_arr, marker=markers[mode_key],
                    color=colors[mode_key], label=labels[mode_key],
                    linewidth=2, markersize=8)

        # R-dependent PK
        cap_colors = {1.5: 'green', 2.0: 'purple'}
        for pk_cap in [1.5, 2.0]:
            cap_key = f'r_dependent_cap{pk_cap}'
            if cap_key in results:
                lm_arr = []
                req_arr = []
                for l_max in l_max_values:
                    if l_max in results[cap_key]:
                        lm_arr.append(l_max)
                        req_arr.append(results[cap_key][l_max]['R_eq'])
                ax.plot(lm_arr, req_arr, marker='D',
                        color=cap_colors[pk_cap],
                        label=f'R-dependent (cap={pk_cap})',
                        linewidth=2, markersize=8, linestyle='--')

        ax.axhline(R_EXP, color='k', linestyle=':', linewidth=2,
                    label=f'Experiment ({R_EXP})')
        ax.set_xlabel('l_max')
        ax.set_ylabel('R_eq (bohr)')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(plot_dir / 'selfconsistent_pk_req_vs_lmax.png', dpi=150)
        plt.close()
        print(f"Plot C saved: {plot_dir}/selfconsistent_pk_req_vs_lmax.png")

        # --- Plot D: PES curves comparison at l_max=4 ---
        if 4 in pes_results:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.set_title('PES Curves at l_max=4: PK Mode Comparison',
                          fontsize=14)

            for mode_key, label, color in [
                ('l_dependent', 'l-dependent', 'blue'),
                ('selfconsistent', 'self-consistent', 'red'),
            ]:
                data = pes_results[4][mode_key]
                ax.plot(data['R'], data['E_total'], '-o', color=color,
                        label=label, markersize=4)

            # R-dependent
            for pk_cap, color in [(1.5, 'green'), (2.0, 'purple')]:
                cap_key = f'r_dependent_cap{pk_cap}'
                if cap_key in results and 4 in results[cap_key]:
                    data = results[cap_key][4]
                    ax.plot(data['R'], data['E_total'], '--D', color=color,
                            label=f'R-dep (cap={pk_cap})', markersize=4)

            ax.axvline(R_EXP, color='k', linestyle=':', alpha=0.5)
            ax.set_xlabel('R (bohr)')
            ax.set_ylabel('E_total (Ha)')
            ax.legend()
            ax.grid(True, alpha=0.3)

            plt.tight_layout()
            plt.savefig(plot_dir / 'selfconsistent_pk_pes_lmax4.png',
                        dpi=150)
            plt.close()
            print(f"Plot D saved: "
                  f"{plot_dir}/selfconsistent_pk_pes_lmax4.png")

    except ImportError:
        print("matplotlib not available — skipping plots")

    print("\nDone.")


if __name__ == '__main__':
    main()
