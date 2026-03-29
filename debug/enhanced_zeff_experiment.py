"""
Track A, Task 5: PK-free composed geometry with enhanced Z_eff.

Hypothesis: encoding core exclusion within Z_eff(r) rather than as a
separate PK barrier preserves angular channel symmetry and reduces
l_max divergence.

Enhanced Z_eff:
    Z_eff_enhanced(r) = Z_eff_standard(r) - A_rep * exp(-r^2 / (2 sigma^2))

Two A_rep values are tested:
  - "energy_matched": A_rep set so integrated repulsive potential energy
    matches the PK barrier (A_rep~40.9, Z_eff(0)~-37.9). Very strong.
  - "moderate": A_rep=4 giving Z_eff(0)=-1. Mild repulsion.

Experiment:
  1. Solve core, derive PK parameters
  2. Compute enhanced Z_eff curves
  3. Run composed LiH at l_max = 2, 3, 4 with enhanced Z_eff, NO PK
  4. Extract channel weights at l_max=4
  5. Compare drift rates with standard PK results
"""

import numpy as np
import json
import time
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import curve_fit

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.level4_multichannel import (
    build_angular_hamiltonian,
    _channel_list,
)
from geovac.composed_diatomic import _v_cross_nuc_1s


# ---------------------------------------------------------------------------
# Enhanced Z_eff construction
# ---------------------------------------------------------------------------

def make_enhanced_z_eff(core: CoreScreening, sigma: float, A_rep: float):
    """Build enhanced Z_eff(r) callable."""
    def z_eff_enhanced(r_val):
        r_arr = np.atleast_1d(np.asarray(r_val, dtype=float))
        z_std = core.z_eff(np.maximum(r_arr, 1e-6))
        repulsion = A_rep * np.exp(-r_arr**2 / (2 * sigma**2))
        result = z_std - repulsion
        if np.ndim(r_val) == 0:
            return float(result[0])
        return result
    return z_eff_enhanced


def compute_energy_matched_A_rep(core: CoreScreening, pk: AbInitioPK) -> float:
    """Compute A_rep that matches the PK barrier's integrated energy."""
    sigma = pk.r_core
    r = np.linspace(0.002, 15.0, 50000)
    psi2s_r2 = (1.0 / 8.0) * (2.0 - r)**2 * np.exp(-r) * r**2

    # PK energy: integral V_PK(r) |psi|^2 r^2 dr
    V_pk = pk.A * np.exp(-pk.B * r**2) / r**2
    E_pk = np.trapezoid(V_pk * psi2s_r2, r)

    # Enhanced Z_eff repulsion potential: A_rep * exp(-r^2/(2s^2)) / r
    # Energy: A_rep * integral [exp(-r^2/(2s^2)) / r] |psi|^2 r^2 dr
    integrand_rep = np.exp(-r**2 / (2 * sigma**2)) * psi2s_r2 / r
    E_rep_per_Arep = np.trapezoid(integrand_rep, r)

    return E_pk / E_rep_per_Arep, E_pk


# ---------------------------------------------------------------------------
# Adiabatic curve + radial solve with Z_A_func
# ---------------------------------------------------------------------------

def solve_valence_with_zeff(
    R: float, Z_A_func: callable, Z_B: float, l_max: int,
    n_alpha: int = 100, n_Re: int = 300,
    R_e_min: float = 0.1, R_e_max: float = 15.0,
    z0: float = 0.0, n_theta: int = 64,
) -> float:
    """Solve Level 4 valence problem with Z_A_func (no PK)."""
    Z_A_bare = float(Z_A_func(0.0))

    R_e_angular = np.concatenate([
        np.linspace(R_e_min, 1.0, 40),
        np.linspace(1.0, 3.0, 40),
        np.linspace(3.0, 6.0, 30),
        np.linspace(6.0, R_e_max, 20),
    ])
    R_e_angular = np.unique(R_e_angular)

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    n_Re_ang = len(R_e_angular)
    mu_vals = np.zeros(n_Re_ang)

    for i, R_e in enumerate(R_e_angular):
        rho = R / (2.0 * R_e)
        H = build_angular_hamiltonian(
            alpha, rho, R_e, l_max, Z=1.0,
            m_max=0, Z_A=Z_A_bare, Z_B=Z_B, z0=z0,
            Z_A_func=Z_A_func, n_theta=n_theta,
            pk_potentials=None, pk_projector=None,
        )
        evals = eigh(H, eigvals_only=True)
        mu_vals[i] = evals[0]

    U = (mu_vals + 15.0 / 8.0) / R_e_angular**2

    U_spline = CubicSpline(R_e_angular, U, extrapolate=True)
    h_Re = (R_e_max - R_e_min) / (n_Re + 1)
    R_e_radial = R_e_min + (np.arange(n_Re) + 1) * h_Re
    V_radial = U_spline(R_e_radial)

    diag = np.ones(n_Re) / h_Re**2 + V_radial
    off_diag = -0.5 * np.ones(n_Re - 1) / h_Re**2

    evals_rad, _ = eigh_tridiagonal(
        diag, off_diag, select='i', select_range=(0, 0),
    )
    return float(evals_rad[0])


# ---------------------------------------------------------------------------
# Channel weight extraction
# ---------------------------------------------------------------------------

def get_channel_weights(
    R: float, R_e_target: float, Z_A_func: callable,
    Z_B: float, l_max: int, n_alpha: int = 100, z0: float = 0.0,
) -> tuple:
    """Extract channel weight distribution at (R, R_e)."""
    Z_A_bare = float(Z_A_func(0.0))
    homonuclear = abs(Z_A_bare - Z_B) < 1e-10

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h
    rho = R / (2.0 * R_e_target)

    H = build_angular_hamiltonian(
        alpha, rho, R_e_target, l_max, Z=1.0,
        m_max=0, Z_A=Z_A_bare, Z_B=Z_B, z0=z0,
        Z_A_func=Z_A_func,
        pk_potentials=None, pk_projector=None,
    )

    _, evecs = eigh(H)
    vec = evecs[:, 0]

    channels = _channel_list(l_max, homonuclear=homonuclear)

    weights = {}
    for ic, ch in enumerate(channels):
        w = np.sum(vec[ic * n_alpha:(ic + 1) * n_alpha]**2) * h
        weights[str(ch)] = float(w)

    total = sum(weights.values())
    if total > 0:
        weights = {k: v / total for k, v in weights.items()}

    return weights, channels


# ---------------------------------------------------------------------------
# Morse fit
# ---------------------------------------------------------------------------

def _morse(R, E_min, D_e, a, R_eq):
    return E_min + D_e * (1.0 - np.exp(-a * (R - R_eq)))**2


def fit_morse(R_grid, E_grid):
    valid = ~np.isnan(E_grid)
    R_v, E_v = R_grid[valid], E_grid[valid]
    if len(R_v) < 4:
        return {'R_eq': np.nan, 'E_min': np.nan, 'D_e': np.nan}

    i_min = np.argmin(E_v)
    R_eq_g, E_min_g = R_v[i_min], E_v[i_min]
    D_e_g = max(E_v[-1] - E_min_g, 0.005)

    try:
        popt, _ = curve_fit(
            _morse, R_v, E_v,
            p0=[E_min_g, D_e_g, 1.0, R_eq_g], maxfev=10000,
        )
        return {'R_eq': float(popt[3]), 'E_min': float(popt[0]),
                'D_e': float(abs(popt[1])), 'a': float(abs(popt[2]))}
    except RuntimeError:
        return {'R_eq': float(R_eq_g), 'E_min': float(E_min_g),
                'D_e': float(D_e_g)}


# ---------------------------------------------------------------------------
# PES scan for a single (mode, l_max)
# ---------------------------------------------------------------------------

def scan_pes(
    z_eff_func, l_max, E_core, R_grid,
    Z_A_bare_nuclear=3.0, Z_B=1.0, n_core=2,
    n_alpha=100, n_Re=300, label="",
):
    """Scan PES for a single enhanced Z_eff mode and l_max."""
    Z_A_bare_eff = float(z_eff_func(0.0))
    homonuclear = abs(Z_A_bare_eff - Z_B) < 1e-10
    channels = _channel_list(l_max, homonuclear=homonuclear)
    n_ch = len(channels)

    print(f"\n  [{label}] l_max={l_max}, Z_A_bare={Z_A_bare_eff:.2f}, "
          f"n_ch={n_ch}, homo={homonuclear}")
    print(f"  {'R':>6s} {'E_elec':>10s} {'E_comp':>12s} {'time':>6s}")

    n_R = len(R_grid)
    E_comp = np.zeros(n_R)
    E_elec_arr = np.zeros(n_R)

    for i, R in enumerate(R_grid):
        t0 = time.time()
        try:
            E_el = solve_valence_with_zeff(
                R, z_eff_func, Z_B, l_max,
                n_alpha=n_alpha, n_Re=n_Re,
            )
            V_NN = Z_A_bare_nuclear * Z_B / R
            V_cross = _v_cross_nuc_1s(Z_A_bare_nuclear, n_core, Z_B, R)
            E_c = E_core + E_el + V_NN + V_cross
            E_comp[i] = E_c
            E_elec_arr[i] = E_el
            dt = time.time() - t0
            print(f"  {R:6.3f} {E_el:10.4f} {E_c:12.6f} {dt:6.1f}s")
        except Exception as e:
            E_comp[i] = np.nan
            dt = time.time() - t0
            print(f"  {R:6.3f}  FAILED: {e}  ({dt:.1f}s)")

    morse = fit_morse(R_grid, E_comp)
    R_eq = morse['R_eq']
    err = (R_eq - 3.015) / 3.015 * 100 if not np.isnan(R_eq) else np.nan

    # Check for monotonic PES
    valid = ~np.isnan(E_comp)
    i_min = np.argmin(E_comp[valid]) if np.any(valid) else 0
    is_mono = (i_min == 0 or i_min == np.sum(valid) - 1)

    print(f"  R_eq={R_eq:.3f} (err={err:+.1f}%), "
          f"E_min={morse['E_min']:.6f}, D_e={morse['D_e']:.6f}")
    if is_mono:
        print(f"  WARNING: PES minimum at boundary (monotonic)")

    return {
        'R_grid': R_grid.tolist(),
        'E_composed': E_comp.tolist(),
        'E_elec': E_elec_arr.tolist(),
        'R_eq': float(R_eq),
        'E_min': float(morse['E_min']),
        'D_e': float(morse['D_e']),
        'R_eq_error_pct': float(err) if not np.isnan(err) else None,
        'n_channels': n_ch,
        'homonuclear': homonuclear,
        'is_monotonic': bool(is_mono),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 64)
    print("Track A, Task 5: Enhanced Z_eff Experiment")
    print("=" * 64)
    t_total = time.time()

    # --- 1. Solve core ---
    print("\n--- Step 1: Solve Li 1s^2 core ---")
    core = CoreScreening(Z=3, l_max=2, n_alpha=200)
    core.solve(verbose=False)
    E_core = core.energy
    pk = AbInitioPK(core, n_core=2)
    sigma = pk.r_core
    print(f"  E_core = {E_core:.6f} Ha")
    print(f"  PK: A={pk.A:.4f}, B={pk.B:.4f}, r_core={sigma:.4f}")

    # --- 2. Compute A_rep values ---
    print("\n--- Step 2: Enhanced Z_eff parameters ---")
    A_rep_matched, E_pk = compute_energy_matched_A_rep(core, pk)
    A_rep_moderate = 4.0  # Z_eff(0) = 3 - 4 = -1

    print(f"  Energy-matched: A_rep={A_rep_matched:.2f}, "
          f"Z_eff(0)={core.z_eff(0.001)-A_rep_matched:.2f}")
    print(f"  Moderate:       A_rep={A_rep_moderate:.2f}, "
          f"Z_eff(0)={core.z_eff(0.001)-A_rep_moderate:.2f}")
    print(f"  E_pk per electron = {E_pk:.6f} Ha")

    # Tabulate both curves
    modes = {
        'moderate': make_enhanced_z_eff(core, sigma, A_rep_moderate),
        'energy_matched': make_enhanced_z_eff(core, sigma, A_rep_matched),
    }

    r_tab = [0.0, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0]
    print(f"\n  {'r':>6s} {'Z_std':>8s} {'Z_mod':>8s} {'Z_ematch':>8s}")
    for rv in r_tab:
        rv_s = max(rv, 1e-6)
        z_s = core.z_eff(rv_s)
        z_m = modes['moderate'](rv_s)
        z_e = modes['energy_matched'](rv_s)
        print(f"  {rv:6.2f} {z_s:8.4f} {z_m:8.4f} {z_e:8.4f}")

    # --- 3. PES scans ---
    R_grid = np.concatenate([
        np.linspace(2.0, 2.5, 3),
        np.linspace(2.7, 4.0, 10),
        np.linspace(4.5, 7.0, 5),
    ])

    results = {
        'E_core': float(E_core),
        'sigma': float(sigma),
        'A_rep_moderate': float(A_rep_moderate),
        'A_rep_energy_matched': float(A_rep_matched),
        'E_pk_per_electron': float(E_pk),
        'reference_R_eq': 3.015,
    }

    for mode_name, z_eff_func in modes.items():
        print(f"\n{'='*64}")
        print(f"Mode: {mode_name}")
        print(f"{'='*64}")

        for l_max in [2, 3, 4]:
            key = f'{mode_name}_l{l_max}'
            results[key] = scan_pes(
                z_eff_func, l_max, E_core, R_grid,
                label=mode_name, n_alpha=100, n_Re=300,
            )

    # --- 4. Channel weights at l_max=4 ---
    print(f"\n{'='*64}")
    print("Channel weights at l_max=4")
    print(f"{'='*64}")

    for mode_name, z_eff_func in modes.items():
        key = f'{mode_name}_l4'
        R_eq = results[key]['R_eq']
        if np.isnan(R_eq) or R_eq < 1.5 or R_eq > 8.0:
            R_eq = 3.0  # fallback
        R_e_target = 1.5

        weights, ch_list = get_channel_weights(
            R_eq, R_e_target, z_eff_func, 1.0, 4, n_alpha=100,
        )

        print(f"\n  [{mode_name}] R={R_eq:.3f}, R_e={R_e_target}")
        sorted_w = sorted(weights.items(), key=lambda x: -x[1])
        odd_total = 0.0
        for ch_str, w in sorted_w[:10]:
            if w > 0.001:
                print(f"    {ch_str}: {w*100:.1f}%")
            ch = eval(ch_str)
            if (ch[0] + ch[1]) % 2 != 0:
                odd_total += w

        print(f"    Odd-l total: {odd_total*100:.1f}%")
        print(f"    (0,0) weight: {weights.get('(0, 0)', 0.0)*100:.1f}%")

        results[f'{mode_name}_weights_l4'] = {
            'R': float(R_eq), 'R_e': R_e_target,
            'weights': weights, 'odd_l_total': float(odd_total),
        }

    # --- 5. Save ---
    with open('debug/data/enhanced_zeff_experiment.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nData saved to debug/data/enhanced_zeff_experiment.json")

    # --- 6. Summary ---
    print(f"\n{'='*64}")
    print("SUMMARY")
    print(f"{'='*64}")

    print(f"\nR_eq vs l_max:")
    print(f"  {'Mode':>20s} {'l=2':>8s} {'l=3':>8s} {'l=4':>8s} "
          f"{'drift':>12s}")
    print(f"  {'-'*20} {'-'*8} {'-'*8} {'-'*8} {'-'*12}")

    for mode_name in ['moderate', 'energy_matched']:
        r_vals = []
        for l in [2, 3, 4]:
            key = f'{mode_name}_l{l}'
            r = results[key]['R_eq']
            r_vals.append(r)
        if not any(np.isnan(r_vals)):
            drift = (r_vals[2] - r_vals[0]) / 2.0
        else:
            drift = np.nan
        print(f"  {mode_name:>20s} {r_vals[0]:8.3f} {r_vals[1]:8.3f} "
              f"{r_vals[2]:8.3f} {drift:+11.3f}")

    print(f"  {'l_dep PK (ref)':>20s} {'3.176':>8s} {'3.488':>8s} "
          f"{'3.783':>8s} {'+0.303':>12s}")
    print(f"  {'ch_blind PK (ref)':>20s} {'3.222':>8s} {'3.583':>8s} "
          f"{'3.956':>8s} {'+0.367':>12s}")
    print(f"  {'HeH+ bare (ref)':>20s} {'1.464':>8s} {'1.651':>8s} "
          f"{'1.989':>8s} {'+0.262':>12s}")
    print(f"  {'Experiment':>20s} {'':>8s} {'':>8s} "
          f"{'':>8s} {'3.015':>12s}")

    elapsed = time.time() - t_total
    print(f"\nTotal: {elapsed:.0f}s ({elapsed/60:.1f} min)")


if __name__ == '__main__':
    main()
