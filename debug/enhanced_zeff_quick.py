"""
Quick enhanced Z_eff experiment — focused on key results.

Tests two A_rep values (moderate=4, energy_matched=40.88) at l_max=2,4
with reduced R grid (8 points) and n_alpha=80 for speed.
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


def make_z_eff(core, sigma, A_rep):
    def z_eff(r_val):
        r_arr = np.atleast_1d(np.asarray(r_val, dtype=float))
        z_std = core.z_eff(np.maximum(r_arr, 1e-6))
        rep = A_rep * np.exp(-r_arr**2 / (2 * sigma**2))
        result = z_std - rep
        if np.ndim(r_val) == 0:
            return float(result[0])
        return result
    return z_eff


def solve_valence(R, Z_A_func, Z_B, l_max, n_alpha=80, n_Re=250,
                   R_e_min=0.1, R_e_max=15.0, z0=0.0):
    Z_A_bare = float(Z_A_func(0.0))

    R_e_ang = np.concatenate([
        np.linspace(R_e_min, 1.0, 30),
        np.linspace(1.0, 3.0, 30),
        np.linspace(3.0, 6.0, 20),
        np.linspace(6.0, R_e_max, 15),
    ])
    R_e_ang = np.unique(R_e_ang)

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    mu_vals = np.zeros(len(R_e_ang))
    for i, R_e in enumerate(R_e_ang):
        rho = R / (2.0 * R_e)
        H = build_angular_hamiltonian(
            alpha, rho, R_e, l_max, Z=1.0,
            m_max=0, Z_A=Z_A_bare, Z_B=Z_B, z0=z0,
            Z_A_func=Z_A_func, n_theta=48,
            pk_potentials=None, pk_projector=None,
        )
        evals = eigh(H, eigvals_only=True)
        mu_vals[i] = evals[0]

    U = (mu_vals + 15.0 / 8.0) / R_e_ang**2
    U_spline = CubicSpline(R_e_ang, U, extrapolate=True)

    h_Re = (R_e_max - R_e_min) / (n_Re + 1)
    R_e_rad = R_e_min + (np.arange(n_Re) + 1) * h_Re
    V_rad = U_spline(R_e_rad)

    diag = np.ones(n_Re) / h_Re**2 + V_rad
    off = -0.5 * np.ones(n_Re - 1) / h_Re**2
    ev, _ = eigh_tridiagonal(diag, off, select='i', select_range=(0, 0))
    return float(ev[0])


def get_weights(R, R_e, Z_A_func, Z_B, l_max, n_alpha=80):
    Z_A_bare = float(Z_A_func(0.0))
    homo = abs(Z_A_bare - Z_B) < 1e-10

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h
    rho = R / (2.0 * R_e)

    H = build_angular_hamiltonian(
        alpha, rho, R_e, l_max, Z=1.0,
        m_max=0, Z_A=Z_A_bare, Z_B=Z_B,
        Z_A_func=Z_A_func,
        pk_potentials=None, pk_projector=None,
    )
    _, evecs = eigh(H)
    vec = evecs[:, 0]

    channels = _channel_list(l_max, homonuclear=homo)
    weights = {}
    for ic, ch in enumerate(channels):
        w = np.sum(vec[ic * n_alpha:(ic + 1) * n_alpha]**2) * h
        weights[str(ch)] = float(w)

    total = sum(weights.values())
    if total > 0:
        weights = {k: v / total for k, v in weights.items()}
    return weights, channels


def _morse(R, E_min, D_e, a, R_eq):
    return E_min + D_e * (1.0 - np.exp(-a * (R - R_eq)))**2


def main():
    print("=" * 64)
    print("Enhanced Z_eff Quick Experiment")
    print("=" * 64)
    t_total = time.time()

    # Core solve
    print("\nSolving core...")
    core = CoreScreening(Z=3, l_max=2, n_alpha=200)
    core.solve(verbose=False)
    E_core = core.energy
    pk = AbInitioPK(core, n_core=2)
    sigma = pk.r_core
    print(f"  E_core={E_core:.6f}, PK: A={pk.A:.4f}, B={pk.B:.4f}, "
          f"r_core={sigma:.4f}")

    # Energy-matched A_rep
    r = np.linspace(0.002, 15.0, 50000)
    psi2s_r2 = (1.0/8.0) * (2.0-r)**2 * np.exp(-r) * r**2
    V_pk = pk.A * np.exp(-pk.B * r**2) / r**2
    E_pk = np.trapezoid(V_pk * psi2s_r2, r)
    integ = np.trapezoid(np.exp(-r**2/(2*sigma**2)) * psi2s_r2 / r, r)
    A_rep_em = E_pk / integ

    print(f"\n  A_rep (energy-matched) = {A_rep_em:.2f} -> Z_eff(0)={3-A_rep_em:.2f}")
    print(f"  A_rep (moderate) = 4.00 -> Z_eff(0)=-1.00")

    # R grid (reduced)
    R_grid = np.array([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.5, 7.0])
    Z_B = 1.0
    n_alpha = 80
    n_Re = 250

    results = {'E_core': float(E_core), 'sigma': float(sigma)}

    modes = {
        'moderate': (4.0, make_z_eff(core, sigma, 4.0)),
        'energy_matched': (A_rep_em, make_z_eff(core, sigma, A_rep_em)),
    }

    for mode_name, (A_rep_val, z_func) in modes.items():
        Z_A_bare = float(z_func(0.0))
        print(f"\n{'='*64}")
        print(f"Mode: {mode_name} (A_rep={A_rep_val:.2f}, "
              f"Z_A_bare={Z_A_bare:.2f})")
        print(f"{'='*64}")

        for l_max in [2, 4]:
            homo = abs(Z_A_bare - Z_B) < 1e-10
            channels = _channel_list(l_max, homonuclear=homo)
            n_ch = len(channels)

            print(f"\n  l_max={l_max}, n_ch={n_ch}, homo={homo}")
            print(f"  {'R':>6s} {'E_elec':>10s} {'E_comp':>12s} {'dt':>6s}")

            E_comp = np.zeros(len(R_grid))
            for i, R in enumerate(R_grid):
                t0 = time.time()
                try:
                    E_el = solve_valence(R, z_func, Z_B, l_max,
                                          n_alpha=n_alpha, n_Re=n_Re)
                    V_NN = 3.0 * Z_B / R
                    V_cross = _v_cross_nuc_1s(3.0, 2, Z_B, R)
                    E_comp[i] = E_core + E_el + V_NN + V_cross
                    dt = time.time() - t0
                    print(f"  {R:6.3f} {E_el:10.4f} {E_comp[i]:12.6f} "
                          f"{dt:6.1f}s")
                except Exception as e:
                    E_comp[i] = np.nan
                    print(f"  {R:6.3f}  FAILED: {e}")

            # Morse fit
            valid = ~np.isnan(E_comp)
            R_v, E_v = R_grid[valid], E_comp[valid]
            i_min = np.argmin(E_v)
            is_mono = (i_min == 0 or i_min == len(E_v) - 1)

            R_eq = R_v[i_min]
            if not is_mono and len(E_v) >= 4:
                try:
                    popt, _ = curve_fit(
                        _morse, R_v, E_v,
                        p0=[E_v[i_min], max(E_v[-1]-E_v[i_min], 0.005),
                            1.0, R_v[i_min]], maxfev=10000)
                    R_eq = float(popt[3])
                except RuntimeError:
                    pass

            err = (R_eq - 3.015) / 3.015 * 100
            print(f"  -> R_eq={R_eq:.3f} (err={err:+.1f}%)"
                  f"{'  MONOTONIC' if is_mono else ''}")

            key = f'{mode_name}_l{l_max}'
            results[key] = {
                'R_grid': R_grid.tolist(),
                'E_composed': E_comp.tolist(),
                'R_eq': float(R_eq),
                'R_eq_error_pct': float(err),
                'n_channels': n_ch,
                'is_monotonic': bool(is_mono),
            }

    # Channel weights at l_max=4
    print(f"\n{'='*64}")
    print("Channel weights (l_max=4)")
    print(f"{'='*64}")

    for mode_name, (_, z_func) in modes.items():
        key = f'{mode_name}_l4'
        R_eq = results[key]['R_eq']
        if R_eq < 1.5 or R_eq > 8.0:
            R_eq = 3.0

        weights, _ = get_weights(R_eq, 1.5, z_func, Z_B, 4, n_alpha)
        sorted_w = sorted(weights.items(), key=lambda x: -x[1])
        odd = sum(v for k, v in weights.items()
                  if (eval(k)[0] + eval(k)[1]) % 2 != 0)

        print(f"\n  [{mode_name}] R={R_eq:.3f}")
        for ch, w in sorted_w[:8]:
            if w > 0.001:
                print(f"    {ch}: {w*100:.1f}%")
        print(f"    Odd-l total: {odd*100:.1f}%")

        results[f'{mode_name}_weights_l4'] = {
            'weights': weights, 'odd_l_total': float(odd),
        }

    # Save
    with open('debug/data/enhanced_zeff_experiment.json', 'w') as f:
        json.dump(results, f, indent=2)

    # Summary
    print(f"\n{'='*64}")
    print("SUMMARY")
    print(f"{'='*64}")
    print(f"\n{'Mode':>22s} {'l=2':>8s} {'l=4':>8s} {'drift/l':>10s}")
    print(f"{'-'*22} {'-'*8} {'-'*8} {'-'*10}")

    for mn in ['moderate', 'energy_matched']:
        r2 = results[f'{mn}_l2']['R_eq']
        r4 = results[f'{mn}_l4']['R_eq']
        drift = (r4 - r2) / 2.0
        m2 = '*' if results[f'{mn}_l2']['is_monotonic'] else ''
        m4 = '*' if results[f'{mn}_l4']['is_monotonic'] else ''
        print(f"  {mn:>20s} {r2:7.3f}{m2} {r4:7.3f}{m4} {drift:+9.3f}")

    print(f"  {'l_dep PK (ref)':>20s} {'3.176':>8s} {'3.783':>8s} "
          f"{'+0.303':>10s}")
    print(f"  {'ch_blind PK (ref)':>20s} {'3.222':>8s} {'3.956':>8s} "
          f"{'+0.367':>10s}")
    print(f"  {'HeH+ bare (ref)':>20s} {'1.464':>8s} {'1.989':>8s} "
          f"{'+0.262':>10s}")
    print(f"\n  * = monotonic PES (minimum at boundary)")
    print(f"  Reference R_eq = 3.015 bohr")

    elapsed = time.time() - t_total
    print(f"\nTotal: {elapsed:.0f}s ({elapsed/60:.1f} min)")


if __name__ == '__main__':
    main()
