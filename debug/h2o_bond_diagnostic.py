"""
Diagnostic sweep: Level 4 bond pair accuracy vs charge asymmetry.

Tests whether the Level 4 solver degrades at high Z_A/Z_B ratios,
as relevant to H2O (Z_eff=6, Z_B=1).

Sections:
  1. Bare bond pair scan: Z_A = 1,2,3,4,6,8 with Z_B=1, no PK
  2. PES shape comparison: Z_A=2 vs Z_A=6
  3. Z_eff(r) profile for O (Z=8)
  4. Bond pair at R=1.81 with Z_eff injection vs constant Z_A=6
  5. Bond pair at R=1.81 with PK pseudopotential for Z=8

Output: debug/data/h2o_bond_diagnostic.json
"""

import json
import time
import numpy as np
import sys
import os

# Ensure project root is on path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.level4_multichannel import (
    solve_level4_h2_multichannel,
    build_angular_hamiltonian,
    _channel_list,
    compute_adiabatic_curve_mc,
)
from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal


# Solver parameters (match H2O solver)
L_MAX = 2
N_ALPHA = 100
N_RE = 300

# Reference R_eq values (bohr)
REFS = {
    1: {'name': 'H2', 'R_eq': 1.401, 'E_exact': -1.17447},
    2: {'name': 'HeH+', 'R_eq': 1.463, 'E_exact': -2.9787},
    3: {'name': 'LiH2+-like'},
    4: {'name': 'BeH3+-like'},
    6: {'name': 'OH5+-like (O valence)'},
    8: {'name': 'NeH7+-like'},
}


def solve_bond_pair_pes(Z_A: float, Z_B: float = 1.0,
                        R_grid: np.ndarray = None,
                        pk_potentials=None,
                        label: str = '') -> dict:
    """Solve Level 4 bond pair PES and extract R_eq."""
    if R_grid is None:
        # Adaptive range: higher Z_A -> tighter bond expected
        R_min = 0.3
        R_max = 6.0
        R_grid = np.concatenate([
            np.arange(R_min, 1.0, 0.1),
            np.arange(1.0, 3.0, 0.1),
            np.arange(3.0, R_max + 0.01, 0.2),
        ])

    energies = []
    print(f"\n{'='*60}")
    print(f"PES scan: Z_A={Z_A}, Z_B={Z_B}  {label}")
    print(f"  l_max={L_MAX}, n_alpha={N_ALPHA}, n_Re={N_RE}")
    print(f"  R range: {R_grid[0]:.2f} to {R_grid[-1]:.2f}, {len(R_grid)} points")
    print(f"  {'R':>6s}  {'E_total':>12s}  {'E_elec':>12s}")
    print(f"  {'-'*6}  {'-'*12}  {'-'*12}")

    for R in R_grid:
        try:
            result = solve_level4_h2_multichannel(
                R=R, Z_A=Z_A, Z_B=Z_B,
                l_max=L_MAX, n_alpha=N_ALPHA, n_Re=N_RE,
                verbose=False,
                pk_potentials=pk_potentials,
            )
            E_total = result['E_total']
            E_elec = result['E_elec']
            energies.append({
                'R': float(R),
                'E_total': float(E_total),
                'E_elec': float(E_elec),
            })
            print(f"  {R:6.3f}  {E_total:12.6f}  {E_elec:12.6f}")
        except Exception as e:
            print(f"  {R:6.3f}  FAILED: {e}")
            energies.append({'R': float(R), 'E_total': None, 'E_elec': None})

    # Find minimum
    valid = [e for e in energies if e['E_total'] is not None]
    if not valid:
        return {'Z_A': Z_A, 'Z_B': Z_B, 'error': 'all points failed'}

    E_arr = np.array([e['E_total'] for e in valid])
    R_arr = np.array([e['R'] for e in valid])
    i_min = np.argmin(E_arr)
    R_eq = R_arr[i_min]
    E_min = E_arr[i_min]

    # Rough D_e from dissociation limit
    E_dissoc = E_arr[-1]
    D_e = E_dissoc - E_min

    print(f"\n  R_eq = {R_eq:.3f} bohr")
    print(f"  E_min = {E_min:.6f} Ha")
    print(f"  D_e = {D_e:.6f} Ha")
    print(f"  Bound: {'YES' if D_e > 0 else 'NO'}")

    return {
        'Z_A': float(Z_A),
        'Z_B': float(Z_B),
        'R_eq': float(R_eq),
        'E_min': float(E_min),
        'D_e': float(D_e),
        'energies': energies,
        'label': label,
    }


def section1_charge_asymmetry_scan() -> list:
    """Section 1: Bare bond pair scan across Z_A values."""
    print("\n" + "=" * 70)
    print("SECTION 1: Bare bond pair scan (no PK, no screening)")
    print("=" * 70)

    results = []
    for Z_A in [1, 2, 3, 4, 6, 8]:
        t0 = time.time()
        res = solve_bond_pair_pes(Z_A, Z_B=1.0, label=REFS[Z_A].get('name', ''))
        res['time_s'] = time.time() - t0

        ref = REFS[Z_A]
        if 'R_eq' in ref:
            err = abs(res['R_eq'] - ref['R_eq']) / ref['R_eq'] * 100
            res['R_eq_ref'] = ref['R_eq']
            res['R_eq_error_pct'] = err
            print(f"  R_eq error: {err:.1f}% (ref: {ref['R_eq']:.3f})")
        print(f"  Time: {res['time_s']:.1f}s")
        results.append(res)

    # Summary table
    print("\n" + "=" * 70)
    print("SUMMARY: R_eq vs charge asymmetry")
    print("=" * 70)
    print(f"  {'Z_A':>4s}  {'Z_A/Z_B':>7s}  {'R_eq':>8s}  {'R_eq_ref':>8s}"
          f"  {'Err%':>6s}  {'D_e':>8s}  {'Bound':>5s}")
    print(f"  {'-'*4}  {'-'*7}  {'-'*8}  {'-'*8}  {'-'*6}  {'-'*8}  {'-'*5}")
    for r in results:
        ref_str = f"{r.get('R_eq_ref', '---'):>8.3f}" if 'R_eq_ref' in r else '     ---'
        err_str = f"{r.get('R_eq_error_pct', 0):>6.1f}" if 'R_eq_error_pct' in r else '   ---'
        bound = 'YES' if r['D_e'] > 0 else 'NO'
        print(f"  {r['Z_A']:4.0f}  {r['Z_A']/r['Z_B']:7.1f}  {r['R_eq']:8.3f}"
              f"  {ref_str}  {err_str}  {r['D_e']:8.4f}  {bound:>5s}")

    return results


def section2_pes_shape(scan_results: list) -> dict:
    """Section 2: PES shape comparison for Z_A=2 and Z_A=6."""
    print("\n" + "=" * 70)
    print("SECTION 2: PES shape comparison (Z_A=2 vs Z_A=6)")
    print("=" * 70)

    shape_data = {}
    for Z_A_target in [2, 6]:
        res = next((r for r in scan_results if r['Z_A'] == Z_A_target), None)
        if res is None:
            print(f"  No data for Z_A={Z_A_target}")
            continue

        valid = [e for e in res['energies'] if e['E_total'] is not None]
        R_arr = np.array([e['R'] for e in valid])
        E_arr = np.array([e['E_total'] for e in valid])

        # Shift so minimum = 0
        E_shifted = E_arr - np.min(E_arr)
        R_eq = res['R_eq']

        # Curvature at minimum (finite difference)
        i_min = np.argmin(E_arr)
        if 0 < i_min < len(E_arr) - 1:
            dR = R_arr[i_min + 1] - R_arr[i_min]
            k = (E_arr[i_min - 1] - 2 * E_arr[i_min] + E_arr[i_min + 1]) / dR**2
        else:
            k = float('nan')

        shape_data[f'Z_A={Z_A_target}'] = {
            'R_eq': float(R_eq),
            'curvature_k': float(k),
            'D_e': float(res['D_e']),
            'E_min': float(res['E_min']),
        }

        print(f"\n  Z_A={Z_A_target}:")
        print(f"    R_eq = {R_eq:.3f} bohr")
        print(f"    Curvature k = {k:.4f} Ha/bohr^2")
        print(f"    D_e = {res['D_e']:.6f} Ha")
        print(f"    PES depth/curvature ratio = {res['D_e']/k:.4f} bohr^2"
              if k > 0 else "    Curvature invalid")

    if 'Z_A=2' in shape_data and 'Z_A=6' in shape_data:
        k2 = shape_data['Z_A=2']['curvature_k']
        k6 = shape_data['Z_A=6']['curvature_k']
        if k2 > 0 and k6 > 0:
            print(f"\n  Curvature ratio k(6)/k(2) = {k6/k2:.2f}")
            print(f"  Stiffer PES -> R_eq less sensitive to corrections"
                  if k6 > k2 else
                  f"  Flatter PES -> R_eq MORE sensitive to corrections")

    return shape_data


def section3_zeff_profile() -> dict:
    """Section 3: Z_eff(r) profile for oxygen (Z=8)."""
    print("\n" + "=" * 70)
    print("SECTION 3: Z_eff(r) profile for O (Z=8)")
    print("=" * 70)

    t0 = time.time()
    core = CoreScreening(Z=8, l_max=L_MAX, n_alpha=200)
    core.solve(verbose=False)
    t_core = time.time() - t0
    print(f"  Core solved in {t_core:.1f}s, E_core = {core.energy:.6f} Ha")

    # Evaluate Z_eff at key distances
    r_points = [0.0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 1.81, 2.0, 3.0, 5.0]
    r_arr = np.array(r_points)
    z_eff_arr = core.z_eff(r_arr)

    print(f"\n  {'r (bohr)':>10s}  {'Z_eff(r)':>10s}  {'N_screen':>10s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}")
    zeff_data = []
    for r, z in zip(r_points, z_eff_arr):
        n_screen = 8.0 - z
        print(f"  {r:10.3f}  {z:10.4f}  {n_screen:10.4f}")
        zeff_data.append({'r': r, 'Z_eff': float(z), 'N_screen': float(n_screen)})

    # Key diagnostic: Z_eff at R_eq
    z_at_req = float(core.z_eff(1.81))
    deviation = z_at_req - 6.0
    print(f"\n  At R_eq=1.81 bohr (expt O-H distance):")
    print(f"    Z_eff(1.81) = {z_at_req:.4f}")
    print(f"    Deviation from asymptotic (6.0): {deviation:+.4f}")
    print(f"    Relative deviation: {abs(deviation)/6.0*100:.1f}%")

    if abs(deviation) > 0.1:
        print(f"    ** SIGNIFICANT: Z_eff deviates from 6 by > 0.1 at R_eq **")
        print(f"    ** Using constant Z_A=6 overbinds the O side **")
    else:
        print(f"    Z_eff ~ 6 at R_eq -- constant Z_A approximation is valid")

    return {
        'E_core': float(core.energy),
        'z_eff_profile': zeff_data,
        'z_eff_at_Req': float(z_at_req),
        'deviation_from_6': float(deviation),
        'core_object': core,  # for reuse (not serialized)
    }


def _solve_single_R_with_zeff_func(
    R: float, Z_A_const: float, Z_B: float,
    Z_A_func=None,
) -> float:
    """
    Solve Level 4 at a single R with optional Z_eff(r) injection.

    Uses the same adiabatic approach as solve_level4_h2_multichannel
    but passes Z_A_func to build_angular_hamiltonian.
    """
    homonuclear = (Z_A_const == Z_B)
    channels = _channel_list(L_MAX, homonuclear=homonuclear)
    n_ch = len(channels)

    # Same R_e grid as the main solver
    R_e_min, R_e_max = 0.3, 15.0
    R_e_angular = np.concatenate([
        np.linspace(R_e_min, 1.0, 40),
        np.linspace(1.0, 3.0, 40),
        np.linspace(3.0, 6.0, 30),
        np.linspace(6.0, R_e_max, 20),
    ])
    R_e_angular = np.unique(R_e_angular)

    # Origin at charge center for heteronuclear
    z0 = R * (Z_A_const - Z_B) / (2.0 * (Z_A_const + Z_B))

    h = (np.pi / 2) / (N_ALPHA + 1)
    alpha = (np.arange(N_ALPHA) + 1) * h

    mu_vals = np.zeros(len(R_e_angular))
    for i, R_e in enumerate(R_e_angular):
        rho = R / (2.0 * R_e)

        H = build_angular_hamiltonian(
            alpha, rho, R_e, l_max=L_MAX, Z=Z_A_const,
            m_max=0, Z_A=Z_A_const, Z_B=Z_B,
            z0=z0,
            Z_A_func=Z_A_func,
        )
        evals, _ = eigh(H)
        mu_vals[i] = evals[0]

    U = (mu_vals + 15.0 / 8.0) / R_e_angular**2

    # Radial solve
    U_spline = CubicSpline(R_e_angular, U, extrapolate=True)
    h_Re = (R_e_max - R_e_min) / (N_RE + 1)
    R_e_radial = R_e_min + (np.arange(N_RE) + 1) * h_Re
    V_radial = U_spline(R_e_radial)

    diag = np.ones(N_RE) / h_Re**2 + V_radial
    off_diag = -0.5 * np.ones(N_RE - 1) / h_Re**2
    evals_rad, _ = eigh_tridiagonal(diag, off_diag,
                                     select='i', select_range=(0, 0))
    E_elec = evals_rad[0]
    E_total = E_elec + Z_A_const * Z_B / R
    return float(E_total)


def section4_zeff_injection(core: CoreScreening) -> dict:
    """Section 4: Bond pair with Z_eff injection at R=1.81."""
    print("\n" + "=" * 70)
    print("SECTION 4: Bond pair with Z_eff(r) injection at R=1.81")
    print("=" * 70)

    R = 1.81
    Z_A = 6.0
    Z_B = 1.0

    # Constant Z_A=6
    print(f"\n  Solving with constant Z_A={Z_A}...")
    t0 = time.time()
    E_const = _solve_single_R_with_zeff_func(R, Z_A, Z_B, Z_A_func=None)
    t1 = time.time()
    print(f"    E_total(const Z_A=6) = {E_const:.6f} Ha  ({t1-t0:.1f}s)")

    # Z_eff(r) injection
    print(f"  Solving with Z_eff(r) injection...")
    t0 = time.time()
    E_zeff = _solve_single_R_with_zeff_func(R, Z_A, Z_B, Z_A_func=core.z_eff)
    t1 = time.time()
    print(f"    E_total(Z_eff(r))    = {E_zeff:.6f} Ha  ({t1-t0:.1f}s)")

    dE = E_zeff - E_const
    print(f"\n    Delta E = {dE:+.6f} Ha ({dE*27.211:.1f} meV)")
    print(f"    Z_eff injection {'raises' if dE > 0 else 'lowers'} energy"
          f" by {abs(dE):.6f} Ha")

    if abs(dE) > 0.01:
        print(f"    ** SIGNIFICANT: |dE| > 10 mHa **")
        print(f"    ** The constant-Z_A approximation introduces "
              f"meaningful error **")
    else:
        print(f"    Constant Z_A=6 is a reasonable approximation at this R")

    # Also try at a few more R values for trend
    print(f"\n  R-dependent comparison:")
    print(f"    {'R':>6s}  {'E(const)':>12s}  {'E(Z_eff)':>12s}  {'dE':>10s}")
    print(f"    {'-'*6}  {'-'*12}  {'-'*12}  {'-'*10}")

    injection_data = []
    for R_test in [1.0, 1.5, 1.81, 2.5, 3.5]:
        E_c = _solve_single_R_with_zeff_func(R_test, Z_A, Z_B, Z_A_func=None)
        E_z = _solve_single_R_with_zeff_func(R_test, Z_A, Z_B,
                                              Z_A_func=core.z_eff)
        d = E_z - E_c
        print(f"    {R_test:6.2f}  {E_c:12.6f}  {E_z:12.6f}  {d:+10.6f}")
        injection_data.append({
            'R': R_test, 'E_const': E_c, 'E_zeff': E_z, 'dE': d,
        })

    return {
        'R': R,
        'E_const_ZA6': E_const,
        'E_zeff_injection': E_zeff,
        'dE': dE,
        'R_sweep': injection_data,
    }


def section5_pk_effect(core: CoreScreening) -> dict:
    """Section 5: Bond pair with PK pseudopotential for Z=8."""
    print("\n" + "=" * 70)
    print("SECTION 5: Bond pair with PK pseudopotential (Z=8)")
    print("=" * 70)

    pk = AbInitioPK(core, n_core=2)
    print(f"\n  {pk.summary()}")

    pk_potentials = [pk.pk_dict(atom='A')]
    R = 1.81
    Z_A = 6.0
    Z_B = 1.0

    # Without PK
    print(f"\n  Solving at R={R} WITHOUT PK...")
    result_nopk = solve_level4_h2_multichannel(
        R=R, Z_A=Z_A, Z_B=Z_B,
        l_max=L_MAX, n_alpha=N_ALPHA, n_Re=N_RE,
        verbose=False,
    )
    E_nopk = result_nopk['E_total']
    print(f"    E_total(no PK) = {E_nopk:.6f} Ha")

    # With PK
    print(f"  Solving at R={R} WITH PK...")
    result_pk = solve_level4_h2_multichannel(
        R=R, Z_A=Z_A, Z_B=Z_B,
        l_max=L_MAX, n_alpha=N_ALPHA, n_Re=N_RE,
        verbose=False,
        pk_potentials=pk_potentials,
    )
    E_pk = result_pk['E_total']
    print(f"    E_total(PK)    = {E_pk:.6f} Ha")

    dE_pk = E_pk - E_nopk
    print(f"\n    PK effect: {dE_pk:+.6f} Ha ({dE_pk*27.211:.1f} meV)")
    print(f"    PK A = {pk.A:.4f}, B = {pk.B:.4f}")

    # PES comparison: with and without PK
    print(f"\n  Short PES scan with/without PK:")
    print(f"    {'R':>6s}  {'E(no PK)':>12s}  {'E(PK)':>12s}  {'dE':>10s}")
    print(f"    {'-'*6}  {'-'*12}  {'-'*12}  {'-'*10}")

    pk_pes_data = []
    R_scan = np.arange(1.0, 4.1, 0.25)
    E_nopk_arr = []
    E_pk_arr = []

    for R_val in R_scan:
        r1 = solve_level4_h2_multichannel(
            R=R_val, Z_A=Z_A, Z_B=Z_B,
            l_max=L_MAX, n_alpha=N_ALPHA, n_Re=N_RE,
            verbose=False,
        )
        r2 = solve_level4_h2_multichannel(
            R=R_val, Z_A=Z_A, Z_B=Z_B,
            l_max=L_MAX, n_alpha=N_ALPHA, n_Re=N_RE,
            verbose=False,
            pk_potentials=pk_potentials,
        )
        e1 = r1['E_total']
        e2 = r2['E_total']
        d = e2 - e1
        print(f"    {R_val:6.2f}  {e1:12.6f}  {e2:12.6f}  {d:+10.6f}")
        E_nopk_arr.append(e1)
        E_pk_arr.append(e2)
        pk_pes_data.append({
            'R': float(R_val), 'E_nopk': float(e1),
            'E_pk': float(e2), 'dE': float(d),
        })

    E_nopk_arr = np.array(E_nopk_arr)
    E_pk_arr = np.array(E_pk_arr)
    R_eq_nopk = R_scan[np.argmin(E_nopk_arr)]
    R_eq_pk = R_scan[np.argmin(E_pk_arr)]

    print(f"\n    R_eq(no PK) = {R_eq_nopk:.2f}")
    print(f"    R_eq(PK)    = {R_eq_pk:.2f}")
    print(f"    PK shifts R_eq by {R_eq_pk - R_eq_nopk:+.2f} bohr")

    return {
        'pk_A': float(pk.A),
        'pk_B': float(pk.B),
        'pk_r_core': float(pk.r_core),
        'E_nopk_at_181': float(E_nopk),
        'E_pk_at_181': float(E_pk),
        'dE_pk_at_181': float(dE_pk),
        'R_eq_nopk': float(R_eq_nopk),
        'R_eq_pk': float(R_eq_pk),
        'pk_pes': pk_pes_data,
    }


def main():
    t_total = time.time()
    results = {}

    # Section 1: Charge asymmetry scan
    scan_results = section1_charge_asymmetry_scan()
    results['section1_scan'] = []
    for r in scan_results:
        # Remove non-serializable items
        entry = {k: v for k, v in r.items() if k != 'core_object'}
        results['section1_scan'].append(entry)

    # Section 2: PES shape
    results['section2_shape'] = section2_pes_shape(scan_results)

    # Section 3: Z_eff profile
    zeff_result = section3_zeff_profile()
    results['section3_zeff'] = {
        k: v for k, v in zeff_result.items() if k != 'core_object'
    }
    core = zeff_result['core_object']

    # Section 4: Z_eff injection
    results['section4_injection'] = section4_zeff_injection(core)

    # Section 5: PK effect
    results['section5_pk'] = section5_pk_effect(core)

    # Final report
    t_elapsed = time.time() - t_total
    results['total_time_s'] = t_elapsed

    print("\n" + "=" * 70)
    print("FINAL DIAGNOSTIC REPORT")
    print("=" * 70)

    # Summary table
    print(f"\n  {'Z_A':>4s}  {'Z_A/Z_B':>7s}  {'R_eq':>8s}  {'R_eq_ref':>8s}"
          f"  {'Err%':>6s}  {'D_e':>8s}")
    print(f"  {'-'*4}  {'-'*7}  {'-'*8}  {'-'*8}  {'-'*6}  {'-'*8}")
    for r in scan_results:
        ref_str = (f"{r['R_eq_ref']:8.3f}" if 'R_eq_ref' in r
                   else '     ---')
        err_str = (f"{r['R_eq_error_pct']:6.1f}" if 'R_eq_error_pct' in r
                   else '   ---')
        print(f"  {r['Z_A']:4.0f}  {r['Z_A']/r['Z_B']:7.1f}"
              f"  {r['R_eq']:8.3f}  {ref_str}  {err_str}  {r['D_e']:8.4f}")

    # Key findings
    print(f"\n  Key findings:")

    z_at_req = results['section3_zeff']['z_eff_at_Req']
    dev = results['section3_zeff']['deviation_from_6']
    print(f"    Z_eff(r=1.81) = {z_at_req:.4f}  (deviation from 6: {dev:+.4f})")

    dE_inj = results['section4_injection']['dE']
    print(f"    Z_eff injection effect at R=1.81: dE = {dE_inj:+.6f} Ha")

    dE_pk = results['section5_pk']['dE_pk_at_181']
    print(f"    PK effect at R=1.81: dE = {dE_pk:+.6f} Ha")

    R_eq_nopk = results['section5_pk']['R_eq_nopk']
    R_eq_pk = results['section5_pk']['R_eq_pk']
    print(f"    R_eq shift from PK: {R_eq_nopk:.2f} -> {R_eq_pk:.2f}")

    shape = results['section2_shape']
    if 'Z_A=2' in shape and 'Z_A=6' in shape:
        k2 = shape['Z_A=2']['curvature_k']
        k6 = shape['Z_A=6']['curvature_k']
        if k2 > 0 and k6 > 0:
            print(f"    PES curvature ratio k(6)/k(2) = {k6/k2:.2f}"
                  f" ({'stiffer' if k6 > k2 else 'flatter'} at high Z_A)")

    print(f"\n  Total diagnostic time: {t_elapsed:.0f}s")

    # Save results
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'data', 'h2o_bond_diagnostic.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved to {outpath}")


if __name__ == '__main__':
    main()
