"""
Systematic Stress Testing of the Prolate Spheroidal Lattice
============================================================
Phase 1: Excited states, limiting cases, grid convergence, asymmetric charges.
Phase 2: V_ee integrals and H2 CI.

Date: 2026-03-13
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import time
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_spheroidal_lattice import (
    ProlateSpheroidalLattice,
    scan_h2plus_pes,
    fit_spectroscopic_constants,
)


# ============================================================
# Reference data
# ============================================================
# H2+ exact: Bates, Ledsham & Stewart (1953)
H2PLUS_REF = {
    '1sigma_g': {
        'R_eq': 1.997,
        'E_min': -0.6026,  # Total energy at equilibrium
        'D_e': 0.1026,     # Dissociation energy (relative to H + H+)
        'dissoc_limit': -0.500,  # H(1s) + H+
    },
    '1sigma_u': {
        # Purely repulsive (no minimum in total energy)
        'dissoc_limit': -0.500,  # H(1s) + H+
        'E_R2': -0.1675,        # approximate E_total at R=2.0
    },
    '1pi_u': {
        # Shallow minimum around R=8 bohr
        'dissoc_limit': -0.125,  # H(2p) + H+
    },
}

# HeH2+ (one electron, Z_A=2, Z_B=1)
# United atom limit: Li2+ (Z=3) -> E = -Z^2/2 = -4.5 Ha
HEH2PLUS_REF = {
    'dissoc_limit': -2.000,  # He+(1s) + H+ -> E = -Z^2/2 = -2.0 Ha
}


def run_test(name: str, func, *args, **kwargs):
    """Run a test and return (status, result_dict)."""
    t0 = time.time()
    try:
        result = func(*args, **kwargs)
        dt = time.time() - t0
        result['time'] = dt
        return result
    except Exception as e:
        dt = time.time() - t0
        return {'status': 'FAIL', 'error': str(e), 'time': dt}


# ============================================================
# Category 1: Excited States of H2+
# ============================================================
def test_sigma_g_baseline(N_xi: int = 5000) -> dict:
    """1sigma_g ground state at R=2.0 (baseline)."""
    lat = ProlateSpheroidalLattice(R=2.0, N_xi=N_xi)
    E_tot = lat.total_energy()
    E_exact = -0.6026
    err = abs(E_tot - E_exact) / abs(E_exact) * 100
    return {
        'status': 'PASS' if err < 2.0 else 'FAIL',
        'E_total': E_tot,
        'E_exact': E_exact,
        'error_pct': err,
    }


def test_sigma_u(N_xi: int = 5000) -> dict:
    """1sigma_u first excited state at R=2.0."""
    lat = ProlateSpheroidalLattice(R=2.0, N_xi=N_xi, n_angular=1)
    E_tot = lat.total_energy()
    E_ref = -0.1675
    err = abs(E_tot - E_ref) / abs(E_ref) * 100
    # Also check it's above 1sigma_g
    lat_g = ProlateSpheroidalLattice(R=2.0, N_xi=N_xi)
    E_g = lat_g.total_energy()
    ordering_ok = E_tot > E_g
    return {
        'status': 'PASS' if err < 5.0 and ordering_ok else 'FAIL',
        'E_total': E_tot,
        'E_ref': E_ref,
        'error_pct': err,
        'ordering_ok': ordering_ok,
        'E_ground': E_g,
    }


def test_pi_u(N_xi: int = 5000) -> dict:
    """1pi_u state: m=1, shallow minimum around R=8."""
    # Check dissociation limit
    lat_R15 = ProlateSpheroidalLattice(R=15.0, N_xi=N_xi, xi_max=60.0, m=1)
    E_R15 = lat_R15.total_energy()
    dissoc_err = abs(E_R15 - (-0.125)) / 0.125 * 100

    # Find approximate minimum
    R_vals = np.array([4.0, 6.0, 8.0, 10.0, 12.0])
    pes = scan_h2plus_pes(R_vals, N_xi=N_xi, m=1, verbose=False)
    idx_min = np.argmin(pes['E_total'])
    has_min = 0 < idx_min < len(R_vals) - 1
    E_min = pes['E_total'][idx_min]
    R_min = R_vals[idx_min]

    return {
        'status': 'PASS' if dissoc_err < 5.0 else 'MARGINAL',
        'E_R15': E_R15,
        'dissoc_limit': -0.125,
        'dissoc_err_pct': dissoc_err,
        'has_minimum': has_min,
        'R_min_approx': R_min,
        'E_min_approx': E_min,
    }


def test_state_ordering(N_xi: int = 5000) -> dict:
    """Verify correct energy ordering of states at R=2.0."""
    states = {}
    configs = [
        ('1sigma_g', {'n_angular': 0, 'm': 0}),
        ('1sigma_u', {'n_angular': 1, 'm': 0}),
        ('1pi_u', {'n_angular': 0, 'm': 1}),
    ]
    for name, kw in configs:
        lat = ProlateSpheroidalLattice(R=2.0, N_xi=N_xi, **kw)
        states[name] = lat.total_energy()

    # Ground state should be lowest
    correct = states['1sigma_g'] < states['1sigma_u']
    return {
        'status': 'PASS' if correct else 'FAIL',
        'energies': states,
        'ordering_correct': correct,
    }


# ============================================================
# Category 2: Limiting Cases
# ============================================================
def test_united_atom_limit(N_xi: int = 5000) -> dict:
    """As R -> 0, H2+ -> He+ (E_elec -> -2.0 Ha)."""
    results = {}
    for R in [0.2, 0.5, 1.0]:
        lat = ProlateSpheroidalLattice(R=R, N_xi=N_xi)
        E_el, c2, A = lat.solve()
        E_tot = E_el + 1.0 / R
        results[R] = {'E_elec': E_el, 'E_total': E_tot, 'c2': c2}

    # At R=0.2, E_elec should approach -2.0 (He+)
    E_el_small = results[0.2]['E_elec']
    err = abs(E_el_small - (-2.0)) / 2.0 * 100
    # Total energy diverges (V_NN -> inf), but E_elec should approach united atom
    return {
        'status': 'PASS' if err < 15 else 'MARGINAL',
        'results': results,
        'E_elec_R02': E_el_small,
        'He_plus_limit': -2.0,
        'error_pct': err,
    }


def test_dissociation_limit(N_xi: int = 8000) -> dict:
    """As R -> inf, E_total -> -0.5 Ha (H + H+)."""
    results = {}
    for R in [10.0, 20.0, 50.0]:
        lat = ProlateSpheroidalLattice(
            R=R, N_xi=N_xi, xi_max=max(25.0, R * 3)
        )
        E_tot = lat.total_energy()
        results[R] = E_tot

    err_50 = abs(results[50.0] - (-0.5)) / 0.5 * 100
    return {
        'status': 'PASS' if err_50 < 5.0 else 'MARGINAL',
        'results': results,
        'limit': -0.500,
        'error_R50_pct': err_50,
    }


def test_sigma_u_dissociation(N_xi: int = 8000) -> dict:
    """sigma_u dissociation: E -> -0.5 Ha (same limit as sigma_g)."""
    results = {}
    for R in [10.0, 20.0]:
        lat = ProlateSpheroidalLattice(
            R=R, N_xi=N_xi, xi_max=max(25.0, R * 3),
            n_angular=1,
        )
        E_tot = lat.total_energy()
        results[R] = E_tot

    err_20 = abs(results[20.0] - (-0.5)) / 0.5 * 100
    return {
        'status': 'PASS' if err_20 < 15.0 else 'MARGINAL',
        'results': results,
        'limit': -0.500,
        'error_R20_pct': err_20,
    }


# ============================================================
# Category 3: Asymmetric Charges
# ============================================================
def test_heh2plus(N_xi: int = 5000) -> dict:
    """HeH2+ (Z_A=2, Z_B=1): one-electron heteronuclear."""
    # PES scan
    R_vals = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0])
    pes = scan_h2plus_pes(R_vals, Z_A=2, Z_B=1, N_xi=N_xi, verbose=False)
    fit = fit_spectroscopic_constants(pes['R'], pes['E_total'])

    # Dissociation limit: He+(1s) + H+ -> -2.0 Ha
    lat_R10 = ProlateSpheroidalLattice(
        R=10.0, Z_A=2, Z_B=1, N_xi=N_xi, xi_max=40.0
    )
    E_R10 = lat_R10.total_energy()

    return {
        'status': 'PASS' if fit['E_min'] < -2.0 else 'FAIL',
        'R_eq': fit['R_eq'],
        'E_min': fit['E_min'],
        'D_e': -2.0 - fit['E_min'],  # Relative to He+ + H+
        'E_R10': E_R10,
        'dissoc_limit': -2.0,
        'boundary': fit['boundary'],
    }


def test_extreme_asymmetry(N_xi: int = 5000) -> dict:
    """LiHe3+ (Z_A=3, Z_B=2) and CH5+ (Z_A=6, Z_B=1)."""
    results = {}
    for name, Z_A, Z_B in [('LiHe3+', 3, 2), ('CH5+', 6, 1)]:
        try:
            lat = ProlateSpheroidalLattice(
                R=2.0, Z_A=Z_A, Z_B=Z_B, N_xi=N_xi
            )
            E_tot = lat.total_energy()
            # United atom limit: -(Z_A+Z_B)^2/2
            E_united = -(Z_A + Z_B)**2 / 2.0
            # Dissociation limit: -max(Z_A,Z_B)^2/2
            E_dissoc = -max(Z_A, Z_B)**2 / 2.0
            results[name] = {
                'E_total': E_tot,
                'E_united_limit': E_united,
                'E_dissoc_limit': E_dissoc,
                'bound': E_tot < E_dissoc,
            }
        except Exception as e:
            results[name] = {'error': str(e)}

    all_ok = all('E_total' in v for v in results.values())
    return {
        'status': 'PASS' if all_ok else 'FAIL',
        'results': results,
    }


# ============================================================
# Category 4: Grid Convergence & Scaling
# ============================================================
def test_grid_convergence(R: float = 2.0) -> dict:
    """Grid convergence study at fixed R."""
    results = []
    for N_xi in [200, 500, 1000, 2000, 5000, 10000]:
        t0 = time.time()
        lat = ProlateSpheroidalLattice(R=R, N_xi=N_xi)
        E_tot = lat.total_energy()
        dt = time.time() - t0
        results.append({
            'N_xi': N_xi,
            'E_total': E_tot,
            'time': dt,
        })

    # Richardson extrapolation from two finest
    E_5k = results[-2]['E_total']
    E_10k = results[-1]['E_total']
    E_extrap = (4 * E_10k - E_5k) / 3  # Assuming O(h^2)
    err_5k = abs(E_5k - E_extrap) / abs(E_extrap) * 100
    err_10k = abs(E_10k - E_extrap) / abs(E_extrap) * 100

    return {
        'status': 'PASS',
        'results': results,
        'E_extrapolated': E_extrap,
        'err_5k_pct': err_5k,
        'err_10k_pct': err_10k,
    }


def test_z_scaling() -> dict:
    """How does required N_xi scale with Z?"""
    results = {}
    target_err = 1.0  # Target < 1% error
    for name, Z_A, Z_B, R in [
        ('H2+', 1, 1, 2.0),
        ('HeH2+', 2, 1, 1.5),
        ('LiH3+', 3, 1, 2.0),
    ]:
        # Find minimum N_xi for < 1% error (use N=10000 as reference)
        lat_ref = ProlateSpheroidalLattice(
            R=R, Z_A=Z_A, Z_B=Z_B, N_xi=10000
        )
        E_ref = lat_ref.total_energy()

        min_N = None
        for N_xi in [500, 1000, 2000, 3000, 5000, 8000]:
            lat = ProlateSpheroidalLattice(
                R=R, Z_A=Z_A, Z_B=Z_B, N_xi=N_xi
            )
            E = lat.total_energy()
            err = abs(E - E_ref) / abs(E_ref) * 100
            if err < target_err and min_N is None:
                min_N = N_xi

        results[name] = {
            'Z_A': Z_A, 'Z_B': Z_B,
            'E_ref': E_ref,
            'min_N_for_1pct': min_N,
        }

    return {
        'status': 'PASS',
        'results': results,
    }


# ============================================================
# Category 5: Spectroscopic Constants (Full PES)
# ============================================================
def test_h2plus_pes(N_xi: int = 8000) -> dict:
    """Full H2+ PES with spectroscopic constants."""
    R_vals = np.arange(0.5, 8.01, 0.25)
    pes = scan_h2plus_pes(R_vals, N_xi=N_xi, verbose=False)
    fit = fit_spectroscopic_constants(pes['R'], pes['E_total'])

    R_err = abs(fit['R_eq'] - 1.997) / 1.997 * 100
    E_err = abs(fit['E_min'] - (-0.6026)) / 0.6026 * 100
    D_err = abs(fit['D_e'] - 0.1026) / 0.1026 * 100

    return {
        'status': 'PASS' if R_err < 1.0 and E_err < 1.0 else 'MARGINAL',
        'R_eq': fit['R_eq'],
        'E_min': fit['E_min'],
        'D_e': fit['D_e'],
        'k': fit['k'],
        'R_err_pct': R_err,
        'E_err_pct': E_err,
        'D_err_pct': D_err,
    }


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("  PROLATE SPHEROIDAL LATTICE — SYSTEMATIC STRESS TESTS")
    print("  Date: 2026-03-13")
    print("=" * 70)

    all_results = {}

    # --- Category 1: Excited States ---
    print("\n--- Category 1: Excited States of H2+ ---")

    print("\n[1.1] 1sigma_g baseline...")
    r = run_test('1sigma_g', test_sigma_g_baseline)
    all_results['1sigma_g'] = r
    print(f"  {r['status']}  E_total={r.get('E_total', 'N/A'):.6f}  "
          f"err={r.get('error_pct', 'N/A'):.2f}%  [{r['time']:.1f}s]")

    print("\n[1.2] 1sigma_u (antibonding)...")
    r = run_test('1sigma_u', test_sigma_u)
    all_results['1sigma_u'] = r
    print(f"  {r['status']}  E_total={r.get('E_total', 'N/A'):.6f}  "
          f"err={r.get('error_pct', 'N/A'):.2f}%  [{r['time']:.1f}s]")
    if 'ordering_ok' in r:
        print(f"  Ordering: 1sigma_g ({r.get('E_ground', 0):.4f}) < "
              f"1sigma_u ({r.get('E_total', 0):.4f}) = {r['ordering_ok']}")

    print("\n[1.3] 1pi_u (m=1)...")
    r = run_test('1pi_u', test_pi_u)
    all_results['1pi_u'] = r
    print(f"  {r['status']}  E_R15={r.get('E_R15', 'N/A'):.6f}  "
          f"dissoc_err={r.get('dissoc_err_pct', 'N/A'):.2f}%  [{r['time']:.1f}s]")
    if 'has_minimum' in r:
        print(f"  Has minimum: {r['has_minimum']}  "
              f"R_min~{r.get('R_min_approx', 'N/A'):.1f}  "
              f"E_min~{r.get('E_min_approx', 'N/A'):.4f}")

    print("\n[1.4] State ordering check...")
    r = run_test('ordering', test_state_ordering)
    all_results['ordering'] = r
    print(f"  {r['status']}  [{r['time']:.1f}s]")
    if 'energies' in r:
        for st, E in r['energies'].items():
            print(f"    {st:12s}: {E:.6f} Ha")

    # --- Category 2: Limiting Cases ---
    print("\n--- Category 2: Limiting Cases ---")

    print("\n[2.1] United atom limit (R->0)...")
    r = run_test('united_atom', test_united_atom_limit)
    all_results['united_atom'] = r
    print(f"  {r['status']}  E_elec(R=0.2)={r.get('E_elec_R02', 'N/A'):.4f}  "
          f"(He+ limit: -2.0)  err={r.get('error_pct', 'N/A'):.1f}%  [{r['time']:.1f}s]")
    if 'results' in r:
        for R_val, data in r['results'].items():
            if isinstance(data, dict) and 'E_elec' in data:
                print(f"    R={R_val:.1f}: E_elec={data['E_elec']:.4f}  "
                      f"E_total={data['E_total']:.4f}")

    print("\n[2.2] Dissociation limit (R->inf, sigma_g)...")
    r = run_test('dissociation', test_dissociation_limit)
    all_results['dissociation'] = r
    print(f"  {r['status']}  err(R=50)={r.get('error_R50_pct', 'N/A'):.2f}%  [{r['time']:.1f}s]")
    if 'results' in r:
        for R_val, E in r['results'].items():
            print(f"    R={R_val:.0f}: E_total={E:.6f}")

    print("\n[2.3] sigma_u dissociation limit...")
    r = run_test('dissoc_sigma_u', test_sigma_u_dissociation)
    all_results['dissoc_sigma_u'] = r
    print(f"  {r['status']}  err(R=20)={r.get('error_R20_pct', 'N/A'):.2f}%  [{r['time']:.1f}s]")
    if 'results' in r:
        for R_val, E in r['results'].items():
            print(f"    R={R_val:.0f}: E_total={E:.6f}")

    # --- Category 3: Asymmetric Charges ---
    print("\n--- Category 3: Asymmetric Charges ---")

    print("\n[3.1] HeH2+ (Z_A=2, Z_B=1)...")
    r = run_test('HeH2+', test_heh2plus)
    all_results['HeH2+'] = r
    print(f"  {r['status']}  R_eq={r.get('R_eq', 'N/A'):.3f}  "
          f"E_min={r.get('E_min', 'N/A'):.4f}  D_e={r.get('D_e', 'N/A'):.4f}  [{r['time']:.1f}s]")

    print("\n[3.2] Extreme asymmetry (LiHe3+, CH5+)...")
    r = run_test('extreme', test_extreme_asymmetry)
    all_results['extreme'] = r
    print(f"  {r['status']}  [{r['time']:.1f}s]")
    if 'results' in r:
        for name, data in r['results'].items():
            if 'E_total' in data:
                print(f"    {name}: E_total={data['E_total']:.4f}  "
                      f"bound={data['bound']}")
            else:
                print(f"    {name}: ERROR {data.get('error', 'unknown')}")

    # --- Category 4: Grid Convergence ---
    print("\n--- Category 4: Grid Convergence & Scaling ---")

    print("\n[4.1] Grid convergence at R=2.0...")
    r = run_test('grid_conv', test_grid_convergence)
    all_results['grid_conv'] = r
    print(f"  {r['status']}  [{r['time']:.1f}s]")
    if 'results' in r:
        print(f"  {'N_xi':>8s} {'E_total':>12s} {'time':>8s}")
        for pt in r['results']:
            print(f"  {pt['N_xi']:8d} {pt['E_total']:12.8f} {pt['time']:8.2f}s")
        print(f"  E_extrapolated = {r.get('E_extrapolated', 'N/A'):.8f}")
        print(f"  err(N=5k) = {r.get('err_5k_pct', 'N/A'):.4f}%")
        print(f"  err(N=10k) = {r.get('err_10k_pct', 'N/A'):.4f}%")

    print("\n[4.2] Z-scaling: minimum grid for 1% accuracy...")
    r = run_test('z_scaling', test_z_scaling)
    all_results['z_scaling'] = r
    print(f"  {r['status']}  [{r['time']:.1f}s]")
    if 'results' in r:
        for name, data in r['results'].items():
            print(f"    {name}: E_ref={data['E_ref']:.4f}  "
                  f"min_N={data.get('min_N_for_1pct', 'N/A')}")

    # --- Category 5: Full PES ---
    print("\n--- Category 5: Spectroscopic Constants ---")

    print("\n[5.1] H2+ full PES (N_xi=8000)...")
    r = run_test('h2plus_pes', test_h2plus_pes)
    all_results['h2plus_pes'] = r
    print(f"  {r['status']}  [{r['time']:.1f}s]")
    if 'R_eq' in r:
        print(f"    R_eq = {r['R_eq']:.4f} bohr  (exact 1.997, err {r.get('R_err_pct', 'N/A'):.2f}%)")
        print(f"    E_min = {r['E_min']:.6f} Ha  (exact -0.6026, err {r.get('E_err_pct', 'N/A'):.2f}%)")
        print(f"    D_e = {r['D_e']:.6f} Ha  (exact 0.1026, err {r.get('D_err_pct', 'N/A'):.2f}%)")

    # ============================================================
    # Summary Table
    # ============================================================
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"  {'Test':<25s} {'Status':<10s} {'Error':<15s} {'Time':>8s}")
    print("-" * 65)

    summary_rows = [
        ('1sigma_g baseline', '1sigma_g', 'error_pct'),
        ('1sigma_u antibonding', '1sigma_u', 'error_pct'),
        ('1pi_u (m=1)', '1pi_u', 'dissoc_err_pct'),
        ('State ordering', 'ordering', None),
        ('United atom (R->0)', 'united_atom', 'error_pct'),
        ('Dissociation (R->inf)', 'dissociation', 'error_R50_pct'),
        ('sigma_u dissociation', 'dissoc_sigma_u', 'error_R20_pct'),
        ('HeH2+', 'HeH2+', None),
        ('Extreme asymmetry', 'extreme', None),
        ('Grid convergence', 'grid_conv', 'err_10k_pct'),
        ('Z-scaling', 'z_scaling', None),
        ('H2+ PES', 'h2plus_pes', 'R_err_pct'),
    ]

    n_pass = n_marginal = n_fail = 0
    for label, key, err_key in summary_rows:
        r = all_results.get(key, {})
        status = r.get('status', '???')
        err_str = ''
        if err_key and err_key in r:
            err_str = f"{r[err_key]:.2f}%"
        t = r.get('time', 0)

        if status == 'PASS':
            n_pass += 1
        elif status == 'MARGINAL':
            n_marginal += 1
        else:
            n_fail += 1

        print(f"  {label:<25s} {status:<10s} {err_str:<15s} {t:7.1f}s")

    print("-" * 65)
    print(f"  PASS: {n_pass}  MARGINAL: {n_marginal}  FAIL: {n_fail}")
    total_time = sum(r.get('time', 0) for r in all_results.values())
    print(f"  Total time: {total_time:.0f}s")

    # Save results
    outpath = os.path.join(os.path.dirname(__file__), 'data', 'prolate_stress_tests.txt')
    with open(outpath, 'w') as f:
        f.write("# Prolate Spheroidal Lattice Stress Tests\n")
        f.write(f"# Date: 2026-03-13\n")
        f.write(f"# PASS: {n_pass}  MARGINAL: {n_marginal}  FAIL: {n_fail}\n\n")
        for label, key, err_key in summary_rows:
            r = all_results.get(key, {})
            f.write(f"{label}: {r.get('status', '???')}")
            if err_key and err_key in r:
                f.write(f"  err={r[err_key]:.4f}%")
            f.write(f"  time={r.get('time', 0):.1f}s\n")
            # Write detail
            for k, v in r.items():
                if k not in ('status', 'time', 'results', 'error'):
                    f.write(f"  {k} = {v}\n")
            f.write("\n")
    print(f"\n  Results saved to {outpath}")

    return all_results


if __name__ == '__main__':
    results = main()
