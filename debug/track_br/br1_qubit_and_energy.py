"""
Track BR-1: Qubit metrics at varying max_n + single-point energy at R=3.015.

Part A: Qubit Hamiltonian metrics for LiH at max_n = 1, 2, 3
Part B: Composed total energy at R=3.015 bohr for l_max = 2, 3, 4

Output: debug/track_br/br1_results.json
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Part A: Qubit metrics
# ---------------------------------------------------------------------------

def run_qubit_metrics():
    """Compute qubit Hamiltonian metrics for LiH at max_n = 1, 2, 3."""
    from geovac.composed_qubit import build_composed_hamiltonian, lih_spec, GAUSSIAN_LIH_PUBLISHED
    from geovac.trotter_bounds import pauli_1norm
    from geovac.measurement_grouping import qwc_groups

    results = []
    for max_n in [1, 2, 3]:
        print(f"\n{'='*60}")
        print(f"Part A: LiH qubit metrics at max_n = {max_n}")
        print(f"{'='*60}")

        spec = lih_spec(max_n_core=max_n, max_n_val=max_n)

        # Full Hamiltonian (with PK)
        res_full = build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=True)
        norm_full = pauli_1norm(res_full['qubit_op'])
        groups_full = qwc_groups(res_full['qubit_op'])
        n_qwc_full = len(groups_full)

        # Electronic-only (without PK)
        res_elec = build_composed_hamiltonian(spec, pk_in_hamiltonian=False, verbose=False)
        norm_elec = pauli_1norm(res_elec['qubit_op'])
        groups_elec = qwc_groups(res_elec['qubit_op'])
        n_qwc_elec = len(groups_elec)

        entry = {
            'max_n': max_n,
            'Q': res_full['Q'],
            'M': res_full['M'],
            'N_pauli_full': res_full['N_pauli'],
            'N_pauli_elec': res_elec['N_pauli'],
            '1_norm_full_Ha': round(norm_full, 4),
            '1_norm_elec_Ha': round(norm_elec, 4),
            'n_qwc_full': n_qwc_full,
            'n_qwc_elec': n_qwc_elec,
            'nuclear_repulsion': round(res_full['nuclear_repulsion'], 6),
            'wall_time_s': round(res_full['wall_time_s'], 2),
        }
        results.append(entry)

        print(f"  Q = {entry['Q']}, M = {entry['M']}")
        print(f"  N_Pauli (full) = {entry['N_pauli_full']}")
        print(f"  N_Pauli (elec) = {entry['N_pauli_elec']}")
        print(f"  1-norm  (full) = {entry['1_norm_full_Ha']:.4f} Ha")
        print(f"  1-norm  (elec) = {entry['1_norm_elec_Ha']:.4f} Ha")
        print(f"  QWC groups (full) = {entry['n_qwc_full']}")
        print(f"  QWC groups (elec) = {entry['n_qwc_elec']}")

    # Gaussian baselines
    print(f"\n{'='*60}")
    print("Gaussian baselines (Trenev et al.):")
    for basis, data in GAUSSIAN_LIH_PUBLISHED.items():
        print(f"  {basis}: Q={data['Q']}, N_Pauli={data['N_pauli']}")

    return results


# ---------------------------------------------------------------------------
# Part B: Single-point energy at R=3.015 bohr
# ---------------------------------------------------------------------------

def run_single_point_energy():
    """Compute composed total energy at R=3.015 for l_max = 2, 3, 4."""
    from geovac.composed_diatomic import ComposedDiatomicSolver

    R_target = 3.015
    R_grid = np.array([R_target])

    results = []
    for l_max in [2, 3, 4]:
        print(f"\n{'='*60}")
        print(f"Part B: LiH composed energy at R={R_target}, l_max={l_max}")
        print(f"{'='*60}")

        t0 = time.time()
        solver = ComposedDiatomicSolver.LiH_ab_initio(
            l_max=l_max,
            level4_method='variational_2d',
            pk_channel_mode='l_dependent',
            verbose=True,
        )
        solver.solve_core()
        E_core = solver.E_core
        pk_A = solver.pk_A
        pk_B = solver.pk_B
        print(f"  E_core = {E_core:.6f} Ha, PK: A={pk_A:.4f}, B={pk_B:.4f}")

        # Single-point solve (bypass scan_pes which needs >= 3 points)
        from geovac.composed_diatomic import _v_cross_nuc_1s
        E_elec = solver._solve_valence_at_R(R_target, n_Re=200)
        V_NN = solver.Z_A_bare * solver.Z_B / R_target
        V_cross = _v_cross_nuc_1s(solver.Z_A_bare, solver.n_core, solver.Z_B, R_target)
        E_composed = E_core + V_cross + E_elec + V_NN
        wall = time.time() - t0

        # Reference: exact nonrelativistic LiH at R=3.015
        E_exact = -8.0706  # Ha (Cencek & Rychlewski, approximate)
        error_ha = E_composed - E_exact
        error_pct = (E_composed - E_exact) / abs(E_exact) * 100

        entry = {
            'l_max': l_max,
            'R': R_target,
            'E_composed': round(E_composed, 6),
            'E_elec': round(E_elec, 6),
            'E_core': round(E_core, 6),
            'V_NN': round(V_NN, 6),
            'V_cross': round(V_cross, 6),
            'E_exact_ref': E_exact,
            'error_Ha': round(error_ha, 6),
            'error_pct': round(error_pct, 3),
            'pk_A': round(pk_A, 4),
            'pk_B': round(pk_B, 4),
            'wall_time_s': round(wall, 1),
        }
        results.append(entry)

        print(f"\n  E_composed = {E_composed:.6f} Ha")
        print(f"  E_exact    = {E_exact:.4f} Ha")
        print(f"  Error      = {error_ha:+.6f} Ha ({error_pct:+.3f}%)")
        print(f"  Wall time  = {wall:.1f}s")

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    out_dir = Path(__file__).parent
    out_file = out_dir / 'br1_results.json'

    print("=" * 70)
    print("Track BR-1: LiH Qubit Metrics + Single-Point Energy")
    print("=" * 70)

    # Part A
    qubit_results = run_qubit_metrics()

    # Part B
    energy_results = run_single_point_energy()

    # Summary table
    print(f"\n{'='*70}")
    print("SUMMARY: Part A — Qubit Metrics")
    print(f"{'='*70}")
    print(f"{'max_n':>6} {'Q':>4} {'N_Pauli':>8} {'1-norm':>10} {'QWC':>5}")
    for r in qubit_results:
        print(f"{r['max_n']:>6} {r['Q']:>4} {r['N_pauli_full']:>8} "
              f"{r['1_norm_full_Ha']:>10.2f} {r['n_qwc_full']:>5}")

    print(f"\n{'='*70}")
    print("SUMMARY: Part B — Single-Point Energy at R=3.015")
    print(f"{'='*70}")
    print(f"{'l_max':>6} {'E_composed':>12} {'Error(Ha)':>10} {'Error(%)':>10}")
    for r in energy_results:
        print(f"{r['l_max']:>6} {r['E_composed']:>12.6f} "
              f"{r['error_Ha']:>+10.6f} {r['error_pct']:>+10.3f}%")

    # Check convergence direction
    if len(energy_results) >= 2:
        energies = [r['E_composed'] for r in energy_results]
        diffs = [energies[i+1] - energies[i] for i in range(len(energies)-1)]
        monotonic_down = all(d < 0 for d in diffs)
        monotonic_up = all(d > 0 for d in diffs)
        print(f"\n  Convergence: {'monotonic decreasing' if monotonic_down else 'monotonic increasing' if monotonic_up else 'non-monotonic'}")
        for i, d in enumerate(diffs):
            print(f"  l_max {energy_results[i]['l_max']}->{energy_results[i+1]['l_max']}: dE = {d:+.6f} Ha")

    # Save
    output = {
        'description': 'Track BR-1: LiH qubit metrics + single-point energy convergence',
        'qubit_metrics': qubit_results,
        'energy_convergence': energy_results,
    }
    with open(out_file, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_file}")
