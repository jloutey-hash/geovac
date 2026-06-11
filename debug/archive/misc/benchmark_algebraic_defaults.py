"""
Benchmark algebraic vs quadrature defaults for Z_eff and exchange.

Compares single-point energies at R_eq for LiH, BeH₂, and H₂O
using both old (spline/numerical) and new (spectral_laguerre/algebraic_laguerre)
defaults.
"""

import time
import numpy as np
import json
from pathlib import Path

# ======================================================================
# LiH comparison
# ======================================================================

def benchmark_lih():
    """Compare LiH single-point energies at several R values."""
    from geovac.core_screening import CoreScreening
    from geovac.ab_initio_pk import AbInitioPK
    from geovac.level4_multichannel import solve_level4_h2_multichannel
    from geovac.composed_diatomic import ComposedDiatomicSolver, _v_cross_nuc_1s

    results = {}

    # Test at R_eq and a few other R values
    R_values = [2.5, 3.015, 4.0, 5.0]

    for zeff_method in ['spline', 'spectral_laguerre']:
        print(f"\n{'='*60}")
        print(f"LiH: zeff_method='{zeff_method}'")
        print(f"{'='*60}")

        t0 = time.time()

        # Solve core with this zeff_method
        core = CoreScreening(Z=3, l_max=2, n_alpha=200,
                            zeff_method=zeff_method)
        core.solve(verbose=True)
        E_core = core.energy
        core_time = time.time() - t0

        # Derive PK
        pk = AbInitioPK(core, n_core=2)
        pk_d = pk.pk_dict(atom='A')
        pk_d['channel_mode'] = 'l_dependent'
        pk_potentials = [pk_d]

        # Compare Z_eff values at key r-points
        r_test = np.array([0.5, 1.0, 2.0, 3.0, 5.0])
        z_eff_vals = core.z_eff(r_test)
        print(f"  Z_eff at r={r_test}: {z_eff_vals}")

        # Scan energies
        for R in R_values:
            t1 = time.time()
            result = solve_level4_h2_multichannel(
                R=R, Z_A=1.0, Z_B=1.0, l_max=2,
                n_alpha=100, n_Re=300, verbose=False,
                pk_potentials=pk_potentials,
            )
            E_elec = result['E_elec']
            V_NN = 3.0 * 1.0 / R
            V_cross = _v_cross_nuc_1s(3.0, 2, 1.0, R)
            E_composed = E_core + V_cross + E_elec + V_NN
            dt = time.time() - t1

            key = f"LiH_{zeff_method}_R{R:.3f}"
            results[key] = {
                'E_core': E_core,
                'E_elec': E_elec,
                'E_composed': E_composed,
                'time': dt,
                'z_eff_method': zeff_method,
            }
            print(f"  R={R:.3f}: E_composed={E_composed:.6f} Ha ({dt:.1f}s)")

        results[f"LiH_{zeff_method}_core_time"] = core_time

    return results


def benchmark_lih_pes():
    """Full LiH PES comparison — old vs new defaults."""
    from geovac.composed_diatomic import ComposedDiatomicSolver

    results = {}

    for label, zeff_kw, slater_kw in [
        ('quadrature', {'zeff_method': 'spline'}, {}),
        ('algebraic', {}, {}),  # Uses new defaults
    ]:
        print(f"\n{'='*60}")
        print(f"LiH full PES: {label}")
        print(f"{'='*60}")

        t0 = time.time()

        # Create solver — need to patch CoreScreening call
        solver = ComposedDiatomicSolver.LiH_ab_initio(
            l_max=2,
            pk_channel_mode='l_dependent',
            verbose=True,
        )

        # Override zeff_method on core if needed
        if 'zeff_method' in zeff_kw:
            # Solve core manually with specified method
            from geovac.core_screening import CoreScreening
            solver.core = CoreScreening(
                Z=3, l_max=2, n_alpha=200,
                zeff_method=zeff_kw['zeff_method'],
            )
            solver.core.solve(verbose=True)
            solver.E_core = solver.core.energy

            # Re-derive PK
            from geovac.ab_initio_pk import AbInitioPK
            solver.ab_initio_pk = AbInitioPK(solver.core, n_core=2)
            solver.pk_A = solver.ab_initio_pk.A
            solver.pk_B = solver.ab_initio_pk.B
            pk_d = solver.ab_initio_pk.pk_dict(atom='A')
            pk_d['channel_mode'] = 'l_dependent'
            solver.pk_potentials = [pk_d]
        else:
            solver.solve_core()

        # PES scan
        R_grid = np.concatenate([
            np.linspace(2.0, 2.5, 3),
            np.linspace(2.7, 4.0, 10),
            np.linspace(4.5, 7.0, 5),
        ])
        pes = solver.scan_pes(R_grid=R_grid)
        spectro = solver.fit_spectroscopic_constants()
        total_time = time.time() - t0

        results[f'LiH_{label}'] = {
            'R_eq': spectro['R_eq'],
            'E_min': spectro['E_min'],
            'D_e': spectro['D_e'],
            'omega_e': spectro['omega_e'],
            'total_time': total_time,
            'E_core': solver.E_core,
        }

        print(f"\n  R_eq = {spectro['R_eq']:.4f} bohr")
        print(f"  D_e = {spectro['D_e']:.6f} Ha")
        print(f"  omega_e = {spectro['omega_e']:.1f} cm-1")
        print(f"  Total time: {total_time:.1f}s")

    return results


def benchmark_beh2():
    """BeH₂ single-point comparison at R_eq."""
    from geovac.composed_triatomic import ComposedTriatomicSolver

    results = {}

    for label in ['quadrature', 'algebraic']:
        print(f"\n{'='*60}")
        print(f"BeH₂: {label}")
        print(f"{'='*60}")

        t0 = time.time()

        solver = ComposedTriatomicSolver.BeH2(
            l_max=2,
            coupling_mode='full_exchange',
            verbose=True,
        )

        if label == 'quadrature':
            from geovac.core_screening import CoreScreening
            solver.core = CoreScreening(
                Z=4, l_max=2, n_alpha=200,
                zeff_method='spline',
            )
            solver.core.solve(verbose=True)
            solver.E_core = solver.core.energy
            from geovac.ab_initio_pk import AbInitioPK
            solver.ab_initio_pk = AbInitioPK(solver.core, n_core=2)
            solver.pk_A = solver.ab_initio_pk.A
            solver.pk_B = solver.ab_initio_pk.B
            pk_d = solver.ab_initio_pk.pk_dict(atom='A')
            pk_d['channel_mode'] = solver.pk_channel_mode
            solver.pk_potentials = [pk_d]
        else:
            solver.solve_core()

        # Scan
        R_grid = np.concatenate([
            np.linspace(1.8, 2.2, 3),
            np.linspace(2.3, 3.5, 8),
            np.linspace(4.0, 6.0, 4),
        ])
        pes = solver.scan_pes(R_grid=R_grid)
        spectro = solver.fit_spectroscopic_constants()
        total_time = time.time() - t0

        results[f'BeH2_{label}'] = {
            'R_eq': spectro['R_eq'],
            'E_min': spectro['E_min'],
            'D_e': spectro['D_e'],
            'total_time': total_time,
        }

        print(f"\n  R_eq = {spectro['R_eq']:.4f} bohr")
        print(f"  Total time: {total_time:.1f}s")

    return results


def print_comparison(results):
    """Print formatted comparison table."""
    print("\n" + "=" * 80)
    print("COMPARISON TABLE: Quadrature vs Algebraic Defaults")
    print("=" * 80)

    for mol in ['LiH', 'BeH2']:
        q_key = f'{mol}_quadrature'
        a_key = f'{mol}_algebraic'
        if q_key in results and a_key in results:
            q = results[q_key]
            a = results[a_key]

            print(f"\n{mol}:")
            print(f"  {'Metric':<15} {'Quadrature':>12} {'Algebraic':>12} {'Delta':>10} {'Delta%':>8}")
            print(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*10} {'-'*8}")

            for metric in ['R_eq', 'E_min', 'D_e', 'total_time']:
                qv = q.get(metric, float('nan'))
                av = a.get(metric, float('nan'))
                delta = av - qv
                pct = delta / abs(qv) * 100 if abs(qv) > 1e-15 else 0.0
                unit = 's' if metric == 'total_time' else ('Ha' if 'E' in metric or 'D' in metric else 'bohr')
                print(f"  {metric:<15} {qv:>12.6f} {av:>12.6f} {delta:>10.6f} {pct:>7.3f}%")


if __name__ == '__main__':
    print("=" * 80)
    print("Track R: Algebraic Default Benchmark")
    print("=" * 80)

    # Start with LiH full PES (most important benchmark)
    lih_results = benchmark_lih_pes()
    print_comparison(lih_results)

    # Save results
    output_dir = Path("debug/data")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "track_r_benchmark.json"

    # Convert numpy types for JSON serialization
    def to_json_safe(obj):
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    safe_results = {k: {kk: to_json_safe(vv) for kk, vv in v.items()}
                    for k, v in lih_results.items()}
    with open(output_file, 'w') as f:
        json.dump(safe_results, f, indent=2)
    print(f"\nResults saved to {output_file}")
