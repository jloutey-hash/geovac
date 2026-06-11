"""
Convergence study: Level 3 coupled-channel He solver at l_max=0..5 with q_mode='exact'.

Tests both n_channels=3 and n_channels=5 for l_max>=3.
Saves results to debug/data/lmax_convergence_exact_v2.json.
"""

import json
import os
import sys
import time
import signal
import traceback

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

EXACT_HE = -2.903724377  # Pekeris value

TIMEOUT = 600  # seconds per run


def run_single(l_max, n_channels, timeout=TIMEOUT):
    """Run solver for given l_max and n_channels, with wall-clock timing."""
    label = f"l_max={l_max}, n_ch={n_channels}"
    print(f"\n{'='*60}")
    print(f"Running: {label}")
    print(f"{'='*60}")

    t0 = time.time()
    try:
        result = solve_hyperspherical_algebraic_coupled(
            Z=2.0,
            n_basis=15,
            n_channels=n_channels,
            l_max=l_max,
            n_R=200,
            N_R_radial=3000,
            q_mode='exact',
            verbose=True,
        )
        dt = time.time() - t0

        energy = float(result['energy'])
        error_pct = abs((energy - EXACT_HE) / EXACT_HE) * 100
        total_dim = int(result.get('total_dim', result.get('n_basis', 0) * (l_max + 1)))

        entry = {
            'l_max': l_max,
            'n_channels': n_channels,
            'energy': energy,
            'error_pct': error_pct,
            'wall_time_s': round(dt, 2),
            'total_dim': total_dim,
            'status': 'ok',
        }
        print(f"  Energy: {energy:.6f} Ha, Error: {error_pct:.4f}%, Time: {dt:.1f}s, Dim: {total_dim}")
        return entry

    except Exception as e:
        dt = time.time() - t0
        print(f"  FAILED after {dt:.1f}s: {e}")
        traceback.print_exc()
        return {
            'l_max': l_max,
            'n_channels': n_channels,
            'energy': None,
            'error_pct': None,
            'wall_time_s': round(dt, 2),
            'total_dim': None,
            'status': f'error: {str(e)[:100]}',
        }


def main():
    results = []

    # Part 1: l_max = 0..5 with n_channels=3
    for lm in range(6):
        entry = run_single(lm, n_channels=3)
        results.append(entry)

    # Part 2: l_max = 3,4,5 with n_channels=5
    for lm in [3, 4, 5]:
        entry = run_single(lm, n_channels=5)
        results.append(entry)

    # Save JSON
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'lmax_convergence_exact_v2.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")

    # Print summary table
    print(f"\n{'='*80}")
    print(f"{'l_max':>5} {'n_ch':>5} {'Energy (Ha)':>14} {'Error%':>10} {'Time(s)':>10} {'Dim':>6} {'Status':>10}")
    print(f"{'-'*80}")
    for r in results:
        e_str = f"{r['energy']:.6f}" if r['energy'] is not None else "N/A"
        err_str = f"{r['error_pct']:.4f}" if r['error_pct'] is not None else "N/A"
        dim_str = str(r['total_dim']) if r['total_dim'] is not None else "N/A"
        print(f"{r['l_max']:>5} {r['n_channels']:>5} {e_str:>14} {err_str:>10} {r['wall_time_s']:>10.1f} {dim_str:>6} {r['status']:>10}")
    print(f"{'='*80}")


if __name__ == '__main__':
    main()
