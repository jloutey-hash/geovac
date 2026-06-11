"""
Track DI Sprint 4: He excited states and atomic spectrum.

Computes He singlet + triplet S, P, D states from graph-native CI
and the 2D variational solver. Compares to Drake's NR reference values.
"""

import json
import numpy as np
from pathlib import Path

from geovac.casimir_ci import compute_he_spectrum, HE_NR_REFERENCE
from geovac.level3_variational import solve_he_variational_2d


def run_graph_native_spectrum() -> dict:
    """Compute He spectrum from graph-native CI at multiple n_max."""
    results = {}

    for n_max in range(3, 6):
        print(f"\n{'='*60}")
        print(f"Graph-native CI: n_max = {n_max}")
        print(f"{'='*60}")

        result = compute_he_spectrum(n_max=n_max, n_states=3,
                                     max_L=min(2, n_max - 1))

        print(f"\nSector dimensions: {result['sector_dims']}")

        print(f"\n{'State':<10} {'Energy (Ha)':>14} {'Reference':>14} {'Error %':>10}")
        print("-" * 50)
        for label in sorted(result['states'].keys()):
            s = result['states'][label]
            ref_str = f"{s['reference']:.6f}" if s['reference'] else "N/A"
            err_str = f"{s['error_pct']:.3f}" if s['error_pct'] is not None else "N/A"
            print(f"{label:<10} {s['energy']:>14.6f} {ref_str:>14} {err_str:>10}")

        if result['transitions']:
            print(f"\n{'Transition':<20} {'dE comp':>10} {'dE exact':>10} {'Err %':>8}")
            print("-" * 50)
            for label, t in sorted(result['transitions'].items()):
                print(f"{label:<20} {t['dE_computed']:>10.6f} "
                      f"{t['dE_exact']:>10.6f} {t['error_pct']:>8.3f}")

        # Serialize for JSON
        states_ser = {}
        for label, s in result['states'].items():
            states_ser[label] = {
                'energy': float(s['energy']),
                'reference': float(s['reference']) if s['reference'] is not None else None,
                'error_pct': float(s['error_pct']) if s['error_pct'] is not None else None,
            }
        transitions_ser = {}
        for label, t in result['transitions'].items():
            transitions_ser[label] = {
                'dE_computed': float(t['dE_computed']),
                'dE_exact': float(t['dE_exact']),
                'error_pct': float(t['error_pct']),
            }

        results[f'n_max_{n_max}'] = {
            'states': states_ser,
            'transitions': transitions_ser,
            'sector_dims': result['sector_dims'],
        }

    return results


def run_variational_2d_spectrum() -> dict:
    """Compute 2^3S from 2D variational solver."""
    results = {}

    print(f"\n{'='*60}")
    print("2D Variational Solver: Singlet vs Triplet")
    print(f"{'='*60}")

    for symmetry in ['singlet', 'triplet']:
        for l_max in [0, 1, 2]:
            print(f"\n  symmetry={symmetry}, l_max={l_max}:")
            res = solve_he_variational_2d(
                Z=2, n_basis_R=20, n_basis_alpha=20, l_max=l_max,
                symmetry=symmetry, n_states=3,
            )
            energies = res['energies'][:3]
            key = f'{symmetry}_lmax{l_max}'
            results[key] = {
                'energies': [float(e) for e in energies],
                'dim_total': res['dim_total'],
            }

            ref_key = '1_1S' if symmetry == 'singlet' else '2_3S'
            ref = HE_NR_REFERENCE.get(ref_key, None)

            for i, E in enumerate(energies):
                ref_str = ""
                if i == 0 and ref is not None:
                    err = abs((E - ref) / ref) * 100.0
                    ref_str = f"  (ref={ref:.6f}, err={err:.3f}%)"
                print(f"    E[{i}] = {E:.6f}{ref_str}")

    return results


def main() -> None:
    """Run full spectrum analysis and save results."""
    all_results = {}

    all_results['graph_native'] = run_graph_native_spectrum()
    all_results['variational_2d'] = run_variational_2d_spectrum()

    # Save
    out_path = Path(__file__).parent / 'data' / 'track_di_he_spectrum.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
