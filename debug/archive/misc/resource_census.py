"""
GeoVac Full Library Resource Census
====================================


Builds all supported molecules via the ecosystem export API and
collects quantum resource metrics:
  - Q (qubits)
  - N_pauli (non-identity Pauli terms)
  - 1-norm (Ha)
  - Pauli/Q ratio

Outputs a formatted markdown table and JSON data file.

Usage:
    python debug/resource_census.py

Author: GeoVac Development Team
Date: April 2026
"""

import json
import os
import sys
import time
from pathlib import Path

# Ensure the local project root is first on sys.path so that
# the local geovac package is imported, not an editable install
# pointing elsewhere (e.g. OneDrive).
_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from geovac.ecosystem_export import hamiltonian


# All supported systems in display order
SYSTEMS = [
    # Atomic / diatomic
    ('He',   'atomic'),
    ('H2',   'diatomic'),
    # First-row main-group hydrides
    ('LiH',  '1st-row'),
    ('BeH2', '1st-row'),
    ('CH4',  '1st-row'),
    ('NH3',  '1st-row'),
    ('H2O',  '1st-row'),
    ('HF',   '1st-row'),
    # Second-row ([Ne] frozen core)
    ('NaH',  '2nd-row'),
    ('MgH2', '2nd-row'),
    ('SiH4', '2nd-row'),
    ('PH3',  '2nd-row'),
    ('H2S',  '2nd-row'),
    ('HCl',  '2nd-row'),
    # Third-row s-block ([Ar] frozen core)
    ('KH',   '3rd-row-s'),
    ('CaH2', '3rd-row-s'),
    # Third-row p-block ([Ar]3d10 frozen core)
    ('GeH4', '3rd-row-p'),
    ('AsH3', '3rd-row-p'),
    ('H2Se', '3rd-row-p'),
    ('HBr',  '3rd-row-p'),
    # Transition metal hydrides (Z=21-30)
    ('ScH',  'TM'),
    ('TiH',  'TM'),
    ('VH',   'TM'),
    ('CrH',  'TM'),
    ('MnH',  'TM'),
    ('FeH',  'TM'),
    ('CoH',  'TM'),
    ('NiH',  'TM'),
    ('CuH',  'TM'),
    ('ZnH',  'TM'),
]


def run_census() -> list:
    """Build all molecules and collect resource metrics."""
    results = []
    total_t0 = time.perf_counter()

    for system, category in SYSTEMS:
        t0 = time.perf_counter()
        try:
            H = hamiltonian(system, verbose=False)
            wall = time.perf_counter() - t0

            Q = H.n_qubits
            n_terms = H.n_terms
            # Non-identity count
            n_pauli = n_terms - 1  # subtract identity term
            one_norm = H.one_norm
            ratio = n_pauli / Q if Q > 0 else 0.0

            results.append({
                'system': system,
                'category': category,
                'Q': Q,
                'N_pauli': n_pauli,
                'N_terms': n_terms,
                'one_norm': round(one_norm, 2),
                'pauli_per_Q': round(ratio, 2),
                'wall_s': round(wall, 1),
                'status': 'OK',
            })
            print(f'  {system:6s}  Q={Q:3d}  N_pauli={n_pauli:5d}  '
                  f'1-norm={one_norm:10.2f}  ratio={ratio:.2f}  ({wall:.1f}s)')
        except Exception as e:
            wall = time.perf_counter() - t0
            results.append({
                'system': system,
                'category': category,
                'status': 'FAIL',
                'error': str(e),
                'wall_s': round(wall, 1),
            })
            print(f'  {system:6s}  FAIL: {e}')

    total_wall = time.perf_counter() - total_t0
    print(f'\nTotal wall time: {total_wall:.1f}s')
    print(f'Systems built: {sum(1 for r in results if r["status"] == "OK")}/{len(results)}')
    return results


def format_markdown_table(results: list) -> str:
    """Format results as a markdown table."""
    lines = [
        '| System | Category | Q | N_Pauli | 1-Norm (Ha) | Pauli/Q |',
        '|:-------|:---------|--:|--------:|------------:|--------:|',
    ]
    for r in results:
        if r['status'] != 'OK':
            lines.append(f'| {r["system"]} | {r["category"]} | — | FAIL | — | — |')
            continue
        lines.append(
            f'| {r["system"]:6s} | {r["category"]:9s} | '
            f'{r["Q"]:3d} | {r["N_pauli"]:5d} | '
            f'{r["one_norm"]:11.2f} | {r["pauli_per_Q"]:6.2f} |'
        )
    return '\n'.join(lines)


def main() -> None:
    print('GeoVac Full Library Resource Census')
    print('=' * 60)

    results = run_census()

    # Save JSON
    out_dir = Path('debug/data')
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / 'resource_census.json'
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nJSON saved to {json_path}')

    # Print markdown table
    table = format_markdown_table(results)
    print('\n' + table)

    # Save markdown
    md_path = out_dir / 'resource_census.md'
    with open(md_path, 'w') as f:
        f.write('# GeoVac Resource Census\n\n')
        f.write(f'**{sum(1 for r in results if r["status"] == "OK")}** molecules at max_n=2\n\n')
        f.write(table + '\n')
    print(f'Markdown saved to {md_path}')

    # Summary statistics
    ok = [r for r in results if r['status'] == 'OK']
    if ok:
        main_group = [r for r in ok if r['category'] not in ('atomic', 'diatomic', 'TM')]
        tm = [r for r in ok if r['category'] == 'TM']
        if main_group:
            ratios = [r['pauli_per_Q'] for r in main_group]
            print(f'\nMain-group Pauli/Q: {min(ratios):.2f} - {max(ratios):.2f} '
                  f'(mean {sum(ratios)/len(ratios):.2f})')
        if tm:
            ratios = [r['pauli_per_Q'] for r in tm]
            print(f'TM hydride Pauli/Q: {min(ratios):.2f} - {max(ratios):.2f} '
                  f'(mean {sum(ratios)/len(ratios):.2f})')
        q_vals = [r['Q'] for r in ok]
        print(f'Qubit range: {min(q_vals)} - {max(q_vals)}')


if __name__ == '__main__':
    main()
