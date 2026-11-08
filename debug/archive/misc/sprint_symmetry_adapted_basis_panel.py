"""Verification panel for symmetry-adapted basis hidden-Z_2 detection.

Runs the combined-rotation sector decomposition + hidden-Z_2 detection on
the representative molecule panel (LiH, BeH2, H2O, NH3, CH4, HF) and writes
results to ``debug/data/sprint_symmetry_adapted_basis_panel.json``.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from geovac.molecular_spec import (
    lih_spec, beh2_spec, h2o_spec, nh3_spec, ch4_spec, hf_spec,
)
from geovac.balanced_coupled import (
    _get_nuclei_for_beh2, _get_nuclei_for_h2o,
)
from geovac.extended_tapering import extended_tapered_from_spec
from geovac.symmetry_adapted_basis import (
    decompose_hamiltonian_by_sector,
    find_hidden_z2_in_sector,
    extended_plus_hidden_tapered_from_spec,
)


def main():
    panel = [
        ('LiH', lih_spec(), None),
        ('HF', hf_spec(), None),
        ('BeH2', beh2_spec(), _get_nuclei_for_beh2(2.5)),
        ('H2O', h2o_spec(), _get_nuclei_for_h2o()),
        ('NH3', nh3_spec(), None),
        ('CH4', ch4_spec(), None),
    ]

    results = []
    print(f"{'mol':<6} {'M':<4} {'n_sect':<7} "
          f"{'max_off_h1':<12} {'max_off_eri':<12} "
          f"{'n_cand':<8} {'n_valid':<8} "
          f"{'Q_ext':<6} {'Q_eph':<6} {'save_hid':<10} {'time(s)':<8}")
    print('-' * 110)

    for name, spec, nuc in panel:
        t0 = time.time()
        sectors = decompose_hamiltonian_by_sector(
            spec, nuclei=nuc, builder='composed',
        )
        meta = sectors.pop(('__meta__',))
        cands = find_hidden_z2_in_sector(
            spec, nuclei=nuc, builder='composed',
        )
        valid = [c for c in cands if c['is_valid']]
        ext = extended_tapered_from_spec(
            spec, nuclei=nuc, use_hopf=True, use_ell_parity=True,
            use_atom_swap=bool(nuc), use_inversion=False,
        )
        eph = extended_plus_hidden_tapered_from_spec(
            spec, nuclei=nuc, use_hopf=True, use_ell_parity=True,
            use_atom_swap=bool(nuc), use_inversion=False,
        )
        sector_summary = [
            {'sector': list(k), 'dim': v['dim'],
             'orbitals': v['orbital_indices']}
            for k, v in sectors.items()
        ]
        n_pauli_ext = len(ext['qubit_op_tapered'].terms) - (
            1 if () in ext['qubit_op_tapered'].terms else 0
        )
        n_pauli_eph = len(eph['qubit_op_tapered_plus_hidden'].terms) - (
            1 if () in eph['qubit_op_tapered_plus_hidden'].terms else 0
        )
        m = sum(s['dim'] for s in sectors.values())
        wall = time.time() - t0
        row = {
            'molecule': name,
            'use_atom_swap': bool(nuc),
            'M': m,
            'Q_naive': ext['Q_naive'],
            'Q_extended': ext['Q_tapered'],
            'Q_extended_plus_hidden': eph['Q_tapered_hidden'],
            'delta_Q_extended': ext['delta_Q'],
            'delta_Q_hidden': ext['Q_tapered'] - eph['Q_tapered_hidden'],
            'n_sectors': len(sectors),
            'n_hidden_candidates': len(cands),
            'n_hidden_valid': len(valid),
            'n_hidden_kept_after_dedup_audit': len(eph['hidden_z2_kept']),
            'max_offsector_h1': meta['max_offsector_h1'],
            'max_offsector_eri': meta['max_offsector_eri'],
            'commute_hopf_swap': meta['commute_hopf_swap'],
            'pauli_extended': n_pauli_ext,
            'pauli_ext_plus_hidden': n_pauli_eph,
            'pauli_reduction_pct': 100.0 * (n_pauli_ext - n_pauli_eph) / max(n_pauli_ext, 1),
            'sectors': sector_summary,
            'wall_time_s': wall,
        }
        results.append(row)
        print(f"{name:<6} {m:<4} {len(sectors):<7} "
              f"{meta['max_offsector_h1']:<12.2e} "
              f"{meta['max_offsector_eri']:<12.2e} "
              f"{len(cands):<8} {len(valid):<8} "
              f"{ext['Q_tapered']:<6} {eph['Q_tapered_hidden']:<6} "
              f"{ext['Q_tapered'] - eph['Q_tapered_hidden']:<10} "
              f"{wall:<8.2f}")

    out = Path(__file__).parent / 'data' / 'sprint_symmetry_adapted_basis_panel.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2))
    print(f"\nWritten: {out}")

    # Summary
    print("\n=== SUMMARY ===")
    total_save = sum(r['delta_Q_hidden'] for r in results)
    total_extended = sum(r['Q_extended'] for r in results)
    total_eph = sum(r['Q_extended_plus_hidden'] for r in results)
    print(f"Total Q saved by hidden Z_2 across panel: {total_save}")
    print(f"Total Q (extended-only):    {total_extended}")
    print(f"Total Q (extended+hidden):  {total_eph}")
    print(f"Pauli reduction range: "
          f"{min(r['pauli_reduction_pct'] for r in results):.1f}%"
          f" – {max(r['pauli_reduction_pct'] for r in results):.1f}%")


if __name__ == '__main__':
    main()
