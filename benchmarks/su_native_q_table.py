"""SU-B / SU-C native-Q resource table for single-bond alkali / alkaline-earth
   hydrides across the GeoVac library + Sunaga 2025 comparison.

Part of Sprint 3 Track SU (Sunaga head-to-head market test).

Builds scalar and relativistic composed qubit Hamiltonians at n_max=2 for:
  - LiH (Z=3):      1s^2 core + 2s + bond block ("structure B", Li-like core)
  - NaH (Z=11):     [Ne] frozen core + bond block (alkali monohydride)
  - KH  (Z=19):     [Ar] frozen core + bond block
  - CaH (Z=20):     [Ar] frozen core + bond block (alkaline-earth monohydride)
  - SrH (Z=38):     [Kr] frozen core + bond block (Sunaga target)
  - BaH (Z=56):     [Xe] frozen core + bond block (Sunaga target)

LiH is structure-different from the alkali/alkaline-earth monohydrides:
  Li has 1s^2 explicit core + 2s, making LiH a three-block composed spec
  (core + lone + bond), Q=30 at n_max=2.  The alkali/alkaline-earth
  monohydrides (NaH, KH, CaH, SrH, BaH) share identical block topology
  (single frozen-core + single bond block), Q=20 at n_max=2.

Comparison baseline (CORRECTED from arXiv:2406.04992 HTML, 2026-04-15):
  Sunaga et al. 2025 (PRA 111, 022817), RaH at Q=18:
    - Relativistic Hamiltonian: 12,556 Pauli strings (47,099 integrals)
    - Non-relativistic Hamiltonian:  2,740 Pauli strings (4,249 integrals)
  The previous T4 memo conflated Pauli strings with integral counts.

Writes:
    debug/data/su_native_q_table.json
    debug/data/su_resource_tables.json (superset for SU-D)
"""
from __future__ import annotations

import json
import os
import time
from typing import Dict, Any, Tuple

from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import (
    lih_spec, nah_spec, kh_spec,
    lih_spec_relativistic,
    srh_spec_relativistic, bah_spec_relativistic, cah_spec_relativistic,
    _alkaline_earth_monohydride_spec,
)
from geovac.measurement_grouping import qwc_groups as _qwc_groups


# Corrected Sunaga 2025 reference values (from arXiv:2406.04992 HTML, Section Theory/Methodology)
SUNAGA_RAH_18Q = {
    'citation': 'Sunaga et al., PRA 111, 022817 (2025); arXiv:2406.04992',
    'molecule': 'RaH',
    'Q': 18,
    'active_space': '3 occ + 15 unocc spinorbitals',
    'basis': 'cv2z Dyall, DC Hamiltonian, JW mapping',
    # CORRECTED values:
    'pauli_rel': 12556,            # Pauli strings in rel. Hamiltonian
    'integrals_rel': 47099,        # One- and two-electron integrals (rel)
    'pauli_nonrel': 2740,          # Pauli strings in NR Hamiltonian
    'integrals_nonrel': 4249,      # One- and two-electron integrals (NR)
    # PDM operator (NOT Hamiltonian):
    'pauli_pdm_rel': 107,
    'integrals_pdm_rel': 162,
    'pauli_pdm_nonrel': 67,
    'integrals_pdm_nonrel': 66,
    # Previous T4 memo mistake (record for provenance):
    '_prior_t4_error': (
        'The Tier 2 T4 memo and benchmarks/relativistic_comparison.py '
        'had 47099 labeled as Pauli (rel) and 12556 as Pauli (nonrel). '
        'These were integrals (rel) and Pauli (rel) respectively. '
        'Verified 2026-04-15 from arXiv HTML; see '
        'debug/su_sunaga_comparison.md.'
    ),
}


def _ala_spec_scalar(Z: int, name: str, max_n: int = 2):
    """Scalar alkaline-earth monohydride spec (wrapper for _alkaline_earth_monohydride_spec)."""
    return _alkaline_earth_monohydride_spec(
        Z, name, max_n=max_n, relativistic=False,
    )


def _qubit_op_metrics(qubit_op) -> Dict[str, Any]:
    """Non-identity Pauli count, 1-norm (total/ni), and QWC groups."""
    n_total = 0
    n_ni = 0
    lam_total = 0.0
    lam_ni = 0.0
    for term, coef in qubit_op.terms.items():
        c = abs(complex(coef))
        n_total += 1
        lam_total += c
        if term:  # non-identity
            n_ni += 1
            lam_ni += c
    try:
        qwc = len(_qwc_groups(qubit_op))
    except Exception:
        qwc = -1
    return {
        'N_pauli_total': int(n_total),
        'N_pauli_ni': int(n_ni),
        'lambda_total': float(lam_total),
        'lambda_ni': float(lam_ni),
        'qwc': int(qwc),
    }


def _build_and_measure(spec) -> Dict[str, Any]:
    t0 = time.time()
    res = build_composed_hamiltonian(spec)
    dt = time.time() - t0
    m = _qubit_op_metrics(res['qubit_op'])
    m['Q'] = int(res['Q'])
    m['M'] = int(res['M'])
    m['wall_s'] = float(dt)
    m['name'] = spec.name
    return m


def scalar_spec_for(name: str, Z: int, max_n: int = 2):
    """Return the scalar spec factory result for given molecule."""
    if name == 'LiH':
        return lih_spec(max_n=max_n)
    elif name == 'NaH':
        return nah_spec(max_n=max_n)
    elif name == 'KH':
        return kh_spec(max_n=max_n)
    elif name in ('CaH', 'SrH', 'BaH'):
        return _ala_spec_scalar(Z, name, max_n=max_n)
    raise ValueError(name)


def rel_spec_for(name: str, Z: int, max_n: int = 2):
    """Return the relativistic spec factory result, or None if not available."""
    if name == 'LiH':
        return lih_spec_relativistic(max_n=max_n)
    elif name == 'CaH':
        return cah_spec_relativistic(max_n=max_n)
    elif name == 'SrH':
        return srh_spec_relativistic(max_n=max_n)
    elif name == 'BaH':
        return bah_spec_relativistic(max_n=max_n)
    # NaH, KH have no production _relativistic variant yet
    return None


def main():
    outdir = os.path.join(os.path.dirname(__file__), '..', 'debug', 'data')
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    # Single-bond alkali / alkaline-earth library at n_max=2
    targets: list[Tuple[str, int]] = [
        ('LiH',  3),   # structure-different (has explicit Li 1s^2 core + 2s)
        ('NaH', 11),
        ('KH',  19),
        ('CaH', 20),
        ('SrH', 38),
        ('BaH', 56),
    ]

    rows = []
    for name, Z in targets:
        print(f'Building {name} (Z={Z})...')
        scalar = _build_and_measure(scalar_spec_for(name, Z))
        rel_spec_obj = rel_spec_for(name, Z)
        rel = _build_and_measure(rel_spec_obj) if rel_spec_obj is not None else None
        row = {
            'molecule': name,
            'Z': Z,
            'n_max': 2,
            'scalar': scalar,
            'relativistic': rel,
        }
        if rel is not None:
            row['rel_over_scalar_pauli'] = rel['N_pauli_ni'] / max(scalar['N_pauli_ni'], 1)
        # Ratios to Sunaga RaH-18q relativistic (CORRECTED):
        sunaga_pauli = SUNAGA_RAH_18Q['pauli_rel']
        row['native_q_ratio_scalar_vs_sunaga'] = scalar['N_pauli_ni'] / sunaga_pauli
        if rel is not None:
            row['native_q_ratio_rel_vs_sunaga'] = rel['N_pauli_ni'] / sunaga_pauli
        rows.append(row)

    out_json = os.path.join(outdir, 'su_native_q_table.json')
    with open(out_json, 'w') as f:
        json.dump({
            'sprint': 'Sprint 3 Track SU Part 2 (SU-B)',
            'date': '2026-04-15',
            'description': 'Native-Q resource table for single-bond '
                           'alkali / alkaline-earth hydrides in the GeoVac '
                           'library at n_max=2, scalar and relativistic '
                           'composed pipelines.',
            'sunaga_reference': SUNAGA_RAH_18Q,
            'rows': rows,
        }, f, indent=2)
    print(f'\nWrote {out_json}')

    # Emit a concise markdown table
    md_lines = []
    md_lines.append('# SU-B native-Q resource table')
    md_lines.append('')
    md_lines.append(
        '| Molecule | Z | Q | N_Pauli (sca) | N_Pauli (rel) | '
        'lambda_ni (sca, Ha) | lambda_ni (rel, Ha) | QWC (rel) | '
        'rel/sca | vs Sunaga (rel) |'
    )
    md_lines.append(
        '|----------|---|---|---|---|---|---|---|---|---|'
    )
    for row in rows:
        sca = row['scalar']
        rel = row['relativistic']
        if rel is not None:
            md_lines.append(
                f"| {row['molecule']} | {row['Z']} | "
                f"{sca['Q']} | {sca['N_pauli_ni']} | {rel['N_pauli_ni']} | "
                f"{sca['lambda_ni']:.2f} | {rel['lambda_ni']:.2f} | "
                f"{rel['qwc']} | "
                f"{row['rel_over_scalar_pauli']:.2f}x | "
                f"{row['native_q_ratio_rel_vs_sunaga']:.4f}x |"
            )
        else:
            md_lines.append(
                f"| {row['molecule']} | {row['Z']} | "
                f"{sca['Q']} | {sca['N_pauli_ni']} | --- | "
                f"{sca['lambda_ni']:.2f} | --- | --- | --- | "
                f"sca {row['native_q_ratio_scalar_vs_sunaga']:.4f}x |"
            )
    md_lines.append('')
    md_lines.append(
        f"Sunaga RaH Q=18 rel = {SUNAGA_RAH_18Q['pauli_rel']} Pauli strings "
        f"(CORRECTED: previous T4 memo had 47,099 which is the integral count, "
        f"not Pauli string count)."
    )

    md_text = '\n'.join(md_lines) + '\n'
    print('\n' + md_text)

    # Also save superset resource table for SU-D (Paper 20 update)
    out_super = os.path.join(outdir, 'su_resource_tables.json')
    with open(out_super, 'w') as f:
        json.dump({
            'sprint': 'Sprint 3 Track SU',
            'date': '2026-04-15',
            'native_q_table': rows,
            'sunaga_reference_corrected': SUNAGA_RAH_18Q,
            'notes': {
                'sunaga_correction': (
                    'The reference values for Sunaga 2025 (RaH Q=18) '
                    'have been CORRECTED. Previous value of 47,099 for '
                    'relativistic Pauli count was actually the integral '
                    'count. Correct values: 12,556 Pauli (rel), 47,099 '
                    'integrals (rel), 2,740 Pauli (nonrel), 4,249 integrals '
                    '(nonrel). Verified 2026-04-15 from arXiv:2406.04992 '
                    'HTML, Theory/Methodology section.'
                ),
                'supplemental_status': (
                    'Sunaga Tables S1-S3 contain ground-state energies '
                    '(S1), hardware PDM measurements (S2), and clique-wise '
                    'PDM contributions (S3).  They do NOT contain per-molecule '
                    'Pauli term counts, 1-norms, or QWC groups for the '
                    '18-qubit Hamiltonians of BeH/MgH/CaH/SrH/BaH.  The '
                    'previously-flagged DEFERRED matched-Q comparison is '
                    'therefore not fetchable from the supplement — the '
                    'public paper+SI does not publish per-molecule '
                    'Hamiltonian resource counts.'
                ),
                'matched_q_feasibility': (
                    'GeoVac natively operates at Q = 2*M where M is the total '
                    'spatial-orbital count across composed blocks.  At '
                    'n_max=2 for a single-bond alkali/alkaline-earth '
                    'monohydride with one frozen-core block and one bond '
                    'block, M = 10 giving Q = 20.  To reach Q=18 would '
                    'require either dropping one spatial orbital from the '
                    'closed 5-orbital (1s, 2s, 2p_x, 2p_y, 2p_z) shell of '
                    'a block (breaking the natural Fock s+p angular shell) '
                    'or reducing a block to n_max=1 (giving Q<20 but in '
                    'a non-natural way).  Matched-Q=18 is a structural '
                    'mismatch, not a GeoVac limitation; the native-Q '
                    'point is the intended comparison.'
                ),
            },
        }, f, indent=2)
    print(f'Wrote {out_super}')


if __name__ == '__main__':
    main()
