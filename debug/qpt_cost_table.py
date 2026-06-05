"""Burkat-Fitzpatrick QPT Toffoli-cost table for the GeoVac library.

Cost model: O(d^3) Toffoli for the unitary-group Paldus transform on d
spatial orbitals (Burkat-Fitzpatrick arXiv:2506.09151).

For comparison context, also report:
  - GeoVac composed Pauli count (already in the library benchmarks)
  - Standard Bravyi 2017 + per-block Hopf-Z2 reductions
  - Joint sector-dim reductions (LiH-style, computed for each molecule)

The QPT circuit operates on the qubit register AFTER the Hopf-Z2 tapering;
the d in the cost formula is the number of spatial orbitals retained,
which is M (the composed builder's spatial orbital count) before tapering
or M - Delta_Q_spatial after.

For the SCOPING memo: report all 28 GeoVac library molecules.
"""

from __future__ import annotations

import json
import os
from typing import Dict, List

import numpy as np

from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import (
    lih_spec, beh2_spec, h2o_spec, ch4_spec, nh3_spec, hf_spec,
    nah_spec, mgh2_spec, sih4_spec, ph3_spec, h2s_spec, hcl_spec,
    kh_spec, cah2_spec, geh4_spec, ash3_spec, h2se_spec, hbr_spec,
)


# Spec functions to test (subset of library; the rest are multi-center and
# need spec construction through the lif_spec / co_spec / etc. helpers).
_LIBRARY = [
    ('LiH', lih_spec),
    ('BeH2', beh2_spec),
    ('CH4', ch4_spec),
    ('NH3', nh3_spec),
    ('H2O', h2o_spec),
    ('HF', hf_spec),
    ('NaH', nah_spec),
    ('MgH2', mgh2_spec),
    ('SiH4', sih4_spec),
    ('PH3', ph3_spec),
    ('H2S', h2s_spec),
    ('HCl', hcl_spec),
    ('KH', kh_spec),
    ('CaH2', cah2_spec),
    ('GeH4', geh4_spec),
    ('AsH3', ash3_spec),
    ('H2Se', h2se_spec),
    ('HBr', hbr_spec),
]


def main():
    results = []

    for label, spec_fn in _LIBRARY:
        try:
            spec = spec_fn()
            res = build_composed_hamiltonian(spec, verbose=False)
            M = int(res['h1'].shape[0])
            Q_full = 2 * M
            N_pauli = int(res['N_pauli'])

            # Pre-tapering qubit count
            Q = Q_full

            # Standard Bravyi 2017 tapering: -2 qubits (N, Sz parity)
            Q_after_bravyi = Q - 2

            # Hopf-Z2 per-sub-block: ΔQ = 2 + n_sub_blocks (one P per sub-block
            # with at least one l>=1 orbital; estimate from spec).
            # For first-row hydrides, n_sb is typically the number of blocks.
            # For frozen-core systems, n_sb depends on the active space.
            n_blocks = len(spec.blocks)
            # Per-sub-block estimate: count blocks with l_min < l_max (i.e.,
            # blocks containing at least one l>=1 orbital).
            n_sb_estimate = sum(
                1 for blk in spec.blocks if blk.max_n >= 2
            )
            # Hopf-Z2 ΔQ
            delta_Q_hopf = 2 + n_sb_estimate
            Q_after_hopf = Q - delta_Q_hopf

            # QPT Toffoli cost: O(M^3) operating on spatial orbitals
            # Lit-search agent estimate: ~M^3 with implicit constant ~1
            qpt_toffoli_naive = M ** 3

            # For reference: typical QPE cost for a chemistry Hamiltonian
            # ~lambda * t_QPE_per_unit (lambda = 1-norm of H).
            lam_total = float(res.get('lambda_total', 0.0))
            lam_ni = float(res.get('lambda_ni', 0.0))

            results.append({
                'molecule': label,
                'M': M,
                'Q_full': Q_full,
                'Q_after_bravyi': Q_after_bravyi,
                'Q_after_hopf_per_block': Q_after_hopf,
                'delta_Q_hopf': delta_Q_hopf,
                'N_pauli_composed': N_pauli,
                'lambda_total': lam_total,
                'lambda_ni': lam_ni,
                'qpt_toffoli_M_cubed': qpt_toffoli_naive,
                'n_blocks': n_blocks,
                'n_sub_blocks_hopf': n_sb_estimate,
            })

            print(f"  {label:6s}: M={M:3d}, Q={Q:3d}, "
                  f"N_pauli={N_pauli:5d}, "
                  f"dQ_hopf={delta_Q_hopf:2d}, "
                  f"QPT~M^3={qpt_toffoli_naive:8d}")
        except Exception as e:
            print(f"  {label:6s}: ERROR {e}")
            results.append({'molecule': label, 'error': str(e)})

    # Summary statistics
    successful = [r for r in results if 'error' not in r]
    if successful:
        print(f"\n=== Summary across {len(successful)} molecules ===")
        Q_vals = [r['Q_full'] for r in successful]
        qpt_costs = [r['qpt_toffoli_M_cubed'] for r in successful]
        n_pauli = [r['N_pauli_composed'] for r in successful]

        print(f"  Q range:                 {min(Q_vals)} - {max(Q_vals)}")
        print(f"  QPT M^3 Toffoli range:   {min(qpt_costs)} - {max(qpt_costs)}")
        print(f"  Composed Pauli range:    {min(n_pauli)} - {max(n_pauli)}")
        print(f"  dQ_hopf range:           {min(r['delta_Q_hopf'] for r in successful)} - "
              f"{max(r['delta_Q_hopf'] for r in successful)}")

        # Total QPT cost across library
        total_qpt = sum(qpt_costs)
        print(f"  Total QPT cost (sum M^3): {total_qpt}")

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'qpt_cost_table.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump({
            'library_panel': results,
            'qpt_scaling_model': 'O(M^3) Toffoli, Burkat-Fitzpatrick arXiv:2506.09151',
            'notes': (
                "M^3 is the leading-order Toffoli cost from the lit-search "
                "summary. The Burkat-Fitzpatrick paper has explicit constants "
                "not captured here; this is a scoping estimate. The qubit count "
                "Q is reported pre-tapering; QPT operates on the M spatial "
                "orbital register (== Q/2 in standard JW). After Hopf-Z2 "
                "tapering, the effective M may be smaller."
            ),
        }, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
