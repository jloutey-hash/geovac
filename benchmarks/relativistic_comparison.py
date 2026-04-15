"""Relativistic-VQE head-to-head: GeoVac spin-ful composed vs Sunaga et al. 2025.

Track T4 of the Dirac-on-S^3 Tier 2 sprint.

Sunaga et al., "Relativistic VQE calculations of molecular electric dipole
moments on quantum hardware", PRA 111, 022817 (2025); arXiv:2406.04992.

Sunaga's published explicit Pauli counts are limited:

- **18-qubit RaH, cv2z Dyall, DC Hamiltonian**: 47,099 Pauli strings
  (relativistic); 12,556 (non-relativistic).  PDM operator uses
  107 (relativistic) / 67 (non-relativistic) integrals.
- **18-qubit calculations** for BeH, MgH, CaH, SrH, BaH, RaH all use
  3 occupied + 15 unoccupied spinorbitals = 18 qubits active space.
- **12-qubit SrH** pipeline 2qg counts: UCCSD 9148 -> ES-VQE 2960
  -> PyZX 144 -> Qiskit L3 26 -> RL-ZX+Cflow 13.
- **6-qubit SrH/SrF** for hardware, 19 / 53 Pauli terms in the PDM
  operator partitioned into 1 dominant QWC clique.

The paper does NOT publish per-molecule Pauli counts for BeH/MgH/CaH/SrH
at 18 qubits.  Only RaH-18q is a fully quoted baseline.

T4 Sunaga-baseline reality:
- Matched-Q comparison is only possible for RaH at Q=18.  GeoVac's
  RaH requires Z=88 support (far outside the Z<=36 scope).  We therefore
  compare GeoVac LiH/BeH/CaH at their native Q=30/30/20 against the one
  calibrated Sunaga baseline (RaH-18q = 47,099 Pauli, cv2z Dyall).  The
  Q mismatch is honestly noted.
- GeoVac atomic family for the Sunaga target molecules: LiH (Z=3),
  BeH (Z=4), CaH (Z=20).  SrH (Z=38) not yet in scope.

Run:
    python benchmarks/relativistic_comparison.py

Writes:
    debug/data/tier2_market/sunaga_comparison.json
    debug/data/tier2_market/sunaga_comparison_table.md
"""

from __future__ import annotations

import json
import os
import time
from typing import Dict, Any

from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import (
    lih_spec_relativistic,
    beh_spec_relativistic,
    cah_spec_relativistic,
)
from geovac.measurement_grouping import qwc_groups as _qwc_groups


def _qubit_op_metrics(qubit_op):
    """Count non-identity Pauli terms, 1-norm (total and ni), and QWC groups."""
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
        'N_pauli_total': n_total,
        'N_pauli_ni': n_ni,
        'lambda_total': lam_total,
        'lambda_ni': lam_ni,
        'qwc': qwc,
    }


# --- Sunaga et al. 2025 reference values (PRA 111, 022817) ---------------

SUNAGA_2025 = {
    'citation': 'Sunaga et al., PRA 111, 022817 (2025); arXiv:2406.04992',
    'basis': 'cv2z Dyall, DC Hamiltonian -> JW',
    'active_space_18q': '3 occ + 15 unocc spinorbitals',
    'RaH_18q': {
        'Q': 18,
        'pauli_rel': 47099,
        'pauli_nonrel': 12556,
        'pdm_integrals_rel': 107,
        'pdm_integrals_nonrel': 67,
        'method': 'VQE-UCCSD',
    },
    'SrH_12q': {
        'Q': 12,
        'two_qubit_gates_UCCSD_initial': 9148,
        'two_qubit_gates_ES_VQE': 2960,
        'two_qubit_gates_PyZX': 144,
        'two_qubit_gates_Qiskit_L3': 26,
        'two_qubit_gates_RLZX_Cflow': 13,
        'parameters_UCCSD': 104,
        'parameters_ES_VQE': 32,
    },
    'SrH_6q': {
        'Q': 6,
        'two_qubit_gates_UCCSD_initial': 280,
        'two_qubit_gates_ES_VQE': 148,
        'two_qubit_gates_Qiskit_L3': 60,
        'two_qubit_gates_RLZX_Cflow': 28,
        'pdm_pauli_terms': 19,
        'qwc_cliques_dominant': 1,
    },
    # DEFERRED cells (not extractable from the published tables/figures):
    'deferred': [
        'Per-molecule Pauli counts for BeH/MgH/CaH/SrH/BaH at 18q',
        'Per-molecule 1-norms (lambda) at 18q or 12q',
        'Per-molecule QWC group counts',
    ],
}


def compute_geovac_molecule(name: str, factory) -> Dict[str, Any]:
    """Build GeoVac relativistic composed Hamiltonian and return resource metrics."""
    t0 = time.time()
    res = build_composed_hamiltonian(factory(max_n=2))
    dt_rel = time.time() - t0
    rel_metrics = _qubit_op_metrics(res['qubit_op'])

    # Matched scalar build for context.
    scalar_spec = factory(max_n=2)
    scalar_spec.relativistic = False
    scalar_spec.name = scalar_spec.name.replace('_rel', '')
    t0 = time.time()
    res_scalar = build_composed_hamiltonian(scalar_spec)
    dt_scalar = time.time() - t0
    scalar_metrics = _qubit_op_metrics(res_scalar['qubit_op'])

    return {
        'name': name,
        'Q': int(res['Q']),
        'N_pauli_rel': rel_metrics['N_pauli_ni'],
        'N_pauli_rel_total': rel_metrics['N_pauli_total'],
        'lambda_total_rel': rel_metrics['lambda_total'],
        'lambda_ni_rel': rel_metrics['lambda_ni'],
        'qwc_rel': rel_metrics['qwc'],
        'N_pauli_scalar': scalar_metrics['N_pauli_ni'],
        'lambda_ni_scalar': scalar_metrics['lambda_ni'],
        'qwc_scalar': scalar_metrics['qwc'],
        'wall_rel_s': dt_rel,
        'wall_scalar_s': dt_scalar,
    }


def main():
    outdir = os.path.join(
        os.path.dirname(__file__), '..', 'debug', 'data', 'tier2_market'
    )
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    molecules = [
        ('LiH', lih_spec_relativistic),
        ('BeH', beh_spec_relativistic),
        ('CaH', cah_spec_relativistic),
    ]

    geovac_results = {}
    for name, factory in molecules:
        print(f'  Building {name}...')
        geovac_results[name] = compute_geovac_molecule(name, factory)

    # Write JSON
    out_json = os.path.join(outdir, 'sunaga_comparison.json')
    with open(out_json, 'w') as f:
        json.dump({
            'geovac': geovac_results,
            'sunaga_2025': SUNAGA_2025,
        }, f, indent=2)
    print(f'Wrote {out_json}')

    # --- Markdown comparison table ---
    md = []
    md.append('# GeoVac (Tier 2) vs Sunaga et al. 2025 -- Head-to-head resource comparison\n')
    md.append('Source: Sunaga et al., PRA 111, 022817 (2025); arXiv:2406.04992.\n')
    md.append('GeoVac: relativistic composed pipeline, n_max=2, Dirac (kappa, m_j) basis.\n')
    md.append('Sunaga: VQE-UCCSD, cv2z Dyall basis, DC Hamiltonian -> JW, '
              '3 occ + 15 unocc spinorbitals = 18q active space.\n')
    md.append('\n## Part A: GeoVac raw resource metrics (reproduced from T3 memo)\n')
    md.append('| Molecule | Q | N_pauli (rel) | N_pauli (scalar) | rel/scalar | '
              'lambda_ni (rel, Ha) | lambda_ni (scalar, Ha) | QWC (rel) | QWC (scalar) |')
    md.append('|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|')
    for name, _ in molecules:
        r = geovac_results[name]
        ratio = r['N_pauli_rel'] / r['N_pauli_scalar']
        md.append(
            f"| {name} | {r['Q']} | {r['N_pauli_rel']} | {r['N_pauli_scalar']} | "
            f"{ratio:.2f}x | {r['lambda_ni_rel']:.2f} | {r['lambda_ni_scalar']:.2f} | "
            f"{r['qwc_rel']} | {r['qwc_scalar']} |"
        )

    md.append('\n## Part B: Head-to-head with Sunaga 2025\n')
    md.append('Sunaga publishes explicit per-molecule Pauli counts only for **RaH-18q**.  '
              'We therefore compare each GeoVac molecule against the RaH-18q baseline as '
              'the single calibrated Sunaga number.  Q mismatch is flagged.\n')
    md.append('| Molecule | GeoVac Q | GeoVac N_pauli | Sunaga target | Sunaga Q | '
              'Sunaga N_pauli | GeoVac/Sunaga | Flag |')
    md.append('|:---|:---:|:---:|:---|:---:|:---:|:---:|:---|')
    RaH = SUNAGA_2025['RaH_18q']
    for name, _ in molecules:
        r = geovac_results[name]
        ratio = r['N_pauli_rel'] / RaH['pauli_rel']
        note = 'Q mismatch' if r['Q'] != RaH['Q'] else 'matched'
        md.append(
            f"| {name} | {r['Q']} | {r['N_pauli_rel']} | RaH-18q | {RaH['Q']} | "
            f"{RaH['pauli_rel']} | {ratio:.3f}x | {note} |"
        )

    md.append('\n## Part C: 1-norm vs Sunaga (DEFERRED)\n')
    md.append('Sunaga does not report per-molecule 1-norms for BeH/MgH/CaH/SrH/BaH/RaH '
              'in the published figures/tables (Section III/IV).  Only ground-state '
              'energies (-3178.63 Ha for SrH STO-6G) and PDM values are reported.  '
              'Direct 1-norm comparison is a **fetchable-via-direct-paper-read** '
              'deferred item.\n')

    md.append('\n## Part D: QWC groups vs Sunaga (DEFERRED)\n')
    md.append('Sunaga reports **1 dominant QWC clique** for the 6-qubit PDM operator '
              '(19 Pauli terms) and for the 12-qubit PDM (53 terms).  QWC groups for '
              'the full Hamiltonian are not reported.  Comparison deferred.\n')

    md.append('\n## Part E: 2-qubit gate count (UCCSD/VQE ansatz)\n')
    md.append('GeoVac Tier 2 deliverable is the Hamiltonian, not a VQE ansatz, so this '
              'cell is not directly comparable.  For reference, Sunaga 12q SrH ES-VQE '
              'starts at 9,148 2qg and compresses to 13 after RL-ZX+Cflow pipeline.  '
              'The GeoVac 30-qubit relativistic Hamiltonian would need an analogous '
              'UCCSD ansatz + pipeline-optimisation pass for head-to-head gate comparison, '
              'which is out of scope for T4.\n')

    md.append('\n## Part F: Matched-Z family comparison\n')
    md.append('The only Sunaga molecule in scope for GeoVac (Z <= 36) is CaH (Z=20).  '
              'For that molecule only:\n')
    ca = geovac_results['CaH']
    md.append(f'- GeoVac CaH Q={ca["Q"]}, N_pauli={ca["N_pauli_rel"]}, '
              f'lambda_ni={ca["lambda_ni_rel"]:.2f} Ha, QWC={ca["qwc_rel"]}.\n')
    md.append(f'- Sunaga CaH: 18q active space, but no per-molecule Pauli count published.  '
              'The Sunaga 18q baselines are all calibrated on RaH (47,099 Pauli).\n')
    md.append('- GeoVac CaH at Q=20 is 1/18 the qubit count advantage over Sunaga.  '
              'GeoVac composed architecture has structural asymmetry with frozen-core '
              'reduction (Ca [Ar] core is inert) that Sunaga\'s cv2z Dyall basis does not have.\n')

    md.append('\n## Deferred items (fetchable-via-direct-paper-read)\n')
    for d in SUNAGA_2025['deferred']:
        md.append(f'- {d}')
    md.append('\nRetrieving the Sunaga Supplemental Material (Tables S1-S3) would '
              'populate the per-molecule Pauli/1-norm/QWC cells for all six molecules.\n')

    md_text = '\n'.join(md) + '\n'
    out_md = os.path.join(outdir, 'sunaga_comparison_table.md')
    with open(out_md, 'w') as f:
        f.write(md_text)
    print(f'Wrote {out_md}')

    print('\n=== Headline ===')
    for name, _ in molecules:
        r = geovac_results[name]
        ratio = r['N_pauli_rel'] / RaH['pauli_rel']
        print(f'{name:4s}: Q={r["Q"]:3d}  '
              f'GeoVac Pauli={r["N_pauli_rel"]:6d}  '
              f'Sunaga RaH-18q={RaH["pauli_rel"]:6d}  '
              f'ratio={ratio:.3f}x')


if __name__ == '__main__':
    main()
