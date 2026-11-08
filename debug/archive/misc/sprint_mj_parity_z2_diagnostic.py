"""Sprint m_j-parity Z₂ — direct commutation audit on relativistic chemistry.

Follow-on from v3.92.0 κ-parity NEGATIVE-STRUCTURAL closure. The Dirac-Coulomb
operator conserves $j$ and $m_j$ but NOT $\\kappa$. This sprint tests the
analogous parity stabilizer P_{m_j} = prod_{q: m_j(q) < 0} Z_q in the
ORIGINAL basis (each Dirac orbital is one JW qubit; no rotation).

Verdict gate:
  POSITIVE: residual <1e-10 (per-block or global) → ship as production tapering
  NEGATIVE: residual >1e-6 → document mechanism; record structural reason

The relativistic chemistry ERI has selection rule m_j_a + m_j_b = m_j_c + m_j_d
(total M_J conservation). This conserves M_J SUM but DOES NOT conserve the
parity of the count $N_{m_j < 0}$ in general — e.g., process
(+3/2)+(-1/2) → (+1/2)+(+1/2) conserves M_J = 1 but flips N_{m_j<0} by 1.

The question: does this happen at non-zero amplitude in the chemistry
construction? Or do additional constraints kill it?
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass

from openfermion import QubitOperator
from openfermion.utils import commutator

from geovac.molecular_spec import (
    lih_spec_relativistic, beh_spec_relativistic, cah_spec_relativistic,
)
from geovac.composed_qubit_relativistic import (
    build_composed_hamiltonian_relativistic,
)
from geovac.relativistic_tapering import (
    enumerate_relativistic_orbital_table,
)


def _z_string_from_indices(qubit_indices) -> QubitOperator:
    op = QubitOperator(())
    for q in qubit_indices:
        op *= QubitOperator(((int(q), 'Z'),))
    return op


def build_mj_parity_stabilizers(
    orbital_table: List[Tuple[Any, int, int, int]],
    mode: str = 'per_block',
) -> List[QubitOperator]:
    """For each (sub_block) or globally, emit
       P_{m_j} = prod_{q: two_m_j(q) < 0} Z_q
    over the JW qubits of orbitals with negative m_j.

    Note: this is the ORIGINAL-basis Z-string (no rotation).
    Each Dirac orbital is one JW qubit.
    """
    M = len(orbital_table)
    if mode == 'global':
        neg_qubits = [
            i for i, (_, _, _, two_mj) in enumerate(orbital_table) if two_mj < 0
        ]
        return [_z_string_from_indices(neg_qubits)] if neg_qubits else []

    from collections import OrderedDict
    sb_to_neg: 'OrderedDict[Any, List[int]]' = OrderedDict()
    for q, (sb, _, _, two_mj) in enumerate(orbital_table):
        sb_to_neg.setdefault(sb, [])
        if two_mj < 0:
            sb_to_neg[sb].append(q)
    return [
        _z_string_from_indices(qs) for qs in sb_to_neg.values() if qs
    ]


def _commutes(op_a: QubitOperator, op_b: QubitOperator,
              atol: float = 1e-10) -> Tuple[bool, float]:
    c = commutator(op_a, op_b)
    max_coef = max((abs(v) for v in c.terms.values()), default=0.0)
    return (max_coef < atol), float(max_coef)


def diagnose_one(name: str, spec) -> Dict[str, Any]:
    print(f"\n{'='*72}")
    print(f"System: {name}")
    print('=' * 72)

    t0 = time.time()
    result = build_composed_hamiltonian_relativistic(spec, verbose=False)
    qop = result['qubit_op']
    Q = result['Q']
    print(f"  Q = {Q}, N_pauli = {result['N_pauli']}  "
          f"(build wall time {time.time()-t0:.1f}s)")

    orbital_table = enumerate_relativistic_orbital_table(spec)
    n_neg = sum(1 for (_, _, _, two_mj) in orbital_table if two_mj < 0)
    print(f"  orbital_table: {len(orbital_table)} orbitals, "
          f"{n_neg} with two_m_j < 0")

    # Global
    stabs_global = build_mj_parity_stabilizers(orbital_table, mode='global')
    stabs_perblock = build_mj_parity_stabilizers(orbital_table, mode='per_block')

    audit: Dict[str, Any] = {}
    if stabs_global:
        t0 = time.time()
        ok, res = _commutes(qop, stabs_global[0])
        print(f"  global  P_{{m_j}}: residual = {res:.4e}  "
              f"({'PASS' if ok else 'FAIL'})  [{time.time()-t0:.1f}s]")
        audit['global'] = {'residual': res, 'passes': ok}

    audit['per_block'] = []
    for k, op in enumerate(stabs_perblock):
        t0 = time.time()
        ok, res = _commutes(qop, op)
        print(f"  block[{k}] P_{{m_j}}: residual = {res:.4e}  "
              f"({'PASS' if ok else 'FAIL'})  [{time.time()-t0:.1f}s]")
        audit['per_block'].append({'block': k, 'residual': res, 'passes': ok})

    # Per-block all-pass?
    pb_all = all(c['passes'] for c in audit['per_block'])
    gl_pass = audit.get('global', {}).get('passes', False)
    if pb_all and gl_pass:
        verdict = 'POSITIVE (per_block)'
    elif gl_pass:
        verdict = 'PARTIAL (global only)'
    else:
        verdict = 'NEGATIVE'
    print(f"  VERDICT: {verdict}")

    return {
        'name': name,
        'Q': Q,
        'N_pauli': result['N_pauli'],
        'n_negative_mj_orbitals': n_neg,
        'audit': audit,
        'verdict': verdict,
    }


def main():
    print('=' * 72)
    print("Sprint m_j-parity Z₂ — relativistic chemistry direct-basis audit")
    print('=' * 72)

    panel = [
        ('LiH_rel', lih_spec_relativistic(R=3.015, max_n=2)),
        ('BeH_rel', beh_spec_relativistic(R=2.538, max_n=2)),
        ('CaH_rel', cah_spec_relativistic(R=3.78, max_n=2)),
    ]
    results = []
    for name, spec in panel:
        try:
            results.append(diagnose_one(name, spec))
        except Exception as e:
            print(f"  ERROR on {name}: {e}")
            results.append({'name': name, 'error': str(e)})

    # Save
    out_path = Path(__file__).parent / 'data' / 'sprint_mj_parity_z2.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump({'panel': results}, f, indent=2, default=str)
    print(f"\nData written to {out_path}")


if __name__ == '__main__':
    main()
