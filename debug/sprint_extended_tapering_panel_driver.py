"""Sprint M-vS Symmetries — verification panel for extended tapering.

Tests `geovac/extended_tapering.py` across representative molecules.
For each system, compares:
  - Q_naive (no tapering)
  - Q_hopf (Hopf-only via existing z2_tapering.hopf_tapered_from_spec)
  - Q_extended (Hopf + ℓ-parity + atom-swap + inversion via extended)
  - N_pauli at each level
  - FCI spectrum bit-exact preservation (where feasible)

Decision gate:
  - All extended tapering candidates: Q_extended ≤ Q_hopf
  - FCI ground state preserved at <1e-10 between tapered/untapered (small systems)
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import scipy.linalg as la

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass

from openfermion.linalg import get_sparse_operator
from scipy.sparse.linalg import eigsh

from geovac.molecular_spec import (
    lih_spec, beh2_spec, h2o_spec, nh3_spec, ch4_spec, hf_spec,
)
from geovac.z2_tapering import hopf_tapered_from_spec
from geovac.extended_tapering import extended_tapered_from_spec


def count_pauli(qop) -> int:
    return len(qop.terms) - (1 if () in qop.terms else 0)


def smallest_eig(qop, n_qubits: int, k: int = 1):
    """Get the lowest k eigenvalues via sparse JW + eigsh.
    Only feasible for small systems (n_qubits ≤ ~16)."""
    sparse = get_sparse_operator(qop, n_qubits=n_qubits)
    eigs, _ = eigsh(sparse, k=k, which='SA')
    return sorted(float(e) for e in eigs)


def get_default_nuclei(canonical: str, R: Optional[float] = None) -> Optional[List[Dict]]:
    """Get default nuclei list for representative panel."""
    from geovac.balanced_coupled import (
        _get_nuclei_for_lih, _get_nuclei_for_beh2, _get_nuclei_for_h2o,
    )
    if canonical == 'LiH':
        return _get_nuclei_for_lih(None, R or 3.015)
    if canonical == 'BeH2':
        return _get_nuclei_for_beh2(R or 2.5)
    if canonical == 'H2O':
        return _get_nuclei_for_h2o()
    if canonical == 'HF':
        return [
            {'Z': 9.0, 'position': (0.0, 0.0, 0.0), 'label': 'F'},
            {'Z': 1.0, 'position': (0.0, 0.0, 1.733), 'label': 'H'},
        ]
    if canonical == 'NH3':
        # Trigonal pyramid; H at C3v positions
        import math
        R_NH = 1.913  # bohr
        beta = math.radians(106.7 / 2.0)  # half HNH
        # H positions: trigonal pyramid with C3 axis along z
        # N at origin, H_1 in xz plane, H_2, H_3 rotated by 120° around z
        return [
            {'Z': 7.0, 'position': (0.0, 0.0, 0.0), 'label': 'N'},
            {
                'Z': 1.0,
                'position': (R_NH * math.sin(beta), 0.0, R_NH * math.cos(beta)),
                'label': 'H1',
            },
            {
                'Z': 1.0,
                'position': (
                    R_NH * math.sin(beta) * math.cos(2.0 * math.pi / 3.0),
                    R_NH * math.sin(beta) * math.sin(2.0 * math.pi / 3.0),
                    R_NH * math.cos(beta),
                ),
                'label': 'H2',
            },
            {
                'Z': 1.0,
                'position': (
                    R_NH * math.sin(beta) * math.cos(4.0 * math.pi / 3.0),
                    R_NH * math.sin(beta) * math.sin(4.0 * math.pi / 3.0),
                    R_NH * math.cos(beta),
                ),
                'label': 'H3',
            },
        ]
    if canonical == 'CH4':
        # Tetrahedral; H at T_d positions
        import math
        R_CH = 2.063  # bohr
        a = R_CH / math.sqrt(3.0)
        return [
            {'Z': 6.0, 'position': (0.0, 0.0, 0.0), 'label': 'C'},
            {'Z': 1.0, 'position': ( a,  a,  a), 'label': 'H1'},
            {'Z': 1.0, 'position': ( a, -a, -a), 'label': 'H2'},
            {'Z': 1.0, 'position': (-a,  a, -a), 'label': 'H3'},
            {'Z': 1.0, 'position': (-a, -a,  a), 'label': 'H4'},
        ]
    return None


def run_one(name: str, spec, verify_fci: bool = True):
    print(f"\n{'=' * 72}")
    print(f"System: {name}")
    print(f"{'=' * 72}")

    # Hopf-only baseline (production reference)
    t0 = time.time()
    hopf_out = hopf_tapered_from_spec(spec, mode='per_block', verbose=False)
    qop_naive = hopf_out['qubit_op_naive']
    qop_hopf = hopf_out['qubit_op_tapered']
    Q_naive = hopf_out['Q_naive']
    Q_hopf = hopf_out['Q_tapered']
    n_p_naive = count_pauli(qop_naive)
    n_p_hopf = count_pauli(qop_hopf)
    print(f"  Hopf-only: Q {Q_naive} → {Q_hopf} (ΔQ = {Q_naive - Q_hopf}), "
          f"N_pauli {n_p_naive} → {n_p_hopf} "
          f"[{time.time()-t0:.1f}s]")

    # Extended tapering (Hopf + ℓ-parity + atom-swap + inversion)
    nuclei = get_default_nuclei(name)
    t0 = time.time()
    try:
        ext_out = extended_tapered_from_spec(
            spec,
            use_hopf=True, use_ell_parity=True,
            use_atom_swap=True, use_inversion=True,
            nuclei=nuclei,
            verbose=False,
        )
        qop_ext = ext_out['qubit_op_tapered']
        Q_ext = ext_out['Q_tapered']
        n_p_ext = count_pauli(qop_ext)
        kinds = ext_out['kinds_kept']
        dropped = ext_out['dropped']
        ell_kept = sum(1 for k in kinds if k.startswith('ell_parity'))
        swap_kept = sum(1 for k in kinds if k.startswith('atom_swap'))
        inv_kept = sum(1 for k in kinds if k.startswith('inversion'))
        hopf_kept = sum(1 for k in kinds if k.startswith('hopf'))
        print(f"  Extended:  Q {Q_naive} → {Q_ext} (ΔQ = {Q_naive - Q_ext}), "
              f"N_pauli {n_p_naive} → {n_p_ext} "
              f"[{time.time()-t0:.1f}s]")
        print(f"    Stabilizers kept: hopf={hopf_kept}, ell_parity={ell_kept}, "
              f"atom_swap={swap_kept}, inversion={inv_kept}")
        if dropped:
            print(f"    Stabilizers dropped (commutator audit): "
                  f"{[(k, r) for k, _, r in dropped]}")

        delta_extended_vs_hopf = Q_hopf - Q_ext
        print(f"    Extended saves {delta_extended_vs_hopf} additional qubits "
              f"vs Hopf-only baseline.")

        # FCI cross-check (small systems only)
        if verify_fci and Q_ext <= 16 and Q_naive <= 20:
            t0 = time.time()
            E_ext = smallest_eig(qop_ext, n_qubits=Q_ext, k=1)[0]
            E_naive = smallest_eig(qop_naive, n_qubits=Q_naive, k=1)[0]
            print(f"    Spectrum check: E_min(naive) = {E_naive:.10f}, "
                  f"E_min(ext) = {E_ext:.10f}, |diff| = {abs(E_ext - E_naive):.2e} "
                  f"[{time.time()-t0:.1f}s]")
        elif verify_fci:
            print(f"    Spectrum check skipped (Q_ext={Q_ext} or "
                  f"Q_naive={Q_naive} too large for sparse eigsh)")

        return {
            'name': name,
            'Q_naive': Q_naive,
            'Q_hopf': Q_hopf,
            'Q_extended': Q_ext,
            'N_pauli_naive': n_p_naive,
            'N_pauli_hopf': n_p_hopf,
            'N_pauli_extended': n_p_ext,
            'delta_Q_hopf': Q_naive - Q_hopf,
            'delta_Q_extended': Q_naive - Q_ext,
            'delta_Q_extended_vs_hopf': Q_hopf - Q_ext,
            'kinds_kept': kinds,
            'dropped': dropped,
        }
    except Exception as e:
        print(f"  Extended FAILED: {e}")
        import traceback
        traceback.print_exc()
        return {'name': name, 'error': str(e)}


def main():
    print("=" * 72)
    print("Sprint M-vS Symmetries — Extended Tapering Verification Panel")
    print("Date: 2026-06-07 (session continuation)")
    print("=" * 72)

    panel = [
        ('LiH', lih_spec()),  # 3 sub-blocks, no equivalent atoms
        ('BeH2', beh2_spec()),  # 5 sub-blocks, H_1 ↔ H_2 + inversion
        ('H2O', h2o_spec()),  # 7 sub-blocks, H_1 ↔ H_2
        ('NH3', nh3_spec()),  # 7 sub-blocks (1 core + 3 bonds × 2), 3-fold
        ('CH4', ch4_spec()),  # 9 sub-blocks, S_4 Klein V_4
        ('HF', hf_spec()),  # heteronuclear, no equivalent atoms
    ]

    results = []
    for name, spec in panel:
        try:
            results.append(run_one(name, spec))
        except Exception as e:
            print(f"  {name} FAILED at spec construction: {e}")
            import traceback
            traceback.print_exc()
            results.append({'name': name, 'error': str(e)})

    # Summary table
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"{'System':6s}  {'Q_naive':>8s}  {'Q_hopf':>8s}  {'Q_ext':>8s}  "
          f"{'ΔQ_hopf':>8s}  {'ΔQ_ext':>8s}  {'extra':>6s}  "
          f"{'N_pauli_naive':>14s}  {'N_pauli_hopf':>14s}  {'N_pauli_ext':>13s}")
    for r in results:
        if 'error' in r:
            print(f"{r['name']:6s}  ERROR: {r['error']}")
            continue
        print(f"{r['name']:6s}  {r['Q_naive']:>8d}  {r['Q_hopf']:>8d}  "
              f"{r['Q_extended']:>8d}  {r['delta_Q_hopf']:>8d}  "
              f"{r['delta_Q_extended']:>8d}  "
              f"{r['delta_Q_extended_vs_hopf']:>6d}  "
              f"{r['N_pauli_naive']:>14d}  {r['N_pauli_hopf']:>14d}  "
              f"{r['N_pauli_extended']:>13d}")

    # Save
    import json
    out_path = Path(__file__).parent / 'data' / 'sprint_extended_tapering_panel.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)

    def _jsonable(obj):
        if isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, dict):
            return {k: _jsonable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_jsonable(v) for v in obj]
        return obj

    with open(out_path, 'w') as f:
        json.dump({'panel': [_jsonable(r) for r in results]}, f, indent=2)
    print(f"\nData written to {out_path}")


if __name__ == '__main__':
    main()
