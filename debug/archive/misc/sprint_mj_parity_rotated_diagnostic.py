"""Sprint m_j-parity Z₂ rotated-basis — relativistic chemistry.

Follow-on from v3.94.0 direct-basis NEGATIVE: rotate to the
m_j → −m_j (sym, antisym) basis before constructing the Z-string,
analogous to the non-relativistic Hopf m_l → −m_l Z₂ from
`geovac/z2_tapering.py`.

Structural argument: under all-orbital m_j → −m_j flip, the four
3-j sign factors in the jj-coupled X_k angular coefficient products
multiply to +1, so the ERI tensor is symmetric under the full
reflection. The rotation to (sym, antisym) basis then isolates an
antisym sector closed under H, and the corresponding Z-string commutes.

Verdict gate:
  POSITIVE: residual < 1e-10 on rotated H_rel for LiH_rel/BeH_rel/CaH_rel.
  NEGATIVE: residual > 1e-6. Document mechanism.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import scipy.linalg as la

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass

from openfermion import FermionOperator, QubitOperator, jordan_wigner
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


def build_mj_pm_rotation(
    orbital_table: List[Tuple[Any, int, int, int]],
) -> Tuple[np.ndarray, np.ndarray]:
    """Build the rotation to the m_j → −m_j (sym, antisym) eigenbasis.

    For each (sb_key, n_fock, kappa) shell, pair orbitals
    (two_m_j = +x, two_m_j = -x). Rotate to:
        φ_+ = (φ_{+x} + φ_{-x}) / √2  (parity +1)
        φ_- = (φ_{+x} - φ_{-x}) / √2  (parity −1)
    """
    M = len(orbital_table)
    U = np.zeros((M, M))
    parity = np.zeros(M, dtype=int)

    idx_of = {orb: i for i, orb in enumerate(orbital_table)}
    handled = [False] * M

    for i, (sb, n_fock, kappa, two_m_j) in enumerate(orbital_table):
        if handled[i]:
            continue
        if two_m_j == 0:
            U[i, i] = 1.0
            parity[i] = +1
            handled[i] = True
            continue

        partner_key = (sb, n_fock, kappa, -two_m_j)
        if partner_key not in idx_of:
            raise RuntimeError(
                f"orbital {orbital_table[i]} has no m_j → −m_j partner"
            )
        j = idx_of[partner_key]
        if j == i:
            U[i, i] = 1.0
            parity[i] = +1
            handled[i] = True
            continue

        if two_m_j > 0:
            i_plus, i_minus = i, j
        else:
            i_plus, i_minus = j, i

        inv_sqrt2 = 1.0 / np.sqrt(2.0)
        U[i_plus, i_plus] = inv_sqrt2
        U[i_plus, i_minus] = inv_sqrt2
        U[i_minus, i_plus] = inv_sqrt2
        U[i_minus, i_minus] = -inv_sqrt2
        parity[i_plus] = +1
        parity[i_minus] = -1
        handled[i_plus] = True
        handled[i_minus] = True

    err = float(np.max(np.abs(U @ U.T - np.eye(M))))
    if err > 1e-12:
        raise RuntimeError(f"m_j rotation not orthogonal: err={err:.2e}")

    return U, parity


def densify_eri(eri_sparse: Dict, Q: int) -> np.ndarray:
    """Convert sparse ERI dict to dense Q⁴ tensor."""
    eri = np.zeros((Q, Q, Q, Q))
    for (a, b, c, d), val in eri_sparse.items():
        eri[a, b, c, d] = val
    return eri


def symmetrize_eri_for_hermiticity(eri: np.ndarray) -> np.ndarray:
    """Mirror the symmetrization in build_composed_hamiltonian_relativistic:
    sym_eri[a,b,c,d] = 0.5 * (eri[a,b,c,d] + eri[c,d,a,b])
    """
    return 0.5 * (eri + eri.transpose(2, 3, 0, 1))


def rotate_h1_eri(
    h1: np.ndarray, eri: np.ndarray, U: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """h1' = U h1 U^T; eri' physicist <ab|V|cd> covariantly transformed."""
    h1_rot = U @ h1 @ U.T
    # Physicist convention: ⟨ab|V|cd⟩.  Rotation acts on all four orbital
    # slots: pre-rotation orbitals → post-rotation orbitals via U_{p,a}.
    eri_rot = np.einsum(
        'pa,qb,rc,sd,abcd->pqrs', U, U, U, U, eri, optimize='optimal',
    )
    return h1_rot, eri_rot


def build_fermion_op_from_rotated_integrals(
    h1: np.ndarray, eri_sym: np.ndarray, nuc_rep: float,
) -> FermionOperator:
    """Mirror the FermionOperator assembly in
    build_composed_hamiltonian_relativistic.

    One-body: h1[p,q] a_p^† a_q
    Two-body: physicist <a b | V | c d> → 0.5 a^† b^† d c
    """
    Q = h1.shape[0]
    fop = FermionOperator((), float(nuc_rep))
    for p in range(Q):
        for q in range(Q):
            hv = h1[p, q]
            if abs(hv) < 1e-12:
                continue
            fop += FermionOperator(((p, 1), (q, 0)), hv)

    for a in range(Q):
        for b in range(Q):
            if a == b:
                continue
            for c in range(Q):
                for d in range(Q):
                    if c == d:
                        continue
                    val = eri_sym[a, b, c, d]
                    if abs(val) < 1e-14:
                        continue
                    fop += FermionOperator(
                        ((a, 1), (b, 1), (d, 0), (c, 0)),
                        0.5 * val,
                    )
    return fop


def build_mj_parity_stabilizers_rotated(
    parity: np.ndarray, orbital_table, mode: str = 'per_block',
) -> List[QubitOperator]:
    """Z-string on antisym orbitals (parity = -1)."""
    M = len(parity)
    if mode == 'global':
        anti = [i for i in range(M) if parity[i] == -1]
        return [_z_string_from_indices(anti)] if anti else []

    from collections import OrderedDict
    sb_to_anti: 'OrderedDict[Any, List[int]]' = OrderedDict()
    for p, (sb, _, _, _) in enumerate(orbital_table):
        sb_to_anti.setdefault(sb, [])
        if parity[p] == -1:
            sb_to_anti[sb].append(p)
    return [_z_string_from_indices(anti) for anti in sb_to_anti.values() if anti]


def _commutes(op_a, op_b, atol: float = 1e-10):
    c = commutator(op_a, op_b)
    max_coef = max((abs(v) for v in c.terms.values()), default=0.0)
    return (max_coef < atol), float(max_coef)


def diagnose_one(name: str, spec) -> Dict[str, Any]:
    print(f"\n{'='*72}")
    print(f"System: {name}")
    print('=' * 72)

    t0 = time.time()
    result = build_composed_hamiltonian_relativistic(spec, verbose=False)
    Q = result['Q']
    h1 = result['h1']
    eri_sparse = result['eri_sparse']
    nuc_rep = result['nuclear_repulsion']
    print(f"  Q = {Q}, original N_pauli = {result['N_pauli']}  "
          f"(builder {time.time()-t0:.1f}s)")

    orbital_table = enumerate_relativistic_orbital_table(spec)
    assert len(orbital_table) == Q

    # Densify ERI and symmetrize (mirroring builder)
    print(f"  densifying ERI ({len(eri_sparse)} sparse → {Q**4} dense)...")
    t0 = time.time()
    eri = densify_eri(eri_sparse, Q)
    eri_sym = symmetrize_eri_for_hermiticity(eri)
    print(f"  ({time.time()-t0:.1f}s)")

    # Verify: rebuild fermion_op from dense, JW, compare to builder's qop
    print(f"  rebuilding fermion_op from dense h1+eri_sym + JW...")
    t0 = time.time()
    fop_rebuilt = build_fermion_op_from_rotated_integrals(h1, eri_sym, nuc_rep)
    qop_rebuilt = jordan_wigner(fop_rebuilt)
    print(f"  ({time.time()-t0:.1f}s)")
    qop_orig = result['qubit_op']
    diff_op = qop_orig - qop_rebuilt
    max_diff = max(abs(v) for v in diff_op.terms.values()) if diff_op.terms else 0.0
    print(f"  ||qop_orig - qop_rebuilt||_max = {max_diff:.4e}  "
          f"({'PASS' if max_diff < 1e-8 else 'FAIL — densification convention mismatch'})")

    # Build rotation
    print(f"  building m_j → −m_j rotation U_mj...")
    U, parity = build_mj_pm_rotation(orbital_table)
    n_plus = int(np.sum(parity == +1))
    n_minus = int(np.sum(parity == -1))
    print(f"  parity counts: +1={n_plus}, −1={n_minus}")

    # Rotate
    t0 = time.time()
    h1_rot, eri_rot = rotate_h1_eri(h1, eri_sym, U)
    eri_rot_sym = symmetrize_eri_for_hermiticity(eri_rot)
    print(f"  rotation wall time {time.time()-t0:.1f}s")

    # Build rotated fermion op + JW
    t0 = time.time()
    fop_rot = build_fermion_op_from_rotated_integrals(
        h1_rot, eri_rot_sym, nuc_rep,
    )
    qop_rot = jordan_wigner(fop_rot)
    N_pauli_rot = len(qop_rot.terms) - (1 if () in qop_rot.terms else 0)
    print(f"  N_pauli_rotated = {N_pauli_rot}  ({time.time()-t0:.1f}s)")

    # Build stabilizers
    stabs_global = build_mj_parity_stabilizers_rotated(
        parity, orbital_table, mode='global',
    )
    stabs_perblock = build_mj_parity_stabilizers_rotated(
        parity, orbital_table, mode='per_block',
    )
    print(f"  stabilizers: global={len(stabs_global)}, "
          f"per_block={len(stabs_perblock)}")

    # Audit
    audit: Dict[str, Any] = {'per_block': []}
    if stabs_global:
        ok, res = _commutes(qop_rot, stabs_global[0])
        print(f"  global  P_{{m_j}}^rot: residual = {res:.4e}  "
              f"({'PASS' if ok else 'FAIL'})")
        audit['global'] = {'residual': res, 'passes': ok}
    for k, op in enumerate(stabs_perblock):
        ok, res = _commutes(qop_rot, op)
        print(f"  block[{k}] P_{{m_j}}^rot: residual = {res:.4e}  "
              f"({'PASS' if ok else 'FAIL'})")
        audit['per_block'].append({
            'block': k, 'residual': res, 'passes': ok,
        })

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
        'N_pauli_orig': result['N_pauli'],
        'N_pauli_rotated': N_pauli_rot,
        'densification_check_residual': max_diff,
        'parity_plus': n_plus, 'parity_minus': n_minus,
        'audit': audit,
        'verdict': verdict,
    }


def main():
    print('=' * 72)
    print("Sprint m_j-parity Z₂ ROTATED-BASIS — relativistic chemistry audit")
    print('Date: 2026-06-08')
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
            import traceback
            traceback.print_exc()
            results.append({'name': name, 'error': str(e)})

    out_path = Path(__file__).parent / 'data' / 'sprint_mj_parity_rotated.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump({'panel': results}, f, indent=2, default=str)
    print(f"\nData written to {out_path}")


if __name__ == '__main__':
    main()
