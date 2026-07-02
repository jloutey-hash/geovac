"""Sprint S2-v2 — Closed form for chi_k^{h1} via cross-cut hopping matrix rank.

We claim:
    chi_k^{h1}  =  2 * rank( (h_{pq})_{p <= k, q > k} )
       + epsilon_pure_left + epsilon_pure_right
where epsilon_pure_X are 0/1 indicators for whether pure-X-side h1 terms exist.

Reasoning (formal):
- Under JW with the standard {alpha, beta} interleaving, the one-body operator
  h1 = sum_{pq} h_{pq} a^dag_p a_q (with p, q spin-orbital indices) becomes
  a sum of Pauli strings. For p < q the term a^dag_p a_q + h.c. expands as
      (h_{pq} / 2) * [X_p Z_{p+1..q-1} X_q  +  Y_p Z_{p+1..q-1} Y_q]
  i.e. two Pauli words sharing the same Z-trail.
- At a cut k between qubits k and k+1:
    * p, q <= k:  pure-left Pauli string => contributes to M[P_L, I_R] column.
    * p, q > k:  pure-right Pauli => contributes to M[I_L, P_R] row.
    * p <= k < q (cross-cut): contributes to a NON-trivial (P_L, P_R) pair,
      with P_L = X_p Z_{p+1..k} or Y_p Z_{p+1..k}, and P_R = Z_{k+1..q-1} X_q
      or Y_q Z_{k+1..q-1} -- two left-right pairs per (p, q).
- The CROSS-CUT rank contribution: per (p, q) cross-cut pair, h1 contributes
  TWO rank-1 outer products to M (XX channel and YY channel), and these are
  independent. So the cross-cut contribution to chi_k^{h1} is
      2 * (number of independent cross-cut h_{pq} channels)
    = 2 * rank( h_{cross-cut} )
  where h_{cross-cut} is the (k+1) x (Q-k-1) sub-matrix of h with p <= k,
  q > k.
- Pure-left contributions: H_left_h1 = sum_{p,q <= k} h_{pq} a^dag_p a_q.
  This contributes a rank-1 outer product (H_left_h1, I_R) IFF H_left_h1 != 0.
- Pure-right contribution: symmetric, IFF H_right_h1 != 0.
- Constant identity: contributes 1 mode (I_L, I_R) if there is any constant.
  For our chemistry h1 with no scalar term, this is absent.

Total: chi_k^{h1} = 2 * rank(h_{cross-cut}) + 1_{LL} + 1_{RR}

Empirical test:
1. Build LiH composed Hamiltonian h1 + V_ee.
2. Extract h1 only as a QubitOperator.
3. Measure chi_k^{h1} profile.
4. Compute 2 * rank(h_{cross-cut}) for each cut k.
5. Verify the formula.
"""

from __future__ import annotations

import json
import os
import numpy as np

from openfermion import QubitOperator
from openfermion.transforms import jordan_wigner
from openfermion.ops import FermionOperator

from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import lih_spec


def operator_schmidt_rank(qop: QubitOperator, cut: int, n_qubits: int,
                          rel_thr: float = 1e-12):
    """Compute operator Schmidt rank of qop at cut, plus the M matrix."""
    coef = {}
    for pauli_tuple, c in qop.terms.items():
        left = tuple((q, op) for (q, op) in pauli_tuple if q < cut)
        right = tuple((q, op) for (q, op) in pauli_tuple if q >= cut)
        coef[(left, right)] = c
    left_keys = sorted(set(k[0] for k in coef))
    right_keys = sorted(set(k[1] for k in coef))
    if not left_keys or not right_keys:
        return 0, np.array([])
    left_idx = {k: i for i, k in enumerate(left_keys)}
    right_idx = {k: i for i, k in enumerate(right_keys)}
    M = np.zeros((len(left_keys), len(right_keys)))
    for (l, r), v in coef.items():
        M[left_idx[l], right_idx[r]] = float(v.real if hasattr(v, 'real') else v)
    sv = np.linalg.svd(M, compute_uv=False)
    rank = int(np.sum(sv > rel_thr * max(sv[0], 1e-30)))
    return rank, sv


def h1_to_qubit(h1: np.ndarray) -> QubitOperator:
    """Convert a one-body (spatial) Hamiltonian to a QubitOperator under JW
    with the standard alpha/beta interleaved ordering used by composed_qubit.
    h1 is M x M (spatial) -- expanded to 2M x 2M spin-orbital and then JW.
    """
    M = h1.shape[0]
    fop = FermionOperator()
    for p in range(M):
        for q in range(M):
            if abs(h1[p, q]) < 1e-14:
                continue
            # alpha/beta interleaved spin-orbital index: 2p, 2p+1 are p_alpha,
            # p_beta
            for s in (0, 1):
                fop += FermionOperator(f"{2*p + s}^ {2*q + s}", h1[p, q])
    qop = jordan_wigner(fop)
    qop.compress()
    return qop


def h1_spinorb_cross_cut_rank(h1: np.ndarray, cut: int) -> int:
    """Return rank of the cross-cut h1 sub-matrix in spin-orbital indexing.

    Spin-orbital index p_so = 2 * p_spat + s (s in {0, 1}). Cut k corresponds
    to spin-orbitals [0, k) on left, [k, 2M) on right.
    """
    M = h1.shape[0]
    Q = 2 * M
    # Expand h1 to spin-orbital basis: h1_so[2p+s, 2q+s] = h1[p, q] (diagonal in
    # spin). Other entries zero.
    h1_so = np.zeros((Q, Q))
    for p in range(M):
        for q in range(M):
            for s in (0, 1):
                h1_so[2*p + s, 2*q + s] = h1[p, q]
    if cut <= 0 or cut >= Q:
        return 0
    cross = h1_so[:cut, cut:]
    return int(np.linalg.matrix_rank(cross, tol=1e-12))


def main():
    spec = lih_spec()
    res = build_composed_hamiltonian(spec, verbose=False)
    h1 = res['h1']
    h1_pk = res.get('h1_pk', None)
    M = h1.shape[0]
    Q = 2 * M
    if h1_pk is not None:
        # h1 in the returned dict already includes PK if pk_in_hamiltonian=True
        h1_total = h1
    else:
        h1_total = h1
    print(f"LiH: M = {M} spatial orbitals, Q = {Q} qubits")
    print(f"h1 nonzero entries: {int(np.count_nonzero(np.abs(h1_total) > 1e-12))}")

    h1_qop = h1_to_qubit(h1_total)
    print(f"h1 as QubitOperator: {len(h1_qop.terms)} Pauli terms\n")

    # Compute chi_k^{h1} profile and compare to formula
    print(f"{'cut':>4s}  {'chi^h1':>8s}  {'cross_rk':>9s}  "
          f"{'2*cross':>8s}  {'LL':>3s}  {'RR':>3s}  {'pred':>5s}  match?")
    rows = []
    for cut in range(1, Q):
        chi_h1, sv = operator_schmidt_rank(h1_qop, cut, Q)
        cross_rk = h1_spinorb_cross_cut_rank(h1_total, cut)
        # Check pure-left / pure-right blocks for nonzero content
        h1_so = np.zeros((Q, Q))
        for p in range(M):
            for q in range(M):
                for s in (0, 1):
                    h1_so[2*p + s, 2*q + s] = h1_total[p, q]
        LL_nonzero = int(np.any(np.abs(h1_so[:cut, :cut]) > 1e-12))
        RR_nonzero = int(np.any(np.abs(h1_so[cut:, cut:]) > 1e-12))
        pred = 2 * cross_rk + LL_nonzero + RR_nonzero
        match = "OK" if pred == chi_h1 else "FAIL"
        print(f"{cut:>4d}  {chi_h1:>8d}  {cross_rk:>9d}  "
              f"{2*cross_rk:>8d}  {LL_nonzero:>3d}  {RR_nonzero:>3d}  "
              f"{pred:>5d}  {match}")
        rows.append({'cut': cut, 'chi_h1': chi_h1, 'cross_rk': cross_rk,
                     'LL_nonzero': LL_nonzero, 'RR_nonzero': RR_nonzero,
                     'pred': pred, 'match': pred == chi_h1})

    n_match = sum(1 for r in rows if r['match'])
    print(f"\nFormula match: {n_match}/{len(rows)}")

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'sprint_s2_v2_h1_contribution.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump({'M': M, 'Q': Q,
                   'h1_nonzero': int(np.count_nonzero(np.abs(h1_total) > 1e-12)),
                   'rows': rows}, f, indent=2)
    print(f"Saved to {out_path}")


if __name__ == '__main__':
    main()
