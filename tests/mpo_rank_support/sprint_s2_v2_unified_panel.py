"""Sprint S2-v2 — Unified verification panel for Theorem 3.2.A.unified.

For every cut k in LiH, report:
  chi^H       -- empirical MPO bond rank of the full Hamiltonian
  chi^h1      -- bond rank of the one-body piece alone
  chi^Vee     -- bond rank of the two-body piece alone
  chi^L=0,1,2 -- per-multipole-L pieces
  bound_h1    -- 2 * rank(h_cross-cut) + 1_LL + 1_RR  (closed form, Thm A)
  bound_sub   -- chi^h1 + chi^Vee (subadditivity)
  sum_L       -- sum over L of chi^L (per-L subadditivity)

Three rank inheritance corollaries this panel verifies:
  (1) chi^h1 == bound_h1 at every cut  (EXACT closed form, Theorem 3.2.A.A)
  (2) chi^H <= chi^h1 + chi^Vee at every cut (subadditivity, generically slack)
  (3) chi^Vee <= sum_L chi^L (per-L subadditivity, generically slack from L-mixing)
  (4) Profile dominated by chi^L=0 when 1s/2s pair is cross-cut active,
      then transitions to chi^L=1 when 2s closes -- structural inheritance
      from the multipole density-operator construction.
"""

from __future__ import annotations

import json
import os
from collections import defaultdict
import numpy as np

from openfermion import QubitOperator
from openfermion.transforms import jordan_wigner
from openfermion.ops import FermionOperator

from geovac.composed_qubit import (
    build_composed_hamiltonian, _ck_coefficient, _enumerate_states,
)
from geovac.molecular_spec import lih_spec


def operator_schmidt_rank(qop: QubitOperator, cut: int, n_qubits: int,
                          rel_thr: float = 1e-12) -> int:
    coef = {}
    for pauli_tuple, c in qop.terms.items():
        left = tuple((q, op) for (q, op) in pauli_tuple if q < cut)
        right = tuple((q, op) for (q, op) in pauli_tuple if q >= cut)
        coef[(left, right)] = c
    left_keys = sorted(set(k[0] for k in coef))
    right_keys = sorted(set(k[1] for k in coef))
    if not left_keys or not right_keys:
        return 0
    left_idx = {k: i for i, k in enumerate(left_keys)}
    right_idx = {k: i for i, k in enumerate(right_keys)}
    M = np.zeros((len(left_keys), len(right_keys)))
    for (l, r), v in coef.items():
        M[left_idx[l], right_idx[r]] = float(v.real if hasattr(v, 'real') else v)
    sv = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(sv > rel_thr * max(sv[0], 1e-30)))


def h1_to_qubit(h1: np.ndarray) -> QubitOperator:
    M = h1.shape[0]
    fop = FermionOperator()
    for p in range(M):
        for q in range(M):
            if abs(h1[p, q]) < 1e-14:
                continue
            for s in (0, 1):
                fop += FermionOperator(f"{2*p+s}^ {2*q+s}", h1[p, q])
    qop = jordan_wigner(fop)
    qop.compress()
    return qop


def h1_cross_cut_rank(h1: np.ndarray, cut: int):
    M = h1.shape[0]
    Q = 2 * M
    h1_so = np.zeros((Q, Q))
    for p in range(M):
        for q in range(M):
            for s in (0, 1):
                h1_so[2*p+s, 2*q+s] = h1[p, q]
    cross = h1_so[:cut, cut:]
    LL = int(np.any(np.abs(h1_so[:cut, :cut]) > 1e-12))
    RR = int(np.any(np.abs(h1_so[cut:, cut:]) > 1e-12))
    return int(np.linalg.matrix_rank(cross, tol=1e-12)), LL, RR


def build_Vee_total_qubit(eri: np.ndarray) -> QubitOperator:
    M = eri.shape[0]
    fop = FermionOperator()
    for p in range(M):
        for q in range(M):
            for r in range(M):
                for s in range(M):
                    v = eri[p, r, q, s]
                    if abs(v) < 1e-15:
                        continue
                    for s1 in (0, 1):
                        for s2 in (0, 1):
                            fop += FermionOperator(
                                f"{2*p+s1}^ {2*q+s2}^ {2*s+s2} {2*r+s1}",
                                0.5 * v)
    qop = jordan_wigner(fop)
    qop.compress()
    return qop


def build_Vee_per_L_qubit(eri, states_by_block, block_offsets):
    Vee_L = defaultdict(FermionOperator)
    for blk_idx, (states, offset) in enumerate(
            zip(states_by_block, block_offsets)):
        n_sp = len(states)
        ck = {}
        for a in range(n_sp):
            la, ma = states[a][1], states[a][2]
            for c in range(n_sp):
                lc, mc = states[c][1], states[c][2]
                for k in range(0, la + lc + 1):
                    if (la + lc + k) % 2 != 0:
                        continue
                    v = _ck_coefficient(la, ma, lc, mc, k)
                    if abs(v) > 1e-15:
                        ck[(a, c, k)] = v
        ac_k_map = defaultdict(list)
        for (a, c, k), v in ck.items():
            ac_k_map[(a, c)].append((k, v))
        for a in range(n_sp):
            for b in range(n_sp):
                for c in range(n_sp):
                    for d in range(n_sp):
                        if states[a][2] + states[b][2] != states[c][2] + states[d][2]:
                            continue
                        ap, bp = offset + a, offset + b
                        cp, dp = offset + c, offset + d
                        v_full = float(eri[ap, cp, bp, dp])
                        if abs(v_full) < 1e-15:
                            continue
                        ks_ac = set(k for k, _ in ac_k_map.get((a, c), []))
                        ks_bd = set(k for k, _ in ac_k_map.get((b, d), []))
                        ks = sorted(ks_ac & ks_bd)
                        if not ks:
                            continue
                        if len(ks) == 1:
                            L = ks[0]
                            for s1 in (0, 1):
                                for s2 in (0, 1):
                                    Vee_L[L] += FermionOperator(
                                        f"{2*ap+s1}^ {2*bp+s2}^ "
                                        f"{2*dp+s2} {2*cp+s1}", 0.5 * v_full)
                        else:
                            ck_pairs = {k: (dict(ac_k_map[(a,c)])[k]
                                            * dict(ac_k_map[(b,d)])[k])
                                        for k in ks}
                            total = sum(v*v for v in ck_pairs.values())
                            for k in ks:
                                frac = ck_pairs[k]**2 / total
                                for s1 in (0, 1):
                                    for s2 in (0, 1):
                                        Vee_L[k] += FermionOperator(
                                            f"{2*ap+s1}^ {2*bp+s2}^ "
                                            f"{2*dp+s2} {2*cp+s1}",
                                            0.5 * v_full * frac)
    out = {}
    for L, fop in Vee_L.items():
        qop = jordan_wigner(fop)
        qop.compress()
        out[L] = qop
    return out


def main():
    spec = lih_spec()
    res = build_composed_hamiltonian(spec, verbose=False)
    h1 = res['h1']
    eri = res['eri']
    H_qop = res['qubit_op']
    M = h1.shape[0]
    Q = 2 * M

    states_by_block = []
    block_offsets = []
    offset = 0
    for blk in spec.blocks:
        l_min = getattr(blk, 'l_min', 0)
        cs = _enumerate_states(blk.max_n, l_min=l_min)
        states_by_block.append(cs); block_offsets.append(offset); offset += len(cs)
        if blk.has_h_partner:
            pn = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            ps = _enumerate_states(pn)
            states_by_block.append(ps); block_offsets.append(offset); offset += len(ps)

    h1_qop = h1_to_qubit(h1)
    Vee_qop = build_Vee_total_qubit(eri)
    Vee_L = build_Vee_per_L_qubit(eri, states_by_block, block_offsets)

    print(f"M={M}, Q={Q}; sub-block offsets={block_offsets}")
    print(f"Pauli counts: H={len(H_qop.terms)} h1={len(h1_qop.terms)} "
          f"Vee={len(Vee_qop.terms)} {[f'L={L}:{len(Vee_L[L].terms)}' for L in sorted(Vee_L)]}")
    print()

    print(f"{'cut':>4} {'chi^H':>5} {'chi^h1':>6} {'bnd_h1':>6} "
          f"{'chi^V':>6} {'L=0':>4} {'L=1':>4} {'L=2':>4} "
          f"{'sumL':>5} {'sub':>5} {'h1+V':>5} {'lhs<=rhs':>9}")
    print("-" * 90)
    rows = []
    for cut in range(1, Q):
        chi_H = operator_schmidt_rank(H_qop, cut, Q)
        chi_h1 = operator_schmidt_rank(h1_qop, cut, Q)
        chi_Vee = operator_schmidt_rank(Vee_qop, cut, Q)
        rk, LL, RR = h1_cross_cut_rank(h1, cut)
        bnd_h1 = 2 * rk + LL + RR
        chi_per_L = {L: operator_schmidt_rank(Vee_L[L], cut, Q) for L in Vee_L}
        sum_L = sum(chi_per_L.values())
        sub_h1V = chi_h1 + chi_Vee

        # Validate three inequalities
        eq_h1 = chi_h1 == bnd_h1
        ineq_sub = chi_H <= sub_h1V
        ineq_perL = chi_Vee <= sum_L
        ok_chr = "*" if (eq_h1 and ineq_sub and ineq_perL) else "X"

        print(f"{cut:>4} {chi_H:>5} {chi_h1:>6} {bnd_h1:>6} "
              f"{chi_Vee:>6} {chi_per_L.get(0, 0):>4} {chi_per_L.get(1, 0):>4} "
              f"{chi_per_L.get(2, 0):>4} {sum_L:>5} {sub_h1V:>5} "
              f"{sub_h1V:>5} {ok_chr:>9}")
        rows.append({
            'cut': cut, 'chi_H': chi_H, 'chi_h1': chi_h1, 'bnd_h1': bnd_h1,
            'chi_Vee': chi_Vee, 'chi_per_L': chi_per_L,
            'sum_L': sum_L, 'sub_h1V': sub_h1V,
            'eq_h1': eq_h1, 'ineq_sub': ineq_sub, 'ineq_perL': ineq_perL,
        })

    n_eq_h1 = sum(1 for r in rows if r['eq_h1'])
    n_ineq_sub = sum(1 for r in rows if r['ineq_sub'])
    n_ineq_perL = sum(1 for r in rows if r['ineq_perL'])
    print(f"\nTheorem 3.2.A.unified verification:")
    print(f"  (A) chi^h1 == 2*rank(h_cross) + 1_LL + 1_RR: {n_eq_h1}/{len(rows)}")
    print(f"  (B) chi^H <= chi^h1 + chi^Vee (subadditivity): {n_ineq_sub}/{len(rows)}")
    print(f"  (C) chi^Vee <= sum_L chi^L (per-L subadditivity): {n_ineq_perL}/{len(rows)}")

    # Compute average slack for (B) and (C)
    slack_B = [r['sub_h1V'] - r['chi_H'] for r in rows]
    slack_C = [r['sum_L'] - r['chi_Vee'] for r in rows]
    print(f"  (B) average slack: {np.mean(slack_B):.2f} (max {max(slack_B)})")
    print(f"  (C) average slack: {np.mean(slack_C):.2f} (max {max(slack_C)})")

    # Universal interior profile check
    print(f"\nUniversal interior profile within each sub-block:")
    print(f"  Expected: [4, 16, 16, 9, 9, 9, 6, 3, 3, 2]  (cuts 1..10)")
    interior_cuts = [(1,10), (11,20), (21,30)]
    for (lo, hi) in interior_cuts:
        prof = [r['chi_H'] for r in rows if lo <= r['cut'] < hi]
        print(f"  sub-block cuts [{lo}, {hi}): chi_H profile = {prof}")

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'sprint_s2_v2_unified_panel.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump({
            'M': M, 'Q': Q,
            'sub_block_offsets': block_offsets,
            'rows': rows,
            'verification': {
                'eq_h1_count': f'{n_eq_h1}/{len(rows)}',
                'ineq_sub_count': f'{n_ineq_sub}/{len(rows)}',
                'ineq_perL_count': f'{n_ineq_perL}/{len(rows)}',
                'mean_slack_B': float(np.mean(slack_B)),
                'mean_slack_C': float(np.mean(slack_C)),
            },
        }, f, indent=2, default=str)
    print(f"\nSaved to {out_path}")


if __name__ == '__main__':
    main()
