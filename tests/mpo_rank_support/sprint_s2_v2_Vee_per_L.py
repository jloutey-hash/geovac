"""Sprint S2-v2 — Per-L decomposition of chi_k^{V_ee}.

We build V_ee = sum_L V_ee^L using the Slater-Gaunt structure of the
composed builder, then compute chi_k^{V_ee, L} (the operator Schmidt rank
of V_ee^L alone) at every cut k. Compare to:

  chi_k^{V_ee, total}  vs  sum_L chi_k^{V_ee, L}

The bound chi(A+B) <= chi(A) + chi(B) means the SUM provides an upper bound
on chi^{V_ee, total}. Where the bound is tight, the per-L decomposition is
the right basis for the closed form. Where the bound has slack, the
multi-L mixing creates rank reduction (analogous to F3 graded compression).

We then compare per-cut to chi^H total and to the universal {4, 16, 16, 9,
9, 9, 6, 3, 3, 2} profile to attribute each row to a specific cause.
"""

from __future__ import annotations

import json
import os
from collections import defaultdict
from typing import Dict, Tuple, List

import numpy as np

from openfermion import QubitOperator
from openfermion.transforms import jordan_wigner
from openfermion.ops import FermionOperator

from geovac.composed_qubit import (
    build_composed_hamiltonian,
    _ck_coefficient,
    _enumerate_states,
)
from geovac.molecular_spec import lih_spec


def operator_schmidt_rank(qop: QubitOperator, cut: int, n_qubits: int,
                          rel_thr: float = 1e-12):
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


def build_Vee_per_L_qubit(eri_full: np.ndarray, M_sp: int,
                          states_by_block: List[List[Tuple[int, int, int]]],
                          block_offsets: List[int],
                          ) -> Dict[int, QubitOperator]:
    """Split V_ee tensor by multipole order L using c^L(la,ma;lc,mc).

    Strategy: for each block (sub-block in the composed builder), use the
    pre-existing Slater-Condon decomposition

        V^L_{abcd}  =  sum_{k = L} c^L(la,ma; lc,mc) c^L(lb,mb; ld,md) R^L

    We can't easily re-extract R^L from the ERI tensor alone (would need to
    invert the Gaunt structure), but we can EMPIRICALLY decompose V_ee by L
    by FILTERING the ERI tensor to those (a, b, c, d) where the dominant
    L-channel is L_max(a,c) = max(|la-lc|, |lb-ld|, ...).

    Easier approach: just enumerate the V_ee tensor and PROJECT onto each L
    via the Slater structure: V_ee^L_{abcd} keeps only terms where
    L = k_ac = k_bd in the original Slater sum.

    For composed builder block-diagonal V_ee, each (a, b, c, d) lies in a
    single block. Within a block, V^k_{abcd} can be extracted by reading
    c^k coefficients directly.
    """
    # We'll reconstruct V_ee^L block-by-block, since the composed builder
    # only has within-block ERIs (F4 from DF memo).
    Vee_L: Dict[int, FermionOperator] = defaultdict(FermionOperator)

    for blk_idx, (states, offset) in enumerate(
            zip(states_by_block, block_offsets)):
        n_sp = len(states)

        # Per-block c^k table
        ck = {}
        for a in range(n_sp):
            la, ma = states[a][1], states[a][2]
            for c in range(n_sp):
                lc, mc = states[c][1], states[c][2]
                k_max = la + lc
                for k in range(0, k_max + 1):
                    if (la + lc + k) % 2 != 0:
                        continue
                    v = _ck_coefficient(la, ma, lc, mc, k)
                    if abs(v) > 1e-15:
                        ck[(a, c, k)] = v

        # For each k, build the L=k component of the ERI from the c^k
        # decomposition; we need R^k itself which lives in eri_full only via
        # its sum. We'll use the LINEAR algebra trick:
        #
        #   eri_full_{abcd} = sum_k c^k(a,c) c^k(b,d) R^k_{a,c,b,d}
        #
        # So if at fixed (a, b, c, d) there's only ONE allowed k_ac=k_bd, the
        # whole entry IS V^k. If multiple k are allowed, we split using the
        # known c^k weights (the radial part R^k is the same).
        #
        # Practical method: for entries with multiple allowed L, distribute
        # by least-squares solving the c^k mixing. For our atomic {1s, 2s,
        # 2p} basis, at most 1-2 L channels mix per quartet.
        #
        # Even simpler: rebuild V^k DIRECTLY from a re-solve with R^k
        # available. But R^k requires re-running the integrals. Skip that --
        # use the unique-L identification per quartet when possible.
        ac_k_map = defaultdict(list)
        for (a, c, k), v in ck.items():
            ac_k_map[(a, c)].append((k, v))

        # For each (a, b, c, d) in the within-block ERI:
        ks_seen = set()
        # Reconstruct V^k by inverting Slater per-quartet:
        # eri[a, b, c, d] = sum_k_allowed c^k(a,c) c^k(b,d) R^k
        # If only one k is allowed for (a, c) AND (b, d), then
        #   V^k_{abcd} = eri[a, b, c, d]
        # If multiple k allowed, the same eri entry has contributions from
        # different k. We can decompose by SOLVING a system over k.
        #
        # We collect all (a, b, c, d) quartets that share the same k-list
        # and then look at how they couple.
        for a in range(n_sp):
            for b in range(n_sp):
                for c in range(n_sp):
                    for d in range(n_sp):
                        if states[a][2] + states[b][2] != \
                                states[c][2] + states[d][2]:
                            continue
                        # Spatial index in the full builder offset
                        ap = offset + a
                        bp = offset + b
                        cp = offset + c
                        dp = offset + d
                        val_full = float(eri_full[ap, cp, bp, dp])
                        if abs(val_full) < 1e-15:
                            continue
                        # Allowed k: intersection of (a, c) k's and (b, d) k's
                        ks_ac = set(k for k, _ in ac_k_map.get((a, c), []))
                        ks_bd = set(k for k, _ in ac_k_map.get((b, d), []))
                        ks = sorted(ks_ac & ks_bd)
                        if len(ks) == 0:
                            continue
                        if len(ks) == 1:
                            # Single-L quartet: all goes to one L
                            L = ks[0]
                            ks_seen.add(L)
                            # FermionOperator entry, with spin sum and 1/2
                            for s1 in (0, 1):
                                for s2 in (0, 1):
                                    fo = FermionOperator(
                                        f"{2*ap + s1}^ {2*bp + s2}^ "
                                        f"{2*dp + s2} {2*cp + s1}",
                                        0.5 * val_full,
                                    )
                                    Vee_L[L] += fo
                        else:
                            # Multi-L quartet: need to split the val by
                            # c^k(a,c) c^k(b,d) weights. We can solve a
                            # small linear system, but at this point each
                            # quartet gives ONE number from the ERI tensor;
                            # to split we'd need an R^k cache.
                            #
                            # For now: attribute proportionally to the
                            # |c^k(a,c) c^k(b,d)|^2 weight as an approximation.
                            # This is an ANSATZ that overestimates per-L
                            # rank but is structurally consistent.
                            ck_pairs = {
                                k: (
                                    dict(ac_k_map[(a, c)])[k]
                                    * dict(ac_k_map[(b, d)])[k]
                                )
                                for k in ks
                            }
                            total_w2 = sum(v*v for v in ck_pairs.values())
                            for k in ks:
                                w = ck_pairs[k]
                                frac = w * w / total_w2
                                ks_seen.add(k)
                                for s1 in (0, 1):
                                    for s2 in (0, 1):
                                        fo = FermionOperator(
                                            f"{2*ap + s1}^ {2*bp + s2}^ "
                                            f"{2*dp + s2} {2*cp + s1}",
                                            0.5 * val_full * frac,
                                        )
                                        Vee_L[k] += fo
        print(f"  Block {blk_idx}: Ls seen = {sorted(ks_seen)}")

    # JW each
    Vee_L_qubit = {}
    for L, fop in Vee_L.items():
        qop = jordan_wigner(fop)
        qop.compress()
        Vee_L_qubit[L] = qop
    return Vee_L_qubit


def main():
    spec = lih_spec()
    res = build_composed_hamiltonian(spec, verbose=False)
    h1 = res['h1']
    eri = res['eri']
    M_sp = h1.shape[0]
    Q = 2 * M_sp

    # Enumerate states per sub-block (need to know which orbitals belong to
    # which block to apply per-block c^k decomposition).
    states_by_block = []
    block_offsets = []
    offset = 0
    for blk in spec.blocks:
        l_min = getattr(blk, 'l_min', 0)
        center_states = _enumerate_states(blk.max_n, l_min=l_min)
        states_by_block.append(center_states)
        block_offsets.append(offset)
        offset += len(center_states)
        if blk.has_h_partner:
            partner_max_n = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_states = _enumerate_states(partner_max_n)
            states_by_block.append(partner_states)
            block_offsets.append(offset)
            offset += len(partner_states)

    print(f"Sub-blocks: {len(states_by_block)}; offsets: {block_offsets}; "
          f"M = {M_sp}, Q = {Q}")
    print(f"ERI nonzero entries: {int(np.count_nonzero(np.abs(eri) > 1e-15))}")

    # Build V_ee total as qubit op
    Vee_fop = FermionOperator()
    for p in range(M_sp):
        for q in range(M_sp):
            for r in range(M_sp):
                for s in range(M_sp):
                    v = eri[p, r, q, s]  # physicist convention <pq|rs>
                    if abs(v) < 1e-15:
                        continue
                    for s1 in (0, 1):
                        for s2 in (0, 1):
                            Vee_fop += FermionOperator(
                                f"{2*p + s1}^ {2*q + s2}^ "
                                f"{2*s + s2} {2*r + s1}",
                                0.5 * v,
                            )
    Vee_qop = jordan_wigner(Vee_fop)
    Vee_qop.compress()
    print(f"V_ee total Pauli terms: {len(Vee_qop.terms)}")

    print(f"\nBuilding per-L V_ee...")
    Vee_L = build_Vee_per_L_qubit(eri, M_sp, states_by_block, block_offsets)
    for L, qop in sorted(Vee_L.items()):
        print(f"  V_ee^L={L}: {len(qop.terms)} Pauli terms")

    # Compute chi profiles
    print(f"\n{'cut':>4s}  {'chi^Vee':>8s}  ", end='')
    for L in sorted(Vee_L.keys()):
        print(f"chi^L={L:1d} ", end='')
    print(f" {'sum_L':>6s}  match?")

    rows = []
    for cut in range(1, Q):
        chi_Vee, _ = operator_schmidt_rank(Vee_qop, cut, Q)
        chi_per_L = {L: operator_schmidt_rank(Vee_L[L], cut, Q)[0]
                     for L in Vee_L}
        sum_L = sum(chi_per_L.values())
        match = "tight" if chi_Vee == sum_L else (
            "slack" if chi_Vee < sum_L else "FAIL")
        print(f"{cut:>4d}  {chi_Vee:>8d}  ", end='')
        for L in sorted(Vee_L.keys()):
            print(f"  {chi_per_L[L]:>4d}  ", end='')
        print(f" {sum_L:>6d}  {match}")
        rows.append({'cut': cut, 'chi_Vee': chi_Vee,
                     **{f'chi_L{L}': chi_per_L[L] for L in Vee_L},
                     'sum_L': sum_L, 'slack': sum_L - chi_Vee})

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'sprint_s2_v2_Vee_per_L.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump({'M_sp': M_sp, 'Q': Q, 'Ls': sorted(Vee_L.keys()),
                   'n_pauli_Vee': len(Vee_qop.terms),
                   'n_pauli_per_L': {str(L): len(qop.terms)
                                     for L, qop in Vee_L.items()},
                   'rows': rows}, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == '__main__':
    main()
