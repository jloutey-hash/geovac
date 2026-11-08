"""S2 closed-form: derive chi_k from chemistry-operator counting.

Goal: find a structural formula chi_k = f(orbital_basis, m_closed) that
matches the empirical post-orbital sequence (1, 16, 9, 9, 3, 2) for the
{1s, 2s, 2p_{-1,0,+1}} valence basis.

Approach: at each post-orbital cut, enumerate the *fermionic operators*
(not Pauli strings) that contribute to the Hamiltonian with mixed
left-right support. The matrix rank chi_k counts independent
left-fermionic-operators (with spin) modulo linear combinations.

Strategy:
  1. For each orbital m_closed (= number of spatial orbitals on left),
     enumerate all (p, q, r, s) ERI quartets with Gaunt nonzero such that
     some of {p, q, r, s} are on left and some on right.
  2. Classify by ERI type (1L3R, 2L2R, 3L1R) and count.
  3. Map to "independent left fermionic operators."
  4. Add 2 (additive baseline) and compare to chi_k.
"""

from __future__ import annotations

from itertools import combinations, product
from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np

from geovac.composed_qubit import _ck_coefficient, _enumerate_states


def enumerate_orbitals_5(n_max: int = 2, l_min: int = 0):
    """Return the 5 orbitals (1s, 2s, 2p_-1, 2p_0, 2p_+1) for the standard sub-block."""
    return _enumerate_states(n_max, l_min=l_min)


def count_cross_cut_eri(orbitals, m_closed: int) -> Dict:
    """For each ERI quartet (p, q, r, s) with Gaunt-allowed coupling,
    classify by left/right split.

    Convention: orbital index < m_closed -> LEFT, otherwise -> RIGHT.

    Each quartet (p, q, r, s) represents the chemist ERI (pr|qs) with
    g[p, r, q, s] potentially nonzero if Gaunt-allowed at some L.
    """
    M = len(orbitals)
    counts = {
        'total_gaunt_allowed': 0,
        'all_left': 0,
        'all_right': 0,
        '1L3R': 0,
        '2L2R': 0,
        '3L1R': 0,
        'by_split': defaultdict(int),
        'cross_cut_quartets': [],
    }

    def is_left(p):
        return p < m_closed

    for p in range(M):
        for q in range(M):
            for r in range(M):
                for s in range(M):
                    # Gaunt-allowed: |l_p - l_r| <= L <= l_p + l_r,
                    # parity (l_p + l_r + L) even, etc.
                    # Check if any L gives nonzero ck_pr and ck_qs.
                    l_p, m_p_q = orbitals[p][1], orbitals[p][2]
                    l_r, m_r_q = orbitals[r][1], orbitals[r][2]
                    l_q, m_q_q = orbitals[q][1], orbitals[q][2]
                    l_s, m_s_q = orbitals[s][1], orbitals[s][2]

                    if m_p_q + m_q_q != m_r_q + m_s_q:
                        continue

                    L_min1 = abs(l_p - l_r)
                    L_max1 = l_p + l_r
                    L_min2 = abs(l_q - l_s)
                    L_max2 = l_q + l_s
                    L_overlap = range(max(L_min1, L_min2), min(L_max1, L_max2) + 1)

                    found_L = False
                    for L in L_overlap:
                        if (l_p + l_r + L) % 2 != 0:
                            continue
                        if (l_q + l_s + L) % 2 != 0:
                            continue
                        # Check Gaunt nonzero for some choice
                        c_pr = _ck_coefficient(l_p, m_p_q, l_r, m_r_q, L)
                        c_qs = _ck_coefficient(l_q, m_q_q, l_s, m_s_q, L)
                        if abs(c_pr) > 1e-12 and abs(c_qs) > 1e-12:
                            found_L = True
                            break
                    if not found_L:
                        continue

                    counts['total_gaunt_allowed'] += 1
                    split = sum(is_left(idx) for idx in [p, q, r, s])
                    if split == 4:
                        counts['all_left'] += 1
                    elif split == 0:
                        counts['all_right'] += 1
                    elif split == 1:
                        counts['1L3R'] += 1
                        counts['by_split'][f'1L3R'] += 1
                    elif split == 2:
                        counts['2L2R'] += 1
                        counts['by_split'][f'2L2R'] += 1
                    elif split == 3:
                        counts['3L1R'] += 1
                        counts['by_split'][f'3L1R'] += 1
                    if 0 < split < 4:
                        counts['cross_cut_quartets'].append(
                            (p, q, r, s, split)
                        )
    return counts


def count_independent_left_operators(orbitals, m_closed: int) -> Dict:
    """At post-orbital cut m_closed, count independent left fermionic
    *operators* arising from the Hamiltonian, modulo spin.

    Each cross-cut ERI term g[p,r,q,s] decomposes into a tensor product:
       (left fermionic operator) (x) (right fermionic operator)

    Examples (with L = on left, R = on right):
      - g[L, R, L, R]: a_{p_L}^† a_{q_L} on left * a_{r_R}^† a_{s_R} on right
        (number-like 2L2R, density-density)
      - g[L, R, R, L]: a_{p_L}^† a_{s_L} on left * a_{r_R}^† a_{q_R} on right
        (exchange-like 2L2R, transposed)
      - g[L, R, R, R]: a_{p_L}^† on left * complex 3R operator on right
        (1L3R, hopping)

    Each LEFT fermionic operator pattern contributes 1 bond-rank dim, but
    with spin doubling.

    Returns a count and a list of distinct (left operator pattern, spatial).
    """
    M = len(orbitals)
    distinct_left_ops = set()
    by_type = defaultdict(set)

    def is_left(p):
        return p < m_closed

    for p in range(M):
        for q in range(M):
            for r in range(M):
                for s in range(M):
                    # Gaunt check
                    l_p, m_p_q = orbitals[p][1], orbitals[p][2]
                    l_r, m_r_q = orbitals[r][1], orbitals[r][2]
                    l_q, m_q_q = orbitals[q][1], orbitals[q][2]
                    l_s, m_s_q = orbitals[s][1], orbitals[s][2]
                    if m_p_q + m_q_q != m_r_q + m_s_q:
                        continue
                    L_overlap = range(max(abs(l_p - l_r), abs(l_q - l_s)),
                                      min(l_p + l_r, l_q + l_s) + 1)
                    found = False
                    for L in L_overlap:
                        if (l_p + l_r + L) % 2 != 0 or (l_q + l_s + L) % 2 != 0:
                            continue
                        c1 = _ck_coefficient(l_p, m_p_q, l_r, m_r_q, L)
                        c2 = _ck_coefficient(l_q, m_q_q, l_s, m_s_q, L)
                        if abs(c1) > 1e-12 and abs(c2) > 1e-12:
                            found = True
                            break
                    if not found:
                        continue

                    split = (is_left(p), is_left(q), is_left(r), is_left(s))
                    if all(split) or not any(split):
                        continue

                    # The chemist ERI a_p^† a_q^† a_s a_r in physicist
                    # The fermionic operator is a_p^† a_q^† a_s a_r.
                    # Left part: orbitals that are on left.
                    # For the bond rank, we need the "left" operator as a
                    # specific fermionic-operator pattern.
                    # Split the dagger/no-dagger:
                    # p: dagger, q: dagger, s: no, r: no (chemist convention)
                    left_indices_with_dag = []
                    if is_left(p):
                        left_indices_with_dag.append((p, '†'))
                    if is_left(q):
                        left_indices_with_dag.append((q, '†'))
                    if is_left(s):
                        left_indices_with_dag.append((s, ''))
                    if is_left(r):
                        left_indices_with_dag.append((r, ''))

                    # Canonical form: sort by orbital index, separated by dag
                    left_op_key = tuple(sorted(left_indices_with_dag))
                    n_left_dag = sum(1 for _, d in left_indices_with_dag if d == '†')
                    n_left_undag = sum(1 for _, d in left_indices_with_dag if d == '')

                    distinct_left_ops.add(left_op_key)
                    by_type[(n_left_dag, n_left_undag)].add(left_op_key)

    return {
        'n_distinct_left_ops_spatial': len(distinct_left_ops),
        'by_dag_undag_count': {str(k): len(v) for k, v in by_type.items()},
    }


def main():
    orbitals = enumerate_orbitals_5(n_max=2, l_min=0)
    print(f"5-orbital sub-block: {[str(o) for o in orbitals]}")

    print(f"\n{'='*70}")
    print(f"Cross-cut ERI quartet counts at each post-orbital cut")
    print(f"{'='*70}")

    chi_measured = {0: 1, 1: 16, 2: 9, 3: 9, 4: 3, 5: 2}

    for m_closed in [0, 1, 2, 3, 4, 5]:
        counts = count_cross_cut_eri(orbitals, m_closed)
        n_total = counts['total_gaunt_allowed']
        n_cross = counts['1L3R'] + counts['2L2R'] + counts['3L1R']
        print(f"\nm_closed = {m_closed}:")
        print(f"  Total Gaunt-allowed ERI quartets:  {n_total}")
        print(f"  All-left + all-right:               "
              f"{counts['all_left']} + {counts['all_right']}")
        print(f"  Cross-cut quartets (1L3R+2L2R+3L1R): {n_cross}")
        print(f"    1L3R: {counts['1L3R']}, 2L2R: {counts['2L2R']}, "
              f"3L1R: {counts['3L1R']}")

        op_counts = count_independent_left_operators(orbitals, m_closed)
        print(f"  Distinct LEFT fermionic operators (spatial only): "
              f"{op_counts['n_distinct_left_ops_spatial']}")
        print(f"    by (n_dag, n_undag): {op_counts['by_dag_undag_count']}")
        print(f"  PREDICTED chi_k = 2 + n_distinct_spatial = "
              f"{2 + op_counts['n_distinct_left_ops_spatial']}")
        print(f"  MEASURED chi_k = {chi_measured.get(m_closed, '?')}")

    # Hypothesis check: chi_k - 2 ratio
    print(f"\n{'='*70}")
    print(f"Hypothesis check")
    print(f"{'='*70}")
    print(f"{'m_closed':>10s}  {'chi_meas':>10s}  {'chi-2':>8s}  "
          f"{'n_dist_left_spatial':>22s}  {'spin_factor':>12s}")
    for m_closed in [0, 1, 2, 3, 4, 5]:
        op_counts = count_independent_left_operators(orbitals, m_closed)
        n_d = op_counts['n_distinct_left_ops_spatial']
        chi = chi_measured.get(m_closed, 0)
        chi_minus_2 = chi - 2
        ratio = chi_minus_2 / max(n_d, 1)
        print(f"  {m_closed:>10d}  {chi:>10d}  {chi_minus_2:>8d}  "
              f"{n_d:>22d}  {ratio:>12.4f}")


if __name__ == '__main__':
    main()
