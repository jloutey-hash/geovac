"""
Track DF Sprint 2: Analytical ERI sparsity analysis for nested H-set basis.

Computes the number of nonzero two-electron integrals (ERIs) in the
H-set coupled angular momentum basis for 4 electrons at L=0,
and compares to the uncoupled (flat) basis and composed architecture.

The H-set coupling scheme:
  Individual: l1, l2, l3, l4
  Pair-coupled: l12 in |l1-l2|..l1+l2,  l34 in |l3-l4|..l3+l4
  Total: L in |l12-l34|..l12+l34,  we fix L=0 => l12 = l34

ERI selection rules arise from:
  1. Gaunt integrals: |l - l'| <= k <= l + l', parity constraint
  2. 6j recoupling: additional vanishing from Racah coefficients
  3. Triangle inequality on coupled angular momenta

We count nonzero ERIs for:
  (a) Intra-pair: V_ee(1,2) and V_ee(3,4) — block-diagonal in partner pair
  (b) Inter-pair: V_ee(1,3), V_ee(1,4), V_ee(2,3), V_ee(2,4) — couples both pairs
"""

import numpy as np
from itertools import product
from sympy.physics.wigner import wigner_3j, wigner_6j
from sympy import N as sympify_N
from typing import List, Tuple, Dict
import json


# ============================================================================
# Gaunt integral selection rules
# ============================================================================

def gaunt_nonzero(l1: int, l2: int, k: int) -> bool:
    """Check if Gaunt integral G(l1, k, l2) can be nonzero.

    Selection rules:
    1. Triangle inequality: |l1 - l2| <= k <= l1 + l2
    2. Parity: l1 + k + l2 must be even
    """
    if k < abs(l1 - l2) or k > l1 + l2:
        return False
    if (l1 + k + l2) % 2 != 0:
        return False
    return True


def allowed_k_values(l1: int, l2: int, l_max: int) -> List[int]:
    """Return all allowed multipole orders k for Gaunt integral G(l1, k, l2)."""
    ks = []
    for k in range(0, 2 * l_max + 1):
        if gaunt_nonzero(l1, l2, k):
            ks.append(k)
    return ks


# ============================================================================
# Basis enumeration
# ============================================================================

def enumerate_uncoupled_basis_M0(l_max: int) -> List[Tuple[int, int, int, int]]:
    """Enumerate uncoupled (l1, l2, l3, l4) channels with M=0 (all m_i=0)."""
    basis = []
    for l1, l2, l3, l4 in product(range(l_max + 1), repeat=4):
        basis.append((l1, l2, l3, l4))
    return basis


def enumerate_hset_basis_L0(l_max: int) -> List[Tuple[int, int, int, int, int]]:
    """Enumerate H-set coupled basis states at L=0.

    Each state is (l1, l2, l12, l3, l4) with l34 = l12 (since L=0).

    Returns list of (l1, l2, l12, l3, l4) tuples.
    """
    basis = []
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            for l12 in range(abs(l1 - l2), l1 + l2 + 1):
                for l3 in range(l_max + 1):
                    for l4 in range(l_max + 1):
                        l34 = l12  # L=0 requires l12 = l34
                        if l34 < abs(l3 - l4) or l34 > l3 + l4:
                            continue
                        basis.append((l1, l2, l12, l3, l4))
    return basis


# ============================================================================
# 6j symbol evaluation (cached)
# ============================================================================

_6j_cache = {}

def sixj(j1, j2, j3, j4, j5, j6) -> float:
    """Evaluate Wigner 6j symbol with caching."""
    key = (j1, j2, j3, j4, j5, j6)
    if key not in _6j_cache:
        val = float(sympify_N(wigner_6j(j1, j2, j3, j4, j5, j6)))
        _6j_cache[key] = val
    return _6j_cache[key]


# ============================================================================
# Intra-pair ERI: V_ee(1,2) in coupled basis
# ============================================================================

def count_intrapair_12_eri(basis: List[Tuple], l_max: int) -> Tuple[int, int]:
    """Count nonzero ERI matrix elements for V_ee(1,2) in H-set basis.

    V_ee(1,2) acts on electrons 1 and 2 only.
    In the H-set basis, l3 and l4 are spectators, so l34 is conserved.
    Since L=0 requires l12=l34, l12 is also conserved (l12 = l34 = l12').

    The coupling matrix element is:
    <l1' l2' l12' | V_k(1,2) | l1 l2 l12>
    = sum_k R_k(alpha) * <l1'|P_k|l1> * <l2'|P_k|l2> * (recoupling)

    Selection rules:
    - l3' = l3, l4' = l4 (spectators)
    - l12' = l12 (from l34 conservation + L=0)
    - Gaunt: |l1-l1'| <= k <= l1+l1', |l2-l2'| <= k <= l2+l2'
    - 6j recoupling: {l1' l2' l12; k k' ???} — but since l12 is conserved
      and we sum over m's, the recoupling factor is:
      (-1)^{l1'+l2'+l12} * {l1 l2 l12; l2' l1' k} (Racah W)

    Returns (nonzero, total) pair counts.
    """
    n = len(basis)
    nonzero = 0
    total = 0

    for i, (l1, l2, l12, l3, l4) in enumerate(basis):
        for j, (l1p, l2p, l12p, l3p, l4p) in enumerate(basis):
            total += 1

            # Spectator conservation
            if l3p != l3 or l4p != l4:
                continue
            # l12 conservation (since l34=l12 and l34'=l12', and l3'=l3, l4'=l4)
            if l12p != l12:
                continue

            # Check if any multipole k allows nonzero coupling
            any_k = False
            for k in range(0, 2 * l_max + 1):
                if not gaunt_nonzero(l1, l1p, k):
                    continue
                if not gaunt_nonzero(l2, l2p, k):
                    continue
                # 6j recoupling factor: {l1 l2 l12; l2' l1' k}
                # This can vanish even when individual Gaunt integrals don't
                val = sixj(l1, l2, l12, l2p, l1p, k)
                if abs(val) > 1e-15:
                    any_k = True
                    break

            if any_k:
                nonzero += 1

    return nonzero, total


def count_intrapair_34_eri(basis: List[Tuple], l_max: int) -> Tuple[int, int]:
    """Count nonzero ERI for V_ee(3,4). Symmetric to V_ee(1,2)."""
    n = len(basis)
    nonzero = 0
    total = 0

    for i, (l1, l2, l12, l3, l4) in enumerate(basis):
        for j, (l1p, l2p, l12p, l3p, l4p) in enumerate(basis):
            total += 1

            # Spectator conservation
            if l1p != l1 or l2p != l2:
                continue
            if l12p != l12:
                continue

            # l34 = l12, l34' = l12' = l12
            l34 = l12

            any_k = False
            for k in range(0, 2 * l_max + 1):
                if not gaunt_nonzero(l3, l3p, k):
                    continue
                if not gaunt_nonzero(l4, l4p, k):
                    continue
                val = sixj(l3, l4, l34, l4p, l3p, k)
                if abs(val) > 1e-15:
                    any_k = True
                    break

            if any_k:
                nonzero += 1

    return nonzero, total


# ============================================================================
# Inter-pair ERI: V_ee(1,3) in coupled basis
# ============================================================================

def count_interpair_13_eri(basis: List[Tuple], l_max: int) -> Tuple[int, int]:
    """Count nonzero ERI for V_ee(1,3) in H-set basis.

    V_ee(1,3) couples electrons from different pairs. Both l12 and l34
    can change (since both l1 and l3 are active).

    The coupling involves:
    1. Gaunt on electron 1: |l1-l1'| <= k <= l1+l1'
    2. Gaunt on electron 3: |l3-l3'| <= k <= l3+l3' (same k for given multipole)
    3. Spectator: l2' = l2, l4' = l4
    4. Recoupling: 9j symbol connecting (l1,l2,l12)(l3,l4,l34)(k,0,k)
       to (l1',l2',l12')(l3',l4',l34')(k,0,k)

    The matrix element factorizes as:
    <l1'l2'l12'; l3'l4'l34'; L=0 | V_k(1,3) | l1 l2 l12; l3 l4 l34; L=0>

    With l2'=l2, l4'=l4, this involves a recoupling from
    [(l1,l2)l12, (l3,l4)l34]_L=0  to  [(l1,l3)_k, (l2,l4)_?]
    which is a 9j symbol.

    Selection rules:
    - l2' = l2, l4' = l4 (spectators)
    - l12' in |l1'-l2|..l1'+l2 (new pair coupling)
    - l34' = l12' (L=0 constraint: l34' = l12')
    - l34' in |l3'-l4|..l3'+l4 (must be valid)
    - Gaunt: same k for both electrons 1 and 3
    - 9j recoupling may introduce additional zeros
    """
    n = len(basis)
    nonzero = 0
    total = 0

    for i, (l1, l2, l12, l3, l4) in enumerate(basis):
        l34 = l12  # L=0
        for j, (l1p, l2p, l12p, l3p, l4p) in enumerate(basis):
            l34p = l12p  # L=0
            total += 1

            # Spectator conservation
            if l2p != l2 or l4p != l4:
                continue

            # Check if any multipole k allows nonzero coupling
            any_k = False
            for k in range(0, 2 * l_max + 1):
                if not gaunt_nonzero(l1, l1p, k):
                    continue
                if not gaunt_nonzero(l3, l3p, k):
                    continue

                # 9j recoupling:
                # { l1  l2  l12  }
                # { l3  l4  l34  }
                # { k   0   k    }  (since L=0, the bottom-right is k when coupled)
                #
                # But actually, for V_ee(1,3) with multipole P_k,
                # the recoupling connects:
                # [(l1,l2)l12, (l3,l4)l34]_0 <-> [(l1,l3)_X, (l2,l4)_Y]_0
                #
                # The relevant coefficient is the 9j symbol:
                # { l1   l2   l12  }
                # { l3   l4   l34  }
                # { l1'  l2'  l12' }  -- but this isn't quite right either
                #
                # More precisely, the matrix element of T^k(1) * T^k(3)
                # between coupled states involves:
                # sum over intermediate couplings, with the
                # Wigner-Eckart theorem applied twice.
                #
                # For the specific case of a scalar product T^k(1).T^k(3),
                # the result in the coupled scheme is:
                #
                # <(l1'l2)l12'; (l3'l4)l34'; 0| T^k(1).T^k(3) |(l1 l2)l12; (l3 l4)l34; 0>
                #
                # = (-1)^{l12'+l34+k} * {l12' l34' 0; l34 l12 k} (this is 0 unless l12'=l34')
                #   * (-1)^{l2+l1'+l12} * {l1 l12 l2; l12' l1' k}  [pair 12 recoupling]
                #   * (-1)^{l4+l3'+l34} * {l3 l34 l4; l34' l3' k}  [pair 34 recoupling]
                #   * <l1'||T^k||l1> * <l3'||T^k||l3>  [reduced matrix elements]
                #
                # The key 6j symbol {l12' l34' 0; l34 l12 k} with L=0 gives:
                # = delta(l12', l34') * (-1)^{l12+k} / sqrt(2*l12+1) * delta(l12', l34')
                #   ... actually for {a b 0; b' a' c}:
                # {a b 0; b' a' c} = delta(a,a') * delta(b,b') * (-1)^{a+b+c} / sqrt((2a+1)(2b+1))
                #   NO — that's only when the third column is 0.
                #
                # Let me use the actual 6j identity:
                # {j1 j2 0; j3 j4 j5} = delta(j1,j4)*delta(j2,j3) * (-1)^{j1+j2+j5} / sqrt((2j1+1)(2j2+1))
                #   when j3=j2, j4=j1.
                # But here the third column is (0, k, k) not (0, 0, 0).
                # So we have {l12' l34' 0; l34 l12 k}.
                # The top-right is 0, so: delta(l12',l12)*delta(l34',l34) ... NO.
                # Actually {a b J; d e k} with J=0: requires a=b and d=e,
                # and then = (-1)^{a+d+k} * delta(a,b)*delta(d,e) / sqrt((2a+1)(2d+1))
                #
                # Wait: {a b 0; d e k}: triangle rules require:
                #   (a, b, 0) => a = b
                #   (d, e, k) => |d-e| <= k <= d+e
                #   (a, e, k) => ... hmm, wrong. The 6j has rows/columns:
                #   {j1 j2 j3; j4 j5 j6}:
                #     triads: (j1,j2,j3), (j1,j5,j6), (j4,j2,j6), (j4,j5,j3)
                #   So {l12' l34' 0; l34 l12 k}:
                #     triads: (l12', l34', 0) => l12' = l34' [already enforced by L=0]
                #             (l12', l12, k) => |l12'-l12| <= k <= l12'+l12
                #             (l34, l34', k) => |l34-l34'| <= k <= l34+l34'
                #             (l34, l12, 0) => l34 = l12 [already enforced]
                #
                #   Since l12'=l34' and l34=l12, the conditions become:
                #     |l12'-l12| <= k <= l12'+l12
                #   This is a SELECTION RULE on how much l12 can change!
                #
                # So the full selection rules for V_ee(1,3) are:
                # 1. l2' = l2, l4' = l4 (spectators)
                # 2. Gaunt on l1: |l1-l1'| <= k <= l1+l1', parity
                # 3. Gaunt on l3: |l3-l3'| <= k <= l3+l3', parity
                # 4. |l12'-l12| <= k <= l12'+l12 (from top-level 6j)
                # 5. l12' in |l1'-l2|..l1'+l2 and l12' in |l3'-l4|..l3'+l4
                # 6. Two 6j factors for pair recoupling must be nonzero

                # Check condition 4: |l12'-l12| <= k
                if k < abs(l12p - l12) or k > l12p + l12:
                    continue

                # 6j for pair 12 recoupling: {l1 l2 l12; l12' l1' k}
                # Triad check: (l1,l2,l12), (l1,l1',k) -- wait, 6j triads:
                # {a b c; d e f}: (a,b,c), (a,e,f), (d,b,f), (d,e,c)
                # {l1 l2 l12; l12' l1' k}:
                #   (l1,l2,l12) ok by construction
                #   (l1,l1',k) ok by Gaunt check above
                #   (l12',l2,k) need: |l12'-l2| <= k <= l12'+l2
                #   (l12',l1',l12) need: |l12'-l1'| <= l12 <= l12'+l1'

                # Check 6j triad conditions
                if k < abs(l12p - l2) or k > l12p + l2:
                    continue
                if l12 < abs(l12p - l1p) or l12 > l12p + l1p:
                    continue

                val_12 = sixj(l1, l2, l12, l12p, l1p, k)
                if abs(val_12) < 1e-15:
                    continue

                # 6j for pair 34 recoupling: {l3 l4 l34; l34' l3' k}
                # = {l3 l4 l12; l12' l3' k} since l34=l12, l34'=l12'
                if k < abs(l12p - l4) or k > l12p + l4:
                    continue
                if l12 < abs(l12p - l3p) or l12 > l12p + l3p:
                    continue

                val_34 = sixj(l3, l4, l12, l12p, l3p, k)
                if abs(val_34) < 1e-15:
                    continue

                # Top-level 6j: {l12' l12' 0; l12 l12 k}
                # = (-1)^{2*l12+k} / (2*l12+1) ... only if triangle holds
                val_top = sixj(l12p, l12p, 0, l12, l12, k)
                if abs(val_top) < 1e-15:
                    continue

                any_k = True
                break

            if any_k:
                nonzero += 1

    return nonzero, total


# ============================================================================
# Uncoupled (flat) basis ERI counting
# ============================================================================

def count_uncoupled_eri_pair(basis_unc: List[Tuple], pair: Tuple[int, int],
                             l_max: int) -> Tuple[int, int]:
    """Count nonzero ERI for V_ee(i,j) in uncoupled basis.

    pair = (i, j) with i,j in {0,1,2,3} indexing electrons.
    Spectators: electrons not in pair must have same l.
    Active: Gaunt selection rule on each active electron.
    """
    i, j = pair
    spectators = [x for x in range(4) if x not in (i, j)]

    n = len(basis_unc)
    nonzero = 0
    total = 0

    for a, state in enumerate(basis_unc):
        for b, statep in enumerate(basis_unc):
            total += 1

            # Spectator check
            if any(state[s] != statep[s] for s in spectators):
                continue

            li, lj = state[i], state[j]
            lip, ljp = statep[i], statep[j]

            any_k = False
            for k in range(0, 2 * l_max + 1):
                if gaunt_nonzero(li, lip, k) and gaunt_nonzero(lj, ljp, k):
                    any_k = True
                    break

            if any_k:
                nonzero += 1

    return nonzero, total


# ============================================================================
# S_4 [2,2] projection
# ============================================================================

def s4_22_projected_dimension(l_max: int) -> int:
    """Compute dimension of S4 [2,2] subspace in uncoupled basis.

    Uses character formula: dim = (2/24) * sum_g chi(g) * Tr(D(g))
    """
    from itertools import permutations

    basis = enumerate_uncoupled_basis_M0(l_max)
    n = len(basis)
    ch_index = {ch: i for i, ch in enumerate(basis)}

    # S4 characters for [2,2]
    def cycle_type(perm):
        visited = [False] * 4
        cycles = []
        for i in range(4):
            if visited[i]:
                continue
            clen = 0
            j = i
            while not visited[j]:
                visited[j] = True
                j = perm[j]
                clen += 1
            cycles.append(clen)
        return tuple(sorted(cycles, reverse=True))

    chi_map = {
        (1,1,1,1): 2, (2,1,1): 0, (3,1): -1, (4,): 0, (2,2): 2
    }

    dim_proj = 0.0
    for perm in permutations(range(4)):
        chi = chi_map[cycle_type(perm)]
        if chi == 0:
            continue
        # Trace of D(perm)
        tr = 0
        for ch in basis:
            permuted = tuple(ch[perm[k]] for k in range(4))
            if permuted == ch:
                tr += 1
        dim_proj += chi * tr

    dim_proj *= 2 / 24  # dim_irrep / |G|
    return int(round(dim_proj))


# ============================================================================
# Main analysis
# ============================================================================

def run_analysis(l_max: int) -> Dict:
    """Run full ERI sparsity analysis at given l_max."""
    print(f"\n{'='*70}")
    print(f"  ERI SPARSITY ANALYSIS: l_max = {l_max}")
    print(f"{'='*70}")

    # Enumerate bases
    basis_unc = enumerate_uncoupled_basis_M0(l_max)
    basis_hset = enumerate_hset_basis_L0(l_max)
    n_unc = len(basis_unc)
    n_hset = len(basis_hset)

    print(f"\nUncoupled basis (M=0): {n_unc} states")
    print(f"H-set coupled basis (L=0): {n_hset} states")

    # S4 projected dimensions
    dim_22 = s4_22_projected_dimension(l_max)
    print(f"S4 [2,2] projected dimension (uncoupled): {dim_22}")

    # --- Uncoupled ERI counting ---
    print(f"\n--- Uncoupled basis ERI counts ---")

    # Intra-pair (0,1)
    nz_unc_01, tot_unc_01 = count_uncoupled_eri_pair(basis_unc, (0, 1), l_max)
    dens_unc_01 = nz_unc_01 / tot_unc_01 * 100
    print(f"V_ee(1,2) uncoupled: {nz_unc_01}/{tot_unc_01} = {dens_unc_01:.1f}%")

    # Inter-pair (0,2)
    nz_unc_02, tot_unc_02 = count_uncoupled_eri_pair(basis_unc, (0, 2), l_max)
    dens_unc_02 = nz_unc_02 / tot_unc_02 * 100
    print(f"V_ee(1,3) uncoupled: {nz_unc_02}/{tot_unc_02} = {dens_unc_02:.1f}%")

    # Total ERI (all 6 pairs)
    total_eri_unc = n_unc * n_unc  # total possible per pair
    total_nz_unc = 0
    for pair in [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]:
        nz, _ = count_uncoupled_eri_pair(basis_unc, pair, l_max)
        total_nz_unc += nz
    dens_unc_total = total_nz_unc / (6 * total_eri_unc) * 100
    print(f"Total V_ee (all 6 pairs) uncoupled: {total_nz_unc}/{6*total_eri_unc} = {dens_unc_total:.1f}%")

    # --- H-set coupled ERI counting ---
    print(f"\n--- H-set coupled basis ERI counts ---")

    # Intra-pair V_ee(1,2)
    nz_h12, tot_h12 = count_intrapair_12_eri(basis_hset, l_max)
    dens_h12 = nz_h12 / tot_h12 * 100
    print(f"V_ee(1,2) H-set: {nz_h12}/{tot_h12} = {dens_h12:.1f}%")

    # Intra-pair V_ee(3,4)
    nz_h34, tot_h34 = count_intrapair_34_eri(basis_hset, l_max)
    dens_h34 = nz_h34 / tot_h34 * 100
    print(f"V_ee(3,4) H-set: {nz_h34}/{tot_h34} = {dens_h34:.1f}%")

    # Inter-pair V_ee(1,3)
    nz_h13, tot_h13 = count_interpair_13_eri(basis_hset, l_max)
    dens_h13 = nz_h13 / tot_h13 * 100
    print(f"V_ee(1,3) H-set: {nz_h13}/{tot_h13} = {dens_h13:.1f}%")

    # By symmetry, V_ee(2,4) has same density as V_ee(1,3)
    # V_ee(1,4) and V_ee(2,3) also have same density (cross-pair, one from each)
    # For the total, we need all 6 pairs. Intra: (1,2)+(3,4). Inter: (1,3)+(1,4)+(2,3)+(2,4).
    # By pair symmetry: (1,3)=(2,4) and (1,4)=(2,3)
    # For now approximate inter-pair density as same for all 4 cross terms.
    total_nz_hset = nz_h12 + nz_h34 + 4 * nz_h13
    total_possible_hset = tot_h12 + tot_h34 + 4 * tot_h13
    dens_hset_total = total_nz_hset / total_possible_hset * 100
    print(f"Total V_ee (est, 2 intra + 4 inter) H-set: {total_nz_hset}/{total_possible_hset} = {dens_hset_total:.1f}%")

    # --- Summary ---
    print(f"\n--- Summary ---")
    print(f"{'Metric':<40} {'Uncoupled':>12} {'H-set':>12}")
    print(f"{'Basis states':<40} {n_unc:>12} {n_hset:>12}")
    print(f"{'Intra-pair V_ee(1,2) density':<40} {dens_unc_01:>11.1f}% {dens_h12:>11.1f}%")
    print(f"{'Inter-pair V_ee(1,3) density':<40} {dens_unc_02:>11.1f}% {dens_h13:>11.1f}%")
    print(f"{'Total ERI density (all 6 pairs)':<40} {dens_unc_total:>11.1f}% {dens_hset_total:>11.1f}%")

    results = {
        'l_max': l_max,
        'uncoupled': {
            'n_states': n_unc,
            's4_22_dim': dim_22,
            'vee_12': {'nonzero': nz_unc_01, 'total': tot_unc_01, 'density_pct': dens_unc_01},
            'vee_13': {'nonzero': nz_unc_02, 'total': tot_unc_02, 'density_pct': dens_unc_02},
            'total_density_pct': dens_unc_total,
        },
        'hset': {
            'n_states': n_hset,
            'vee_12': {'nonzero': nz_h12, 'total': tot_h12, 'density_pct': dens_h12},
            'vee_34': {'nonzero': nz_h34, 'total': tot_h34, 'density_pct': dens_h34},
            'vee_13': {'nonzero': nz_h13, 'total': tot_h13, 'density_pct': dens_h13},
            'total_density_pct': dens_hset_total,
        },
    }

    return results


if __name__ == '__main__':
    all_results = {}

    # l_max = 1
    r1 = run_analysis(l_max=1)
    all_results['l_max_1'] = r1

    # l_max = 2 (may be slow due to 6j evaluations)
    print("\n\nStarting l_max=2 analysis (may take a few minutes)...")
    r2 = run_analysis(l_max=2)
    all_results['l_max_2'] = r2

    # Save results
    output_path = 'debug/data/track_df_eri_sparsity.json'
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {output_path}")

    # Go/no-go assessment
    print(f"\n{'='*70}")
    print("  GO/NO-GO ASSESSMENT")
    print(f"{'='*70}")
    gate = 15.0  # percent
    for key, res in all_results.items():
        dens = res['hset']['total_density_pct']
        status = "PASS" if dens <= gate else "FAIL"
        print(f"{key}: H-set total ERI density = {dens:.1f}% — {status} (gate: {gate}%)")
