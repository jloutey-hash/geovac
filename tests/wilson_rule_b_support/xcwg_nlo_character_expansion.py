"""
XCWG-E: Next-to-leading-order strong-coupling character expansion on Rule B.

Date: 2026-05-16
Follows XCWG-D (leading-order area law, sigma>0 at all beta) and tests
whether the NLO-corrected string tension sigma_NLO(beta) stays positive at
all beta (3D compact U(1) signature) or crosses zero at finite beta_c
(4D-style deconfinement transition).

Theoretical setup
=================

For compact U(1) lattice gauge theory with Wilson action
    S_W = -beta sum_P cos(theta_P)
the character expansion gives, for a Wilson loop W with edge indicator w in Z^E,
    <W>_beta = N(w; beta) / Z(beta)
where
    N(w; beta) = sum_{n in Z^P : d_1^T n = w} prod_P (I_{n_p}(beta) / I_0(beta))
    Z(beta)    = sum_{n in Z^P : d_1^T n = 0} prod_P (I_{n_p}(beta) / I_0(beta))

LO:
  - N includes the trivial surface n = n_LO (e.g. n_p=1 on A_min plaquettes), giving (I_1/I_0)^A_min
  - Z includes only n=0, giving 1
  - <W>_LO = (I_1/I_0)^A_min

NLO:
  - N additionally includes n_LO + z for each closed 2-cycle z (size k_delta and up)
  - Z additionally includes pure closed 2-cycles z alone
  - <W>_NLO = N_NLO / Z_NLO

We define
  sigma_NLO(beta) := -log <W>_NLO / A_min

If sigma_NLO(beta) > 0 at all beta: confining at NLO (3D compact U(1)).
If sigma_NLO(beta) crosses zero at some beta_c: deconfinement transition (4D-like).

Methodology
===========

1. Build d_1: P x E plaquette-boundary operator (each row is a signed L=4 plaquette)
2. Enumerate small closed 2-cycles (integer solutions of d_1^T z = 0)
3. For a representative Wilson loop W with A_min plaquettes:
   - Compute LO surface assignment n_LO
   - Generate N_NLO by adding closed 2-cycles to n_LO (with bounded total |n|_1)
   - Generate Z_NLO independently by enumerating closed 2-cycles
   - Compute <W>_NLO = N_NLO / Z_NLO
   - Define sigma_NLO = -log <W>_NLO / A_min

Output
======
    - tests/wilson_rule_b_support/data/xcwg_nlo_character_expansion.json
    - debug/xcwg_nlo_character_expansion_memo.md (separately)

Notes
=====
- Uses scipy.special.iv for Bessel functions; supports I_0 and I_1 (LO) plus
  I_{>=2} for higher-charge configurations (subleading in beta)
- Surface enumeration over Z is exponential; we cap at |n|_1 <= A_min + k_delta + 2
- At n_max=3 we restrict to a single Wilson loop and a small cycle search
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import numpy as np
from scipy.special import iv

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))  # repo root (tests/wilson_rule_b_support/ -> two levels up)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from geovac.ihara_zeta_dirac import build_dirac_s3_graph  # noqa: E402

from xcwg_wilson_loop_scaling import (  # noqa: E402
    signed_incidence, adjacency_list, enumerate_primitive_closed_walks,
)
from xcwg_strong_coupling_wilson import (  # noqa: E402
    edge_set_of_walk, adjacent_plaquettes,
    build_area_A_loop, is_valid_closed_loop, i1_over_i0, sigma_beta,
)

sys.stdout.reconfigure(line_buffering=True)


# =============================================================================
# Bessel ratio I_n(beta) / I_0(beta) for any n
# =============================================================================

def In_over_I0(n: int, beta: float) -> float:
    """Compute I_n(beta) / I_0(beta) for integer n (positive or negative)."""
    n_abs = abs(n)
    if n_abs == 0:
        return 1.0
    if beta <= 0.0:
        return 0.0
    if beta < 100.0:
        num = float(iv(n_abs, beta))
        den = float(iv(0, beta))
        if den <= 0.0:
            return 0.0
        return num / den
    # Large-beta asymptotic: I_n/I_0 -> 1 - (4n^2-1)/(8 beta) - ...
    correction = -(4 * n_abs * n_abs - 1) / (8.0 * beta)
    return 1.0 + correction


# =============================================================================
# Build plaquettes and d_1
# =============================================================================

def build_d1_and_plaquettes(n_max: int, plaq_cap: int = 500) -> Dict:
    """Build d_1 and plaquette list."""
    t0 = time.time()
    A, labels, deg, desc = build_dirac_s3_graph(n_max, "B")
    V = A.shape[0]
    E_tot = int(A.sum() // 2)
    B_inc, edges, edge_idx = signed_incidence(A)
    adj = adjacency_list(A)

    plaquettes: List[Tuple[int, ...]] = []
    for walk in enumerate_primitive_closed_walks(adj, 4):
        plaquettes.append(walk)
        if len(plaquettes) >= plaq_cap:
            break
    P = len(plaquettes)
    capped = len(plaquettes) >= plaq_cap

    d_1 = np.zeros((P, E_tot), dtype=np.int8)
    plaq_edge_sets: List[frozenset] = []
    plaq_signed_vecs: List[np.ndarray] = []
    for p, walk in enumerate(plaquettes):
        L = len(walk)
        for i in range(L):
            u = walk[i]
            v = walk[(i + 1) % L]
            ei = edge_idx[(u, v)]
            e_canon = (min(u, v), max(u, v))
            sign = +1 if (u, v) == e_canon else -1
            d_1[p, ei] += sign
        plaq_edge_sets.append(edge_set_of_walk(walk, edge_idx))
        plaq_signed_vecs.append(d_1[p].astype(np.float64).copy())

    plaq_adj = adjacent_plaquettes(plaq_edge_sets)

    return {
        "n_max": n_max,
        "V": V, "E": E_tot, "P": P,
        "beta_1_graph": E_tot - V + 1,
        "d_1": d_1, "plaquettes": plaquettes,
        "plaq_edge_sets": plaq_edge_sets,
        "plaq_signed_vecs": plaq_signed_vecs,
        "plaq_adj": plaq_adj,
        "edges": edges, "edge_idx": edge_idx, "B_inc": B_inc, "adj": adj,
        "labels": labels,
        "plaq_cap_applied": capped,
        "build_time_sec": time.time() - t0,
    }


# =============================================================================
# Closed 2-cycle homology and small-cycle enumeration
# =============================================================================

def compute_2cycle_dimension(d_1: np.ndarray) -> Dict:
    """Q-dimension of ker(d_1^T): closed 2-cycles in Q^P."""
    d_1_float = d_1.astype(np.float64)
    P = d_1.shape[0]
    rank = int(np.linalg.matrix_rank(d_1_float))
    b_2 = P - rank
    return {"rank_d1": rank, "P": P, "b_2": b_2}


def find_smallest_closed_2cycle(
    d_1: np.ndarray,
    plaq_adj: Dict[int, List[int]],
    max_size: int = 8,
    time_budget_sec: float = 60.0,
) -> Optional[Dict]:
    """Find ONE smallest integer closed 2-cycle z (with z_p in {-1, 0, +1}).

    BFS over connected plaquette subsets, increasing size; at each size try
    sign assignments. Returns first cycle found with smallest size.
    """
    P = d_1.shape[0]
    t0 = time.time()

    def signs_for_cycle(subset: List[int]) -> Optional[List[int]]:
        k = len(subset)
        d_sub = d_1[subset, :].astype(np.int64)
        # 2^{k-1} patterns
        for sign_bits in range(2 ** (k - 1)):
            signs = [1] + [1 if (sign_bits >> i) & 1 == 0 else -1
                          for i in range(k - 1)]
            sa = np.asarray(signs, dtype=np.int64)
            bd = sa @ d_sub
            if np.all(bd == 0):
                return signs
        return None

    for size in range(2, max_size + 1):
        if time.time() - t0 > time_budget_sec:
            return None
        # BFS to enumerate connected subsets of size = current
        seen: set = set()
        for start in range(P):
            if time.time() - t0 > time_budget_sec:
                break
            stack: List[Tuple[int, ...]] = [(start,)]
            local_seen = {(start,)}
            while stack:
                if time.time() - t0 > time_budget_sec:
                    break
                cluster = stack.pop()
                k = len(cluster)
                if k == size:
                    sig = tuple(sorted(cluster))
                    if sig in seen:
                        continue
                    seen.add(sig)
                    signs = signs_for_cycle(list(sig))
                    if signs is not None:
                        return {
                            "plaq_set": list(sig),
                            "signs": signs,
                            "size": size,
                        }
                    continue
                # Extend
                candidates: set = set()
                for c in cluster:
                    candidates.update(plaq_adj.get(c, []))
                candidates = candidates - set(cluster)
                for nxt in candidates:
                    new_cluster = tuple(sorted(cluster + (nxt,)))
                    if new_cluster in local_seen:
                        continue
                    local_seen.add(new_cluster)
                    stack.append(cluster + (nxt,))
    return None


def enumerate_closed_2cycles_up_to_size(
    d_1: np.ndarray,
    plaq_adj: Dict[int, List[int]],
    size_target: int,
    max_count: int = 200,
    time_budget_sec: float = 120.0,
) -> List[Dict]:
    """Enumerate closed 2-cycles of size exactly size_target (or smallest available)."""
    P = d_1.shape[0]
    t0 = time.time()
    cycles: List[Dict] = []

    def signs_for_cycle(subset: List[int]) -> Optional[List[int]]:
        k = len(subset)
        d_sub = d_1[subset, :].astype(np.int64)
        for sign_bits in range(2 ** (k - 1)):
            signs = [1] + [1 if (sign_bits >> i) & 1 == 0 else -1
                          for i in range(k - 1)]
            sa = np.asarray(signs, dtype=np.int64)
            bd = sa @ d_sub
            if np.all(bd == 0):
                return signs
        return None

    seen: set = set()
    for start in range(P):
        if len(cycles) >= max_count or time.time() - t0 > time_budget_sec:
            break
        stack: List[Tuple[int, ...]] = [(start,)]
        local_seen = {(start,)}
        while stack and len(cycles) < max_count:
            if time.time() - t0 > time_budget_sec:
                break
            cluster = stack.pop()
            k = len(cluster)
            if k == size_target:
                sig = tuple(sorted(cluster))
                if sig in seen:
                    continue
                seen.add(sig)
                signs = signs_for_cycle(list(sig))
                if signs is not None:
                    cycles.append({
                        "plaq_set": list(sig),
                        "signs": signs,
                        "size": k,
                    })
                continue
            if k > size_target:
                continue
            candidates: set = set()
            for c in cluster:
                candidates.update(plaq_adj.get(c, []))
            candidates = candidates - set(cluster)
            for nxt in candidates:
                new_cluster = tuple(sorted(cluster + (nxt,)))
                if new_cluster in local_seen:
                    continue
                local_seen.add(new_cluster)
                stack.append(cluster + (nxt,))
    return cycles


# =============================================================================
# Normalized NLO Wilson loop: <W>_NLO = N_NLO(W) / Z_NLO
# =============================================================================

def compute_W_normalized(
    d_1: np.ndarray,
    n_LO: np.ndarray,  # length-P int array, LO surface
    w: np.ndarray,     # length-E Wilson loop indicator
    cycles: List[Dict],
    beta: float,
    truncation_n1: int,
    n_max_per_plaq: int = 2,
) -> Dict:
    """Compute <W>_NLO = N(w)/Z properly normalized.

    N(w; beta) = sum_{n: d_1^T n = w, |n|_1 <= truncation} prod_p (I_{n_p}/I_0)
    Z(beta)    = sum_{n: d_1^T n = 0, |n|_1 <= truncation} prod_p (I_{n_p}/I_0)

    Restricts to assignments built as n_LO + integer combinations of cycles for N
    and integer combinations of cycles for Z.
    Single-cycle additions with alpha in {-1, +1, +2, -2} are considered.
    """
    P = d_1.shape[0]

    def contribution(n_arr: np.ndarray) -> float:
        prod = 1.0
        for p in range(P):
            np_charge = int(n_arr[p])
            if np_charge != 0:
                prod *= In_over_I0(np_charge, beta)
        return prod

    def is_valid(n_arr: np.ndarray, target_w: np.ndarray) -> bool:
        check = (n_arr.astype(np.int64) @ d_1).astype(np.int64)
        return np.array_equal(check, target_w.astype(np.int64))

    # ---- Numerator: N(w; beta) ----
    N_total = 0.0
    n_LO_n1 = int(np.sum(np.abs(n_LO)))
    # LO term
    if is_valid(n_LO, w):
        N_total += contribution(n_LO)
    N_LO_term = N_total
    seen_n: set = set()
    seen_n.add(tuple(int(x) for x in n_LO))

    for z_data in cycles:
        z = np.zeros(P, dtype=np.int64)
        for p, s in zip(z_data["plaq_set"], z_data["signs"]):
            z[p] = s
        for alpha in [1, -1, 2, -2]:
            n_new = n_LO + alpha * z
            n1 = int(np.sum(np.abs(n_new)))
            if n1 > truncation_n1:
                continue
            if int(np.abs(n_new).max()) > n_max_per_plaq:
                continue
            sig = tuple(int(x) for x in n_new)
            if sig in seen_n:
                continue
            seen_n.add(sig)
            if is_valid(n_new, w):
                N_total += contribution(n_new)

    # ---- Denominator: Z(beta) ----
    Z_total = 1.0  # trivial n=0 assignment
    seen_z: set = set()
    seen_z.add(tuple([0] * P))
    for z_data in cycles:
        z = np.zeros(P, dtype=np.int64)
        for p, s in zip(z_data["plaq_set"], z_data["signs"]):
            z[p] = s
        for alpha in [1, -1, 2, -2]:
            n_new = alpha * z
            n1 = int(np.sum(np.abs(n_new)))
            if n1 > truncation_n1:
                continue
            if int(np.abs(n_new).max()) > n_max_per_plaq:
                continue
            sig = tuple(int(x) for x in n_new)
            if sig in seen_z:
                continue
            seen_z.add(sig)
            if is_valid(n_new, np.zeros_like(w)):
                Z_total += contribution(n_new)

    W_LO_value = N_LO_term  # = (I_1/I_0)^|n_LO|_1
    W_NLO_value = N_total / Z_total if Z_total > 0 else float("nan")
    return {
        "beta": beta,
        "N_total": N_total,
        "Z_total": Z_total,
        "W_LO": W_LO_value,
        "W_NLO": W_NLO_value,
        "W_NLO_minus_LO": W_NLO_value - W_LO_value,
        "n_N_terms": len(seen_n),
        "n_Z_terms": len(seen_z),
    }


# =============================================================================
# Driver
# =============================================================================

def run_nlo(
    n_max: int,
    plaq_cap: int = 200,
    cycle_size_target: Optional[int] = None,
    cycle_max_count: int = 50,
    truncation_extra: int = 8,
    beta_grid: List[float] = None,
    cycle_time_budget_sec: float = 60.0,
) -> Dict:
    """Run normalized NLO Wilson loop analysis at one n_max."""
    if beta_grid is None:
        beta_grid = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0]

    print(f"\n{'='*72}")
    print(f"Normalized NLO at n_max = {n_max}")
    print(f"{'='*72}")
    g = build_d1_and_plaquettes(n_max, plaq_cap=plaq_cap)
    print(f"  V={g['V']}, E={g['E']}, P={g['P']}, beta_1={g['beta_1_graph']}")

    cyc_info = compute_2cycle_dimension(g["d_1"])
    print(f"  rank(d_1)={cyc_info['rank_d1']}, b_2={cyc_info['b_2']}")

    # Find smallest closed 2-cycle
    print(f"\n  Searching for smallest closed 2-cycle (max size 8, budget {cycle_time_budget_sec}s) ...")
    t1 = time.time()
    smallest = find_smallest_closed_2cycle(g["d_1"], g["plaq_adj"],
                                           max_size=8,
                                           time_budget_sec=cycle_time_budget_sec)
    print(f"    Smallest cycle search took {time.time()-t1:.1f}s")
    if smallest is None:
        print("    ⚠ No closed 2-cycle found within budget")
        k_delta = None
        cycles: List[Dict] = []
    else:
        k_delta = smallest["size"]
        print(f"    k_delta = {k_delta} (sample: plaq_set={smallest['plaq_set']})")
        # Enumerate cycles at this size
        if cycle_size_target is None:
            cycle_size_target = k_delta
        t2 = time.time()
        cycles = enumerate_closed_2cycles_up_to_size(
            g["d_1"], g["plaq_adj"],
            size_target=cycle_size_target,
            max_count=cycle_max_count,
            time_budget_sec=cycle_time_budget_sec,
        )
        print(f"    Enumerated {len(cycles)} cycles of size {cycle_size_target} in {time.time()-t2:.1f}s")

    # Select Wilson loop = plaquette 0 with A_min=1
    w = g["d_1"][0].astype(np.int64)
    n_LO = np.zeros(g["P"], dtype=np.int64)
    n_LO[0] = 1
    A_min = 1
    truncation_n1 = A_min + (k_delta if k_delta is not None else truncation_extra) + 1

    print(f"\n  Wilson loop: plaquette 0, A_min = {A_min}")
    print(f"  Truncation |n|_1 <= {truncation_n1}")

    # Beta scan
    results = []
    print(f"\n  beta scan ({len(beta_grid)} points):")
    print(f"  {'beta':>10}  {'sigma_LO':>14}  {'sigma_NLO':>14}  {'N':>14}  {'Z':>14}  {'<W>_NLO':>14}")
    for beta in beta_grid:
        r = compute_W_normalized(g["d_1"], n_LO, w, cycles, beta, truncation_n1)
        s_NLO = (-np.log(r["W_NLO"]) / A_min
                if r["W_NLO"] > 0 else float("nan"))
        r["sigma_LO"] = float(sigma_beta(beta))
        r["sigma_NLO"] = float(s_NLO)
        results.append(r)
        print(f"  {beta:>10.4g}  {r['sigma_LO']:>14.6e}  {s_NLO:>14.6e}  "
              f"{r['N_total']:>14.6e}  {r['Z_total']:>14.6e}  "
              f"{r['W_NLO']:>14.6e}")

    # Verdict
    all_pos = all(r["sigma_NLO"] > 0 for r in results if not np.isnan(r["sigma_NLO"]))
    sigmas = [r["sigma_NLO"] for r in results if not np.isnan(r["sigma_NLO"])]
    min_sigma = min(sigmas) if sigmas else float("nan")
    min_idx = int(np.argmin([r["sigma_NLO"] for r in results]))
    min_beta = results[min_idx]["beta"]
    # Zero crossing
    signs = [np.sign(r["sigma_NLO"]) if not np.isnan(r["sigma_NLO"]) else 0
             for r in results]
    has_crossing = any(signs[i] * signs[i+1] < 0
                       for i in range(len(signs) - 1)
                       if signs[i] != 0 and signs[i+1] != 0)
    print(f"\n  All sigma_NLO > 0: {all_pos}")
    print(f"  Minimum sigma_NLO: {min_sigma:.6e} at beta = {min_beta:.4g}")
    print(f"  Zero crossing detected: {has_crossing}")

    return {
        "n_max": n_max,
        "V": g["V"], "E": g["E"], "P": g["P"],
        "beta_1_graph": g["beta_1_graph"],
        "plaq_cap_applied": g["plaq_cap_applied"],
        "rank_d1": cyc_info["rank_d1"],
        "b_2": cyc_info["b_2"],
        "k_delta": k_delta,
        "smallest_cycle_sample": smallest,
        "n_cycles_enumerated": len(cycles),
        "cycles_sample": cycles[:5],  # Sample
        "A_min": A_min,
        "truncation_n1": truncation_n1,
        "beta_grid": beta_grid,
        "scan_results": [
            {
                "beta": float(r["beta"]),
                "sigma_LO": float(r["sigma_LO"]),
                "sigma_NLO": float(r["sigma_NLO"]) if not np.isnan(r["sigma_NLO"]) else None,
                "N_total": float(r["N_total"]),
                "Z_total": float(r["Z_total"]),
                "W_LO": float(r["W_LO"]),
                "W_NLO": float(r["W_NLO"]) if not np.isnan(r["W_NLO"]) else None,
                "delta_W": float(r["W_NLO_minus_LO"]) if not np.isnan(r["W_NLO_minus_LO"]) else None,
            }
            for r in results
        ],
        "verdict": {
            "all_positive": all_pos,
            "min_sigma_NLO": float(min_sigma) if not np.isnan(min_sigma) else None,
            "min_beta": float(min_beta),
            "has_zero_crossing": has_crossing,
        },
    }


# =============================================================================
# Creutz-ratio cross-check
# =============================================================================

def creutz_ratios(
    g: Dict,
    cycles: List[Dict],
    A_values: List[int],
    beta_grid: List[float],
    truncation_extra: int = 7,
) -> Dict:
    """For each A in A_values, build a Wilson loop with minimal area A,
    compute <W>_NLO, and report sigma_eff(A1, A2) := -log[W(A2)/W(A1)] / (A2-A1).

    At LO, sigma_eff is constant = sigma_LO(beta) for all A pairs.
    At NLO, sigma_eff varies with A; the asymptote at large A is the true sigma.
    """
    d_1 = g["d_1"]
    P = d_1.shape[0]
    plaq_adj = g["plaq_adj"]
    plaq_signed_vecs = g["plaq_signed_vecs"]
    B_inc = g["B_inc"]

    cluster_data: Dict[int, Dict] = {}
    for A in A_values:
        if A == 1:
            cluster = [0]
        else:
            visited = {0}
            queue = [0]
            cluster = [0]
            while queue and len(cluster) < A:
                p = queue.pop(0)
                for nbr in plaq_adj.get(p, []):
                    if nbr not in visited:
                        visited.add(nbr)
                        cluster.append(nbr)
                        queue.append(nbr)
                        if len(cluster) >= A:
                            break
            if len(cluster) < A:
                continue
        # Build LO surface with proper orientations
        n_LO = np.zeros(P, dtype=np.int64)
        n_LO[cluster[0]] = 1
        running_boundary = d_1[cluster[0]].astype(np.float64).copy()
        used_edges = set(np.where(d_1[cluster[0]] != 0)[0].tolist())
        for p in cluster[1:]:
            cand = d_1[p].astype(np.float64).copy()
            cand_edges = set(np.where(cand != 0)[0].tolist())
            shared = used_edges & cand_edges
            sign = 1
            if shared:
                e = next(iter(shared))
                if running_boundary[e] * cand[e] > 0:
                    sign = -1
            n_LO[p] = sign
            running_boundary = running_boundary + sign * cand
            used_edges = set(np.where(running_boundary != 0)[0].tolist())
        boundary = running_boundary.astype(np.int64)
        # Validate
        if not is_valid_closed_loop(boundary, B_inc):
            continue
        if not np.all(np.abs(boundary) <= 1):
            continue
        cluster_data[A] = {
            "cluster": cluster,
            "n_LO": n_LO,
            "w": boundary,
        }

    # For each A, scan beta
    rows = []
    truncation = max(A_values) + truncation_extra + 1
    for beta in beta_grid:
        W_by_A: Dict[int, float] = {}
        for A in sorted(cluster_data.keys()):
            cd = cluster_data[A]
            r = compute_W_normalized(d_1, cd["n_LO"], cd["w"], cycles, beta, truncation)
            W_by_A[A] = r["W_NLO"]
        sigmas_eff: List[Dict] = []
        As = sorted(W_by_A.keys())
        for i in range(len(As) - 1):
            A1, A2 = As[i], As[i + 1]
            W1, W2 = W_by_A[A1], W_by_A[A2]
            if W1 > 0 and W2 > 0 and A2 > A1:
                s_eff = -np.log(W2 / W1) / (A2 - A1)
            else:
                s_eff = float("nan")
            sigmas_eff.append({
                "A_pair": [A1, A2],
                "W_A1": float(W1),
                "W_A2": float(W2),
                "sigma_eff": float(s_eff) if not np.isnan(s_eff) else None,
            })
        rows.append({
            "beta": float(beta),
            "sigma_LO": float(sigma_beta(beta)),
            "W_by_A": {int(A): float(W) for A, W in W_by_A.items()},
            "sigmas_eff": sigmas_eff,
        })
    return {
        "A_values": sorted(cluster_data.keys()),
        "rows": rows,
    }


# =============================================================================
# Main
# =============================================================================

def main():
    out = {
        "sprint": "XCWG-E normalized NLO character expansion (2026-05-16)",
        "formula": "<W>_NLO = N_NLO(w; beta) / Z_NLO(beta)",
        "sigma_NLO_def": "sigma_NLO := -log <W>_NLO / A_min",
        "rationale": (
            "LO predicts <W>=(I_1/I_0)^A for all beta. NLO adds (I_1/I_0)^k_delta "
            "corrections from closed 2-cycles to BOTH N and Z. The proper ratio "
            "remains <1 at all beta if 3D-like (positive sigma), or hits 1 at "
            "some beta_c if 4D-like (deconfinement transition)."
        ),
    }
    beta_grid = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0]
    out["beta_grid"] = beta_grid

    # n_max=2
    print("\n" + "#" * 72)
    print("# Normalized NLO at n_max = 2")
    print("#" * 72)
    out["n_max_2"] = run_nlo(
        n_max=2, plaq_cap=200,
        cycle_max_count=30,
        beta_grid=beta_grid,
        cycle_time_budget_sec=120.0,
    )

    # Creutz-ratio at n_max=2 to cross-check
    print("\n" + "-" * 72)
    print("Creutz-ratio cross-check at n_max=2")
    print("-" * 72)
    g2 = build_d1_and_plaquettes(2, plaq_cap=200)
    # Use the cycles we already found in the run_nlo call (they were enumerated)
    # but rebuild a focused set
    cycles_n2 = enumerate_closed_2cycles_up_to_size(g2["d_1"], g2["plaq_adj"],
                                                    size_target=6,
                                                    max_count=20,
                                                    time_budget_sec=60.0)
    print(f"  Using {len(cycles_n2)} cycles of size 6 for Creutz analysis")
    cr = creutz_ratios(g2, cycles_n2, A_values=[1, 2, 3, 4],
                       beta_grid=beta_grid, truncation_extra=7)
    out["n_max_2_creutz"] = cr
    print(f"  A values: {cr['A_values']}")
    print(f"  {'beta':>8}  {'sig_LO':>11}  ", end="")
    for i in range(len(cr["A_values"]) - 1):
        A1 = cr["A_values"][i]
        A2 = cr["A_values"][i + 1]
        print(f"{'sig_eff('+str(A1)+'->'+str(A2)+')':>16}  ", end="")
    print()
    for row in cr["rows"]:
        beta = row["beta"]
        s_lo = row["sigma_LO"]
        print(f"  {beta:>8.4g}  {s_lo:>11.4e}  ", end="")
        for s_eff_data in row["sigmas_eff"]:
            s_eff = s_eff_data["sigma_eff"]
            if s_eff is not None:
                print(f"{s_eff:>16.4e}  ", end="")
            else:
                print(f"{'nan':>16}  ", end="")
        print()

    # n_max=3 (limited)
    print("\n" + "#" * 72)
    print("# Normalized NLO at n_max = 3 (limited)")
    print("#" * 72)
    try:
        out["n_max_3"] = run_nlo(
            n_max=3, plaq_cap=200,
            cycle_max_count=10,
            beta_grid=beta_grid,
            cycle_time_budget_sec=90.0,
        )
    except Exception as e:
        out["n_max_3"] = {"error": str(e)}
        print(f"  n_max=3 failed: {e}")

    # Final verdict
    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)
    if "verdict" in out["n_max_2"]:
        v2 = out["n_max_2"]["verdict"]
        print(f"  n_max=2: k_delta = {out['n_max_2']['k_delta']}")
        print(f"  n_max=2: all_positive = {v2['all_positive']}, "
              f"min_sigma = {v2['min_sigma_NLO']:.4e} at beta = {v2['min_beta']:.4g}")
        print(f"  n_max=2: zero_crossing = {v2['has_zero_crossing']}")
    if "n_max_3" in out and "verdict" in out["n_max_3"]:
        v3 = out["n_max_3"]["verdict"]
        print(f"  n_max=3: k_delta = {out['n_max_3'].get('k_delta')}")
        print(f"  n_max=3: all_positive = {v3['all_positive']}, "
              f"min_sigma = {v3['min_sigma_NLO']}")

    # Save
    out_path = os.path.join(_HERE, "data", "xcwg_nlo_character_expansion.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    def _sanitize(obj):
        if isinstance(obj, dict):
            return {str(k): _sanitize(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_sanitize(x) for x in obj]
        if isinstance(obj, (np.integer, np.int64, np.int32, np.int8)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    out_clean = _sanitize(out)
    with open(out_path, "w") as f:
        json.dump(out_clean, f, indent=2, default=str)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
