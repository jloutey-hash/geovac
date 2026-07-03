"""
XCWG Wilson-loop scaling sprint (2026-05-15).

Extends the XCWG-B Track 3 Wilson-loop pilot from n_max=2 to n_max in {3, 4, 5}
and fits the scaling exponent alpha from <S(L)> ~ L^alpha at L in {4, 6, 8}.

Second witness for "Rule B Wilson U(1) is 3D compact U(1) on a non-cubic graph":
    alpha ~ 2  -> AREA LAW (confined, expected for 3D compact U(1))
    alpha ~ 1  -> PERIMETER LAW (deconfined, anomalous 4D-like physics)
    1 < a < 2  -> intermediate / finite-size crossover

First witness (MK beta_eff -> 0 from Track B1) already passed.

Uses K = d_1^T d_1 (Track B3 correction; NOT L_1 = B^T B which gives a contact
zero on closed loops since closed loops live in ker B).

Output:
    - tests/wilson_rule_b_support/data/xcwg_wilson_loop_scaling.json
    - debug/xcwg_wilson_loop_scaling_memo.md (separately written)
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections import defaultdict
from typing import Dict, List, Tuple, Iterable

import numpy as np

# Make geovac importable without installing the package.
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))  # repo root (tests/wilson_rule_b_support/ -> two levels up)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from geovac.ihara_zeta_dirac import build_dirac_s3_graph  # noqa: E402


# =============================================================================
# Graph utilities
# =============================================================================

def signed_incidence(A: np.ndarray) -> Tuple[np.ndarray, List[Tuple[int, int]],
                                              Dict[Tuple[int, int], int]]:
    """Build signed incidence matrix B (V x E) with canonical orientation u<v."""
    V = A.shape[0]
    edges: List[Tuple[int, int]] = []
    for u in range(V):
        for v in range(u + 1, V):
            if A[u, v]:
                edges.append((u, v))
    E = len(edges)
    B = np.zeros((V, E), dtype=np.int8)
    edge_idx: Dict[Tuple[int, int], int] = {}
    for e_idx, (u, v) in enumerate(edges):
        B[u, e_idx] = -1
        B[v, e_idx] = +1
        edge_idx[(u, v)] = e_idx
        edge_idx[(v, u)] = e_idx
    return B, edges, edge_idx


def adjacency_list(A: np.ndarray) -> List[List[int]]:
    """Return adjacency list (list of neighbor index lists)."""
    V = A.shape[0]
    return [list(map(int, np.nonzero(A[u])[0])) for u in range(V)]


# =============================================================================
# Primitive closed walk enumeration
# =============================================================================
#
# A primitive closed walk of length L on an undirected graph is a sequence
# of vertices v0, v1, ..., v_{L-1}, v_L = v_0 such that:
#   * (v_i, v_{i+1}) is an edge for every i
#   * v_{i+1} != v_{i-1} (non-backtracking)
#   * v_{L-1} != v_1 (closure is non-backtracking)
#   * the walk cannot be written as the repetition of a shorter walk
#
# We enumerate by DFS from each oriented starting edge and keep one canonical
# representative per cycle (lex-min over rotations + orientation reversal).
#
# At L=8 the count grows substantially. We aim to be memory-efficient by NOT
# materializing all walks in memory; instead we accumulate the action sum S
# in-stream and only keep canonical-representative counts.

def _canonical_walk(walk: Tuple[int, ...]) -> Tuple[int, ...]:
    """Lex-min rotation+reversal of a closed vertex walk (length L+1, v_L = v_0)."""
    # Cyclic walk represented as tuple of L vertices (drop the closing repeat).
    L = len(walk) - 1
    base = walk[:L]
    best = base
    for k in range(L):
        rot = base[k:] + base[:k]
        if rot < best:
            best = rot
    rev = tuple(reversed(base))
    for k in range(L):
        rot = rev[k:] + rev[:k]
        if rot < best:
            best = rot
    return best


def _is_primitive(walk: Tuple[int, ...]) -> bool:
    """Check primitivity of a closed walk represented by L vertices."""
    L = len(walk)
    for d in range(1, L):
        if L % d == 0 and d < L:
            period = walk[:d]
            if all(walk[k*d:(k+1)*d] == period for k in range(L // d)):
                return False
    return True


def enumerate_primitive_closed_walks(
    adj: List[List[int]],
    L: int,
    vertex_order: List[int] = None,
) -> Iterable[Tuple[int, ...]]:
    """Generator of canonical-representative primitive closed walks of length L.

    Yields each canonical walk (tuple of L vertices, no repeated start) once.
    Uses DFS from each starting vertex v0; closes when path returns to v0 and
    length == L; filters by non-backtracking at every step including closure.

    Memory: O(walk depth) = O(L). Output: O(number of canonical walks).

    vertex_order: optional permutation of [0..V-1] specifying DFS start order.
        Useful for randomized sampling when capping (DFS-from-0 systematically
        explores walks involving low-l shells first; shuffling avoids this bias).
    """
    V = len(adj)
    seen = set()
    order = vertex_order if vertex_order is not None else list(range(V))
    for v0 in order:
        # Stack of partial walks: each is list of vertices [v0, v1, ..., v_k]
        # plus the "previous" vertex for non-backtracking check.
        # We do iterative DFS to avoid recursion depth issues.
        stack = [(v0,)]
        while stack:
            path = stack.pop()
            k = len(path)
            cur = path[-1]
            prev = path[-2] if k >= 2 else -1
            if k == L:
                # Closure check: need edge cur -> v0, AND non-backtracking at
                # BOTH closure boundaries:
                #   (a) v0 != path[-2]  : closure edge cur->v0 doesn't reverse
                #                         the previous edge path[-2]->cur
                #   (b) path[-1] != path[1] : closure cur->v0 followed by v0->path[1]
                #                             doesn't backtrack
                if v0 in adj[cur]:
                    if k < 2:
                        continue
                    if v0 == path[-2]:
                        continue  # closure edge backtracks immediately
                    if path[-1] == path[1]:
                        continue  # closure step backtracks the first step
                    closed = tuple(path)
                    if _is_primitive(closed):
                        canon = _canonical_walk(closed + (v0,))
                        if canon not in seen:
                            seen.add(canon)
                            yield canon
                continue
            # Extend
            for w in adj[cur]:
                if w == prev:
                    continue  # backtrack
                stack.append(path + (w,))


def vertex_walk_to_edge_indicator(
    walk: Tuple[int, ...],
    edge_idx: Dict[Tuple[int, int], int],
    n_edges: int,
) -> np.ndarray:
    """Convert a canonical closed walk (L vertices, implicit closure) to a
    signed edge-indicator vector C in R^E.

    The walk's edges are (walk[i], walk[i+1]) for i in 0..L-1, with closure
    edge (walk[L-1], walk[0]). The sign convention: +1 if edge traversed in
    canonical orientation (u < v), -1 if reversed.
    """
    L = len(walk)
    C = np.zeros(n_edges, dtype=np.float64)
    for i in range(L):
        u = walk[i]
        v = walk[(i + 1) % L]
        ei = edge_idx[(u, v)]
        e_canon = min(u, v), max(u, v)
        sign = +1 if (u, v) == e_canon else -1
        C[ei] += sign
    return C


# =============================================================================
# Wilson kinetic operator K = d_1^T d_1
# =============================================================================

def build_K_from_L4_plaquettes(
    plaquettes: List[Tuple[int, ...]],
    edge_idx: Dict[Tuple[int, int], int],
    n_edges: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Build plaquette boundary d_1 (P x E) and Wilson kinetic K = d_1^T d_1.

    plaquettes: list of canonical L=4 closed walks (each a tuple of 4 vertices).
    Returns (d_1, K).
    """
    P = len(plaquettes)
    d_1 = np.zeros((P, n_edges), dtype=np.int8)
    for p, walk in enumerate(plaquettes):
        L = len(walk)
        for i in range(L):
            u = walk[i]
            v = walk[(i + 1) % L]
            ei = edge_idx[(u, v)]
            e_canon = (min(u, v), max(u, v))
            sign = +1 if (u, v) == e_canon else -1
            d_1[p, ei] += sign
    K = (d_1.astype(np.float64).T @ d_1.astype(np.float64))
    return d_1, K


# =============================================================================
# Wilson loop action <S(L)>
# =============================================================================

def average_wilson_action(
    adj: List[List[int]],
    edges: List[Tuple[int, int]],
    edge_idx: Dict[Tuple[int, int], int],
    K_pinv: np.ndarray,
    L: int,
    max_walks_for_stats: int = 200_000,
    report_every: int = 25_000,
) -> Dict:
    """Compute <S(L)> = <C^T K_pinv C> averaged over all primitive closed walks
    of length L on adj.

    Streams the enumeration: accumulates sum_S, sum_S2, count without storing
    walks in memory. Caps at max_walks_for_stats for very large counts (in
    which case the mean is over a representative sample — the canonical
    representative ordering is deterministic so the cap is reproducible).

    Returns dict with mean, std, count, sample_size, walltime.
    """
    n_edges = len(edges)
    t0 = time.time()
    sum_S = 0.0
    sum_S2 = 0.0
    count = 0
    S_min = +np.inf
    S_max = -np.inf
    capped = False
    for walk in enumerate_primitive_closed_walks(adj, L):
        C = vertex_walk_to_edge_indicator(walk, edge_idx, n_edges)
        S = float(C @ K_pinv @ C)
        sum_S += S
        sum_S2 += S * S
        if S < S_min:
            S_min = S
        if S > S_max:
            S_max = S
        count += 1
        if count >= max_walks_for_stats:
            capped = True
            break
        if count % report_every == 0:
            dt = time.time() - t0
            print(f"    L={L}: {count} walks, <S>={sum_S/count:.4f}, elapsed {dt:.1f}s")
    dt = time.time() - t0
    if count == 0:
        return {
            "L": L,
            "count": 0,
            "mean": float("nan"),
            "std": float("nan"),
            "min": float("nan"),
            "max": float("nan"),
            "walltime_sec": dt,
            "capped": False,
        }
    mean = sum_S / count
    var = max(0.0, sum_S2 / count - mean * mean)
    std = float(np.sqrt(var))
    return {
        "L": L,
        "count": count,
        "mean": mean,
        "std": std,
        "min": S_min,
        "max": S_max,
        "walltime_sec": dt,
        "capped": capped,
    }


# =============================================================================
# Scaling-exponent fit
# =============================================================================

def fit_scaling_exponent(L_values: List[int], S_means: List[float]) -> Dict:
    """Fit <S(L)> = A * L^alpha by log-log linear regression.

    Returns dict with alpha, log_A, residual, predicted_S_at_each_L,
    and the area-law/perimeter-law verdicts.

    With n_data points, returns:
        - alpha and std error from least-squares
        - residual sum of squares (log space)
        - r^2 = 1 - SS_res / SS_tot
        - predicted S at each L for cross-check
        - successive ratios S(L_{i+1})/S(L_i)
    """
    L_arr = np.array(L_values, dtype=np.float64)
    S_arr = np.array(S_means, dtype=np.float64)
    if np.any(S_arr <= 0):
        return {
            "alpha": float("nan"),
            "log_A": float("nan"),
            "r_squared": float("nan"),
            "comment": "Non-positive S; cannot log-fit.",
        }
    logL = np.log(L_arr)
    logS = np.log(S_arr)
    # Standard least-squares: logS = log_A + alpha * logL
    # Use numpy.polyfit (degree 1) for slope+intercept, then estimate stderr
    n = len(L_arr)
    if n < 2:
        return {"alpha": float("nan"), "comment": "Need >= 2 points."}
    # Manual OLS for stderr
    L_mean = np.mean(logL)
    S_mean = np.mean(logS)
    SS_xx = np.sum((logL - L_mean) ** 2)
    SS_xy = np.sum((logL - L_mean) * (logS - S_mean))
    alpha = SS_xy / SS_xx
    log_A = S_mean - alpha * L_mean
    # Residuals & R^2
    predicted_logS = log_A + alpha * logL
    SS_res = float(np.sum((logS - predicted_logS) ** 2))
    SS_tot = float(np.sum((logS - S_mean) ** 2))
    r_squared = 1.0 - SS_res / SS_tot if SS_tot > 0 else 1.0
    # Stderr of alpha: sigma_residual / sqrt(SS_xx); sigma_residual = sqrt(SS_res / (n-2))
    if n > 2:
        sigma_resid = float(np.sqrt(SS_res / (n - 2)))
        alpha_stderr = sigma_resid / float(np.sqrt(SS_xx))
    else:
        alpha_stderr = float("nan")  # Can't estimate stderr from 2 points
    predicted_S = np.exp(predicted_logS).tolist()
    ratios = []
    for i in range(1, len(L_arr)):
        ratios.append({
            "L_pair": [int(L_arr[i-1]), int(L_arr[i])],
            "S_ratio": float(S_arr[i] / S_arr[i-1]),
            "L_ratio_to_alpha_2": float((L_arr[i] / L_arr[i-1]) ** 2),
            "L_ratio_to_alpha_1": float(L_arr[i] / L_arr[i-1]),
        })
    return {
        "alpha": float(alpha),
        "alpha_stderr": float(alpha_stderr) if not np.isnan(alpha_stderr) else None,
        "log_A": float(log_A),
        "A": float(np.exp(log_A)),
        "r_squared": float(r_squared),
        "ss_residual_log": SS_res,
        "predicted_S": predicted_S,
        "ratios": ratios,
        "verdict": (
            "AREA LAW" if abs(alpha - 2.0) < 0.3
            else "PERIMETER LAW" if abs(alpha - 1.0) < 0.3
            else "INTERMEDIATE"
        ),
    }


# =============================================================================
# Main driver
# =============================================================================

def run_n_max(n_max: int, target_lengths: List[int],
              caps: Dict[int, int]) -> Dict:
    """Run the Wilson-loop scaling at one n_max.

    caps[L] = max walks to enumerate at length L (stream and cap).
    """
    print(f"\n{'='*70}")
    print(f"n_max = {n_max}")
    print(f"{'='*70}")

    t_build_0 = time.time()
    A, labels, deg, desc = build_dirac_s3_graph(n_max, "B")
    V = A.shape[0]
    E_tot = int(A.sum() // 2)
    print(f"  Rule B graph: V={V}, E={E_tot}, beta_1={E_tot - V + 1}")

    B_inc, edges, edge_idx = signed_incidence(A)
    adj = adjacency_list(A)

    # Step 1: Enumerate L=4 plaquettes (this defines the 2-cells of the
    # Hodge 2-complex; K = d_1^T d_1 is built from them).
    print(f"  enumerating L=4 plaquettes (cap = {caps[4]})...")
    t0 = time.time()
    plaquettes = []
    for walk in enumerate_primitive_closed_walks(adj, 4):
        plaquettes.append(walk)
        if len(plaquettes) >= caps[4]:
            break
    print(f"    found {len(plaquettes)} plaquettes in {time.time()-t0:.1f}s")

    # Step 2: Build d_1 and K = d_1^T d_1
    t0 = time.time()
    d_1, K = build_K_from_L4_plaquettes(plaquettes, edge_idx, E_tot)
    K_rank = int(np.linalg.matrix_rank(K))
    d_1_rank = int(np.linalg.matrix_rank(d_1))
    beta_1 = E_tot - V + 1
    print(f"    d_1 shape={d_1.shape}, rank={d_1_rank}; "
          f"K shape={K.shape}, rank={K_rank}; beta_1={beta_1}")
    if K_rank != beta_1:
        print(f"    WARNING: rank(K)={K_rank} != beta_1={beta_1}")

    # Step 3: K_pinv (Moore-Penrose; works on cokernel of K)
    t0 = time.time()
    K_pinv = np.linalg.pinv(K)
    print(f"    K_pinv computed in {time.time()-t0:.1f}s")

    # Step 4: For each target length L, compute <S(L)>
    by_length = {}
    for L in target_lengths:
        print(f"  computing <S(L={L})> (cap = {caps[L]})...")
        result = average_wilson_action(
            adj, edges, edge_idx, K_pinv, L,
            max_walks_for_stats=caps[L],
        )
        print(f"    L={L}: {result['count']} walks, <S>={result['mean']:.4f} "
              f"+/- {result['std']:.4f}, range [{result['min']:.4f}, {result['max']:.4f}]"
              f"{', CAPPED' if result['capped'] else ''}")
        by_length[L] = result

    # Step 5: Scaling fit
    L_used = [L for L in target_lengths if by_length[L]["count"] > 0]
    S_means = [by_length[L]["mean"] for L in L_used]
    fit = fit_scaling_exponent(L_used, S_means)
    stderr_str = ""
    if fit.get("alpha_stderr") is not None:
        stderr_str = f" +/- {fit['alpha_stderr']:.4f}"
    print(f"  fit alpha = {fit['alpha']:.4f}{stderr_str}, "
          f"R^2 = {fit['r_squared']:.4f}, verdict: {fit['verdict']}")

    return {
        "n_max": n_max,
        "graph": {"V": V, "E": E_tot, "beta_1": beta_1},
        "K": {"shape": list(K.shape), "rank": K_rank, "d_1_rank": d_1_rank,
              "d_1_shape": list(d_1.shape)},
        "plaquettes": {"count": len(plaquettes), "cap_applied": len(plaquettes) >= caps[4]},
        "wilson_action_by_length": by_length,
        "scaling_fit": fit,
        "walltime_total_sec": time.time() - t_build_0,
    }


def main() -> None:
    results = {
        "sprint": "XCWG Wilson-loop scaling (extension of Track 3 pilot)",
        "date": "2026-05-15",
        "Wilson_kinetic_operator": "K = d_1^T d_1 (Track B3 correction)",
        "scaling_form": "<S(L)> ~ A * L^alpha; alpha~2 area, alpha~1 perimeter",
    }

    # First reproduce n_max=2 for the record (already in pilot, but include L=8).
    # Caps chosen to balance compute time:
    #   n_max=2: graph is tiny, enumerate everything
    #   n_max=3: L=8 may give ~10^4-10^5 walks; cap at 100k
    #   n_max=4: L=8 likely ~10^6-10^7; cap at 200k for representative sample
    #   n_max=5: only attempt L=4 and L=6; cap aggressively

    # n_max = 2
    results["n_max_2"] = run_n_max(
        2,
        target_lengths=[4, 6, 8],
        caps={4: 1_000_000, 6: 1_000_000, 8: 1_000_000},
    )

    # n_max = 3
    results["n_max_3"] = run_n_max(
        3,
        target_lengths=[4, 6, 8],
        caps={4: 1_000_000, 6: 1_000_000, 8: 200_000},
    )

    # n_max = 4
    results["n_max_4"] = run_n_max(
        4,
        target_lengths=[4, 6, 8],
        caps={4: 1_000_000, 6: 500_000, 8: 200_000},
    )

    # n_max = 5 (tentative; L=8 likely intractable)
    try:
        results["n_max_5"] = run_n_max(
            5,
            target_lengths=[4, 6],  # skip L=8 by design
            caps={4: 500_000, 6: 200_000, 8: 0},
        )
        results["n_max_5"]["L8_skipped"] = (
            "L=8 skipped at n_max=5: estimated intractable; "
            "L=4 and L=6 only."
        )
    except Exception as e:
        results["n_max_5"] = {"error": str(e), "comment": "n_max=5 failed; partial data only"}

    # Cross-n_max summary
    alphas = {}
    for k in ["n_max_2", "n_max_3", "n_max_4", "n_max_5"]:
        if k in results and "scaling_fit" in results[k]:
            alphas[k] = {
                "alpha": results[k]["scaling_fit"]["alpha"],
                "alpha_stderr": results[k]["scaling_fit"].get("alpha_stderr"),
                "r_squared": results[k]["scaling_fit"]["r_squared"],
                "verdict": results[k]["scaling_fit"]["verdict"],
            }
    results["cross_n_max_alpha_trend"] = alphas

    # Witness assessment
    alphas_clean = [
        v["alpha"] for v in alphas.values()
        if isinstance(v.get("alpha"), float) and not np.isnan(v["alpha"])
    ]
    if alphas_clean:
        # Drift toward 2 (area) or 1 (perimeter)?
        first = alphas_clean[0]
        last = alphas_clean[-1]
        drift = last - first
        witness_passes = abs(last - 2.0) < 0.4
        results["witness_assessment"] = {
            "first_alpha": first,
            "last_alpha": last,
            "drift": drift,
            "witness_passes_area_law": witness_passes,
            "comment": (
                "Witness PASSES area-law verdict" if witness_passes
                else "Witness FAILS area-law verdict (alpha not close to 2)"
            ),
        }

    out_path = os.path.join(_HERE, "data", "xcwg_wilson_loop_scaling.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {out_path}")

    # Final summary
    print("\n" + "=" * 70)
    print("FINAL SCALING TABLE")
    print("=" * 70)
    print(f"{'n_max':<8} {'alpha':<12} {'stderr':<12} {'R^2':<8} {'verdict':<20}")
    for k, v in alphas.items():
        nm = k.replace("n_max_", "")
        a = v.get("alpha", float("nan"))
        se = v.get("alpha_stderr")
        se_str = f"{se:.4f}" if se is not None else "n/a"
        r2 = v.get("r_squared", float("nan"))
        vd = v.get("verdict", "")
        print(f"{nm:<8} {a:<12.4f} {se_str:<12} {r2:<8.4f} {vd:<20}")


if __name__ == "__main__":
    main()
