"""
XCWG-H Polyakov-loop analog on time-extended Rule B graph G_B x C_{N_t} (2026-05-16).

Sixth-witness diagnostic for the 3D compact U(1) verdict in Paper 41 v3.

Theoretical setup
=================

The Polyakov loop in finite-temperature lattice gauge theory is
   P(x) = prod_{t=0}^{N_t-1} U_{(x,t) -> (x,t+1)}
the path-ordered product of *temporal* link variables around the
Euclidean-time circle of length N_t.  Its expectation value is the
canonical deconfinement order parameter:
   <P> = 0    confined  (center Z_N unbroken)
   <P> != 0   deconfined (center Z_N spontaneously broken)

For 3D compact U(1):
   Polyakov 1977 PLB 59, 82 / 1975 PLB 72, 477:
   permanent confinement at all beta, all T.  <P> = 0 universally.
   No deconfinement transition: the would-be Kaluza-Klein reduction to
   2D U(1) is also confined (Wilson 1974), so 3D U(1) at finite T stays
   confined down to T -> infinity.

For 4D compact U(1):
   BMK 1977 / Guth 1980: transition at beta_c ~ 1.01 at T=0.  At finite T,
   the deconfinement temperature T_c depends on beta; the theory enters
   the Coulomb / deconfined phase above some N_t-dependent threshold.

This script builds the time-extended product graph G_B x C_{N_t},
defines the Wilson U(1) action on its (spatial + mixed temporal-spatial)
plaquettes, and measures <P> via Metropolis-Hastings MC at multiple beta
and N_t.

Expected outcome (if Rule B Wilson U(1) is indeed 3D compact U(1)):
   <P>(beta, N_t) = 0 within MC noise at ALL tested (beta, N_t).
   Two-point correlator <P(v) P*(w)> exponentially decaying with
   spatial graph distance.

If <P> != 0 at some (beta, N_t), that would CONTRADICT the 3D compact
U(1) verdict and indicate a Kaluza-Klein deconfinement transition.

Output
======
    debug/data/xcwg_polyakov_loop.json
    debug/xcwg_polyakov_loop_memo.md          (separately)
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter, defaultdict
from typing import Dict, List, Tuple, Optional

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from geovac.ihara_zeta_dirac import build_dirac_s3_graph  # noqa: E402

sys.stdout.reconfigure(line_buffering=True)

TWO_PI = 2.0 * np.pi


# =============================================================================
# Product graph G_B x C_{N_t} construction
# =============================================================================

def build_product_graph(n_max: int, N_t: int) -> Dict:
    """Construct the product graph G_B x C_{N_t}.

    Vertices are pairs (v, t) for v in V_B, t in {0, ..., N_t-1}, indexed
    as i = v * N_t + t (so v = i // N_t, t = i % N_t).

    Edges come in two types:
      SPATIAL: (v, t) - (w, t) for (v,w) in E_B, each t.  E_B * N_t edges.
      TEMPORAL: (v, t) - (v, t+1 mod N_t) for each v in V_B, each t.
        V_B * N_t edges.

    All edges are oriented lex-min in the global index (so each undirected
    edge appears once); we tag whether it's spatial or temporal.

    Plaquettes come in two types:
      SPATIAL: inherit each 4-plaquette of G_B, one copy at each time t.
      TEMPORAL: rectangular plaquettes (v,t) - (w,t) - (w,t+1) - (v,t+1)
        for each spatial edge (v,w) and each t.  E_B * N_t plaquettes.
    """
    # Build Rule B
    A_B, labels_B, deg_B, desc = build_dirac_s3_graph(n_max, "B")
    V_B = int(A_B.shape[0])
    # Spatial edges in V_B (undirected, lex-min)
    spatial_edges_B: List[Tuple[int, int]] = []
    for u in range(V_B):
        for v in range(u + 1, V_B):
            if A_B[u, v]:
                spatial_edges_B.append((u, v))
    E_B = len(spatial_edges_B)

    # Enumerate Rule B 4-plaquettes (primitive non-backtracking)
    plaqs_B = enumerate_4_plaquettes(A_B)
    P_B = len(plaqs_B)

    # Total product-graph dimensions
    V = V_B * N_t
    E_spatial = E_B * N_t
    E_temporal = V_B * N_t
    E = E_spatial + E_temporal
    P_spatial = P_B * N_t
    P_temporal = E_B * N_t
    P = P_spatial + P_temporal

    # Vertex index: (v, t) -> v * N_t + t
    def vi(v: int, t: int) -> int:
        return v * N_t + (t % N_t)

    # Build edge list and edge type
    # edges[e] = (a, b) with a < b in the product-graph index
    # edge_type[e] = 's' or 't'
    edges: List[Tuple[int, int]] = []
    edge_types: List[str] = []
    # Map from edge (a,b) with a<b -> e_idx
    edge_index: Dict[Tuple[int, int], int] = {}
    # For spatial edges, also record which Rule-B edge and time slice
    # so we can build the boundary maps efficiently
    spatial_edge_info: List[Tuple[int, int]] = []  # (rule_b_edge_idx, t)
    temporal_edge_info: List[Tuple[int, int]] = []  # (v_in_B, t)  meaning t -> t+1

    # Spatial edges
    for t in range(N_t):
        for k, (u, v) in enumerate(spatial_edges_B):
            a, b = vi(u, t), vi(v, t)
            if a > b:
                a, b = b, a
            e_idx = len(edges)
            edges.append((a, b))
            edge_types.append('s')
            edge_index[(a, b)] = e_idx
            spatial_edge_info.append((k, t))

    # Temporal edges: (v, t) - (v, t+1 mod N_t)
    for v in range(V_B):
        for t in range(N_t):
            a, b = vi(v, t), vi(v, (t + 1) % N_t)
            # For N_t=2, both spatial vertices v have temporal edge (v,0)-(v,1)
            # which is the same edge as (v,1)-(v,0); we'd double-count.
            # The cycle C_2 has 2 vertices and only 1 edge (it's the same edge
            # in both directions).  We handle this by enforcing a < b and
            # checking if already in edge_index.
            if a > b:
                a, b = b, a
            if (a, b) in edge_index:
                # Already counted; this is the N_t=2 case
                continue
            e_idx = len(edges)
            edges.append((a, b))
            edge_types.append('t')
            edge_index[(a, b)] = e_idx
            # Temporal info: this edge goes from (v, t) to (v, t+1 mod N_t)
            # We store t such that vi(v,t) is the smaller endpoint.
            # For N_t >= 3, this is straightforward; for N_t = 2, the edge
            # connects (v,0) and (v,1) directly.
            temporal_edge_info.append((v, t))

    # Recompute counts in case N_t=2 collapsed temporal edges
    E = len(edges)
    E_temporal_actual = len(temporal_edge_info)
    E_spatial_actual = len(spatial_edge_info)
    assert E == E_spatial_actual + E_temporal_actual

    # Build plaquette list with signed boundaries (d_1: P x E)
    # Plaquette boundary signs are determined by traversing the plaquette
    # in a consistent orientation; for each edge in the boundary we record
    # its sign relative to the canonical (a, b) orientation with a < b.
    plaquettes: List[List[Tuple[int, int]]] = []  # each: list of (e_idx, sign)
    plaquette_types: List[str] = []

    # SPATIAL plaquettes: each Rule B 4-plaquette at each time slice.
    # plaq_B is a sequence of vertices (v0, v1, v2, v3) forming a 4-cycle in G_B.
    # The boundary in the product graph at time t is the 4 spatial edges:
    # (v0,t)-(v1,t), (v1,t)-(v2,t), (v2,t)-(v3,t), (v3,t)-(v0,t).
    for t in range(N_t):
        for cyc in plaqs_B:
            # cyc is a tuple of 4 vertex indices in V_B
            boundary: List[Tuple[int, int]] = []
            ok = True
            for k in range(4):
                u, w = cyc[k], cyc[(k + 1) % 4]
                a, b = vi(u, t), vi(w, t)
                if a < b:
                    sign = +1
                else:
                    a, b = b, a
                    sign = -1
                e_idx = edge_index.get((a, b))
                if e_idx is None:
                    ok = False
                    break
                boundary.append((e_idx, sign))
            if ok:
                plaquettes.append(boundary)
                plaquette_types.append('s')

    # TEMPORAL plaquettes: for each spatial edge (u,w) in G_B and each t,
    # the rectangular plaquette
    #   (u,t) -> (w,t) -> (w,t+1) -> (u,t+1) -> (u,t)
    # uses 2 spatial edges (one at time t, one at time t+1) and 2 temporal edges.
    # For N_t=2, the rectangular plaquette has temporal edges (u,0)-(u,1) and
    # (w,0)-(w,1), which are still distinct edges (different spatial endpoints).
    # The two spatial edges (u,w) at t=0 and t=1 are distinct.  So even for
    # N_t=2 we get genuine rectangular plaquettes; the only thing that
    # collapses is one of them (the "front" and "back" rectangles between
    # times 0,1 and 1,0 are the same rectangle traversed in opposite order).
    # We enforce uniqueness by only taking t in {0, ..., N_t-1} but if N_t=2,
    # the t=0 and t=1 rectangles coincide as undirected plaquettes.
    seen_temporal_plaqs: set = set()
    for k, (u, w) in enumerate(spatial_edges_B):
        for t in range(N_t):
            t_next = (t + 1) % N_t
            # 4 edges in order around the rectangle
            # (u,t) -> (w,t):  spatial at time t
            a, b = vi(u, t), vi(w, t)
            if a < b:
                e1, s1 = edge_index[(a, b)], +1
            else:
                e1, s1 = edge_index[(b, a)], -1
            # (w,t) -> (w,t_next): temporal at vertex w
            a, b = vi(w, t), vi(w, t_next)
            if a < b:
                e2, s2 = edge_index[(a, b)], +1
            else:
                e2, s2 = edge_index[(b, a)], -1
            # (w,t_next) -> (u,t_next): spatial at time t_next, reversed
            a, b = vi(w, t_next), vi(u, t_next)
            if a < b:
                e3, s3 = edge_index[(a, b)], +1
            else:
                e3, s3 = edge_index[(b, a)], -1
            # (u,t_next) -> (u,t): temporal at vertex u, reversed
            a, b = vi(u, t_next), vi(u, t)
            if a < b:
                e4, s4 = edge_index[(a, b)], +1
            else:
                e4, s4 = edge_index[(b, a)], -1
            # Canonical key: sorted set of 4 edges
            key = tuple(sorted([e1, e2, e3, e4]))
            if key in seen_temporal_plaqs:
                continue
            seen_temporal_plaqs.add(key)
            plaquettes.append([(e1, s1), (e2, s2), (e3, s3), (e4, s4)])
            plaquette_types.append('t')

    P = len(plaquettes)

    # Build d_1: P x E signed boundary matrix
    d_1 = np.zeros((P, E), dtype=np.int8)
    for p_idx, boundary in enumerate(plaquettes):
        for e_idx, sign in boundary:
            d_1[p_idx, e_idx] += sign

    # Build temporal-edge lookup at each spatial vertex
    # For Polyakov loop at spatial vertex v: collect the temporal edges
    # in cyclic order around the C_{N_t} cycle anchored at v.
    # polyakov_path[v] = list of (e_idx, sign) such that following these
    # in order yields the loop (v,0) -> (v,1) -> ... -> (v, N_t-1) -> (v,0).
    polyakov_path: List[List[Tuple[int, int]]] = []
    for v in range(V_B):
        path: List[Tuple[int, int]] = []
        for t in range(N_t):
            a_orig, b_orig = vi(v, t), vi(v, (t + 1) % N_t)
            a, b = (a_orig, b_orig) if a_orig < b_orig else (b_orig, a_orig)
            sign = +1 if a_orig < b_orig else -1
            e_idx = edge_index[(a, b)]
            path.append((e_idx, sign))
        polyakov_path.append(path)

    # Build spatial-distance matrix on G_B (used for two-point Polyakov correlator)
    spatial_dist = bfs_distance_matrix(A_B)

    return {
        "n_max": n_max,
        "N_t": N_t,
        "V_B": V_B,
        "E_B": E_B,
        "P_B": P_B,
        "V": V,
        "E": E,
        "E_spatial": E_spatial_actual,
        "E_temporal": E_temporal_actual,
        "P": P,
        "P_spatial": sum(1 for t in plaquette_types if t == 's'),
        "P_temporal": sum(1 for t in plaquette_types if t == 't'),
        "edges": edges,
        "edge_types": edge_types,
        "edge_index": edge_index,
        "spatial_edge_info": spatial_edge_info,
        "temporal_edge_info": temporal_edge_info,
        "plaquettes": plaquettes,
        "plaquette_types": plaquette_types,
        "d_1": d_1,
        "polyakov_path": polyakov_path,
        "spatial_dist": spatial_dist,
        "rule_b_desc": desc,
    }


def bfs_distance_matrix(A: np.ndarray) -> np.ndarray:
    """Compute V x V matrix of graph distances by BFS."""
    V = A.shape[0]
    D = np.full((V, V), -1, dtype=np.int32)
    for s in range(V):
        D[s, s] = 0
        queue = [s]
        while queue:
            next_q = []
            for u in queue:
                for w in np.nonzero(A[u])[0]:
                    if D[s, w] < 0:
                        D[s, w] = D[s, u] + 1
                        next_q.append(int(w))
            queue = next_q
    return D


def enumerate_4_plaquettes(A: np.ndarray) -> List[Tuple[int, int, int, int]]:
    """Enumerate primitive 4-cycles in undirected graph A.

    Returns list of 4-tuples (v0, v1, v2, v3) with v0 < min(v1, v2, v3)
    and v0 < v2 (canonicalization); each 4-cycle appears exactly once.
    """
    V = A.shape[0]
    plaqs: List[Tuple[int, int, int, int]] = []
    seen: set = set()
    for v0 in range(V):
        neighbors0 = np.nonzero(A[v0])[0]
        for v1 in neighbors0:
            if v1 <= v0:
                continue
            neighbors1 = np.nonzero(A[v1])[0]
            for v2 in neighbors1:
                if v2 == v0:
                    continue
                neighbors2 = np.nonzero(A[v2])[0]
                for v3 in neighbors2:
                    if v3 == v1 or v3 == v0:
                        continue
                    if A[v3, v0] == 0:
                        continue
                    # Canonical: 4-cycle as set {v0, v1, v2, v3} with vertex order.
                    # Use lexmin rotation/reflection of (v0, v1, v2, v3).
                    cyc = (v0, int(v1), int(v2), int(v3))
                    canon = canonical_4cycle(cyc)
                    if canon not in seen:
                        seen.add(canon)
                        plaqs.append(canon)
    return plaqs


def canonical_4cycle(cyc: Tuple[int, int, int, int]) -> Tuple[int, int, int, int]:
    """Canonicalize 4-cycle under rotation and reflection: lex-min."""
    candidates = []
    for k in range(4):
        rot = cyc[k:] + cyc[:k]
        candidates.append(rot)
        candidates.append(rot[::-1])
    return min(candidates)


# =============================================================================
# Wilson action and Metropolis sweeps
# =============================================================================

def plaquette_angles(d_1: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """theta_P = sum_e sign(P,e) * theta_e for each plaquette."""
    return d_1.astype(np.float64) @ theta


def plaquettes_containing_edge(d_1: np.ndarray) -> List[np.ndarray]:
    """For each edge, list of (plaq_idx, sign) entries."""
    P, E = d_1.shape
    out: List[np.ndarray] = []
    for e in range(E):
        col = d_1[:, e]
        nz = np.nonzero(col)[0]
        out.append(np.stack([nz, col[nz]], axis=1).astype(np.int64))
    return out


def metropolis_sweep(
    d_1: np.ndarray,
    theta: np.ndarray,
    theta_P: np.ndarray,
    edge_data: List[np.ndarray],
    beta: float,
    delta_max: float,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """One Metropolis sweep over all links."""
    E = theta.shape[0]
    n_accept = 0
    for e in range(E):
        delta = rng.uniform(-delta_max, delta_max)
        ed = edge_data[e]
        if ed.size == 0:
            theta[e] = (theta[e] + delta + np.pi) % TWO_PI - np.pi
            n_accept += 1
            continue
        plaq_indices = ed[:, 0]
        signs = ed[:, 1].astype(np.float64)
        theta_P_old = theta_P[plaq_indices]
        theta_P_new = theta_P_old + signs * delta
        dS = beta * np.sum(np.cos(theta_P_old) - np.cos(theta_P_new))
        if dS <= 0.0 or rng.random() < np.exp(-dS):
            theta[e] = (theta[e] + delta + np.pi) % TWO_PI - np.pi
            theta_P[plaq_indices] = theta_P_new
            n_accept += 1
    return theta, theta_P, n_accept


def tune_delta_max(
    d_1: np.ndarray,
    edge_data: List[np.ndarray],
    beta: float,
    rng: np.random.Generator,
    target_accept: float = 0.5,
    n_tune_sweeps: int = 20,
) -> float:
    """Auto-tune delta_max for ~target_accept acceptance."""
    E = d_1.shape[1]
    delta_max = np.pi
    theta = np.zeros(E)
    theta_P = plaquette_angles(d_1, theta)
    for it in range(8):
        accepts = 0
        attempts = 0
        for _ in range(n_tune_sweeps):
            for e in range(E):
                delta = rng.uniform(-delta_max, delta_max)
                ed = edge_data[e]
                if ed.size == 0:
                    theta[e] = (theta[e] + delta + np.pi) % TWO_PI - np.pi
                    accepts += 1
                    attempts += 1
                    continue
                plaq_indices = ed[:, 0]
                signs = ed[:, 1].astype(np.float64)
                theta_P_old = theta_P[plaq_indices]
                theta_P_new = theta_P_old + signs * delta
                dS = beta * np.sum(np.cos(theta_P_old) - np.cos(theta_P_new))
                attempts += 1
                if dS <= 0.0 or rng.random() < np.exp(-dS):
                    theta[e] = (theta[e] + delta + np.pi) % TWO_PI - np.pi
                    theta_P[plaq_indices] = theta_P_new
                    accepts += 1
        rate = accepts / max(attempts, 1)
        if abs(rate - target_accept) < 0.05:
            return delta_max
        if rate > target_accept:
            delta_max = min(np.pi, delta_max * 1.4)
        else:
            delta_max = delta_max * 0.7
    return delta_max


# =============================================================================
# Polyakov loop measurement
# =============================================================================

def measure_polyakov(
    theta: np.ndarray,
    polyakov_path: List[List[Tuple[int, int]]],
) -> np.ndarray:
    """For each spatial vertex v in V_B, compute the Polyakov loop:
        P(v) = exp(i * sum_t sign(P,e_t) * theta(e_t))
    where the sum is over the temporal edges around the loop at vertex v.

    Returns complex array of shape (V_B,) of Polyakov-loop values.
    """
    V_B = len(polyakov_path)
    P_vals = np.zeros(V_B, dtype=np.complex128)
    for v in range(V_B):
        phase = 0.0
        for e_idx, sign in polyakov_path[v]:
            phase += sign * theta[e_idx]
        P_vals[v] = np.exp(1j * phase)
    return P_vals


# =============================================================================
# MC driver at one (beta, N_t)
# =============================================================================

def mc_run_one_point(
    g: Dict,
    beta: float,
    n_therm: int,
    n_sample: int,
    sample_interval: int,
    rng_seed: int,
    verbose: bool = True,
) -> Dict:
    """Run MC at one (beta, N_t), measure <P>, <|P|>, <P(v)P*(w)>."""
    rng = np.random.default_rng(rng_seed)
    d_1 = g["d_1"]
    E = d_1.shape[1]
    V_B = g["V_B"]
    polyakov_path = g["polyakov_path"]
    spatial_dist = g["spatial_dist"]
    edge_data = plaquettes_containing_edge(d_1)

    # Tune
    delta_max = tune_delta_max(d_1, edge_data, beta, np.random.default_rng(rng_seed + 1000),
                                target_accept=0.5, n_tune_sweeps=15)

    # Cold start
    theta = np.zeros(E)
    theta_P = plaquette_angles(d_1, theta)

    # Thermalize
    t0 = time.time()
    n_accept_therm = 0
    for sweep in range(n_therm):
        theta, theta_P, n_acc = metropolis_sweep(
            d_1, theta, theta_P, edge_data, beta, delta_max, rng
        )
        n_accept_therm += n_acc
    t_therm = time.time() - t0
    therm_accept_rate = n_accept_therm / (n_therm * E)

    # Sample
    t0 = time.time()
    # Per-config measurements: P(v) for each v (complex), magnitude, etc.
    P_per_config: List[np.ndarray] = []     # each (V_B,) complex
    # Two-point correlator accumulators by spatial distance
    corr_sum_by_d: Dict[int, complex] = defaultdict(complex)
    corr_count_by_d: Dict[int, int] = defaultdict(int)
    abs_corr_sum_by_d: Dict[int, float] = defaultdict(float)
    n_accept_samp = 0
    for k in range(n_sample):
        for _ in range(sample_interval):
            theta, theta_P, n_acc = metropolis_sweep(
                d_1, theta, theta_P, edge_data, beta, delta_max, rng
            )
            n_accept_samp += n_acc
        # Measure
        P_v = measure_polyakov(theta, polyakov_path)
        P_per_config.append(P_v)
        # Two-point correlator: <P(v) P*(w)> as function of d(v,w)
        for v in range(V_B):
            for w in range(V_B):
                d = int(spatial_dist[v, w])
                if d < 0:
                    continue
                term = P_v[v] * np.conj(P_v[w])
                corr_sum_by_d[d] += term
                corr_count_by_d[d] += 1
                abs_corr_sum_by_d[d] += abs(term)
    t_samp = time.time() - t0
    samp_accept_rate = n_accept_samp / (n_sample * sample_interval * E)

    # Aggregate
    P_arr = np.stack(P_per_config, axis=0)  # (n_sample, V_B)
    # <P> averaged over lattice and config
    P_mean = complex(np.mean(P_arr))
    re_P_mean = float(P_mean.real)
    im_P_mean = float(P_mean.imag)
    # |<P>| (gauge-invariant magnitude of the mean)
    abs_P_mean = abs(P_mean)
    # <|P|> averaged: each P(v) is a phase on unit circle, |P(v)|=1 exactly
    # since it's exp(i*phase).  So <|P|> = 1 by construction; what's
    # informative is <|<P>|_lattice>: the magnitude of the spatial average
    # within each configuration.
    abs_lat_avg_per_config = np.abs(np.mean(P_arr, axis=1))   # (n_sample,)
    abs_lat_avg_mean = float(np.mean(abs_lat_avg_per_config))
    abs_lat_avg_se = float(np.std(abs_lat_avg_per_config, ddof=1)
                             / np.sqrt(n_sample)) if n_sample > 1 else 0.0
    # Per-config <P>_lat
    lat_avg_per_config = np.mean(P_arr, axis=1)   # (n_sample,) complex
    re_lat_avg_mean = float(lat_avg_per_config.real.mean())
    im_lat_avg_mean = float(lat_avg_per_config.imag.mean())
    re_lat_avg_se = float(lat_avg_per_config.real.std(ddof=1)
                            / np.sqrt(n_sample)) if n_sample > 1 else 0.0
    im_lat_avg_se = float(lat_avg_per_config.imag.std(ddof=1)
                            / np.sqrt(n_sample)) if n_sample > 1 else 0.0
    # Re P_mean / Im P_mean SE
    # Block standard error: 10 blocks
    n_blocks = min(10, n_sample)
    block_size = max(n_sample // n_blocks, 1)
    block_re_P = np.array([
        P_arr[i*block_size:(i+1)*block_size].real.mean()
        for i in range(n_blocks)
    ])
    block_im_P = np.array([
        P_arr[i*block_size:(i+1)*block_size].imag.mean()
        for i in range(n_blocks)
    ])
    re_P_se = float(block_re_P.std(ddof=1) / np.sqrt(n_blocks)) if n_blocks > 1 else 0.0
    im_P_se = float(block_im_P.std(ddof=1) / np.sqrt(n_blocks)) if n_blocks > 1 else 0.0
    # |<P>| SE: standard error of the magnitude of the spatially+temporally averaged P
    block_lat_P = np.array([
        complex(P_arr[i*block_size:(i+1)*block_size].mean())
        for i in range(n_blocks)
    ])
    abs_P_se = float(np.abs(block_lat_P).std(ddof=1) / np.sqrt(n_blocks)) if n_blocks > 1 else 0.0

    # Correlator by distance
    correlator: Dict[int, Dict] = {}
    distances = sorted(corr_count_by_d.keys())
    for d in distances:
        c_avg = corr_sum_by_d[d] / corr_count_by_d[d] / n_sample
        abs_c_avg = abs_corr_sum_by_d[d] / corr_count_by_d[d] / n_sample
        correlator[d] = {
            "n_pairs_per_config": corr_count_by_d[d] // n_sample,
            "n_total_pairs": corr_count_by_d[d],
            "re_corr": float(c_avg.real),
            "im_corr": float(c_avg.imag),
            "abs_corr": float(abs(c_avg)),
            "mean_abs_term": float(abs_c_avg),
        }

    if verbose:
        print(f"    [beta={beta:.3g}, N_t={g['N_t']}]: "
              f"<P>={re_P_mean:+.4f}+{im_P_mean:+.4f}j +/-({re_P_se:.4f},{im_P_se:.4f}); "
              f"|<P>|={abs_P_mean:.4f}+/-{abs_P_se:.4f}; "
              f"<|<P>_lat|>={abs_lat_avg_mean:.4f}+/-{abs_lat_avg_se:.4f}; "
              f"therm/samp={t_therm:.1f}/{t_samp:.1f}s; accept={samp_accept_rate:.3f}")

    return {
        "beta": float(beta),
        "N_t": int(g["N_t"]),
        "n_max": int(g["n_max"]),
        "re_P_mean": re_P_mean,
        "im_P_mean": im_P_mean,
        "re_P_se": re_P_se,
        "im_P_se": im_P_se,
        "abs_P_mean": float(abs_P_mean),
        "abs_P_se": abs_P_se,
        "abs_lat_avg_mean": abs_lat_avg_mean,
        "abs_lat_avg_se": abs_lat_avg_se,
        "re_lat_avg_mean": re_lat_avg_mean,
        "im_lat_avg_mean": im_lat_avg_mean,
        "re_lat_avg_se": re_lat_avg_se,
        "im_lat_avg_se": im_lat_avg_se,
        "delta_max": float(delta_max),
        "therm_accept_rate": float(therm_accept_rate),
        "samp_accept_rate": float(samp_accept_rate),
        "n_therm": int(n_therm),
        "n_sample": int(n_sample),
        "sample_interval": int(sample_interval),
        "therm_time_sec": float(t_therm),
        "samp_time_sec": float(t_samp),
        "correlator_by_spatial_distance": correlator,
    }


# =============================================================================
# Main driver
# =============================================================================

def main():
    out: Dict = {
        "sprint": "XCWG-H Polyakov loop on G_B x C_{N_t} (2026-05-16)",
        "method": (
            "Construct product graph G_B x C_{N_t} (Rule B at n_max=2 "
            "cross-multiplied with cyclic chain of length N_t). "
            "Wilson U(1) action on (spatial + mixed temporal) plaquettes. "
            "Metropolis-Hastings MC, measure Polyakov loop "
            "P(v) = prod_t U_{(v,t)->(v,t+1)} at each spatial vertex, "
            "and two-point correlator <P(v) P*(w)> vs spatial distance."
        ),
        "physics_target": (
            "3D compact U(1) is confined at all beta, all T (Polyakov 1977). "
            "If Rule B Wilson U(1) is indeed 3D compact U(1), <P>=0 within "
            "MC noise at ALL (beta, N_t).  Deviation would indicate "
            "Kaluza-Klein deconfinement and contradict the 3D verdict."
        ),
    }

    # Build product graph at multiple N_t values; n_max=2 is the primary
    # operating point given XCWG-F machinery.
    n_max = 2
    N_t_values = [2, 4, 8]
    beta_values = [0.3, 1.0, 3.0]

    out["n_max"] = n_max
    out["N_t_values"] = N_t_values
    out["beta_values"] = beta_values

    all_results: List[Dict] = []
    graph_summaries: List[Dict] = []

    for N_t in N_t_values:
        print(f"\n{'='*72}\nBuilding product graph G_B x C_{{N_t}} at n_max={n_max}, N_t={N_t}\n{'='*72}")
        t_build = time.time()
        g = build_product_graph(n_max, N_t)
        t_build = time.time() - t_build
        print(f"  V={g['V']}, E={g['E']} (spatial={g['E_spatial']}, temporal={g['E_temporal']}), "
              f"P={g['P']} (spatial={g['P_spatial']}, temporal={g['P_temporal']})")
        print(f"  V_B={g['V_B']}, E_B={g['E_B']}, P_B={g['P_B']}; build time={t_build:.2f}s")

        graph_summaries.append({
            "n_max": n_max,
            "N_t": N_t,
            "V": g["V"],
            "E": g["E"],
            "E_spatial": g["E_spatial"],
            "E_temporal": g["E_temporal"],
            "P": g["P"],
            "P_spatial": g["P_spatial"],
            "P_temporal": g["P_temporal"],
            "V_B": g["V_B"],
            "E_B": g["E_B"],
            "P_B": g["P_B"],
            "build_time_sec": t_build,
        })

        # Sanity: spectral structure
        d_1 = g["d_1"]
        rank_d_1 = int(np.linalg.matrix_rank(d_1.astype(np.float64)))
        print(f"  rank(d_1) = {rank_d_1}")

        for beta in beta_values:
            r = mc_run_one_point(
                g, beta=beta,
                n_therm=2000, n_sample=800, sample_interval=20,
                rng_seed=10_000 + 1000 * N_t + int(100 * beta),
                verbose=True,
            )
            all_results.append(r)

    out["graph_summaries"] = graph_summaries
    out["all_results"] = all_results

    # ----------------------------------------------------------------
    # Verdict: persistent confinement (3D U(1) signature)?
    # ----------------------------------------------------------------
    # NOTE on N_t=2 degeneracy:
    #   The cyclic chain C_2 has 2 vertices but only 1 undirected edge per
    #   spatial vertex.  The Polyakov loop at vertex v is then
    #       P(v) = U_{(v,0)->(v,1)} * U_{(v,1)->(v,0)} = exp(i*theta_e)*exp(-i*theta_e) = 1
    #   identically for any gauge configuration.  This is a kinematic identity,
    #   NOT a deconfinement signal: the "loop" backtracks on a single edge so
    #   it carries no net winding.  N_t=2 is structurally degenerate and we
    #   exclude it from the confinement verdict.  The first nontrivial N_t is 3,
    #   and we use N_t in {4, 8} as the operating points.
    print(f"\n{'='*72}\nVERDICT (seventh witness: Polyakov-loop deconfinement test)\n{'='*72}")
    verdict_summary: List[Dict] = []
    nondegenerate_confined = True
    for r in all_results:
        N_t = r["N_t"]; beta = r["beta"]
        re_P, im_P = r["re_P_mean"], r["im_P_mean"]
        re_se, im_se = r["re_P_se"], r["im_P_se"]
        abs_P = r["abs_P_mean"]
        abs_se = r["abs_P_se"]
        # Confined: Re P and Im P both within ~3-sigma of zero, AND
        # |<P>| within ~3-sigma of zero.
        re_zero = abs(re_P) <= 3 * max(re_se, 1e-6)
        im_zero = abs(im_P) <= 3 * max(im_se, 1e-6)
        abs_zero = abs_P <= 3 * max(abs_se, 1e-6)
        # N_t=2 is degenerate (trivial Polyakov loop = 1 by construction)
        degenerate = (N_t == 2)
        if degenerate:
            confined_status = "DEGENERATE"
            confined = None
        else:
            confined = re_zero and im_zero and abs_zero
            confined_status = "CONFINED" if confined else "DECONFINED"
            if not confined:
                nondegenerate_confined = False
        verdict_summary.append({
            "N_t": N_t, "beta": beta,
            "re_P_mean": re_P, "re_P_se": re_se, "re_zero_within_3sigma": re_zero,
            "im_P_mean": im_P, "im_P_se": im_se, "im_zero_within_3sigma": im_zero,
            "abs_P_mean": abs_P, "abs_P_se": abs_se, "abs_zero_within_3sigma": abs_zero,
            "degenerate_N_t_2": bool(degenerate),
            "confined_verdict": confined,
            "verdict_label": confined_status,
        })
        print(f"  N_t={N_t}, beta={beta}: <P>=({re_P:+.4f}+/-{re_se:.4f}, "
              f"{im_P:+.4f}+/-{im_se:.4f}); |<P>|={abs_P:.4f}+/-{abs_se:.4f}; "
              f"verdict={confined_status}")

    out["verdict_summary"] = verdict_summary
    out["seventh_witness_pass"] = bool(nondegenerate_confined)
    out["seventh_witness_note"] = (
        "PASS criterion: <P>=0 within 3-sigma at all NON-DEGENERATE (beta, N_t). "
        "N_t=2 is structurally degenerate (Polyakov loop is identity by "
        "single-edge backtracking) and excluded from the verdict."
    )
    print(f"\n  Seventh witness pass (non-degenerate persistent confinement): {nondegenerate_confined}")
    print(f"  (N_t=2 excluded as degenerate; N_t in {{4, 8}} are the nontrivial tests)")

    # Save
    def _sanitize(obj):
        if isinstance(obj, dict):
            # Don't include numpy arrays or non-serializable items
            return {str(k): _sanitize(v) for k, v in obj.items()
                    if not (k in ("edges", "edge_types", "edge_index",
                                  "spatial_edge_info", "temporal_edge_info",
                                  "plaquettes", "plaquette_types", "d_1",
                                  "polyakov_path", "spatial_dist"))}
        if isinstance(obj, (list, tuple)):
            return [_sanitize(x) for x in obj]
        if isinstance(obj, (np.integer, np.int64, np.int32, np.int8)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, complex):
            return {"real": float(obj.real), "imag": float(obj.imag)}
        return obj

    out_path = os.path.join(_HERE, "data", "xcwg_polyakov_loop.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(_sanitize(out), f, indent=2, default=str)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
