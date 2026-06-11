"""
Track alpha-D (Phase 4B): Hopf Projection as Graph Morphism

Tests whether the spectral invariants of the GeoVac S^3 graph at n_max=3,
together with its S^1 -> S^3 -> S^2 Hopf quotient (m-fiber collapse), encode
the fine-structure constant K = 137.036 or its Paper 2 decomposition
K = pi*(B + F - Delta) = pi*(42 + pi^2/6 - 1/40) ~ pi*43.6199.

Construction:
  - S^3 graph: GeoVac lattice(Z=1, max_n=3). Nodes = (n,l,m), 14 states.
    Edges from Paper 7 selection rules:
        angular: (n,l,m) <-> (n,l,m+/-1)    [Delta m = +/-1 within (n,l)]
        radial:  (n,l,m) <-> (n+/-1,l,m)    [Delta n = +/-1, same (l,m)]
    (these are the GeometricLattice.adjacency entries.)
  - S^2 quotient graph: collapse the 2l+1 m-substates of each (n,l) into one
    node labeled (n,l). At n_max=3 this gives 6 nodes:
        (1,0) w=1, (2,0) w=1, (2,1) w=3, (3,0) w=1, (3,1) w=3, (3,2) w=5.
    Quotient adjacency:
        A_S2[a,b] = sum of S^3 edges from any m in sector a to any m in sector b.
    The only S^3 edges crossing sectors are the radial T+/- transitions
    (Delta n=+/-1, same (l,m)). Angular L+/- transitions stay inside a sector
    and become self-loops of weight 2*l (counted twice by the symmetric sum,
    one in each direction); we drop those to form the unweighted quotient
    Laplacian and separately record them as the "fiber" contribution.
  - S^1 fiber graphs: for each (n,l) with l>=1, the 2l+1 m-states form a PATH
    graph P_{2l+1} (not a cycle), because the angular transitions are only
    m <-> m+/-1 with no wrap (|m|<=l). Paper 7's selection rules do NOT close
    the m-ring. We therefore compute both path-graph Laplacians (production)
    and cycle-graph Laplacians (alternative) for completeness.

All eigenvalues of these Laplacians are algebraic integers over Q (actually
integers for unweighted path/cycle graphs, and small integers for the S^3
and quotient Laplacians here). The spectral sums use mpmath at 50 dps.

Author: Track alpha-D sub-agent, 2026-04-10
"""

from __future__ import annotations

import json
import sys
import os
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# mpmath for high-precision spectral sums
try:
    import mpmath as mp
    mp.mp.dps = 50
    HAVE_MP = True
except ImportError:
    HAVE_MP = False

# Import GeoVac lattice
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from geovac.lattice import GeometricLattice  # noqa: E402


# ---------- Targets (Paper 2, Sec VIII and Eq. for K) ----------
K_ALPHA = 137.035999084          # reciprocal fine-structure constant
K_OVER_PI = K_ALPHA / np.pi      # ~ 43.619934
B_PAPER = 42.0                   # integer Casimir trace (B = 42)
F_PAPER = (np.pi ** 2) / 6.0     # ~ 1.6449 (zeta(2))
DELTA_PAPER = 1.0 / 40.0         # ~ 0.025
B_PLUS_F = B_PAPER + F_PAPER                     # ~ 43.6449
B_PLUS_F_MINUS_DELTA = B_PLUS_F - DELTA_PAPER    # ~ 43.6199 ~ K/pi
K_FROM_PAPER2 = np.pi * B_PLUS_F_MINUS_DELTA     # ~ 137.0357

TARGETS = {
    "K_alpha_137.036": K_ALPHA,
    "K/pi_43.6199": K_OVER_PI,
    "B_42": B_PAPER,
    "F_pi2/6": F_PAPER,
    "Delta_1/40": DELTA_PAPER,
    "B+F_43.6449": B_PLUS_F,
    "B+F-Delta_43.6199": B_PLUS_F_MINUS_DELTA,
    "pi*(B+F-Delta)_137.0357": K_FROM_PAPER2,
}


# ---------- Utility: graph invariants ----------

def eigvals_sym(L: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Eigenvalues of a symmetric matrix, sorted ascending, tiny parts zeroed."""
    w = np.linalg.eigvalsh(L.astype(np.float64))
    w[np.abs(w) < tol] = 0.0
    return np.sort(w)


def spectral_invariants(name: str, w: np.ndarray) -> dict:
    """
    Compute: zero modes, zeta(1,2,3), log det', det', algebraic connectivity,
    spectral entropy. Uses mpmath for spectral sums.
    """
    w = np.asarray(w, dtype=np.float64)
    zero_modes = int(np.sum(w < 1e-10))
    pos = w[w > 1e-10]
    n_nonzero = len(pos)

    invariants: dict = {
        "name": name,
        "n_nodes": int(len(w)),
        "zero_modes": zero_modes,
        "n_nonzero_eigs": n_nonzero,
        "eigenvalues": [float(x) for x in w],
        "trace": float(np.sum(w)),
    }

    if n_nonzero == 0:
        invariants.update({
            "zeta_1": None, "zeta_2": None, "zeta_3": None,
            "logdet_prime": None, "det_prime": None,
            "alg_connectivity": None,
            "spectral_entropy": None,
            "largest_eig": None,
        })
        return invariants

    if HAVE_MP:
        pos_mp = [mp.mpf(str(float(x))) for x in pos]
        z1 = mp.fsum(1 / lam for lam in pos_mp)
        z2 = mp.fsum(1 / (lam ** 2) for lam in pos_mp)
        z3 = mp.fsum(1 / (lam ** 3) for lam in pos_mp)
        logdet = mp.fsum(mp.log(lam) for lam in pos_mp)
        detp = mp.exp(logdet)
        invariants.update({
            "zeta_1": float(z1),
            "zeta_2": float(z2),
            "zeta_3": float(z3),
            "zeta_1_str": mp.nstr(z1, 20),
            "logdet_prime": float(logdet),
            "logdet_prime_str": mp.nstr(logdet, 20),
            "det_prime": float(detp),
        })
    else:
        invariants.update({
            "zeta_1": float(np.sum(1.0 / pos)),
            "zeta_2": float(np.sum(1.0 / pos ** 2)),
            "zeta_3": float(np.sum(1.0 / pos ** 3)),
            "logdet_prime": float(np.sum(np.log(pos))),
            "det_prime": float(np.exp(np.sum(np.log(pos)))),
        })

    invariants["alg_connectivity"] = float(np.min(pos))
    invariants["largest_eig"] = float(np.max(w))

    # Spectral entropy from p_i = lam_i / sum(lam_j), skipping zeros
    total = float(np.sum(pos))
    p = pos / total
    p = p[p > 0]
    invariants["spectral_entropy"] = float(-np.sum(p * np.log(p)))
    return invariants


def match_targets(value: float, tol_rel: float = 0.02) -> list:
    """Find all targets within tol_rel relative error."""
    hits = []
    for name, t in TARGETS.items():
        if abs(t) < 1e-12:
            continue
        rel = abs(value - t) / abs(t)
        if rel < tol_rel:
            hits.append({"target": name, "target_value": t, "rel_error": rel})
    return hits


# ---------- Build S^3 graph ----------

def build_s3_graph(max_n: int = 3):
    """Build S^3 graph via GeometricLattice. Returns states, A, L, edge list."""
    lat = GeometricLattice(max_n=max_n, nuclear_charge=1)
    A = lat.adjacency.toarray()
    # Ensure symmetric binary adjacency (weights are 1.0 by default)
    A = ((A + A.T) / 2 > 0).astype(np.float64)
    deg = np.sum(A, axis=1)
    L = np.diag(deg) - A
    edges = []
    for i in range(A.shape[0]):
        for j in range(i + 1, A.shape[0]):
            if A[i, j] > 0:
                edges.append([i, j])
    return lat.states, A, L, edges


# ---------- Build S^2 quotient graph ----------

def build_s2_quotient(states, A):
    """
    Collapse m-substates: sectors labeled by (n,l).
    A_S2[a,b] = # of S^3 edges with one endpoint in sector a, the other in b.
    Self-loop contributions (angular L+/-) are stored separately for bookkeeping.
    Standard graph Laplacian of the quotient is computed from the off-diagonal
    intersector count (standard unweighted quotient).
    """
    sector_of = {}
    sector_list = []  # ordered list of sectors
    sector_weight = {}  # multiplicity 2l+1

    for i, (n, l, m) in enumerate(states):
        key = (n, l)
        if key not in sector_of.values():
            sector_list.append(key)
        sector_of[i] = key

    # Ensure sector_list is unique (preserves insertion order)
    seen = set()
    ordered = []
    for i in range(len(states)):
        key = sector_of[i]
        if key not in seen:
            ordered.append(key)
            seen.add(key)
    sector_list = ordered
    sector_index = {s: k for k, s in enumerate(sector_list)}

    for s in sector_list:
        n, l = s
        sector_weight[s] = 2 * l + 1

    n_sec = len(sector_list)
    # Intersector edge count (off-diagonal) and within-sector edge count (diagonal)
    inter = np.zeros((n_sec, n_sec), dtype=np.float64)
    intra = np.zeros(n_sec, dtype=np.float64)  # angular edges within a sector (fiber)
    N = A.shape[0]
    for i in range(N):
        a = sector_index[sector_of[i]]
        for j in range(i + 1, N):
            if A[i, j] > 0:
                b = sector_index[sector_of[j]]
                if a == b:
                    intra[a] += 1.0
                else:
                    inter[a, b] += 1.0
                    inter[b, a] += 1.0

    # Unweighted quotient Laplacian from intersector adjacency
    A_S2 = inter.copy()
    deg = np.sum(A_S2, axis=1)
    L_S2 = np.diag(deg) - A_S2

    return sector_list, sector_weight, A_S2, L_S2, intra


# ---------- Fiber graphs ----------

def path_graph_laplacian(n: int) -> np.ndarray:
    """Laplacian of the path graph P_n (n nodes, n-1 edges)."""
    if n <= 0:
        return np.zeros((0, 0))
    if n == 1:
        return np.zeros((1, 1))
    A = np.zeros((n, n))
    for i in range(n - 1):
        A[i, i + 1] = 1.0
        A[i + 1, i] = 1.0
    D = np.diag(np.sum(A, axis=1))
    return D - A


def cycle_graph_laplacian(n: int) -> np.ndarray:
    """Laplacian of the cycle graph C_n."""
    if n <= 0:
        return np.zeros((0, 0))
    if n == 1:
        return np.zeros((1, 1))
    if n == 2:
        return np.array([[1.0, -1.0], [-1.0, 1.0]])
    A = np.zeros((n, n))
    for i in range(n):
        A[i, (i + 1) % n] = 1.0
        A[(i + 1) % n, i] = 1.0
    D = np.diag(np.sum(A, axis=1))
    return D - A


def fiber_analysis(sector_list):
    """Build path and cycle fiber Laplacians for each (n,l), l>=1."""
    fibers_path = {}
    fibers_cycle = {}
    for (n, l) in sector_list:
        if l >= 1:
            sz = 2 * l + 1
            fibers_path[(n, l)] = path_graph_laplacian(sz)
            fibers_cycle[(n, l)] = cycle_graph_laplacian(sz)
    return fibers_path, fibers_cycle


# ---------- Main analysis ----------

def run():
    out = {
        "meta": {
            "track": "alpha-D",
            "phase": "4B",
            "max_n": 3,
            "targets": TARGETS,
            "mpmath_dps": 50 if HAVE_MP else None,
        }
    }

    # --- Subtask 1: S^3 graph ---
    states, A_S3, L_S3, edges = build_s3_graph(max_n=3)
    out["s3_graph"] = {
        "nodes": [list(s) for s in states],
        "n_nodes": len(states),
        "edges": edges,
        "n_edges": len(edges),
        "degree_sequence": [int(x) for x in np.sum(A_S3, axis=1)],
    }
    w_S3 = eigvals_sym(L_S3)
    inv_S3 = spectral_invariants("S3 (n_max=3)", w_S3)
    out["s3_invariants"] = inv_S3

    # Verify topological identity: graph eigenvalues should include 0
    # (the constant mode) if graph is connected.
    assert inv_S3["zero_modes"] >= 1

    # --- Subtask 2: S^2 quotient ---
    sector_list, sector_weight, A_S2, L_S2, intra = build_s2_quotient(states, A_S3)
    out["s2_quotient"] = {
        "sectors": [list(s) for s in sector_list],
        "sector_weights": {str(list(s)): int(sector_weight[s]) for s in sector_list},
        "A_S2": A_S2.tolist(),
        "intra_sector_edge_counts": intra.tolist(),
        "description": (
            "Nodes are (n,l). Off-diagonal entries = # of S^3 edges crossing sectors. "
            "Intra_sector_edge_counts records the # of angular L+/- edges inside each "
            "sector (which become the S^1 fiber)."
        ),
    }
    w_S2 = eigvals_sym(L_S2)
    inv_S2 = spectral_invariants("S2 quotient (n_max=3)", w_S2)
    out["s2_invariants"] = inv_S2

    # --- Subtask 3: S^1 fibers ---
    fibers_path, fibers_cycle = fiber_analysis(sector_list)
    fiber_results_path = {}
    fiber_results_cycle = {}
    sum_logdet_path = 0.0
    sum_logdet_cycle = 0.0

    for key, L_f in fibers_path.items():
        w = eigvals_sym(L_f)
        inv = spectral_invariants(f"fiber_path_{key}", w)
        fiber_results_path[str(list(key))] = inv
        if inv["logdet_prime"] is not None:
            sum_logdet_path += inv["logdet_prime"]

    for key, L_f in fibers_cycle.items():
        w = eigvals_sym(L_f)
        inv = spectral_invariants(f"fiber_cycle_{key}", w)
        fiber_results_cycle[str(list(key))] = inv
        if inv["logdet_prime"] is not None:
            sum_logdet_cycle += inv["logdet_prime"]

    out["fibers_path"] = fiber_results_path
    out["fibers_cycle"] = fiber_results_cycle
    out["fiber_logdet_sum_path"] = sum_logdet_path
    out["fiber_logdet_sum_cycle"] = sum_logdet_cycle

    # --- Subtask 4: comparison / combinations ---
    # Named scalars to test against targets.
    candidates = {}

    def add(name, value):
        try:
            candidates[name] = float(value)
        except Exception:
            candidates[name] = None

    # Raw invariants
    add("S3_trace", inv_S3["trace"])
    add("S3_logdet", inv_S3["logdet_prime"])
    add("S3_detp", inv_S3["det_prime"])
    add("S3_zeta1", inv_S3["zeta_1"])
    add("S3_zeta2", inv_S3["zeta_2"])
    add("S3_zeta3", inv_S3["zeta_3"])
    add("S3_entropy", inv_S3["spectral_entropy"])
    add("S3_alg_conn", inv_S3["alg_connectivity"])
    add("S3_largest_eig", inv_S3["largest_eig"])

    add("S2_trace", inv_S2["trace"])
    add("S2_logdet", inv_S2["logdet_prime"])
    add("S2_detp", inv_S2["det_prime"])
    add("S2_zeta1", inv_S2["zeta_1"])
    add("S2_zeta2", inv_S2["zeta_2"])
    add("S2_zeta3", inv_S2["zeta_3"])
    add("S2_entropy", inv_S2["spectral_entropy"])
    add("S2_alg_conn", inv_S2["alg_connectivity"])
    add("S2_largest_eig", inv_S2["largest_eig"])

    add("fiber_logdet_sum_path", sum_logdet_path)
    add("fiber_logdet_sum_cycle", sum_logdet_cycle)

    # Ratios, differences, products
    if inv_S2["logdet_prime"] and inv_S3["logdet_prime"]:
        add("S3_over_S2_logdet", inv_S3["logdet_prime"] / inv_S2["logdet_prime"])
        add("S3_minus_S2_logdet", inv_S3["logdet_prime"] - inv_S2["logdet_prime"])
    if inv_S2["zeta_1"] and inv_S3["zeta_1"]:
        add("S3_over_S2_zeta1", inv_S3["zeta_1"] / inv_S2["zeta_1"])
        add("S3_minus_S2_zeta1", inv_S3["zeta_1"] - inv_S2["zeta_1"])
    if inv_S2["det_prime"] and inv_S3["det_prime"]:
        add("S3_detp_over_S2_detp", inv_S3["det_prime"] / inv_S2["det_prime"])

    # Weighted quotient trace: use 2l+1 multiplicities
    weighted_s2_trace = sum(
        sector_weight[s] * v
        for s, v in zip(sector_list, np.sum(A_S2, axis=1))
    )
    add("S2_weighted_degree_trace", weighted_s2_trace)

    # Paper 2 style combinations
    if inv_S2["logdet_prime"] is not None:
        x = inv_S2["logdet_prime"]
        add("S2_logdet_times_pi", x * np.pi)
        add("S2_logdet_times_pi2_6", x * F_PAPER)
        add("S2_logdet_times_42", x * 42.0)
        add("S2_logdet_plus_42", x + 42.0)
        add("S2_logdet_plus_pi2_6", x + F_PAPER)
    if inv_S3["logdet_prime"] is not None:
        x = inv_S3["logdet_prime"]
        add("S3_logdet_times_pi", x * np.pi)
        add("S3_logdet_times_42", x * 42.0)
        add("S3_logdet_plus_42", x + 42.0)
    if inv_S3["trace"] is not None:
        x = inv_S3["trace"]
        add("S3_trace_times_pi", x * np.pi)
    if inv_S2["trace"] is not None:
        x = inv_S2["trace"]
        add("S2_trace_times_pi", x * np.pi)

    # Degeneracy-weighted Casimir trace (the Paper 2 "B" object)
    # B_formal = sum_{(n,l)} (2l+1) * l*(l+1), truncated at n_max=3
    # The Casimir of SO(3) on orbital l is l(l+1); the fiber degeneracy is 2l+1.
    casimir_sum = 0
    casimir_breakdown = []
    for (n, l) in sector_list:
        c = (2 * l + 1) * l * (l + 1)
        casimir_sum += c
        casimir_breakdown.append({"sector": [n, l], "2l+1": 2 * l + 1,
                                   "l(l+1)": l * (l + 1), "contribution": c})
    out["casimir_weighted_trace"] = {
        "value": casimir_sum,
        "target_B_42": B_PAPER,
        "match": bool(casimir_sum == 42),
        "breakdown": casimir_breakdown,
        "description": (
            "Sum over all (n,l) sectors at n_max=3 of (2l+1)*l(l+1). "
            "This is the degeneracy-weighted Casimir trace of the SO(3) orbital "
            "angular momentum operator, truncated at the n_max=3 shell. "
            "It equals exactly 42 -- the 'B = 42' of Paper 2."
        ),
    }
    add("casimir_weighted_trace", casimir_sum)
    add("casimir_times_pi", casimir_sum * np.pi)
    add("casimir_times_pi_plus_F", casimir_sum * np.pi + F_PAPER * np.pi)
    add("pi_times_(casimir+F-Delta)",
        np.pi * (casimir_sum + F_PAPER - DELTA_PAPER))

    # Multiplicity-weighted (normalized) Laplacian on the S^2 quotient:
    # L_w = M^{-1/2} L_S2 M^{-1/2} with M = diag(2l+1)
    mult = np.array([sector_weight[s] for s in sector_list], dtype=float)
    Minv12 = np.diag(1.0 / np.sqrt(mult))
    L_S2_weighted = Minv12 @ L_S2 @ Minv12
    w_S2_w = eigvals_sym(L_S2_weighted)
    inv_S2_w = spectral_invariants("S2 quotient (2l+1)-weighted", w_S2_w)
    out["s2_weighted_invariants"] = inv_S2_w
    add("S2w_trace", inv_S2_w["trace"])
    add("S2w_zeta1", inv_S2_w["zeta_1"])
    add("S2w_logdet", inv_S2_w["logdet_prime"])

    # Graph exchange constant (several definitions)
    gec = {}
    if inv_S3["zeta_1"] is not None and inv_S2["zeta_1"] is not None:
        gec["zeta1_diff"] = inv_S3["zeta_1"] - inv_S2["zeta_1"]
    if inv_S3["logdet_prime"] is not None and inv_S2["logdet_prime"] is not None:
        gec["logdet_diff"] = inv_S3["logdet_prime"] - inv_S2["logdet_prime"]
    # Von Neumann entropy diff
    gec["vn_entropy_diff"] = inv_S3["spectral_entropy"] - inv_S2["spectral_entropy"]
    gec["fiber_logdet_sum_path"] = sum_logdet_path
    gec["fiber_logdet_sum_cycle"] = sum_logdet_cycle
    out["graph_exchange_constant"] = gec
    for k, v in gec.items():
        add(f"gec_{k}", v)

    # Near-miss table
    near_misses = []
    for name, v in candidates.items():
        if v is None:
            continue
        hits = match_targets(v, tol_rel=0.05)
        for h in hits:
            near_misses.append({
                "candidate": name,
                "candidate_value": v,
                **h,
            })
    near_misses.sort(key=lambda r: r["rel_error"])
    out["candidates"] = candidates
    out["near_misses"] = near_misses

    # --- Write JSON ---
    outdir = PROJECT_ROOT / "debug" / "data" / "track_alpha_phase4b"
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "track_d_graph_morphism.json", "w") as f:
        json.dump(out, f, indent=2, default=str)

    # --- Write markdown analysis ---
    write_markdown(out, outdir / "track_d_analysis.md")

    print_summary(out)
    return out


def write_markdown(out, path):
    lines = []
    lines.append("# Track alpha-D: Hopf Projection as Graph Morphism")
    lines.append("")
    lines.append("## Construction")
    lines.append("")
    lines.append(
        "- **S^3 graph**: GeoVac `GeometricLattice(max_n=3, Z=1)`. Nodes are "
        "(n,l,m) for 1<=n<=3, 0<=l<n, -l<=m<=l (14 nodes). Edges from Paper 7 "
        "selection rules: angular Delta m=+/-1 within (n,l), radial Delta n=+/-1 "
        "with same (l,m)."
    )
    lines.append(
        "- **S^2 quotient**: collapse 2l+1 m-substates of each (n,l) into one "
        "node labeled (n,l); 6 nodes with multiplicities 1,1,3,1,3,5 (sum=14). "
        "Intersector adjacency counts S^3 edges crossing sectors (these are "
        "exactly the radial Delta n=+/-1 edges). Intra-sector edges are angular "
        "L+/- transitions and become the fiber."
    )
    lines.append(
        "- **S^1 fibers**: each (n,l) sector with l>=1 has 2l+1 m-states "
        "connected by angular L+/- transitions. Paper 7 selection rules give a "
        "PATH graph P_{2l+1}, not a cycle. Cycle alternative computed for "
        "completeness."
    )
    lines.append("")

    lines.append("## S^3 graph at n_max = 3")
    lines.append("")
    lines.append(f"Nodes: {out['s3_graph']['n_nodes']}, edges: {out['s3_graph']['n_edges']}")
    lines.append("")
    lines.append(f"Laplacian eigenvalues (S^3): {out['s3_invariants']['eigenvalues']}")
    lines.append("")
    lines.append("## S^2 quotient")
    lines.append("")
    lines.append("Intersector adjacency A_S2:")
    lines.append("```")
    for row in out["s2_quotient"]["A_S2"]:
        lines.append("  " + " ".join(f"{x:4.0f}" for x in row))
    lines.append("```")
    lines.append(f"Laplacian eigenvalues (S^2): {out['s2_invariants']['eigenvalues']}")
    lines.append("")

    lines.append("## Fiber Laplacian (path-graph) eigenvalue sums")
    lines.append("")
    for k, v in out["fibers_path"].items():
        lines.append(f"- {k}: eigs={v['eigenvalues']}, logdet'={v['logdet_prime']}")
    lines.append(f"Sum of fiber log det' (path): {out['fiber_logdet_sum_path']}")
    lines.append("")

    lines.append("## Spectral invariants comparison")
    lines.append("")
    lines.append("| Invariant | S^3 | S^2 |")
    lines.append("|---|---|---|")
    for key in ["n_nodes", "zero_modes", "trace", "zeta_1", "zeta_2", "zeta_3",
                "logdet_prime", "det_prime", "alg_connectivity",
                "largest_eig", "spectral_entropy"]:
        a = out["s3_invariants"].get(key)
        b = out["s2_invariants"].get(key)
        lines.append(f"| {key} | {a} | {b} |")
    lines.append("")

    lines.append("## Graph exchange constant candidates")
    lines.append("")
    for k, v in out["graph_exchange_constant"].items():
        lines.append(f"- {k} = {v}")
    lines.append("")

    lines.append("## Candidates vs targets (rel < 5%)")
    lines.append("")
    if out["near_misses"]:
        lines.append("| Candidate | Value | Target | Target value | Rel err |")
        lines.append("|---|---|---|---|---|")
        for row in out["near_misses"]:
            lines.append(
                f"| {row['candidate']} | {row['candidate_value']:.6g} | "
                f"{row['target']} | {row['target_value']:.6g} | "
                f"{row['rel_error']:.3e} |"
            )
    else:
        lines.append("No candidates matched any target within 5% relative error.")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(
        "AMBIGUOUS / TAUTOLOGICAL. Interpretation of the results:"
    )
    lines.append("")
    lines.append(
        "1. **The raw graph Laplacian spectral invariants of L_{S^3} and L_{S^2} "
        "do NOT encode K.** Eigenvalues are small integers and golden ratios "
        "(Fibonacci values from the n-chain paths). logdet'(L_{S^3}) = 7.208, "
        "logdet'(L_{S^2}) = 2.890, zeta_1(L_{S^3}) = 7.7, zeta_1(L_{S^2}) = 1.5. "
        "None of these hit K, K/pi, B=42, or B+F-Delta to better than 3% in any "
        "natural combination."
    )
    lines.append("")
    lines.append(
        "2. **The 'graph exchange constant' (logdet/zeta/entropy differences "
        "between S^3 and S^2) does not hit any target.** logdet_diff = 4.317, "
        "zeta1_diff = 6.2, VN entropy diff = 1.338, fiber logdet sum = 3.807. "
        "None are close to 42, 43.62, pi^2/6, or 137."
    )
    lines.append("")
    lines.append(
        "3. **The integer B = 42 IS recovered exactly**, but as the "
        "degeneracy-weighted Casimir sum Sum_{(n,l)} (2l+1)*l(l+1) truncated at "
        "n_max=3, summing contributions 0+0+6+0+6+30 = 42. This is the "
        "quotient-space Casimir trace, not a graph Laplacian spectral "
        "invariant. It reproduces Paper 2's B = 42 construction exactly and "
        "therefore gives K = pi*(B+F-Delta) to 4.8e-7 relative error. **But "
        "this is Paper 2's existing formula restated in graph language, not a "
        "new derivation.** The tautology is: B=42 comes from SO(3) "
        "representation theory of the S^2 base, which is identical to "
        "weighting (n,l) sectors by (2l+1)*l(l+1). Calling the sectors 'nodes "
        "of the quotient graph' does not change the computation."
    )
    lines.append("")
    lines.append(
        "4. **Conclusion**: the Hopf-quotient spectral route does not produce "
        "an independent derivation of K at n_max=3. The graph morphism is "
        "well-defined, but its spectral invariants live in a different "
        "algebraic world (small integers, Fibonacci) from the Paper 2 targets "
        "(42, pi^2/6, 1/40, 137). The only combination that hits is the one "
        "that bypasses the graph Laplacian entirely and uses the SO(3) Casimir "
        "directly -- which is Paper 2's original construction."
    )
    lines.append("")
    lines.append(
        "Recommendation: **CLEAN NEGATIVE for the 'spectral invariants of the "
        "graph quotient encode K' hypothesis.** The honest partial result is "
        "that the S^2-quotient node set carries the SO(3) Casimir weighting, "
        "so B=42 is naturally expressed as a weighted trace on the quotient; "
        "this is a re-labeling of Paper 2, not new physics. No independent "
        "prediction of K or K/pi from Laplacian spectra was obtained."
    )
    with open(path, "w") as f:
        f.write("\n".join(lines))


def print_summary(out):
    print("=" * 70)
    print("Track alpha-D: Hopf Projection as Graph Morphism")
    print("=" * 70)
    print()
    print(f"S^3 graph: {out['s3_graph']['n_nodes']} nodes, "
          f"{out['s3_graph']['n_edges']} edges")
    print(f"S^3 eigenvalues: {out['s3_invariants']['eigenvalues']}")
    print(f"S^3 logdet' = {out['s3_invariants']['logdet_prime']}")
    print(f"S^3 zeta(1) = {out['s3_invariants']['zeta_1']}")
    print()
    print(f"S^2 quotient: {out['s2_invariants']['n_nodes']} nodes")
    print(f"S^2 eigenvalues: {out['s2_invariants']['eigenvalues']}")
    print(f"S^2 logdet' = {out['s2_invariants']['logdet_prime']}")
    print(f"S^2 zeta(1) = {out['s2_invariants']['zeta_1']}")
    print()
    print(f"Fiber logdet sum (path): {out['fiber_logdet_sum_path']}")
    print(f"Fiber logdet sum (cycle): {out['fiber_logdet_sum_cycle']}")
    print()
    print("Near-misses (rel < 5%):")
    for row in out["near_misses"][:15]:
        print(f"  {row['candidate']:35s} = {row['candidate_value']:12.6f}  -> "
              f"{row['target']:25s} rel={row['rel_error']:.3e}")
    if not out["near_misses"]:
        print("  (none)")
    print()


if __name__ == "__main__":
    run()
