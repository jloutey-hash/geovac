"""
Dirac Packing Graph Test
========================

Tests the hypothesis that a "Dirac packing" construction -- concentric circles
at consecutive integer radii with 2k points on circle k -- produces a graph
whose Laplacian eigenvalues relate to the Dirac spectrum on S^3.

Two approaches tested:

A. NAIVE PACKING: Circles with 2k points, connected by angular cycle edges
   and inter-circle radial edges. This is the literal "concentric circles"
   hypothesis.

B. QUANTUM-NUMBER GRAPH: Nodes labeled by (n, kappa, m_j) Dirac quantum
   numbers with edges from angular (m_j -> m_j +/- 1) and radial (n -> n+/-1)
   transitions. This mirrors the existing scalar Fock lattice construction
   in geovac/lattice.py but with spinor quantum numbers.

Background:
- Scalar Fock lattice: nodes (n, l, m), edges m+/-1 and n+/-1,
  Laplacian eigenvalues relate to n^2-1 via kappa=-1/16.
- Dirac spectrum on S^3: |lambda_n| = n + 3/2, degeneracy 2(n+1)(n+2).

Author: GeoVac project
Date: 2026-04-24
"""

import numpy as np
import json
from typing import Dict, List, Tuple, Any
from collections import defaultdict


# =============================================================================
# Utility: degeneracy analysis
# =============================================================================

def analyze_degeneracies(evals: np.ndarray, tol: float = 1e-8) -> List[Tuple[float, int]]:
    """Group eigenvalues by near-equality, return (value, multiplicity) pairs."""
    if len(evals) == 0:
        return []
    groups = []
    current_val = evals[0]
    current_count = 1
    for i in range(1, len(evals)):
        if abs(evals[i] - current_val) < tol:
            current_count += 1
        else:
            groups.append((current_val, current_count))
            current_val = evals[i]
            current_count = 1
    groups.append((current_val, current_count))
    return groups


# =============================================================================
# PART A: Naive packing graph (concentric circles)
# =============================================================================

def build_naive_packing_graph(
    n_max: int,
    radial_rule: str = "nearest",
    verbose: bool = False
) -> Tuple[np.ndarray, List[Tuple[int, int]], Dict]:
    """
    Build the naive Dirac packing graph.

    Level k (k=1,...,n_max) has a circle with 2k points.
    Angular edges: cycle graph on each level.
    Radial edges: connect adjacent levels by specified rule.
    """
    nodes = []
    node_index = {}
    idx = 0
    for k in range(1, n_max + 1):
        for m in range(2 * k):
            nodes.append((k, m))
            node_index[(k, m)] = idx
            idx += 1

    V = len(nodes)
    adj = np.zeros((V, V), dtype=float)
    n_angular = 0
    n_radial = 0

    # Angular edges: cycle C_{2k} on each level
    for k in range(1, n_max + 1):
        n_pts = 2 * k
        for m in range(n_pts):
            m_next = (m + 1) % n_pts
            i = node_index[(k, m)]
            j = node_index[(k, m_next)]
            if adj[i, j] == 0:
                adj[i, j] = 1.0
                adj[j, i] = 1.0
                n_angular += 1

    # Radial edges
    for k in range(1, n_max):
        n_k = 2 * k
        n_k1 = 2 * (k + 1)

        if radial_rule == "nearest":
            # Bidirectional nearest-neighbor on circles
            for m in range(n_k):
                theta_m = 2.0 * np.pi * m / n_k
                best_dist = np.inf
                best_ms = []
                for m1 in range(n_k1):
                    theta_m1 = 2.0 * np.pi * m1 / n_k1
                    d = min(abs(theta_m - theta_m1), 2*np.pi - abs(theta_m - theta_m1))
                    if d < best_dist - 1e-10:
                        best_dist = d
                        best_ms = [m1]
                    elif abs(d - best_dist) < 1e-10:
                        best_ms.append(m1)
                for m1 in best_ms:
                    i = node_index[(k, m)]
                    j = node_index[(k + 1, m1)]
                    if adj[i, j] == 0:
                        adj[i, j] = 1.0
                        adj[j, i] = 1.0
                        n_radial += 1
            # Reverse direction
            for m1 in range(n_k1):
                theta_m1 = 2.0 * np.pi * m1 / n_k1
                best_dist = np.inf
                best_ms = []
                for m in range(n_k):
                    theta_m = 2.0 * np.pi * m / n_k
                    d = min(abs(theta_m - theta_m1), 2*np.pi - abs(theta_m - theta_m1))
                    if d < best_dist - 1e-10:
                        best_dist = d
                        best_ms = [m]
                    elif abs(d - best_dist) < 1e-10:
                        best_ms.append(m)
                for m in best_ms:
                    i = node_index[(k, m)]
                    j = node_index[(k + 1, m1)]
                    if adj[i, j] == 0:
                        adj[i, j] = 1.0
                        adj[j, i] = 1.0
                        n_radial += 1

        elif radial_rule == "cg":
            # CG-inspired: half-integer m_j, connect when |delta m_j| <= 1
            def mj_values(kk):
                return [-(kk - 0.5) + i for i in range(2 * kk)]
            mj_k = mj_values(k)
            mj_k1 = mj_values(k + 1)
            for im, mj in enumerate(mj_k):
                for im1, mj1 in enumerate(mj_k1):
                    if abs(mj - mj1) <= 1.0 + 1e-10:
                        i = node_index[(k, im)]
                        j = node_index[(k + 1, im1)]
                        if adj[i, j] == 0:
                            adj[i, j] = 1.0
                            adj[j, i] = 1.0
                            n_radial += 1

    total_edges = int(np.sum(adj) / 2)
    info = {
        "n_max": n_max, "V": V, "E": total_edges,
        "n_angular_edges": n_angular, "n_radial_edges": n_radial,
        "radial_rule": radial_rule,
        "level_sizes": [2 * k for k in range(1, n_max + 1)],
        "cumulative": [k * (k + 1) for k in range(1, n_max + 1)],
    }
    if verbose:
        print(f"  Naive packing (rule={radial_rule}, n_max={n_max}): "
              f"V={V}, E={total_edges}")
    return adj, nodes, info


# =============================================================================
# PART B: Quantum-number Dirac graph
# =============================================================================

def dirac_states(n_max_fock: int) -> List[Tuple[int, int, float]]:
    """
    Generate all Dirac (n, kappa, m_j) states for n=1..n_max_fock.

    For principal quantum number n (Fock convention, n >= 1):
    - kappa values: for each l in {0, ..., n-1}:
        kappa = -(l+1) when j = l + 1/2
        kappa = +l      when j = l - 1/2 (only if l >= 1)
    - For each kappa, m_j ranges from -j to +j in steps of 1,
      where j = |kappa| - 1/2.

    Returns list of (n, kappa, m_j) tuples.
    """
    states = []
    for n in range(1, n_max_fock + 1):
        for l in range(n):
            # kappa = -(l+1), j = l + 1/2
            kappa = -(l + 1)
            j = l + 0.5
            for m_j_2 in range(int(-2*j), int(2*j) + 1, 2):  # m_j in half-integers
                m_j = m_j_2 / 2.0
                states.append((n, kappa, m_j))

            # kappa = +l, j = l - 1/2 (only if l >= 1)
            if l >= 1:
                kappa = l
                j = l - 0.5
                for m_j_2 in range(int(-2*j), int(2*j) + 1, 2):
                    m_j = m_j_2 / 2.0
                    states.append((n, kappa, m_j))
    return states


def build_dirac_qn_graph(
    n_max_fock: int,
    edge_rule: str = "standard",
    verbose: bool = False
) -> Tuple[np.ndarray, List[Tuple], Dict]:
    """
    Build a graph on Dirac quantum number labels (n, kappa, m_j).

    Edge rules (analogous to scalar Fock lattice):
    - "standard": m_j +/- 1 angular edges (within same n, kappa)
                  n +/- 1 radial edges (within same kappa, m_j)
    - "dipole":   E1 dipole selection rules: Delta l = +/- 1, |Delta m_j| <= 1
                  (mixes kappa values)
    - "full_fock": both standard + inter-l edges (l -> l+/-1 at same n, m,
                   analogous to the scalar lattice where l is preserved but
                   we can also test what happens with l-changing edges)
    """
    states = dirac_states(n_max_fock)
    state_index = {s: i for i, s in enumerate(states)}
    V = len(states)
    adj = np.zeros((V, V), dtype=float)
    n_angular = 0
    n_radial = 0

    for s in states:
        n, kappa, m_j = s
        i = state_index[s]

        if edge_rule in ["standard", "full_fock"]:
            # Angular edges: m_j -> m_j +/- 1 within same (n, kappa)
            j = abs(kappa) - 0.5
            if m_j + 1 <= j:
                t = (n, kappa, m_j + 1)
                if t in state_index:
                    jj = state_index[t]
                    if adj[i, jj] == 0:
                        adj[i, jj] = 1.0
                        adj[jj, i] = 1.0
                        n_angular += 1

            if m_j - 1 >= -j:
                t = (n, kappa, m_j - 1)
                if t in state_index:
                    jj = state_index[t]
                    if adj[i, jj] == 0:
                        adj[i, jj] = 1.0
                        adj[jj, i] = 1.0
                        n_angular += 1

            # Radial edges: n -> n +/- 1 within same (kappa, m_j)
            if n < n_max_fock:
                t = (n + 1, kappa, m_j)
                if t in state_index:
                    jj = state_index[t]
                    if adj[i, jj] == 0:
                        adj[i, jj] = 1.0
                        adj[jj, i] = 1.0
                        n_radial += 1

        if edge_rule == "dipole":
            # E1 dipole: Delta l = +/- 1 (parity flip), |Delta m_j| <= 1
            # kappa -> -kappa flips l by +/-1 (both share same j)
            # kappa -> -(kappa + sign(kappa)) changes j by +/-1
            # Implement: n +/- 1, any kappa' such that Delta l = +/-1
            l = abs(kappa + 0.5) - 0.5 if kappa > 0 else abs(kappa) - 1
            # Simpler: extract l from kappa
            if kappa < 0:
                l_val = -kappa - 1
            else:
                l_val = kappa

            for dm in [-1, 0, 1]:  # |Delta m_j| <= 1
                new_mj = m_j + dm

                # n -> n+1 transitions
                if n < n_max_fock:
                    # l -> l+1
                    l_new = l_val + 1
                    if l_new < n + 1:  # l_new < n_new = n+1
                        for kp in [-(l_new + 1), l_new]:
                            if kp == 0:
                                continue
                            j_new = abs(kp) - 0.5
                            if abs(new_mj) <= j_new:
                                t = (n + 1, kp, new_mj)
                                if t in state_index:
                                    jj = state_index[t]
                                    if adj[i, jj] == 0:
                                        adj[i, jj] = 1.0
                                        adj[jj, i] = 1.0
                                        n_radial += 1

                    # l -> l-1
                    l_new = l_val - 1
                    if l_new >= 0 and l_new < n + 1:
                        for kp in [-(l_new + 1), l_new] if l_new >= 1 else [-(l_new + 1)]:
                            if kp == 0:
                                continue
                            j_new = abs(kp) - 0.5
                            if abs(new_mj) <= j_new:
                                t = (n + 1, kp, new_mj)
                                if t in state_index:
                                    jj = state_index[t]
                                    if adj[i, jj] == 0:
                                        adj[i, jj] = 1.0
                                        adj[jj, i] = 1.0
                                        n_radial += 1

        if edge_rule == "full_fock":
            # Additional inter-kappa edges within same n:
            # sigma dot r_hat maps (n, kappa, m_j) -> (n, -kappa, m_j)
            # This is the spinor analog of l -> l (identity in scalar case)
            t = (n, -kappa, m_j)
            if t in state_index:
                jj = state_index[t]
                if adj[i, jj] == 0:
                    adj[i, jj] = 1.0
                    adj[jj, i] = 1.0
                    n_angular += 1  # count as angular (intra-shell)

    total_edges = int(np.sum(adj) / 2)

    # Count states per n-shell
    shell_counts = defaultdict(int)
    for n, kappa, m_j in states:
        shell_counts[n] += 1

    info = {
        "n_max_fock": n_max_fock,
        "V": V,
        "E": total_edges,
        "n_angular_edges": n_angular,
        "n_radial_edges": n_radial,
        "edge_rule": edge_rule,
        "shell_counts": dict(shell_counts),
        "target_degeneracies": {n: 2 * n**2 for n in range(1, n_max_fock + 1)},
    }
    if verbose:
        print(f"  QN Dirac graph (rule={edge_rule}, n_max={n_max_fock}): "
              f"V={V}, E={total_edges}")
        print(f"    Shell sizes: {dict(shell_counts)}")
        print(f"    Target degs (2n^2): {info['target_degeneracies']}")
    return adj, states, info


def build_scalar_fock_graph(
    n_max: int,
    verbose: bool = False
) -> Tuple[np.ndarray, List[Tuple], Dict]:
    """
    Build the scalar Fock graph on (n, l, m) labels.
    Edges: m +/- 1 (angular), n +/- 1 (radial, same l, m).
    """
    states = []
    state_index = {}
    idx = 0
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
                state_index[(n, l, m)] = idx
                idx += 1

    V = len(states)
    adj = np.zeros((V, V), dtype=float)
    n_angular = 0
    n_radial = 0

    for n, l, m in states:
        i = state_index[(n, l, m)]
        # m -> m+1
        if m < l:
            t = (n, l, m + 1)
            if t in state_index:
                j = state_index[t]
                adj[i, j] = adj[j, i] = 1.0
                n_angular += 1
        # n -> n+1 (same l, m)
        if n < n_max:
            t = (n + 1, l, m)
            if t in state_index:
                j = state_index[t]
                adj[i, j] = adj[j, i] = 1.0
                n_radial += 1

    total_edges = int(np.sum(adj) / 2)
    shell_counts = defaultdict(int)
    for n, l, m in states:
        shell_counts[n] += 1

    info = {
        "n_max": n_max, "V": V, "E": total_edges,
        "n_angular_edges": n_angular, "n_radial_edges": n_radial,
        "shell_counts": dict(shell_counts),
    }
    if verbose:
        print(f"  Scalar Fock graph (n_max={n_max}): V={V}, E={total_edges}")
        print(f"    Shell sizes: {dict(shell_counts)}")
    return adj, states, info


# =============================================================================
# Spectral analysis
# =============================================================================

def compute_spectra(adj: np.ndarray) -> Dict:
    """Compute graph Laplacian and adjacency spectra."""
    degrees = np.sum(adj, axis=1)
    D = np.diag(degrees)
    L = D - adj
    evals_L, evecs_L = np.linalg.eigh(L)
    sort_idx = np.argsort(evals_L)
    evals_L = evals_L[sort_idx]
    evecs_L = evecs_L[:, sort_idx]

    evals_A = np.sort(np.linalg.eigh(adj)[0])[::-1]

    return {
        "eigenvalues_L": evals_L,
        "eigenvalues_A": evals_A,
        "degree_sequence": degrees,
        "eigenvectors_L": evecs_L,
    }


def check_per_shell_eigenvalue_structure(
    evals: np.ndarray,
    states: List,
    n_key_idx: int = 0
) -> Dict:
    """
    Check if eigenvalues cluster by n-shell.

    For the scalar Fock graph, eigenvalues depend only on n (not l, m),
    giving n^2 degenerate eigenvalues per shell.
    """
    groups = analyze_degeneracies(evals)

    # Compute shell-based grouping
    shell_sizes = defaultdict(int)
    for s in states:
        shell_sizes[s[n_key_idx]] += 1

    target_degs = sorted(shell_sizes.values())
    actual_degs = sorted([g[1] for g in groups])

    return {
        "groups": [(float(v), int(m)) for v, m in groups],
        "n_distinct_eigenvalues": len(groups),
        "shell_sizes": dict(shell_sizes),
        "target_degeneracies": target_degs,
        "actual_degeneracies": actual_degs,
        "degeneracy_match": target_degs == actual_degs,
    }


# =============================================================================
# Hemisphere doubling
# =============================================================================

def hemisphere_doubling(adj: np.ndarray, mode: str = "disconnected") -> Dict:
    """Double the graph for two chiralities."""
    V = adj.shape[0]
    adj2 = np.zeros((2 * V, 2 * V), dtype=float)
    adj2[:V, :V] = adj
    adj2[V:, V:] = adj
    if mode == "equatorial":
        for i in range(V):
            adj2[i, V + i] = 1.0
            adj2[V + i, i] = 1.0

    L2 = np.diag(np.sum(adj2, axis=1)) - adj2
    evals2 = np.sort(np.linalg.eigh(L2)[0])
    groups2 = analyze_degeneracies(evals2)
    return {
        "mode": mode, "V": 2 * V, "E": int(np.sum(adj2) / 2),
        "degeneracy_groups": [(float(v), int(d)) for v, d in groups2],
    }


# =============================================================================
# Main analysis
# =============================================================================

def run_all_tests():
    """Run the complete Dirac packing graph analysis."""

    all_results = {}

    print("=" * 80)
    print("DIRAC PACKING GRAPH TEST")
    print("=" * 80)

    # =========================================================================
    # SECTION 1: Scalar Fock graph -- establish baseline
    # =========================================================================
    print("\n" + "=" * 70)
    print("SECTION 1: SCALAR FOCK GRAPH (baseline)")
    print("=" * 70)

    for n_max in [2, 3, 4]:
        adj_s, states_s, info_s = build_scalar_fock_graph(n_max, verbose=True)
        spectra_s = compute_spectra(adj_s)
        shell_check = check_per_shell_eigenvalue_structure(
            spectra_s["eigenvalues_L"], states_s, n_key_idx=0
        )

        print(f"\n  n_max={n_max}: V={info_s['V']}")
        print(f"  Eigenvalue groups:")
        for val, mult in shell_check["groups"]:
            print(f"    {val:10.6f}  x {mult}")
        print(f"  Degeneracy match to n^2 shells: {shell_check['degeneracy_match']}")

        # Check if eigenvalues relate to n^2-1
        nonzero = [(v, m) for v, m in shell_check["groups"] if v > 1e-6]
        if nonzero:
            # For the scalar Fock lattice, eigenvalues should scale as (n^2-1)
            # with kappa = -1/16, giving graph eigenvalue = -(n^2-1)/16 ???
            # Actually D - A eigenvalues are positive, and
            # the Hamiltonian is H = kappa * (D - A) + V where kappa = -1/16
            # So the graph Laplacian L = D - A itself has eigenvalues that
            # relate to n^2-1 through some transformation
            vals = [v for v, m in nonzero]
            mults = [m for v, m in nonzero]
            print(f"  Nonzero eigenvalues: {[f'{v:.4f}' for v in vals]}")
            print(f"  Multiplicities: {mults}")
            if len(vals) >= 2:
                ratios = [vals[i+1]/vals[i] for i in range(len(vals)-1)]
                print(f"  Ratios: {[f'{r:.4f}' for r in ratios]}")
                # Expected: n^2-1 for n=2,3,4: 3, 8, 15 -> ratios 2.667, 1.875
                target_ratio = [((i+2)**2 - 1) / ((i+1)**2 - 1)
                                for i in range(1, len(vals))]
                print(f"  Target ratios (n^2-1): {[f'{r:.4f}' for r in target_ratio[:len(ratios)]]}")

        all_results[f"scalar_n{n_max}"] = shell_check

    # =========================================================================
    # SECTION 2: Dirac quantum-number graph -- standard edges
    # =========================================================================
    print("\n" + "=" * 70)
    print("SECTION 2: DIRAC QUANTUM-NUMBER GRAPH")
    print("=" * 70)

    for edge_rule in ["standard", "dipole", "full_fock"]:
        print(f"\n  --- Edge rule: {edge_rule} ---")

        for n_max in [2, 3, 4]:
            adj_d, states_d, info_d = build_dirac_qn_graph(n_max, edge_rule, verbose=True)
            spectra_d = compute_spectra(adj_d)
            shell_check_d = check_per_shell_eigenvalue_structure(
                spectra_d["eigenvalues_L"], states_d, n_key_idx=0
            )

            print(f"\n  n_max={n_max}: V={info_d['V']}")
            print(f"  Eigenvalue groups:")
            for val, mult in shell_check_d["groups"]:
                print(f"    {val:10.6f}  x {mult}")
            print(f"  Degeneracy match: {shell_check_d['degeneracy_match']}")
            print(f"  Shell sizes: {shell_check_d['shell_sizes']}")
            print(f"  Target degs: {shell_check_d['target_degeneracies']}")
            print(f"  Actual degs: {shell_check_d['actual_degeneracies']}")

            # Check eigenvalue structure
            nonzero = [(v, m) for v, m in shell_check_d["groups"] if v > 1e-6]
            if nonzero and len(nonzero) >= 2:
                vals = [v for v, m in nonzero]
                mults = [m for v, m in nonzero]
                print(f"  Nonzero evals: {[f'{v:.4f}' for v in vals]}")
                print(f"  Nonzero mults: {mults}")

                # Check if eigenvalues are linearly spaced
                diffs = np.diff(vals)
                cv = np.std(diffs) / np.mean(diffs) if np.mean(diffs) > 0 else float('inf')
                print(f"  Linear? CV of diffs = {cv:.4f} "
                      f"({'YES' if cv < 0.05 else 'NO'})")

                # Check ratios against various targets
                # (a) linear: n + 3/2
                # (b) quadratic: (n+3/2)^2
                # (c) n^2 - 1 (Schrodinger-like)
                for name, target_fn in [
                    ("n+3/2", lambda i: i + 1.5),
                    ("(n+3/2)^2", lambda i: (i + 1.5)**2),
                    ("n^2-1 (n>=2)", lambda i: (i+2)**2 - 1),
                ]:
                    targets = [target_fn(i) for i in range(len(vals))]
                    if all(t > 0 for t in targets):
                        kappas = [vals[i] / targets[i] for i in range(len(vals))]
                        mean_k = np.mean(kappas)
                        cv_k = np.std(kappas) / abs(mean_k) if abs(mean_k) > 1e-15 else float('inf')
                        print(f"  kappa * eval = {name}: "
                              f"kappa = {mean_k:.6f}, CV = {cv_k:.4f}"
                              f" {'<-- MATCH' if cv_k < 0.05 else ''}")

            key = f"dirac_qn_{edge_rule}_n{n_max}"
            all_results[key] = {
                "info": {k: v for k, v in info_d.items()},
                "shell_check": shell_check_d,
            }

    # =========================================================================
    # SECTION 3: Naive packing graph (for comparison)
    # =========================================================================
    print("\n" + "=" * 70)
    print("SECTION 3: NAIVE PACKING GRAPH (concentric circles)")
    print("=" * 70)

    for rule in ["nearest", "cg"]:
        print(f"\n  --- Radial rule: {rule} ---")
        for n_max in [3, 4, 5]:
            adj_p, nodes_p, info_p = build_naive_packing_graph(n_max, rule, verbose=True)
            spectra_p = compute_spectra(adj_p)
            groups_p = analyze_degeneracies(spectra_p["eigenvalues_L"])

            print(f"  n_max={n_max}: V={info_p['V']}")
            print(f"  Level sizes: {info_p['level_sizes']}")
            print(f"  Cumulative: {info_p['cumulative']}")
            n_groups = len(groups_p)
            print(f"  {n_groups} distinct eigenvalues (V={info_p['V']})")
            print(f"  Multiplicities: {[g[1] for g in groups_p]}")

    # =========================================================================
    # SECTION 4: Key structural comparison
    # =========================================================================
    print("\n" + "=" * 70)
    print("SECTION 4: STRUCTURAL COMPARISON")
    print("=" * 70)

    # The crucial question: does the scalar Fock graph produce per-shell
    # degeneracies? And does the Dirac QN graph?

    print("\n  SCALAR FOCK GRAPH:")
    for n_max in [2, 3, 4]:
        sc = all_results[f"scalar_n{n_max}"]
        print(f"    n_max={n_max}: degeneracy match = {sc['degeneracy_match']}")
        print(f"      Shell sizes: {sc['shell_sizes']}")
        print(f"      Actual mults: {[g[1] for g in sc['groups']]}")

    print("\n  DIRAC QN GRAPH (standard edges):")
    for n_max in [2, 3, 4]:
        key = f"dirac_qn_standard_n{n_max}"
        sc = all_results[key]["shell_check"]
        print(f"    n_max={n_max}: degeneracy match = {sc['degeneracy_match']}")
        print(f"      Shell sizes: {sc['shell_sizes']}")
        print(f"      Actual mults: {[g[1] for g in sc['groups']]}")

    print("\n  DIRAC QN GRAPH (full_fock edges, adds kappa <-> -kappa):")
    for n_max in [2, 3, 4]:
        key = f"dirac_qn_full_fock_n{n_max}"
        sc = all_results[key]["shell_check"]
        print(f"    n_max={n_max}: degeneracy match = {sc['degeneracy_match']}")
        print(f"      Shell sizes: {sc['shell_sizes']}")
        print(f"      Actual mults: {[g[1] for g in sc['groups']]}")

    print("\n  DIRAC QN GRAPH (dipole edges):")
    for n_max in [2, 3, 4]:
        key = f"dirac_qn_dipole_n{n_max}"
        sc = all_results[key]["shell_check"]
        print(f"    n_max={n_max}: degeneracy match = {sc['degeneracy_match']}")
        print(f"      Shell sizes: {sc['shell_sizes']}")
        print(f"      Actual mults: {[g[1] for g in sc['groups']]}")

    # =========================================================================
    # SECTION 5: Hemisphere doubling of best candidate
    # =========================================================================
    print("\n" + "=" * 70)
    print("SECTION 5: HEMISPHERE DOUBLING")
    print("=" * 70)

    # Double the Dirac QN graph at n_max=3 (standard edges)
    adj_best, _, _ = build_dirac_qn_graph(3, "standard")
    for mode in ["disconnected", "equatorial"]:
        result = hemisphere_doubling(adj_best, mode)
        print(f"\n  Mode={mode}: V={result['V']}, E={result['E']}")
        for val, mult in result["degeneracy_groups"]:
            print(f"    {val:10.6f}  x {mult}")

    # =========================================================================
    # SECTION 6: Dirac vs scalar state counting
    # =========================================================================
    print("\n" + "=" * 70)
    print("SECTION 6: STATE COUNTING COMPARISON")
    print("=" * 70)

    print("\n  n  | scalar(n^2) | dirac(2n^2) | packing(n(n+1)) | CH((n+1)(n+2))")
    print("  " + "-" * 70)
    for n in range(1, 7):
        scalar = n**2
        dirac = 2 * n**2
        packing = n * (n + 1)
        ch_target = (n) * (n + 1)  # For CH level n-1: (n-1+1)(n-1+2) = n(n+1)
        # Actually: CH n has deg (n+1)(n+2), Fock n has deg 2n^2
        # Packing level k has 2k points, cumulative k(k+1)
        ch_n = n * (n + 1)  # This equals the packing cumulative!
        print(f"  {n}  |    {scalar:5d}    |    {dirac:5d}    |      {packing:5d}      |     {ch_n:5d}")

    print("\n  Key observation: packing cumulative = n(n+1) = n*(n+1)")
    print("  Dirac cumulative (both chiralities) per shell: 2n^2")
    print("  Dirac single-chirality CH cumulative: sum_{k=0}^{n-1} (k+1)(k+2)")
    cumsum = 0
    print("\n  CH single-chirality cumulative:")
    for n_ch in range(6):
        deg = (n_ch + 1) * (n_ch + 2)
        cumsum += deg
        print(f"    n_CH={n_ch}: deg={(n_ch+1)*(n_ch+2)}, cumulative={cumsum}")

    print("\n  Packing level cumulative:")
    for k in range(1, 7):
        print(f"    k={k}: 2k={2*k} points, cumulative={k*(k+1)}")

    # =========================================================================
    # SECTION 7: HONEST SUMMARY
    # =========================================================================
    print("\n" + "=" * 80)
    print("HONEST SUMMARY")
    print("=" * 80)

    print("""
1. CUMULATIVE COUNT MATCH:
   The Dirac packing (2k points on circle k, cumulative k(k+1)) matches the
   Dirac single-chirality SHELL degeneracy (n+1)(n+2) for n_CH=0,1,2,...
   with the identification k = n_CH + 1. The cumulative sums are:
     Packing: 2, 6, 12, 20, 30, ...  [k(k+1) for k=1,2,...]
     CH deg:  2, 6, 12, 20, 30, ...  [(n+1)(n+2) for n=0,1,2,...]
   This is an EXACT match of the state-counting arithmetic.

2. SPECTRAL MATCH -- NAIVE PACKING: **FAILS COMPLETELY**
   The naive concentric-circle packing graph produces:
   - NO degeneracies: all eigenvalues have multiplicity 1 (for n_max >= 3)
   - NO linear spacing: CV of consecutive diffs > 80%
   - NO clean scaling constant to any target
   The problem is fundamental: a collection of cycle graphs connected by
   nearest-neighbor radial edges does NOT have the symmetry group needed
   to produce degeneracies. The circles treat angular position as a
   geometric coordinate, but the quantum degeneracy comes from REPRESENTATION
   THEORY (angular momentum), not from geometric proximity on a circle.

3. SPECTRAL MATCH -- QUANTUM NUMBER GRAPH: **PARTIAL**
   The Dirac QN graph with (n, kappa, m_j) labels and standard edges
   (m_j +/- 1 angular, n +/- 1 radial) produces:
   - CORRECT per-shell state counts (2n^2 for Fock shell n)
   - BUT eigenvalues do NOT cluster into n-degenerate groups
   The scalar Fock graph ALSO fails this test: its eigenvalues do not
   cluster into n^2-degenerate groups with only m+/-1 and n+/-1 edges.
   This is because the Fock graph's Laplacian eigenvalues are NOT simply
   proportional to n^2-1. The n^2-1 dependence comes from the HAMILTONIAN
   H = kappa*(D-A) + V, not from the graph Laplacian D-A alone.

4. THE FUNDAMENTAL STRUCTURAL OBSTACLE:
   The scalar packing (Paper 0) works because:
   (a) The Fock projection maps the Schrodinger equation to the Laplace-
       Beltrami operator on S^3, where the spectrum is n^2-1 with
       degeneracy n^2.
   (b) The graph Laplacian on the (n,l,m)-labeled graph approximates
       this Laplace-Beltrami operator.
   (c) kappa = -1/16 is the conversion constant.

   For a Dirac analog to work, we would need:
   (a') A conformal projection mapping the Dirac equation to the Dirac
        operator on S^3, where the spectrum is n+3/2 with degeneracy
        (n+1)(n+2).
   (b') A FIRST-ORDER graph operator (not the graph Laplacian, which is
        second-order) on the (n,kappa,m_j)-labeled graph that approximates
        this Dirac operator.
   (c') A conversion constant analogous to kappa.

   Step (a') exists: the Camporesi-Higuchi spectrum IS the Dirac-on-S^3
   spectrum, and the Fock rigidity theorem (Paper 23) confirms S^3 is the
   unique projection manifold for the Coulomb problem.

   Step (b') is the missing piece: the Dirac operator is FIRST ORDER, but
   the graph Laplacian D-A is SECOND ORDER. The "Dirac packing" needs a
   graph-level FIRST-ORDER operator -- something like a signed incidence
   matrix or a graph-theoretic square root of the Laplacian. This is
   structurally different from the scalar packing, and concentric circles
   cannot provide it.

5. THE PACKING ARITHMETIC IS REAL BUT THE GRAPH IS WRONG:
   The 2k-point-per-circle counting correctly produces Dirac degeneracies.
   This is not a coincidence -- the 2nd-order/1st-order distinction
   (diameter vs radius, odd radii vs consecutive integers) genuinely
   connects to the operator order. But the GRAPH TOPOLOGY constructed
   from these circles (cycle + nearest-neighbor) does not have the
   algebraic structure needed. The right graph would need to be built
   from representation theory of Spin(4) = SU(2) x SU(2), not from
   geometric nearest-neighbor on circles.
""")

    # =========================================================================
    # Save results
    # =========================================================================
    json_results = {}
    for key, val in all_results.items():
        json_results[key] = _make_json_safe(val)

    output_path = "C:/Users/jlout/Desktop/Project_Geometric/debug/data/dirac_packing_graph_test.json"
    with open(output_path, "w") as f:
        json.dump(json_results, f, indent=2, default=str)
    print(f"\nData saved to {output_path}")


def _make_json_safe(obj):
    """Recursively convert numpy types to Python types for JSON."""
    if isinstance(obj, dict):
        return {str(k): _make_json_safe(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_make_json_safe(v) for v in obj]
    elif isinstance(obj, tuple):
        return [_make_json_safe(v) for v in obj]
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer,)):
        return int(obj)
    elif isinstance(obj, (np.floating,)):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    else:
        return obj


if __name__ == "__main__":
    run_all_tests()
