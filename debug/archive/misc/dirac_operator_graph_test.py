"""
Dirac operator on packing graphs: D_graph = [[0, B^T], [B, 0]]

Tests whether the signed incidence matrix of packing graphs produces
eigenvalues related to the Dirac spectrum on S^3: |lambda_n| = n + 3/2.

Graph A: Dirac packing (2k nodes on circle k)
Graph B: Scalar Fock packing (2(2k-1) nodes on circle k) -- control

The graph-theoretic Dirac operator D_graph acts on the (V+E)-dimensional
combined node+edge space. Its eigenvalues are +/-sqrt(eigenvalues of L_0)
plus zero modes. Its square is block-diagonal:
  D^2_graph = diag(B^T B, B B^T) = diag(L_1, L_0)
"""

import numpy as np
from numpy.linalg import eigh, svd
import json
from collections import Counter


# ============================================================
# Part 1: Build graphs
# ============================================================

def build_dirac_packing_graph(N_max):
    """
    Dirac packing: level k has 2k nodes on a circle (k=1,...,N_max).
    Angular edges: cycle within each level.
    Radial edges: (k, m) -> (k+1, m) and (k+1, m+1) for m=0,...,2k-1.

    Returns: nodes (list of (k,m)), edges (list of ((k1,m1),(k2,m2))), adjacency dict
    """
    nodes = []
    node_index = {}
    idx = 0
    for k in range(1, N_max + 1):
        for m in range(2 * k):
            nodes.append((k, m))
            node_index[(k, m)] = idx
            idx += 1

    V = len(nodes)
    edges = []

    # Angular edges: cycle within each level
    # Each level k has 2k nodes forming a cycle: 0-1-2-...(2k-1)-0
    # A cycle on n nodes has exactly n edges (including the closing edge).
    # We orient each edge m -> (m+1) mod n_k.
    edge_set = set()
    for k in range(1, N_max + 1):
        n_k = 2 * k
        for m in range(n_k):
            m_next = (m + 1) % n_k
            # Canonical form: use the pair with smaller first index
            # to avoid duplicates for n_k=2 case
            n1 = (k, min(m, m_next))
            n2 = (k, max(m, m_next))
            if (n1, n2) not in edge_set:
                edge_set.add((n1, n2))
                # Orient: tail=m, head=m_next
                edges.append(((k, m), (k, m_next)))

    # Radial edges: (k, m) -> (k+1, m) and (k+1, m+1)
    for k in range(1, N_max):
        n_k = 2 * k
        for m in range(n_k):
            edges.append(((k, m), (k + 1, m)))
            edges.append(((k, m), (k + 1, m + 1)))

    return nodes, edges, node_index


def build_scalar_fock_packing_graph(N_max):
    """
    Scalar Fock packing: level k has 2(2k-1) nodes on a circle (k=1,...,N_max).
    Cumulative counts: 2, 8, 18, 32, ... = 2k^2 (matching n^2 Fock degeneracy).
    Angular edges: cycle within each level.
    Radial edges: nearest angular neighbors between adjacent levels.
    """
    nodes = []
    node_index = {}
    idx = 0
    for k in range(1, N_max + 1):
        n_k = 2 * (2 * k - 1)  # 2, 6, 10, 14, ...
        for m in range(n_k):
            nodes.append((k, m))
            node_index[(k, m)] = idx
            idx += 1

    V = len(nodes)
    edges = []

    # Angular edges: cycle within each level
    edge_set = set()
    for k in range(1, N_max + 1):
        n_k = 2 * (2 * k - 1)
        for m in range(n_k):
            m_next = (m + 1) % n_k
            canon = ((k, min(m, m_next)), (k, max(m, m_next)))
            if canon not in edge_set:
                edge_set.add(canon)
                edges.append(((k, m), (k, m_next)))

    # Radial edges: connect by angular proximity
    for k in range(1, N_max):
        n_k = 2 * (2 * k - 1)
        n_k1 = 2 * (2 * (k + 1) - 1)
        # Place nodes at angles
        angles_k = [2 * np.pi * m / n_k for m in range(n_k)]
        angles_k1 = [2 * np.pi * m / n_k1 for m in range(n_k1)]

        # For each node on level k, connect to nearest neighbors on level k+1
        for m in range(n_k):
            theta = angles_k[m]
            # Find the 2 or 3 nearest nodes on level k+1
            dists = []
            for m1 in range(n_k1):
                d = abs(angles_k1[m1] - theta)
                d = min(d, 2 * np.pi - d)
                dists.append((d, m1))
            dists.sort()
            # Connect to closest 2 (or 3 if very close)
            threshold = 2 * np.pi / n_k1 + 1e-10
            for d, m1 in dists:
                if d <= threshold:
                    edges.append(((k, m), (k + 1, m1)))

    return nodes, edges, node_index


# ============================================================
# Part 2: Signed incidence matrix B
# ============================================================

def build_incidence_matrix(nodes, edges, node_index):
    """
    Build V x E signed incidence matrix B.
    Orientation: tail -> head.
    For radial: lower level = tail (+1), upper level = head (-1).
    For angular: lower m = tail (+1), higher m = head (-1).

    B[v, e] = +1 if v is tail of edge e
    B[v, e] = -1 if v is head of edge e
    """
    V = len(nodes)
    E = len(edges)
    B = np.zeros((V, E), dtype=float)

    for e_idx, (n1, n2) in enumerate(edges):
        i1 = node_index[n1]
        i2 = node_index[n2]
        # n1 is tail (+1), n2 is head (-1)
        B[i1, e_idx] = +1.0
        B[i2, e_idx] = -1.0

    return B


def verify_laplacian(B, nodes, edges, node_index):
    """Verify L_0 = B B^T equals D - A (combinatorial graph Laplacian)."""
    V = len(nodes)
    L0 = B @ B.T

    # Build D - A directly
    A = np.zeros((V, V))
    for n1, n2 in edges:
        i1 = node_index[n1]
        i2 = node_index[n2]
        A[i1, i2] = 1.0
        A[i2, i1] = 1.0
    D = np.diag(A.sum(axis=1))
    L_DA = D - A

    diff = np.max(np.abs(L0 - L_DA))
    return diff, L0, A


# ============================================================
# Part 3: D_graph construction and diagonalization
# ============================================================

def build_dirac_operator(B):
    """
    D_graph = | 0    B   |     B is V x E
              | B^T  0   |     B^T is E x V

    Acts on (node_vec, edge_vec) in V+E space:
      D(v, e) = (B @ e, B^T @ v)

    D^2 = diag(B B^T, B^T B) = diag(L_0, L_1)

    Size (V+E) x (V+E)
    """
    V, E = B.shape
    D = np.zeros((V + E, V + E))
    D[:V, V:] = B         # V x E block (top-right)
    D[V:, :V] = B.T       # E x V block (bottom-left)
    return D


# ============================================================
# Part 4: Analysis functions
# ============================================================

def analyze_eigenvalues(evals, label, N_max, V, E):
    """Detailed analysis of D_graph eigenvalues."""
    results = {}

    # Sort by absolute value
    abs_sorted = np.sort(np.abs(evals))
    signed_sorted = np.sort(evals)

    # Identify zero modes (|lambda| < 1e-10)
    zero_mask = np.abs(evals) < 1e-10
    n_zero = np.sum(zero_mask)
    nonzero_evals = evals[~zero_mask]

    # Positive eigenvalues (sorted)
    pos_evals = np.sort(nonzero_evals[nonzero_evals > 0])
    neg_evals = np.sort(nonzero_evals[nonzero_evals < 0])

    # Check +/- pairing
    if len(pos_evals) == len(neg_evals):
        pair_diff = np.max(np.abs(pos_evals + neg_evals[::-1]))
        results['pm_pairing_error'] = float(pair_diff)
    else:
        results['pm_pairing_error'] = None

    # Multiplicities (group eigenvalues within tolerance)
    tol = 1e-8
    unique_pos = []
    multiplicities = []
    for ev in pos_evals:
        found = False
        for i, u in enumerate(unique_pos):
            if abs(ev - u) < tol:
                multiplicities[i] += 1
                found = True
                break
        if not found:
            unique_pos.append(ev)
            multiplicities.append(1)

    results['n_zero_modes'] = int(n_zero)
    results['n_positive'] = len(pos_evals)
    results['unique_positive_eigenvalues'] = [float(x) for x in unique_pos]
    results['multiplicities'] = multiplicities

    # Check linearity: are eigenvalues evenly spaced?
    if len(unique_pos) >= 2:
        spacings = np.diff(unique_pos)
        if len(spacings) >= 2:
            spacing_cv = np.std(spacings) / np.mean(spacings) if np.mean(spacings) > 0 else float('inf')
            results['spacing_cv'] = float(spacing_cv)
            results['mean_spacing'] = float(np.mean(spacings))
        else:
            results['spacing_cv'] = 0.0
            results['mean_spacing'] = float(spacings[0]) if len(spacings) > 0 else None

    # Try scaling to Dirac spectrum |lambda_n| = n + 3/2
    if len(unique_pos) >= 2:
        # Best linear fit: eigenvalue = a * n + b
        n_vals = np.arange(len(unique_pos))
        A_mat = np.column_stack([n_vals, np.ones(len(n_vals))])
        coeffs, residual, _, _ = np.linalg.lstsq(A_mat, unique_pos, rcond=None)
        results['linear_fit_slope'] = float(coeffs[0])
        results['linear_fit_intercept'] = float(coeffs[1])
        predicted = coeffs[0] * n_vals + coeffs[1]
        rel_errors = np.abs(predicted - unique_pos) / np.abs(unique_pos)
        results['linear_fit_max_rel_error'] = float(np.max(rel_errors))

        # Specific check: does kappa * eigenvalue = n + 3/2?
        target_dirac = n_vals + 1.5
        if len(unique_pos) >= 1:
            kappa_fit = np.mean(target_dirac / np.array(unique_pos))
            scaled = kappa_fit * np.array(unique_pos)
            dirac_residuals = np.abs(scaled - target_dirac) / target_dirac
            results['dirac_kappa'] = float(kappa_fit)
            results['dirac_max_rel_error'] = float(np.max(dirac_residuals))

        # Specific check: does kappa * eigenvalue = sqrt(n^2 - 1) for n=1,2,...?
        n_schrod = np.arange(1, len(unique_pos) + 1)
        target_schrod = np.sqrt(n_schrod**2 - 1.0)
        # Skip n=1 (target = 0) for scaling
        if len(unique_pos) >= 2 and target_schrod[0] == 0:
            kappa_s = np.mean(target_schrod[1:] / np.array(unique_pos[1:]))
            scaled_s = kappa_s * np.array(unique_pos)
            # exclude first for residual if target is 0
            schrod_res = np.abs(scaled_s[1:] - target_schrod[1:]) / target_schrod[1:]
            results['schrod_kappa'] = float(kappa_s)
            results['schrod_max_rel_error'] = float(np.max(schrod_res))
        elif len(unique_pos) >= 2:
            kappa_s = np.mean(target_schrod / np.array(unique_pos))
            scaled_s = kappa_s * np.array(unique_pos)
            schrod_res = np.abs(scaled_s - target_schrod) / np.maximum(target_schrod, 1e-15)
            results['schrod_kappa'] = float(kappa_s)
            results['schrod_max_rel_error'] = float(np.max(schrod_res))

    return results, unique_pos, multiplicities


def analyze_eigenvectors(D_mat, evals, evecs, V, E):
    """Examine eigenvector node/edge localization."""
    results = []
    # Find smallest nonzero positive eigenvalue
    pos_mask = evals > 1e-10
    if not np.any(pos_mask):
        return results

    pos_idx = np.where(pos_mask)[0]
    sorted_pos_idx = pos_idx[np.argsort(evals[pos_idx])]

    # Examine first 3 smallest positive eigenvalues
    for i, idx in enumerate(sorted_pos_idx[:min(5, len(sorted_pos_idx))]):
        vec = evecs[:, idx]
        node_part = vec[:V]
        edge_part = vec[V:]

        node_norm_sq = np.sum(node_part**2)
        edge_norm_sq = np.sum(edge_part**2)

        results.append({
            'eigenvalue': float(evals[idx]),
            'node_weight': float(node_norm_sq),
            'edge_weight': float(edge_norm_sq),
            'node_over_edge_ratio': float(node_norm_sq / max(edge_norm_sq, 1e-15)),
        })

    return results


# ============================================================
# Part 5: Weighted and normalized variants
# ============================================================

def build_weighted_incidence(nodes, edges, node_index, weight_fn):
    """Build weighted incidence matrix."""
    V = len(nodes)
    E = len(edges)
    B = np.zeros((V, E), dtype=float)

    for e_idx, (n1, n2) in enumerate(edges):
        i1 = node_index[n1]
        i2 = node_index[n2]
        w = weight_fn(n1, n2)
        B[i1, e_idx] = +w
        B[i2, e_idx] = -w

    return B


# ============================================================
# Main execution
# ============================================================

def run_analysis():
    all_results = {}

    print("=" * 80)
    print("DIRAC OPERATOR ON PACKING GRAPHS")
    print("=" * 80)

    for graph_type in ['dirac', 'scalar']:
        all_results[graph_type] = {}

        print(f"\n{'=' * 70}")
        print(f"GRAPH TYPE: {graph_type.upper()} PACKING")
        print(f"{'=' * 70}")

        if graph_type == 'dirac':
            print("Level k has 2k nodes. Cumulative: 2, 6, 12, 20, 30, ...")
            print("Dirac degeneracy g_n = (n+1)(n+2): 2, 6, 12, 20, 30, ...")
            print("Target spectrum: |lambda_n| = n + 3/2 = 3/2, 5/2, 7/2, ...")
        else:
            print("Level k has 2(2k-1) nodes. Cumulative: 2, 8, 18, 32, 50, ...")
            print("Fock degeneracy n^2: 1, 4, 9, 16, 25, ... (cumulative 2*sum = 2n^2)")
            print("Target spectrum: sqrt(n^2 - 1) = 0, sqrt(3), sqrt(8), ...")

        for N_max in [2, 3, 4, 5]:
            print(f"\n{'-' * 60}")
            print(f"  N_max = {N_max}")
            print(f"{'-' * 60}")

            # Build graph
            if graph_type == 'dirac':
                nodes, edges, node_index = build_dirac_packing_graph(N_max)
            else:
                nodes, edges, node_index = build_scalar_fock_packing_graph(N_max)

            V = len(nodes)
            E = len(edges)
            print(f"  V = {V}, E = {E}, V+E = {V+E}")

            # Level structure
            level_counts = Counter(k for k, m in nodes)
            print(f"  Level counts: {dict(sorted(level_counts.items()))}")

            # Cumulative counts
            cum = 0
            print("  Cumulative state counts: ", end="")
            for k in sorted(level_counts.keys()):
                cum += level_counts[k]
                print(f"N({k})={cum}", end="  ")
            print()

            if graph_type == 'dirac':
                print("  Dirac degeneracies (n+1)(n+2): ", end="")
                for n in range(N_max):
                    print(f"g_{n}={(n+1)*(n+2)}", end="  ")
                print()

            # Build incidence matrix
            B = build_incidence_matrix(nodes, edges, node_index)

            # Verify Laplacian
            diff, L0, A = verify_laplacian(B, nodes, edges, node_index)
            print(f"  L_0 = B B^T verification: max|L0 - (D-A)| = {diff:.2e}")

            # Compute node Laplacian eigenvalues
            L0_evals = np.sort(np.linalg.eigvalsh(L0))
            print(f"  L_0 eigenvalues (first 10): {np.round(L0_evals[:10], 6)}")
            print(f"  L_0 eigenvalues (unique, tol=1e-6):")
            unique_L0 = []
            L0_mults = []
            for ev in L0_evals:
                found = False
                for i, u in enumerate(unique_L0):
                    if abs(ev - u) < 1e-6:
                        L0_mults[i] += 1
                        found = True
                        break
                if not found:
                    unique_L0.append(ev)
                    L0_mults.append(1)
            for u, mult in zip(unique_L0, L0_mults):
                print(f"    lambda = {u:12.6f}  mult = {mult}")

            # Edge Laplacian
            L1 = B.T @ B
            L1_evals = np.sort(np.linalg.eigvalsh(L1))

            # Build D_graph
            D_mat = build_dirac_operator(B)
            evals, evecs = eigh(D_mat)

            print(f"\n  D_graph eigenvalues (all {len(evals)}):")
            # Group by value
            unique_D = []
            D_mults = []
            for ev in sorted(evals):
                found = False
                for i, u in enumerate(unique_D):
                    if abs(ev - u) < 1e-8:
                        D_mults[i] += 1
                        found = True
                        break
                if not found:
                    unique_D.append(ev)
                    D_mults.append(1)
            for u, mult in zip(unique_D, D_mults):
                print(f"    lambda = {u:12.8f}  mult = {mult}")

            # SVD verification
            U_svd, sigma, Vt_svd = svd(B, full_matrices=False)
            sigma_sorted = np.sort(sigma)[::-1]
            pos_D = np.sort(evals[evals > 1e-10])
            print(f"\n  SVD verification:")
            print(f"    Singular values of B: {np.round(sigma_sorted, 8)}")
            if len(pos_D) > 0:
                print(f"    Positive D_graph eigenvalues: {np.round(pos_D, 8)}")
                # Match singular values to positive eigenvalues
                if len(sigma_sorted[sigma_sorted > 1e-10]) == len(pos_D):
                    svd_match = np.max(np.abs(np.sort(sigma_sorted[sigma_sorted > 1e-10]) - np.sort(pos_D)))
                    print(f"    Max |sigma - lambda_+| = {svd_match:.2e}")
                else:
                    print(f"    Count mismatch: {len(sigma_sorted[sigma_sorted > 1e-10])} SVs vs {len(pos_D)} pos evals")

            # D^2 verification
            D2_evals = evals**2
            D2_sorted = np.sort(D2_evals)
            # Should be union of L0 and L1 eigenvalues
            L_union = np.sort(np.concatenate([L0_evals, L1_evals]))
            if len(D2_sorted) == len(L_union):
                d2_check = np.max(np.abs(D2_sorted - L_union))
                print(f"    D^2 vs L0 U L1 check: max error = {d2_check:.2e}")

            # Betti numbers
            beta_0 = np.sum(np.abs(L0_evals) < 1e-10)
            beta_1 = E - V + beta_0
            n_zero = np.sum(np.abs(evals) < 1e-10)
            print(f"\n  Topology:")
            print(f"    beta_0 (components) = {beta_0}")
            print(f"    beta_1 (cycles) = {beta_1}")
            print(f"    dim(ker D_graph) = {n_zero} (should be beta_0 + beta_1 = {beta_0 + beta_1})")

            # Detailed eigenvalue analysis
            analysis, unique_pos, mults = analyze_eigenvalues(
                evals, graph_type, N_max, V, E
            )

            print(f"\n  Eigenvalue analysis:")
            print(f"    +/- pairing error: {analysis.get('pm_pairing_error', 'N/A')}")
            print(f"    Unique positive eigenvalues: {len(unique_pos)}")

            if 'spacing_cv' in analysis:
                print(f"    Spacing CV: {analysis['spacing_cv']:.6f} (0 = exactly linear)")

            if 'linear_fit_slope' in analysis:
                print(f"    Linear fit: lambda = {analysis['linear_fit_slope']:.6f} * n + {analysis['linear_fit_intercept']:.6f}")
                print(f"    Linear fit max rel error: {analysis['linear_fit_max_rel_error']:.6f}")

            if 'dirac_kappa' in analysis:
                print(f"\n    Dirac comparison (target: n + 3/2):")
                print(f"      kappa_D = {analysis['dirac_kappa']:.8f}")
                print(f"      Max rel error: {analysis['dirac_max_rel_error']:.6f}")
                # Print the scaled values vs targets
                n_vals = np.arange(len(unique_pos))
                target_dirac = n_vals + 1.5
                scaled = analysis['dirac_kappa'] * np.array(unique_pos)
                print(f"      {'n':>4s} {'eigenvalue':>12s} {'scaled':>12s} {'target':>12s} {'rel_err':>10s}")
                for n, ev, sc, tg in zip(n_vals, unique_pos, scaled, target_dirac):
                    re = abs(sc - tg) / tg
                    print(f"      {n:4d} {ev:12.6f} {sc:12.6f} {tg:12.6f} {re:10.6f}")

            if 'schrod_kappa' in analysis:
                print(f"\n    Schrodinger comparison (target: sqrt(n^2 - 1)):")
                print(f"      kappa_S = {analysis['schrod_kappa']:.8f}")
                print(f"      Max rel error: {analysis['schrod_max_rel_error']:.6f}")

            # Degeneracy comparison
            print(f"\n    Degeneracy structure (positive eigenvalues):")
            print(f"      {'level':>5s} {'eigenvalue':>12s} {'mult':>6s}", end="")
            if graph_type == 'dirac':
                print(f" {'(n+1)(n+2)':>12s}")
            else:
                print(f" {'n^2':>12s}")
            for i, (ev, mult) in enumerate(zip(unique_pos, mults)):
                if graph_type == 'dirac':
                    target_deg = (i + 1) * (i + 2)
                    print(f"      {i:5d} {ev:12.6f} {mult:6d} {target_deg:12d}")
                else:
                    target_deg = (i + 1)**2
                    print(f"      {i:5d} {ev:12.6f} {mult:6d} {target_deg:12d}")

            # Eigenvector analysis
            evec_analysis = analyze_eigenvectors(D_mat, evals, evecs, V, E)
            if evec_analysis:
                print(f"\n    Eigenvector structure (smallest positive eigenvalues):")
                print(f"      {'eigenvalue':>12s} {'node_wt':>10s} {'edge_wt':>10s} {'node/edge':>10s}")
                for ea in evec_analysis:
                    print(f"      {ea['eigenvalue']:12.6f} {ea['node_weight']:10.6f} {ea['edge_weight']:10.6f} {ea['node_over_edge_ratio']:10.4f}")

            # Store results
            all_results[graph_type][N_max] = {
                'V': V, 'E': E,
                'beta_0': int(beta_0),
                'beta_1': int(beta_1),
                'n_zero_modes': int(n_zero),
                'unique_positive_eigenvalues': [float(x) for x in unique_pos],
                'positive_multiplicities': mults,
                'L0_eigenvalues': [float(x) for x in unique_L0],
                'L0_multiplicities': L0_mults,
                'analysis': analysis,
                'eigenvector_analysis': evec_analysis,
            }

    # ============================================================
    # Part 5: Weighted variant (Dirac packing only)
    # ============================================================

    print(f"\n{'=' * 70}")
    print("WEIGHTED INCIDENCE MATRIX (Dirac packing, weight = 1/k for radial)")
    print(f"{'=' * 70}")

    all_results['dirac_weighted'] = {}

    for N_max in [3, 4, 5]:
        nodes, edges, node_index = build_dirac_packing_graph(N_max)
        V = len(nodes)
        E = len(edges)

        def weight_fn(n1, n2):
            k1, m1 = n1
            k2, m2 = n2
            if k1 == k2:
                return 1.0  # angular
            else:
                return 1.0 / max(k1, k2)  # radial: weight by 1/radius

        Bw = build_weighted_incidence(nodes, edges, node_index, weight_fn)
        Dw = build_dirac_operator(Bw)
        evals_w, evecs_w = eigh(Dw)

        print(f"\n  N_max = {N_max}, V={V}, E={E}")

        # Group eigenvalues
        unique_Dw = []
        Dw_mults = []
        for ev in sorted(evals_w):
            found = False
            for i, u in enumerate(unique_Dw):
                if abs(ev - u) < 1e-8:
                    Dw_mults[i] += 1
                    found = True
                    break
            if not found:
                unique_Dw.append(ev)
                Dw_mults.append(1)

        pos_w = [(u, m) for u, m in zip(unique_Dw, Dw_mults) if u > 1e-10]
        print(f"  Positive eigenvalues (weighted):")
        for u, m in pos_w:
            print(f"    lambda = {u:12.8f}  mult = {m}")

        # Compare to Dirac target
        if pos_w:
            unique_pos_w = [u for u, m in pos_w]
            n_vals = np.arange(len(unique_pos_w))
            target_dirac = n_vals + 1.5
            kappa_w = np.mean(target_dirac / np.array(unique_pos_w))
            scaled_w = kappa_w * np.array(unique_pos_w)
            print(f"  Dirac fit: kappa = {kappa_w:.8f}")
            print(f"    {'n':>4s} {'eigenvalue':>12s} {'scaled':>12s} {'target':>12s} {'rel_err':>10s}")
            for n, ev, sc, tg in zip(n_vals, unique_pos_w, scaled_w, target_dirac):
                re = abs(sc - tg) / tg
                print(f"    {n:4d} {ev:12.6f} {sc:12.6f} {tg:12.6f} {re:10.6f}")

        all_results['dirac_weighted'][N_max] = {
            'V': V, 'E': E,
            'unique_positive_eigenvalues': [float(u) for u, m in pos_w],
            'positive_multiplicities': [m for u, m in pos_w],
        }

    # ============================================================
    # Part 5b: Alternative weight: sqrt(k) for radial edges
    # ============================================================

    print(f"\n{'=' * 70}")
    print("WEIGHTED INCIDENCE: sqrt(k) radial weight (Dirac packing)")
    print(f"{'=' * 70}")

    all_results['dirac_weighted_sqrt'] = {}

    for N_max in [3, 4, 5]:
        nodes, edges, node_index = build_dirac_packing_graph(N_max)
        V = len(nodes)
        E = len(edges)

        def weight_fn_sqrt(n1, n2):
            k1, m1 = n1
            k2, m2 = n2
            if k1 == k2:
                return 1.0
            else:
                return np.sqrt(max(k1, k2))

        Bw2 = build_weighted_incidence(nodes, edges, node_index, weight_fn_sqrt)
        Dw2 = build_dirac_operator(Bw2)
        evals_w2, _ = eigh(Dw2)

        unique_Dw2 = []
        Dw2_mults = []
        for ev in sorted(evals_w2):
            found = False
            for i, u in enumerate(unique_Dw2):
                if abs(ev - u) < 1e-8:
                    Dw2_mults[i] += 1
                    found = True
                    break
            if not found:
                unique_Dw2.append(ev)
                Dw2_mults.append(1)

        pos_w2 = [(u, m) for u, m in zip(unique_Dw2, Dw2_mults) if u > 1e-10]

        print(f"\n  N_max = {N_max}")
        for u, m in pos_w2:
            print(f"    lambda = {u:12.8f}  mult = {m}")

        all_results['dirac_weighted_sqrt'][N_max] = {
            'unique_positive_eigenvalues': [float(u) for u, m in pos_w2],
            'positive_multiplicities': [m for u, m in pos_w2],
        }

    # ============================================================
    # Part 6: Summary table and verdict
    # ============================================================

    print(f"\n{'=' * 70}")
    print("SUMMARY TABLE")
    print(f"{'=' * 70}")

    print(f"\n{'Graph':>10s} {'N_max':>5s} {'V':>5s} {'E':>5s} {'dim(D)':>6s} "
          f"{'#zero':>6s} {'b0':>4s} {'b1':>4s} {'#pos_uniq':>9s} "
          f"{'spacing_CV':>10s} {'linear?':>8s}")
    print("-" * 90)

    for gtype in ['dirac', 'scalar']:
        for N_max in sorted(all_results[gtype].keys()):
            r = all_results[gtype][N_max]
            a = r['analysis']
            cv = a.get('spacing_cv', None)
            cv_str = f"{cv:.4f}" if cv is not None else "N/A"
            linear = "YES" if (cv is not None and cv < 0.05) else "NO" if cv is not None else "?"
            print(f"{gtype:>10s} {N_max:>5d} {r['V']:>5d} {r['E']:>5d} {r['V']+r['E']:>6d} "
                  f"{r['n_zero_modes']:>6d} {r['beta_0']:>4d} {r['beta_1']:>4d} "
                  f"{len(r['unique_positive_eigenvalues']):>9d} "
                  f"{cv_str:>10s} {linear:>8s}")

    # ============================================================
    # Additional analysis: eigenvalue ratios
    # ============================================================

    print(f"\n{'=' * 70}")
    print("EIGENVALUE RATIO ANALYSIS")
    print(f"{'=' * 70}")

    for gtype in ['dirac', 'scalar']:
        print(f"\n  {gtype.upper()} packing:")
        for N_max in sorted(all_results[gtype].keys()):
            r = all_results[gtype][N_max]
            up = r['unique_positive_eigenvalues']
            if len(up) >= 2:
                ratios = [up[i+1]/up[i] for i in range(len(up)-1)]
                print(f"    N_max={N_max}: ratios lambda_{i+1}/lambda_i = {[f'{r:.4f}' for r in ratios]}")

                # Check if eigenvalues go as n + offset (Dirac-like)
                if len(up) >= 3:
                    # Fit lambda = a*n + b
                    n_arr = np.arange(len(up), dtype=float)
                    A_fit = np.column_stack([n_arr, np.ones(len(n_arr))])
                    coeffs = np.linalg.lstsq(A_fit, up, rcond=None)[0]
                    print(f"           Linear fit: {coeffs[0]:.6f} * n + {coeffs[1]:.6f}")
                    residuals = np.array(up) - (coeffs[0] * n_arr + coeffs[1])
                    print(f"           Residuals: {[f'{r:.6f}' for r in residuals]}")

                # Check if eigenvalues go as sqrt(something)
                sq_vals = [u**2 for u in up]
                sq_diffs = [sq_vals[i+1] - sq_vals[i] for i in range(len(sq_vals)-1)]
                print(f"           lambda^2 values: {[f'{s:.4f}' for s in sq_vals]}")
                print(f"           lambda^2 diffs: {[f'{d:.4f}' for d in sq_diffs]}")

    # ============================================================
    # Key comparison: Dirac packing eigenvalue^2 vs L0 eigenvalues
    # ============================================================

    print(f"\n{'=' * 70}")
    print("D_graph^2 vs LAPLACIAN EIGENVALUES")
    print(f"{'=' * 70}")

    for gtype in ['dirac', 'scalar']:
        print(f"\n  {gtype.upper()} packing:")
        for N_max in sorted(all_results[gtype].keys()):
            r = all_results[gtype][N_max]
            print(f"    N_max={N_max}:")
            print(f"      L_0 eigenvalues: {[f'{e:.4f}' for e in r['L0_eigenvalues']]}")
            print(f"      L_0 multiplicities: {r['L0_multiplicities']}")

            # The sqrt of nonzero L0 eigenvalues should appear in D_graph
            L0_nonzero = [e for e in r['L0_eigenvalues'] if e > 1e-10]
            sqrt_L0 = [np.sqrt(e) for e in L0_nonzero]
            print(f"      sqrt(L_0 nonzero): {[f'{s:.6f}' for s in sqrt_L0]}")
            print(f"      D_graph positive: {[f'{e:.6f}' for e in r['unique_positive_eigenvalues']]}")

    # ============================================================
    # Verdict
    # ============================================================

    print(f"\n{'=' * 70}")
    print("VERDICT")
    print(f"{'=' * 70}")

    print("""
    1. D_graph eigenvalues are +/-sqrt(L_0 eigenvalues) U {zero modes}.
       This is a mathematical identity, not a surprise.

    2. The question reduces to: does L_0 (the combinatorial graph Laplacian
       of the packing graph) have eigenvalues that go as (n + 3/2)^2 for
       the Dirac packing?

    3. Since L_0 = D - A is the standard graph Laplacian, and we already know
       from Paper 7 that the scalar Fock graph's Laplacian gives n^2 - 1
       (with kappa = -1/16), the question is whether the Dirac packing
       graph's Laplacian gives (n + 3/2)^2 - c for some constant c.
    """)

    # Direct comparison for Dirac packing
    print("  DIRAC PACKING L_0 eigenvalue analysis:")
    for N_max in sorted(all_results['dirac'].keys()):
        r = all_results['dirac'][N_max]
        L0_nonzero = [e for e in r['L0_eigenvalues'] if e > 1e-10]
        if L0_nonzero:
            # Check if L0 eigenvalues go as (n + 3/2)^2
            n_vals = np.arange(len(L0_nonzero))
            target = [(n + 1.5)**2 for n in n_vals]

            if len(L0_nonzero) >= 2:
                kappa_L = target[1] / L0_nonzero[1] if L0_nonzero[1] > 0 else None
                if kappa_L:
                    scaled = [kappa_L * e for e in L0_nonzero]
                    print(f"    N_max={N_max}: L0 nonzero = {[f'{e:.4f}' for e in L0_nonzero]}")
                    print(f"              Target (n+3/2)^2 = {[f'{t:.4f}' for t in target[:len(L0_nonzero)]]}")
                    print(f"              kappa_L = {kappa_L:.6f}")
                    print(f"              kappa_L * L0 = {[f'{s:.4f}' for s in scaled]}")

            # Also check: do L0 eigenvalues go as n*(n+3)?
            # Because (n+3/2)^2 = n^2 + 3n + 9/4, and (n+3/2)^2 - 9/4 = n(n+3)
            target2 = [n * (n + 3) for n in n_vals]
            print(f"              Target n(n+3) = {target2[:len(L0_nonzero)]}")

    # Check: does the raw L0 on the Dirac graph match n(n+d) for any d?
    print("\n  DIRAC PACKING: probing L_0 eigenvalue structure")
    for N_max in [4, 5]:
        r = all_results['dirac'][N_max]
        L0_all = r['L0_eigenvalues']
        L0_mults = r['L0_multiplicities']
        print(f"    N_max={N_max}:")
        for ev, mult in zip(L0_all, L0_mults):
            if ev > 1e-10:
                # Solve n(n+d) = ev for d, assuming n = integer
                # Try matching with sequential indices
                pass

        # Just print the raw nonzero eigenvalues with their square roots
        nonzero = [(ev, mult) for ev, mult in zip(L0_all, L0_mults) if ev > 1e-10]
        print(f"    Nonzero L_0: ", end="")
        for ev, mult in nonzero:
            print(f"({ev:.4f}, mult={mult})", end="  ")
        print()
        print(f"    sqrt(L_0):   ", end="")
        for ev, mult in nonzero:
            print(f"({np.sqrt(ev):.4f}, mult={mult})", end="  ")
        print()

    # Save results
    output_path = 'C:/Users/jlout/Desktop/Project_Geometric/debug/data/dirac_operator_graph_test.json'

    # Convert numpy types for JSON serialization
    def convert(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert(v) for v in obj]
        return obj

    with open(output_path, 'w') as f:
        json.dump(convert(all_results), f, indent=2)

    print(f"\n  Results saved to {output_path}")

    return all_results


if __name__ == '__main__':
    results = run_analysis()
