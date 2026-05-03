"""
Test whether photon angular momentum quantum numbers commute with L₁.
=====================================================================

The Fock graph has nodes (n, l, m) and edges of two types:
  - Angular (L±): (n, l, m) <-> (n, l, m±1)  -- carries Δm = ±1
  - Radial  (T±): (n, l, m) <-> (n±1, l, m)  -- carries Δm = 0

No edge changes l (the graph is disconnected by l-sector).

We define edge-space operators:
  M_edge[e] = m_target - m_source  (signed, on directed edges)
  |Δm|[e] = |m_target - m_source|  (unsigned)
  Δn[e] = |n_target - n_source|    (radial vs angular classifier)

Question: does L₁ = B^T B block-diagonalize by these quantum numbers?
If [L₁, diag(|Δm|)] = 0, then L₁ eigenmodes have definite |Δm| and
the photon naturally carries angular momentum quantum numbers.

Since L₁ acts on UNDIRECTED edges (each edge appears once), we use
|Δm| as the natural quantum number (the sign depends on orientation
convention, but the magnitude is intrinsic to the edge).

Author: Claude Code
Date: 2026-05-01
"""

import json
import sys
from pathlib import Path

import numpy as np

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from geovac.lattice import GeometricLattice


def build_fock_graph_data(n_max: int):
    """Build graph topology and return nodes, edges, incidence, L1."""
    lat = GeometricLattice(max_n=n_max, topological_weights=False)
    states = lat.states  # list of (n, l, m)
    V = len(states)

    # Build state index
    state_index = {s: i for i, s in enumerate(states)}

    # Extract undirected edges (i < j) from adjacency
    adj = lat.adjacency
    rows, cols = adj.nonzero()
    edge_set = set()
    for r, c in zip(rows, cols):
        if r < c:
            edge_set.add((int(r), int(c)))
    edges = sorted(edge_set)
    E = len(edges)

    # Build signed incidence matrix B (V x E)
    B = np.zeros((V, E), dtype=np.float64)
    for k, (i, j) in enumerate(edges):
        B[i, k] = 1.0
        B[j, k] = -1.0

    # Edge Laplacian L1 = B^T B
    L1 = B.T @ B

    # Node Laplacian L0 = B B^T (for reference)
    L0 = B @ B.T

    return {
        'states': states,
        'V': V,
        'edges': edges,
        'E': E,
        'B': B,
        'L1': L1,
        'L0': L0,
    }


def classify_edges(states, edges):
    """Classify each edge by its quantum number transfer.

    Returns arrays indexed by edge index:
      delta_m[e] = m_j - m_i  (signed, based on canonical i<j ordering)
      abs_delta_m[e] = |m_j - m_i|
      delta_n[e] = |n_j - n_i|  (0 for angular, 1 for radial)
      delta_l[e] = |l_j - l_i|  (should be 0 for all edges on Fock graph)
      edge_type[e] = 'angular' or 'radial'
    """
    E = len(edges)
    delta_m = np.zeros(E, dtype=int)
    abs_delta_m = np.zeros(E, dtype=int)
    delta_n = np.zeros(E, dtype=int)
    delta_l = np.zeros(E, dtype=int)
    edge_types = []

    for k, (i, j) in enumerate(edges):
        n_i, l_i, m_i = states[i]
        n_j, l_j, m_j = states[j]

        delta_m[k] = m_j - m_i
        abs_delta_m[k] = abs(m_j - m_i)
        delta_n[k] = abs(n_j - n_i)
        delta_l[k] = abs(l_j - l_i)

        if n_i == n_j and l_i == l_j:
            edge_types.append('angular')
        elif l_i == l_j and m_i == m_j:
            edge_types.append('radial')
        else:
            edge_types.append('unknown')

    return {
        'delta_m': delta_m,
        'abs_delta_m': abs_delta_m,
        'delta_n': delta_n,
        'delta_l': delta_l,
        'edge_types': edge_types,
    }


def check_commutator(L1, quantum_number_array, name):
    """Check if L1 commutes with diag(quantum_number_array).

    Returns:
      - Frobenius norm of [L1, M]
      - Relative norm: ||[L1, M]|| / (||L1|| * ||M||)
      - Block diagonal fraction: what fraction of L1's nonzero entries
        connect edges with the SAME quantum number value
    """
    E = L1.shape[0]
    M = np.diag(quantum_number_array.astype(float))

    commutator = L1 @ M - M @ L1
    frob_comm = np.linalg.norm(commutator, 'fro')
    frob_L1 = np.linalg.norm(L1, 'fro')
    frob_M = np.linalg.norm(M, 'fro')

    if frob_L1 * frob_M > 0:
        rel_norm = frob_comm / (frob_L1 * frob_M)
    else:
        rel_norm = 0.0

    # Block diagonal analysis: for each nonzero L1[i,j] (i != j),
    # check if quantum_number[i] == quantum_number[j]
    n_off_diag_nonzero = 0
    n_block_diagonal = 0
    n_off_block = 0

    for i in range(E):
        for j in range(E):
            if i == j:
                continue
            if abs(L1[i, j]) > 1e-12:
                n_off_diag_nonzero += 1
                if quantum_number_array[i] == quantum_number_array[j]:
                    n_block_diagonal += 1
                else:
                    n_off_block += 1

    if n_off_diag_nonzero > 0:
        block_diag_fraction = n_block_diagonal / n_off_diag_nonzero
    else:
        block_diag_fraction = 1.0

    return {
        'name': name,
        'frobenius_norm_commutator': frob_comm,
        'frobenius_norm_L1': frob_L1,
        'frobenius_norm_M': frob_M,
        'relative_norm': rel_norm,
        'n_off_diag_nonzero_L1': n_off_diag_nonzero,
        'n_block_diagonal': n_block_diagonal,
        'n_off_block': n_off_block,
        'block_diagonal_fraction': block_diag_fraction,
        'commutes_exactly': frob_comm < 1e-10,
    }


def analyze_L1_block_structure(L1, quantum_number_array, name, states, edges):
    """Detailed block structure analysis.

    Groups edges by their quantum number value and checks whether L1
    is block-diagonal in that grouping.
    """
    E = L1.shape[0]
    unique_vals = sorted(set(quantum_number_array.tolist()))

    blocks = {}
    for val in unique_vals:
        mask = quantum_number_array == val
        indices = np.where(mask)[0]
        blocks[val] = indices

    # For each pair of distinct blocks, compute the cross-block norm
    cross_block_info = []
    total_cross_block_norm_sq = 0.0
    total_within_block_norm_sq = 0.0

    for val in unique_vals:
        idx = blocks[val]
        sub_block = L1[np.ix_(idx, idx)]
        total_within_block_norm_sq += np.linalg.norm(sub_block, 'fro')**2

    for i_val, val1 in enumerate(unique_vals):
        for val2 in unique_vals[i_val+1:]:
            idx1 = blocks[val1]
            idx2 = blocks[val2]
            cross = L1[np.ix_(idx1, idx2)]
            cross_norm = np.linalg.norm(cross, 'fro')
            total_cross_block_norm_sq += 2 * cross_norm**2  # factor 2 for symmetry
            if cross_norm > 1e-12:
                cross_block_info.append({
                    'block1_value': int(val1),
                    'block2_value': int(val2),
                    'cross_block_frobenius': cross_norm,
                    'n_edges_block1': len(idx1),
                    'n_edges_block2': len(idx2),
                })

    total_norm_sq = np.linalg.norm(L1, 'fro')**2

    # Diagonal contribution
    diag_norm_sq = np.sum(np.diag(L1)**2)
    off_diag_norm_sq = total_norm_sq - diag_norm_sq

    return {
        'name': name,
        'unique_values': [int(v) for v in unique_vals],
        'block_sizes': {int(v): len(blocks[v]) for v in unique_vals},
        'total_frobenius_sq': total_norm_sq,
        'within_block_frobenius_sq': total_within_block_norm_sq,
        'cross_block_frobenius_sq': total_cross_block_norm_sq,
        'fraction_within_block': total_within_block_norm_sq / total_norm_sq if total_norm_sq > 0 else 1.0,
        'fraction_cross_block': total_cross_block_norm_sq / total_norm_sq if total_norm_sq > 0 else 0.0,
        'is_block_diagonal': total_cross_block_norm_sq < 1e-10,
        'nonzero_cross_blocks': cross_block_info,
    }


def eigenvalue_decomposition_by_block(L1, quantum_number_array, name):
    """Compute L1 eigenvalues and check if eigenvectors are localized in blocks."""
    E = L1.shape[0]
    if E == 0:
        return {'name': name, 'note': 'empty edge space'}

    eigenvalues, eigenvectors = np.linalg.eigh(L1)

    unique_vals = sorted(set(quantum_number_array.tolist()))
    blocks = {val: np.where(quantum_number_array == val)[0] for val in unique_vals}

    # For each eigenvector, compute the fraction of its weight in each block
    eigvec_info = []
    for k in range(E):
        ev = eigenvalues[k]
        vec = eigenvectors[:, k]
        weight_by_block = {}
        total_weight = np.sum(vec**2)
        for val in unique_vals:
            idx = blocks[val]
            block_weight = np.sum(vec[idx]**2)
            weight_by_block[int(val)] = block_weight / total_weight if total_weight > 0 else 0.0

        # Dominant block
        dominant_val = max(weight_by_block, key=weight_by_block.get)
        dominant_fraction = weight_by_block[dominant_val]

        eigvec_info.append({
            'eigenvalue': float(ev),
            'dominant_block': dominant_val,
            'dominant_fraction': dominant_fraction,
            'is_pure': dominant_fraction > 0.999,
        })

    n_pure = sum(1 for e in eigvec_info if e['is_pure'])

    return {
        'name': name,
        'n_eigenvectors': E,
        'n_pure_block_eigenvectors': n_pure,
        'fraction_pure': n_pure / E if E > 0 else 1.0,
        'eigenvector_details': eigvec_info,
    }


def run_analysis(n_max: int):
    """Run the full commutator analysis at a given n_max."""
    print(f"\n{'='*70}")
    print(f"  PHOTON M_L COMMUTATOR ANALYSIS AT n_max = {n_max}")
    print(f"{'='*70}")

    # Build graph
    data = build_fock_graph_data(n_max)
    states = data['states']
    edges = data['edges']
    L1 = data['L1']
    V, E = data['V'], data['E']

    print(f"\n  Graph: V = {V} nodes, E = {E} edges")

    # Classify edges
    edge_class = classify_edges(states, edges)

    # Edge type census
    n_angular = edge_class['edge_types'].count('angular')
    n_radial = edge_class['edge_types'].count('radial')
    n_unknown = edge_class['edge_types'].count('unknown')
    print(f"  Edge types: {n_angular} angular (dm=+-1), {n_radial} radial (dm=0)")
    if n_unknown > 0:
        print(f"  WARNING: {n_unknown} edges of unknown type!")

    # Verify Δl = 0 for all edges
    assert np.all(edge_class['delta_l'] == 0), "Expected Δl=0 for all edges!"
    print(f"  Δl = 0 confirmed for ALL edges (graph decomposes by l)")

    # Values of |Δm|
    unique_abs_dm = sorted(set(edge_class['abs_delta_m'].tolist()))
    print(f"  |Δm| values on edges: {unique_abs_dm}")
    for val in unique_abs_dm:
        count = np.sum(edge_class['abs_delta_m'] == val)
        print(f"    |Δm| = {val}: {count} edges")

    results = {'n_max': n_max, 'V': V, 'E': E}
    results['n_angular_edges'] = n_angular
    results['n_radial_edges'] = n_radial

    # =====================================================================
    # Test 1: Does L1 commute with diag(|Δm|)?
    # =====================================================================
    print(f"\n  --- Test 1: [L₁, diag(|Δm|)] ---")
    comm_abs_dm = check_commutator(L1, edge_class['abs_delta_m'], '|delta_m|')
    print(f"    ||[L₁, |Δm|]||_F = {comm_abs_dm['frobenius_norm_commutator']:.6e}")
    print(f"    Relative norm    = {comm_abs_dm['relative_norm']:.6e}")
    print(f"    Block-diagonal fraction of L₁ entries: {comm_abs_dm['block_diagonal_fraction']:.4f}")
    print(f"    Commutes exactly: {comm_abs_dm['commutes_exactly']}")
    results['commutator_abs_delta_m'] = comm_abs_dm

    # =====================================================================
    # Test 2: Block structure analysis for |Δm|
    # =====================================================================
    print(f"\n  --- Test 2: Block structure of L₁ by |Δm| ---")
    block_abs_dm = analyze_L1_block_structure(
        L1, edge_class['abs_delta_m'], '|delta_m|', states, edges
    )
    print(f"    Block sizes: {block_abs_dm['block_sizes']}")
    print(f"    Fraction of ||L₁||²_F within blocks: {block_abs_dm['fraction_within_block']:.6f}")
    print(f"    Fraction of ||L₁||²_F cross-block:   {block_abs_dm['fraction_cross_block']:.6f}")
    print(f"    Is block-diagonal: {block_abs_dm['is_block_diagonal']}")
    if block_abs_dm['nonzero_cross_blocks']:
        print(f"    Cross-block couplings:")
        for cb in block_abs_dm['nonzero_cross_blocks']:
            print(f"      |Δm|={cb['block1_value']} <-> |Δm|={cb['block2_value']}: "
                  f"||cross||_F = {cb['cross_block_frobenius']:.4f}")
    results['block_structure_abs_delta_m'] = block_abs_dm

    # =====================================================================
    # Test 3: Does L1 commute with diag(Δn)?
    # (Δn classifies edges as angular vs radial)
    # =====================================================================
    print(f"\n  --- Test 3: [L₁, diag(Δn)] (angular vs radial edge classification) ---")
    comm_dn = check_commutator(L1, edge_class['delta_n'], 'delta_n')
    print(f"    ||[L₁, Δn]||_F = {comm_dn['frobenius_norm_commutator']:.6e}")
    print(f"    Relative norm   = {comm_dn['relative_norm']:.6e}")
    print(f"    Block-diagonal fraction: {comm_dn['block_diagonal_fraction']:.4f}")
    print(f"    Commutes exactly: {comm_dn['commutes_exactly']}")
    results['commutator_delta_n'] = comm_dn

    # =====================================================================
    # Test 4: Block structure for Δn (angular vs radial)
    # =====================================================================
    print(f"\n  --- Test 4: Block structure of L₁ by Δn ---")
    block_dn = analyze_L1_block_structure(
        L1, edge_class['delta_n'], 'delta_n', states, edges
    )
    print(f"    Block sizes: {block_dn['block_sizes']}")
    print(f"    Fraction within blocks: {block_dn['fraction_within_block']:.6f}")
    print(f"    Fraction cross-block:   {block_dn['fraction_cross_block']:.6f}")
    print(f"    Is block-diagonal: {block_dn['is_block_diagonal']}")
    results['block_structure_delta_n'] = block_dn

    # =====================================================================
    # Test 5: Signed Δm -- does L1 commute with diag(Δm_signed)?
    # For undirected edges with canonical orientation i<j, this is m_j - m_i
    # =====================================================================
    print(f"\n  --- Test 5: [L₁, diag(Δm_signed)] ---")
    comm_signed_dm = check_commutator(L1, edge_class['delta_m'], 'delta_m_signed')
    print(f"    ||[L₁, Δm_signed]||_F = {comm_signed_dm['frobenius_norm_commutator']:.6e}")
    print(f"    Relative norm          = {comm_signed_dm['relative_norm']:.6e}")
    print(f"    Block-diagonal fraction: {comm_signed_dm['block_diagonal_fraction']:.4f}")
    print(f"    Commutes exactly: {comm_signed_dm['commutes_exactly']}")
    results['commutator_delta_m_signed'] = comm_signed_dm

    # =====================================================================
    # Test 6: Combined quantum number (|Δm|, within_l_sector)
    # Since l is conserved, we can also label edges by the l-sector they're in
    # =====================================================================
    print(f"\n  --- Test 6: L₁ block structure by l-sector ---")
    edge_l_sector = np.array([states[edges[k][0]][1] for k in range(E)], dtype=int)
    block_l = analyze_L1_block_structure(L1, edge_l_sector, 'l_sector', states, edges)
    print(f"    Block sizes: {block_l['block_sizes']}")
    print(f"    Fraction within blocks: {block_l['fraction_within_block']:.6f}")
    print(f"    Fraction cross-block:   {block_l['fraction_cross_block']:.6f}")
    print(f"    Is block-diagonal: {block_l['is_block_diagonal']}")
    results['block_structure_l_sector'] = block_l

    # =====================================================================
    # Test 7: Eigenvector localization
    # =====================================================================
    print(f"\n  --- Test 7: L₁ eigenvector localization by |Δm| ---")
    eigvec_abs_dm = eigenvalue_decomposition_by_block(
        L1, edge_class['abs_delta_m'], '|delta_m|'
    )
    print(f"    Total eigenvectors: {eigvec_abs_dm['n_eigenvectors']}")
    print(f"    Pure-block (>99.9%): {eigvec_abs_dm['n_pure_block_eigenvectors']}")
    print(f"    Fraction pure: {eigvec_abs_dm['fraction_pure']:.4f}")
    results['eigenvector_localization_abs_dm'] = {
        'n_eigenvectors': eigvec_abs_dm['n_eigenvectors'],
        'n_pure': eigvec_abs_dm['n_pure_block_eigenvectors'],
        'fraction_pure': eigvec_abs_dm['fraction_pure'],
    }

    # Show a few mixed eigenvectors if any
    mixed = [e for e in eigvec_abs_dm['eigenvector_details'] if not e['is_pure']]
    if mixed:
        print(f"\n    Mixed eigenvectors (not >99.9% in one block):")
        for m in mixed[:10]:
            print(f"      λ = {m['eigenvalue']:.4f}, dominant |Δm|={m['dominant_block']}, "
                  f"fraction={m['dominant_fraction']:.4f}")

    # =====================================================================
    # Test 8: Eigenvector localization by l-sector
    # =====================================================================
    print(f"\n  --- Test 8: L₁ eigenvector localization by l-sector ---")
    eigvec_l = eigenvalue_decomposition_by_block(L1, edge_l_sector, 'l_sector')
    print(f"    Total eigenvectors: {eigvec_l['n_eigenvectors']}")
    print(f"    Pure-block (>99.9%): {eigvec_l['n_pure_block_eigenvectors']}")
    print(f"    Fraction pure: {eigvec_l['fraction_pure']:.4f}")
    results['eigenvector_localization_l_sector'] = {
        'n_eigenvectors': eigvec_l['n_eigenvectors'],
        'n_pure': eigvec_l['n_pure_block_eigenvectors'],
        'fraction_pure': eigvec_l['fraction_pure'],
    }

    return results


def within_l_sector_analysis(n_max: int):
    """Analyze the within-l-sector structure at n_max=3 to understand mixing.

    Within an l-sector, edges are either:
      - Angular: (n, l, m) <-> (n, l, m+1)  -- carries |dm|=1
      - Radial:  (n, l, m) <-> (n+1, l, m)  -- carries |dm|=0

    L1 couples edges that SHARE A NODE. An angular edge and a radial edge
    share a node iff they meet at the same (n,l,m) vertex. This is the
    mechanism for the |dm| mixing.

    We explore: within a single l-sector, how does L1 decompose?
    """
    print(f"\n{'='*70}")
    print(f"  WITHIN-L-SECTOR ANALYSIS AT n_max = {n_max}")
    print(f"{'='*70}")

    data = build_fock_graph_data(n_max)
    states = data['states']
    edges = data['edges']
    L1 = data['L1']
    E = data['E']

    # Classify edges by l-sector
    edge_l_sector = np.array([states[edges[k][0]][1] for k in range(E)], dtype=int)
    edge_class = classify_edges(states, edges)

    unique_l = sorted(set(edge_l_sector.tolist()))

    for l_val in unique_l:
        mask = edge_l_sector == l_val
        indices = np.where(mask)[0]
        n_edges_l = len(indices)
        if n_edges_l < 2:
            continue

        # Extract sub-block of L1 for this l-sector
        L1_sub = L1[np.ix_(indices, indices)]

        # Classify edges within this sector
        sub_abs_dm = edge_class['abs_delta_m'][indices]
        n_angular = np.sum(sub_abs_dm == 1)
        n_radial = np.sum(sub_abs_dm == 0)

        print(f"\n  l = {l_val} sector: {n_edges_l} edges ({n_angular} angular, {n_radial} radial)")

        # Check within-sector |dm| commutator
        if n_angular > 0 and n_radial > 0:
            M_sub = np.diag(sub_abs_dm.astype(float))
            comm_sub = L1_sub @ M_sub - M_sub @ L1_sub
            frob_comm = np.linalg.norm(comm_sub, 'fro')
            frob_L1_sub = np.linalg.norm(L1_sub, 'fro')
            print(f"    ||[L1_sub, |dm|]||_F = {frob_comm:.4f} (||L1_sub||_F = {frob_L1_sub:.4f})")
            print(f"    Relative: {frob_comm/frob_L1_sub:.4f}" if frob_L1_sub > 0 else "")

            # Show the actual cross-block entries
            angular_idx = np.where(sub_abs_dm == 1)[0]
            radial_idx = np.where(sub_abs_dm == 0)[0]
            cross = L1_sub[np.ix_(angular_idx, radial_idx)]
            print(f"    Cross-block (angular x radial) sub-matrix ({len(angular_idx)}x{len(radial_idx)}):")
            if cross.size <= 100:  # only print if small
                for row in cross:
                    print(f"      [{', '.join(f'{x:5.1f}' for x in row)}]")

            # Check what L1 eigenvalues look like in this sector
            evals = np.sort(np.linalg.eigvalsh(L1_sub))
            print(f"    L1 eigenvalues in l={l_val}: {[f'{e:.4f}' for e in evals]}")

        # Also label edges by their m-value (using the average m of the two nodes)
        # For angular edges: m_avg = (m + m+1)/2 = m + 0.5
        # For radial edges: m_avg = m
        m_labels = []
        for k_local, k_global in enumerate(indices):
            i, j = edges[k_global]
            n_i, l_i, m_i = states[i]
            n_j, l_j, m_j = states[j]
            if m_i == m_j:  # radial
                m_labels.append(f"R:m={m_i}")
            else:  # angular
                m_labels.append(f"A:m={min(m_i,m_j)}->{max(m_i,m_j)}")
        print(f"    Edge labels: {m_labels}")


def main():
    """Run the full analysis at n_max = 2, 3, 4."""
    all_results = {}

    for n_max in [2, 3, 4]:
        results = run_analysis(n_max)
        all_results[f'n_max_{n_max}'] = results

    # =====================================================================
    # Summary
    # =====================================================================
    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"{'='*70}")

    print("\n  Does L₁ block-diagonalize by |Δm| (photon angular momentum)?")
    print(f"  {'n_max':<6} {'E':<6} {'[L₁,|Δm|]=0?':<15} {'Block-diag frac':<18} {'Eigvec pure frac'}")
    print(f"  {'-'*60}")
    for n_max in [2, 3, 4]:
        r = all_results[f'n_max_{n_max}']
        comm = r['commutator_abs_delta_m']
        eigv = r['eigenvector_localization_abs_dm']
        blk = r['block_structure_abs_delta_m']
        print(f"  {n_max:<6} {r['E']:<6} {str(comm['commutes_exactly']):<15} "
              f"{blk['fraction_within_block']:<18.6f} {eigv['fraction_pure']:.4f}")

    print("\n  Does L₁ block-diagonalize by l-sector?")
    print(f"  {'n_max':<6} {'E':<6} {'Block-diag?':<15} {'Within-block frac'}")
    print(f"  {'-'*50}")
    for n_max in [2, 3, 4]:
        r = all_results[f'n_max_{n_max}']
        blk = r['block_structure_l_sector']
        print(f"  {n_max:<6} {r['E']:<6} {str(blk['is_block_diagonal']):<15} "
              f"{blk['fraction_within_block']:.6f}")

    print("\n  Does L₁ block-diagonalize by Δn (angular vs radial)?")
    print(f"  {'n_max':<6} {'E':<6} {'[L₁,Δn]=0?':<15} {'Block-diag frac':<18}")
    print(f"  {'-'*50}")
    for n_max in [2, 3, 4]:
        r = all_results[f'n_max_{n_max}']
        comm = r['commutator_delta_n']
        blk = r['block_structure_delta_n']
        print(f"  {n_max:<6} {r['E']:<6} {str(comm['commutes_exactly']):<15} "
              f"{blk['fraction_within_block']:<18.6f}")

    # Within-l-sector analysis for n_max=3 and 4
    within_l_sector_analysis(3)
    within_l_sector_analysis(4)

    # Save results
    output_path = Path(__file__).parent / "data" / "photon_ml_commutator.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Clean results for JSON (remove numpy arrays)
    def clean_for_json(obj):
        if isinstance(obj, dict):
            return {k: clean_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [clean_for_json(v) for v in obj]
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        return obj

    with output_path.open('w') as f:
        json.dump(clean_for_json(all_results), f, indent=2)

    print(f"\n  Results saved to: {output_path}")

    # Final interpretation
    print(f"\n{'='*70}")
    print(f"  INTERPRETATION")
    print(f"{'='*70}")

    # Check if l-sector is exact
    l_exact = all(all_results[f'n_max_{n}']['block_structure_l_sector']['is_block_diagonal']
                  for n in [2, 3, 4])

    # Check if |Δm| is exact
    dm_exact = all(all_results[f'n_max_{n}']['commutator_abs_delta_m']['commutes_exactly']
                   for n in [2, 3, 4])

    print(f"\n  L₁ exactly block-diagonal by l-sector: {l_exact}")
    print(f"  L₁ exactly commutes with |Δm|:         {dm_exact}")

    if l_exact:
        print(f"\n  The l-sector decomposition is EXACT. This is expected because")
        print(f"  no edge connects different l values -- the graph is disconnected by l.")
        print(f"  Within each l-sector, L₁ acts independently.")

    if dm_exact:
        print(f"\n  |Δm| commutes exactly with L₁! This means photon modes on the")
        print(f"  Fock graph have well-defined transferred angular momentum |Δm|.")
        print(f"  The photon naturally decomposes into |Δm|=0 (radial) and |Δm|=1 (angular) sectors.")
    else:
        print(f"\n  |Δm| does NOT commute exactly with L₁. The edge Laplacian mixes")
        print(f"  angular (|Δm|=1) and radial (|Δm|=0) edges.")

        # Report the strength of the mixing
        for n_max in [2, 3, 4]:
            blk = all_results[f'n_max_{n_max}']['block_structure_abs_delta_m']
            print(f"    n_max={n_max}: {blk['fraction_cross_block']*100:.2f}% of L₁ Frobenius mass "
                  f"is cross-block (angular<->radial)")


if __name__ == '__main__':
    main()
