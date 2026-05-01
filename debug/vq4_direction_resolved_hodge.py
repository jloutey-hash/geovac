"""
VQ-4: Direction-resolved Hodge decomposition of the Fock graph.
================================================================

Classifies edges of the S^3 Fock scalar graph by their (Dn, Dm)
quantum-number change and builds per-channel Hodge decompositions.

Edge types:
  - T-edges:  Dn = +/-1, Dm = 0  (inter-shell "radial" transitions)
  - L-edges:  Dn = 0,    Dm = +1 (intra-shell "angular" transitions)

In the canonical edge orientation (lower-index -> higher-index), ALL
L-edges have Dm = +1 because the state ordering places m in ascending
order within each (n, l) sub-shell.  There are therefore only TWO
non-empty direction channels for undirected edges:

  Channel T (radial):  Dm = 0, |Dn| >= 1
  Channel L (angular): Dm = +1, Dn = 0

The "L-" channel (Dm = -1) is EMPTY by construction -- each undirected
L-edge is traversed in the m -> m+1 direction by the canonical
orientation.  This is not a physics limitation but a consequence of
the canonical edge ordering.  We document this and focus the analysis
on the T vs L decomposition, which IS physically meaningful.

For each channel mu in {T, L}:
  1. Extract the sub-incidence matrix B_mu (V x E_mu)
  2. Build the sub-node Laplacian L0_mu = B_mu . B_mu^T
  3. Build the sub-edge Laplacian L1_mu = B_mu^T . B_mu
  4. Compute spectra, Betti numbers, photon propagator G_gamma_mu

Key physics questions:
  - Is the photon propagation isotropic or anisotropic?
  - Does G_gamma^{T} + G_gamma^{L} = G_gamma^{full}?
  - Are the T and L channels individually connected?

Author: VQ-4 sprint
Date: 2026-05-01
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# Add project root to path
project_root = Path(__file__).resolve().parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from geovac.lattice import GeometricLattice


# ---------------------------------------------------------------------------
# JSON serialization helper
# ---------------------------------------------------------------------------

def _json_default(obj):
    """Convert numpy types to Python types for JSON serialization."""
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def build_fock_graph_data(n_max: int) -> Dict:
    """Build the Fock graph and extract states, edges, and classifications.

    Returns dict with:
      states: list of (n, l, m) tuples
      edges: list of (v1, v2) with v1 < v2
      edge_info: list of dicts with dn, dl, dm, channel
      V, E: counts
    """
    lat = GeometricLattice(max_n=n_max, topological_weights=False)
    states = lat.states
    V = len(states)

    # Extract undirected edges from adjacency
    adj = lat.adjacency
    rows, cols = adj.nonzero()
    edge_set: set = set()
    for r, c in zip(rows, cols):
        if int(r) < int(c):
            edge_set.add((int(r), int(c)))
    edges = sorted(edge_set)
    E = len(edges)

    # Classify each edge
    edge_info = []
    for e_idx, (v1, v2) in enumerate(edges):
        s1 = states[v1]
        s2 = states[v2]
        dn = s2[0] - s1[0]
        dl = s2[1] - s1[1]
        dm = s2[2] - s1[2]

        if dn != 0 and dm == 0:
            channel = 'T'
        elif dn == 0 and dm == 1:
            channel = 'L+'
        elif dn == 0 and dm == -1:
            channel = 'L-'
        else:
            channel = 'unknown'

        edge_info.append({
            'edge_idx': e_idx,
            'v1': v1, 'v2': v2,
            's1': list(s1), 's2': list(s2),
            'dn': dn, 'dl': dl, 'dm': dm,
            'channel': channel,
        })

    return {
        'states': states,
        'edges': edges,
        'edge_info': edge_info,
        'V': V,
        'E': E,
    }


def build_full_incidence(V: int, edges: List[Tuple[int, int]]) -> np.ndarray:
    """Build the full V x E signed incidence matrix B.

    Convention: edge k = (i, j) with i < j ->
      B[i, k] = +1 (tail/source)
      B[j, k] = -1 (head/sink)
    """
    E = len(edges)
    B = np.zeros((V, E), dtype=np.float64)
    for k, (i, j) in enumerate(edges):
        B[i, k] = 1.0
        B[j, k] = -1.0
    return B


def build_sub_incidence(V: int,
                        all_edges: List[Tuple[int, int]],
                        edge_indices: List[int]) -> np.ndarray:
    """Build the V x E_mu sub-incidence matrix for a channel.

    Parameters
    ----------
    V : int
        Total number of nodes.
    all_edges : list of (v1, v2) pairs
        Full edge list.
    edge_indices : list of int
        Indices into all_edges for this channel.

    Returns
    -------
    B_mu : numpy array (V x E_mu)
    """
    E_mu = len(edge_indices)
    if E_mu == 0:
        return np.zeros((V, 0), dtype=np.float64)
    B = np.zeros((V, E_mu), dtype=np.float64)
    for k, e_idx in enumerate(edge_indices):
        i, j = all_edges[e_idx]
        B[i, k] = 1.0
        B[j, k] = -1.0
    return B


def compute_hodge_data(B_mu: np.ndarray, V: int, channel_name: str) -> Dict:
    """Compute the full Hodge decomposition for a sub-incidence matrix.

    Returns a dict with:
      L0_spectrum, L1_spectrum, beta_0, beta_1,
      G_gamma (photon propagator as pseudoinverse of L1),
      node_coverage, connectivity info.
    """
    E_mu = B_mu.shape[1]

    if E_mu == 0:
        return {
            'channel': channel_name,
            'E_mu': 0,
            'L0_spectrum': [],
            'L1_spectrum': [],
            'L0_nonzero_spectrum': [],
            'L1_nonzero_spectrum': [],
            'beta_0_full_V': V,
            'beta_1_formula': -V,
            'beta_1_kernel': 0,
            'sub_graph_V': 0,
            'sub_graph_beta_0': 0,
            'is_connected_sub': False,
            'G_gamma': [],
            'G_gamma_trace': 0.0,
            'G_gamma_frobenius': 0.0,
            'node_coverage': 0,
            'nodes_touched': [],
            'notes': 'Empty channel (no edges).',
        }

    # Node Laplacian L0 = B . B^T (V x V)
    L0 = B_mu @ B_mu.T

    # Edge Laplacian L1 = B^T . B (E_mu x E_mu)
    L1 = B_mu.T @ B_mu

    # Spectra
    ev_L0 = np.sort(np.linalg.eigvalsh(L0))
    ev_L1 = np.sort(np.linalg.eigvalsh(L1))

    # beta_0 = number of zero eigenvalues of L0 = connected components
    beta_0 = int(np.sum(np.abs(ev_L0) < 1e-10))

    # beta_1 = E_mu - V + beta_0  (but also = number of zero evals of L1)
    beta_1_formula = E_mu - V + beta_0
    beta_1_kernel = int(np.sum(np.abs(ev_L1) < 1e-10))
    # Use the kernel count (more reliable for sub-graphs)
    beta_1 = beta_1_kernel

    # Photon propagator G_gamma = L1^+ (Moore-Penrose pseudoinverse)
    if E_mu > 0:
        G_gamma = np.linalg.pinv(L1)
    else:
        G_gamma = np.zeros((0, 0))

    # Node coverage: which nodes are incident to at least one edge in this channel
    nodes_touched = set()
    for k in range(E_mu):
        for v in range(V):
            if abs(B_mu[v, k]) > 0.5:
                nodes_touched.add(v)
    nodes_touched = sorted(nodes_touched)

    # Check if the sub-graph is connected (beta_0 relative to touched nodes)
    # For the sub-graph restricted to nodes that appear:
    V_touched = len(nodes_touched)
    if V_touched == 0:
        is_connected = False
        sub_beta_0 = 0
    else:
        # Build the adjacency of the sub-graph restricted to touched nodes
        node_map = {v: i for i, v in enumerate(nodes_touched)}
        sub_adj = np.zeros((V_touched, V_touched))
        for k in range(E_mu):
            endpoints = [v for v in range(V) if abs(B_mu[v, k]) > 0.5]
            if len(endpoints) == 2:
                i_sub = node_map[endpoints[0]]
                j_sub = node_map[endpoints[1]]
                sub_adj[i_sub, j_sub] = 1
                sub_adj[j_sub, i_sub] = 1
        sub_L0 = np.diag(sub_adj.sum(axis=1)) - sub_adj
        sub_ev = np.linalg.eigvalsh(sub_L0)
        sub_beta_0 = int(np.sum(np.abs(sub_ev) < 1e-10))
        is_connected = (sub_beta_0 == 1)

    return {
        'channel': channel_name,
        'E_mu': E_mu,
        'L0_spectrum': ev_L0.tolist(),
        'L1_spectrum': ev_L1.tolist(),
        'L0_nonzero_spectrum': [float(x) for x in ev_L0 if abs(x) > 1e-10],
        'L1_nonzero_spectrum': [float(x) for x in ev_L1 if abs(x) > 1e-10],
        'beta_0_full_V': beta_0,
        'beta_1_formula': beta_1_formula,
        'beta_1_kernel': beta_1,
        'sub_graph_V': V_touched,
        'sub_graph_beta_0': sub_beta_0,
        'is_connected_sub': is_connected,
        'G_gamma': G_gamma.tolist(),
        'G_gamma_trace': float(np.trace(G_gamma)),
        'G_gamma_frobenius': float(np.linalg.norm(G_gamma, 'fro')),
        'node_coverage': V_touched,
        'nodes_touched': nodes_touched,
    }


def analyze_nmax(n_max: int) -> Dict:
    """Run the full direction-resolved Hodge analysis at a given n_max.

    Returns a comprehensive results dict.
    """
    print(f"\n{'='*60}")
    print(f"  VQ-4: Direction-resolved Hodge at n_max = {n_max}")
    print(f"{'='*60}")

    # Build graph and classify edges
    gdata = build_fock_graph_data(n_max)
    V = gdata['V']
    E = gdata['E']
    states = gdata['states']
    edges = gdata['edges']
    edge_info = gdata['edge_info']

    print(f"\nGraph: V = {V}, E = {E}")

    # Classify edges by channel
    channels = {'T': [], 'L+': [], 'L-': []}
    for ei in edge_info:
        ch = ei['channel']
        if ch in channels:
            channels[ch].append(ei['edge_idx'])
        else:
            print(f"  WARNING: unknown channel for edge {ei}")

    for ch_name in ['T', 'L+', 'L-']:
        print(f"  Channel {ch_name}: {len(channels[ch_name])} edges")

    # Print edge details
    print("\nEdge classification:")
    for ei in edge_info:
        print(f"  edge {ei['edge_idx']}: {tuple(ei['s1'])} -> {tuple(ei['s2'])}  "
              f"Dn={ei['dn']:+d}  Dm={ei['dm']:+d}  channel={ei['channel']}")

    # Build full incidence matrix and Hodge data
    B_full = build_full_incidence(V, edges)
    L0_full = B_full @ B_full.T
    L1_full = B_full.T @ B_full
    G_full = np.linalg.pinv(L1_full) if E > 0 else np.zeros((0, 0))

    ev_L0_full = np.sort(np.linalg.eigvalsh(L0_full))
    ev_L1_full = np.sort(np.linalg.eigvalsh(L1_full))

    beta_0_full = int(np.sum(np.abs(ev_L0_full) < 1e-10))
    beta_1_full = int(np.sum(np.abs(ev_L1_full) < 1e-10))

    print(f"\nFull graph:")
    print(f"  L0 spectrum: {np.round(ev_L0_full, 6).tolist()}")
    print(f"  L1 spectrum: {np.round(ev_L1_full, 6).tolist()}")
    print(f"  beta_0 = {beta_0_full}, beta_1 = {beta_1_full}")

    # Per-channel Hodge decomposition
    channel_results = {}
    for ch_name in ['T', 'L+', 'L-']:
        ch_edges = channels[ch_name]
        B_mu = build_sub_incidence(V, edges, ch_edges)
        hodge = compute_hodge_data(B_mu, V, ch_name)
        channel_results[ch_name] = hodge

        print(f"\nChannel {ch_name} (E_mu = {hodge['E_mu']}):")
        if hodge['E_mu'] == 0:
            print(f"  EMPTY -- no edges in this channel")
        else:
            print(f"  Sub-graph nodes: {hodge['sub_graph_V']} of {V}")
            print(f"  Connected (sub): {hodge['is_connected_sub']} "
                  f"(sub_beta_0 = {hodge['sub_graph_beta_0']})")
            print(f"  L1 spectrum: {[round(x, 6) for x in hodge['L1_spectrum']]}")
            print(f"  L1 nonzero: {[round(x, 6) for x in hodge['L1_nonzero_spectrum']]}")
            print(f"  beta_1 (kernel): {hodge['beta_1_kernel']}")
            print(f"  G_gamma trace: {hodge['G_gamma_trace']:.6f}")

    # --- Isotropy check ---
    print("\n--- Isotropy check ---")
    nz_T = sorted(channel_results['T']['L1_nonzero_spectrum'])
    nz_Lp = sorted(channel_results['L+']['L1_nonzero_spectrum'])
    nz_Lm = sorted(channel_results['L-']['L1_nonzero_spectrum'])

    if nz_T and nz_Lp:
        print(f"  T  nonzero L1: {[round(x, 6) for x in nz_T]}")
        print(f"  L+ nonzero L1: {[round(x, 6) for x in nz_Lp]}")
        if nz_T == nz_Lp:
            isotropy = 'ISOTROPIC (T = L+ spectra identical)'
        else:
            isotropy = 'ANISOTROPIC (T and L+ spectra differ)'
        print(f"  Verdict: {isotropy}")
    else:
        isotropy = 'N/A (one or both channels empty)'
        print(f"  Verdict: {isotropy}")

    # --- Propagator sum check ---
    print("\n--- Propagator sum check: G_T + G_L+ + G_L- =? G_full ---")
    if E > 0:
        # Build per-channel propagators in the FULL edge space
        # The per-channel sub-incidence matrices produce propagators in their
        # own E_mu-dimensional space. To compare with G_full, we need to
        # embed them back into the full E-dimensional edge space.
        G_sum = np.zeros((E, E))
        for ch_name in ['T', 'L+', 'L-']:
            ch_edges = channels[ch_name]
            if len(ch_edges) == 0:
                continue
            # Build embedding: the sub-propagator in E_mu space maps to
            # a submatrix of E x E via the edge indices
            hodge = channel_results[ch_name]
            G_ch = np.array(hodge['G_gamma'])
            for i_local, e_i in enumerate(ch_edges):
                for j_local, e_j in enumerate(ch_edges):
                    G_sum[e_i, e_j] += G_ch[i_local, j_local]

        diff = np.linalg.norm(G_sum - G_full, 'fro')
        rel_diff = diff / max(np.linalg.norm(G_full, 'fro'), 1e-30)
        propagator_sum_matches = diff < 1e-10

        print(f"  ||G_sum - G_full||_F = {diff:.2e}")
        print(f"  Relative diff = {rel_diff:.2e}")
        print(f"  Match: {propagator_sum_matches}")

        # Diagonal comparison
        diag_full = np.diag(G_full)
        diag_sum = np.diag(G_sum)
        print(f"\n  Per-edge diagonal comparison (G_full vs G_sum):")
        for e_idx in range(E):
            ei = edge_info[e_idx]
            ch = ei['channel']
            print(f"    edge {e_idx} ({ch}): G_full={diag_full[e_idx]:.6f}  "
                  f"G_sum={diag_sum[e_idx]:.6f}  "
                  f"diff={abs(diag_full[e_idx]-diag_sum[e_idx]):.2e}")
    else:
        propagator_sum_matches = True
        diff = 0.0
        rel_diff = 0.0

    # --- Cross-channel coupling in full propagator ---
    print("\n--- Cross-channel coupling in G_full ---")
    if E > 0:
        T_edges = channels['T']
        L_edges = channels['L+']  # L- is empty

        # G_full[T, L] block
        cross_block_norm = 0.0
        if T_edges and L_edges:
            cross_block = G_full[np.ix_(T_edges, L_edges)]
            cross_block_norm = float(np.linalg.norm(cross_block, 'fro'))
            T_block = G_full[np.ix_(T_edges, T_edges)]
            L_block = G_full[np.ix_(L_edges, L_edges)]
            T_block_norm = float(np.linalg.norm(T_block, 'fro'))
            L_block_norm = float(np.linalg.norm(L_block, 'fro'))
            total_norm = float(np.linalg.norm(G_full, 'fro'))

            print(f"  G_full[T,T] Frobenius: {T_block_norm:.6f}")
            print(f"  G_full[L,L] Frobenius: {L_block_norm:.6f}")
            print(f"  G_full[T,L] Frobenius: {cross_block_norm:.6f}")
            print(f"  G_full total Frobenius: {total_norm:.6f}")
            print(f"  Cross/Total ratio: {cross_block_norm/max(total_norm, 1e-30):.6f}")

            # Is the full propagator block-diagonal in the T/L decomposition?
            is_block_diagonal = cross_block_norm < 1e-10
            print(f"  Block-diagonal: {is_block_diagonal}")
        else:
            is_block_diagonal = True
            T_block_norm = 0.0
            L_block_norm = 0.0
            print(f"  Only one non-empty channel -- trivially block-diagonal")
    else:
        is_block_diagonal = True
        cross_block_norm = 0.0
        T_block_norm = 0.0
        L_block_norm = 0.0

    # --- Compile results ---
    result = {
        'n_max': n_max,
        'V': V,
        'E': E,
        'edge_classification': [
            {
                'edge_idx': ei['edge_idx'],
                's1': ei['s1'], 's2': ei['s2'],
                'dn': ei['dn'], 'dm': ei['dm'],
                'channel': ei['channel'],
            }
            for ei in edge_info
        ],
        'channel_counts': {ch: len(channels[ch]) for ch in ['T', 'L+', 'L-']},
        'full_graph': {
            'L0_spectrum': ev_L0_full.tolist(),
            'L1_spectrum': ev_L1_full.tolist(),
            'beta_0': beta_0_full,
            'beta_1': beta_1_full,
            'G_gamma_trace': float(np.trace(G_full)) if E > 0 else 0.0,
            'G_gamma_frobenius': float(np.linalg.norm(G_full, 'fro')) if E > 0 else 0.0,
        },
        'channels': {},
        'isotropy_check': {
            'T_nonzero_L1': nz_T,
            'Lplus_nonzero_L1': nz_Lp,
            'Lminus_nonzero_L1': nz_Lm,
            'verdict': isotropy,
        },
        'propagator_sum_check': {
            'frobenius_diff': float(diff),
            'relative_diff': float(rel_diff),
            'matches': propagator_sum_matches,
        },
        'cross_channel_coupling': {
            'G_TT_frobenius': float(T_block_norm),
            'G_LL_frobenius': float(L_block_norm),
            'G_TL_frobenius': float(cross_block_norm),
            'is_block_diagonal': is_block_diagonal,
        },
    }

    # Add per-channel results (omit large matrices for JSON)
    for ch_name in ['T', 'L+', 'L-']:
        ch = channel_results[ch_name]
        result['channels'][ch_name] = {
            'E_mu': ch['E_mu'],
            'sub_graph_V': ch['sub_graph_V'],
            'sub_graph_beta_0': ch['sub_graph_beta_0'],
            'is_connected_sub': ch['is_connected_sub'],
            'L0_spectrum': ch['L0_spectrum'],
            'L1_spectrum': ch['L1_spectrum'],
            'L0_nonzero_spectrum': ch['L0_nonzero_spectrum'],
            'L1_nonzero_spectrum': ch['L1_nonzero_spectrum'],
            'beta_0_full_V': ch['beta_0_full_V'],
            'beta_1_kernel': ch['beta_1_kernel'],
            'G_gamma_trace': ch['G_gamma_trace'],
            'G_gamma_frobenius': ch['G_gamma_frobenius'],
            'nodes_touched': ch['nodes_touched'],
        }

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Run the VQ-4 direction-resolved Hodge analysis."""

    results = {
        'description': (
            'VQ-4: Direction-resolved Hodge decomposition of the Fock graph. '
            'Edges classified by (Dn, Dm) into T (radial) and L (angular) '
            'channels. Per-channel sub-incidence, sub-Laplacians, and '
            'photon propagators computed.'
        ),
    }

    for n_max in [2, 3]:
        key = f'n_max_{n_max}'
        results[key] = analyze_nmax(n_max)

    # -----------------------------------------------------------------------
    # Cross-n_max comparison
    # -----------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"  Cross-n_max comparison")
    print(f"{'='*60}")

    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        print(f"\nn_max = {n_max}:")
        print(f"  V = {r['V']}, E = {r['E']}")
        cc = r['channel_counts']
        print(f"  T edges: {cc['T']}, L+ edges: {cc['L+']}, L- edges: {cc['L-']}")
        print(f"  Full beta_0 = {r['full_graph']['beta_0']}, "
              f"beta_1 = {r['full_graph']['beta_1']}")

        for ch_name in ['T', 'L+']:
            ch = r['channels'][ch_name]
            if ch['E_mu'] > 0:
                print(f"  Channel {ch_name}: E_mu={ch['E_mu']}, "
                      f"sub_V={ch['sub_graph_V']}, "
                      f"connected={ch['is_connected_sub']}, "
                      f"beta_1={ch['beta_1_kernel']}, "
                      f"G_trace={ch['G_gamma_trace']:.6f}")

    # -----------------------------------------------------------------------
    # Summary findings
    # -----------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"  SUMMARY")
    print(f"{'='*60}")

    # Key finding 1: L- channel is always empty
    print("\n1. L- channel is EMPTY at both n_max=2 and n_max=3.")
    print("   Reason: canonical edge orientation (i < j) always gives Dm = +1")
    print("   for angular edges because states with higher m have higher index.")
    print("   The graph has TWO physically distinct edge types: T and L.")

    # Key finding 2: Isotropy
    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        print(f"\n2. Isotropy at n_max={n_max}: {r['isotropy_check']['verdict']}")

    # Key finding 3: Propagator sum
    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        ps = r['propagator_sum_check']
        print(f"\n3. Propagator sum at n_max={n_max}: "
              f"matches={ps['matches']}, diff={ps['frobenius_diff']:.2e}")

    # Key finding 4: Cross-channel coupling
    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        cc = r['cross_channel_coupling']
        print(f"\n4. Cross-channel in G_full at n_max={n_max}: "
              f"block_diag={cc['is_block_diagonal']}, "
              f"G_TL={cc['G_TL_frobenius']:.6f}")

    # Save results (convert numpy types for JSON compatibility)
    output_path = Path(__file__).resolve().parent / 'data' / 'vq4_direction_resolved_hodge.json'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=_json_default)
    print(f"\nResults saved to: {output_path}")

    # -----------------------------------------------------------------------
    # Write memo
    # -----------------------------------------------------------------------
    write_memo(results)


def write_memo(results: Dict) -> None:
    """Write the analysis memo."""
    memo_path = Path(__file__).resolve().parent / 'vq4_direction_resolved_hodge_memo.md'

    r2 = results['n_max_2']
    r3 = results['n_max_3']

    lines = []
    lines.append("# VQ-4: Direction-Resolved Hodge Decomposition of the Fock Graph")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append("The Fock scalar graph edges are classified by their quantum-number")
    lines.append("change (Dn, Dm) into two physically distinct channels:")
    lines.append("")
    lines.append("- **T-channel** (radial): Dn = +/-1, Dm = 0 (inter-shell transitions)")
    lines.append("- **L-channel** (angular): Dn = 0, Dm = +/-1 (intra-shell transitions)")
    lines.append("")
    lines.append("Per-channel Hodge decompositions (sub-incidence B_mu, sub-Laplacians")
    lines.append("L0_mu and L1_mu, photon propagator G_mu = L1_mu^+) are computed at")
    lines.append("n_max = 2 and n_max = 3.")
    lines.append("")

    # Finding 0: Three channels reduce to two
    lines.append("## Finding 0: Three Channels Reduce to Two")
    lines.append("")
    lines.append("The naive classification into three channels (L+, L-, T) yields only")
    lines.append("two non-empty channels. The L- channel (Dm = -1) is **empty** at both")
    lines.append("n_max = 2 and n_max = 3. This is structural: the canonical edge")
    lines.append("orientation (lower-index node -> higher-index node) always assigns")
    lines.append("Dm = +1 to angular edges because the state ordering places m in")
    lines.append("ascending order within each (n, l) sub-shell. Each undirected L-edge")
    lines.append("represents BOTH the m -> m+1 and m+1 -> m transitions; the")
    lines.append("\"direction\" is a labeling convention, not a physical distinction.")
    lines.append("")
    lines.append("The physically meaningful decomposition is **T vs L** (radial vs angular).")
    lines.append("")

    # Finding 1: Edge counts
    lines.append("## Finding 1: Edge Counts and Asymmetry")
    lines.append("")
    lines.append("| n_max | V | E | E_T | E_L | E_T/E_L |")
    lines.append("|:------|:--|:--|:----|:----|:--------|")
    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        E_T = r['channel_counts']['T']
        E_L = r['channel_counts']['L+']
        ratio = f"{E_T/E_L:.3f}" if E_L > 0 else "inf"
        lines.append(f"| {n_max} | {r['V']} | {r['E']} | {E_T} | {E_L} | {ratio} |")
    lines.append("")
    lines.append("At n_max = 2: 1 T-edge and 2 L-edges. The T-channel is the minority.")
    lines.append("At n_max = 3: 5 T-edges and 8 L-edges. L-edges dominate because the")
    lines.append("l = 2 shell contributes 4 angular edges (a path graph of length 5).")
    lines.append("")

    # Finding 2: Connectivity
    lines.append("## Finding 2: Per-Channel Connectivity")
    lines.append("")
    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        lines.append(f"**n_max = {n_max}:**")
        lines.append("")
        for ch_name in ['T', 'L+']:
            ch = r['channels'][ch_name]
            if ch['E_mu'] > 0:
                lines.append(f"- Channel {ch_name}: {ch['sub_graph_V']} nodes touched, "
                             f"{ch['sub_graph_beta_0']} connected components, "
                             f"connected = {ch['is_connected_sub']}")
            else:
                lines.append(f"- Channel {ch_name}: empty")
        lines.append("")

    lines.append("The T-channel connects nodes across different n-shells (same l, m).")
    lines.append("The L-channel connects nodes within the same (n, l) sub-shell across m.")
    lines.append("Neither channel alone connects the full graph -- photon propagation")
    lines.append("requires BOTH radial and angular hops.")
    lines.append("")

    # Finding 3: Spectral anisotropy
    lines.append("## Finding 3: Spectral Anisotropy")
    lines.append("")
    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        iso = r['isotropy_check']
        lines.append(f"**n_max = {n_max}:** {iso['verdict']}")
        lines.append("")
        lines.append(f"- T  nonzero L1 eigenvalues: "
                     f"{[round(x, 6) for x in iso['T_nonzero_L1']]}")
        lines.append(f"- L+ nonzero L1 eigenvalues: "
                     f"{[round(x, 6) for x in iso['Lplus_nonzero_L1']]}")
        lines.append("")

    lines.append("The T and L channels have **different spectra**, confirming the")
    lines.append("photon propagation is inherently anisotropic on the Fock graph.")
    lines.append("The T-channel (radial) carries different eigenvalues than the")
    lines.append("L-channel (angular). This is a structural feature of the graph")
    lines.append("topology, not a dynamical effect.")
    lines.append("")

    # Finding 4: Propagator sum
    lines.append("## Finding 4: Propagator Additivity")
    lines.append("")
    lines.append("Does the sum of per-channel propagators (embedded in the full edge")
    lines.append("space) equal the full propagator?")
    lines.append("")
    lines.append("| n_max | ||G_T + G_L - G_full||_F | Relative diff | Match? |")
    lines.append("|:------|:------------------------|:--------------|:-------|")
    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        ps = r['propagator_sum_check']
        lines.append(f"| {n_max} | {ps['frobenius_diff']:.2e} | "
                     f"{ps['relative_diff']:.2e} | {ps['matches']} |")
    lines.append("")

    # Interpret the propagator sum result
    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        ps = r['propagator_sum_check']
        if ps['matches']:
            lines.append(f"At n_max = {n_max}: **YES** -- the propagator decomposes additively.")
        else:
            lines.append(f"At n_max = {n_max}: **NO** -- the propagator does NOT decompose "
                         f"additively (diff = {ps['frobenius_diff']:.2e}).")
    lines.append("")

    # Finding 5: Cross-channel coupling
    lines.append("## Finding 5: Cross-Channel Coupling in G_full")
    lines.append("")
    lines.append("The full propagator G_full can be partitioned into T-T, L-L, and")
    lines.append("T-L blocks. If the T-L cross-block is nonzero, the photon")
    lines.append("propagation mixes radial and angular channels.")
    lines.append("")
    lines.append("| n_max | ||G_{TT}|| | ||G_{LL}|| | ||G_{TL}|| | Block-diag? |")
    lines.append("|:------|:-----------|:-----------|:-----------|:------------|")
    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        cc = r['cross_channel_coupling']
        lines.append(f"| {n_max} | {cc['G_TT_frobenius']:.6f} | "
                     f"{cc['G_LL_frobenius']:.6f} | "
                     f"{cc['G_TL_frobenius']:.6f} | "
                     f"{cc['is_block_diagonal']} |")
    lines.append("")

    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        cc = r['cross_channel_coupling']
        if cc['is_block_diagonal']:
            lines.append(f"At n_max = {n_max}: G_full is **block-diagonal** in T/L channels.")
            lines.append("The photon propagates independently in each direction channel.")
        else:
            lines.append(f"At n_max = {n_max}: G_full has **nonzero cross-channel coupling**.")
            lines.append("The photon propagation mixes radial and angular channels.")
    lines.append("")

    # Finding 6: Per-channel propagator traces
    lines.append("## Finding 6: Per-Channel Propagator Traces")
    lines.append("")
    lines.append("| n_max | Tr(G_T) | Tr(G_L) | Tr(G_full) | Tr(G_T)+Tr(G_L) |")
    lines.append("|:------|:--------|:--------|:-----------|:----------------|")
    for n_max in [2, 3]:
        r = results[f'n_max_{n_max}']
        tr_T = r['channels']['T']['G_gamma_trace']
        tr_L = r['channels']['L+']['G_gamma_trace']
        tr_full = r['full_graph']['G_gamma_trace']
        lines.append(f"| {n_max} | {tr_T:.6f} | {tr_L:.6f} | "
                     f"{tr_full:.6f} | {tr_T + tr_L:.6f} |")
    lines.append("")

    # Structural interpretation
    lines.append("## Structural Interpretation")
    lines.append("")
    lines.append("The direction-resolved Hodge decomposition reveals that the Fock graph")
    lines.append("photon has two structurally distinct propagation channels:")
    lines.append("")
    lines.append("1. **Radial (T)**: photon hops between shells (Dn = +/-1, same l and m).")
    lines.append("   These are the edges that connect different energy levels.")
    lines.append("")
    lines.append("2. **Angular (L)**: photon hops within a shell (Dm = +/-1, same n and l).")
    lines.append("   These are the edges that rotate the magnetic quantum number.")
    lines.append("")
    lines.append("The spectra are generically different (anisotropic photon). The full")
    lines.append("propagator's T-L cross-block determines whether radial and angular")
    lines.append("photon modes couple or propagate independently.")
    lines.append("")
    lines.append("This gives the scalar graph photon effective \"polarization\" structure")
    lines.append("from pure topology, without introducing explicit vector labels. The")
    lines.append("anisotropy is a consequence of the Fock graph having two geometrically")
    lines.append("distinct edge types built into its quantum-number lattice.")
    lines.append("")

    # Data files
    lines.append("## Data Files")
    lines.append("")
    lines.append("- `debug/data/vq4_direction_resolved_hodge.json` -- full numerical results")
    lines.append("- `debug/vq4_direction_resolved_hodge.py` -- computation script")
    lines.append("")

    memo_text = "\n".join(lines) + "\n"
    with open(memo_path, 'w') as f:
        f.write(memo_text)
    print(f"Memo saved to: {memo_path}")


if __name__ == '__main__':
    main()
