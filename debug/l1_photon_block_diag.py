"""
L₁ photon block-diagonalization diagnostic.
============================================

Investigate whether the edge Laplacian L₁ = B^T B on the Fock/Dirac graph
block-diagonalizes by photon quantum numbers (q, m_q) derived from the
edge endpoints.

Background:
- The Fock graph has nodes labeled by (n, l, m) quantum numbers.
- Edges connect adjacent shells (Δn = ±1) with selection rules, and
  magnetic transitions (Δm = ±1) within shells.
- Each edge implicitly carries photon quantum numbers:
  * T-edges (Δn = ±1, Δl = 0, Δm = 0): q = 0, m_q = 0 (monopole)
  * L-edges (Δn = 0, Δl = 0, Δm = ±1): q = 0, m_q = ±1 (but l unchanged)

  Actually on the scalar Fock graph, edges have Δl = 0, so all edges
  preserve l. The photon carried by an edge (n₁,l₁,m₁)→(n₂,l₂,m₂)
  has m_q = m₂ - m₁ and q is constrained by the angular momentum
  triangle |l₁-l₂| ≤ q ≤ l₁+l₂. Since l₁ = l₂ on the scalar graph,
  q ∈ {0, 1, ..., 2l} in principle, but the dominant/natural channel
  for same-l transitions is q = 0 (monopole), while m-ladder edges
  have Δm = ±1 which requires q ≥ 1.

  For T-edges: Δn=±1, Δl=0, Δm=0 → m_q = 0. Natural q = 0.
  For L-edges: Δn=0, Δl=0, Δm=±1 → m_q = ±1. Requires q ≥ 1. Natural q = 1.

For the Dirac Rule B graph, edges have Δl = ±1 (dipole), so q = 1
is the dominant channel with m_q = m_j₂ - m_j₁.

The question: does L₁ block-diagonalize by these edge quantum numbers?

Author: Investigation sprint, May 2026.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from geovac.lattice import GeometricLattice
from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels, kappa_to_l
from geovac.ihara_zeta_dirac import build_dirac_s3_graph


# ---------------------------------------------------------------------------
# Helper: build signed incidence and L1 for the Fock graph (numpy path)
# ---------------------------------------------------------------------------

def build_fock_graph_numpy(n_max: int):
    """Build Fock graph topology and return edges, states, B, L1."""
    lat = GeometricLattice(max_n=n_max, topological_weights=False)
    states = lat.states
    V = len(states)
    state_index = {s: i for i, s in enumerate(states)}

    adj = lat.adjacency
    rows, cols = adj.nonzero()
    edge_set = set()
    for r, c in zip(rows, cols):
        if r < c:
            edge_set.add((int(r), int(c)))
    edges = sorted(edge_set)
    E = len(edges)

    # Build signed incidence B (V x E)
    B = np.zeros((V, E), dtype=np.float64)
    for k, (i, j) in enumerate(edges):
        B[i, k] = 1
        B[j, k] = -1

    # Edge Laplacian
    L1 = B.T @ B

    return states, edges, B, L1, V, E


def build_dirac_graph_numpy(n_max: int, rule: str = "B"):
    """Build Dirac graph and return edges, labels, B, L1."""
    A_adj, labels, deg, desc = build_dirac_s3_graph(n_max, rule)
    V = len(labels)

    # Extract edges from adjacency
    edges = []
    for i in range(V):
        for j in range(i + 1, V):
            if A_adj[i, j]:
                edges.append((i, j))
    E = len(edges)

    # Build signed incidence B
    B = np.zeros((V, E), dtype=np.float64)
    for k, (i, j) in enumerate(edges):
        B[i, k] = 1
        B[j, k] = -1

    L1 = B.T @ B
    return labels, edges, B, L1, V, E


# ---------------------------------------------------------------------------
# Edge quantum number labeling
# ---------------------------------------------------------------------------

def label_fock_edges(states, edges):
    """Label each edge of the Fock graph with photon quantum numbers.

    For edge (i, j) connecting state (n₁, l₁, m₁) to (n₂, l₂, m₂):
    - Δn = n₂ - n₁, Δl = l₂ - l₁, Δm = m₂ - m₁
    - m_q = m₂ - m₁ (magnetic quantum number transferred)
    - edge_type: 'T' if Δn ≠ 0, 'L' if Δm ≠ 0
    - q: for same-l edges, if Δm = 0 then q=0 (monopole/scalar),
      if |Δm| = 1 then q = 1 (dipole-like). More precisely, the
      "natural" q is |Δm| (since q ≥ |m_q|).
    """
    edge_labels = []
    for e_idx, (i, j) in enumerate(edges):
        s1 = states[i]
        s2 = states[j]
        n1, l1, m1 = s1
        n2, l2, m2 = s2

        dn = n2 - n1
        dl = l2 - l1
        dm = m2 - m1  # = m_q

        # Classify edge type
        if dn != 0:
            edge_type = 'T'  # radial / shell-to-shell
        else:
            edge_type = 'L'  # angular / magnetic

        # q: minimal angular momentum the "photon" must carry
        # q ≥ |m_q| is the triangle rule constraint.
        # On the scalar Fock graph, all edges have Δl = 0.
        # T-edges: q=0, m_q=0 (monopole transfer)
        # L-edges: q=1, m_q=±1 (need at least q=1 to carry Δm=±1)
        m_q = dm
        if dm == 0 and dl == 0:
            q = 0
        else:
            q = max(abs(dl), abs(dm))  # minimal q from triangle

        edge_labels.append({
            'edge_idx': e_idx,
            'nodes': (i, j),
            'states': (s1, s2),
            'dn': dn, 'dl': dl, 'dm': dm,
            'edge_type': edge_type,
            'q': q,
            'm_q': m_q,
        })

    return edge_labels


def label_dirac_edges(labels, edges, rule="B"):
    """Label Dirac graph edges with photon quantum numbers.

    For Rule B (E1 dipole): Δl = ±1, so q = 1 and m_q derived from
    the angular momentum coupling.
    For Rule A (κ-preserving): Δl = 0, so edges are monopole (q=0)
    for T-type and q=1 for L-type (Δm_j = ±1).
    """
    edge_labels = []
    for e_idx, (i, j) in enumerate(edges):
        a = labels[i]
        b = labels[j]

        dn = b.n_fock - a.n_fock
        la = kappa_to_l(a.kappa)
        lb = kappa_to_l(b.kappa)
        dl = lb - la
        dmj2 = b.two_m_j - a.two_m_j  # 2 * Δm_j

        # m_q = m_j2 - m_j1 (half-integer in general)
        # For analysis, use 2*m_q = 2*m_j2 - 2*m_j1
        two_m_q = dmj2

        if rule == "A":
            # Rule A: Δκ = 0, so Δl = 0
            if dn != 0 and dmj2 == 0:
                edge_type = 'T'
                q = 0
            else:
                edge_type = 'L'
                q = 1  # Δm_j = ±1 requires q ≥ 1
        else:
            # Rule B: Δl = ±1 (E1 dipole), q = 1
            q = 1
            if dn != 0:
                edge_type = 'T'
            else:
                edge_type = 'L'

        edge_labels.append({
            'edge_idx': e_idx,
            'nodes': (i, j),
            'dn': dn, 'dl': dl,
            'two_dmj': dmj2,
            'edge_type': edge_type,
            'q': q,
            'two_m_q': two_m_q,
        })

    return edge_labels


# ---------------------------------------------------------------------------
# Block-diagonalization analysis
# ---------------------------------------------------------------------------

def analyze_block_diagonal(L1: np.ndarray, edge_labels: list,
                            sector_key_fn, sector_name: str):
    """Check whether L1 is block-diagonal under a given sector labeling.

    Parameters
    ----------
    L1 : ndarray (E x E)
    edge_labels : list of dicts with edge quantum numbers
    sector_key_fn : callable, edge_label_dict -> hashable sector key
    sector_name : str, description of the sector decomposition

    Returns
    -------
    dict with analysis results
    """
    E = L1.shape[0]

    # Assign each edge to a sector
    edge_sectors = [sector_key_fn(el) for el in edge_labels]
    unique_sectors = sorted(set(edge_sectors), key=str)

    # Build sector index mapping
    sector_indices = {}
    for sec in unique_sectors:
        idx = [i for i, s in enumerate(edge_sectors) if s == sec]
        sector_indices[str(sec)] = idx

    # Check block diagonality: are there nonzero entries L1[e_i, e_j]
    # where e_i and e_j are in different sectors?
    total_frob_sq = np.sum(L1**2)
    diag_block_frob_sq = 0.0
    cross_block_frob_sq = 0.0

    # Count nonzero cross-sector entries
    n_cross_nonzero = 0
    n_diag_nonzero = 0
    cross_entries = []

    for i in range(E):
        for j in range(E):
            val = L1[i, j]
            if abs(val) < 1e-14:
                continue
            si = edge_sectors[i]
            sj = edge_sectors[j]
            if si == sj:
                diag_block_frob_sq += val**2
                n_diag_nonzero += 1
            else:
                cross_block_frob_sq += val**2
                n_cross_nonzero += 1
                if len(cross_entries) < 20:  # limit output
                    cross_entries.append({
                        'i': i, 'j': j,
                        'sector_i': str(si), 'sector_j': str(sj),
                        'value': float(val),
                    })

    is_block_diagonal = (n_cross_nonzero == 0)

    cross_fraction = cross_block_frob_sq / total_frob_sq if total_frob_sq > 0 else 0.0
    diag_fraction = diag_block_frob_sq / total_frob_sq if total_frob_sq > 0 else 0.0

    # Sector dimensions
    sector_dims = {str(sec): len(idx) for sec, idx in sector_indices.items()}

    return {
        'sector_name': sector_name,
        'is_block_diagonal': bool(is_block_diagonal),
        'n_sectors': len(unique_sectors),
        'sector_dims': sector_dims,
        'total_frobenius_sq': float(total_frob_sq),
        'diagonal_block_frobenius_sq': float(diag_block_frob_sq),
        'cross_block_frobenius_sq': float(cross_block_frob_sq),
        'cross_fraction': float(cross_fraction),
        'diag_fraction': float(diag_fraction),
        'n_cross_nonzero': n_cross_nonzero,
        'n_diag_nonzero': n_diag_nonzero,
        'sample_cross_entries': cross_entries,
    }


def analyze_fock_graph(n_max: int) -> dict:
    """Full block-diagonal analysis for the scalar Fock graph at n_max."""
    print(f"\n{'='*60}")
    print(f"Scalar Fock graph, n_max = {n_max}")
    print(f"{'='*60}")

    states, edges, B, L1, V, E = build_fock_graph_numpy(n_max)
    edge_labels = label_fock_edges(states, edges)

    print(f"  V = {V}, E = {E}")

    # Print edge classification summary
    q_vals = [el['q'] for el in edge_labels]
    mq_vals = [el['m_q'] for el in edge_labels]
    types = [el['edge_type'] for el in edge_labels]
    print(f"  Edge types: {dict((t, types.count(t)) for t in set(types))}")
    print(f"  q values: {dict((q, q_vals.count(q)) for q in sorted(set(q_vals)))}")
    print(f"  m_q values: {dict((m, mq_vals.count(m)) for m in sorted(set(mq_vals)))}")

    results = {
        'graph_type': 'fock_scalar',
        'n_max': n_max,
        'V': V,
        'E': E,
        'edge_type_counts': dict((t, types.count(t)) for t in sorted(set(types))),
        'q_counts': dict((q, q_vals.count(q)) for q in sorted(set(q_vals))),
        'm_q_counts': dict((m, mq_vals.count(m)) for m in sorted(set(mq_vals))),
    }

    # Analysis 1: block-diagonal by (q, m_q)
    print(f"\n  Analysis 1: block-diagonal by (q, m_q)?")
    res_qmq = analyze_block_diagonal(
        L1, edge_labels,
        sector_key_fn=lambda el: (el['q'], el['m_q']),
        sector_name="(q, m_q)"
    )
    print(f"    Block-diagonal: {res_qmq['is_block_diagonal']}")
    print(f"    Sectors: {res_qmq['n_sectors']}, dims: {res_qmq['sector_dims']}")
    print(f"    Cross-block fraction: {res_qmq['cross_fraction']:.6f}")
    if not res_qmq['is_block_diagonal']:
        print(f"    Cross-block nonzero entries: {res_qmq['n_cross_nonzero']}")
        for ce in res_qmq['sample_cross_entries'][:5]:
            print(f"      L1[{ce['i']},{ce['j']}] = {ce['value']:.4f}, "
                  f"sectors {ce['sector_i']} vs {ce['sector_j']}")
    results['by_q_mq'] = res_qmq

    # Analysis 2: block-diagonal by q alone
    print(f"\n  Analysis 2: block-diagonal by q alone?")
    res_q = analyze_block_diagonal(
        L1, edge_labels,
        sector_key_fn=lambda el: el['q'],
        sector_name="q"
    )
    print(f"    Block-diagonal: {res_q['is_block_diagonal']}")
    print(f"    Sectors: {res_q['n_sectors']}, dims: {res_q['sector_dims']}")
    print(f"    Cross-block fraction: {res_q['cross_fraction']:.6f}")
    if not res_q['is_block_diagonal']:
        print(f"    Cross-block nonzero entries: {res_q['n_cross_nonzero']}")
    results['by_q'] = res_q

    # Analysis 3: block-diagonal by m_q alone
    print(f"\n  Analysis 3: block-diagonal by m_q alone?")
    res_mq = analyze_block_diagonal(
        L1, edge_labels,
        sector_key_fn=lambda el: el['m_q'],
        sector_name="m_q"
    )
    print(f"    Block-diagonal: {res_mq['is_block_diagonal']}")
    print(f"    Sectors: {res_mq['n_sectors']}, dims: {res_mq['sector_dims']}")
    print(f"    Cross-block fraction: {res_mq['cross_fraction']:.6f}")
    if not res_mq['is_block_diagonal']:
        print(f"    Cross-block nonzero entries: {res_mq['n_cross_nonzero']}")
    results['by_mq'] = res_mq

    # Analysis 4: block-diagonal by edge_type (T vs L)
    print(f"\n  Analysis 4: block-diagonal by edge type (T vs L)?")
    res_type = analyze_block_diagonal(
        L1, edge_labels,
        sector_key_fn=lambda el: el['edge_type'],
        sector_name="edge_type"
    )
    print(f"    Block-diagonal: {res_type['is_block_diagonal']}")
    print(f"    Sectors: {res_type['n_sectors']}, dims: {res_type['sector_dims']}")
    print(f"    Cross-block fraction: {res_type['cross_fraction']:.6f}")
    results['by_edge_type'] = res_type

    # Analysis 5: block-diagonal by l-sector of source node
    # (the Fock graph splits into l-sectors since T±/L± preserve l)
    print(f"\n  Analysis 5: block-diagonal by l-sector?")
    def l_sector(el):
        s1 = el['states'][0]
        return s1[1]  # l of source node
    res_l = analyze_block_diagonal(
        L1, edge_labels,
        sector_key_fn=l_sector,
        sector_name="l_sector"
    )
    print(f"    Block-diagonal: {res_l['is_block_diagonal']}")
    print(f"    Sectors: {res_l['n_sectors']}, dims: {res_l['sector_dims']}")
    print(f"    Cross-block fraction: {res_l['cross_fraction']:.6f}")
    results['by_l_sector'] = res_l

    # Also print the L1 eigenvalues for reference
    evals = np.sort(np.linalg.eigvalsh(L1))
    n_zero = np.sum(np.abs(evals) < 1e-10)
    print(f"\n  L1 spectrum: {n_zero} zero eigenvalues, "
          f"nonzero: {np.round(evals[evals > 0.5], 4)[:10]}...")
    results['L1_n_zero'] = int(n_zero)
    results['L1_nonzero_evals'] = sorted([float(e) for e in evals if e > 0.5])

    return results


def analyze_dirac_graph(n_max: int, rule: str = "B") -> dict:
    """Full block-diagonal analysis for the Dirac graph."""
    print(f"\n{'='*60}")
    print(f"Dirac graph Rule {rule}, n_max = {n_max}")
    print(f"{'='*60}")

    labels, edges, B, L1, V, E = build_dirac_graph_numpy(n_max, rule)
    edge_labels = label_dirac_edges(labels, edges, rule)

    print(f"  V = {V}, E = {E}")

    # Edge classification summary
    q_vals = [el['q'] for el in edge_labels]
    two_mq_vals = [el['two_m_q'] for el in edge_labels]
    types = [el['edge_type'] for el in edge_labels]
    print(f"  Edge types: {dict((t, types.count(t)) for t in sorted(set(types)))}")
    print(f"  q values: {dict((q, q_vals.count(q)) for q in sorted(set(q_vals)))}")
    print(f"  2*m_q values: {dict((m, two_mq_vals.count(m)) for m in sorted(set(two_mq_vals)))}")

    results = {
        'graph_type': f'dirac_rule_{rule}',
        'n_max': n_max,
        'V': V,
        'E': E,
        'edge_type_counts': dict((t, types.count(t)) for t in sorted(set(types))),
        'q_counts': dict((q, q_vals.count(q)) for q in sorted(set(q_vals))),
        'two_m_q_counts': dict((m, two_mq_vals.count(m)) for m in sorted(set(two_mq_vals))),
    }

    # Analysis 1: block-diagonal by (q, 2*m_q)
    print(f"\n  Analysis 1: block-diagonal by (q, 2*m_q)?")
    res_qmq = analyze_block_diagonal(
        L1, edge_labels,
        sector_key_fn=lambda el: (el['q'], el['two_m_q']),
        sector_name="(q, 2*m_q)"
    )
    print(f"    Block-diagonal: {res_qmq['is_block_diagonal']}")
    print(f"    Sectors: {res_qmq['n_sectors']}, dims: {res_qmq['sector_dims']}")
    print(f"    Cross-block fraction: {res_qmq['cross_fraction']:.6f}")
    if not res_qmq['is_block_diagonal']:
        print(f"    Cross-block nonzero entries: {res_qmq['n_cross_nonzero']}")
        for ce in res_qmq['sample_cross_entries'][:5]:
            print(f"      L1[{ce['i']},{ce['j']}] = {ce['value']:.4f}, "
                  f"sectors {ce['sector_i']} vs {ce['sector_j']}")
    results['by_q_2mq'] = res_qmq

    # Analysis 2: block-diagonal by q alone
    print(f"\n  Analysis 2: block-diagonal by q alone?")
    res_q = analyze_block_diagonal(
        L1, edge_labels,
        sector_key_fn=lambda el: el['q'],
        sector_name="q"
    )
    print(f"    Block-diagonal: {res_q['is_block_diagonal']}")
    print(f"    Sectors: {res_q['n_sectors']}, dims: {res_q['sector_dims']}")
    print(f"    Cross-block fraction: {res_q['cross_fraction']:.6f}")
    results['by_q'] = res_q

    # Analysis 3: block-diagonal by 2*m_q alone
    print(f"\n  Analysis 3: block-diagonal by 2*m_q alone?")
    res_mq = analyze_block_diagonal(
        L1, edge_labels,
        sector_key_fn=lambda el: el['two_m_q'],
        sector_name="2*m_q"
    )
    print(f"    Block-diagonal: {res_mq['is_block_diagonal']}")
    print(f"    Sectors: {res_mq['n_sectors']}, dims: {res_mq['sector_dims']}")
    print(f"    Cross-block fraction: {res_mq['cross_fraction']:.6f}")
    if not res_mq['is_block_diagonal']:
        print(f"    Cross-block nonzero entries: {res_mq['n_cross_nonzero']}")
    results['by_2mq'] = res_mq

    # Analysis 4: block-diagonal by edge type (T vs L)
    print(f"\n  Analysis 4: block-diagonal by edge type (T vs L)?")
    res_type = analyze_block_diagonal(
        L1, edge_labels,
        sector_key_fn=lambda el: el['edge_type'],
        sector_name="edge_type"
    )
    print(f"    Block-diagonal: {res_type['is_block_diagonal']}")
    print(f"    Sectors: {res_type['n_sectors']}, dims: {res_type['sector_dims']}")
    print(f"    Cross-block fraction: {res_type['cross_fraction']:.6f}")
    results['by_edge_type'] = res_type

    # Analysis 5: block-diagonal by Delta-l
    print(f"\n  Analysis 5: block-diagonal by Delta-l?")
    res_dl = analyze_block_diagonal(
        L1, edge_labels,
        sector_key_fn=lambda el: el['dl'],
        sector_name="dl"
    )
    print(f"    Block-diagonal: {res_dl['is_block_diagonal']}")
    print(f"    Sectors: {res_dl['n_sectors']}, dims: {res_dl['sector_dims']}")
    print(f"    Cross-block fraction: {res_dl['cross_fraction']:.6f}")
    results['by_dl'] = res_dl

    # L1 eigenvalues
    evals = np.sort(np.linalg.eigvalsh(L1))
    n_zero = np.sum(np.abs(evals) < 1e-10)
    print(f"\n  L1 spectrum: {n_zero} zero eigenvalues, "
          f"nonzero: {np.round(evals[evals > 0.5], 4)[:10]}...")
    results['L1_n_zero'] = int(n_zero)
    results['L1_nonzero_evals'] = sorted([float(e) for e in evals if e > 0.5])

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    all_results = {}

    # Scalar Fock graph at n_max = 2 and 3
    for n_max in [2, 3]:
        key = f'fock_nmax{n_max}'
        all_results[key] = analyze_fock_graph(n_max)

    # Dirac graph Rule B at n_max = 2 and 3
    for n_max in [2, 3]:
        key = f'dirac_B_nmax{n_max}'
        all_results[key] = analyze_dirac_graph(n_max, rule="B")

    # Dirac graph Rule A at n_max = 2 and 3
    for n_max in [2, 3]:
        key = f'dirac_A_nmax{n_max}'
        all_results[key] = analyze_dirac_graph(n_max, rule="A")

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")

    for key, res in all_results.items():
        print(f"\n{key}:")
        graph_type = res['graph_type']
        nmax = res['n_max']

        # Find the best block-diagonalization
        analyses = []
        for akey in ['by_q_mq', 'by_q_2mq', 'by_q', 'by_mq', 'by_2mq',
                      'by_edge_type', 'by_l_sector', 'by_dl']:
            if akey in res:
                a = res[akey]
                analyses.append((akey, a['is_block_diagonal'],
                                  a['cross_fraction'], a['n_sectors']))

        for akey, is_bd, cross_frac, n_sec in analyses:
            status = "YES" if is_bd else f"NO (cross={cross_frac:.4f})"
            print(f"  {akey:20s}: {status:30s} [{n_sec} sectors]")

    # Save results
    output_path = project_root / "debug" / "data" / "l1_photon_block_diag.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert any non-serializable types
    def clean_for_json(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, tuple):
            return list(obj)
        return obj

    def deep_clean(d):
        if isinstance(d, dict):
            return {str(k): deep_clean(v) for k, v in d.items()}
        if isinstance(d, list):
            return [deep_clean(v) for v in d]
        return clean_for_json(d)

    with output_path.open("w") as f:
        json.dump(deep_clean(all_results), f, indent=2)

    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
