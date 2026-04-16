"""
D6-A/B: Enumerate three Z_2 symmetries on the S^3 graph and test Fock projection compatibility.

Track D6 of the Dirac-on-S^3 program.

Three candidate Z_2 symmetries:
  sigma: m -> -m (azimuthal reflection / complex conjugation of Y_l^m)
  P:     kappa -> -kappa (Dirac kappa-parity, spinor-only)
  C:     chirality gamma^5 flip (spinor-only, no scalar action)

Author: GeoVac D6 sub-agent
Date: April 2026
"""

import json
import numpy as np
from fractions import Fraction
from collections import OrderedDict


# ============================================================================
# D6-A.1: Build the n_max=3 scalar graph explicitly
# ============================================================================

def build_scalar_graph(n_max: int = 3):
    """Build the scalar (n,l,m) graph. Returns nodes, edges, adjacency matrix."""
    # Generate states in canonical order: n=1..n_max, l=0..n-1, m=-l..l
    nodes = []
    node_index = {}
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                node_index[(n, l, m)] = len(nodes)
                nodes.append((n, l, m))

    N = len(nodes)
    adj = np.zeros((N, N), dtype=int)
    edges = []

    for i, (n, l, m) in enumerate(nodes):
        # Angular: m <-> m+1
        if m < l:
            target = (n, l, m + 1)
            j = node_index[target]
            if adj[i, j] == 0:
                adj[i, j] = 1
                adj[j, i] = 1
                edges.append((i, j))

        # Radial: n <-> n+1 (same l, m)
        if n < n_max:
            target = (n + 1, l, m)
            if target in node_index:
                j = node_index[target]
                if adj[i, j] == 0:
                    adj[i, j] = 1
                    adj[j, i] = 1
                    edges.append((i, j))

    return nodes, edges, adj, node_index


# ============================================================================
# D6-A.2: sigma (m-reflection) Z_2 action on scalar graph
# ============================================================================

def sigma_permutation(nodes, node_index):
    """Build the permutation matrix P_sigma for m -> -m."""
    N = len(nodes)
    P = np.zeros((N, N), dtype=int)
    for i, (n, l, m) in enumerate(nodes):
        target = (n, l, -m)
        j = node_index[target]
        P[i, j] = 1
    return P


def analyze_sigma(nodes, node_index, adj):
    """Full analysis of the sigma (m -> -m) symmetry."""
    P = sigma_permutation(nodes, node_index)
    N = len(nodes)

    # Check P^2 = I
    P_sq = P @ P
    is_involution = np.array_equal(P_sq, np.eye(N, dtype=int))

    # Check commutation with adjacency: P A P^T = A
    P_A_Pt = P @ adj @ P.T
    commutes_with_adj = np.array_equal(P_A_Pt, adj)

    # Fixed points (m=0 states) and 2-cycles
    fixed_points = []
    two_cycles = []
    visited = set()
    for i, (n, l, m) in enumerate(nodes):
        if i in visited:
            continue
        target = (n, l, -m)
        j = node_index[target]
        if i == j:
            fixed_points.append(nodes[i])
            visited.add(i)
        else:
            two_cycles.append((nodes[i], nodes[j]))
            visited.add(i)
            visited.add(j)

    # Eigenvalues of P_sigma (should be +1 and -1)
    eigenvalues = np.linalg.eigvalsh(P.astype(float))
    n_plus = int(np.sum(np.abs(eigenvalues - 1.0) < 1e-10))
    n_minus = int(np.sum(np.abs(eigenvalues + 1.0) < 1e-10))

    return {
        "is_involution": bool(is_involution),
        "commutes_with_adj": bool(commutes_with_adj),
        "fixed_points": fixed_points,
        "n_fixed": len(fixed_points),
        "two_cycles": two_cycles,
        "n_two_cycles": len(two_cycles),
        "eigenvalue_decomposition": {"n_plus1": n_plus, "n_minus1": n_minus},
        "permutation_matrix": P.tolist(),
    }


# ============================================================================
# D6-A.3: P (kappa-parity) Z_2 action on Dirac labels
# ============================================================================

def build_dirac_labels(n_max: int = 3):
    """Enumerate all Dirac labels (n_fock, kappa, two_m_j) at n_max.

    For each n_fock = 1..n_max:
      For each l = 0..n_fock-1:
        kappa = -(l+1) gives j = l+1/2
        kappa = +l (only if l >= 1) gives j = l-1/2
      For each kappa: two_m_j in range(-2j, 2j+1, step=2)
    """
    labels = []
    for n_fock in range(1, n_max + 1):
        for l in range(n_fock):
            kappa_options = [-(l + 1)]
            if l >= 1:
                kappa_options.append(l)
            for kappa in kappa_options:
                j2 = 2 * abs(kappa) - 1  # 2*j
                for two_mj in range(-j2, j2 + 1, 2):
                    labels.append((n_fock, kappa, two_mj))
    return labels


def kappa_to_l(kappa):
    """l from kappa."""
    if kappa < 0:
        return -kappa - 1
    return kappa


def kappa_to_j2(kappa):
    """2*j from kappa."""
    return 2 * abs(kappa) - 1


def analyze_kappa_parity(dirac_labels):
    """Analyze the kappa -> -kappa map on Dirac labels.

    Key question: does kappa -> -kappa map the label set to itself?
    """
    label_set = set(dirac_labels)
    label_index = {lab: i for i, lab in enumerate(dirac_labels)}

    N = len(dirac_labels)

    # Try kappa -> -kappa: (n, kappa, two_mj) -> (n, -kappa, two_mj)
    # This changes l and j:
    #   kappa < 0: l = -kappa-1, j = l+1/2. -kappa > 0: l' = -kappa = l+1, j' = l+1-1/2 = l+1/2.
    #     So j is preserved but l changes: l -> l+1.
    #   kappa > 0: l = kappa, j = l-1/2. -kappa < 0: l' = kappa-1 = l-1, j' = l-1+1/2 = l-1/2.
    #     So j is preserved but l changes: l -> l-1.
    # The two_m_j range depends on j = |kappa|-1/2. Since |-kappa| = |kappa|, j is unchanged.
    # So two_mj is always valid after the flip.
    # BUT: l' must satisfy l' < n_fock. Check this.

    maps_to_self = True
    kappa_perm = np.zeros((N, N), dtype=int)
    unmapped = []
    mapped_pairs = []
    fixed = []

    for i, (n, kappa, two_mj) in enumerate(dirac_labels):
        target = (n, -kappa, two_mj)
        if target in label_set:
            j_idx = label_index[target]
            kappa_perm[i, j_idx] = 1
            if i == j_idx:
                fixed.append(dirac_labels[i])
            elif i < j_idx:
                mapped_pairs.append((dirac_labels[i], dirac_labels[j_idx]))
        else:
            maps_to_self = False
            unmapped.append(dirac_labels[i])

    # Check if it's a valid permutation
    is_permutation = (np.sum(kappa_perm, axis=1) == 1).all() and (np.sum(kappa_perm, axis=0) == 1).all()
    is_involution = np.array_equal(kappa_perm @ kappa_perm, np.eye(N, dtype=int)) if is_permutation else False

    # Detailed analysis of what kappa -> -kappa does
    kappa_analysis = []
    for n, kappa, two_mj in dirac_labels:
        l_orig = kappa_to_l(kappa)
        j2_orig = kappa_to_j2(kappa)
        new_kappa = -kappa
        l_new = kappa_to_l(new_kappa) if new_kappa != 0 else None
        j2_new = kappa_to_j2(new_kappa) if new_kappa != 0 else None
        kappa_analysis.append({
            "label": (n, kappa, two_mj),
            "l": l_orig, "j2": j2_orig,
            "target_kappa": new_kappa,
            "target_l": l_new, "target_j2": j2_new,
            "target_in_set": (n, new_kappa, two_mj) in label_set,
            "l_change": l_new - l_orig if l_new is not None else None,
        })

    return {
        "n_dirac_labels": N,
        "maps_to_self": bool(maps_to_self),
        "is_permutation": bool(is_permutation),
        "is_involution": bool(is_involution),
        "n_unmapped": len(unmapped),
        "unmapped_labels": [str(u) for u in unmapped],
        "n_fixed_points": len(fixed),
        "fixed_points": [str(f) for f in fixed],
        "n_two_cycles": len(mapped_pairs),
        "kappa_analysis_summary": kappa_analysis[:10],  # first 10 for readability
    }


# ============================================================================
# D6-A.4: C (chirality) -- scalar graph has no action
# ============================================================================

def analyze_chirality(n_max: int = 3):
    """Document that chirality has no action on the scalar graph.

    On the Dirac-on-S^3 spectrum, each CH level n has eigenvalues
    +|lambda_n| and -|lambda_n|, corresponding to the two Weyl sectors.
    The chirality flip swaps these sectors.

    On the scalar graph (14 nodes at n_max=3), there is no spinor
    structure and chirality is undefined.
    """
    # Count states in each sector for the Dirac case
    dirac_counts = {}
    for n_ch in range(n_max):
        n_fock = n_ch + 1
        g_weyl = (n_ch + 1) * (n_ch + 2)  # per chirality sector
        g_dirac = 2 * g_weyl
        dirac_counts[n_fock] = {
            "g_weyl_plus": g_weyl,
            "g_weyl_minus": g_weyl,
            "g_dirac_total": g_dirac,
        }

    total_dirac = sum(v["g_dirac_total"] for v in dirac_counts.values())
    total_weyl = sum(v["g_weyl_plus"] for v in dirac_counts.values())

    return {
        "acts_on_scalar_graph": False,
        "reason": "Chirality (gamma^5) acts on spinor sectors. The scalar graph has no spinor structure.",
        "dirac_level_counts": dirac_counts,
        "total_dirac_states_nmax3": total_dirac,
        "total_weyl_states_nmax3": total_weyl,
        "chirality_splits_as": "Each Weyl sector has identical state count. "
                               "At n_max=3 (CH n=0,1,2): total 40 Dirac = 20+20 Weyl.",
    }


# ============================================================================
# D6-B.5: Fock projection compatibility of sigma
# ============================================================================

def fock_compatibility_sigma(nodes, node_index, adj):
    """Check if sigma commutes with the adjacency matrix (graph automorphism).

    If P_sigma A = A P_sigma, then sigma is a graph automorphism and
    automatically compatible with the Fock projection (eigenvalues preserved).
    """
    P = sigma_permutation(nodes, node_index)
    PA = P @ adj
    AP = adj @ P
    commutes = np.array_equal(PA, AP)

    # Also check: does sigma preserve the degree matrix D?
    degrees = np.sum(adj, axis=1)
    sigma_degrees = P @ degrees
    preserves_degree = np.array_equal(degrees, sigma_degrees)

    # Check: does sigma preserve the full graph Laplacian L = D - A?
    D_mat = np.diag(degrees)
    L = D_mat - adj
    P_L_Pt = P @ L @ P.T
    preserves_laplacian = np.array_equal(P_L_Pt, L)

    return {
        "sigma_commutes_with_A": bool(commutes),
        "sigma_preserves_degree": bool(preserves_degree),
        "sigma_preserves_laplacian": bool(preserves_laplacian),
        "is_graph_automorphism": bool(commutes),
        "fock_compatible": bool(commutes),
        "reason": "sigma (m->-m) is a graph automorphism iff P_sigma A = A P_sigma. "
                  "Graph automorphisms commute with the Laplacian and hence preserve "
                  "the Fock eigenvalue spectrum. Fock projection maps eigenvalues to "
                  "hydrogenic energies, so sigma-compatibility is automatic."
    }


# ============================================================================
# D6-B.6: Dirac accidental degeneracy and kappa-parity
# ============================================================================

def analyze_accidental_degeneracy(n_max: int = 3):
    """Identify Dirac accidental degeneracy pairs and their relation to kappa -> -kappa.

    The Dirac formula E(n,j) depends on n and j only, not on l separately.
    States with same (n, j) but different l are accidentally degenerate.

    In kappa language: |kappa| = j + 1/2, so same j means same |kappa|.
    The two l values for a given j are:
      kappa = -(j+1/2) = -(l+1) with l = j-1/2  (kappa < 0)
      kappa = +(j+1/2) = +l'    with l' = j+1/2  (kappa > 0)

    So the accidental degeneracy pair swaps kappa -> -kappa!
    Specifically: kappa and -kappa have the same |kappa|, hence same j,
    and the l values differ by 1 (l_new = l_old + 1 when kappa < 0 -> kappa > 0).

    The map IS kappa -> -kappa, not kappa -> -kappa - sgn(kappa).
    """
    pairs = []
    for n_fock in range(1, n_max + 1):
        for l in range(n_fock):
            kappa_neg = -(l + 1)  # j = l + 1/2
            kappa_pos = l + 1     # j = (l+1) - 1/2 = l + 1/2 (same j!)
            # But kappa_pos = l+1 requires l' = kappa_pos = l+1 < n_fock
            if l + 1 < n_fock:
                # Verify: kappa_pos = -(kappa_neg) = l+1, and kappa_to_l(l+1) = l+1
                # kappa_to_j(kappa_neg) = |-(l+1)| - 1/2 = l + 1/2
                # kappa_to_j(kappa_pos) = |l+1| - 1/2 = l + 1/2  ✓ same j
                j_val = Fraction(2 * (l + 1) - 1, 2)
                l_neg = l   # l from kappa_neg
                l_pos = l + 1  # l from kappa_pos

                # Standard spectroscopic notation
                l_names = "spdfghiklm"
                l_name_neg = l_names[l_neg] if l_neg < len(l_names) else f"l={l_neg}"
                l_name_pos = l_names[l_pos] if l_pos < len(l_names) else f"l={l_pos}"

                pairs.append({
                    "n_fock": n_fock,
                    "j": str(j_val),
                    "kappa_1": kappa_neg,
                    "l_1": l_neg,
                    "name_1": f"{n_fock}{l_name_neg}_{j_val}",
                    "kappa_2": kappa_pos,
                    "l_2": l_pos,
                    "name_2": f"{n_fock}{l_name_pos}_{j_val}",
                    "is_kappa_negation": kappa_pos == -kappa_neg,
                })

    # Extend to n=4 for the full 6 pairs
    pairs_n4 = []
    for n_fock in range(1, 5):
        for l in range(n_fock):
            kappa_neg = -(l + 1)
            kappa_pos = l + 1
            if l + 1 < n_fock:
                j_val = Fraction(2 * (l + 1) - 1, 2)
                l_names = "spdfghiklm"
                l_neg = l
                l_pos = l + 1
                l_name_neg = l_names[l_neg] if l_neg < len(l_names) else f"l={l_neg}"
                l_name_pos = l_names[l_pos] if l_pos < len(l_names) else f"l={l_pos}"
                pairs_n4.append({
                    "n_fock": n_fock,
                    "j": str(j_val),
                    "name_1": f"{n_fock}{l_name_neg}_{j_val}",
                    "name_2": f"{n_fock}{l_name_pos}_{j_val}",
                    "is_kappa_negation": True,
                })

    return {
        "pairs_within_nmax3": pairs,
        "n_pairs_nmax3": len(pairs),
        "pairs_through_n4": pairs_n4,
        "n_pairs_n4": len(pairs_n4),
        "structural_result": (
            "The Dirac accidental degeneracy map IS exactly kappa -> -kappa. "
            "For each degenerate pair (n, kappa) and (n, -kappa): "
            "|kappa| is the same so j is the same, but l differs by 1. "
            "This is NOT kappa -> -kappa - sgn(kappa); it is the simple negation."
        ),
    }


# ============================================================================
# D6-B Key question: relation between sigma and P
# ============================================================================

def analyze_sigma_vs_kappa_parity():
    """Analyze whether sigma (m -> -m) and P (kappa -> -kappa) are related.

    sigma acts on the scalar graph: (n, l, m) -> (n, l, -m)
    P acts on the Dirac labels: (n, kappa, two_mj) -> (n, -kappa, two_mj)

    These are operations on DIFFERENT spaces. The question is whether
    there is a natural embedding where they become the same operation.
    """
    analysis = {
        "sigma_space": "scalar (n, l, m) — 14 nodes at n_max=3",
        "P_space": "Dirac (n, kappa, two_mj) — varies with n_max",
        "sigma_action": "m -> -m (preserves n, l)",
        "P_action": "kappa -> -kappa (preserves n, two_mj; changes l by ±1, preserves j)",

        "are_same_operation": False,
        "reason": (
            "sigma and P act on different quantum numbers and different spaces. "
            "sigma flips m (magnetic quantum number) while preserving l. "
            "P flips kappa (spin-orbit alignment) while preserving j but changing l. "
            "In the scalar graph, P has no meaning (no spin/kappa structure). "
            "In the Dirac graph, sigma generalizes to m_j -> -m_j (time reversal), "
            "which is independent of the kappa -> -kappa map."
        ),

        "embedding_analysis": (
            "The natural two-to-one lift from scalar (n,l,m) to Dirac (n,kappa,m_j) is: "
            "(n, l, m) lifts to (n, -(l+1), 2m±1) and (n, +l, 2m±1) for l>=1. "
            "Under this lift: "
            "  sigma (m -> -m) lifts to m_j -> -m_j (time reversal on spinors). "
            "  P (kappa -> -kappa) has no scalar pre-image (it mixes different l). "
            "These are genuinely independent Z_2 symmetries."
        ),

        "group_structure": (
            "On the full Dirac label space, the three Z_2s generate: "
            "  sigma_spinor: m_j -> -m_j (time reversal) "
            "  P: kappa -> -kappa (accidental degeneracy / parity) "
            "  C: chirality flip (Weyl sector swap) "
            "sigma_spinor and P commute (act on different quantum numbers). "
            "C acts on a separate index (chirality). "
            "Together they generate Z_2 x Z_2 x Z_2 on the full Dirac-on-S^3 space. "
            "On the scalar graph, only sigma (= sigma_spinor restricted to m) acts."
        ),
    }
    return analysis


# ============================================================================
# Main: run all analyses and save
# ============================================================================

def main():
    print("=" * 70)
    print("D6-A/B: Z_2 symmetry enumeration on the S^3 graph")
    print("=" * 70)

    # D6-A.1: Build graph
    print("\n--- D6-A.1: n_max=3 scalar graph ---")
    nodes, edges, adj, node_index = build_scalar_graph(n_max=3)
    print(f"Nodes: {len(nodes)}")
    print(f"Edges: {len(edges)}")
    for i, (n, l, m) in enumerate(nodes):
        print(f"  [{i:2d}] (n={n}, l={l}, m={m:+d})")
    print(f"\nEdge list:")
    for i, j in edges:
        print(f"  {nodes[i]} -- {nodes[j]}")

    # Verify expected counts
    assert len(nodes) == 14, f"Expected 14 nodes, got {len(nodes)}"
    assert len(edges) == 13, f"Expected 13 edges, got {len(edges)}"

    # D6-A.2: sigma analysis
    print("\n--- D6-A.2: sigma (m -> -m) ---")
    sigma_result = analyze_sigma(nodes, node_index, adj)
    print(f"Is involution (P^2 = I): {sigma_result['is_involution']}")
    print(f"Commutes with adjacency: {sigma_result['commutes_with_adj']}")
    print(f"Fixed points (m=0): {sigma_result['n_fixed']}")
    for fp in sigma_result['fixed_points']:
        print(f"  {fp}")
    print(f"2-cycles: {sigma_result['n_two_cycles']}")
    for c in sigma_result['two_cycles']:
        print(f"  {c[0]} <-> {c[1]}")
    print(f"Eigenvalue decomposition: {sigma_result['eigenvalue_decomposition']}")

    # D6-A.3: kappa-parity analysis
    print("\n--- D6-A.3: P (kappa -> -kappa) on Dirac labels ---")
    dirac_labels = build_dirac_labels(n_max=3)
    print(f"Total Dirac labels at n_max=3: {len(dirac_labels)}")
    for i, (n, kappa, two_mj) in enumerate(dirac_labels):
        l = kappa_to_l(kappa)
        j2 = kappa_to_j2(kappa)
        print(f"  [{i:2d}] n={n}, kappa={kappa:+d}, 2*m_j={two_mj:+d}  "
              f"(l={l}, j={j2}/2)")

    kappa_result = analyze_kappa_parity(dirac_labels)
    print(f"\nMaps label set to itself: {kappa_result['maps_to_self']}")
    print(f"Is valid permutation: {kappa_result['is_permutation']}")
    print(f"Is involution: {kappa_result['is_involution']}")
    print(f"Unmapped labels: {kappa_result['n_unmapped']}")
    for u in kappa_result['unmapped_labels']:
        print(f"  {u}")
    print(f"Fixed points: {kappa_result['n_fixed_points']}")
    for f in kappa_result['fixed_points']:
        print(f"  {f}")
    print(f"2-cycles: {kappa_result['n_two_cycles']}")

    # D6-A.4: chirality
    print("\n--- D6-A.4: C (chirality) ---")
    chirality_result = analyze_chirality(n_max=3)
    print(f"Acts on scalar graph: {chirality_result['acts_on_scalar_graph']}")
    print(f"Reason: {chirality_result['reason']}")
    print(f"Chirality split: {chirality_result['chirality_splits_as']}")

    # D6-B.5: Fock compatibility
    print("\n--- D6-B.5: Fock projection compatibility ---")
    fock_result = fock_compatibility_sigma(nodes, node_index, adj)
    print(f"sigma commutes with A: {fock_result['sigma_commutes_with_A']}")
    print(f"sigma preserves degree: {fock_result['sigma_preserves_degree']}")
    print(f"sigma preserves Laplacian: {fock_result['sigma_preserves_laplacian']}")
    print(f"sigma is graph automorphism: {fock_result['is_graph_automorphism']}")
    print(f"sigma is Fock-compatible: {fock_result['fock_compatible']}")

    # D6-B.6: Accidental degeneracy
    print("\n--- D6-B.6: Dirac accidental degeneracy pairs ---")
    degen_result = analyze_accidental_degeneracy(n_max=3)
    print(f"Pairs within n_max=3: {degen_result['n_pairs_nmax3']}")
    for p in degen_result['pairs_within_nmax3']:
        print(f"  {p['name_1']} <-> {p['name_2']}  "
              f"(kappa {p['kappa_1']:+d} <-> {p['kappa_2']:+d}, "
              f"is negation: {p['is_kappa_negation']})")
    print(f"\nAll 6 pairs through n=4:")
    for p in degen_result['pairs_through_n4']:
        print(f"  {p['name_1']} <-> {p['name_2']}  (kappa negation: {p['is_kappa_negation']})")
    print(f"\n{degen_result['structural_result']}")

    # Key question: sigma vs P
    print("\n--- Key question: sigma vs kappa-parity ---")
    relation = analyze_sigma_vs_kappa_parity()
    print(f"Are same operation: {relation['are_same_operation']}")
    print(f"Reason: {relation['reason']}")
    print(f"\nEmbedding analysis: {relation['embedding_analysis']}")
    print(f"\nGroup structure: {relation['group_structure']}")

    # ========================================================================
    # Save results
    # ========================================================================
    results = {
        "graph": {
            "n_max": 3,
            "n_nodes": len(nodes),
            "n_edges": len(edges),
            "nodes": [list(n) for n in nodes],
            "edges": [[int(i), int(j)] for i, j in edges],
            "adjacency_matrix": adj.tolist(),
        },
        "sigma_m_reflection": {
            "is_involution": sigma_result["is_involution"],
            "commutes_with_adj": sigma_result["commutes_with_adj"],
            "n_fixed_points": sigma_result["n_fixed"],
            "fixed_points": [list(fp) for fp in sigma_result["fixed_points"]],
            "n_two_cycles": sigma_result["n_two_cycles"],
            "two_cycles": [[list(c[0]), list(c[1])] for c in sigma_result["two_cycles"]],
            "eigenvalue_decomposition": sigma_result["eigenvalue_decomposition"],
        },
        "kappa_parity": {
            "n_dirac_labels": kappa_result["n_dirac_labels"],
            "maps_to_self": kappa_result["maps_to_self"],
            "is_permutation": kappa_result["is_permutation"],
            "is_involution": kappa_result["is_involution"],
            "n_unmapped": kappa_result["n_unmapped"],
            "unmapped_labels": kappa_result["unmapped_labels"],
            "n_fixed_points": kappa_result["n_fixed_points"],
            "n_two_cycles": kappa_result["n_two_cycles"],
        },
        "chirality": {
            "acts_on_scalar_graph": chirality_result["acts_on_scalar_graph"],
            "total_dirac_states": chirality_result["total_dirac_states_nmax3"],
            "total_weyl_states": chirality_result["total_weyl_states_nmax3"],
        },
        "fock_compatibility": {
            "sigma_commutes_with_A": fock_result["sigma_commutes_with_A"],
            "sigma_preserves_laplacian": fock_result["sigma_preserves_laplacian"],
            "sigma_is_graph_automorphism": fock_result["is_graph_automorphism"],
            "sigma_fock_compatible": fock_result["fock_compatible"],
        },
        "accidental_degeneracy": {
            "n_pairs_nmax3": degen_result["n_pairs_nmax3"],
            "pairs_nmax3": degen_result["pairs_within_nmax3"],
            "n_pairs_n4": degen_result["n_pairs_n4"],
            "is_kappa_negation": True,
            "structural_result": degen_result["structural_result"],
        },
        "sigma_vs_P_relation": {
            "are_same_operation": relation["are_same_operation"],
            "sigma_space": relation["sigma_space"],
            "P_space": relation["P_space"],
            "group_structure": relation["group_structure"],
        },
    }

    output_path = "C:/Users/jlout/Desktop/Project_Geometric/debug/data/d6_z2_enumeration.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
