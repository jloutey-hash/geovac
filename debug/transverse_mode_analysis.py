"""
Transverse photon mode analysis: can co-exact eigenvectors carry
angular momentum labels?

The 3 remaining selection rules (SO(4) channel count, Ward identity,
charge conjugation) require the photon to carry quantum numbers (L, M_L).
On the continuum S^3, co-exact 1-forms at eigenvalue n(n+2) carry
angular momentum q = n.

This script:
1. Diagonalizes d_1 d_1^T (co-exact Laplacian) and gets eigenvectors
2. Analyzes each eigenvector's distribution across shell-pair edges
3. Tests for rotational symmetry (m-quantum-number structure)
4. Checks whether continuum eigenvalue matching n(n+2) gives usable q labels
5. Builds a mode-resolved self-energy and tests SO(4) channel count
"""

import numpy as np
from pathlib import Path
import json
import time
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from geovac.lattice import GeometricLattice
from geovac.graph_qed_vertex import (
    build_projection_matrix, build_vertex_tensor,
    vertex_tensor_to_matrices
)
from geovac.graph_qed_photon import build_fock_graph
from geovac.su2_wilson_gauge import enumerate_plaquettes


def build_shell_edge_map(lat, edges):
    """Map each edge to its shell-pair (n_src, n_tgt).

    Returns dict: edge_index -> (n_low, n_high) and inverse map.
    """
    states = lat.states

    edge_shells = {}
    shell_pair_edges = {}  # (n_low, n_high) -> list of edge indices

    for e_idx, (i, j) in enumerate(edges):
        n_i = states[i][0]
        n_j = states[j][0]
        n_low, n_high = min(n_i, n_j), max(n_i, n_j)
        edge_shells[e_idx] = (n_low, n_high)

        key = (n_low, n_high)
        if key not in shell_pair_edges:
            shell_pair_edges[key] = []
        shell_pair_edges[key].append(e_idx)

    return edge_shells, shell_pair_edges, states


def analyze_coexact_modes(n_max):
    """Full analysis of co-exact eigenvectors."""
    print(f"\n{'='*70}")
    print(f"CO-EXACT MODE ANALYSIS at n_max={n_max}")
    print(f"{'='*70}")

    # Build graph
    lat = GeometricLattice(max_n=n_max, topological_weights=False)
    V = lat.num_states
    adj_sparse = lat.adjacency
    adj = np.array(adj_sparse.todense(), dtype=float)

    # Get edges
    rows, cols = adj_sparse.nonzero()
    edge_set = set()
    for r, c in zip(rows, cols):
        if r < c:
            edge_set.add((int(r), int(c)))
    edges = sorted(edge_set)
    E = len(edges)
    edge_index = {e: i for i, e in enumerate(edges)}

    # Incidence matrix
    B = np.zeros((V, E), dtype=float)
    for k, (i, j) in enumerate(edges):
        B[i, k] = 1.0
        B[j, k] = -1.0

    # Plaquettes and d_1
    plaquettes = enumerate_plaquettes(adj, max_length=8, both_orientations=False)
    F = len(plaquettes)

    d1 = np.zeros((E, F), dtype=float)
    for f_idx, plaq in enumerate(plaquettes):
        for oe in plaq:
            src, tgt = oe.source, oe.target
            if src < tgt:
                e_key = (src, tgt)
                sign = +1.0
            else:
                e_key = (tgt, src)
                sign = -1.0
            if e_key in edge_index:
                e_idx = edge_index[e_key]
                d1[e_idx, f_idx] += sign

    # Co-exact Laplacian
    L1_coexact = d1 @ d1.T
    evals, evecs = np.linalg.eigh(L1_coexact)

    # Identify nonzero modes
    tol = 1e-8
    nz_mask = evals > tol
    nz_evals = evals[nz_mask]
    nz_evecs = evecs[:, nz_mask]  # each column is an eigenvector
    n_modes = len(nz_evals)

    print(f"\n  Graph: V={V}, E={E}, F={F}")
    print(f"  Co-exact modes: {n_modes}")
    print(f"  Eigenvalues: {[f'{x:.4f}' for x in nz_evals]}")

    # Continuum comparison: co-exact 1-forms on S^3 have eigenvalues n(n+2)
    continuum_evals = [n*(n+2) for n in range(1, 10)]
    print(f"\n  Continuum n(n+2): {continuum_evals[:8]}")

    # Try to assign q labels by nearest continuum eigenvalue
    print(f"\n  Mode-to-continuum matching:")
    q_assignments = []
    for i, ev in enumerate(nz_evals):
        best_q = min(range(1, 10), key=lambda q: abs(q*(q+2) - ev))
        residual = ev - best_q*(best_q+2)
        q_assignments.append(best_q)
        print(f"    Mode {i}: eigenvalue {ev:.4f}, nearest q={best_q} "
              f"(continuum {best_q*(best_q+2)}), residual {residual:+.4f}")

    # Shell-pair distribution of each eigenvector
    edge_shells, shell_pair_edges, states = build_shell_edge_map(lat, edges)

    print(f"\n  Shell-pair edge counts:")
    for pair, eidxs in sorted(shell_pair_edges.items()):
        print(f"    shells ({pair[0]},{pair[1]}): {len(eidxs)} edges")

    print(f"\n  Eigenvector shell-pair weight distribution:")
    print(f"  {'Mode':<6} {'Eigenval':<10}", end="")
    for pair in sorted(shell_pair_edges.keys()):
        print(f" ({pair[0]},{pair[1]})", end="")
    print()

    mode_profiles = []
    for i in range(n_modes):
        vec = nz_evecs[:, i]
        profile = {}
        print(f"  {i:<6} {nz_evals[i]:<10.4f}", end="")
        for pair in sorted(shell_pair_edges.keys()):
            eidxs = shell_pair_edges[pair]
            weight = np.sum(vec[eidxs]**2)
            profile[pair] = weight
            print(f" {weight:>6.3f}", end="")
        print()
        mode_profiles.append(profile)

    # Key question: do modes localize on specific shell pairs?
    # A photon with q labels should connect shells (n, n') with |n-n'| related to q
    # On the Fock graph, edges only connect adjacent shells (Delta n = 1)
    # So ALL co-exact modes must live on nearest-neighbor shell pairs
    # The q label must come from WITHIN-shell structure (m quantum numbers)

    print(f"\n  === WITHIN-SHELL-PAIR STRUCTURE ===")
    print(f"  (Edges connect (n, n+1) pairs. q must be encoded in m-structure.)")

    # For each shell pair, decompose eigenvectors by the m-quantum numbers
    # of their endpoints
    for pair in sorted(shell_pair_edges.keys()):
        eidxs = shell_pair_edges[pair]
        if len(eidxs) < 2:
            continue

        n_low, n_high = pair
        print(f"\n  Shell pair ({n_low}, {n_high}): {len(eidxs)} edges")

        # Get the (l, m) labels for each endpoint
        for e_idx in eidxs[:min(5, len(eidxs))]:
            i, j = edges[e_idx]
            n_i, l_i, m_i = states[i]
            n_j, l_j, m_j = states[j]
            if n_i > n_j:
                n_i, l_i, m_i, n_j, l_j, m_j = n_j, l_j, m_j, n_i, l_i, m_i
            print(f"    e{e_idx}: ({n_i},{l_i},{m_i}) -- ({n_j},{l_j},{m_j})"
                  f"  Delta_l={l_j-l_i}, Delta_m={m_j-m_i}")

    # === MODE-RESOLVED SELF-ENERGY ===
    # Build vertex tensor
    print(f"\n  === MODE-RESOLVED SELF-ENERGY ===")
    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    fock_data = build_fock_graph(n_max)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats_sympy = vertex_tensor_to_matrices(entries, N_dirac, E_fock)
    V_mats = [np.array(vm.tolist(), dtype=float) for vm in V_mats_sympy]

    assert E_fock == E, f"Edge mismatch: {E_fock} vs {E}"

    # Compute self-energy from each individual co-exact mode
    print(f"\n  Per-mode self-energy contributions:")
    print(f"  {'Mode':<6} {'q_assign':<9} {'Tr(Sigma_k)':<14} {'GS block':<12} {'Diag frac':<10}")

    mode_sigmas = []
    for k in range(n_modes):
        vec = nz_evecs[:, k]
        ev = nz_evals[k]
        # Mode-k propagator: (1/lambda_k) * |v_k><v_k|
        G_k = (1.0 / ev) * np.outer(vec, vec)

        # Self-energy from this mode alone
        Sigma_k = np.zeros((N_dirac, N_dirac))
        for e1 in range(E_fock):
            for e2 in range(E_fock):
                g_ee = G_k[e1, e2]
                if abs(g_ee) < 1e-15:
                    continue
                Sigma_k += g_ee * (V_mats[e1] @ V_mats[e2].T)

        mode_sigmas.append(Sigma_k)

        # Diagnostics
        tr_k = np.trace(Sigma_k)
        # GS indices (n_fock=1, kappa=-1)
        gs_indices = [i for i, dl in enumerate(dirac_labels)
                      if dl.n_fock == 1 and dl.kappa == -1]
        gs_block = Sigma_k[np.ix_(gs_indices, gs_indices)]
        gs_trace = np.trace(gs_block)
        diag_frac = np.sum(np.diag(Sigma_k)**2) / np.sum(Sigma_k**2) if np.sum(Sigma_k**2) > 1e-20 else 0

        print(f"  {k:<6} q={q_assignments[k]:<6} {tr_k:<14.6f} {gs_trace:<12.6e} {diag_frac:<10.4f}")

    # Which Dirac labels does each mode couple?
    # Analyze the pattern: does mode k only couple specific (n1, n2) electron pairs?
    print(f"\n  === ELECTRON-STATE COUPLING PATTERN PER MODE ===")
    print(f"  (Which electron shell-pairs does each photon mode connect?)")

    for k in range(n_modes):
        Sigma_k = mode_sigmas[k]
        # For each pair of Dirac states (a, b), check if Sigma_k[a,b] is nonzero
        # Group by shell pair of the electrons
        shell_coupling = {}
        for a in range(N_dirac):
            n_a = dirac_labels[a].n_fock
            for b in range(N_dirac):
                n_b = dirac_labels[b].n_fock
                if abs(Sigma_k[a, b]) > 1e-10:
                    key = (min(n_a, n_b), max(n_a, n_b))
                    if key not in shell_coupling:
                        shell_coupling[key] = 0.0
                    shell_coupling[key] += Sigma_k[a, b]**2

        if shell_coupling:
            total = sum(shell_coupling.values())
            top_pairs = sorted(shell_coupling.items(), key=lambda x: -x[1])[:4]
            pairs_str = ", ".join(f"({p[0]},{p[1]}):{v/total:.2f}" for p, v in top_pairs)
            print(f"    Mode {k} (q={q_assignments[k]}): {pairs_str}")

    # === KEY TEST: Does q-resolved propagator enforce triangle inequality? ===
    # In the continuum, SO(4) channel count W(n1, n2, q) = 0 unless
    # |n1 - n2| <= q <= n1 + n2 - 2 (the triangle condition)
    #
    # If we assign q to each co-exact mode, does Sigma restricted to
    # modes with q violating the triangle inequality vanish?

    print(f"\n  === TRIANGLE INEQUALITY TEST ===")
    print(f"  Testing: does photon mode with assigned q only couple")
    print(f"  electron states (n1, n2) satisfying |n1-n2| <= q <= n1+n2-2?")

    for k in range(n_modes):
        q = q_assignments[k]
        Sigma_k = mode_sigmas[k]
        violations = 0
        total_nonzero = 0

        for a in range(N_dirac):
            n1 = dirac_labels[a].n_fock
            for b in range(a, N_dirac):
                n2 = dirac_labels[b].n_fock
                if abs(Sigma_k[a, b]) > 1e-10:
                    total_nonzero += 1
                    # Triangle condition (using Fock convention n_fock = n_CH + 1)
                    # Continuum uses n_CH; our n_fock = n_CH + 1
                    n1_ch = n1 - 1
                    n2_ch = n2 - 1
                    if n1_ch + n2_ch < q or abs(n1_ch - n2_ch) > q:
                        violations += 1

        viol_frac = violations / total_nonzero if total_nonzero > 0 else 0
        status = "PASS" if violations == 0 else f"FAIL ({violations}/{total_nonzero} = {viol_frac:.1%})"
        print(f"    Mode {k} (q={q}): {status}")

    # === WARD IDENTITY TEST ===
    # Ward identity: sum over photon modes of (vertex * propagator * vertex)
    # should relate to electron propagator.
    # Simple test: does Sigma_T commute with h1 (free Dirac operator)?
    # In continuum QED, [Sigma, H_0] has specific structure from Ward identity.

    # Build h1 (free Dirac on graph)
    # h1 diagonal = eigenvalues lambda_n = n + 3/2 (CH convention)
    h1 = np.zeros((N_dirac, N_dirac))
    for a in range(N_dirac):
        n_fock = dirac_labels[a].n_fock
        n_ch = n_fock - 1  # Camporesi-Higuchi convention
        h1[a, a] = n_ch + 1.5  # |lambda_n| = n + 3/2

    # Full transverse self-energy
    Sigma_T = sum(mode_sigmas)

    comm = Sigma_T @ h1 - h1 @ Sigma_T
    comm_norm = np.linalg.norm(comm)
    sigma_norm = np.linalg.norm(Sigma_T)

    print(f"\n  === WARD IDENTITY DIAGNOSTIC ===")
    print(f"  ||[Sigma_T, H_0]|| / ||Sigma_T|| = {comm_norm/sigma_norm:.6f}")
    print(f"  (Zero would indicate Sigma_T is diagonal in H_0 eigenbasis)")

    # Check if Sigma_T is diagonal in n-shell blocks
    n_shells = sorted(set(dl.n_fock for dl in dirac_labels))
    print(f"\n  Sigma_T block structure by n-shell:")
    for n1 in n_shells:
        for n2 in n_shells:
            if n2 < n1:
                continue
            idx1 = [i for i, dl in enumerate(dirac_labels) if dl.n_fock == n1]
            idx2 = [i for i, dl in enumerate(dirac_labels) if dl.n_fock == n2]
            block = Sigma_T[np.ix_(idx1, idx2)]
            block_norm = np.linalg.norm(block)
            if block_norm > 1e-10:
                print(f"    ({n1},{n2}): ||block|| = {block_norm:.6f}")

    return {
        "n_max": n_max,
        "n_modes": n_modes,
        "eigenvalues": nz_evals.tolist(),
        "q_assignments": q_assignments,
        "mode_profiles": mode_profiles,
    }


def main():
    t0 = time.time()

    print("="*70)
    print("TRANSVERSE PHOTON MODE ANALYSIS")
    print("Goal: identify photon quantum numbers in co-exact sector")
    print("="*70)

    results = {}
    results["n_max_3"] = analyze_coexact_modes(n_max=3)
    results["n_max_4"] = analyze_coexact_modes(n_max=4)

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total elapsed: {elapsed:.1f}s")

    # Save
    output_path = Path(__file__).parent / "data" / "transverse_mode_analysis.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def serializer(x):
        if isinstance(x, np.floating):
            return float(x)
        if isinstance(x, np.integer):
            return int(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, tuple):
            return list(x)
        return x

    # Convert tuple keys to strings for JSON
    json_results = {}
    for key, val in results.items():
        if isinstance(val, dict):
            clean = {}
            for k, v in val.items():
                if k == "mode_profiles":
                    # Convert tuple keys
                    clean[k] = [{str(kk): vv for kk, vv in prof.items()} for prof in v]
                else:
                    clean[k] = v
            json_results[key] = clean
        else:
            json_results[key] = val

    with open(output_path, "w") as f:
        json.dump(json_results, f, indent=2, default=serializer)

    print(f"Results saved: {output_path}")


if __name__ == "__main__":
    main()
