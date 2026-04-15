"""
V_ee spectral analysis in the graph eigenbasis for He.

Investigates the structure of V_ee when transformed to the basis where
h1 (the graph Laplacian one-body Hamiltonian) is diagonal. In this basis,
h1 is diagonal and V_ee carries ALL the correlation physics.

Key questions:
1. What is the effective rank of V_ee in the graph eigenbasis?
2. Where is the cusp singularity concentrated in matrix space?
3. How much correlation comes from off-diagonal V_ee elements?
4. What is the spectrum of V_ee (one-body and many-body)?

Track DI analysis script.
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.casimir_ci import (
    _build_graph_h1, two_electron_integral, build_graph_native_fci
)


def build_vee_one_body_matrix(Z: int, n_max: int, orbitals: list) -> np.ndarray:
    """Build the one-body V_ee matrix: V_ee[i,j] = sum_k <ij|1/r12|kj> averaged.

    Actually, V_ee is a TWO-body operator. For the one-body representation,
    we build the full 4-index tensor and contract/reshape as needed.

    Returns the full 4-index tensor g[i,j,k,l] = <ij|kl> at k_orb=Z.
    """
    n_spatial = len(orbitals)
    k_orb = float(Z)
    g = np.zeros((n_spatial, n_spatial, n_spatial, n_spatial))
    for i in range(n_spatial):
        ni, li, mi = orbitals[i]
        for j in range(n_spatial):
            nj, lj, mj = orbitals[j]
            for k in range(n_spatial):
                nk, lk, mk = orbitals[k]
                for l_idx in range(n_spatial):
                    nl, ll, ml = orbitals[l_idx]
                    g[i, j, k, l_idx] = two_electron_integral(
                        ni, li, mi, nj, lj, mj,
                        nk, lk, mk, nl, ll, ml, k_orb)
    return g


def analyze_vee_graph_eigenbasis(Z: int, n_max: int) -> dict:
    """Full spectral analysis of V_ee in the graph eigenbasis."""
    print(f"\n{'='*70}")
    print(f"  n_max = {n_max}, Z = {Z}")
    print(f"{'='*70}")

    # Step 1: Build h1 and get orbital list
    h1_mat, orbitals = _build_graph_h1(Z, n_max)
    n_spatial = len(orbitals)
    print(f"  n_spatial = {n_spatial} orbitals")
    print(f"  Orbitals: {orbitals}")

    # Step 2: Diagonalize h1 by (l,m) sector (same as build_graph_consistent_fci)
    from typing import Dict, Tuple, List
    sectors: Dict[Tuple[int, int], List[int]] = {}
    for idx, (n, l, m) in enumerate(orbitals):
        key = (l, m)
        if key not in sectors:
            sectors[key] = []
        sectors[key].append(idx)

    U = np.eye(n_spatial)
    evals = np.zeros(n_spatial)
    orbital_labels_eigen = [''] * n_spatial

    for (l, m), indices in sectors.items():
        if len(indices) == 1:
            idx = indices[0]
            evals[idx] = h1_mat[idx, idx]
            orbital_labels_eigen[idx] = f"({orbitals[idx][0]},{l},{m})"
        else:
            block = h1_mat[np.ix_(indices, indices)]
            block_evals, block_U = np.linalg.eigh(block)
            for i_local, i_global in enumerate(indices):
                evals[i_global] = block_evals[i_local]
                orbital_labels_eigen[i_global] = f"e{i_local}({l},{m})"
                for j_local, j_global in enumerate(indices):
                    U[j_global, i_global] = block_U[j_local, i_local]

    # Verify U is orthogonal
    orth_err = np.max(np.abs(U.T @ U - np.eye(n_spatial)))
    print(f"  U orthogonality error: {orth_err:.2e}")

    # h1 eigenvalues
    print(f"\n  h1 eigenvalues (graph eigenbasis):")
    sort_idx = np.argsort(evals)
    for i in sort_idx:
        print(f"    {orbital_labels_eigen[i]:>12s}: {evals[i]:.6f} Ha")

    h1_spread = evals.max() - evals.min()

    # Step 3: Build V_ee 4-index tensor in (n,l,m) basis
    print(f"\n  Building V_ee tensor ({n_spatial}^4 = {n_spatial**4} elements)...")
    g_nlm = build_vee_one_body_matrix(Z, n_max, orbitals)
    print(f"  Done. Nonzero elements: {np.count_nonzero(np.abs(g_nlm) > 1e-15)}/{n_spatial**4}")

    # Step 4: Transform to graph eigenbasis
    g_graph = np.einsum('ia,jb,ijkl,kc,ld->abcd', U, U, g_nlm, U, U, optimize=True)

    # ======================================================================
    # ANALYSIS A: One-body-like projections of V_ee
    # ======================================================================

    # Contract V_ee to a "Coulomb" one-body matrix: J[a,c] = sum_b g[a,b,c,b]
    # This is the mean-field part of V_ee
    J_nlm = np.einsum('aibj,bj->ai', g_nlm.reshape(n_spatial, n_spatial, n_spatial, n_spatial), np.eye(n_spatial))
    # Actually: J[a,c] = sum_b <ab|cb> (summing over occupied orbital b in mean-field)
    # For a proper analysis, let's look at the RAW tensor structure

    # Step 5: Reshape to 2-body matrix for spectral analysis
    # V_ee as a matrix on the 2-particle space: V[(i,j),(k,l)] = g[i,j,k,l]
    # This is the "supermatrix" form
    V_super_nlm = g_nlm.reshape(n_spatial**2, n_spatial**2)
    V_super_graph = g_graph.reshape(n_spatial**2, n_spatial**2)

    print(f"\n  --- V_ee SUPERMATRIX (one-particle-pair space, {n_spatial**2}x{n_spatial**2}) ---")

    # Eigenvalues of V_ee supermatrix in graph eigenbasis
    evals_V = np.linalg.eigvalsh(V_super_graph)
    evals_V_sorted = np.sort(evals_V)[::-1]  # descending
    print(f"  V_ee supermatrix eigenvalues (top 10 / bottom 5):")
    for i in range(min(10, len(evals_V_sorted))):
        print(f"    lam_{i+1} = {evals_V_sorted[i]:.6f}")
    if len(evals_V_sorted) > 10:
        print(f"    ...")
        for i in range(max(10, len(evals_V_sorted)-5), len(evals_V_sorted)):
            print(f"    lam_{i+1} = {evals_V_sorted[i]:.6f}")

    # Singular values
    sv = np.linalg.svd(V_super_graph, compute_uv=False)
    sv_sorted = np.sort(sv)[::-1]
    print(f"\n  Singular values of V_ee supermatrix (top 10):")
    for i in range(min(10, len(sv_sorted))):
        print(f"    sv_{i+1} = {sv_sorted[i]:.6f}")

    # Effective rank
    sv_max = sv_sorted[0]
    thresholds = [0.01, 0.001, 1e-4, 1e-6]
    print(f"\n  Effective rank (singular values > threshold * sv_max):")
    for t in thresholds:
        rank = np.sum(sv_sorted > t * sv_max)
        print(f"    threshold {t:.0e}: rank {rank}/{len(sv_sorted)}")

    # Frobenius norm
    frob = np.linalg.norm(V_super_graph, 'fro')
    print(f"\n  Frobenius norm of V_ee supermatrix: {frob:.6f}")
    print(f"  h1 eigenvalue spread: {h1_spread:.6f}")
    print(f"  Ratio Frob(V_ee) / |h1 spread|: {frob / abs(h1_spread):.4f}")

    # Sparsity
    total = V_super_graph.size
    nonzero_1e6 = np.count_nonzero(np.abs(V_super_graph) > 1e-6)
    nonzero_1e3 = np.count_nonzero(np.abs(V_super_graph) > 1e-3)
    nonzero_1e1 = np.count_nonzero(np.abs(V_super_graph) > 1e-1)
    print(f"\n  Sparsity of V_ee supermatrix:")
    print(f"    |V| > 1e-1: {nonzero_1e1}/{total} ({100*nonzero_1e1/total:.1f}%)")
    print(f"    |V| > 1e-3: {nonzero_1e3}/{total} ({100*nonzero_1e3/total:.1f}%)")
    print(f"    |V| > 1e-6: {nonzero_1e6}/{total} ({100*nonzero_1e6/total:.1f}%)")

    # ======================================================================
    # ANALYSIS B: Diagonal vs off-diagonal V_ee in FCI
    # ======================================================================
    print(f"\n  --- FCI DIAGONAL VS FULL V_ee ---")

    # Build full FCI
    H_fci = build_graph_native_fci(Z, n_max, m_total=0, spin='singlet')
    n_configs = H_fci.shape[0]
    fci_evals = np.linalg.eigvalsh(H_fci)
    E_fci = fci_evals[0]
    print(f"  FCI ground state: {E_fci:.6f} Ha ({n_configs} configs)")

    # Build "diagonal V_ee" FCI: keep only diagonal elements of V_ee in graph eigenbasis
    # We need to build the FCI matrix using graph-eigenbasis h1 (diagonal) and
    # only diagonal V_ee
    # Reconstruct configs in graph eigenbasis
    eigen_m = [0] * n_spatial
    for (l, m), indices in sectors.items():
        for idx in indices:
            eigen_m[idx] = m

    configs = []
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if eigen_m[i] + eigen_m[j] == 0:
                configs.append((i, j))

    H_diag_vee = np.zeros((n_configs, n_configs))
    for I in range(n_configs):
        i, j = configs[I]
        for J in range(I, n_configs):
            p, q = configs[J]

            bra_perms = [(i, j)]
            if i != j:
                bra_perms.append((j, i))
            ket_perms = [(p, q)]
            if p != q:
                ket_perms.append((q, p))

            N_IJ = np.sqrt(float(len(bra_perms)))
            N_PQ = np.sqrt(float(len(ket_perms)))

            me = 0.0
            for a, b in bra_perms:
                for c, d in ket_perms:
                    # One-body: diagonal in eigenbasis
                    if a == c and b == d:
                        me += evals[a] + evals[b]
                    # Two-body: ONLY diagonal g[a,b,a,b] terms
                    # Diagonal means the bra orbital pair equals the ket orbital pair
                    if a == c and b == d:
                        me += g_graph[a, b, a, b]

            me /= (N_IJ * N_PQ)
            H_diag_vee[I, J] = me
            H_diag_vee[J, I] = me

    diag_evals = np.linalg.eigvalsh(H_diag_vee)
    E_diag = diag_evals[0]

    E_exact = -2.903724377034
    err_fci = (E_fci - E_exact) / abs(E_exact) * 100
    err_diag = (E_diag - E_exact) / abs(E_exact) * 100

    print(f"  Diagonal-V_ee ground state: {E_diag:.6f} Ha")
    print(f"  Full FCI error: {err_fci:.3f}%")
    print(f"  Diagonal-only error: {err_diag:.3f}%")
    print(f"  Off-diagonal V_ee contribution: {E_fci - E_diag:.6f} Ha")
    print(f"  Off-diagonal as % of correlation: {100*(E_fci - E_diag)/(E_fci - E_diag + (E_diag - E_exact)):.1f}%"
          if abs(E_diag - E_exact) > 1e-10 else "")

    # ======================================================================
    # ANALYSIS C: Largest V_ee matrix elements in graph eigenbasis
    # ======================================================================
    print(f"\n  --- LARGEST V_ee ELEMENTS IN GRAPH EIGENBASIS ---")

    # Collect all (i,j,k,l) with large |g_graph[i,j,k,l]|
    flat_indices = np.argsort(np.abs(g_graph).ravel())[::-1]
    print(f"  Top 20 |V_ee[a,b,c,d]| in graph eigenbasis:")
    count = 0
    seen = set()
    for flat_idx in flat_indices:
        if count >= 20:
            break
        idx_tuple = np.unravel_index(flat_idx, g_graph.shape)
        a, b, c, d = idx_tuple
        # Skip permutation-related duplicates
        canon = tuple(sorted([(a, c), (b, d)]))
        if canon in seen:
            continue
        seen.add(canon)
        val = g_graph[a, b, c, d]
        la = orbital_labels_eigen[a]
        lb = orbital_labels_eigen[b]
        lc = orbital_labels_eigen[c]
        ld = orbital_labels_eigen[d]
        # Classify: diagonal (a=c,b=d), exchange (a=d,b=c), or off-diag
        if a == c and b == d:
            typ = "DIAG"
        elif a == d and b == c:
            typ = "EXCH"
        else:
            typ = "OFF"
        print(f"    V[{la},{lb},{lc},{ld}] = {val:+.6f}  [{typ}]")
        count += 1

    # ======================================================================
    # ANALYSIS D: Angular momentum channel analysis
    # ======================================================================
    print(f"\n  --- ANGULAR MOMENTUM CHANNEL COUPLING ---")

    # For each pair of (l_a, l_b) sectors, compute the total V_ee coupling strength
    l_values = sorted(set(orbitals[i][1] for i in range(n_spatial)))
    print(f"  l values present: {l_values}")

    # In the (n,l,m) basis, which l-l' channels have strongest V_ee?
    for la in l_values:
        for lb in l_values:
            if lb < la:
                continue
            total_coupling = 0.0
            count_elem = 0
            for i in range(n_spatial):
                if orbitals[i][1] != la:
                    continue
                for j in range(n_spatial):
                    if orbitals[j][1] != lb:
                        continue
                    for k in range(n_spatial):
                        for l_idx in range(n_spatial):
                            val = g_nlm[i, j, k, l_idx]
                            if abs(val) > 1e-15:
                                total_coupling += val**2
                                count_elem += 1
            if count_elem > 0:
                print(f"  Channel l=({la},{lb}): Frobenius contribution = {np.sqrt(total_coupling):.6f}, "
                      f"nonzero elements = {count_elem}")

    # ======================================================================
    # ANALYSIS E: FCI V_ee spectrum (many-body matrix)
    # ======================================================================
    print(f"\n  --- FCI MANY-BODY V_ee MATRIX ---")

    # Build the FCI matrix with ONLY V_ee (no h1)
    H_vee_only = np.zeros((n_configs, n_configs))
    for I in range(n_configs):
        i_conf, j_conf = configs[I]
        for J in range(I, n_configs):
            p_conf, q_conf = configs[J]

            bra_perms = [(i_conf, j_conf)]
            if i_conf != j_conf:
                bra_perms.append((j_conf, i_conf))
            ket_perms = [(p_conf, q_conf)]
            if p_conf != q_conf:
                ket_perms.append((q_conf, p_conf))

            N_IJ = np.sqrt(float(len(bra_perms)))
            N_PQ = np.sqrt(float(len(ket_perms)))

            me = 0.0
            for a, b in bra_perms:
                for c, d in ket_perms:
                    me += g_graph[a, b, c, d]

            me /= (N_IJ * N_PQ)
            H_vee_only[I, J] = me
            H_vee_only[J, I] = me

    vee_evals = np.linalg.eigvalsh(H_vee_only)
    vee_evals_sorted = np.sort(vee_evals)[::-1]
    print(f"  V_ee FCI matrix eigenvalues ({n_configs}x{n_configs}):")
    for i in range(min(15, n_configs)):
        print(f"    lam_{i+1} = {vee_evals_sorted[i]:.6f}")
    if n_configs > 15:
        print(f"    ... ({n_configs - 15} more)")

    vee_sv = np.linalg.svd(H_vee_only, compute_uv=False)
    vee_sv_sorted = np.sort(vee_sv)[::-1]
    vee_sv_max = vee_sv_sorted[0]
    print(f"\n  V_ee FCI matrix singular values (top 10):")
    for i in range(min(10, len(vee_sv_sorted))):
        print(f"    sv_{i+1} = {vee_sv_sorted[i]:.6f}")

    print(f"\n  V_ee FCI effective rank:")
    for t in [0.01, 0.001, 1e-4, 1e-6]:
        rank = np.sum(vee_sv_sorted > t * vee_sv_max)
        print(f"    threshold {t:.0e}: rank {rank}/{n_configs}")

    # Also build h1-only FCI matrix for comparison
    H_h1_only = np.zeros((n_configs, n_configs))
    for I in range(n_configs):
        i_conf, j_conf = configs[I]
        for J in range(I, n_configs):
            p_conf, q_conf = configs[J]

            bra_perms = [(i_conf, j_conf)]
            if i_conf != j_conf:
                bra_perms.append((j_conf, i_conf))
            ket_perms = [(p_conf, q_conf)]
            if p_conf != q_conf:
                ket_perms.append((q_conf, p_conf))

            N_IJ = np.sqrt(float(len(bra_perms)))
            N_PQ = np.sqrt(float(len(ket_perms)))

            me = 0.0
            for a, b in bra_perms:
                for c, d in ket_perms:
                    if a == c and b == d:
                        me += evals[a] + evals[b]

            me /= (N_IJ * N_PQ)
            H_h1_only[I, J] = me
            H_h1_only[J, I] = me

    # Verify: H_full = H_h1 + H_vee
    H_full = H_h1_only + H_vee_only
    full_evals = np.linalg.eigvalsh(H_full)
    print(f"\n  Verification: H_h1 + H_vee ground state = {full_evals[0]:.6f} Ha")
    print(f"  vs FCI ground state = {E_fci:.6f} Ha")
    print(f"  Difference: {abs(full_evals[0] - E_fci):.2e} Ha")

    return {
        'n_max': n_max,
        'n_spatial': n_spatial,
        'n_configs': n_configs,
        'h1_evals': evals,
        'h1_spread': h1_spread,
        'g_graph': g_graph,
        'g_nlm': g_nlm,
        'U': U,
        'orbitals': orbitals,
        'orbital_labels_eigen': orbital_labels_eigen,
        'E_fci': E_fci,
        'E_diag': E_diag,
        'vee_fci_evals': vee_evals_sorted,
        'vee_super_sv': sv_sorted,
        'vee_fci_sv': vee_sv_sorted,
    }


def make_heatmap(results: dict, save_path: str):
    """Create heatmap of |V_ee| in graph eigenbasis."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    g_graph = results['g_graph']
    n_spatial = results['n_spatial']
    labels = results['orbital_labels_eigen']

    # Reshape to supermatrix for visualization
    V_super = np.abs(g_graph.reshape(n_spatial**2, n_spatial**2))
    # Replace zeros with small value for log scale
    V_super = np.where(V_super > 1e-16, V_super, 1e-16)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Left: supermatrix heatmap
    ax = axes[0]
    vmin = 1e-6
    vmax = V_super.max()
    im = ax.imshow(V_super, norm=LogNorm(vmin=vmin, vmax=vmax),
                   cmap='viridis', aspect='auto')
    ax.set_title(f'|V_ee| supermatrix (graph eigenbasis, n_max={results["n_max"]})')
    ax.set_xlabel('Ket pair index (c,d)')
    ax.set_ylabel('Bra pair index (a,b)')
    plt.colorbar(im, ax=ax, label='|V_ee| (Ha)')

    # Right: diagonal block structure
    # Extract the "direct" part: V_direct[a,b] = g_graph[a,b,a,b]
    V_direct = np.zeros((n_spatial, n_spatial))
    for a in range(n_spatial):
        for b in range(n_spatial):
            V_direct[a, b] = g_graph[a, b, a, b]

    ax = axes[1]
    vabs = np.abs(V_direct)
    vabs = np.where(vabs > 1e-16, vabs, 1e-16)
    im = ax.imshow(vabs, norm=LogNorm(vmin=1e-6, vmax=vabs.max()),
                   cmap='inferno', aspect='auto')
    ax.set_title(f'|V_direct[a,b]| = |<ab|ab>| (graph eigenbasis)')
    ax.set_xlabel('Orbital b')
    ax.set_ylabel('Orbital a')
    if n_spatial <= 14:
        ax.set_xticks(range(n_spatial))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=7)
        ax.set_yticks(range(n_spatial))
        ax.set_yticklabels(labels, fontsize=7)
    plt.colorbar(im, ax=ax, label='|V_direct| (Ha)')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"\n  Heatmap saved to {save_path}")
    plt.close()


def make_singular_value_plot(all_results: list, save_path: str):
    """Plot singular value decay across n_max values."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: V_ee supermatrix singular values
    ax = axes[0]
    for r in all_results:
        sv = r['vee_super_sv']
        sv_norm = sv / sv[0]
        ax.semilogy(range(1, len(sv_norm)+1), sv_norm,
                    'o-', markersize=3, label=f'n_max={r["n_max"]}')
    ax.set_xlabel('Singular value index')
    ax.set_ylabel('σ / σ_max')
    ax.set_title('V_ee supermatrix singular value decay')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Right: V_ee FCI matrix singular values
    ax = axes[1]
    for r in all_results:
        sv = r['vee_fci_sv']
        sv_norm = sv / sv[0]
        ax.semilogy(range(1, len(sv_norm)+1), sv_norm,
                    'o-', markersize=3, label=f'n_max={r["n_max"]}')
    ax.set_xlabel('Singular value index')
    ax.set_ylabel('σ / σ_max')
    ax.set_title('V_ee FCI matrix singular value decay')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"\n  Singular value plot saved to {save_path}")
    plt.close()


def main():
    print("V_ee Spectral Analysis in Graph Eigenbasis for He")
    print("=" * 70)

    all_results = []
    for n_max in [2, 3, 4, 5]:
        result = analyze_vee_graph_eigenbasis(Z=2, n_max=n_max)
        all_results.append(result)

    # Summary table
    print("\n\n" + "=" * 70)
    print("  SUMMARY TABLE")
    print("=" * 70)
    E_exact = -2.903724377034
    print(f"  {'n_max':>5} {'n_sp':>5} {'n_cfg':>6} {'E_FCI':>12} {'E_diag':>12} "
          f"{'err_FCI%':>9} {'err_diag%':>10} {'dE_offdiag':>11} {'rank_1%':>8} {'rank_0.1%':>10}")
    for r in all_results:
        err_fci = (r['E_fci'] - E_exact) / abs(E_exact) * 100
        err_diag = (r['E_diag'] - E_exact) / abs(E_exact) * 100
        delta = r['E_fci'] - r['E_diag']
        sv = r['vee_super_sv']
        rank_1 = np.sum(sv > 0.01 * sv[0])
        rank_01 = np.sum(sv > 0.001 * sv[0])
        print(f"  {r['n_max']:>5} {r['n_spatial']:>5} {r['n_configs']:>6} "
              f"{r['E_fci']:>12.6f} {r['E_diag']:>12.6f} "
              f"{err_fci:>9.3f} {err_diag:>10.3f} {delta:>11.6f} "
              f"{rank_1:>8} {rank_01:>10}")

    # Make plots for n_max=4
    plot_dir = os.path.join(os.path.dirname(__file__), 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    # Use n_max=4 for the detailed heatmap
    r4 = [r for r in all_results if r['n_max'] == 4][0]
    make_heatmap(r4, os.path.join(plot_dir, 'vee_graph_eigenbasis.png'))

    # Singular value decay across all n_max
    make_singular_value_plot(all_results, os.path.join(plot_dir, 'vee_sv_decay.png'))

    print("\nDone.")


if __name__ == '__main__':
    main()
