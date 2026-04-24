"""
Follow-up analysis of Dirac native graph spectrum.
====================================================

Key observations from the first run that need deeper investigation:

1. Rule A per-kappa sectors: each kappa sector at fixed n_max has eigenvalues
   that are PATH GRAPH eigenvalues. The m_j direction is a path P_{2j+1}
   and the n direction is a path P_{n_count}. The Laplacian of the
   Cartesian product P_a x P_b has eigenvalues mu_i + nu_j where mu, nu
   are path Laplacian eigenvalues. Verify this.

2. Weighted Rule B at n_max=2 produced eigenvalue 1.5 with degeneracy 4,
   matching CH Dirac. Check if this persists at larger n_max.

3. The scalar graph eigenvalues at n_max are NOT just n^2-1 either -- the
   finite graph Laplacian eigenvalues only become n^2-1 in a specific scaling
   limit. Need to understand the relationship more carefully.

Author: GeoVac exploration, April 2026.
"""

from __future__ import annotations
import json
import numpy as np
from geovac.ihara_zeta_dirac import build_dirac_s3_graph
from geovac.lattice import GeometricLattice
from geovac.dirac_matrix_elements import DiracLabel, kappa_to_l, kappa_to_j

np.set_printoptions(precision=12, linewidth=120, suppress=True)


def eigenvalues_of_laplacian(A):
    D = np.diag(A.sum(axis=1).astype(float))
    L = D - A.astype(float)
    return np.sort(np.linalg.eigvalsh(L))


def path_graph_laplacian_eigenvalues(n: int):
    """Eigenvalues of the Laplacian of path graph P_n: 2 - 2*cos(pi*k/n) for k=0,...,n-1."""
    return np.sort([2 - 2*np.cos(np.pi*k/n) for k in range(n)])


def group_eigenvalues(evals, tol=1e-8):
    if len(evals) == 0:
        return []
    groups = []
    current_val = evals[0]
    current_count = 1
    for i in range(1, len(evals)):
        if abs(evals[i] - current_val) < tol:
            current_count += 1
        else:
            groups.append((float(np.mean(evals[max(0,i-current_count):i])), current_count))
            current_val = evals[i]
            current_count = 1
    groups.append((float(np.mean(evals[len(evals)-current_count:])), current_count))
    return groups


def main():
    print("="*80)
    print("  FOLLOW-UP ANALYSIS 1: Rule A per-kappa sectors as path products")
    print("="*80)

    for n_max in [2, 3, 4]:
        print(f"\n  n_max = {n_max}:")
        A_A, labels_A, deg_A, desc_A = build_dirac_s3_graph(n_max, "A")

        kappa_values = sorted(set(lab.kappa for lab in labels_A))
        for kappa in kappa_values:
            indices = [i for i, lab in enumerate(labels_A) if lab.kappa == kappa]
            if len(indices) <= 1:
                continue

            l = kappa_to_l(kappa)
            j = kappa_to_j(kappa)
            two_j = int(2 * float(j))

            # The kappa sector has:
            #   - n_fock values from max(1, l+1) to n_max
            #   - m_j values from -j to j in steps of 1 (count = 2j+1)
            n_vals = sorted(set(labels_A[i].n_fock for i in indices))
            mj_vals = sorted(set(labels_A[i].two_m_j for i in indices))

            n_count = len(n_vals)
            mj_count = len(mj_vals)

            # Rule A adjacency: dn=+/-1 at fixed mj, OR dmj=+/-1 at fixed n
            # This is EXACTLY the Cartesian product graph P_{n_count} x P_{mj_count}

            # Laplacian eigenvalues of P_a x P_b:
            # mu_i + nu_j where mu = path eigenvalues of P_a, nu = path eigenvalues of P_b
            path_n_evals = path_graph_laplacian_eigenvalues(n_count)
            path_mj_evals = path_graph_laplacian_eigenvalues(mj_count)

            # Tensor sum
            predicted = np.sort([mu + nu for mu in path_n_evals for nu in path_mj_evals])

            # Actual sub-graph eigenvalues
            sub_A = A_A[np.ix_(indices, indices)]
            actual = eigenvalues_of_laplacian(sub_A)

            # Compare
            max_diff = np.max(np.abs(predicted - actual))
            print(f"\n    kappa={kappa:+d} (l={l}, j={j}): P_{n_count} x P_{mj_count}")
            print(f"      Path n-direction eigenvalues: {path_n_evals}")
            print(f"      Path mj-direction eigenvalues: {path_mj_evals}")
            print(f"      Predicted (tensor sum): {np.sort(predicted)}")
            print(f"      Actual sub-graph:       {actual}")
            print(f"      Max difference: {max_diff:.2e}")
            if max_diff < 1e-10:
                print(f"      >>> EXACT MATCH: kappa sector = P_{n_count} x P_{mj_count}")

    # =====================================================================
    print("\n\n" + "="*80)
    print("  FOLLOW-UP ANALYSIS 2: Scalar graph per-l sectors as path products")
    print("="*80)

    for n_max in [2, 3, 4]:
        print(f"\n  n_max = {n_max}:")
        lattice = GeometricLattice(max_n=n_max, nuclear_charge=1, topological_weights=False)
        A_s = lattice.adjacency.toarray()
        states = lattice.states

        l_values = sorted(set(s[1] for s in states))
        for l in l_values:
            indices = [i for i, s in enumerate(states) if s[1] == l]
            if len(indices) <= 1:
                continue

            n_vals = sorted(set(states[i][0] for i in indices))
            m_vals = sorted(set(states[i][2] for i in indices))

            n_count = len(n_vals)
            m_count = len(m_vals)

            path_n_evals = path_graph_laplacian_eigenvalues(n_count)
            path_m_evals = path_graph_laplacian_eigenvalues(m_count)

            predicted = np.sort([mu + nu for mu in path_n_evals for nu in path_m_evals])

            sub_A = A_s[np.ix_(indices, indices)]
            actual = eigenvalues_of_laplacian(sub_A)

            max_diff = np.max(np.abs(predicted - actual))
            print(f"\n    l={l}: P_{n_count} x P_{m_count}")
            print(f"      Path n-direction eigenvalues: {path_n_evals}")
            print(f"      Path m-direction eigenvalues: {path_m_evals}")
            print(f"      Predicted (tensor sum): {np.sort(predicted)}")
            print(f"      Actual sub-graph:       {actual}")
            print(f"      Max difference: {max_diff:.2e}")
            if max_diff < 1e-10:
                print(f"      >>> EXACT MATCH: l-sector = P_{n_count} x P_{m_count}")

    # =====================================================================
    print("\n\n" + "="*80)
    print("  FOLLOW-UP ANALYSIS 3: What the scalar graph Laplacian spectrum REALLY is")
    print("="*80)
    print("  The scalar graph eigenvalues n^2-1 are from the TOPOLOGICALLY-WEIGHTED")
    print("  graph (w = 1/(n1*n2)), NOT the unweighted graph (w = 1).")

    for n_max in [2, 3, 4, 5]:
        print(f"\n  n_max = {n_max}:")
        # Unweighted
        lattice_uw = GeometricLattice(max_n=n_max, nuclear_charge=1, topological_weights=False)
        A_uw = lattice_uw.adjacency.toarray()
        evals_uw = eigenvalues_of_laplacian(A_uw)
        groups_uw = group_eigenvalues(evals_uw)
        print(f"    Unweighted graph L eigenvalues: {groups_uw}")

        # Weighted (topological)
        lattice_tw = GeometricLattice(max_n=n_max, nuclear_charge=1, topological_weights=True)
        A_tw = lattice_tw.adjacency.toarray()
        evals_tw = eigenvalues_of_laplacian(A_tw)
        groups_tw = group_eigenvalues(evals_tw)
        print(f"    Weighted graph L eigenvalues:   {groups_tw}")

        # Compare with n^2-1 scaled by kappa = -1/16
        scaled_target = [(n**2 - 1) * (-1/16) for n in range(1, n_max+1)]
        print(f"    Target kappa*(n^2-1): {scaled_target}")

        # The actual Hamiltonian in geovac is H = kappa*(D - A + W) where W = -Z/n^2
        # The graph Laplacian D-A has path graph eigenvalues in each l-sector.
        # The FULL Hamiltonian = kappa*(D - A) + kappa*diag(-Z/n^2)
        # Eigenvalues are NOT simply n^2-1.

        # Build single-particle H = kappa*(D - A) + kappa*diag(-Z/n^2) directly
        kappa_val = -1.0/16.0
        Z_val = 1
        states_uw = lattice_uw.states
        D_uw = np.diag(A_uw.sum(axis=1))
        L_uw = D_uw - A_uw
        W_uw = np.diag([float(-Z_val / s[0]**2) for s in states_uw])
        H_uw = kappa_val * (L_uw + W_uw)
        h_evals = np.sort(np.linalg.eigvalsh(H_uw))
        h_groups = group_eigenvalues(h_evals)
        print(f"    Unweighted H = kappa*(L + W) eigenvalues: {h_groups}")

        # Expected: -Z^2/(2n^2) Rydberg energies
        rydberg = [-1.0 / (2.0 * n**2) for n in range(1, n_max+1)]
        print(f"    Expected Rydberg -1/(2n^2): {rydberg}")

    # =====================================================================
    print("\n\n" + "="*80)
    print("  FOLLOW-UP ANALYSIS 4: Does D - A of the topologically-weighted scalar graph")
    print("  give n^2-1 eigenvalues?")
    print("="*80)

    for n_max in [2, 3, 4, 5]:
        print(f"\n  n_max = {n_max}:")
        lattice_tw = GeometricLattice(max_n=n_max, nuclear_charge=1, topological_weights=True)
        A_tw = lattice_tw.adjacency.toarray()
        evals_tw = eigenvalues_of_laplacian(A_tw)
        groups_tw = group_eigenvalues(evals_tw)
        print(f"    Weighted (D-A) groups: {groups_tw}")

        # Rescale by -16 (multiply by -kappa^{-1})
        rescaled = [(-16 * v, d) for v, d in groups_tw]
        print(f"    Rescaled by -16:       {rescaled}")

        # Compare to -(n^2-1) pattern
        target = {float(n**2 - 1): n**2 for n in range(1, n_max+1)}
        print(f"    Target n^2-1 (value: deg): {target}")

    # =====================================================================
    print("\n\n" + "="*80)
    print("  FOLLOW-UP ANALYSIS 5: The scalar adjacency A with w=1/(n1*n2)")
    print("  Check: does scale*(D-A) + diagonal(-Z/n^2) give -Z^2/(2n^2)?")
    print("="*80)

    for n_max in [3, 4]:
        print(f"\n  n_max = {n_max}:")
        lattice_tw = GeometricLattice(max_n=n_max, nuclear_charge=1, topological_weights=True)
        A_tw = lattice_tw.adjacency.toarray()
        states = lattice_tw.states

        # Build H = kappa * (D - A) + kappa * diag(-Z/n^2)
        kappa = -1.0/16.0
        Z = 1
        D_tw = np.diag(A_tw.sum(axis=1))
        L_tw = D_tw - A_tw
        W = np.diag([float(-Z / s[0]**2) for s in states])
        H = kappa * (L_tw + W)
        h_evals = np.sort(np.linalg.eigvalsh(H))
        h_groups = group_eigenvalues(h_evals)
        print(f"    H = kappa*(L + W) eigenvalues: {h_groups}")

        rydberg = [(-1.0 / (2.0 * n**2), n**2) for n in range(1, n_max+1)]
        print(f"    Expected Rydberg: {rydberg}")

    # =====================================================================
    print("\n\n" + "="*80)
    print("  FOLLOW-UP ANALYSIS 6: Weighted Rule B hit at n_max=2")
    print("  Does the 1.5 (deg 4) persist in larger graphs?")
    print("="*80)

    for n_max in [2, 3, 4]:
        A_B, labels_B, _, desc = build_dirac_s3_graph(n_max, "B")
        A_wB = np.zeros_like(A_B, dtype=float)
        V = len(labels_B)
        for i in range(V):
            for j in range(i + 1, V):
                if A_B[i, j] == 1:
                    n1 = labels_B[i].n_fock
                    n2 = labels_B[j].n_fock
                    w = 1.0 / (n1 * n2)
                    A_wB[i, j] = w
                    A_wB[j, i] = w

        evals = eigenvalues_of_laplacian(A_wB)
        groups = group_eigenvalues(evals)
        print(f"\n  n_max={n_max}: Weighted Rule B eigenvalue groups:")
        for v, d in groups:
            marker = " <-- CH match!" if abs(v - 1.5) < 0.01 else ""
            print(f"    {v:12.8f}  deg = {d}{marker}")

    # =====================================================================
    print("\n\n" + "="*80)
    print("  FOLLOW-UP ANALYSIS 7: Total eigenvalue sums and traces")
    print("="*80)

    for n_max in [2, 3, 4]:
        for rule in ["A", "B"]:
            A, labels, deg, desc = build_dirac_s3_graph(n_max, rule)
            L = np.diag(A.sum(axis=1).astype(float)) - A.astype(float)
            tr = np.trace(L)
            total_deg = int(A.sum()) // 2
            evals = eigenvalues_of_laplacian(A)
            eval_sum = np.sum(evals)
            print(f"  Rule {rule} n_max={n_max}: tr(L) = {tr:.4f}, "
                  f"sum(evals) = {eval_sum:.4f}, "
                  f"2*|E| = {int(A.sum())}, V = {len(labels)}")

    # Compare with scalar
    for n_max in [2, 3, 4]:
        lattice = GeometricLattice(max_n=n_max, nuclear_charge=1, topological_weights=False)
        A_s = lattice.adjacency.toarray()
        L_s = np.diag(A_s.sum(axis=1)) - A_s
        tr = np.trace(L_s)
        evals = eigenvalues_of_laplacian(A_s)
        print(f"  Scalar n_max={n_max}: tr(L) = {tr:.4f}, "
              f"sum(evals) = {np.sum(evals):.4f}, "
              f"2*|E| = {int(A_s.sum())}, V = {len(lattice.states)}")

    print("\n\nDone.")


if __name__ == "__main__":
    main()
