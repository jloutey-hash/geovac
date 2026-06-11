"""Test whether the full Hamiltonian on the (n,κ,m_j) Dirac graph
reproduces Dirac-Coulomb fine-structure energies.

Construction mirrors the scalar graph-native h1:
  scalar:  h1[i,i] = -Z²/(2n²),  h1[i,j] = κ·(-A[i,j])
  Dirac:   h1[i,i] = -Z²/(2n²) + ΔE_FS(n,j),  h1[i,j] = κ·(-A_dirac[i,j])

Three diagonal variants tested:
  (A) Rydberg only:  -Z²/(2n²)                         [α=0 limit]
  (B) SO only:       -Z²/(2n²) + H_SO(n,κ,α)           [spin-orbit, no Darwin/MV]
  (C) Full FS:       -Z²/(2n²) + ΔE_FS(n,j,α)          [all α² terms]

κ_scale = -1/16 for off-diagonal, same as scalar graph.
"""

import numpy as np
import json
from pathlib import Path

# Import Dirac graph builder
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from geovac.ihara_zeta_dirac import build_dirac_s3_graph, DiracLabel
from geovac.dirac_matrix_elements import kappa_to_l
from geovac.lattice import GeometricLattice

ALPHA = 7.2973525693e-3  # fine-structure constant
KAPPA = -1.0 / 16.0      # universal topological constant


def dirac_coulomb_energy(n: int, kappa: int, Z: int = 1) -> float:
    """Exact Dirac-Coulomb energy in Hartree (rest mass subtracted)."""
    j = abs(kappa) - 0.5
    n_r = n - abs(kappa)
    if n_r < 0:
        return None
    gamma = np.sqrt(kappa**2 - (Z * ALPHA)**2)
    denom = np.sqrt(1 + (Z * ALPHA)**2 / (n_r + gamma)**2)
    c_au = 1.0 / ALPHA
    return c_au**2 * (1.0 / denom - 1.0)


def fine_structure_shift(n: int, kappa: int, Z: int = 1) -> float:
    """Full α² fine-structure shift: SO + Darwin + MV."""
    j = abs(kappa) - 0.5
    return -(Z**4 * ALPHA**2) / (2.0 * n**4) * (n / (j + 0.5) - 0.75)


def spin_orbit_shift(n: int, kappa: int, Z: int = 1) -> float:
    """Breit-Pauli H_SO diagonal. Zero at l=0 (Kramers)."""
    l = kappa_to_l(kappa)
    if l == 0:
        return 0.0
    return -(Z**4 * ALPHA**2 * (kappa + 1)) / (4.0 * n**3 * l * (l + 0.5) * (l + 1))


def build_hamiltonian(adj_matrix, labels, diagonal_fn, kappa_scale=KAPPA):
    """Build H = kappa_scale*(D-A) + diag(diagonal_fn(label))."""
    N = len(labels)
    A = adj_matrix.copy().astype(float)
    D = np.diag(A.sum(axis=1))
    L = D - A  # graph Laplacian

    H = kappa_scale * L
    for i, lab in enumerate(labels):
        H[i, i] += diagonal_fn(lab)
    return H


def scalar_reference(n_max: int, Z: int = 1):
    """Build scalar graph Hamiltonian and compare to Rydberg."""
    lattice = GeometricLattice(max_n=n_max, nuclear_charge=Z)
    orbitals = list(lattice.states)
    n_states = lattice.num_states

    A = lattice.adjacency
    if hasattr(A, 'toarray'):
        A = A.toarray()
    D = np.diag(A.sum(axis=1))
    L = D - A

    H = KAPPA * L
    for i, (n, l, m) in enumerate(orbitals):
        H[i, i] += -Z**2 / (2.0 * n**2)

    evals = np.sort(np.linalg.eigvalsh(H))
    rydberg = sorted([-Z**2 / (2.0 * n**2) for n in range(1, n_max + 1)
                       for l in range(n) for m in range(-l, l + 1)])
    rydberg = np.array(rydberg)

    return evals, rydberg, orbitals


def dirac_test(n_max: int, Z: int = 1):
    """Build Dirac graph Hamiltonian variants and compare to Dirac-Coulomb."""
    adj, labels, _deg, _desc = build_dirac_s3_graph(n_max, adjacency_rule='A')

    # Exact Dirac-Coulomb energies for each label
    exact_dirac = []
    for lab in labels:
        e = dirac_coulomb_energy(lab.n_fock, lab.kappa, Z)
        exact_dirac.append(e if e is not None else 0.0)
    exact_dirac = np.array(exact_dirac)

    # Variant A: Rydberg diagonal only (α=0 limit)
    def diag_rydberg(lab):
        return -Z**2 / (2.0 * lab.n_fock**2)

    H_A = build_hamiltonian(adj, labels, diag_rydberg)
    evals_A = np.sort(np.linalg.eigvalsh(H_A))

    # Variant B: Rydberg + spin-orbit only
    def diag_so(lab):
        return -Z**2 / (2.0 * lab.n_fock**2) + spin_orbit_shift(lab.n_fock, lab.kappa, Z)

    H_B = build_hamiltonian(adj, labels, diag_so)
    evals_B = np.sort(np.linalg.eigvalsh(H_B))

    # Variant C: Rydberg + full fine structure
    def diag_fs(lab):
        return -Z**2 / (2.0 * lab.n_fock**2) + fine_structure_shift(lab.n_fock, lab.kappa, Z)

    H_C = build_hamiltonian(adj, labels, diag_fs)
    evals_C = np.sort(np.linalg.eigvalsh(H_C))

    # Sorted exact energies for comparison
    exact_sorted = np.sort(exact_dirac)

    # Also compute: what are the DIAGONAL-ONLY eigenvalues (no graph coupling)?
    diag_only_ryd = np.sort(np.array([diag_rydberg(lab) for lab in labels]))
    diag_only_so = np.sort(np.array([diag_so(lab) for lab in labels]))
    diag_only_fs = np.sort(np.array([diag_fs(lab) for lab in labels]))

    return {
        'labels': labels,
        'evals_A': evals_A,
        'evals_B': evals_B,
        'evals_C': evals_C,
        'exact_dirac': exact_sorted,
        'diag_only_ryd': diag_only_ryd,
        'diag_only_so': diag_only_so,
        'diag_only_fs': diag_only_fs,
    }


def analyze_degeneracies(evals, tol=1e-10):
    """Group eigenvalues by degeneracy."""
    groups = []
    i = 0
    while i < len(evals):
        val = evals[i]
        count = 1
        while i + count < len(evals) and abs(evals[i + count] - val) < tol:
            count += 1
        groups.append((val, count))
        i += count
    return groups


def main():
    results = {}
    Z = 1  # hydrogen

    print("=" * 80)
    print("NATIVE DIRAC GRAPH: FULL HAMILTONIAN TEST")
    print("H = kappa*(D-A) + diag(V),  kappa = -1/16")
    print("=" * 80)

    # ---- Scalar reference ----
    print("\n--- SCALAR REFERENCE ---")
    for n_max in [2, 3, 4]:
        evals_s, rydberg_s, orbs = scalar_reference(n_max, Z)
        max_err = np.max(np.abs(evals_s - rydberg_s))
        rms_err = np.sqrt(np.mean((evals_s - rydberg_s)**2))
        print(f"  n_max={n_max}: {len(orbs)} states, max|ΔE|={max_err:.6f} Ha, "
              f"rms={rms_err:.6f} Ha")

        # Degeneracy check
        groups = analyze_degeneracies(evals_s, tol=1e-6)
        deg_str = ", ".join(f"{v:.6f}(×{c})" for v, c in groups[:6])
        print(f"    degens: {deg_str}")

        results[f'scalar_n{n_max}'] = {
            'n_states': len(orbs),
            'evals': evals_s.tolist(),
            'rydberg': rydberg_s.tolist(),
            'max_err': float(max_err),
            'rms_err': float(rms_err),
        }

    # ---- Dirac graph ----
    print("\n--- DIRAC GRAPH (Rule A, κ-preserving) ---")
    for n_max in [2, 3, 4]:
        d = dirac_test(n_max, Z)
        n_states = len(d['labels'])
        print(f"\n  n_max={n_max}: {n_states} Dirac nodes")

        # For each variant, compare to exact Dirac-Coulomb
        for tag, evals_key in [('A_rydberg', 'evals_A'),
                                ('B_so', 'evals_B'),
                                ('C_fullFS', 'evals_C')]:
            evals = d[evals_key]
            exact = d['exact_dirac']

            # Since eigenvalue ordering may differ, compare sorted
            max_err = np.max(np.abs(evals - exact))
            rms_err = np.sqrt(np.mean((evals - exact)**2))

            # Also compare to diagonal-only (no graph coupling)
            diag_key = {
                'A_rydberg': 'diag_only_ryd',
                'B_so': 'diag_only_so',
                'C_fullFS': 'diag_only_fs',
            }[tag]
            diag_evals = d[diag_key]
            diag_err = np.max(np.abs(diag_evals - exact))

            print(f"    {tag}: max|ΔE|={max_err:.2e} Ha (diag-only: {diag_err:.2e}), "
                  f"rms={rms_err:.2e}")

            # Degeneracy pattern
            groups = analyze_degeneracies(evals, tol=1e-8)
            deg_str = ", ".join(f"{v:.8f}(×{c})" for v, c in groups[:8])
            print(f"      degens: {deg_str}")

            results[f'dirac_n{n_max}_{tag}'] = {
                'evals': evals.tolist(),
                'max_err_vs_dirac': float(max_err),
                'rms_err_vs_dirac': float(rms_err),
                'diag_only_err': float(diag_err),
            }

        # Exact Dirac-Coulomb degeneracies for reference
        exact_groups = analyze_degeneracies(d['exact_dirac'], tol=1e-10)
        deg_str = ", ".join(f"{v:.8f}(×{c})" for v, c in exact_groups[:8])
        print(f"    EXACT Dirac-Coulomb degens: {deg_str}")

        results[f'dirac_n{n_max}_exact'] = {
            'exact_dirac': d['exact_dirac'].tolist(),
        }

    # ---- Key analysis: graph perturbation vs fine-structure scale ----
    print("\n--- SCALE COMPARISON ---")
    print(f"  kappa = {KAPPA} = 1/{1/abs(KAPPA):.0f}")
    print(f"  alpha^2 = {ALPHA**2:.4e}")
    print(f"  FS splitting (n=2, Z=1): {abs(fine_structure_shift(2, -2, 1) - fine_structure_shift(2, 1, 1)):.4e} Ha")
    print(f"  Rydberg gap (n=1→2):     {abs(-0.5 + 0.125):.4f} Ha")

    # ---- Per-κ-sector analysis ----
    print("\n--- PER-κ SECTOR ANALYSIS (n_max=3) ---")
    adj, labels, _deg, _desc = build_dirac_s3_graph(3, adjacency_rule='A')

    kappa_values = sorted(set(lab.kappa for lab in labels))
    for kv in kappa_values:
        # Extract sector
        indices = [i for i, lab in enumerate(labels) if lab.kappa == kv]
        sector_labels = [labels[i] for i in indices]
        n_vals = sorted(set(lab.n_fock for lab in sector_labels))
        mj_vals = sorted(set(lab.two_m_j for lab in sector_labels))

        # Sector adjacency
        idx_map = {old: new for new, old in enumerate(indices)}
        sector_adj = np.zeros((len(indices), len(indices)))
        for i_new, i_old in enumerate(indices):
            for j_new, j_old in enumerate(indices):
                sector_adj[i_new, j_new] = adj[i_old, j_old]

        # Sector Hamiltonian (full FS diagonal)
        l_val = kappa_to_l(kv)
        j_val = abs(kv) - 0.5
        def diag_fs_sector(lab):
            return -Z**2 / (2.0 * lab.n_fock**2) + fine_structure_shift(lab.n_fock, lab.kappa, Z)

        H_sec = build_hamiltonian(sector_adj, sector_labels, diag_fs_sector)
        evals_sec = np.sort(np.linalg.eigvalsh(H_sec))

        # Exact Dirac-Coulomb for this sector
        exact_sec = np.sort(np.array([dirac_coulomb_energy(lab.n_fock, lab.kappa, Z)
                                       for lab in sector_labels]))

        max_err_sec = np.max(np.abs(evals_sec - exact_sec))
        print(f"  kappa={kv:+d} (l={l_val}, j={j_val}): {len(indices)} nodes, "
              f"n in {n_vals}, |m_j|<={j_val}, max|dE|={max_err_sec:.2e} Ha")

        # Show eigenvalue comparison
        for i in range(min(6, len(evals_sec))):
            lab = sector_labels[i]  # approximate — not exactly matched
            print(f"    graph: {evals_sec[i]:.10f}  exact: {exact_sec[i]:.10f}  "
                  f"delta={evals_sec[i]-exact_sec[i]:+.2e}")

    # Save results
    out_path = Path(__file__).parent / 'data' / 'dirac_native_graph_hamiltonian.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
