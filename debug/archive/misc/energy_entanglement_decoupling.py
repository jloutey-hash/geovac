"""
Energy-entanglement decoupling in the graph-native CI basis.

Investigates the quantitative relationship between the graph topology's
contribution to correlation energy vs entanglement entropy, and connects
this to the structural sparsity of GeoVac qubit Hamiltonians.

Three parts:
1. Isoelectronic sequence: fraction of E_corr and S from h1 vs V_ee across Z
2. Basis rotation test: does rotating away from graph basis destroy sparsity?
3. Effective rank as sparsity predictor: r_eff vs ERI count vs Z

Author: GeoVac Development Team
Date: April 2026
"""

import sys
import os
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.casimir_ci import (
    build_graph_native_fci, _build_graph_h1,
    build_fci_matrix, build_fci_polynomial, solve_variational_fast,
    _build_orbital_basis, two_electron_integral,
)

from debug.entanglement_geometry import (
    build_1rdm_from_singlet_ci,
    compute_entanglement_measures,
)


# =============================================================================
# Shared helpers
# =============================================================================

def build_decomposed_hamiltonians(Z_float, n_max):
    """Build decomposed Hamiltonians: h1_diag, h1_offdiag, V_ee_diag, V_ee_offdiag.

    Returns config-space matrices and auxiliary data.
    """
    from geovac.lattice import GeometricLattice

    lattice = GeometricLattice(max_n=n_max)
    orbitals = list(lattice.states)
    n_spatial = lattice.num_states

    # Build full h1 in orbital space
    h1_mat = np.zeros((n_spatial, n_spatial))
    for i, (n, l, m) in enumerate(orbitals):
        h1_mat[i, i] = -Z_float * Z_float / (2.0 * n * n)
    kappa = -1.0 / 16.0
    A = lattice.adjacency
    if hasattr(A, 'toarray'):
        A_dense = A.toarray()
    else:
        A_dense = np.array(A)
    for i in range(n_spatial):
        for j in range(n_spatial):
            if i != j and abs(A_dense[i, j]) > 1e-15:
                h1_mat[i, j] = kappa * (-A_dense[i, j])

    h1_diag_orb = np.diag(np.diag(h1_mat))
    h1_offdiag_orb = h1_mat - h1_diag_orb

    # Build m_total=0 singlet configs
    configs = []
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if orbitals[i][2] + orbitals[j][2] == 0:
                configs.append((i, j))
    n_configs = len(configs)
    k_orb = Z_float

    def build_config_matrix(h1_src, include_vee=False, vee_diag_only=False):
        """Build config-space matrix from orbital-space h1 and optionally V_ee."""
        H = np.zeros((n_configs, n_configs))
        for I, (i, j) in enumerate(configs):
            for J in range(I, n_configs):
                p, q = configs[J]
                I_perms = [(i, j)]
                if i != j:
                    I_perms.append((j, i))
                J_perms = [(p, q)]
                if p != q:
                    J_perms.append((q, p))
                N_I = np.sqrt(float(len(I_perms)))
                N_J = np.sqrt(float(len(J_perms)))
                me = 0.0
                for a, b in I_perms:
                    for c, d in J_perms:
                        if h1_src is not None:
                            if b == d:
                                me += h1_src[a, c]
                            if a == c:
                                me += h1_src[b, d]
                        if include_vee:
                            na, la, ma_ = orbitals[a]
                            nb, lb, mb_ = orbitals[b]
                            nc, lc, mc_ = orbitals[c]
                            nd, ld, md_ = orbitals[d]
                            me += two_electron_integral(
                                na, la, ma_, nb, lb, mb_,
                                nc, lc, mc_, nd, ld, md_, k_orb)
                me /= (N_I * N_J)
                H[I, J] = me
                H[J, I] = me
        return H

    H_h1_diag = build_config_matrix(h1_diag_orb, include_vee=False)
    H_h1_offdiag = build_config_matrix(h1_offdiag_orb, include_vee=False)
    H_vee_full = build_config_matrix(None, include_vee=True)
    H_vee_diag = np.diag(np.diag(H_vee_full))
    H_vee_offdiag = H_vee_full - H_vee_diag

    H_full = H_h1_diag + H_h1_offdiag + H_vee_full

    return {
        'H_full': H_full,
        'H_h1_diag': H_h1_diag,
        'H_h1_offdiag': H_h1_offdiag,
        'H_vee_full': H_vee_full,
        'H_vee_diag': H_vee_diag,
        'H_vee_offdiag': H_vee_offdiag,
        'configs': configs,
        'orbitals': orbitals,
        'n_spatial': n_spatial,
        'n_configs': n_configs,
        'h1_mat': h1_mat,
        'h1_diag_orb': h1_diag_orb,
        'h1_offdiag_orb': h1_offdiag_orb,
    }


def solve_and_entangle(H, configs, n_spatial):
    """Solve FCI and compute entanglement measures."""
    evals, evecs = np.linalg.eigh(H)
    E0 = evals[0]
    ci = evecs[:, 0]
    rho = build_1rdm_from_singlet_ci(ci, configs, n_spatial)
    ent = compute_entanglement_measures(rho)
    return E0, ci, rho, ent


# =============================================================================
# Part 1: Isoelectronic sequence — energy-entanglement decoupling
# =============================================================================

def part1_isoelectronic_decoupling(n_max=4):
    """Quantify energy vs entanglement fractions from h1 vs V_ee across Z."""
    print("\n" + "=" * 70)
    print(f"PART 1: ENERGY-ENTANGLEMENT DECOUPLING (n_max={n_max})")
    print("=" * 70)

    Z_values = [2, 3, 4, 5, 6, 8, 10]

    # Exact non-relativistic energies for He-like ions
    E_exact = {
        2: -2.903724377,
        3: -7.279913413,
        4: -13.65556624,
        5: -22.03097159,
        6: -32.40624671,
        8: -59.15696,
        10: -93.9068,
    }

    results = {}

    for Z in Z_values:
        t0 = time.time()
        print(f"\n  Z = {Z}...")

        decomp = build_decomposed_hamiltonians(float(Z), n_max)
        H_full = decomp['H_full']
        configs = decomp['configs']
        n_spatial = decomp['n_spatial']

        # Full solution
        E_full, ci_full, rho_full, ent_full = solve_and_entangle(
            H_full, configs, n_spatial)
        S_full = ent_full['von_neumann_entropy']

        # HF energy = energy of the 1s^2 configuration (first config = (0,0))
        # This is the diagonal element H[0,0] for the (1s,1s) config
        E_HF = H_full[0, 0]
        E_corr_full = E_full - E_HF

        # --- Energy decomposition ---
        # Expectation values of each component in the full ground state
        E_h1_diag = ci_full @ decomp['H_h1_diag'] @ ci_full
        E_h1_offdiag = ci_full @ decomp['H_h1_offdiag'] @ ci_full
        E_vee_diag = ci_full @ decomp['H_vee_diag'] @ ci_full
        E_vee_offdiag = ci_full @ decomp['H_vee_offdiag'] @ ci_full

        # HF has all weight on first config, so HF energies are:
        E_HF_h1_diag = decomp['H_h1_diag'][0, 0]
        E_HF_h1_offdiag = decomp['H_h1_offdiag'][0, 0]  # should be 0 for (1s,1s)
        E_HF_vee_diag = decomp['H_vee_diag'][0, 0]
        E_HF_vee_offdiag = decomp['H_vee_offdiag'][0, 0]  # should be 0

        # Correlation energy from each component
        dE_h1_diag = E_h1_diag - E_HF_h1_diag
        dE_h1_offdiag = E_h1_offdiag - E_HF_h1_offdiag
        dE_vee_diag = E_vee_diag - E_HF_vee_diag
        dE_vee_offdiag = E_vee_offdiag - E_HF_vee_offdiag

        # Fractions
        if abs(E_corr_full) > 1e-15:
            frac_E_h1_offdiag = dE_h1_offdiag / E_corr_full
            frac_E_vee_offdiag = dE_vee_offdiag / E_corr_full
            frac_E_h1_diag = dE_h1_diag / E_corr_full
            frac_E_vee_diag = dE_vee_diag / E_corr_full
        else:
            frac_E_h1_offdiag = frac_E_vee_offdiag = 0.0
            frac_E_h1_diag = frac_E_vee_diag = 0.0

        # --- Entanglement decomposition ---
        # Solve with subsets of the Hamiltonian to get entanglement from each layer

        # Layer A: h1_diag + V_ee_diag only (no off-diagonal anything)
        H_layer_A = decomp['H_h1_diag'] + decomp['H_vee_diag']
        _, _, _, ent_A = solve_and_entangle(H_layer_A, configs, n_spatial)
        S_A = ent_A['von_neumann_entropy']

        # Layer B: h1_diag + h1_offdiag + V_ee_diag (add graph topology)
        H_layer_B = decomp['H_h1_diag'] + decomp['H_h1_offdiag'] + decomp['H_vee_diag']
        _, _, _, ent_B = solve_and_entangle(H_layer_B, configs, n_spatial)
        S_B = ent_B['von_neumann_entropy']

        # Layer C: full (add V_ee off-diagonal)
        S_C = S_full

        # Entanglement fractions (incremental)
        if S_full > 1e-15:
            frac_S_rational = S_A / S_full
            frac_S_h1_offdiag = (S_B - S_A) / S_full
            frac_S_vee_offdiag = (S_C - S_B) / S_full
        else:
            frac_S_rational = frac_S_h1_offdiag = frac_S_vee_offdiag = 0.0

        elapsed = time.time() - t0

        error_pct = (E_full - E_exact[Z]) / abs(E_exact[Z]) * 100

        print(f"    E_full = {E_full:.6f} Ha (error = {error_pct:.2f}%)")
        print(f"    E_HF = {E_HF:.6f}, E_corr = {E_corr_full:.6f} Ha")
        print(f"    S_full = {S_full:.6f}")
        print(f"    Energy fractions: h1_offdiag={frac_E_h1_offdiag:.3f}, "
              f"Vee_offdiag={frac_E_vee_offdiag:.3f}, "
              f"h1_diag={frac_E_h1_diag:.3f}, Vee_diag={frac_E_vee_diag:.3f}")
        print(f"    Entropy fractions: h1_offdiag={frac_S_h1_offdiag:.4f}, "
              f"Vee_offdiag={frac_S_vee_offdiag:.4f}, "
              f"rational={frac_S_rational:.4f}")
        print(f"    Time: {elapsed:.1f}s")

        results[str(Z)] = {
            'Z': Z,
            'E_full': float(E_full),
            'E_HF': float(E_HF),
            'E_corr_full': float(E_corr_full),
            'error_pct': float(error_pct),
            'S_full': float(S_full),
            'S_A_rational': float(S_A),
            'S_B_topological': float(S_B),
            # Energy decomposition (expectation values)
            'E_h1_diag': float(E_h1_diag),
            'E_h1_offdiag': float(E_h1_offdiag),
            'E_vee_diag': float(E_vee_diag),
            'E_vee_offdiag': float(E_vee_offdiag),
            # Correlation energy from each component
            'dE_h1_diag': float(dE_h1_diag),
            'dE_h1_offdiag': float(dE_h1_offdiag),
            'dE_vee_diag': float(dE_vee_diag),
            'dE_vee_offdiag': float(dE_vee_offdiag),
            # Fractions
            'frac_E_h1_offdiag': float(frac_E_h1_offdiag),
            'frac_E_vee_offdiag': float(frac_E_vee_offdiag),
            'frac_E_h1_diag': float(frac_E_h1_diag),
            'frac_E_vee_diag': float(frac_E_vee_diag),
            'frac_S_rational': float(frac_S_rational),
            'frac_S_h1_offdiag': float(frac_S_h1_offdiag),
            'frac_S_vee_offdiag': float(frac_S_vee_offdiag),
            # Occupation numbers
            'occupation_numbers': [float(x) for x in ent_full['occupation_numbers'][:15]],
            'time_s': float(elapsed),
        }

    return results


# =============================================================================
# Part 2: Basis rotation test — does rotating away destroy sparsity?
# =============================================================================

def build_orbital_integrals(Z_float, n_max):
    """Build 1e and 2e integral arrays in the orbital basis.

    Returns h1[p,q] and eri[p,q,r,s] arrays.
    """
    from geovac.lattice import GeometricLattice

    lattice = GeometricLattice(max_n=n_max)
    orbitals = list(lattice.states)
    n_spatial = lattice.num_states

    # h1 matrix (graph-native)
    h1 = np.zeros((n_spatial, n_spatial))
    for i, (n, l, m) in enumerate(orbitals):
        h1[i, i] = -Z_float * Z_float / (2.0 * n * n)
    kappa = -1.0 / 16.0
    A = lattice.adjacency
    if hasattr(A, 'toarray'):
        A_dense = A.toarray()
    else:
        A_dense = np.array(A)
    for i in range(n_spatial):
        for j in range(n_spatial):
            if i != j and abs(A_dense[i, j]) > 1e-15:
                h1[i, j] = kappa * (-A_dense[i, j])

    # 2e integrals: eri[p,q,r,s] = <pq|rs> (physics convention)
    k_orb = Z_float
    eri = np.zeros((n_spatial, n_spatial, n_spatial, n_spatial))
    for p in range(n_spatial):
        for q in range(n_spatial):
            for r in range(n_spatial):
                for s in range(n_spatial):
                    np_, lp, mp = orbitals[p]
                    nq, lq, mq = orbitals[q]
                    nr, lr, mr = orbitals[r]
                    ns, ls, ms = orbitals[s]
                    eri[p, q, r, s] = two_electron_integral(
                        np_, lp, mp, nq, lq, mq,
                        nr, lr, mr, ns, ls, ms, k_orb)

    return h1, eri, orbitals, n_spatial


def rotate_integrals(h1, eri, U):
    """Rotate 1e and 2e integrals by unitary U.

    h1'[p,q] = sum_ab U[p,a] h1[a,b] U*[q,b]
    eri'[p,q,r,s] = sum_abcd U[p,a] U[q,b] eri[a,b,c,d] U*[r,c] U*[s,d]
    """
    # h1 rotation
    h1_rot = U @ h1 @ U.T  # U is real orthogonal

    # eri rotation via einsum (feasible for small n_spatial)
    eri_rot = np.einsum('pa,qb,abcd,rc,sd->pqrs', U, U, eri, U, U)

    return h1_rot, eri_rot


def count_significant(arr, threshold=1e-10):
    """Count elements with |value| > threshold."""
    return int(np.sum(np.abs(arr) > threshold))


def part2_basis_rotation(n_max=2, Z=2, n_rotations=20, n_angles=15):
    """Test how basis rotation affects integral sparsity."""
    print("\n" + "=" * 70)
    print(f"PART 2: BASIS ROTATION SPARSITY TEST (n_max={n_max}, Z={Z})")
    print("=" * 70)

    t0 = time.time()

    h1, eri, orbitals, n_spatial = build_orbital_integrals(float(Z), n_max)

    # Count significant integrals in graph basis
    n_h1_graph = count_significant(h1)
    n_eri_graph = count_significant(eri)
    total_h1 = n_spatial * n_spatial
    total_eri = n_spatial ** 4

    print(f"\n  n_spatial = {n_spatial}")
    print(f"  Graph basis: h1 significant = {n_h1_graph}/{total_h1} "
          f"({n_h1_graph/total_h1*100:.1f}%)")
    print(f"  Graph basis: ERI significant = {n_eri_graph}/{total_eri} "
          f"({n_eri_graph/total_eri*100:.1f}%)")

    # Generate random rotation generators (antisymmetric matrices)
    np.random.seed(42)
    generators = []
    for _ in range(n_rotations):
        G = np.random.randn(n_spatial, n_spatial)
        G = (G - G.T) / 2.0  # antisymmetric
        G /= np.linalg.norm(G)  # normalize
        generators.append(G)

    # Sweep rotation angle theta from 0 to pi/2
    thetas = np.linspace(0, np.pi / 2, n_angles)

    rotation_results = []

    for i_gen, G in enumerate(generators):
        print(f"\n  Rotation {i_gen + 1}/{n_rotations}...", end="", flush=True)
        for theta in thetas:
            U = scipy_matrix_exp(theta * G)
            h1_rot, eri_rot = rotate_integrals(h1, eri, U)
            n_h1_rot = count_significant(h1_rot)
            n_eri_rot = count_significant(eri_rot)

            rotation_results.append({
                'rotation_idx': i_gen,
                'theta': float(theta),
                'n_h1_significant': n_h1_rot,
                'n_eri_significant': n_eri_rot,
                'h1_density': float(n_h1_rot / total_h1),
                'eri_density': float(n_eri_rot / total_eri),
            })
        print(" done.", flush=True)

    # Compute statistics at each angle
    angle_stats = {}
    for theta in thetas:
        entries = [r for r in rotation_results if abs(r['theta'] - theta) < 1e-10]
        h1_counts = [r['n_h1_significant'] for r in entries]
        eri_counts = [r['n_eri_significant'] for r in entries]
        angle_stats[f"{theta:.4f}"] = {
            'theta': float(theta),
            'h1_mean': float(np.mean(h1_counts)),
            'h1_std': float(np.std(h1_counts)),
            'h1_min': int(np.min(h1_counts)),
            'h1_max': int(np.max(h1_counts)),
            'eri_mean': float(np.mean(eri_counts)),
            'eri_std': float(np.std(eri_counts)),
            'eri_min': int(np.min(eri_counts)),
            'eri_max': int(np.max(eri_counts)),
        }

    elapsed = time.time() - t0
    print(f"\n  Total time: {elapsed:.1f}s")

    # Print summary
    print(f"\n  Angle-dependent sparsity (ERI):")
    print(f"    theta=0.00: ERI = {n_eri_graph} (graph basis)")
    for key in sorted(angle_stats.keys(), key=lambda x: float(x)):
        s = angle_stats[key]
        if s['theta'] > 0.01:
            print(f"    theta={s['theta']:.2f}: ERI = {s['eri_mean']:.0f} "
                  f"+/- {s['eri_std']:.0f} "
                  f"(ratio = {s['eri_mean']/n_eri_graph:.2f}x)")

    results = {
        'Z': Z,
        'n_max': n_max,
        'n_spatial': n_spatial,
        'n_h1_graph': n_h1_graph,
        'n_eri_graph': n_eri_graph,
        'total_h1': total_h1,
        'total_eri': total_eri,
        'n_rotations': n_rotations,
        'n_angles': n_angles,
        'angle_stats': angle_stats,
        'time_s': float(elapsed),
    }

    return results


def scipy_matrix_exp(M):
    """Matrix exponential for real matrices."""
    from scipy.linalg import expm
    return expm(M)


# =============================================================================
# Part 3: Effective rank as sparsity predictor
# =============================================================================

def part3_effective_rank(n_max=4):
    """Compute effective rank and active orbital count vs Z."""
    print("\n" + "=" * 70)
    print(f"PART 3: EFFECTIVE RANK AS SPARSITY PREDICTOR (n_max={n_max})")
    print("=" * 70)

    Z_values = [2, 3, 4, 5, 6, 8, 10]

    results = {}

    for Z in Z_values:
        t0 = time.time()
        print(f"\n  Z = {Z}...")

        decomp = build_decomposed_hamiltonians(float(Z), n_max)
        H_full = decomp['H_full']
        configs = decomp['configs']
        n_spatial = decomp['n_spatial']

        # Full solution
        E_full, ci_full, rho_full, ent_full = solve_and_entangle(
            H_full, configs, n_spatial)

        occ = ent_full['occupation_numbers']
        S = ent_full['von_neumann_entropy']

        # Effective rank: exp(S)
        r_eff = np.exp(S)

        # Active orbitals: occupation > threshold
        active_001 = int(np.sum(np.array(occ) > 0.001))
        active_01 = int(np.sum(np.array(occ) > 0.01))
        active_0001 = int(np.sum(np.array(occ) > 0.0001))

        # Count significant ERIs (number of non-zero Slater integrals)
        # We count unique ERIs with |value| > 1e-10
        k_orb = float(Z)
        orbitals = decomp['orbitals']
        n_eri_significant = 0
        n_eri_total = 0
        for p in range(n_spatial):
            for q in range(p, n_spatial):
                for r in range(n_spatial):
                    for s in range(r, n_spatial):
                        np_, lp, mp = orbitals[p]
                        nq, lq, mq = orbitals[q]
                        nr, lr, mr = orbitals[r]
                        ns, ls, ms = orbitals[s]
                        val = two_electron_integral(
                            np_, lp, mp, nq, lq, mq,
                            nr, lr, mr, ns, ls, ms, k_orb)
                        n_eri_total += 1
                        if abs(val) > 1e-10:
                            n_eri_significant += 1

        # Prediction: Pauli count ~ r_eff^2 (pairwise interactions)
        predicted_pauli_reff2 = r_eff ** 2

        elapsed = time.time() - t0

        print(f"    S = {S:.6f}, r_eff = exp(S) = {r_eff:.4f}")
        print(f"    Active orbitals: >0.001: {active_001}, >0.01: {active_01}, "
              f">0.0001: {active_0001}")
        print(f"    Significant ERIs: {n_eri_significant}/{n_eri_total} "
              f"({n_eri_significant/n_eri_total*100:.1f}%)")
        print(f"    Time: {elapsed:.1f}s")

        results[str(Z)] = {
            'Z': Z,
            'von_neumann_entropy': float(S),
            'effective_rank': float(r_eff),
            'active_001': active_001,
            'active_01': active_01,
            'active_0001': active_0001,
            'n_eri_significant': n_eri_significant,
            'n_eri_total': n_eri_total,
            'eri_density_pct': float(n_eri_significant / n_eri_total * 100),
            'predicted_pauli_reff2': float(predicted_pauli_reff2),
            'occupation_numbers': [float(x) for x in occ[:15]],
            'time_s': float(elapsed),
        }

    return results


# =============================================================================
# Plotting
# =============================================================================

def plot_energy_vs_entropy_fractions(results_p1, outdir):
    """Plot energy fraction vs entropy fraction for h1/V_ee across Z."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    Z_vals = []
    frac_E_h1 = []
    frac_E_vee = []
    frac_S_h1 = []
    frac_S_vee = []

    for key in sorted(results_p1.keys(), key=int):
        r = results_p1[key]
        Z_vals.append(r['Z'])
        frac_E_h1.append(r['frac_E_h1_offdiag'])
        frac_E_vee.append(r['frac_E_vee_offdiag'])
        frac_S_h1.append(r['frac_S_h1_offdiag'])
        frac_S_vee.append(r['frac_S_vee_offdiag'])

    Z_vals = np.array(Z_vals)

    # Left panel: fractions vs Z
    ax = axes[0]
    ax.plot(Z_vals, frac_E_h1, 'bo-', label='E_corr from h1_offdiag', linewidth=2, markersize=8)
    ax.plot(Z_vals, frac_E_vee, 'rs-', label='E_corr from V_ee_offdiag', linewidth=2, markersize=8)
    ax.plot(Z_vals, frac_S_h1, 'b^--', label='S from h1_offdiag', linewidth=2, markersize=8)
    ax.plot(Z_vals, frac_S_vee, 'rv--', label='S from V_ee_offdiag', linewidth=2, markersize=8)
    ax.set_xlabel('Nuclear charge Z', fontsize=13)
    ax.set_ylabel('Fraction of total', fontsize=13)
    ax.set_title('Energy-Entanglement Decoupling\nacross Isoelectronic Sequence', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axhline(y=1, color='k', linewidth=0.5, linestyle=':')

    # Right panel: scatter of E fraction vs S fraction (each point is a Z)
    ax = axes[1]
    ax.scatter(frac_E_h1, frac_S_h1, c='blue', s=100, zorder=5,
               label='h1_offdiag', edgecolors='black')
    ax.scatter(frac_E_vee, frac_S_vee, c='red', s=100, zorder=5,
               label='V_ee_offdiag', edgecolors='black')
    for i, Z in enumerate(Z_vals):
        ax.annotate(f'Z={int(Z)}', (frac_E_h1[i], frac_S_h1[i]),
                    textcoords="offset points", xytext=(5, 5), fontsize=9, color='blue')
        ax.annotate(f'Z={int(Z)}', (frac_E_vee[i], frac_S_vee[i]),
                    textcoords="offset points", xytext=(5, -10), fontsize=9, color='red')

    ax.set_xlabel('Fraction of correlation energy', fontsize=13)
    ax.set_ylabel('Fraction of entanglement entropy', fontsize=13)
    ax.set_title('Energy vs Entropy Contribution\n(anti-correlation = decoupling)', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    # Draw diagonal for reference
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='_no_legend_')

    plt.tight_layout()
    path = os.path.join(outdir, 'energy_vs_entropy_fractions.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved: {path}")


def plot_integral_sparsity(results_p2, outdir):
    """Plot integral count vs basis rotation angle."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    angle_stats = results_p2['angle_stats']
    thetas = []
    eri_means = []
    eri_stds = []
    h1_means = []
    h1_stds = []

    for key in sorted(angle_stats.keys(), key=float):
        s = angle_stats[key]
        thetas.append(s['theta'])
        eri_means.append(s['eri_mean'])
        eri_stds.append(s['eri_std'])
        h1_means.append(s['h1_mean'])
        h1_stds.append(s['h1_std'])

    thetas = np.array(thetas)
    eri_means = np.array(eri_means)
    eri_stds = np.array(eri_stds)
    h1_means = np.array(h1_means)
    h1_stds = np.array(h1_stds)

    n_eri_graph = results_p2['n_eri_graph']
    n_h1_graph = results_p2['n_h1_graph']
    total_eri = results_p2['total_eri']
    total_h1 = results_p2['total_h1']

    # Left: ERI count vs angle
    ax = axes[0]
    ax.fill_between(thetas, eri_means - eri_stds, eri_means + eri_stds,
                    alpha=0.3, color='red')
    ax.plot(thetas, eri_means, 'r-o', linewidth=2, markersize=5, label='Rotated basis (mean)')
    ax.axhline(y=n_eri_graph, color='blue', linewidth=2, linestyle='--',
               label=f'Graph basis ({n_eri_graph})')
    ax.axhline(y=total_eri, color='gray', linewidth=1, linestyle=':',
               label=f'Maximum ({total_eri})')
    ax.set_xlabel('Rotation angle (rad)', fontsize=13)
    ax.set_ylabel('Significant ERI count (|val| > 1e-10)', fontsize=13)
    ax.set_title(f'ERI Sparsity vs Basis Rotation\n(n_max={results_p2["n_max"]}, Z={results_p2["Z"]})',
                 fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Right: ratio to graph basis
    ax = axes[1]
    ratio_eri = eri_means / n_eri_graph
    ax.fill_between(thetas,
                    (eri_means - eri_stds) / n_eri_graph,
                    (eri_means + eri_stds) / n_eri_graph,
                    alpha=0.3, color='red')
    ax.plot(thetas, ratio_eri, 'r-o', linewidth=2, markersize=5, label='ERI ratio')

    ratio_h1 = h1_means / max(n_h1_graph, 1)
    ax.plot(thetas, ratio_h1, 'b-s', linewidth=2, markersize=5, label='h1 ratio')

    ax.axhline(y=1.0, color='k', linewidth=1, linestyle='--')
    ax.set_xlabel('Rotation angle (rad)', fontsize=13)
    ax.set_ylabel('Ratio to graph basis', fontsize=13)
    ax.set_title('Integral Count Ratio\n(>1 = graph basis is sparser)', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(outdir, 'integral_sparsity_vs_rotation.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_effective_rank(results_p3, outdir):
    """Plot effective rank and active orbital count vs Z."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    Z_vals = []
    r_effs = []
    active_001 = []
    active_01 = []
    active_0001 = []
    S_vals = []
    eri_densities = []
    n_eri_sig = []

    for key in sorted(results_p3.keys(), key=int):
        r = results_p3[key]
        Z_vals.append(r['Z'])
        r_effs.append(r['effective_rank'])
        active_001.append(r['active_001'])
        active_01.append(r['active_01'])
        active_0001.append(r['active_0001'])
        S_vals.append(r['von_neumann_entropy'])
        eri_densities.append(r['eri_density_pct'])
        n_eri_sig.append(r['n_eri_significant'])

    Z_vals = np.array(Z_vals)

    # Left: effective rank and active orbitals vs Z
    ax = axes[0]
    ax.plot(Z_vals, r_effs, 'ko-', linewidth=2, markersize=8, label='r_eff = exp(S)')
    ax.plot(Z_vals, active_001, 'rs--', linewidth=2, markersize=7, label='Active (n > 0.001)')
    ax.plot(Z_vals, active_0001, 'b^:', linewidth=2, markersize=7, label='Active (n > 0.0001)')
    ax.set_xlabel('Nuclear charge Z', fontsize=13)
    ax.set_ylabel('Count', fontsize=13)
    ax.set_title('Effective Rank vs Z', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Middle: entropy vs Z
    ax = axes[1]
    ax.plot(Z_vals, S_vals, 'go-', linewidth=2, markersize=8)
    ax.set_xlabel('Nuclear charge Z', fontsize=13)
    ax.set_ylabel('Von Neumann entropy S', fontsize=13)
    ax.set_title('Entanglement Entropy vs Z', fontsize=14)
    ax.grid(True, alpha=0.3)

    # Fit power law S ~ Z^alpha
    log_Z = np.log(Z_vals)
    log_S = np.log(S_vals)
    coeffs = np.polyfit(log_Z, log_S, 1)
    ax.text(0.05, 0.95, f'S ~ Z^{{{coeffs[0]:.2f}}}',
            transform=ax.transAxes, fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Right: ERI density vs Z
    ax = axes[2]
    ax.plot(Z_vals, n_eri_sig, 'mo-', linewidth=2, markersize=8, label='Significant ERIs')
    ax.set_xlabel('Nuclear charge Z', fontsize=13)
    ax.set_ylabel('Count', fontsize=13)
    ax.set_title('Significant ERI Count vs Z\n(n_max=4, threshold=1e-10)', fontsize=14)
    ax.grid(True, alpha=0.3)
    # Note: ERI count is Z-independent (Gaunt selection rules)
    ax.text(0.05, 0.95, 'ERI count is Z-independent\n(Gaunt selection rules)',
            transform=ax.transAxes, fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    plt.tight_layout()
    path = os.path.join(outdir, 'effective_rank_vs_Z.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# =============================================================================
# Main
# =============================================================================

def main():
    outdir = os.path.join(os.path.dirname(__file__), 'plots')
    datadir = os.path.join(os.path.dirname(__file__), 'data')
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(datadir, exist_ok=True)

    all_results = {}

    # Part 1: Isoelectronic decoupling
    print("\n" + "#" * 70)
    print("# PART 1: ISOELECTRONIC SEQUENCE")
    print("#" * 70)
    results_p1 = part1_isoelectronic_decoupling(n_max=4)
    all_results['part1_isoelectronic'] = results_p1

    # Part 2: Basis rotation sparsity
    print("\n" + "#" * 70)
    print("# PART 2: BASIS ROTATION TEST")
    print("#" * 70)
    results_p2 = part2_basis_rotation(n_max=2, Z=2, n_rotations=20, n_angles=15)
    all_results['part2_rotation'] = results_p2

    # Part 3: Effective rank
    print("\n" + "#" * 70)
    print("# PART 3: EFFECTIVE RANK")
    print("#" * 70)
    results_p3 = part3_effective_rank(n_max=4)
    all_results['part3_effective_rank'] = results_p3

    # Save results
    data_path = os.path.join(datadir, 'energy_entanglement_decoupling.json')
    with open(data_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\n  Results saved to: {data_path}")

    # Generate plots
    print("\n" + "#" * 70)
    print("# GENERATING PLOTS")
    print("#" * 70)
    plot_energy_vs_entropy_fractions(results_p1, outdir)
    plot_integral_sparsity(results_p2, outdir)
    plot_effective_rank(results_p3, outdir)

    # Final summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print("\n  Part 1 — Energy-Entanglement Decoupling:")
    print("  " + "-" * 60)
    print(f"  {'Z':>3s}  {'frac_E(h1)':>10s}  {'frac_E(Vee)':>11s}  "
          f"{'frac_S(h1)':>10s}  {'frac_S(Vee)':>11s}  {'S':>8s}")
    for key in sorted(results_p1.keys(), key=int):
        r = results_p1[key]
        print(f"  {r['Z']:3d}  {r['frac_E_h1_offdiag']:10.4f}  "
              f"{r['frac_E_vee_offdiag']:11.4f}  "
              f"{r['frac_S_h1_offdiag']:10.4f}  "
              f"{r['frac_S_vee_offdiag']:11.4f}  "
              f"{r['S_full']:8.5f}")

    print(f"\n  Part 2 — Basis Rotation:")
    print(f"    Graph basis ERI count: {results_p2['n_eri_graph']}")
    angle_stats = results_p2['angle_stats']
    max_theta = max(float(k) for k in angle_stats.keys())
    max_stats = angle_stats[f"{max_theta:.4f}"]
    print(f"    Max rotation ERI count: {max_stats['eri_mean']:.0f} "
          f"(ratio = {max_stats['eri_mean']/results_p2['n_eri_graph']:.2f}x)")

    print(f"\n  Part 3 — Effective Rank:")
    print(f"  {'Z':>3s}  {'r_eff':>6s}  {'Active(>0.001)':>14s}  "
          f"{'ERI_sig':>7s}  {'ERI%':>6s}")
    for key in sorted(results_p3.keys(), key=int):
        r = results_p3[key]
        print(f"  {r['Z']:3d}  {r['effective_rank']:6.3f}  "
              f"{r['active_001']:14d}  "
              f"{r['n_eri_significant']:7d}  "
              f"{r['eri_density_pct']:6.1f}")


if __name__ == '__main__':
    main()
