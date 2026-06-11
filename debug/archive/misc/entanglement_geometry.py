"""
Entanglement geometry of the He ground state in the graph-native CI basis.

Computes:
1. Ground state CI coefficients at n_max = 2..7
2. 1-electron reduced density matrix (1-RDM) from the spin-adapted singlet CI
3. Entanglement measures:
   a. Von Neumann entanglement entropy S = -Tr(rho log rho)
   b. Single-orbital entanglement entropy s_i
   c. Entanglement spectrum (occupation numbers)
4. Scaling analysis: S vs n_max, per-shell entanglement

Also compares graph-native CI with variational-k CI.

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


# =============================================================================
# 1-RDM construction from singlet spin-adapted CI
# =============================================================================

def build_1rdm_from_singlet_ci(
    ci_coeffs: np.ndarray,
    configs: list,
    n_spatial: int,
) -> np.ndarray:
    """Build the spatial 1-RDM from singlet spin-adapted CI coefficients.

    The singlet CI wavefunction is:
        |Psi> = sum_I c_I |Phi_I>
    where each |Phi_I> = [phi_i(1)*phi_j(2) + phi_j(1)*phi_i(2)] / N_IJ
    with N_IJ = sqrt(2) for i!=j, N_IJ = 1 for i==j.

    The spatial 1-RDM is:
        rho_{pq} = <Psi| (sum over electron) |p><q| |Psi>
               = 2 * sum_I,J c_I c_J * <Phi_I| |p(1)><q(1)| |Phi_J>
    (factor 2 for tracing over electron 2 in the singlet).

    For singlet CI, the 1-RDM can be computed as:
        rho_{pq} = sum_I,J c_I c_J * Gamma_{IJ,pq}

    where Gamma_{IJ,pq} accounts for the spatial overlaps after tracing
    out one electron.

    Parameters
    ----------
    ci_coeffs : array of shape (n_configs,)
        CI expansion coefficients.
    configs : list of (i, j) pairs
        Configuration list (spatial orbital indices, i <= j).
    n_spatial : int
        Number of spatial orbitals.

    Returns
    -------
    rho : array of shape (n_spatial, n_spatial)
        Spatial 1-RDM. Trace should equal 2 (two electrons).
    """
    rho = np.zeros((n_spatial, n_spatial))

    for I, (i, j) in enumerate(configs):
        for J, (p, q) in enumerate(configs):
            coeff = ci_coeffs[I] * ci_coeffs[J]
            if abs(coeff) < 1e-16:
                continue

            N_I = np.sqrt(2.0) if i != j else 1.0
            N_J = np.sqrt(2.0) if p != q else 1.0

            # The singlet spatial function for config I is:
            # |I> = [|i,j> + |j,i>] / N_I
            # Trace out electron 2:
            # rho1_{pq}^{IJ} = sum_r <I|p,r><q,r|J>
            #   = (1/N_I N_J) sum_r [
            #       <i|p><j|r> <q|p'><r|q'> + permutations
            #     ]
            # where I=(i,j), J=(p',q').

            # Explicitly: |I> has bra terms |i>|j> and |j>|i> (if i!=j)
            # |J> has ket terms |p>|q> and |q>|p> (if p!=q)
            # rho_{alpha,beta} = sum over electron-2 index gamma:
            #   sum_{(a,b) in I_perms} sum_{(c,d) in J_perms}
            #     <a|alpha> <b|gamma> <gamma|d> <beta|c>
            #   = sum_{(a,b),(c,d)} delta(a,alpha) delta(b,d) delta(beta,c)
            #   ... wait, this is wrong. Let me be more careful.

            # <Phi_I| alpha(1)><beta(1) |Phi_J>
            # = (1/N_I N_J) sum_{(a,b) in I} sum_{(c,d) in J}
            #   <a(1)b(2)| |alpha(1)><beta(1)| |c(1)d(2)>
            # = (1/N_I N_J) sum_{a,b,c,d} <a|alpha><beta|c> <b|d>
            # = (1/N_I N_J) sum_{a,b,c,d} delta(a,alpha)*delta(c,beta)*delta(b,d)

            # I permutations
            I_perms = [(i, j)]
            if i != j:
                I_perms.append((j, i))
            J_perms = [(p, q)]
            if p != q:
                J_perms.append((q, p))

            for a, b in I_perms:
                for c, d in J_perms:
                    if b == d:
                        # Contributes to rho[a, c]
                        rho[a, c] += coeff / (N_I * N_J)

    # Multiply by 2 for spin: each spatial orbital is occupied by
    # alpha+beta spin in a singlet, so the spin-traced spatial 1-RDM
    # has Tr = N_electrons = 2.
    rho *= 2.0

    return rho


def compute_entanglement_measures(rho: np.ndarray):
    """Compute entanglement measures from the spatial 1-RDM.

    For a 2-electron singlet, the spin-orbital 1-RDM is:
        rho_spin = rho_spatial / 2  (each spin gets half)
    and the occupation numbers of the spatial 1-RDM sum to 2.

    We define entanglement using the spatial natural orbital occupation
    numbers n_i (eigenvalues of rho, summing to 2). The von Neumann
    entropy of the normalized 1-RDM (divided by N_e) measures orbital
    entanglement.

    Parameters
    ----------
    rho : array of shape (M, M)
        Spatial 1-RDM with Tr(rho) = 2.

    Returns
    -------
    dict with:
        'occupation_numbers': sorted eigenvalues of rho (descending)
        'von_neumann_entropy': S = -sum n_i/2 * log(n_i/2)
        'single_orbital_entropy': array of s_i for each natural orbital
        'total_correlation': sum of s_i
    """
    # Eigenvalues of the spatial 1-RDM
    occ, nat_orbs = np.linalg.eigh(rho)
    # Sort descending
    idx = np.argsort(occ)[::-1]
    occ = occ[idx]
    nat_orbs = nat_orbs[:, idx]

    # Clean up tiny negative eigenvalues from numerical noise
    occ = np.maximum(occ, 0.0)

    # Von Neumann entropy of the normalized 1-RDM
    # For 2 electrons: rho_norm = rho / 2, eigenvalues lambda_i = n_i / 2
    # S = -Tr(rho_norm log rho_norm) = -sum lambda_i log(lambda_i)
    lam = occ / 2.0
    S = 0.0
    for x in lam:
        if x > 1e-15:
            S -= x * np.log(x)

    # Single-orbital entanglement entropy
    # For each spin-orbital with occupation n (0 <= n <= 1):
    #   s = -n log(n) - (1-n) log(1-n)
    # But here we have spatial occupation numbers (0 <= n_i <= 2).
    # In the spin-orbital picture, each spatial orbital gives TWO
    # spin-orbitals with occupation n_i/2 each.
    # s_i = 2 * [-n_i/2 * log(n_i/2) - (1 - n_i/2) * log(1 - n_i/2)]
    s_orbital = np.zeros(len(occ))
    for k, n in enumerate(occ):
        p = n / 2.0  # spin-orbital occupation
        if p > 1e-15 and p < 1.0 - 1e-15:
            s_orbital[k] = 2.0 * (-p * np.log(p) - (1 - p) * np.log(1 - p))
        elif p > 1e-15:
            s_orbital[k] = 2.0 * (-p * np.log(p))
        # else: 0

    return {
        'occupation_numbers': occ,
        'natural_orbitals': nat_orbs,
        'von_neumann_entropy': S,
        'single_orbital_entropy': s_orbital,
        'total_correlation': np.sum(s_orbital),
    }


def assign_quantum_numbers_to_natural_orbitals(
    nat_orbs: np.ndarray,
    orbitals: list,
) -> list:
    """Assign approximate (n, l) labels to natural orbitals.

    Each natural orbital is a linear combination of the original (n,l,m)
    orbitals. We assign the dominant (n,l) character by finding which
    (n,l) shell has the largest total weight.

    Parameters
    ----------
    nat_orbs : array of shape (n_spatial, n_natural)
        Natural orbital coefficients in the (n,l,m) basis.
    orbitals : list of (n, l, m) tuples
        Original orbital labels.

    Returns
    -------
    list of dicts with 'dominant_n', 'dominant_l', 'shell_weights'
    """
    n_nat = nat_orbs.shape[1]
    assignments = []

    for k in range(n_nat):
        coeffs = nat_orbs[:, k]
        # Group weights by (n, l) shell
        shell_weights = {}
        for idx, (n, l, m) in enumerate(orbitals):
            key = (n, l)
            shell_weights[key] = shell_weights.get(key, 0.0) + coeffs[idx] ** 2

        # Find dominant shell
        dominant = max(shell_weights, key=shell_weights.get)
        assignments.append({
            'dominant_n': dominant[0],
            'dominant_l': dominant[1],
            'shell_weights': {f"({n},{l})": float(w)
                              for (n, l), w in sorted(shell_weights.items())
                              if w > 1e-10},
        })

    return assignments


# =============================================================================
# Main computation
# =============================================================================

def run_graph_native_analysis(n_max_values=None):
    """Run entanglement analysis for graph-native CI at multiple n_max."""
    if n_max_values is None:
        n_max_values = [2, 3, 4, 5, 6, 7]

    E_exact = -2.903724377034
    results = {}

    for n_max in n_max_values:
        print(f"\n{'='*60}")
        print(f"Graph-native CI: n_max = {n_max}")
        print(f"{'='*60}")

        t0 = time.time()

        # Build the FCI Hamiltonian
        H = build_graph_native_fci(Z=2, n_max=n_max, m_total=0, spin='singlet')
        n_configs = H.shape[0]

        # Solve for ground state
        evals, evecs = np.linalg.eigh(H)
        E0 = evals[0]
        ci_coeffs = evecs[:, 0]

        error_pct = abs((E0 - E_exact) / E_exact) * 100

        # Get orbital basis and configs
        h1_mat, orbitals = _build_graph_h1(Z=2, n_max=n_max)
        n_spatial = len(orbitals)

        configs = []
        for i in range(n_spatial):
            for j in range(i, n_spatial):
                if orbitals[i][2] + orbitals[j][2] == 0:
                    configs.append((i, j))

        assert len(configs) == n_configs, \
            f"Config count mismatch: {len(configs)} vs {n_configs}"

        # Build 1-RDM
        rho = build_1rdm_from_singlet_ci(ci_coeffs, configs, n_spatial)
        trace = np.trace(rho)
        print(f"  Tr(rho) = {trace:.10f} (should be 2.0)")

        # Compute entanglement measures
        ent = compute_entanglement_measures(rho)

        # Assign quantum numbers to natural orbitals
        assignments = assign_quantum_numbers_to_natural_orbitals(
            ent['natural_orbitals'], orbitals)

        t_elapsed = time.time() - t0

        # Per-shell entanglement
        shell_entropy = {}
        for k, (occ, s_i, assign) in enumerate(zip(
                ent['occupation_numbers'],
                ent['single_orbital_entropy'],
                assignments)):
            key = f"({assign['dominant_n']},{assign['dominant_l']})"
            if key not in shell_entropy:
                shell_entropy[key] = {'occ_sum': 0.0, 'entropy_sum': 0.0, 'count': 0}
            shell_entropy[key]['occ_sum'] += occ
            shell_entropy[key]['entropy_sum'] += s_i
            shell_entropy[key]['count'] += 1

        # Print results
        print(f"  n_spatial = {n_spatial}, n_configs = {n_configs}")
        print(f"  E0 = {E0:.8f} Ha, error = {error_pct:.4f}%")
        print(f"  Von Neumann entropy S = {ent['von_neumann_entropy']:.6f}")
        print(f"  Total correlation = {ent['total_correlation']:.6f}")
        print(f"  Time: {t_elapsed:.1f}s")

        print(f"\n  Top occupation numbers:")
        for k in range(min(10, len(ent['occupation_numbers']))):
            occ = ent['occupation_numbers'][k]
            s_i = ent['single_orbital_entropy'][k]
            a = assignments[k]
            print(f"    n_{k} = {occ:.8f}  s_{k} = {s_i:.6f}  "
                  f"dominant ({a['dominant_n']},{a['dominant_l']})")

        print(f"\n  Per-shell entanglement:")
        for key in sorted(shell_entropy.keys()):
            d = shell_entropy[key]
            print(f"    {key}: total_occ = {d['occ_sum']:.6f}, "
                  f"total_s = {d['entropy_sum']:.6f}, count = {d['count']}")

        # Dominant CI coefficients
        sorted_ci = np.argsort(np.abs(ci_coeffs))[::-1]
        print(f"\n  Top CI coefficients:")
        for rank, idx in enumerate(sorted_ci[:8]):
            i, j = configs[idx]
            oi, oj = orbitals[i], orbitals[j]
            print(f"    |{oi}><{oj}|  c = {ci_coeffs[idx]:.8f}  "
                  f"|c|^2 = {ci_coeffs[idx]**2:.6f}")

        results[n_max] = {
            'n_spatial': n_spatial,
            'n_configs': n_configs,
            'E0': float(E0),
            'error_pct': float(error_pct),
            'trace_rho': float(trace),
            'von_neumann_entropy': float(ent['von_neumann_entropy']),
            'total_correlation': float(ent['total_correlation']),
            'occupation_numbers': [float(x) for x in ent['occupation_numbers']],
            'single_orbital_entropy': [float(x) for x in ent['single_orbital_entropy']],
            'shell_entropy': {k: {kk: float(vv) if isinstance(vv, (float, np.floating))
                                  else vv for kk, vv in v.items()}
                              for k, v in shell_entropy.items()},
            'natural_orbital_assignments': assignments[:min(15, len(assignments))],
            'time_s': float(t_elapsed),
            'top_ci_coefficients': [
                {'config': [int(configs[idx][0]), int(configs[idx][1])],
                 'orbitals': [list(orbitals[configs[idx][0]]),
                              list(orbitals[configs[idx][1]])],
                 'coefficient': float(ci_coeffs[idx]),
                 'weight': float(ci_coeffs[idx]**2)}
                for idx in sorted_ci[:10]
            ],
        }

    return results


def run_variational_k_comparison(n_max_values=None):
    """Run entanglement analysis for variational-k CI for comparison."""
    if n_max_values is None:
        n_max_values = [2, 3, 4, 5]

    E_exact = -2.903724377034
    results = {}

    for n_max in n_max_values:
        print(f"\n{'='*60}")
        print(f"Variational-k CI: n_max = {n_max}")
        print(f"{'='*60}")

        t0 = time.time()

        # Get optimal k
        B, C = build_fci_polynomial(Z=2, n_max=n_max)
        k_opt, E_opt = solve_variational_fast(Z=2, n_max=n_max, B=B, C=C)

        # Build H at optimal k and solve
        H = B * k_opt + C * k_opt**2
        evals, evecs = np.linalg.eigh(H)
        E0 = evals[0]
        ci_coeffs = evecs[:, 0]
        error_pct = abs((E0 - E_exact) / E_exact) * 100

        # Get orbital basis
        orbitals = _build_orbital_basis(n_max)
        n_spatial = len(orbitals)

        configs = []
        for i in range(n_spatial):
            for j in range(i, n_spatial):
                if orbitals[i][2] + orbitals[j][2] == 0:
                    configs.append((i, j))

        # Build 1-RDM
        rho = build_1rdm_from_singlet_ci(ci_coeffs, configs, n_spatial)
        ent = compute_entanglement_measures(rho)
        assignments = assign_quantum_numbers_to_natural_orbitals(
            ent['natural_orbitals'], orbitals)

        t_elapsed = time.time() - t0

        print(f"  k_opt = {k_opt:.6f}")
        print(f"  E0 = {E0:.8f} Ha, error = {error_pct:.4f}%")
        print(f"  Von Neumann entropy S = {ent['von_neumann_entropy']:.6f}")
        print(f"  Tr(rho) = {np.trace(rho):.10f}")

        results[n_max] = {
            'k_opt': float(k_opt),
            'n_spatial': n_spatial,
            'n_configs': len(configs),
            'E0': float(E0),
            'error_pct': float(error_pct),
            'von_neumann_entropy': float(ent['von_neumann_entropy']),
            'total_correlation': float(ent['total_correlation']),
            'occupation_numbers': [float(x) for x in ent['occupation_numbers']],
            'single_orbital_entropy': [float(x) for x in ent['single_orbital_entropy']],
            'time_s': float(t_elapsed),
        }

    return results


def analyze_three_layer_entanglement(n_max=4):
    """Decompose entanglement into the three-layer structure.

    Layer 1 (Rational): diagonal h1 + diagonal V_ee (mean field)
    Layer 2 (Topological): off-diagonal graph h1 (kappa = -1/16)
    Layer 3 (Transcendental): off-diagonal V_ee (correlation)

    For each layer, we compute the 1-RDM and entanglement if only that
    layer's H were used, and also the incremental entropy from adding
    each layer.
    """
    print(f"\n{'='*60}")
    print(f"Three-layer entanglement decomposition: n_max = {n_max}")
    print(f"{'='*60}")

    h1_mat, orbitals = _build_graph_h1(Z=2, n_max=n_max)
    n_spatial = len(orbitals)

    configs = []
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if orbitals[i][2] + orbitals[j][2] == 0:
                configs.append((i, j))
    n_configs = len(configs)

    k_orb = 2.0  # Z=2

    # Build full V_ee matrix in configuration space
    def build_vee_config_matrix():
        """Build V_ee in singlet configuration space."""
        V = np.zeros((n_configs, n_configs))
        for I, (i, j) in enumerate(configs):
            for J, (p, q) in enumerate(configs):
                if J < I:
                    continue
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
                        na, la, ma = orbitals[a]
                        nb, lb, mb = orbitals[b]
                        nc, lc, mc = orbitals[c]
                        nd, ld, md = orbitals[d]
                        me += two_electron_integral(
                            na, la, ma, nb, lb, mb,
                            nc, lc, mc, nd, ld, md, k_orb)
                me /= (N_I * N_J)
                V[I, J] = me
                V[J, I] = me
        return V

    # Build h1 (one-body) in configuration space
    def build_h1_config_matrix(h1_source):
        """Build one-body H in singlet configuration space."""
        H1 = np.zeros((n_configs, n_configs))
        for I, (i, j) in enumerate(configs):
            for J, (p, q) in enumerate(configs):
                if J < I:
                    continue
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
                        if b == d:
                            me += h1_source[a, c]
                        if a == c:
                            me += h1_source[b, d]
                me /= (N_I * N_J)
                H1[I, J] = me
                H1[J, I] = me
        return H1

    # Decompose h1 into diagonal and off-diagonal parts
    h1_diag = np.diag(np.diag(h1_mat))
    h1_offdiag = h1_mat - h1_diag

    # Build configuration-space matrices for each layer
    print("  Building configuration-space matrices...")
    H1_diag = build_h1_config_matrix(h1_diag)
    H1_offdiag = build_h1_config_matrix(h1_offdiag)
    V_ee = build_vee_config_matrix()

    # Split V_ee into diagonal (mean-field) and off-diagonal (correlation)
    V_ee_diag = np.diag(np.diag(V_ee))
    V_ee_offdiag = V_ee - V_ee_diag

    # Full Hamiltonian check
    H_full = H1_diag + H1_offdiag + V_ee
    H_ref = build_graph_native_fci(Z=2, n_max=n_max)
    print(f"  H reconstruction error: {np.max(np.abs(H_full - H_ref)):.2e}")

    layers = {
        'layer1_rational': H1_diag + V_ee_diag,
        'layer1+2_topological': H1_diag + H1_offdiag + V_ee_diag,
        'layer1+2+3_full': H_full,
        'layer2_only': H1_offdiag,
        'layer3_only': V_ee_offdiag,
    }

    layer_results = {}
    for name, H_layer in layers.items():
        if 'only' in name:
            continue  # These are not full Hamiltonians

        evals, evecs = np.linalg.eigh(H_layer)
        E0 = evals[0]
        ci = evecs[:, 0]

        rho = build_1rdm_from_singlet_ci(ci, configs, n_spatial)
        ent = compute_entanglement_measures(rho)

        print(f"\n  {name}:")
        print(f"    E0 = {E0:.8f} Ha")
        print(f"    S_vN = {ent['von_neumann_entropy']:.6f}")
        print(f"    Top occupations: {ent['occupation_numbers'][:5]}")

        layer_results[name] = {
            'E0': float(E0),
            'von_neumann_entropy': float(ent['von_neumann_entropy']),
            'total_correlation': float(ent['total_correlation']),
            'occupation_numbers': [float(x) for x in ent['occupation_numbers'][:10]],
        }

    return layer_results


# =============================================================================
# Plotting
# =============================================================================

def plot_entropy_vs_nmax(results, save_path):
    """Plot von Neumann entropy vs n_max."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    n_vals = sorted(results.keys())
    S_vals = [results[n]['von_neumann_entropy'] for n in n_vals]
    tc_vals = [results[n]['total_correlation'] for n in n_vals]

    # Left panel: entropy scaling
    axes[0].plot(n_vals, S_vals, 'bo-', linewidth=2, markersize=8, label='S (von Neumann)')
    axes[0].set_xlabel('n_max', fontsize=12)
    axes[0].set_ylabel('Entanglement entropy', fontsize=12)
    axes[0].set_title('Von Neumann entropy vs basis size', fontsize=13)
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    # Right panel: entropy per orbital
    n_orb = [results[n]['n_spatial'] for n in n_vals]
    S_per_orb = [s / no for s, no in zip(S_vals, n_orb)]
    axes[1].plot(n_vals, S_per_orb, 'rs-', linewidth=2, markersize=8)
    axes[1].set_xlabel('n_max', fontsize=12)
    axes[1].set_ylabel('S / n_spatial', fontsize=12)
    axes[1].set_title('Entropy per orbital (area vs volume law)', fontsize=13)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_occupation_spectrum(results, save_path):
    """Plot occupation number spectrum at several n_max values."""
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = plt.cm.viridis(np.linspace(0, 0.8, len(results)))
    for idx, n_max in enumerate(sorted(results.keys())):
        occ = results[n_max]['occupation_numbers']
        # Plot only significant occupations
        sig = [x for x in occ if x > 1e-10]
        ax.semilogy(range(len(sig)), sig, 'o-', color=colors[idx],
                    markersize=5, label=f'n_max={n_max}')

    ax.set_xlabel('Natural orbital index', fontsize=12)
    ax.set_ylabel('Occupation number', fontsize=12)
    ax.set_title('Entanglement spectrum (occupation numbers)', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_per_orbital_entropy(results, save_path):
    """Plot per-orbital entanglement entropy colored by (n,l)."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for idx, n_max in enumerate(sorted(results.keys())):
        if idx >= 6:
            break
        ax = axes[idx]
        data = results[n_max]
        s_i = data['single_orbital_entropy']
        assigns = data.get('natural_orbital_assignments', [])

        # Color by dominant n
        if assigns:
            n_vals = [a['dominant_n'] for a in assigns]
            l_vals = [a['dominant_l'] for a in assigns]
            n_max_val = max(n_vals) if n_vals else 1
            colors = [plt.cm.tab10(n - 1) for n in n_vals]
        else:
            colors = ['blue'] * len(s_i)
            n_vals = [0] * len(s_i)
            l_vals = [0] * len(s_i)

        # Only plot significant entries
        n_show = min(len(s_i), 25)
        bars = ax.bar(range(n_show), s_i[:n_show], color=colors[:n_show])
        ax.set_title(f'n_max = {n_max}', fontsize=11)
        ax.set_xlabel('Natural orbital index')
        ax.set_ylabel('s_i')

        # Add (n,l) labels on top bars
        for k in range(min(n_show, 8)):
            if s_i[k] > 0.001 and assigns:
                ax.text(k, s_i[k] * 1.02, f"({assigns[k]['dominant_n']},{assigns[k]['dominant_l']})",
                        ha='center', va='bottom', fontsize=7, rotation=45)

    plt.suptitle('Per-orbital entanglement entropy by (n,l) character', fontsize=14)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_three_layer(layer_results, save_path):
    """Plot three-layer entanglement decomposition."""
    fig, ax = plt.subplots(figsize=(8, 5))

    layers = ['layer1_rational', 'layer1+2_topological', 'layer1+2+3_full']
    labels = ['Rational\n(diag h1 + diag V_ee)',
              'Rational + Topological\n(+ graph off-diag h1)',
              'Full\n(+ correlation V_ee)']
    S_vals = [layer_results[l]['von_neumann_entropy'] for l in layers]

    bars = ax.bar(range(len(layers)), S_vals, color=['#2ecc71', '#3498db', '#e74c3c'],
                  width=0.6)
    ax.set_xticks(range(len(layers)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel('Von Neumann entropy S', fontsize=12)
    ax.set_title('Three-layer entanglement decomposition', fontsize=13)
    ax.grid(True, alpha=0.3, axis='y')

    for bar, val in zip(bars, S_vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.002,
                f'{val:.4f}', ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


# =============================================================================
# Main
# =============================================================================

if __name__ == '__main__':
    debug_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(debug_dir, 'data')
    plot_dir = os.path.join(debug_dir, 'plots')

    print("="*70)
    print("Entanglement Geometry of He Ground State")
    print("="*70)

    # 1. Graph-native CI analysis
    print("\n\n--- GRAPH-NATIVE CI ANALYSIS ---")
    graph_results = run_graph_native_analysis([2, 3, 4, 5, 6, 7])

    # 2. Variational-k comparison
    print("\n\n--- VARIATIONAL-K CI COMPARISON ---")
    vark_results = run_variational_k_comparison([2, 3, 4, 5])

    # 3. Three-layer decomposition
    print("\n\n--- THREE-LAYER DECOMPOSITION ---")
    layer_results = analyze_three_layer_entanglement(n_max=4)

    # 4. Save results
    all_results = {
        'graph_native': {str(k): v for k, v in graph_results.items()},
        'variational_k': {str(k): v for k, v in vark_results.items()},
        'three_layer': layer_results,
    }

    output_path = os.path.join(data_dir, 'entanglement_results.json')
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved: {output_path}")

    # 5. Generate plots
    print("\n\n--- GENERATING PLOTS ---")
    plot_entropy_vs_nmax(graph_results,
                         os.path.join(plot_dir, 'entanglement_entropy_vs_nmax.png'))
    plot_occupation_spectrum(graph_results,
                             os.path.join(plot_dir, 'entanglement_spectrum.png'))
    plot_per_orbital_entropy(graph_results,
                              os.path.join(plot_dir, 'entanglement_per_orbital.png'))
    plot_three_layer(layer_results,
                      os.path.join(plot_dir, 'entanglement_three_layer.png'))

    # 6. Summary
    print("\n\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print("\nGraph-native CI entanglement entropy:")
    for n_max in sorted(graph_results.keys()):
        r = graph_results[n_max]
        print(f"  n_max={n_max}: S={r['von_neumann_entropy']:.6f}, "
              f"n_configs={r['n_configs']}, error={r['error_pct']:.4f}%")

    print("\nVariational-k CI entanglement entropy:")
    for n_max in sorted(vark_results.keys()):
        r = vark_results[n_max]
        print(f"  n_max={n_max}: S={r['von_neumann_entropy']:.6f}, "
              f"k_opt={r['k_opt']:.4f}, error={r['error_pct']:.4f}%")

    print("\nThree-layer decomposition:")
    for name in ['layer1_rational', 'layer1+2_topological', 'layer1+2+3_full']:
        r = layer_results[name]
        print(f"  {name}: S={r['von_neumann_entropy']:.6f}, E={r['E0']:.6f}")

    # Scaling analysis
    n_vals = sorted(graph_results.keys())
    S_vals = [graph_results[n]['von_neumann_entropy'] for n in n_vals]
    n_orb = [graph_results[n]['n_spatial'] for n in n_vals]
    if len(n_vals) >= 3:
        # Check if S scales as log(M) (area law for 0D bipartition)
        # or linearly with M (volume law)
        from numpy.polynomial import polynomial as P
        log_n = np.log(np.array(n_vals, dtype=float))
        # Fit S = a * log(n_max) + b
        coeffs_log = np.polyfit(log_n, S_vals, 1)
        # Fit S = a * n_max + b
        coeffs_lin = np.polyfit(np.array(n_vals, dtype=float), S_vals, 1)

        print(f"\nScaling analysis:")
        print(f"  S ~ {coeffs_log[0]:.4f} * log(n_max) + {coeffs_log[1]:.4f}")
        print(f"  S ~ {coeffs_lin[0]:.6f} * n_max + {coeffs_lin[1]:.4f}")
        print(f"  Log fit is consistent with area-law-like saturation"
              if abs(coeffs_lin[0]) < 0.01 else
              f"  Linear coefficient nonzero: entropy grows with basis")

    print("\nDone.")
