"""
Entanglement geometry across the first row of the periodic table (He-Ne).

Extends the Sprint 3 analysis (He/Li/Be at n_max=3) to the full first row
(B through Ne) at n_max=2. At n_max=2 there are 5 spatial orbitals (10
spin-orbitals), giving determinant counts from C(10,2)=45 (He) to C(10,10)=1
(Ne), all trivially feasible.

For each atom Z=2..10:
1. Full CI ground state via LatticeIndex with slater_full V_ee
2. 1-RDM, occupation spectrum, von Neumann entropy
3. Single-orbital entropies and mutual information network
4. Hub orbital identification (orbital with largest total MI)
5. Core-valence entanglement partition

Author: GeoVac Development Team
Date: April 2026
"""

import sys
import os
import json
import time
import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from collections import defaultdict
from itertools import combinations

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


# =============================================================================
# 1-RDM from N-electron CI vector (spin-orbital basis)
# =============================================================================

def build_1rdm_spinorb(ci_vector, sd_basis, n_spinorb):
    """Build the spin-orbital 1-RDM from an N-electron CI vector."""
    n_sd = len(sd_basis)
    rho1 = np.zeros((n_spinorb, n_spinorb))

    sd_index = {sd: i for i, sd in enumerate(sd_basis)}

    for I, sd_I in enumerate(sd_basis):
        c_I = ci_vector[I]
        if abs(c_I) < 1e-15:
            continue
        occ_set_I = set(sd_I)

        # Diagonal
        for p in sd_I:
            rho1[p, p] += c_I * c_I

        # Off-diagonal
        for kp, p in enumerate(sd_I):
            remaining = list(sd_I[:kp]) + list(sd_I[kp+1:])
            for q in range(n_spinorb):
                if q == p or q in occ_set_I:
                    continue
                new_sd = tuple(sorted(remaining + [q]))
                J_idx = sd_index.get(new_sd)
                if J_idx is None:
                    continue
                c_J = ci_vector[J_idx]
                if abs(c_J) < 1e-15:
                    continue
                n_before_q_in_J = sum(1 for x in new_sd if x < q)
                n_before_p_in_remaining = sum(1 for x in remaining if x < p)
                phase = (-1) ** (n_before_q_in_J + n_before_p_in_remaining)
                rho1[p, q] += c_I * c_J * phase

    return rho1


def build_spatial_1rdm(rho1_spinorb, n_spatial):
    """Convert spin-orbital 1-RDM to spatial 1-RDM by tracing over spin."""
    rho_spatial = np.zeros((n_spatial, n_spatial))
    for i in range(n_spatial):
        for j in range(n_spatial):
            rho_spatial[i, j] = (rho1_spinorb[2*i, 2*j] +
                                 rho1_spinorb[2*i+1, 2*j+1])
    return rho_spatial


# =============================================================================
# Entanglement measures
# =============================================================================

def compute_entanglement_from_1rdm(rho_spatial, n_electrons):
    """Von Neumann entropy and occupation spectrum from spatial 1-RDM."""
    occ_numbers = np.linalg.eigvalsh(rho_spatial)[::-1]
    rho_norm = rho_spatial / n_electrons
    evals_norm = np.linalg.eigvalsh(rho_norm)
    S = 0.0
    for ev in evals_norm:
        if ev > 1e-15:
            S -= ev * np.log(ev)
    per_orb_entropy = np.zeros(rho_spatial.shape[0])
    for i in range(rho_spatial.shape[0]):
        n_i = rho_spatial[i, i]
        p = n_i / n_electrons
        if p > 1e-15 and p < 1 - 1e-15:
            per_orb_entropy[i] = -p * np.log(p) - (1 - p) * np.log(1 - p)
    return {
        'von_neumann_entropy': float(S),
        'occupation_numbers': occ_numbers,
        'per_orbital_entropy': per_orb_entropy,
    }


def compute_single_orbital_entropies(ci_vector, sd_basis, n_spatial, n_spinorb):
    """Compute single-orbital von Neumann entropy s_i for each spatial orbital."""
    s_single = np.zeros(n_spatial)
    for orb_i in range(n_spatial):
        so_alpha = 2 * orb_i
        so_beta = 2 * orb_i + 1
        p_2 = 0.0
        p_alpha = 0.0
        p_beta = 0.0
        p_0 = 0.0
        for I, sd_I in enumerate(sd_basis):
            c_sq = ci_vector[I] ** 2
            has_alpha = so_alpha in sd_I
            has_beta = so_beta in sd_I
            if has_alpha and has_beta:
                p_2 += c_sq
            elif has_alpha:
                p_alpha += c_sq
            elif has_beta:
                p_beta += c_sq
            else:
                p_0 += c_sq
        probs = [p_0, p_alpha, p_beta, p_2]
        s = 0.0
        for p in probs:
            if p > 1e-15:
                s -= p * np.log(p)
        s_single[orb_i] = s
    return s_single


def compute_two_orbital_entropy(ci_vector, sd_basis, orb_i, orb_j, n_spinorb):
    """Compute two-orbital von Neumann entropy s_{ij}."""
    so_map = {
        2 * orb_i: 0, 2 * orb_i + 1: 1,
        2 * orb_j: 2, 2 * orb_j + 1: 3,
    }
    subspace_sos = set(so_map.keys())
    env_vectors = defaultdict(lambda: np.zeros(16))

    for I, sd_I in enumerate(sd_basis):
        c_I = ci_vector[I]
        if abs(c_I) < 1e-15:
            continue
        sub_occ = [0, 0, 0, 0]
        env_orbs = []
        n_transpositions = 0
        n_sub_before = 0
        for k, so in enumerate(sd_I):
            if so in subspace_sos:
                sub_occ[so_map[so]] = 1
                n_transpositions += k - n_sub_before
                n_sub_before += 1
            else:
                env_orbs.append(so)
        phase = (-1) ** n_transpositions
        env_key = tuple(env_orbs)
        sub_idx = sub_occ[0]*8 + sub_occ[1]*4 + sub_occ[2]*2 + sub_occ[3]
        env_vectors[env_key][sub_idx] += c_I * phase

    rho_ij = np.zeros((16, 16))
    for env_key, vec in env_vectors.items():
        rho_ij += np.outer(vec, vec)
    rho_ij = (rho_ij + rho_ij.T) / 2.0

    evals = np.linalg.eigvalsh(rho_ij)
    evals = np.maximum(evals, 0)
    s_ij = 0.0
    for ev in evals:
        if ev > 1e-15:
            s_ij -= ev * np.log(ev)
    return s_ij


def compute_mutual_information_matrix(ci_vector, sd_basis, n_spatial,
                                      n_spinorb, s_single):
    """Compute mutual information I(i,j) = s_i + s_j - s_{ij} for all pairs."""
    MI = np.zeros((n_spatial, n_spatial))
    significant = np.where(s_single > 1e-8)[0]
    for ii, orb_i in enumerate(significant):
        for jj, orb_j in enumerate(significant):
            if orb_j <= orb_i:
                continue
            s_ij = compute_two_orbital_entropy(ci_vector, sd_basis,
                                               orb_i, orb_j, n_spinorb)
            I_ij = max(0.0, s_single[orb_i] + s_single[orb_j] - s_ij)
            MI[orb_i, orb_j] = I_ij
            MI[orb_j, orb_i] = I_ij
    return MI


def compute_subsystem_entropy(ci_vector, sd_basis, subsystem_orbitals, n_spinorb):
    """Compute von Neumann entropy of a subsystem of spatial orbitals."""
    n_sub_spatial = len(subsystem_orbitals)
    n_sub_spinorb = 2 * n_sub_spatial
    dim_sub = 2 ** n_sub_spinorb

    sub_sos = {}
    idx = 0
    for orb in sorted(subsystem_orbitals):
        sub_sos[2 * orb] = idx
        sub_sos[2 * orb + 1] = idx + 1
        idx += 2
    sub_so_set = set(sub_sos.keys())

    env_vectors = defaultdict(lambda: np.zeros(dim_sub))
    for I, sd_I in enumerate(sd_basis):
        c_I = ci_vector[I]
        if abs(c_I) < 1e-15:
            continue
        sub_occ = [0] * n_sub_spinorb
        env_orbs = []
        n_transpositions = 0
        n_sub_before = 0
        for k, so in enumerate(sd_I):
            if so in sub_so_set:
                sub_occ[sub_sos[so]] = 1
                n_transpositions += k - n_sub_before
                n_sub_before += 1
            else:
                env_orbs.append(so)
        phase = (-1) ** n_transpositions
        env_key = tuple(env_orbs)
        sub_idx = 0
        for bit in sub_occ:
            sub_idx = sub_idx * 2 + bit
        env_vectors[env_key][sub_idx] += c_I * phase

    rho_sub = np.zeros((dim_sub, dim_sub))
    for env_key, vec in env_vectors.items():
        rho_sub += np.outer(vec, vec)
    rho_sub = (rho_sub + rho_sub.T) / 2.0

    evals = np.linalg.eigvalsh(rho_sub)
    evals = np.maximum(evals, 0)
    S = 0.0
    for ev in evals:
        if ev > 1e-15:
            S -= ev * np.log(ev)
    return S


# =============================================================================
# Main per-atom computation
# =============================================================================

def run_atom(Z, n_electrons, n_max, label):
    """Run full entanglement analysis for an atom."""
    print(f"\n{'='*70}")
    print(f"  {label}: Z={Z}, N_e={n_electrons}, n_max={n_max}")
    print(f"{'='*70}")

    t0 = time.time()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        li = __import__('geovac.lattice_index', fromlist=['LatticeIndex']).LatticeIndex(
            n_electrons=n_electrons,
            max_n=n_max,
            nuclear_charge=Z,
            vee_method='slater_full',
        )

    n_spatial = li.lattice.num_states
    n_spinorb = li.n_sp
    n_sd = li.n_sd
    orbitals = list(li.lattice.states)
    orbital_labels = [f"({n},{l},{m})" for n, l, m in orbitals]

    print(f"  n_spatial={n_spatial}, n_spinorb={n_spinorb}, n_sd={n_sd}")

    # Solve for ground state
    print("  Solving FCI...")
    if n_sd == 1:
        # Only one determinant — trivially solved
        H = li.assemble_hamiltonian()
        E0 = float(H[0, 0]) if hasattr(H, 'toarray') else float(H[0, 0])
        ci_vector = np.array([1.0])
        eigvals = np.array([E0])
    else:
        eigvals, eigvecs = li.compute_ground_state(n_states=min(3, n_sd))
        E0 = eigvals[0]
        ci_vector = eigvecs[:, 0]

    # Reference energies (non-relativistic exact, Ha)
    E_exact = {
        2: -2.903724377,   # He
        3: -7.478060323,   # Li
        4: -14.66736,      # Be
        5: -24.65391,      # B
        6: -37.8450,       # C
        7: -54.5892,       # N
        8: -75.0673,       # O
        9: -99.7339,       # F
        10: -128.9376,     # Ne
    }
    ref = E_exact.get(Z)
    if ref is not None:
        error_pct = (E0 - ref) / abs(ref) * 100
        print(f"  E0 = {E0:.6f} Ha (exact: {ref:.6f}, error: {error_pct:.2f}%)")
    else:
        error_pct = None
        print(f"  E0 = {E0:.6f} Ha")

    # Top CI coefficients
    sorted_idx = np.argsort(np.abs(ci_vector))[::-1]
    print("  Top CI coefficients:")
    top_sds = []
    for rank in range(min(5, n_sd)):
        idx = sorted_idx[rank]
        sd = li.sd_basis[idx]
        c = ci_vector[idx]
        labels_list = []
        for so in sd:
            sp = so >> 1
            spin = 'a' if so % 2 == 0 else 'b'
            n, l, m = orbitals[sp]
            labels_list.append(f"({n},{l},{m}){spin}")
        print(f"    [{rank}] c={c:+.6f} (w={c**2:.6f}): {' '.join(labels_list)}")
        top_sds.append({
            'coefficient': float(c),
            'weight': float(c**2),
            'orbitals': labels_list,
        })

    # Build 1-RDM
    print("  Building 1-RDM...")
    rho1_spinorb = build_1rdm_spinorb(ci_vector, li.sd_basis, n_spinorb)
    rho_spatial = build_spatial_1rdm(rho1_spinorb, n_spatial)
    tr = np.trace(rho_spatial)
    print(f"  Tr(rho_spatial) = {tr:.6f} (should be {n_electrons})")

    # Entanglement measures
    ent = compute_entanglement_from_1rdm(rho_spatial, n_electrons)
    print(f"  Von Neumann entropy S = {ent['von_neumann_entropy']:.6f}")

    # Per-orbital occupations
    print("  Per-orbital occupation:")
    for i, (n, l, m) in enumerate(orbitals):
        occ = rho_spatial[i, i]
        if occ > 0.001:
            print(f"    ({n},{l},{m}): n_i = {occ:.6f}")

    # Single-orbital entropies
    print("  Computing single-orbital entropies...")
    s_single = compute_single_orbital_entropies(ci_vector, li.sd_basis,
                                                 n_spatial, n_spinorb)

    # Mutual information
    print("  Computing mutual information matrix...")
    MI = compute_mutual_information_matrix(ci_vector, li.sd_basis,
                                           n_spatial, n_spinorb, s_single)

    # Hub orbital: largest total MI
    total_mi = MI.sum(axis=1)
    hub_idx = np.argmax(total_mi)
    hub_orbital = orbitals[hub_idx]
    hub_label = f"({hub_orbital[0]},{hub_orbital[1]},{hub_orbital[2]})"
    total_mi_sum = total_mi.sum() / 2  # each pair counted twice
    hub_mi_pct = (total_mi[hub_idx] / total_mi_sum * 100) if total_mi_sum > 1e-10 else 0.0

    # Second hub
    total_mi_copy = total_mi.copy()
    total_mi_copy[hub_idx] = -1
    hub2_idx = np.argmax(total_mi_copy)
    hub2_orbital = orbitals[hub2_idx]
    hub2_label = f"({hub2_orbital[0]},{hub2_orbital[1]},{hub2_orbital[2]})"
    hub2_mi_pct = (total_mi[hub2_idx] / total_mi_sum * 100) if total_mi_sum > 1e-10 else 0.0

    print(f"  Hub orbital: {hub_label} ({hub_mi_pct:.1f}% of total MI)")
    print(f"  2nd hub: {hub2_label} ({hub2_mi_pct:.1f}%)")

    # Print top MI pairs
    mi_pairs = []
    for i in range(n_spatial):
        for j in range(i + 1, n_spatial):
            if MI[i, j] > 1e-8:
                mi_pairs.append((i, j, MI[i, j]))
    mi_pairs.sort(key=lambda x: x[2], reverse=True)
    print("  Top mutual information pairs:")
    for i, j, val in mi_pairs[:10]:
        oi = orbitals[i]
        oj = orbitals[j]
        print(f"    I(({oi[0]},{oi[1]},{oi[2]}), ({oj[0]},{oj[1]},{oj[2]})) = {val:.6f}")

    # Core-valence entanglement
    core_orbs = [i for i, (n, l, m) in enumerate(orbitals) if n == 1]
    val_orbs = [i for i, (n, l, m) in enumerate(orbitals) if n >= 2]

    if len(core_orbs) > 0 and len(val_orbs) > 0:
        if len(core_orbs) <= len(val_orbs):
            S_core = compute_subsystem_entropy(ci_vector, li.sd_basis,
                                               core_orbs, n_spinorb)
        else:
            S_core = compute_subsystem_entropy(ci_vector, li.sd_basis,
                                               val_orbs, n_spinorb)
        I_core_val = 2 * S_core  # pure state
    else:
        S_core = 0.0
        I_core_val = 0.0

    print(f"  Core-valence: S_core={S_core:.6f}, I(core,val)={I_core_val:.6f}")

    elapsed = time.time() - t0
    print(f"  Time: {elapsed:.1f}s")

    # Classify hub type
    hub_type = 's' if hub_orbital[1] == 0 else ('p' if hub_orbital[1] == 1 else 'd')

    results = {
        'Z': Z,
        'atom': label,
        'n_electrons': n_electrons,
        'n_max': n_max,
        'n_spatial': n_spatial,
        'n_spinorb': n_spinorb,
        'n_sd': n_sd,
        'E0': float(E0),
        'E_exact': ref,
        'error_pct': float(error_pct) if error_pct is not None else None,
        'von_neumann_entropy': float(ent['von_neumann_entropy']),
        'occupation_numbers': [float(x) for x in ent['occupation_numbers']],
        'per_orbital_occupation': {
            orbital_labels[i]: float(rho_spatial[i, i])
            for i in range(n_spatial)
        },
        'single_orbital_entropy': [float(x) for x in s_single],
        'mutual_information_matrix': MI.tolist(),
        'orbital_labels': orbital_labels,
        'hub_orbital': hub_label,
        'hub_mi_pct': float(hub_mi_pct),
        'hub_type': hub_type,
        'hub2_orbital': hub2_label,
        'hub2_mi_pct': float(hub2_mi_pct),
        'S_core': float(S_core),
        'I_core_valence': float(I_core_val),
        'top_ci_coefficients': top_sds,
        'top_mutual_info_pairs': [
            {'orbital_i': list(orbitals[i]),
             'orbital_j': list(orbitals[j]),
             'I_ij': float(val)}
            for i, j, val in mi_pairs[:15]
        ],
        'elapsed_s': elapsed,
    }

    return results, MI, s_single, orbitals


# =============================================================================
# Plotting
# =============================================================================

def plot_entropy_vs_z(all_results, save_path):
    """Plot von Neumann entropy S vs atomic number Z."""
    atoms = sorted(all_results.values(), key=lambda r: r['Z'])
    Z_vals = [r['Z'] for r in atoms]
    S_vals = [r['von_neumann_entropy'] for r in atoms]
    labels = [r['atom'] for r in atoms]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(Z_vals, S_vals, 'ko-', markersize=8, linewidth=2)

    for z, s, lab in zip(Z_vals, S_vals, labels):
        ax.annotate(lab, (z, s), textcoords="offset points",
                    xytext=(0, 10), ha='center', fontsize=10, fontweight='bold')

    ax.set_xlabel('Atomic number Z', fontsize=12)
    ax.set_ylabel('Von Neumann entropy S (nats)', fontsize=12)
    ax.set_title('Entanglement Entropy Across the First Row (n_max=2)', fontsize=14)
    ax.set_xticks(Z_vals)
    ax.grid(True, alpha=0.3)

    # Mark half-filled p-shell (N, Z=7)
    ax.axvline(x=7, color='blue', linestyle='--', alpha=0.3, label='N (half-filled p)')
    ax.legend()

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_hub_vs_z(all_results, save_path):
    """Plot hub orbital identity and MI% vs Z."""
    atoms = sorted(all_results.values(), key=lambda r: r['Z'])
    Z_vals = [r['Z'] for r in atoms]
    hub_pcts = [r['hub_mi_pct'] for r in atoms]
    hub2_pcts = [r['hub2_mi_pct'] for r in atoms]
    hub_types = [r['hub_type'] for r in atoms]
    hub_labels_list = [r['hub_orbital'] for r in atoms]
    atom_labels = [r['atom'] for r in atoms]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Top panel: hub MI%
    colors_hub = {'s': '#2196F3', 'p': '#FF9800', 'd': '#4CAF50'}
    for z, pct, pct2, ht, hl, al in zip(Z_vals, hub_pcts, hub2_pcts, hub_types,
                                          hub_labels_list, atom_labels):
        ax1.bar(z - 0.15, pct, width=0.3, color=colors_hub.get(ht, 'gray'),
                edgecolor='black', linewidth=0.5)
        ax1.bar(z + 0.15, pct2, width=0.3, color='lightgray',
                edgecolor='black', linewidth=0.5)
        ax1.annotate(hl, (z - 0.15, pct), textcoords="offset points",
                     xytext=(0, 5), ha='center', fontsize=7, rotation=45)

    ax1.set_ylabel('MI share (%)', fontsize=12)
    ax1.set_title('Hub Orbital MI Share (blue=s, orange=p, gray=2nd hub)', fontsize=13)
    ax1.grid(True, alpha=0.3)

    # Bottom panel: hub type as categorical
    type_map = {'s': 0, 'p': 1, 'd': 2}
    type_vals = [type_map.get(ht, -1) for ht in hub_types]
    ax2.plot(Z_vals, type_vals, 'ko-', markersize=10, linewidth=2)
    for z, tv, ht, al in zip(Z_vals, type_vals, hub_types, atom_labels):
        color = colors_hub.get(ht, 'gray')
        ax2.plot(z, tv, 'o', markersize=12, color=color, zorder=3)
        ax2.annotate(al, (z, tv), textcoords="offset points",
                     xytext=(0, -18), ha='center', fontsize=9)

    ax2.set_yticks([0, 1, 2])
    ax2.set_yticklabels(['s', 'p', 'd'])
    ax2.set_xlabel('Atomic number Z', fontsize=12)
    ax2.set_ylabel('Hub orbital type', fontsize=12)
    ax2.set_title('Hub Orbital Type vs Z', fontsize=13)
    ax2.set_xticks(Z_vals)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_occupation_spectra(all_results, save_path):
    """Plot occupation spectra for all atoms overlaid."""
    atoms = sorted(all_results.values(), key=lambda r: r['Z'])

    fig, ax = plt.subplots(figsize=(12, 6))

    cmap = plt.cm.viridis
    n_atoms = len(atoms)

    for k, res in enumerate(atoms):
        occ = np.array(res['occupation_numbers'])
        n_e = res['n_electrons']
        color = cmap(k / max(n_atoms - 1, 1))

        # Normalize occupations to [0, 2] range for comparison
        x = np.arange(len(occ))
        ax.plot(x, occ, 'o-', color=color, markersize=5,
                label=f"{res['atom']} (Z={res['Z']}, S={res['von_neumann_entropy']:.4f})")

    ax.set_xlabel('Natural orbital index', fontsize=12)
    ax.set_ylabel('Occupation number', fontsize=12)
    ax.set_title('Natural Orbital Occupation Spectra (n_max=2)', fontsize=14)
    ax.legend(fontsize=9, loc='upper right')
    ax.set_xlim(-0.3, 5.3)
    ax.set_ylim(-0.05, 2.1)
    ax.axhline(y=2.0, color='gray', linestyle='--', alpha=0.3)
    ax.axhline(y=0.0, color='gray', linestyle='--', alpha=0.3)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_mi_network_panel(all_results, mi_data, s_data, orb_data,
                          panel_atoms, save_path):
    """Plot MI heatmaps for selected atoms in a panel."""
    n_panels = len(panel_atoms)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 4.5))
    if n_panels == 1:
        axes = [axes]

    for ax, atom_label in zip(axes, panel_atoms):
        if atom_label not in mi_data:
            ax.set_title(f'{atom_label} (not computed)')
            continue

        MI = mi_data[atom_label]
        orbitals = orb_data[atom_label]
        n_spatial = MI.shape[0]
        olabels = [f"({n},{l},{m})" for n, l, m in orbitals]

        MI_plot = MI.copy()
        MI_plot[MI_plot < 1e-10] = 1e-10

        mi_nonzero = MI_plot[MI_plot > 1e-10]
        if len(mi_nonzero) > 0:
            vmin = max(1e-6, mi_nonzero.min())
            vmax = mi_nonzero.max()
            if vmax <= vmin:
                vmax = vmin * 10
            im = ax.imshow(MI_plot, cmap='YlOrRd',
                           norm=LogNorm(vmin=vmin, vmax=vmax),
                           aspect='equal')
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        else:
            ax.imshow(np.zeros_like(MI_plot), cmap='YlOrRd', aspect='equal')

        ax.set_xticks(range(n_spatial))
        ax.set_yticks(range(n_spatial))
        ax.set_xticklabels(olabels, rotation=45, ha='right', fontsize=7)
        ax.set_yticklabels(olabels, fontsize=7)

        res = all_results[atom_label]
        ax.set_title(f"{atom_label} (Z={res['Z']}, S={res['von_neumann_entropy']:.4f})")

    fig.suptitle('Mutual Information Networks (n_max=2)', fontsize=14)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


# =============================================================================
# Main
# =============================================================================

def main():
    print("=" * 70)
    print("ENTANGLEMENT GEOMETRY ACROSS THE FIRST ROW (He - Ne)")
    print("n_max=2: 5 spatial orbitals, 10 spin-orbitals")
    print("=" * 70)

    n_max = 2

    atom_specs = [
        (2, 2, "He"),
        (3, 3, "Li"),
        (4, 4, "Be"),
        (5, 5, "B"),
        (6, 6, "C"),
        (7, 7, "N"),
        (8, 8, "O"),
        (9, 9, "F"),
        (10, 10, "Ne"),
    ]

    all_results = {}
    mi_data = {}
    s_data = {}
    orb_data = {}

    for Z, n_e, label in atom_specs:
        res, MI, s_single, orbitals = run_atom(Z, n_e, n_max, label)
        all_results[label] = res
        mi_data[label] = MI
        s_data[label] = s_single
        orb_data[label] = orbitals

    # =================================================================
    # Summary table
    # =================================================================
    print("\n" + "=" * 90)
    print("FIRST-ROW ENTANGLEMENT SUMMARY")
    print("=" * 90)

    print(f"\n{'Z':>3} {'Atom':<4} {'N_e':>4} {'N_SD':>6} {'E0 (Ha)':>12} "
          f"{'Err%':>7} {'S':>8} {'Hub':>10} {'Hub%':>6} "
          f"{'2nd Hub':>10} {'2nd%':>6} {'I_cv':>8}")
    print("-" * 90)

    for Z, n_e, label in atom_specs:
        r = all_results[label]
        err_str = f"{r['error_pct']:.2f}" if r['error_pct'] is not None else "N/A"
        print(f"{r['Z']:>3} {label:<4} {r['n_electrons']:>4} {r['n_sd']:>6} "
              f"{r['E0']:>12.6f} {err_str:>7} {r['von_neumann_entropy']:>8.5f} "
              f"{r['hub_orbital']:>10} {r['hub_mi_pct']:>5.1f}% "
              f"{r['hub2_orbital']:>10} {r['hub2_mi_pct']:>5.1f}% "
              f"{r['I_core_valence']:>8.5f}")

    # =================================================================
    # Key questions
    # =================================================================
    print("\n" + "=" * 70)
    print("KEY FINDINGS")
    print("=" * 70)

    # Hub transition
    print("\nHub orbital progression:")
    for Z, n_e, label in atom_specs:
        r = all_results[label]
        print(f"  {label} (Z={Z}): hub={r['hub_orbital']} ({r['hub_type']}-type, "
              f"{r['hub_mi_pct']:.1f}%), 2nd={r['hub2_orbital']} ({r['hub2_mi_pct']:.1f}%)")

    # Find S maximum
    S_vals = [(label, all_results[label]['von_neumann_entropy']) for _, _, label in atom_specs]
    max_S_atom = max(S_vals, key=lambda x: x[1])
    print(f"\nMaximum entropy: {max_S_atom[0]} with S = {max_S_atom[1]:.6f}")

    # Ne closed shell
    if 'Ne' in all_results:
        ne_s = all_results['Ne']['von_neumann_entropy']
        print(f"Ne (closed shell) entropy: S = {ne_s:.6f}")

    # Save results
    out_path = os.path.join(os.path.dirname(__file__), 'data',
                            'entanglement_first_row.json')
    # Convert numpy arrays for JSON serialization
    json_results = {}
    for label, res in all_results.items():
        json_results[label] = res
    with open(out_path, 'w') as f:
        json.dump(json_results, f, indent=2)
    print(f"\nResults saved to: {out_path}")

    # =================================================================
    # Plots
    # =================================================================
    plot_dir = os.path.join(os.path.dirname(__file__), 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    print("\nGenerating plots...")

    plot_entropy_vs_z(all_results,
                      os.path.join(plot_dir, 'entanglement_first_row_entropy.png'))

    plot_hub_vs_z(all_results,
                  os.path.join(plot_dir, 'entanglement_first_row_hubs.png'))

    plot_occupation_spectra(all_results,
                           os.path.join(plot_dir, 'entanglement_first_row_occupation.png'))

    plot_mi_network_panel(all_results, mi_data, s_data, orb_data,
                          ['B', 'N', 'Ne'],
                          os.path.join(plot_dir, 'entanglement_first_row_network.png'))

    print("\nDone!")


if __name__ == '__main__':
    main()
