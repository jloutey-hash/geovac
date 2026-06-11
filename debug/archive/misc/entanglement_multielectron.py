"""
Entanglement geometry of Li (3e) and Be (4e) ground states.

Extends the He entanglement analysis (Sprints 1-2) to multi-electron atoms.
Uses the LatticeIndex infrastructure with slater_full V_ee for full FCI.

Computes:
1. Ground state via full CI (LatticeIndex with slater_full)
2. 1-RDM from the CI vector: (rho1)_{pq} = <Psi|a+_p a_q|Psi>
3. Entanglement measures: von Neumann entropy, occupation spectrum
4. Mutual information I(i,j) = s_i + s_j - s_{ij} for orbital pairs
5. Core-valence entanglement partition
6. Comparison across He/Li/Be

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
from itertools import combinations

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


# =============================================================================
# 1-RDM from N-electron CI vector (spin-orbital basis)
# =============================================================================

def build_1rdm_spinorb(ci_vector: np.ndarray, sd_basis: list,
                       n_spinorb: int) -> np.ndarray:
    """Build the spin-orbital 1-RDM from an N-electron CI vector.

    (rho1)_{pq} = sum_{I,J} c_I c_J <Phi_I| a+_p a_q |Phi_J>

    The matrix element <Phi_I| a+_p a_q |Phi_J> is nonzero only when:
    - I == J and p == q is occupied in I (diagonal)
    - I and J differ by q -> p replacement (single excitation)

    Parameters
    ----------
    ci_vector : array of shape (n_sd,)
    sd_basis : list of tuples of occupied spin-orbital indices
    n_spinorb : int

    Returns
    -------
    rho1 : array of shape (n_spinorb, n_spinorb)
        Spin-orbital 1-RDM. Tr = N_electrons.
    """
    n_sd = len(sd_basis)
    rho1 = np.zeros((n_spinorb, n_spinorb))

    # Build index for fast lookup
    sd_index = {sd: i for i, sd in enumerate(sd_basis)}

    for I, sd_I in enumerate(sd_basis):
        c_I = ci_vector[I]
        if abs(c_I) < 1e-15:
            continue
        occ_set_I = set(sd_I)

        # Diagonal: <I|a+_p a_p|I> = 1 if p occupied
        for p in sd_I:
            rho1[p, p] += c_I * c_I

        # Off-diagonal: a+_p a_q |J> = +/- |I> when J differs from I by q->p
        # We need <I|a+_p a_q|J> nonzero => J has q occupied (not p),
        # and replacing q with p in J gives I.
        # So for each I and each occupied p in I, consider removing p from I
        # and adding q (not in I) to get J.
        for kp, p in enumerate(sd_I):
            # Remove p from I, add q to get J
            remaining = list(sd_I[:kp]) + list(sd_I[kp+1:])
            remaining_set = set(remaining)

            for q in range(n_spinorb):
                if q == p:
                    continue
                if q in occ_set_I:
                    continue  # q must not be in I (it's in J but not I)

                # J = remaining + {q}, sorted
                new_sd = tuple(sorted(remaining + [q]))
                J_idx = sd_index.get(new_sd)
                if J_idx is None:
                    continue

                c_J = ci_vector[J_idx]
                if abs(c_J) < 1e-15:
                    continue

                # Phase: <I|a+_p a_q|J>
                # a_q |J> removes q from J -> gives remaining (same as I without p)
                # then a+_p inserts p -> gives I
                # Phase from a_q: (-1)^(number of occupied orbitals before q in J)
                n_before_q_in_J = sum(1 for x in new_sd if x < q)
                # Phase from a+_p: (-1)^(number of occupied orbitals before p in result)
                # result after a_q|J> is 'remaining', then inserting p
                n_before_p_in_remaining = sum(1 for x in remaining if x < p)
                phase = (-1) ** (n_before_q_in_J + n_before_p_in_remaining)

                rho1[p, q] += c_I * c_J * phase

    return rho1


def build_spatial_1rdm(rho1_spinorb: np.ndarray, n_spatial: int) -> np.ndarray:
    """Convert spin-orbital 1-RDM to spatial 1-RDM by tracing over spin.

    rho_spatial[i,j] = rho1_spinorb[2*i, 2*j] + rho1_spinorb[2*i+1, 2*j+1]

    Tr(rho_spatial) = N_electrons.
    """
    rho_spatial = np.zeros((n_spatial, n_spatial))
    for i in range(n_spatial):
        for j in range(n_spatial):
            rho_spatial[i, j] = (rho1_spinorb[2*i, 2*j] +
                                 rho1_spinorb[2*i+1, 2*j+1])
    return rho_spatial


def compute_entanglement_from_1rdm(rho_spatial: np.ndarray, n_electrons: int):
    """Compute entanglement measures from spatial 1-RDM.

    Returns dict with von_neumann_entropy, occupation_numbers, per_orbital_entropy.
    """
    # Natural orbital occupation numbers
    occ_numbers = np.linalg.eigvalsh(rho_spatial)[::-1]

    # Von Neumann entropy of the normalized 1-RDM
    # S = -Tr(rho_norm * log(rho_norm)) where rho_norm = rho / N_e
    rho_norm = rho_spatial / n_electrons
    evals_norm = np.linalg.eigvalsh(rho_norm)
    S = 0.0
    for ev in evals_norm:
        if ev > 1e-15:
            S -= ev * np.log(ev)

    # Per-orbital entropy: s_i from diagonal occupation
    per_orb_entropy = np.zeros(rho_spatial.shape[0])
    for i in range(rho_spatial.shape[0]):
        n_i = rho_spatial[i, i]
        # Each spatial orbital can hold 0, 1, or 2 electrons
        # Normalized: n_i / N_e
        p = n_i / n_electrons
        if p > 1e-15 and p < 1 - 1e-15:
            per_orb_entropy[i] = -p * np.log(p) - (1 - p) * np.log(1 - p)

    return {
        'von_neumann_entropy': float(S),
        'occupation_numbers': occ_numbers,
        'per_orbital_entropy': per_orb_entropy,
    }


# =============================================================================
# Single-orbital and two-orbital entropies for mutual information
# =============================================================================

def compute_single_orbital_entropies(ci_vector: np.ndarray, sd_basis: list,
                                     n_spatial: int, n_spinorb: int) -> np.ndarray:
    """Compute single-orbital von Neumann entropy s_i for each spatial orbital.

    For spatial orbital i with spin-orbitals (2i, 2i+1), the reduced density
    matrix rho_i is 4x4 in the {|00>, |10>, |01>, |11>} basis.
    For an eigenstate of S_z, the off-diagonal blocks between different
    particle-number sectors vanish, and we get:
        rho_i = diag(p_0, p_up, p_dn, p_2) for a spin-symmetric state
    where p_0 + p_up + p_dn + p_2 = 1.
    """
    s_single = np.zeros(n_spatial)

    for orb_i in range(n_spatial):
        # Spin-orbitals for this spatial orbital
        so_alpha = 2 * orb_i
        so_beta = 2 * orb_i + 1

        # Compute probabilities: p_0 (empty), p_alpha, p_beta, p_2 (doubly occ)
        p_2 = 0.0  # doubly occupied
        p_alpha = 0.0  # only alpha
        p_beta = 0.0  # only beta
        p_0 = 0.0  # empty

        for I, sd_I in enumerate(sd_basis):
            c_sq = ci_vector[I] ** 2
            has_alpha = so_alpha in sd_I
            has_beta = so_beta in sd_I

            if has_alpha and has_beta:
                p_2 += c_sq
            elif has_alpha and not has_beta:
                p_alpha += c_sq
            elif not has_alpha and has_beta:
                p_beta += c_sq
            else:
                p_0 += c_sq

        # Von Neumann entropy
        probs = [p_0, p_alpha, p_beta, p_2]
        s = 0.0
        for p in probs:
            if p > 1e-15:
                s -= p * np.log(p)
        s_single[orb_i] = s

    return s_single


def compute_two_orbital_entropy(ci_vector: np.ndarray, sd_basis: list,
                                orb_i: int, orb_j: int,
                                n_spinorb: int) -> float:
    """Compute two-orbital von Neumann entropy s_{ij}.

    For spatial orbitals i, j, the reduced density matrix rho_{ij} is 16x16
    in the occupation-number basis of the 4 spin-orbitals
    {i_alpha, i_beta, j_alpha, j_beta}.

    We build rho_{ij} by partial trace of |Psi><Psi| over all other orbitals.
    """
    # Map spin-orbitals to subspace indices
    so_map = {
        2 * orb_i: 0,       # i_alpha
        2 * orb_i + 1: 1,   # i_beta
        2 * orb_j: 2,       # j_alpha
        2 * orb_j + 1: 3,   # j_beta
    }
    subspace_sos = set(so_map.keys())

    # For each SD, decompose into (subspace part, environment part)
    # |SD> = |sub_state> x |env_state>  (with appropriate phase)
    # The subspace state is characterized by which of the 4 spin-orbitals are occupied
    # The environment state is the remaining occupied spin-orbitals

    # Group SDs by their environment state
    # rho_{ij} = sum_env |psi_env><psi_env| where |psi_env> is the subspace vector
    # conditioned on environment state 'env'

    from collections import defaultdict
    env_vectors = defaultdict(lambda: np.zeros(16))

    for I, sd_I in enumerate(sd_basis):
        c_I = ci_vector[I]
        if abs(c_I) < 1e-15:
            continue

        # Split into subspace and environment parts
        sub_occ = [0, 0, 0, 0]  # occupations of the 4 subspace spin-orbitals
        env_orbs = []  # environment spin-orbitals (sorted)

        # Compute phase from reordering: move subspace orbitals to the front
        # Original ordering: sd_I is sorted
        # We want: (subspace orbitals) (environment orbitals)
        # Phase = (-1)^(number of transpositions)
        n_transpositions = 0
        sub_indices_in_sd = []
        for k, so in enumerate(sd_I):
            if so in subspace_sos:
                sub_occ[so_map[so]] = 1
                sub_indices_in_sd.append(k)
            else:
                env_orbs.append(so)

        # Phase: count how many environment orbitals are between subspace orbitals
        # when we move subspace orbs to front
        n_sub_before = 0
        for k, so in enumerate(sd_I):
            if so in subspace_sos:
                # Number of env orbitals before this subspace orbital
                n_env_before_here = k - n_sub_before
                n_transpositions += n_env_before_here
                n_sub_before += 1

        phase = (-1) ** n_transpositions
        env_key = tuple(env_orbs)

        # Subspace state index: binary encoding of sub_occ
        sub_idx = sub_occ[0] * 8 + sub_occ[1] * 4 + sub_occ[2] * 2 + sub_occ[3]
        env_vectors[env_key][sub_idx] += c_I * phase

    # Build rho_{ij} = sum_env |psi_env><psi_env|
    rho_ij = np.zeros((16, 16))
    for env_key, vec in env_vectors.items():
        rho_ij += np.outer(vec, vec)

    # Symmetrize for numerical stability
    rho_ij = (rho_ij + rho_ij.T) / 2.0

    # Von Neumann entropy
    evals = np.linalg.eigvalsh(rho_ij)
    evals = np.maximum(evals, 0)
    s_ij = 0.0
    for ev in evals:
        if ev > 1e-15:
            s_ij -= ev * np.log(ev)

    return s_ij


def compute_mutual_information_matrix(ci_vector: np.ndarray, sd_basis: list,
                                      n_spatial: int, n_spinorb: int,
                                      s_single: np.ndarray) -> np.ndarray:
    """Compute mutual information I(i,j) = s_i + s_j - s_{ij} for all pairs."""
    MI = np.zeros((n_spatial, n_spatial))

    # Only compute for orbitals with non-negligible entropy
    significant = np.where(s_single > 1e-8)[0]
    print(f"  Significant orbitals for MI: {len(significant)} / {n_spatial}")

    for ii, orb_i in enumerate(significant):
        for jj, orb_j in enumerate(significant):
            if orb_j <= orb_i:
                continue

            s_ij = compute_two_orbital_entropy(ci_vector, sd_basis,
                                               orb_i, orb_j, n_spinorb)
            I_ij = s_single[orb_i] + s_single[orb_j] - s_ij
            # MI should be non-negative; small negative values are numerical noise
            I_ij = max(0.0, I_ij)
            MI[orb_i, orb_j] = I_ij
            MI[orb_j, orb_i] = I_ij

    return MI


# =============================================================================
# Core-valence entanglement
# =============================================================================

def compute_core_valence_entropy(ci_vector: np.ndarray, sd_basis: list,
                                 subsystem_orbitals: list, n_spinorb: int) -> float:
    """Compute von Neumann entropy of a subsystem of spatial orbitals.

    Traces out all orbitals NOT in subsystem_orbitals from |Psi><Psi|.

    For subsystems with <= 4 spatial orbitals (8 spin-orbitals, 256 dim),
    uses explicit density matrix. For larger subsystems, uses the 1-RDM
    approximation (exact for weakly correlated subsystems).

    Parameters
    ----------
    subsystem_orbitals : list of spatial orbital indices in the subsystem
    """
    n_sub_spatial = len(subsystem_orbitals)
    n_sub_spinorb = 2 * n_sub_spatial

    # Limit: subsystem must have <= 8 spatial orbitals (16 spin-orbs, 65536 dim)
    if n_sub_spinorb > 16:
        # Use 1-RDM block entropy as approximation
        # Build the 1-RDM restricted to the subsystem
        rho1 = build_1rdm_spinorb(ci_vector, sd_basis, n_spinorb)
        # Extract subsystem block
        sub_sos = []
        for orb in sorted(subsystem_orbitals):
            sub_sos.extend([2 * orb, 2 * orb + 1])
        rho_sub = rho1[np.ix_(sub_sos, sub_sos)]
        # Approximate entropy from 1-RDM eigenvalues
        n_e_sub = np.trace(rho_sub)
        if n_e_sub < 1e-10:
            return 0.0
        rho_norm = rho_sub / n_e_sub
        evals = np.linalg.eigvalsh(rho_norm)
        S = 0.0
        for ev in evals:
            if ev > 1e-15:
                S -= ev * np.log(ev)
        return S

    dim_sub = 2 ** n_sub_spinorb

    # Map subsystem spin-orbitals to subspace indices
    sub_sos = {}
    idx = 0
    for orb in sorted(subsystem_orbitals):
        sub_sos[2 * orb] = idx
        sub_sos[2 * orb + 1] = idx + 1
        idx += 2
    sub_so_set = set(sub_sos.keys())

    # Group by environment state
    from collections import defaultdict
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
                n_env_before_here = k - n_sub_before
                n_transpositions += n_env_before_here
                n_sub_before += 1
            else:
                env_orbs.append(so)

        phase = (-1) ** n_transpositions
        env_key = tuple(env_orbs)

        # Subsystem state index: binary encoding
        sub_idx = 0
        for bit in sub_occ:
            sub_idx = sub_idx * 2 + bit

        env_vectors[env_key][sub_idx] += c_I * phase

    # Build subsystem RDM
    rho_sub = np.zeros((dim_sub, dim_sub))
    for env_key, vec in env_vectors.items():
        rho_sub += np.outer(vec, vec)

    rho_sub = (rho_sub + rho_sub.T) / 2.0

    # Von Neumann entropy
    evals = np.linalg.eigvalsh(rho_sub)
    evals = np.maximum(evals, 0)
    S = 0.0
    for ev in evals:
        if ev > 1e-15:
            S -= ev * np.log(ev)

    return S


# =============================================================================
# Main computation
# =============================================================================

def run_atom(Z: int, n_electrons: int, n_max: int, label: str):
    """Run full entanglement analysis for an atom."""
    print(f"\n{'='*70}")
    print(f"  {label}: Z={Z}, N_e={n_electrons}, n_max={n_max}")
    print(f"{'='*70}")

    t0 = time.time()

    # Build FCI
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
    print("  Assembling Hamiltonian...")
    eigvals, eigvecs = li.compute_ground_state(n_states=3)
    E0 = eigvals[0]
    ci_vector = eigvecs[:, 0]

    # Reference energies (non-relativistic)
    E_exact = {
        (2, 2): -2.903724377,  # He
        (3, 3): -7.478060323,  # Li
        (4, 4): -14.66736,     # Be
    }
    ref = E_exact.get((Z, n_electrons))
    if ref is not None:
        error_pct = (E0 - ref) / abs(ref) * 100
        print(f"  E0 = {E0:.8f} Ha (exact: {ref:.8f}, error: {error_pct:.2f}%)")
    else:
        error_pct = None
        print(f"  E0 = {E0:.8f} Ha")

    # Dominant Slater determinants
    sorted_idx = np.argsort(np.abs(ci_vector))[::-1]
    print("  Top CI coefficients:")
    top_sds = []
    for rank in range(min(5, n_sd)):
        idx = sorted_idx[rank]
        sd = li.sd_basis[idx]
        c = ci_vector[idx]
        # Convert to spatial orbital labels
        labels = []
        for so in sd:
            sp = so >> 1
            spin = 'a' if so % 2 == 0 else 'b'
            n, l, m = orbitals[sp]
            labels.append(f"({n},{l},{m}){spin}")
        print(f"    [{rank}] c={c:+.6f} (w={c**2:.6f}): {' '.join(labels)}")
        top_sds.append({
            'coefficient': float(c),
            'weight': float(c**2),
            'orbitals': labels,
        })

    # Build 1-RDM (spin-orbital)
    print("  Building 1-RDM...")
    rho1_spinorb = build_1rdm_spinorb(ci_vector, li.sd_basis, n_spinorb)
    rho_spatial = build_spatial_1rdm(rho1_spinorb, n_spatial)

    # Verify trace
    tr = np.trace(rho_spatial)
    print(f"  Tr(rho_spatial) = {tr:.6f} (should be {n_electrons})")

    # Entanglement measures
    ent = compute_entanglement_from_1rdm(rho_spatial, n_electrons)
    print(f"  Von Neumann entropy S = {ent['von_neumann_entropy']:.6f}")
    print(f"  Top occupation numbers: {ent['occupation_numbers'][:8]}")

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
    print("  Top single-orbital entropies:")
    s_sorted = np.argsort(s_single)[::-1]
    for k in range(min(10, n_spatial)):
        idx = s_sorted[k]
        n, l, m = orbitals[idx]
        print(f"    ({n},{l},{m}): s={s_single[idx]:.6f}")

    # Mutual information
    print("  Computing mutual information matrix...")
    MI = compute_mutual_information_matrix(ci_vector, li.sd_basis,
                                           n_spatial, n_spinorb, s_single)

    # Top MI pairs
    mi_pairs = []
    for i in range(n_spatial):
        for j in range(i + 1, n_spatial):
            if MI[i, j] > 1e-8:
                mi_pairs.append((i, j, MI[i, j]))
    mi_pairs.sort(key=lambda x: x[2], reverse=True)
    print("  Top mutual information pairs:")
    for i, j, val in mi_pairs[:15]:
        oi = orbitals[i]
        oj = orbitals[j]
        print(f"    I(({oi[0]},{oi[1]},{oi[2]}), ({oj[0]},{oj[1]},{oj[2]})) = {val:.6f}")

    # Core-valence entanglement
    # Core = n=1 orbitals, Valence = n>=2
    core_orbs = [i for i, (n, l, m) in enumerate(orbitals) if n == 1]
    val_orbs = [i for i, (n, l, m) in enumerate(orbitals) if n >= 2]

    # For a pure state, S(core) = S(valence) (Schmidt decomposition).
    # Compute the smaller subsystem's entropy exactly.
    if len(core_orbs) <= len(val_orbs):
        S_core = compute_core_valence_entropy(ci_vector, li.sd_basis,
                                               core_orbs, n_spinorb)
        S_valence = S_core  # pure state => S(A) = S(B)
    else:
        S_valence = compute_core_valence_entropy(ci_vector, li.sd_basis,
                                                  val_orbs, n_spinorb)
        S_core = S_valence

    print(f"  Core-valence entanglement:")
    print(f"    S_core = S_valence = {S_core:.6f} (pure state)")

    # For a pure state, I(core, valence) = S_core + S_valence = 2*S_core
    I_core_val = 2 * S_core
    print(f"    I(core, valence) = 2*S_core = {I_core_val:.6f}")

    elapsed = time.time() - t0
    print(f"  Total time: {elapsed:.1f}s")

    results = {
        'Z': Z,
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
        'top_mutual_info_pairs': [
            {'orbital_i': list(orbitals[i]),
             'orbital_j': list(orbitals[j]),
             'I_ij': float(val)}
            for i, j, val in mi_pairs[:20]
        ],
        'top_ci_coefficients': top_sds,
        'core_orbitals': core_orbs,
        'valence_orbitals': val_orbs,
        'S_core': float(S_core),
        'S_valence': float(S_valence),  # = S_core for pure state
        'I_core_valence': float(I_core_val),  # = 2*S_core for pure state
        'elapsed_s': elapsed,
    }

    return results, MI, s_single, orbitals


# =============================================================================
# Plotting
# =============================================================================

def plot_occupation_spectra(results_dict: dict, save_path: str):
    """Plot occupation number spectra for He/Li/Be side by side."""
    fig, axes = plt.subplots(1, len(results_dict), figsize=(5*len(results_dict), 4.5),
                              sharey=True)
    if len(results_dict) == 1:
        axes = [axes]

    for ax, (label, res) in zip(axes, results_dict.items()):
        occ = np.array(res['occupation_numbers'])
        n_e = res['n_electrons']
        n_orb = len(occ)
        x = np.arange(n_orb)

        # Color by magnitude
        colors = ['#2196F3' if o > 0.1 else '#FF9800' if o > 0.01 else '#9E9E9E'
                  for o in occ]

        ax.bar(x, occ, color=colors, edgecolor='black', linewidth=0.5)
        ax.set_xlabel('Natural orbital index')
        ax.set_title(f'{label} (S={res["von_neumann_entropy"]:.4f})')
        ax.set_xlim(-0.5, min(n_orb, 20) - 0.5)

        # Mark N_e line
        ax.axhline(y=n_e/len(occ), color='red', linestyle='--', alpha=0.3,
                    label=f'uniform = {n_e/len(occ):.3f}')
        ax.legend(fontsize=8)

    axes[0].set_ylabel('Occupation number')
    fig.suptitle('Natural Orbital Occupation Spectra (n_max=3)', fontsize=14)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_mutual_information(MI: np.ndarray, orbitals: list, label: str,
                            save_path: str):
    """Plot mutual information heatmap."""
    n_spatial = MI.shape[0]
    orbital_labels = [f"({n},{l},{m})" for n, l, m in orbitals]

    fig, ax = plt.subplots(figsize=(8, 7))

    # Use log scale for better visibility
    MI_plot = MI.copy()
    MI_plot[MI_plot < 1e-10] = 1e-10

    im = ax.imshow(MI_plot, cmap='YlOrRd',
                   norm=LogNorm(vmin=max(1e-6, MI_plot[MI_plot > 1e-10].min()),
                                vmax=MI_plot.max()),
                   aspect='equal')

    ax.set_xticks(range(n_spatial))
    ax.set_yticks(range(n_spatial))
    ax.set_xticklabels(orbital_labels, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(orbital_labels, fontsize=8)
    ax.set_title(f'Mutual Information I(i,j) - {label}')

    plt.colorbar(im, ax=ax, label='I(i,j)')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_core_valence(results_dict: dict, save_path: str):
    """Plot core-valence entanglement comparison."""
    labels = list(results_dict.keys())
    S_core = [results_dict[l]['S_core'] for l in labels]
    S_val = [results_dict[l]['S_valence'] for l in labels]
    I_cv = [results_dict[l]['I_core_valence'] for l in labels]

    x = np.arange(len(labels))
    width = 0.25

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Bar chart
    ax1.bar(x - width, S_core, width, label='S_core', color='#2196F3')
    ax1.bar(x, S_val, width, label='S_valence', color='#FF9800')
    ax1.bar(x + width, I_cv, width, label='I(core,val)', color='#4CAF50')
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_ylabel('Entropy (nats)')
    ax1.set_title('Core-Valence Entanglement')
    ax1.legend()

    # Ratio: I(core,val) / S_total
    S_total = [results_dict[l]['von_neumann_entropy'] for l in labels]
    ratios = [I_cv[i] / S_total[i] if S_total[i] > 1e-10 else 0 for i in range(len(labels))]

    ax2.bar(labels, ratios, color='#9C27B0')
    ax2.set_ylabel('I(core, valence) / S_total')
    ax2.set_title('Core-Valence Entanglement Fraction')
    ax2.axhline(y=1.0, color='gray', linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_mi_network(MI: np.ndarray, orbitals: list, s_single: np.ndarray,
                    label: str, save_path: str):
    """Plot mutual information as a network graph.

    Node size = single-orbital entropy.
    Edge width/color = mutual information.
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    n_spatial = MI.shape[0]

    # Node positions: arrange by (n, l) shells
    pos = {}
    shell_counts = {}
    for i, (n, l, m) in enumerate(orbitals):
        key = (n, l)
        if key not in shell_counts:
            shell_counts[key] = []
        shell_counts[key].append(i)

    # Arrange shells in concentric rings
    for (n, l), indices in shell_counts.items():
        r = n * 1.5
        n_in_shell = len(indices)
        angle_offset = l * 0.5  # offset by l
        for k, idx in enumerate(indices):
            theta = 2 * np.pi * k / max(n_in_shell, 1) + angle_offset + n * 0.3
            pos[idx] = (r * np.cos(theta), r * np.sin(theta))

    # Draw edges (MI > threshold)
    mi_threshold = 1e-4
    max_mi = MI.max()
    for i in range(n_spatial):
        for j in range(i + 1, n_spatial):
            if MI[i, j] > mi_threshold:
                x = [pos[i][0], pos[j][0]]
                y = [pos[i][1], pos[j][1]]
                # Width and alpha proportional to MI
                w = 0.5 + 4.0 * MI[i, j] / max_mi
                alpha = 0.3 + 0.7 * MI[i, j] / max_mi
                ax.plot(x, y, 'r-', linewidth=w, alpha=alpha, zorder=1)

    # Draw nodes
    max_s = s_single.max() if s_single.max() > 0 else 1.0
    for i, (n, l, m) in enumerate(orbitals):
        x, y = pos[i]
        size = 100 + 800 * s_single[i] / max_s
        color = '#2196F3' if n == 1 else '#FF9800' if n == 2 else '#4CAF50'
        ax.scatter(x, y, s=size, c=color, edgecolors='black', linewidth=1.0,
                   zorder=2)
        ax.annotate(f"({n},{l},{m})", (x, y), textcoords="offset points",
                    xytext=(0, 8), ha='center', fontsize=7, zorder=3)

    ax.set_title(f'Mutual Information Network - {label}\n'
                 f'Blue=n=1, Orange=n=2, Green=n=3; Node size ~ s_i')
    ax.set_aspect('equal')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


# =============================================================================
# Main
# =============================================================================

def main():
    print("=" * 70)
    print("MULTI-ELECTRON ENTANGLEMENT GEOMETRY: He / Li / Be")
    print("=" * 70)

    n_max = 3
    all_results = {}
    mi_data = {}
    s_data = {}
    orb_data = {}

    # --- He (Z=2, 2 electrons) ---
    res, MI, s_single, orbitals = run_atom(2, 2, n_max, "He")
    all_results['He'] = res
    mi_data['He'] = MI
    s_data['He'] = s_single
    orb_data['He'] = orbitals

    # --- Li (Z=3, 3 electrons) ---
    res, MI, s_single, orbitals = run_atom(3, 3, n_max, "Li")
    all_results['Li'] = res
    mi_data['Li'] = MI
    s_data['Li'] = s_single
    orb_data['Li'] = orbitals

    # --- Be (Z=4, 4 electrons) ---
    res, MI, s_single, orbitals = run_atom(4, 4, n_max, "Be")
    all_results['Be'] = res
    mi_data['Be'] = MI
    s_data['Be'] = s_single
    orb_data['Be'] = orbitals

    # =================================================================
    # Summary comparison
    # =================================================================
    print("\n" + "=" * 70)
    print("COMPARISON SUMMARY")
    print("=" * 70)

    print(f"\n{'System':<8} {'Z':>3} {'N_e':>4} {'N_SD':>8} {'E0 (Ha)':>14} "
          f"{'Error%':>8} {'S':>8} {'S_core':>8} {'S_val':>8} {'I_cv':>8}")
    print("-" * 90)
    for label in ['He', 'Li', 'Be']:
        r = all_results[label]
        err_str = f"{r['error_pct']:.2f}" if r['error_pct'] is not None else "N/A"
        print(f"{label:<8} {r['Z']:>3} {r['n_electrons']:>4} {r['n_sd']:>8,} "
              f"{r['E0']:>14.8f} {err_str:>8} {r['von_neumann_entropy']:>8.4f} "
              f"{r['S_core']:>8.4f} {r['S_valence']:>8.4f} {r['I_core_valence']:>8.4f}")

    # Star topology check: for each system, identify if MI is star-centered
    print("\n  Mutual information topology analysis:")
    for label in ['He', 'Li', 'Be']:
        MI = mi_data[label]
        orbitals = orb_data[label]
        n_spatial = MI.shape[0]

        # Hub strength: sum of MI for each orbital
        hub_strength = MI.sum(axis=1)
        hub_idx = np.argmax(hub_strength)
        hub_orb = orbitals[hub_idx]
        total_mi = MI.sum() / 2  # divide by 2 since symmetric

        # Fraction of total MI involving the hub
        hub_mi = hub_strength[hub_idx]
        hub_fraction = hub_mi / (2 * total_mi) if total_mi > 1e-10 else 0

        # Check for secondary hubs
        sorted_hubs = np.argsort(hub_strength)[::-1]
        print(f"\n  {label}:")
        print(f"    Primary hub: ({hub_orb[0]},{hub_orb[1]},{hub_orb[2]}) "
              f"with MI_sum={hub_mi:.6f} ({hub_fraction*100:.1f}% of total)")
        for rank in range(min(5, n_spatial)):
            idx = sorted_hubs[rank]
            orb = orbitals[idx]
            frac = hub_strength[idx] / (2 * total_mi) * 100 if total_mi > 1e-10 else 0
            print(f"    [{rank}] ({orb[0]},{orb[1]},{orb[2]}): "
                  f"MI_sum={hub_strength[idx]:.6f} ({frac:.1f}%)")

    # =================================================================
    # Save results
    # =================================================================
    output_dir = os.path.join(os.path.dirname(__file__), 'data')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'entanglement_multielectron.json')

    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\n  Saved results: {output_path}")

    # =================================================================
    # Generate plots
    # =================================================================
    plot_dir = os.path.join(os.path.dirname(__file__), 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    # Occupation spectra
    plot_occupation_spectra(
        all_results,
        os.path.join(plot_dir, 'entanglement_li_be_spectrum.png')
    )

    # Mutual information heatmaps
    for label in ['Li', 'Be']:
        plot_mutual_information(
            mi_data[label], orb_data[label], label,
            os.path.join(plot_dir, f'entanglement_{label.lower()}_mutual_info.png')
        )

    # MI networks
    for label in ['He', 'Li', 'Be']:
        plot_mi_network(
            mi_data[label], orb_data[label], s_data[label], label,
            os.path.join(plot_dir, f'entanglement_{label.lower()}_network.png')
        )

    # Core-valence comparison
    plot_core_valence(
        all_results,
        os.path.join(plot_dir, 'entanglement_core_valence.png')
    )

    print("\n  All plots saved to debug/plots/")
    print("\nDone.")


if __name__ == '__main__':
    main()
