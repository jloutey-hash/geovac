"""
Entanglement decomposition of the He-like ground state in the graph-native CI basis.

Four investigations:
1. Z-dependence: He-like ions Z=1 to Z=10, with Z_c ~ 1.84 signature
2. Configuration-resolved entanglement: which configs drive entanglement
3. Mutual information matrix: orbital-orbital I(i,j) from 2-RDM
4. Graph-native vs standard CI entanglement comparison

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
from matplotlib.colors import LogNorm

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.casimir_ci import (
    build_graph_native_fci, _build_graph_h1,
    build_fci_matrix, build_fci_polynomial, solve_variational_fast,
    _build_orbital_basis, two_electron_integral,
)

# Import 1-RDM and entanglement code from previous sprint
from debug.entanglement_geometry import (
    build_1rdm_from_singlet_ci,
    compute_entanglement_measures,
    assign_quantum_numbers_to_natural_orbitals,
)


# =============================================================================
# Helper: build graph-native FCI for arbitrary float Z
# =============================================================================

def build_graph_native_fci_float_z(Z_float, n_max, m_total=0):
    """Build graph-native FCI for non-integer Z.

    The graph-native CI uses:
    - Diagonal h1: -Z^2/(2n^2)
    - Off-diagonal h1: kappa * (-A) where kappa=-1/16
    - V_ee: Slater integrals at k_orb = Z

    For non-integer Z we build the h1 matrix manually and construct
    the full Hamiltonian in configuration space.
    """
    from geovac.lattice import GeometricLattice

    lattice = GeometricLattice(max_n=n_max)
    orbitals = list(lattice.states)
    n_spatial = lattice.num_states

    # Build h1 with float Z
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
            if i != j and A_dense[i, j] != 0:
                h1_mat[i, j] = kappa * (-A_dense[i, j])

    # Build configs
    configs = []
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if orbitals[i][2] + orbitals[j][2] == m_total:
                configs.append((i, j))

    n_configs = len(configs)
    H = np.zeros((n_configs, n_configs))
    k_orb = Z_float  # V_ee evaluated at k=Z

    for I in range(n_configs):
        i, j = configs[I]
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
                    # One-body
                    if b == d:
                        me += h1_mat[a, c]
                    if a == c:
                        me += h1_mat[b, d]
                    # Two-body (V_ee)
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

    return H, configs, orbitals, n_spatial


# =============================================================================
# Helper: three-layer decomposition for arbitrary Z
# =============================================================================

def three_layer_decomposition(Z_float, n_max):
    """Build three-layer decomposition for arbitrary Z.

    Returns dict with H matrices and entanglement for each layer.
    """
    from geovac.lattice import GeometricLattice

    lattice = GeometricLattice(max_n=n_max)
    orbitals = list(lattice.states)
    n_spatial = lattice.num_states

    # Build full h1
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
            if i != j and A_dense[i, j] != 0:
                h1_mat[i, j] = kappa * (-A_dense[i, j])

    h1_diag = np.diag(np.diag(h1_mat))
    h1_offdiag = h1_mat - h1_diag

    # Build configs
    configs = []
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if orbitals[i][2] + orbitals[j][2] == 0:
                configs.append((i, j))
    n_configs = len(configs)
    k_orb = Z_float

    # Build config-space matrices
    def build_h1_config(h1_src):
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
                            me += h1_src[a, c]
                        if a == c:
                            me += h1_src[b, d]
                me /= (N_I * N_J)
                H1[I, J] = me
                H1[J, I] = me
        return H1

    def build_vee_config():
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
                        na, la, ma_ = orbitals[a]
                        nb, lb, mb_ = orbitals[b]
                        nc, lc, mc_ = orbitals[c]
                        nd, ld, md_ = orbitals[d]
                        me += two_electron_integral(
                            na, la, ma_, nb, lb, mb_,
                            nc, lc, mc_, nd, ld, md_, k_orb)
                me /= (N_I * N_J)
                V[I, J] = me
                V[J, I] = me
        return V

    H1_diag = build_h1_config(h1_diag)
    H1_offdiag = build_h1_config(h1_offdiag)
    V_ee = build_vee_config()
    V_ee_diag = np.diag(np.diag(V_ee))
    V_ee_offdiag = V_ee - V_ee_diag

    layers = {
        'layer1_rational': H1_diag + V_ee_diag,
        'layer1+2_topological': H1_diag + H1_offdiag + V_ee_diag,
        'layer1+2+3_full': H1_diag + H1_offdiag + V_ee,
    }

    results = {}
    for name, H_layer in layers.items():
        evals, evecs = np.linalg.eigh(H_layer)
        E0 = evals[0]
        ci = evecs[:, 0]
        rho = build_1rdm_from_singlet_ci(ci, configs, n_spatial)
        ent = compute_entanglement_measures(rho)
        results[name] = {
            'E0': float(E0),
            'von_neumann_entropy': float(ent['von_neumann_entropy']),
            'occupation_numbers': [float(x) for x in ent['occupation_numbers'][:10]],
        }

    return results


# =============================================================================
# Part 1: Z-dependence of entanglement
# =============================================================================

def part1_z_dependence(n_max=4):
    """Compute entanglement entropy as a function of Z for He-like ions."""
    print("\n" + "=" * 70)
    print("PART 1: Z-DEPENDENCE OF ENTANGLEMENT (n_max={})".format(n_max))
    print("=" * 70)

    # Exact energies for He-like ions (non-relativistic, from NIST/literature)
    E_exact = {
        1: -0.52775101,   # H-
        2: -2.903724377,  # He
        3: -7.279913413,  # Li+
        4: -13.65556624,  # Be2+
        5: -22.03097159,  # B3+
        6: -32.40624671,  # C4+
        8: -59.15696,     # O6+
        10: -93.9068,     # Ne8+
    }

    Z_values = [1, 1.5, 1.84, 2, 3, 4, 5, 6, 8, 10]
    results = {}

    for Z in Z_values:
        t0 = time.time()
        print(f"\n  Z = {Z}...")

        H, configs, orbitals, n_spatial = build_graph_native_fci_float_z(Z, n_max)
        n_configs = H.shape[0]

        evals, evecs = np.linalg.eigh(H)
        E0 = evals[0]
        ci_coeffs = evecs[:, 0]

        # Error relative to exact (if available)
        Z_int = int(Z) if Z == int(Z) else None
        if Z_int and Z_int in E_exact:
            error_pct = (E0 - E_exact[Z_int]) / abs(E_exact[Z_int]) * 100
        else:
            error_pct = None

        # Build 1-RDM and entanglement
        rho = build_1rdm_from_singlet_ci(ci_coeffs, configs, n_spatial)
        ent = compute_entanglement_measures(rho)

        # Three-layer decomposition
        layer = three_layer_decomposition(Z, n_max)

        S_full = layer['layer1+2+3_full']['von_neumann_entropy']
        S_topo = layer['layer1+2_topological']['von_neumann_entropy']
        S_rat = layer['layer1_rational']['von_neumann_entropy']

        # Fractions
        if S_full > 1e-15:
            frac_vee = (S_full - S_topo) / S_full * 100
            frac_h1_offdiag = (S_topo - S_rat) / S_full * 100
            frac_rational = S_rat / S_full * 100
        else:
            frac_vee = frac_h1_offdiag = frac_rational = 0.0

        elapsed = time.time() - t0

        print(f"    E0 = {E0:.8f} Ha", end="")
        if error_pct is not None:
            print(f"  (error = {error_pct:.2f}%)")
        else:
            print()
        print(f"    S = {S_full:.6f}")
        print(f"    Fractions: V_ee={frac_vee:.1f}%, h1_offdiag={frac_h1_offdiag:.1f}%, "
              f"rational={frac_rational:.1f}%")
        print(f"    Top occ: {ent['occupation_numbers'][:5]}")
        print(f"    Time: {elapsed:.1f}s")

        # Dominant CI coefficients
        sorted_ci = np.argsort(np.abs(ci_coeffs))[::-1]
        top_ci = []
        for rank, idx in enumerate(sorted_ci[:5]):
            i, j = configs[idx]
            oi, oj = orbitals[i], orbitals[j]
            top_ci.append({
                'orbitals': [list(oi), list(oj)],
                'coefficient': float(ci_coeffs[idx]),
                'weight': float(ci_coeffs[idx] ** 2),
            })

        results[str(Z)] = {
            'Z': float(Z),
            'E0': float(E0),
            'error_pct': float(error_pct) if error_pct is not None else None,
            'n_configs': n_configs,
            'von_neumann_entropy': float(S_full),
            'S_rational': float(S_rat),
            'S_topological': float(S_topo),
            'S_full': float(S_full),
            'frac_vee_pct': float(frac_vee),
            'frac_h1_offdiag_pct': float(frac_h1_offdiag),
            'frac_rational_pct': float(frac_rational),
            'occupation_numbers': [float(x) for x in ent['occupation_numbers'][:15]],
            'top_ci_coefficients': top_ci,
            'time_s': float(elapsed),
        }

    return results


# =============================================================================
# Part 2: Configuration-resolved entanglement
# =============================================================================

def part2_config_resolved(n_max=4, Z=2):
    """Decompose entanglement by configuration contribution."""
    print("\n" + "=" * 70)
    print(f"PART 2: CONFIGURATION-RESOLVED ENTANGLEMENT (Z={Z}, n_max={n_max})")
    print("=" * 70)

    H, configs, orbitals, n_spatial = build_graph_native_fci_float_z(Z, n_max)
    evals, evecs = np.linalg.eigh(H)
    ci_coeffs = evecs[:, 0]

    n_configs = len(configs)

    # Build 1-RDM
    rho = build_1rdm_from_singlet_ci(ci_coeffs, configs, n_spatial)
    ent = compute_entanglement_measures(rho)
    S_total = ent['von_neumann_entropy']

    # CI coefficient analysis
    print(f"\n  Total configs: {n_configs}")
    print(f"  Total entropy: {S_total:.6f}")

    # Group configs by (n1,l1)x(n2,l2) type
    config_types = {}
    for idx, (i, j) in enumerate(configs):
        ni, li, mi = orbitals[i]
        nj, lj, mj = orbitals[j]
        key = f"({ni},{li})x({nj},{lj})"
        if key not in config_types:
            config_types[key] = {
                'indices': [],
                'total_weight': 0.0,
                'n1': ni, 'l1': li, 'n2': nj, 'l2': lj,
            }
        config_types[key]['indices'].append(idx)
        config_types[key]['total_weight'] += ci_coeffs[idx] ** 2

    # Sort by total weight
    sorted_types = sorted(config_types.items(),
                          key=lambda x: x[1]['total_weight'], reverse=True)

    print(f"\n  Configuration types by weight:")
    for key, data in sorted_types[:15]:
        print(f"    {key}: weight={data['total_weight']:.6f}, "
              f"count={len(data['indices'])}")

    # Compute entanglement contribution of each config type
    # The 1-RDM off-diagonal elements come from cross-terms c_I * c_J
    # where configs I and J share one orbital.
    # We compute the "partial 1-RDM" from each config pair type.

    # For each configuration type, compute its contribution to the 1-RDM
    # off-diagonal elements (which drive entanglement)
    type_contributions = {}
    for key, data in sorted_types:
        # Build partial 1-RDM from this type's configurations cross all others
        rho_partial = np.zeros((n_spatial, n_spatial))
        for idx in data['indices']:
            i, j = configs[idx]
            for J, (p, q) in enumerate(configs):
                coeff = ci_coeffs[idx] * ci_coeffs[J]
                if abs(coeff) < 1e-16:
                    continue
                N_I = np.sqrt(2.0) if i != j else 1.0
                N_J = np.sqrt(2.0) if p != q else 1.0
                I_perms = [(i, j)]
                if i != j:
                    I_perms.append((j, i))
                J_perms = [(p, q)]
                if p != q:
                    J_perms.append((q, p))
                for a, b in I_perms:
                    for c, d in J_perms:
                        if b == d:
                            rho_partial[a, c] += coeff / (N_I * N_J)
        rho_partial *= 2.0

        # Contribution measure: Frobenius norm of off-diagonal part
        rho_offdiag = rho_partial - np.diag(np.diag(rho_partial))
        frob_contribution = np.linalg.norm(rho_offdiag, 'fro')

        type_contributions[key] = {
            'weight': float(data['total_weight']),
            'count': len(data['indices']),
            'rho_offdiag_frobenius': float(frob_contribution),
            'n1': data['n1'], 'l1': data['l1'],
            'n2': data['n2'], 'l2': data['l2'],
        }

    # Dominant cross-configuration pairs
    # Which pairs (I, J) have the largest |c_I * c_J| product
    # where they share an orbital?
    print(f"\n  Configuration pair entanglement channels:")
    cross_pairs = []
    sorted_ci = np.argsort(np.abs(ci_coeffs))[::-1]
    # Take top-20 configs and analyze their cross-terms
    top_indices = sorted_ci[:20]
    for ii, I_idx in enumerate(top_indices):
        for jj, J_idx in enumerate(top_indices):
            if J_idx <= I_idx:
                continue
            i1, j1 = configs[I_idx]
            i2, j2 = configs[J_idx]
            # Check if they share an orbital
            shared = set()
            if i1 == i2 or i1 == j2:
                shared.add(i1)
            if j1 == i2 or j1 == j2:
                shared.add(j1)
            product = abs(ci_coeffs[I_idx] * ci_coeffs[J_idx])
            if product > 1e-6:
                o1 = f"({orbitals[i1][0]},{orbitals[i1][1]},{orbitals[i1][2]})"
                o2 = f"({orbitals[j1][0]},{orbitals[j1][1]},{orbitals[j1][2]})"
                o3 = f"({orbitals[i2][0]},{orbitals[i2][1]},{orbitals[i2][2]})"
                o4 = f"({orbitals[j2][0]},{orbitals[j2][1]},{orbitals[j2][2]})"
                cross_pairs.append({
                    'config_I': f"{o1}x{o2}",
                    'config_J': f"{o3}x{o4}",
                    'c_I': float(ci_coeffs[I_idx]),
                    'c_J': float(ci_coeffs[J_idx]),
                    'product': float(product),
                    'shared_orbitals': len(shared),
                })

    cross_pairs.sort(key=lambda x: x['product'], reverse=True)
    for cp in cross_pairs[:10]:
        print(f"    {cp['config_I']} <-> {cp['config_J']}: "
              f"|c_I*c_J|={cp['product']:.6f}, shared={cp['shared_orbitals']}")

    results = {
        'Z': Z,
        'n_max': n_max,
        'n_configs': n_configs,
        'total_entropy': float(S_total),
        'config_types': type_contributions,
        'cross_pairs': cross_pairs[:20],
        'ci_coefficients': {
            str(idx): {
                'orbitals': [list(orbitals[configs[idx][0]]),
                             list(orbitals[configs[idx][1]])],
                'coefficient': float(ci_coeffs[idx]),
                'weight': float(ci_coeffs[idx] ** 2),
            }
            for idx in sorted_ci[:20]
        },
    }

    return results


# =============================================================================
# Part 3: Mutual information matrix
# =============================================================================

def part3_mutual_information(n_max=4, Z=2):
    """Compute orbital-orbital mutual information I(i,j) from 2-RDM."""
    print("\n" + "=" * 70)
    print(f"PART 3: MUTUAL INFORMATION MATRIX (Z={Z}, n_max={n_max})")
    print("=" * 70)

    H, configs, orbitals, n_spatial = build_graph_native_fci_float_z(Z, n_max)
    evals, evecs = np.linalg.eigh(H)
    ci_coeffs = evecs[:, 0]

    # Build 1-RDM for single-orbital entropies
    rho = build_1rdm_from_singlet_ci(ci_coeffs, configs, n_spatial)
    ent = compute_entanglement_measures(rho)

    # Single-orbital entropies from occupation numbers of 1-RDM
    # For spatial orbital i, occupation n_i comes from rho eigenvalues.
    # But we need s_i for each ORIGINAL orbital, not natural orbital.
    # s_i = -n_i/2 * log(n_i/2) - (1 - n_i/2) * log(1 - n_i/2)
    # where n_i = rho[i,i] (diagonal of 1-RDM in original basis)
    # multiplied by 2 for spin (two spin-orbitals per spatial orbital).

    # For mutual information we work in the spin-orbital picture.
    # Each spatial orbital gives 2 spin-orbitals (alpha, beta).
    # For a singlet, the occupation of alpha and beta spin-orbitals is equal.
    # The spin-orbital occupation is n_i^(spin) = rho[i,i] / 2 (since Tr(rho) = 2
    # and each spatial orbital splits into 2 spin-orbitals with equal occupation).
    # Actually for the single-orbital reduced density matrix, we should
    # trace over everything except orbital i.

    # For a 2-electron singlet, the reduced density matrix for spatial orbital i
    # in the {|0>, |alpha>, |beta>, |alpha beta>} occupation number basis is:
    # rho_i = diag(1 - n_i/2, ...) -- but more carefully:
    #
    # Let n_i = rho_{ii} = sum_I,J c_I c_J * <I| n_hat_i |J> (spatial)
    # The spin-orbital RDM for orbital i is 4x4 in {|00>, |10>, |01>, |11>}:
    # For a singlet, <a+_alpha a_alpha> = <a+_beta a_beta> = n_i / 2
    # and <a+_alpha a+_beta a_beta a_alpha> = D_ii (double occupancy)
    # where D_ii = Prob(both electrons in orbital i)

    # Compute D_ii (double-occupancy probability for each orbital)
    D = np.zeros(n_spatial)
    for I, (i, j) in enumerate(configs):
        if i == j:  # doubly occupied
            D[i] += ci_coeffs[I] ** 2

    # Single-orbital entropy in the {|00>, |up>, |dn>, |ud>} basis
    # For a singlet: rho_i = diag(p_0, p_up, p_dn, p_2)
    # where p_2 = D[i], p_up = p_dn = (n_i - 2*D[i])/2, p_0 = 1 - n_i + D[i]
    # This is because: n_i = 2*D[i] + p_up + p_dn (each singly-occupied contributes 1)
    # and by singlet symmetry p_up = p_dn.

    def single_orbital_entropy(i):
        """Compute s_i for spatial orbital i."""
        n_i = rho[i, i]
        d_i = D[i]
        p_2 = d_i
        p_1 = (n_i - 2 * d_i) / 2.0  # each spin
        p_0 = 1.0 - n_i + d_i

        # Clamp to avoid log(0)
        probs = [p_0, p_1, p_1, p_2]
        s = 0.0
        for p in probs:
            if p > 1e-15:
                s -= p * np.log(p)
        return s

    s_single = np.array([single_orbital_entropy(i) for i in range(n_spatial)])

    print(f"\n  n_spatial = {n_spatial}")
    print(f"  Top single-orbital entropies:")
    s_sorted = np.argsort(s_single)[::-1]
    for k in range(min(10, n_spatial)):
        idx = s_sorted[k]
        n, l, m = orbitals[idx]
        print(f"    orbital ({n},{l},{m}): s={s_single[idx]:.6f}, "
              f"n_i={rho[idx,idx]:.6f}, D_i={D[idx]:.6f}")

    # Compute two-orbital reduced density matrix for each pair (i, j)
    # For a 2-electron system, the 2-RDM IS the full wavefunction.
    # Gamma_{pq,rs} = <Psi| a+_p a+_q a_s a_r |Psi>
    #
    # For a singlet spin-adapted CI with spatial configs:
    # Each config |I> = (phi_a, phi_b) in the singlet spatial combination.
    # The full wavefunction (including spin) is:
    # |Psi> = sum_I c_I * [a+_{a,alpha} a+_{b,beta} - a+_{a,beta} a+_{b,alpha}] / N_I / sqrt(2) |0>
    #
    # For two spatial orbitals i, j, the reduced density matrix in the
    # 16-dimensional spin-orbital basis {|n_i_up, n_i_dn, n_j_up, n_j_dn>}
    # is computed from the 2-RDM elements.
    #
    # For a 2-electron system this simplifies: the total wavefunction has
    # exactly 2 electrons, so the 2-orbital RDM in the occupation-number
    # basis is a 16x16 matrix, but only the 2-electron sector is nonzero
    # (the total state is pure and has exactly 2 electrons).
    #
    # Actually, for the two-orbital entropy s_{ij}, we trace out everything
    # EXCEPT orbitals i and j. The resulting rho_{ij} is a 16x16 matrix
    # (4 states per orbital x 4 states).
    #
    # For 2 electrons total, rho_{ij} lives in the 0, 1, or 2 electron sectors
    # of the (i,j) subspace.

    # Build the full spin-orbital 2-RDM from the singlet CI
    # For a singlet: <a+_{p,s1} a+_{q,s2} a_{r,s2} a_{s,s1}> involves
    # spatial and spin indices.
    # For the two-orbital entropy we need:
    # rho_{ij} = Tr_{k != i,j} |Psi><Psi|

    # Simpler approach for 2-electron systems:
    # The wavefunction in spin-orbital basis is:
    # |Psi> = sum_I c_I |phi_I>
    # where |phi_I> for singlet config (a, b) is:
    #   (a_alpha, b_beta) - (a_beta, b_alpha) for a != b, times 1/sqrt(2*N_I^2)
    #   (a_alpha, a_beta) for a == b, times 1

    # For each pair of spatial orbitals (i, j), the 2-orbital RDM in the
    # occupation-number basis is obtained by projecting the wavefunction
    # onto states where the 2 electrons are distributed among the 4
    # spin-orbitals {i_alpha, i_beta, j_alpha, j_beta}.

    # Enumerate the 2-electron basis states for spin-orbitals of (i, j):
    # With 4 spin-orbitals and 2 electrons, there are C(4,2) = 6 states:
    # |1100>, |1010>, |1001>, |0110>, |0101>, |0011>
    # These are ordered as: i_a i_b, i_a j_a, i_a j_b, i_b j_a, i_b j_b, j_a j_b

    # For the mutual information, we need I(i,j) = s_i + s_j - s_{ij}
    # where s_{ij} is the von Neumann entropy of the two-orbital RDM.

    # For a 2-electron system, each config puts 2 electrons in specific
    # spin-orbitals. The projection onto the (i,j) subspace:
    # - If both electrons are in orbitals i and j: full overlap
    # - If one electron is in i or j and one elsewhere: partial
    # - If neither electron is in i or j: zero sector

    # Build the two-orbital density matrix
    def two_orbital_entropy(orb_i, orb_j):
        """Compute s_{ij} for spatial orbitals orb_i, orb_j.

        We work in the full 16-dim occupation-number basis of the 4
        spin-orbitals {i_alpha, i_beta, j_alpha, j_beta}, but restricted
        to the 2-electron sector (6 states).

        For a singlet state, the amplitude that both electrons land in
        the (i, j) subspace depends on the CI coefficients.
        """
        # The 2-electron states in the (i,j) spin-orbital space:
        # Label spin-orbitals as 0=i_alpha, 1=i_beta, 2=j_alpha, 3=j_beta
        # 2-electron basis: all pairs (p, q) with p < q
        two_e_basis = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        n_basis = len(two_e_basis)

        # Wavefunction amplitudes in this subspace
        # For each CI config (a, b) in singlet:
        # The spin wavefunction is [alpha(1)beta(2) - beta(1)alpha(2)] / sqrt(2)
        # Spatial: [phi_a(1)phi_b(2) + phi_b(1)phi_a(2)] / N_ab
        # Combined: sum over spatial permutations x spin antisymmetrization

        # Map spatial orbital to spin-orbital indices in the (i,j) subspace
        def spatial_to_subspace(s):
            """Return list of (spin_orb_index, spin) for spatial orbital s in (i,j) subspace."""
            if s == orb_i:
                return [(0, 'alpha'), (1, 'beta')]
            elif s == orb_j:
                return [(2, 'alpha'), (3, 'beta')]
            else:
                return []

        # For each config, compute the amplitude in each 2-electron basis state
        # We need to handle the full sector decomposition:
        # - 2 electrons in (i,j): amplitude in the 6-dim subspace
        # - 1 electron in (i,j), 1 outside: contributes to 1-electron sector
        # - 0 electrons in (i,j): contributes to 0-electron sector

        # For the full 16-dim density matrix, we work sector by sector.
        # rho_{ij} = p_0 |vac><vac| + sum_k p_k |k><k| + rho_2
        # where |vac> is the 0-electron state of (i,j), |k> are 1-electron states,
        # and rho_2 is the 2-electron sector.

        # Full 16-dim basis: enumerate all occupation patterns of 4 spin-orbitals
        # State index = binary encoding of (n_{i_alpha}, n_{i_beta}, n_{j_alpha}, n_{j_beta})

        rho_ij = np.zeros((16, 16))

        for I, (a, b) in enumerate(configs):
            for J, (c, d) in enumerate(configs):
                coeff = ci_coeffs[I] * ci_coeffs[J]
                if abs(coeff) < 1e-16:
                    continue

                # The singlet spatial function for config (a, b):
                # |Psi_I> = N_ab^{-1} * sum_{perm} [phi_a(1) phi_b(2) + perm]
                # with singlet spin: [alpha beta - beta alpha] / sqrt(2)

                N_I = np.sqrt(2.0) if a != b else 1.0
                N_J = np.sqrt(2.0) if c != d else 1.0

                # The spin-orbital determinants for config (a, b) singlet:
                # if a != b: |a_alpha b_beta> - |a_beta b_alpha> + |b_alpha a_beta> - |b_beta a_alpha>
                #          = |a_alpha b_beta> - |a_beta b_alpha> + |b_alpha a_beta> - |b_beta a_alpha>
                # simplified: [a_alpha b_beta - a_beta b_alpha] * sqrt(2) / N_ab / sqrt(2)
                # + [b_alpha a_beta - b_beta a_alpha] * sqrt(2) / N_ab / sqrt(2)
                # For a != b, the full state (spatial symmetric x spin antisymmetric):
                # = [phi_a(1)phi_b(2) + phi_b(1)phi_a(2)] / sqrt(2)
                #   x [alpha(1)beta(2) - beta(1)alpha(2)] / sqrt(2)
                # In second quantization:
                # = (a+_{a,alpha} a+_{b,beta} - a+_{a,beta} a+_{b,alpha}) / sqrt(2)  if a == b: just a+_{a,alpha} a+_{a,beta}

                # List of spin-orbital pairs (with signs) for config I
                def get_spinorb_decomp(x, y):
                    """Return list of (spinorb1, spinorb2, amplitude) for singlet config (x, y)."""
                    # spinorb = (spatial_idx, spin): spin=0 for alpha, spin=1 for beta
                    if x == y:
                        # Doubly occupied: a+_{x,alpha} a+_{x,beta} |0>
                        return [((x, 0), (x, 1), 1.0)]
                    else:
                        # Singlet: (a+_{x,alpha} a+_{y,beta} - a+_{x,beta} a+_{y,alpha}) / sqrt(2)
                        # Plus spatial permutation: already included in N_ab normalization
                        # The full normalized state is:
                        # 1/sqrt(2) * (a+_{x,alpha} a+_{y,beta} - a+_{x,beta} a+_{y,alpha}) |0>
                        # (the spatial symmetric part with normalization is handled by N_ab)
                        return [
                            ((x, 0), (y, 1), 1.0 / np.sqrt(2.0)),
                            ((x, 1), (y, 0), -1.0 / np.sqrt(2.0)),
                        ]

                decomp_I = get_spinorb_decomp(a, b)
                decomp_J = get_spinorb_decomp(c, d)

                # For each pair of spin-orbital determinants, trace out
                # everything except orbitals orb_i and orb_j
                for (s1, s2, amp_I) in decomp_I:
                    for (s3, s4, amp_J) in decomp_J:
                        # s1, s2 are (spatial, spin) for electrons 1, 2 in |I>
                        # s3, s4 are (spatial, spin) for electrons 1, 2 in |J>
                        # a+_{s1} a+_{s2} |0>  and  a+_{s3} a+_{s4} |0>

                        total_amp = coeff * amp_I * amp_J / (N_I * N_J)
                        if abs(total_amp) < 1e-16:
                            continue

                        # Classify each electron as "in subspace" or "outside"
                        def in_subspace(so):
                            return so[0] == orb_i or so[0] == orb_j

                        def to_subspace_idx(so):
                            """Map (spatial, spin) to 0-3 index in (i,j) subspace."""
                            if so[0] == orb_i:
                                return so[1]  # 0=alpha, 1=beta
                            elif so[0] == orb_j:
                                return 2 + so[1]  # 2=alpha, 3=beta
                            return -1

                        in1 = in_subspace(s1)
                        in2 = in_subspace(s2)
                        in3 = in_subspace(s3)
                        in4 = in_subspace(s4)

                        # For the partial trace, the "outside" electrons must match
                        # between bra and ket for the contribution to be nonzero.

                        if in1 and in2 and in3 and in4:
                            # Both electrons in subspace
                            idx_bra1 = to_subspace_idx(s1)
                            idx_bra2 = to_subspace_idx(s2)
                            idx_ket1 = to_subspace_idx(s3)
                            idx_ket2 = to_subspace_idx(s4)

                            # Occupation number basis index: binary encoding
                            bra_occ = [0, 0, 0, 0]
                            bra_occ[idx_bra1] = 1
                            bra_occ[idx_bra2] = 1
                            ket_occ = [0, 0, 0, 0]
                            ket_occ[idx_ket1] = 1
                            ket_occ[idx_ket2] = 1

                            bra_idx = bra_occ[0] * 8 + bra_occ[1] * 4 + bra_occ[2] * 2 + bra_occ[3]
                            ket_idx = ket_occ[0] * 8 + ket_occ[1] * 4 + ket_occ[2] * 2 + ket_occ[3]

                            # Fermionic sign from ordering
                            # a+_{s1} a+_{s2} |0> -> need canonical ordering
                            sign_bra = 1 if idx_bra1 < idx_bra2 else -1
                            sign_ket = 1 if idx_ket1 < idx_ket2 else -1

                            rho_ij[bra_idx, ket_idx] += total_amp * sign_bra * sign_ket

                        elif in1 and not in2 and in3 and not in4:
                            # Electron 1 in subspace, electron 2 outside
                            # Partial trace: outside electrons must match
                            if s2 == s4:  # same outside spin-orbital
                                idx_bra = to_subspace_idx(s1)
                                idx_ket = to_subspace_idx(s3)
                                bra_occ = [0, 0, 0, 0]
                                bra_occ[idx_bra] = 1
                                ket_occ = [0, 0, 0, 0]
                                ket_occ[idx_ket] = 1
                                bra_idx = bra_occ[0] * 8 + bra_occ[1] * 4 + bra_occ[2] * 2 + bra_occ[3]
                                ket_idx = ket_occ[0] * 8 + ket_occ[1] * 4 + ket_occ[2] * 2 + ket_occ[3]
                                rho_ij[bra_idx, ket_idx] += total_amp

                        elif not in1 and in2 and not in3 and in4:
                            # Electron 2 in subspace, electron 1 outside
                            if s1 == s3:
                                idx_bra = to_subspace_idx(s2)
                                idx_ket = to_subspace_idx(s4)
                                bra_occ = [0, 0, 0, 0]
                                bra_occ[idx_bra] = 1
                                ket_occ = [0, 0, 0, 0]
                                ket_occ[idx_ket] = 1
                                bra_idx = bra_occ[0] * 8 + bra_occ[1] * 4 + bra_occ[2] * 2 + bra_occ[3]
                                ket_idx = ket_occ[0] * 8 + ket_occ[1] * 4 + ket_occ[2] * 2 + ket_occ[3]
                                rho_ij[bra_idx, ket_idx] += total_amp

                        elif in1 and not in2 and not in3 and in4:
                            # Electron 1 in (bra), electron 2 in (ket), swapped outside
                            if s2 == s3:  # outside electrons match
                                idx_bra = to_subspace_idx(s1)
                                idx_ket = to_subspace_idx(s4)
                                bra_occ = [0, 0, 0, 0]
                                bra_occ[idx_bra] = 1
                                ket_occ = [0, 0, 0, 0]
                                ket_occ[idx_ket] = 1
                                bra_idx = bra_occ[0] * 8 + bra_occ[1] * 4 + bra_occ[2] * 2 + bra_occ[3]
                                ket_idx = ket_occ[0] * 8 + ket_occ[1] * 4 + ket_occ[2] * 2 + ket_occ[3]
                                # Sign from moving electrons: s2 was electron 2 in bra,
                                # s3 was electron 1 in ket -- crossing gives sign -1
                                rho_ij[bra_idx, ket_idx] -= total_amp

                        elif not in1 and in2 and in3 and not in4:
                            # Electron 2 in (bra), electron 1 in (ket), swapped outside
                            if s1 == s4:
                                idx_bra = to_subspace_idx(s2)
                                idx_ket = to_subspace_idx(s3)
                                bra_occ = [0, 0, 0, 0]
                                bra_occ[idx_bra] = 1
                                ket_occ = [0, 0, 0, 0]
                                ket_occ[idx_ket] = 1
                                bra_idx = bra_occ[0] * 8 + bra_occ[1] * 4 + bra_occ[2] * 2 + bra_occ[3]
                                ket_idx = ket_occ[0] * 8 + ket_occ[1] * 4 + ket_occ[2] * 2 + ket_occ[3]
                                rho_ij[bra_idx, ket_idx] -= total_amp

                        elif not in1 and not in2 and not in3 and not in4:
                            # Both outside: contributes to |vac><vac|
                            if s1 == s3 and s2 == s4:
                                rho_ij[0, 0] += total_amp
                            elif s1 == s4 and s2 == s3:
                                rho_ij[0, 0] -= total_amp

        return rho_ij

    # Compute mutual information for all orbital pairs
    print("\n  Computing mutual information matrix...")
    MI = np.zeros((n_spatial, n_spatial))

    # Only compute for orbitals with non-negligible entropy
    significant = np.where(s_single > 1e-8)[0]
    print(f"  Significant orbitals: {len(significant)} / {n_spatial}")

    for ii, orb_i in enumerate(significant):
        for jj, orb_j in enumerate(significant):
            if orb_j <= orb_i:
                continue

            rho_ij = two_orbital_entropy(orb_i, orb_j)

            # Check hermiticity and trace
            rho_ij = (rho_ij + rho_ij.T) / 2  # symmetrize for numerical stability

            # Compute von Neumann entropy
            evals_ij = np.linalg.eigvalsh(rho_ij)
            evals_ij = np.maximum(evals_ij, 0)
            s_ij = 0.0
            for ev in evals_ij:
                if ev > 1e-15:
                    s_ij -= ev * np.log(ev)

            # Mutual information
            I_ij = s_single[orb_i] + s_single[orb_j] - s_ij
            MI[orb_i, orb_j] = I_ij
            MI[orb_j, orb_i] = I_ij

    # Print top mutual information pairs
    print(f"\n  Top mutual information pairs:")
    mi_pairs = []
    for i in range(n_spatial):
        for j in range(i + 1, n_spatial):
            if MI[i, j] > 1e-8:
                mi_pairs.append((i, j, MI[i, j]))
    mi_pairs.sort(key=lambda x: x[2], reverse=True)

    for i, j, val in mi_pairs[:15]:
        oi = orbitals[i]
        oj = orbitals[j]
        print(f"    I(({oi[0]},{oi[1]},{oi[2]}), ({oj[0]},{oj[1]},{oj[2]})) = {val:.6f}")

    # Build orbital labels
    orbital_labels = [f"({n},{l},{m})" for n, l, m in orbitals]

    results = {
        'Z': Z,
        'n_max': n_max,
        'n_spatial': n_spatial,
        'single_orbital_entropy': [float(x) for x in s_single],
        'mutual_information_matrix': MI.tolist(),
        'orbital_labels': orbital_labels,
        'top_mutual_info_pairs': [
            {'orbital_i': list(orbitals[i]),
             'orbital_j': list(orbitals[j]),
             'I_ij': float(val)}
            for i, j, val in mi_pairs[:20]
        ],
    }

    return results


# =============================================================================
# Part 4: Graph-native vs standard CI comparison
# =============================================================================

def part4_graph_vs_standard(n_max=4, Z=2):
    """Compare entanglement between graph-native and standard CI."""
    print("\n" + "=" * 70)
    print(f"PART 4: GRAPH-NATIVE vs STANDARD CI (Z={Z}, n_max={n_max})")
    print("=" * 70)

    E_exact = -2.903724377034

    results = {}

    # (a) Graph-native CI
    print("\n  (a) Graph-native CI...")
    H_graph, configs_g, orbitals_g, n_spatial_g = build_graph_native_fci_float_z(Z, n_max)
    evals_g, evecs_g = np.linalg.eigh(H_graph)
    E0_g = evals_g[0]
    ci_g = evecs_g[:, 0]
    rho_g = build_1rdm_from_singlet_ci(ci_g, configs_g, n_spatial_g)
    ent_g = compute_entanglement_measures(rho_g)
    error_g = abs((E0_g - E_exact) / E_exact) * 100

    print(f"    E0 = {E0_g:.8f} Ha, error = {error_g:.4f}%")
    print(f"    S = {ent_g['von_neumann_entropy']:.6f}")
    print(f"    Top occ: {ent_g['occupation_numbers'][:5]}")

    results['graph_native'] = {
        'E0': float(E0_g),
        'error_pct': float(error_g),
        'von_neumann_entropy': float(ent_g['von_neumann_entropy']),
        'occupation_numbers': [float(x) for x in ent_g['occupation_numbers']],
        'single_orbital_entropy': [float(x) for x in ent_g['single_orbital_entropy']],
    }

    # (b) Standard CI with k=Z (no optimization)
    print("\n  (b) Standard CI with k=Z...")
    orbitals_s = _build_orbital_basis(n_max)
    n_spatial_s = len(orbitals_s)
    configs_s = []
    for i in range(n_spatial_s):
        for j in range(i, n_spatial_s):
            if orbitals_s[i][2] + orbitals_s[j][2] == 0:
                configs_s.append((i, j))

    H_std = build_fci_matrix(Z=Z, n_max=n_max, k_orb=float(Z))
    evals_s, evecs_s = np.linalg.eigh(H_std)
    E0_s = evals_s[0]
    ci_s = evecs_s[:, 0]
    rho_s = build_1rdm_from_singlet_ci(ci_s, configs_s, n_spatial_s)
    ent_s = compute_entanglement_measures(rho_s)
    error_s = abs((E0_s - E_exact) / E_exact) * 100

    print(f"    E0 = {E0_s:.8f} Ha, error = {error_s:.4f}%")
    print(f"    S = {ent_s['von_neumann_entropy']:.6f}")
    print(f"    Top occ: {ent_s['occupation_numbers'][:5]}")

    results['standard_k_Z'] = {
        'E0': float(E0_s),
        'error_pct': float(error_s),
        'k_orb': float(Z),
        'von_neumann_entropy': float(ent_s['von_neumann_entropy']),
        'occupation_numbers': [float(x) for x in ent_s['occupation_numbers']],
        'single_orbital_entropy': [float(x) for x in ent_s['single_orbital_entropy']],
    }

    # (c) Standard CI with k=k_opt (variational)
    print("\n  (c) Standard CI with k=k_opt (variational)...")
    B, C = build_fci_polynomial(Z=Z, n_max=n_max)
    k_opt, E_opt = solve_variational_fast(Z=Z, n_max=n_max, B=B, C=C)

    H_opt = B * k_opt + C * k_opt ** 2
    evals_o, evecs_o = np.linalg.eigh(H_opt)
    E0_o = evals_o[0]
    ci_o = evecs_o[:, 0]
    rho_o = build_1rdm_from_singlet_ci(ci_o, configs_s, n_spatial_s)
    ent_o = compute_entanglement_measures(rho_o)
    error_o = abs((E0_o - E_exact) / E_exact) * 100

    print(f"    k_opt = {k_opt:.6f}")
    print(f"    E0 = {E0_o:.8f} Ha, error = {error_o:.4f}%")
    print(f"    S = {ent_o['von_neumann_entropy']:.6f}")
    print(f"    Top occ: {ent_o['occupation_numbers'][:5]}")

    results['standard_k_opt'] = {
        'E0': float(E0_o),
        'error_pct': float(error_o),
        'k_opt': float(k_opt),
        'von_neumann_entropy': float(ent_o['von_neumann_entropy']),
        'occupation_numbers': [float(x) for x in ent_o['occupation_numbers']],
        'single_orbital_entropy': [float(x) for x in ent_o['single_orbital_entropy']],
    }

    # Compare entanglement spectra
    print(f"\n  Comparison summary:")
    print(f"    Graph-native: S={ent_g['von_neumann_entropy']:.6f}, error={error_g:.4f}%")
    print(f"    Standard k=Z: S={ent_s['von_neumann_entropy']:.6f}, error={error_s:.4f}%")
    print(f"    Standard k_opt: S={ent_o['von_neumann_entropy']:.6f}, error={error_o:.4f}%")
    print(f"    S_graph/S_std_Z = {ent_g['von_neumann_entropy']/max(ent_s['von_neumann_entropy'], 1e-15):.4f}")
    print(f"    S_graph/S_opt = {ent_g['von_neumann_entropy']/max(ent_o['von_neumann_entropy'], 1e-15):.4f}")

    return results


# =============================================================================
# Plotting functions
# =============================================================================

def plot_entanglement_vs_Z(z_results, save_path):
    """Plot S(Z) curve with Z_c marked."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    Z_vals = [z_results[k]['Z'] for k in sorted(z_results.keys(), key=lambda x: float(x))]
    S_vals = [z_results[k]['von_neumann_entropy'] for k in sorted(z_results.keys(), key=lambda x: float(x))]

    # Left: S vs Z
    ax = axes[0]
    ax.plot(Z_vals, S_vals, 'bo-', linewidth=2, markersize=8)
    ax.axvline(x=1.84, color='red', linestyle='--', alpha=0.7, label='$Z_c \\approx 1.84$')
    ax.set_xlabel('Nuclear charge Z', fontsize=12)
    ax.set_ylabel('Von Neumann entropy S', fontsize=12)
    ax.set_title('Entanglement entropy vs Z (n_max=4)', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Annotate key points
    for k in sorted(z_results.keys(), key=lambda x: float(x)):
        Z = z_results[k]['Z']
        S = z_results[k]['von_neumann_entropy']
        if Z in [1, 2, 10]:
            ax.annotate(f'Z={Z}\nS={S:.4f}', (Z, S),
                        textcoords="offset points", xytext=(10, 10),
                        fontsize=9, arrowprops=dict(arrowstyle='->', color='gray'))

    # Right: S*Z^2 vs Z (test 1/Z^2 scaling)
    ax = axes[1]
    SZ2 = [s * z ** 2 for s, z in zip(S_vals, Z_vals)]
    ax.plot(Z_vals, SZ2, 'rs-', linewidth=2, markersize=8)
    ax.axvline(x=1.84, color='red', linestyle='--', alpha=0.7, label='$Z_c \\approx 1.84$')
    ax.set_xlabel('Nuclear charge Z', fontsize=12)
    ax.set_ylabel('$S \\times Z^2$', fontsize=12)
    ax.set_title('Scaled entropy $S \\times Z^2$ (test for $1/Z^2$ decay)', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_z_threelayer(z_results, save_path):
    """Plot three-layer fractions vs Z."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    sorted_keys = sorted(z_results.keys(), key=lambda x: float(x))
    Z_vals = [z_results[k]['Z'] for k in sorted_keys]
    frac_vee = [z_results[k]['frac_vee_pct'] for k in sorted_keys]
    frac_h1 = [z_results[k]['frac_h1_offdiag_pct'] for k in sorted_keys]
    frac_rat = [z_results[k]['frac_rational_pct'] for k in sorted_keys]

    # Left: stacked area
    ax = axes[0]
    ax.fill_between(Z_vals, 0, frac_rat, alpha=0.3, color='green', label='Rational (diag)')
    ax.fill_between(Z_vals, frac_rat,
                    [r + h for r, h in zip(frac_rat, frac_h1)],
                    alpha=0.3, color='blue', label='Topological (graph h1)')
    ax.fill_between(Z_vals,
                    [r + h for r, h in zip(frac_rat, frac_h1)],
                    [r + h + v for r, h, v in zip(frac_rat, frac_h1, frac_vee)],
                    alpha=0.3, color='red', label='Correlation (off-diag V_ee)')
    ax.axvline(x=1.84, color='black', linestyle='--', alpha=0.5, label='$Z_c$')
    ax.set_xlabel('Nuclear charge Z', fontsize=12)
    ax.set_ylabel('Fraction of total entropy (%)', fontsize=12)
    ax.set_title('Three-layer entropy fractions vs Z', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 110)

    # Right: V_ee fraction vs Z
    ax = axes[1]
    ax.plot(Z_vals, frac_vee, 'ro-', linewidth=2, markersize=8, label='V_ee correlation')
    ax.plot(Z_vals, frac_h1, 'bs-', linewidth=2, markersize=8, label='Graph topology')
    ax.axvline(x=1.84, color='black', linestyle='--', alpha=0.5, label='$Z_c$')
    ax.set_xlabel('Nuclear charge Z', fontsize=12)
    ax.set_ylabel('Fraction of total entropy (%)', fontsize=12)
    ax.set_title('V_ee vs graph topology fractions', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_mutual_info(mi_results, save_path):
    """Plot mutual information heatmap."""
    MI = np.array(mi_results['mutual_information_matrix'])
    labels = mi_results['orbital_labels']
    n = len(labels)

    # Find significant orbitals
    row_max = np.max(MI, axis=1)
    significant = np.where(row_max > 1e-6)[0]

    if len(significant) < 2:
        print("  Not enough significant orbitals for MI heatmap")
        return

    MI_sub = MI[np.ix_(significant, significant)]
    labels_sub = [labels[i] for i in significant]

    fig, ax = plt.subplots(figsize=(10, 8))

    # Use log scale for color
    MI_plot = MI_sub.copy()
    MI_plot[MI_plot < 1e-10] = 1e-10

    im = ax.imshow(MI_plot, cmap='YlOrRd',
                   norm=LogNorm(vmin=max(1e-8, MI_plot[MI_plot > 0].min()),
                                vmax=MI_plot.max()),
                   interpolation='nearest')

    ax.set_xticks(range(len(labels_sub)))
    ax.set_xticklabels(labels_sub, rotation=90, fontsize=8)
    ax.set_yticks(range(len(labels_sub)))
    ax.set_yticklabels(labels_sub, fontsize=8)

    plt.colorbar(im, ax=ax, label='Mutual information I(i,j)')
    ax.set_title('Orbital mutual information map (He, n_max=4)', fontsize=13)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_config_contributions(config_results, save_path):
    """Plot configuration-type contributions to entanglement."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    types = config_results['config_types']
    sorted_types = sorted(types.items(), key=lambda x: x[1]['weight'], reverse=True)

    # Left: CI weights by config type
    ax = axes[0]
    names = [k for k, v in sorted_types[:12]]
    weights = [v['weight'] for k, v in sorted_types[:12]]
    ax.barh(range(len(names)), weights, color='steelblue')
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('Total CI weight $\\sum |c_I|^2$', fontsize=12)
    ax.set_title('CI weight by configuration type', fontsize=13)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')

    # Right: 1-RDM off-diagonal Frobenius norm (entanglement contribution)
    ax = axes[1]
    frobs = [v['rho_offdiag_frobenius'] for k, v in sorted_types[:12]]
    ax.barh(range(len(names)), frobs, color='coral')
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('1-RDM off-diag Frobenius norm', fontsize=12)
    ax.set_title('Entanglement contribution by config type', fontsize=13)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_graph_vs_standard(comparison_results, save_path):
    """Plot graph-native vs standard CI entanglement comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    methods = ['graph_native', 'standard_k_Z', 'standard_k_opt']
    labels = ['Graph-native\n($\\kappa=-1/16$)', 'Standard\n($k=Z$)',
              'Standard\n($k=k_{opt}$)']

    # Left: entropy comparison
    ax = axes[0]
    S_vals = [comparison_results[m]['von_neumann_entropy'] for m in methods]
    colors = ['#e74c3c', '#3498db', '#2ecc71']
    bars = ax.bar(range(len(methods)), S_vals, color=colors, width=0.6)
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel('Von Neumann entropy S', fontsize=12)
    ax.set_title('Entanglement entropy comparison', fontsize=13)
    ax.grid(True, alpha=0.3, axis='y')
    for bar, val in zip(bars, S_vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.001,
                f'{val:.5f}', ha='center', va='bottom', fontsize=10)

    # Right: occupation number spectrum
    ax = axes[1]
    for m, label, color in zip(methods, ['Graph', 'Std k=Z', 'Std k_opt'], colors):
        occ = comparison_results[m]['occupation_numbers']
        sig = [x for x in occ if x > 1e-10]
        ax.semilogy(range(len(sig)), sig, 'o-', color=color,
                    markersize=5, linewidth=1.5, label=label)

    ax.set_xlabel('Natural orbital index', fontsize=12)
    ax.set_ylabel('Occupation number', fontsize=12)
    ax.set_title('Entanglement spectrum comparison', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

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
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    print("=" * 70)
    print("ENTANGLEMENT DECOMPOSITION OF He-LIKE GROUND STATES")
    print("=" * 70)

    all_results = {}

    # Part 1: Z-dependence
    z_results = part1_z_dependence(n_max=4)
    all_results['part1_z_dependence'] = z_results

    # Part 2: Configuration-resolved
    config_results = part2_config_resolved(n_max=4, Z=2)
    all_results['part2_config_resolved'] = config_results

    # Part 3: Mutual information
    mi_results = part3_mutual_information(n_max=4, Z=2)
    all_results['part3_mutual_information'] = mi_results

    # Part 4: Graph vs standard
    comparison_results = part4_graph_vs_standard(n_max=4, Z=2)
    all_results['part4_graph_vs_standard'] = comparison_results

    # Save results
    output_path = os.path.join(data_dir, 'entanglement_decomposition.json')
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved: {output_path}")

    # Generate plots
    print("\n\n--- GENERATING PLOTS ---")
    plot_entanglement_vs_Z(z_results,
                            os.path.join(plot_dir, 'entanglement_vs_Z.png'))
    plot_z_threelayer(z_results,
                      os.path.join(plot_dir, 'entanglement_z_threelayer.png'))
    plot_mutual_info(mi_results,
                      os.path.join(plot_dir, 'entanglement_mutual_info.png'))
    plot_config_contributions(config_results,
                               os.path.join(plot_dir, 'entanglement_config_contributions.png'))
    plot_graph_vs_standard(comparison_results,
                            os.path.join(plot_dir, 'entanglement_graph_vs_standard.png'))

    # Final summary
    print("\n\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)

    print("\n--- Part 1: Z-dependence ---")
    print(f"{'Z':>6} {'S':>10} {'frac_Vee':>10} {'frac_h1':>10} {'error%':>10}")
    for k in sorted(z_results.keys(), key=lambda x: float(x)):
        r = z_results[k]
        err = f"{r['error_pct']:.2f}" if r['error_pct'] is not None else "N/A"
        print(f"{r['Z']:>6.2f} {r['von_neumann_entropy']:>10.6f} "
              f"{r['frac_vee_pct']:>9.1f}% {r['frac_h1_offdiag_pct']:>9.1f}% "
              f"{err:>10}")

    # Scaling fit
    Z_vals = [z_results[k]['Z'] for k in sorted(z_results.keys(), key=lambda x: float(x))]
    S_vals = [z_results[k]['von_neumann_entropy'] for k in sorted(z_results.keys(), key=lambda x: float(x))]
    # Fit S = a / Z^b for Z >= 2
    mask = [z >= 2 for z in Z_vals]
    Z_fit = [z for z, m in zip(Z_vals, mask) if m]
    S_fit = [s for s, m in zip(S_vals, mask) if m]
    if len(Z_fit) >= 3:
        log_Z = np.log(np.array(Z_fit))
        log_S = np.log(np.array(S_fit))
        coeffs = np.polyfit(log_Z, log_S, 1)
        print(f"\n  Scaling fit (Z >= 2): S ~ Z^{coeffs[0]:.3f}")
        print(f"  (Expected: ~Z^{-2} if entanglement = relative correlation strength)")

    print("\n--- Part 2: Config-resolved (top-5 types) ---")
    types = config_results['config_types']
    sorted_types = sorted(types.items(), key=lambda x: x[1]['weight'], reverse=True)
    for k, v in sorted_types[:5]:
        print(f"  {k}: weight={v['weight']:.6f}, frob_offdiag={v['rho_offdiag_frobenius']:.6f}")

    print("\n--- Part 3: Top mutual information pairs ---")
    for p in mi_results['top_mutual_info_pairs'][:8]:
        oi = tuple(p['orbital_i'])
        oj = tuple(p['orbital_j'])
        print(f"  I({oi}, {oj}) = {p['I_ij']:.6f}")

    print("\n--- Part 4: Graph vs Standard ---")
    for method in ['graph_native', 'standard_k_Z', 'standard_k_opt']:
        r = comparison_results[method]
        print(f"  {method:20s}: S={r['von_neumann_entropy']:.6f}, error={r['error_pct']:.4f}%")

    print("\nDone.")
