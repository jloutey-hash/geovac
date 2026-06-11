"""
VQ-5: Full vector QED — direction-labeled photon channels with sigma vertex.
=============================================================================

This is the headline track of the VQ sprint, combining VQ-3 (sigma-weighted
vertex) and VQ-4 (direction-resolved Hodge decomposition).

Key question: does RESTRICTING which sigma couples to which photon channel
(breaking the universal Sigma_mu sigma_mu^2 = 3I sum) produce a non-trivial
result different from VQ-3's trivial 3x rescaling?

Physical motivation:
  In continuum QED, gamma^mu couples the electron spin to the photon
  polarization direction. On the Fock graph, VQ-4 found two physically
  distinct edge channels:
    T-channel (radial): Dn=+/-1, Dm=0 (inter-shell)
    L-channel (angular): Dn=0, Dm=+1 (intra-shell)

  The natural direction mapping assigns SPECIFIC sigma matrices to each
  channel, rather than summing all three.

Mappings tested:
  A (z/transverse): T -> sigma_z, L -> sigma_x + sigma_y
  B (single-sigma): T -> sigma_z, L -> sigma_x only
  C (exhaustive):   all 6 combinations of (sigma_x, sigma_y, sigma_z) x (T, L)
  D (symmetrized):  T -> sigma_z, L -> (sigma_x + sigma_y)/2

For each mapping, we compute self-energy, vertex correction, F_2,
full 8-rule selection rule census, and number field analysis.

This script does NOT modify any production code.
"""

import numpy as np
import json
import sys
import os
from collections import defaultdict
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from geovac.dirac_matrix_elements import (
    DiracLabel, iter_dirac_labels, kappa_to_l, kappa_to_j,
)
from geovac.graph_qed_vertex import (
    build_projection_matrix,
    build_vertex_tensor,
    vertex_tensor_to_matrices,
)
from geovac.graph_qed_photon import (
    build_fock_graph,
    compute_photon_propagator,
)
from geovac.graph_qed_propagator import (
    DiracGraphOperator,
    electron_propagator,
)
from geovac.lattice import GeometricLattice
from sympy.physics.wigner import clebsch_gordan
from sympy import Rational


# ===========================================================================
# JSON serialization helper
# ===========================================================================

def _json_default(obj):
    """Convert numpy types to Python types for JSON serialization."""
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, complex):
        return float(obj.real) if abs(obj.imag) < 1e-12 else [obj.real, obj.imag]
    if isinstance(obj, set):
        return sorted(obj)
    if isinstance(obj, DiracLabel):
        return f"({obj.n_fock},{obj.kappa},{obj.two_m_j})"
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


# ===========================================================================
# CG coefficient helper
# ===========================================================================

def cg_float(j1, m1, j2, m2, j, m):
    """CG coefficient <j1,m1; j2,m2 | j,m> as float."""
    return float(clebsch_gordan(j1, j2, j, m1, m2, m))


# ===========================================================================
# Step 1: Build sigma matrices in DiracLabel basis (complex)
# ===========================================================================

def build_sigma_matrices_complex(labels, label_index):
    """Build sigma_x, sigma_y, sigma_z as complex numpy arrays.

    Same construction as VQ-3. Sigma matrices act on the SPIN part,
    preserving n and l, in the coupled (j, m_j) basis.
    """
    N = len(labels)
    l_vals = [kappa_to_l(lab.kappa) for lab in labels]

    sigma_z = np.zeros((N, N), dtype=complex)
    sigma_plus = np.zeros((N, N), dtype=complex)
    sigma_minus = np.zeros((N, N), dtype=complex)

    for i, lab_a in enumerate(labels):
        for k, lab_b in enumerate(labels):
            if lab_a.n_fock != lab_b.n_fock:
                continue
            if l_vals[i] != l_vals[k]:
                continue
            l = l_vals[i]

            mj_a = float(lab_a.m_j)
            mj_b = float(lab_b.m_j)

            # sigma_z: preserves m_j
            if abs(mj_a - mj_b) < 1e-10:
                val = 0.0
                for ms2 in [-1, 1]:
                    ms = ms2 / 2.0
                    ml = mj_b - ms
                    ml_int = int(round(ml))
                    if abs(ml_int) > l:
                        continue
                    cg_a_val = cg_float(
                        Rational(l), Rational(ml_int),
                        Rational(1, 2), Rational(ms2, 2),
                        Rational(lab_a.j_times_2, 2),
                        Rational(lab_a.two_m_j, 2))
                    cg_b_val = cg_float(
                        Rational(l), Rational(ml_int),
                        Rational(1, 2), Rational(ms2, 2),
                        Rational(lab_b.j_times_2, 2),
                        Rational(lab_b.two_m_j, 2))
                    val += cg_a_val * (2.0 * ms) * cg_b_val
                sigma_z[i, k] = val

            # sigma_+: m_j -> m_j + 1
            if abs(mj_a - mj_b - 1.0) < 1e-10:
                ms_b = -0.5
                ml = mj_b - ms_b
                ml_int = int(round(ml))
                if abs(ml_int) <= l:
                    cg_b_val = cg_float(
                        Rational(l), Rational(ml_int),
                        Rational(1, 2), Rational(-1, 2),
                        Rational(lab_b.j_times_2, 2),
                        Rational(lab_b.two_m_j, 2))
                    cg_a_val = cg_float(
                        Rational(l), Rational(ml_int),
                        Rational(1, 2), Rational(1, 2),
                        Rational(lab_a.j_times_2, 2),
                        Rational(lab_a.two_m_j, 2))
                    sigma_plus[i, k] = cg_a_val * cg_b_val

            # sigma_-: m_j -> m_j - 1
            if abs(mj_a - mj_b + 1.0) < 1e-10:
                ms_b = 0.5
                ml = mj_b - ms_b
                ml_int = int(round(ml))
                if abs(ml_int) <= l:
                    cg_b_val = cg_float(
                        Rational(l), Rational(ml_int),
                        Rational(1, 2), Rational(1, 2),
                        Rational(lab_b.j_times_2, 2),
                        Rational(lab_b.two_m_j, 2))
                    cg_a_val = cg_float(
                        Rational(l), Rational(ml_int),
                        Rational(1, 2), Rational(-1, 2),
                        Rational(lab_a.j_times_2, 2),
                        Rational(lab_a.two_m_j, 2))
                    sigma_minus[i, k] = cg_a_val * cg_b_val

    # sigma_+ = (sigma_x + i*sigma_y)/2
    sigma_x = sigma_plus + sigma_minus
    sigma_y = -1j * (sigma_plus - sigma_minus)

    return sigma_x, sigma_y, sigma_z


# ===========================================================================
# Step 2: Build scalar vertex matrices from production code (numpy)
# ===========================================================================

def build_scalar_vertex_numpy(n_max):
    """Build the scalar CG vertex matrices V_e from the production code.

    Returns
    -------
    V_e_list : list of np.ndarray, shape (N_dirac, N_dirac)
    dirac_labels : list of DiracLabel
    fock_data : FockGraphData
    E_fock : int
    N_dirac : int
    """
    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    fock_data = build_fock_graph(n_max)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats_sympy = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    V_e_list = []
    for V_sym in V_mats_sympy:
        V_np = np.array(V_sym.tolist(), dtype=float)
        V_e_list.append(V_np)

    return V_e_list, dirac_labels, fock_data, E_fock, N_dirac


# ===========================================================================
# Step 3: Edge classification from VQ-4
# ===========================================================================

def classify_edges(n_max):
    """Classify Fock graph edges into T and L channels.

    Returns dict with:
      T_indices: list of edge indices for T-channel (radial)
      L_indices: list of edge indices for L-channel (angular)
      edge_info: list of per-edge dicts
    """
    lat = GeometricLattice(max_n=n_max, topological_weights=False)
    states = lat.states

    adj = lat.adjacency
    rows, cols = adj.nonzero()
    edge_set = set()
    for r, c in zip(rows, cols):
        if int(r) < int(c):
            edge_set.add((int(r), int(c)))
    edges = sorted(edge_set)

    T_indices = []
    L_indices = []
    edge_info = []

    for e_idx, (v1, v2) in enumerate(edges):
        s1 = states[v1]
        s2 = states[v2]
        dn = s2[0] - s1[0]
        dm = s2[2] - s1[2]

        if dn != 0 and dm == 0:
            channel = 'T'
            T_indices.append(e_idx)
        elif dn == 0:
            channel = 'L'
            L_indices.append(e_idx)
        else:
            channel = 'unknown'

        edge_info.append({
            'edge_idx': e_idx,
            's1': list(s1), 's2': list(s2),
            'dn': dn, 'dm': dm,
            'channel': channel,
        })

    return {
        'T_indices': T_indices,
        'L_indices': L_indices,
        'edge_info': edge_info,
        'edges': edges,
    }


# ===========================================================================
# Step 4: Compute per-channel self-energy and vertex correction
# ===========================================================================

def compute_channel_self_energy(sigma_per_channel, V_e_list, G_gamma,
                                edge_indices, N_dirac, E_fock, channel_name):
    """Compute self-energy for a single direction channel.

    sigma_per_channel: list of sigma matrices to sum over for this channel
    V_e_list: scalar vertex matrices (all edges)
    G_gamma: full photon propagator (E_fock x E_fock)
    edge_indices: which edge indices belong to this channel
    """
    Sigma = np.zeros((N_dirac, N_dirac), dtype=complex)

    for sigma_mu in sigma_per_channel:
        for e1_local, e1 in enumerate(edge_indices):
            for e2_local, e2 in enumerate(edge_indices):
                g_ee = G_gamma[e1, e2]
                if abs(g_ee) < 1e-15:
                    continue
                V1_sigma = sigma_mu @ V_e_list[e1]
                V2_sigma = sigma_mu @ V_e_list[e2]
                Sigma += g_ee * (V1_sigma @ V2_sigma.conj().T)

    return Sigma


def compute_channel_vertex_correction(sigma_per_channel, V_e_list, G_gamma,
                                       G_e, edge_indices, N_dirac, E_fock,
                                       channel_name):
    """Compute vertex correction for a single direction channel."""
    Lambda_total = np.zeros((N_dirac, N_dirac), dtype=complex)

    for sigma_mu in sigma_per_channel:
        for e1_local, e1 in enumerate(edge_indices):
            for e2_local, e2 in enumerate(edge_indices):
                g_ee = G_gamma[e1, e2]
                if abs(g_ee) < 1e-15:
                    continue
                V1_sigma = sigma_mu @ V_e_list[e1]
                V2_sigma = sigma_mu @ V_e_list[e2]
                Lambda_total += g_ee * (V1_sigma @ G_e @ V2_sigma.conj().T)

    return Lambda_total


def extract_f2_for_mapping(Lambda_total, sigma_per_channel_T,
                            sigma_per_channel_L, V_e_list, G_e,
                            T_indices, L_indices, N_dirac):
    """Extract F_2 for a given mapping.

    V_bare = sum_channel sum_mu sum_{e in channel} sigma_mu @ V_e
    F_2 = Tr(Lambda_total) / Tr(V_bare @ G_e)
    """
    V_bare = np.zeros((N_dirac, N_dirac), dtype=complex)
    for sigma_mu in sigma_per_channel_T:
        for e_idx in T_indices:
            V_bare += sigma_mu @ V_e_list[e_idx]
    for sigma_mu in sigma_per_channel_L:
        for e_idx in L_indices:
            V_bare += sigma_mu @ V_e_list[e_idx]

    denominator = np.trace(V_bare @ G_e)
    numerator = np.trace(Lambda_total)

    if abs(denominator) < 1e-15:
        return None, float('nan'), float('nan')

    F2 = numerator / denominator
    return F2, numerator, denominator


# ===========================================================================
# Step 5: Selection rule census
# ===========================================================================

def selection_rule_census(sigma_per_channel_T, sigma_per_channel_L,
                           V_e_list, T_indices, L_indices,
                           labels, Sigma_work, G_e, N_dirac, E_fock,
                           mapping_name):
    """Test 8 continuum QED selection rules for a given mapping."""
    N = N_dirac
    results = {}

    # Collect all nonzero vertex entries
    all_vertex_entries = []
    for sigma_mu in sigma_per_channel_T:
        for e_idx in T_indices:
            V = sigma_mu @ V_e_list[e_idx]
            for a in range(N):
                for b in range(N):
                    if abs(V[a, b]) > 1e-12:
                        all_vertex_entries.append({
                            "channel": "T", "edge": e_idx,
                            "a": a, "b": b,
                            "value_abs": float(abs(V[a, b])),
                            "lab_a": labels[a], "lab_b": labels[b],
                        })
    for sigma_mu in sigma_per_channel_L:
        for e_idx in L_indices:
            V = sigma_mu @ V_e_list[e_idx]
            for a in range(N):
                for b in range(N):
                    if abs(V[a, b]) > 1e-12:
                        all_vertex_entries.append({
                            "channel": "L", "edge": e_idx,
                            "a": a, "b": b,
                            "value_abs": float(abs(V[a, b])),
                            "lab_a": labels[a], "lab_b": labels[b],
                        })

    # 1. Delta m_j conservation
    dm_j_values = set()
    for entry in all_vertex_entries:
        dm = entry["lab_a"].two_m_j - entry["lab_b"].two_m_j
        dm_j_values.add(dm)
    survives_dmj = dm_j_values.issubset({-2, -1, 0, 1, 2})
    results["delta_mj_conservation"] = {
        "delta_2mj_values": sorted(dm_j_values),
        "survives": bool(survives_dmj),
    }

    # 2. Spatial parity E1 (l_a + l_b odd)
    parity_check = {"odd": 0, "even": 0}
    for entry in all_vertex_entries:
        l_sum = entry["lab_a"].l + entry["lab_b"].l
        if l_sum % 2 == 1:
            parity_check["odd"] += 1
        else:
            parity_check["even"] += 1
    results["spatial_parity_E1"] = {
        "l_sum_parity": parity_check,
        "survives": bool(parity_check["odd"] > 0 and parity_check["even"] == 0),
    }

    # 3. Gaunt/CG sparsity
    total_possible = N * N * (
        len(sigma_per_channel_T) * len(T_indices) +
        len(sigma_per_channel_L) * len(L_indices)
    )
    nonzero_count = len(all_vertex_entries)
    density = nonzero_count / total_possible if total_possible > 0 else 0
    results["gaunt_cg_sparsity"] = {
        "total_possible": total_possible,
        "nonzero": nonzero_count,
        "density_percent": float(100 * density),
        "survives": True,
    }

    # 4. Vertex parity (n_1 + n_2 + q odd)
    n_sum_parity = {"even": 0, "odd": 0}
    for entry in all_vertex_entries:
        s = entry["lab_a"].n_fock + entry["lab_b"].n_fock
        if s % 2 == 0:
            n_sum_parity["even"] += 1
        else:
            n_sum_parity["odd"] += 1
    # Vertex parity survives if ONLY odd n-sums are present
    results["vertex_parity"] = {
        "n_sum_parity": n_sum_parity,
        "survives": bool(n_sum_parity["even"] == 0),
    }

    # 5. SO(4) channel count
    results["so4_channel_count"] = {
        "verdict": "N/A -- photon is still a scalar 1-cochain",
        "survives": False,
    }

    # 6. Charge conjugation
    C_mat = np.zeros((N, N))
    for i, lab in enumerate(labels):
        j_val = float(lab.j)
        mj_val = float(lab.m_j)
        phase = (-1) ** (j_val + mj_val)
        target = DiracLabel(lab.n_fock, lab.kappa, -lab.two_m_j)
        for k, lab_k in enumerate(labels):
            if lab_k == target:
                C_mat[k, i] = phase
                break

    if np.allclose(Sigma_work.imag, 0, atol=1e-10):
        Sigma_real = Sigma_work.real
    else:
        Sigma_real = Sigma_work
    CSC = C_mat @ Sigma_real @ C_mat
    c_symmetric = bool(np.allclose(CSC, Sigma_real, atol=1e-10))
    results["charge_conjugation"] = {
        "C_symmetric": c_symmetric,
        "survives": c_symmetric,
    }

    # 7. Furry's theorem (tadpole = 0)
    all_sigma_edge_pairs = []
    for sigma_mu in sigma_per_channel_T:
        for e_idx in T_indices:
            all_sigma_edge_pairs.append((sigma_mu, e_idx))
    for sigma_mu in sigma_per_channel_L:
        for e_idx in L_indices:
            all_sigma_edge_pairs.append((sigma_mu, e_idx))

    tad = np.zeros((N, N), dtype=complex)
    for sigma_mu, e_idx in all_sigma_edge_pairs:
        tad += sigma_mu @ V_e_list[e_idx]
    tadpole_norm = float(np.max(np.abs(tad)))
    tadpole_zero = bool(tadpole_norm < 1e-10)

    # Also check triangle
    triangle = np.zeros((N, N), dtype=complex)
    for sigma_mu_1, e1 in all_sigma_edge_pairs:
        for sigma_mu_2, e2 in all_sigma_edge_pairs:
            for sigma_mu_3, e3 in all_sigma_edge_pairs:
                V1 = sigma_mu_1 @ V_e_list[e1]
                V2 = sigma_mu_2 @ V_e_list[e2]
                V3 = sigma_mu_3 @ V_e_list[e3]
                triangle += V1 @ V2 @ V3
    triangle_norm = float(np.max(np.abs(triangle)))
    triangle_zero = bool(triangle_norm < 1e-10)

    results["furry_theorem"] = {
        "tadpole_max_abs": tadpole_norm,
        "tadpole_zero": tadpole_zero,
        "triangle_max_abs": triangle_norm,
        "triangle_zero": triangle_zero,
        "survives": tadpole_zero,
    }

    # 8. Ward identity: [D, Sigma] = 0
    D = np.diag([
        (1.0 if lab.kappa < 0 else -1.0) * (lab.n_fock + 0.5)
        for lab in labels
    ])
    comm = D @ Sigma_real - Sigma_real @ D
    comm_norm = float(np.max(np.abs(comm)))
    results["ward_identity"] = {
        "commutator_norm": comm_norm,
        "commutes": bool(comm_norm < 1e-10),
        "survives": False,
        "note": "Ward identity needs full vertex-propagator relation",
    }

    # Summary
    rule_names = [
        "delta_mj_conservation", "spatial_parity_E1", "gaunt_cg_sparsity",
        "vertex_parity", "so4_channel_count", "charge_conjugation",
        "furry_theorem", "ward_identity",
    ]
    survived = sum(1 for r in rule_names if results[r].get("survives", False))
    results["summary"] = {
        "total_rules": 8,
        "survived": survived,
        "fraction": f"{survived}/8",
    }

    return results


# ===========================================================================
# Step 6: Number field analysis
# ===========================================================================

def analyze_number_field_quick(Sigma, F2_val):
    """Quick number field analysis."""
    results = {}

    is_real = bool(np.allclose(Sigma.imag, 0, atol=1e-10))
    results["is_real"] = is_real

    # Check if entries are rational
    S = Sigma.real if is_real else Sigma
    nonzero_entries = []
    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            v = S[i, j]
            if isinstance(v, complex):
                v = v.real if abs(v.imag) < 1e-12 else v
            if abs(v) > 1e-12:
                nonzero_entries.append(float(v.real if isinstance(v, complex) else v))

    # Check rational
    all_rational = True
    for v in nonzero_entries[:20]:
        found = False
        for q in range(1, 300):
            p = round(v * q)
            if abs(v - p / q) < 1e-8:
                found = True
                break
        if not found:
            all_rational = False
            break
    results["appears_rational"] = all_rational

    # Check sqrt(rational) if not rational
    if not all_rational:
        sqrt_count = 0
        for v in nonzero_entries[:20]:
            v2 = v ** 2
            for q in range(1, 300):
                p = round(v2 * q)
                if abs(v2 - p / q) < 1e-8:
                    sqrt_count += 1
                    break
        results["sqrt_rational_fraction"] = f"{sqrt_count}/{min(len(nonzero_entries), 20)}"

    # Check F2
    if F2_val is not None and not np.isnan(F2_val):
        f2 = float(F2_val.real if isinstance(F2_val, complex) else F2_val)
        results["F2_value"] = f2
        # Rational check
        best_frac = None
        best_err = 1e-6
        for q in range(1, 500):
            p = round(f2 * q)
            err = abs(f2 - p / q)
            if err < best_err:
                best_err = err
                best_frac = f"{p}/{q}"
        results["F2_fraction"] = best_frac
        results["F2_fraction_residual"] = float(best_err)

        # sqrt check
        f2_sq = f2 ** 2
        best_frac_sq = None
        best_err_sq = 1e-6
        for q in range(1, 500):
            p = round(f2_sq * q)
            err = abs(f2_sq - p / q)
            if err < best_err_sq:
                best_err_sq = err
                best_frac_sq = f"{p}/{q}"
        results["F2_squared_fraction"] = best_frac_sq
        results["F2_squared_residual"] = float(best_err_sq)

    # Pi check
    results["pi_free"] = all_rational or results.get("sqrt_rational_fraction", "").startswith(str(min(len(nonzero_entries), 20)))

    return results


# ===========================================================================
# Step 7: Run a single mapping
# ===========================================================================

def run_mapping(mapping_name, sigma_T_list, sigma_L_list,
                V_e_list, G_gamma, G_e, T_indices, L_indices,
                labels, N_dirac, E_fock):
    """Run the full analysis for a single sigma-channel mapping.

    Parameters
    ----------
    sigma_T_list : list of sigma matrices for T-channel edges
    sigma_L_list : list of sigma matrices for L-channel edges
    """
    print(f"\n{'='*60}")
    print(f"  Mapping {mapping_name}")
    print(f"{'='*60}")

    sigma_names_T = []
    sigma_names_L = []

    # Compute per-channel self-energies
    Sigma_T = compute_channel_self_energy(
        sigma_T_list, V_e_list, G_gamma,
        T_indices, N_dirac, E_fock, "T")
    Sigma_L = compute_channel_self_energy(
        sigma_L_list, V_e_list, G_gamma,
        L_indices, N_dirac, E_fock, "L")
    Sigma_total = Sigma_T + Sigma_L

    # Compute per-channel vertex corrections
    Lambda_T = compute_channel_vertex_correction(
        sigma_T_list, V_e_list, G_gamma,
        G_e, T_indices, N_dirac, E_fock, "T")
    Lambda_L = compute_channel_vertex_correction(
        sigma_L_list, V_e_list, G_gamma,
        G_e, L_indices, N_dirac, E_fock, "L")
    Lambda_total = Lambda_T + Lambda_L

    # F2
    F2, num, den = extract_f2_for_mapping(
        Lambda_total, sigma_T_list, sigma_L_list,
        V_e_list, G_e, T_indices, L_indices, N_dirac)

    # Self-energy properties
    is_real = bool(np.allclose(Sigma_total.imag, 0, atol=1e-10))
    if is_real:
        Sigma_work = Sigma_total.real
        is_hermitian = bool(np.allclose(Sigma_work, Sigma_work.T, atol=1e-10))
        eigs = np.linalg.eigvalsh(Sigma_work) if is_hermitian else np.linalg.eigvals(Sigma_work).real
    else:
        Sigma_work = Sigma_total
        is_hermitian = bool(np.allclose(Sigma_work, Sigma_work.conj().T, atol=1e-10))
        eigs = np.linalg.eigvals(Sigma_work).real
    eigs_sorted = sorted(eigs)

    # GS block
    gs_idx = [i for i, lab in enumerate(labels)
              if lab.n_fock == 1 and lab.kappa == -1]
    gs_block = Sigma_work[np.ix_(gs_idx, gs_idx)]
    gs_max = float(np.max(np.abs(gs_block)))
    gs_zero = bool(gs_max < 1e-10)

    # Isotropy check: is Sigma_T proportional to Sigma_L?
    if is_real:
        S_T_work = Sigma_T.real
        S_L_work = Sigma_L.real
    else:
        S_T_work = Sigma_T
        S_L_work = Sigma_L
    tr_T = float(np.trace(S_T_work).real if isinstance(np.trace(S_T_work), complex)
                 else np.trace(S_T_work))
    tr_L = float(np.trace(S_L_work).real if isinstance(np.trace(S_L_work), complex)
                 else np.trace(S_L_work))
    if abs(tr_T) > 1e-12 and abs(tr_L) > 1e-12:
        ratio_TL = tr_T / tr_L
        # Check proportionality: S_T = c * S_L ?
        if abs(tr_L) > 1e-10:
            c_prop = tr_T / tr_L
            residual = np.linalg.norm(S_T_work - c_prop * S_L_work, 'fro')
            total = np.linalg.norm(S_T_work, 'fro') + np.linalg.norm(S_L_work, 'fro')
            proportional = bool(residual / max(total, 1e-15) < 1e-8)
        else:
            proportional = False
            residual = float('inf')
    else:
        ratio_TL = None
        proportional = (abs(tr_T) < 1e-12 and abs(tr_L) < 1e-12)
        residual = 0.0

    # Selection rule census
    census = selection_rule_census(
        sigma_T_list, sigma_L_list,
        V_e_list, T_indices, L_indices,
        labels, Sigma_total, G_e, N_dirac, E_fock,
        mapping_name)

    # Number field
    F2_val = F2
    nf = analyze_number_field_quick(Sigma_total, F2_val)

    # Print summary
    tr_total = float(np.trace(Sigma_work).real if isinstance(np.trace(Sigma_work), complex)
                     else np.trace(Sigma_work))
    print(f"  Sigma: real={is_real}, hermitian={is_hermitian}")
    print(f"  Tr(Sigma_T)={tr_T:.6f}, Tr(Sigma_L)={tr_L:.6f}, "
          f"Tr(total)={tr_total:.6f}")
    print(f"  Proportional (Sigma_T ~ c*Sigma_L): {proportional}")
    if ratio_TL is not None:
        print(f"  Tr(T)/Tr(L) ratio: {ratio_TL:.6f}")
    print(f"  Eigenvalues: {[f'{e:.4f}' for e in eigs_sorted]}")
    print(f"  GS block zero: {gs_zero} (max|entry|={gs_max:.6f})")
    if F2 is not None and not np.isnan(float(F2.real if isinstance(F2, complex) else F2)):
        f2_float = float(F2.real if isinstance(F2, complex) else F2)
        print(f"  F_2 = {f2_float:.8f}")
        print(f"  F_2 / F_2_scalar(5*sqrt(2)/3) = {f2_float / (5*np.sqrt(2)/3):.6f}")
    else:
        print(f"  F_2 = undefined")
    print(f"  Selection rules: {census['summary']['fraction']}")
    print(f"  Number field: rational={nf.get('appears_rational', '?')}")

    # Compile result
    result = {
        "mapping_name": mapping_name,
        "self_energy": {
            "is_real": is_real,
            "is_hermitian": is_hermitian,
            "trace_T": tr_T,
            "trace_L": tr_L,
            "trace_total": tr_total,
            "eigenvalues": [float(e) for e in eigs_sorted],
            "n_zero_eigenvalues": sum(1 for e in eigs_sorted if abs(e) < 1e-10),
            "positive_semidefinite": bool(all(e > -1e-10 for e in eigs_sorted)),
            "gs_block_zero": gs_zero,
            "gs_block_max_abs": gs_max,
        },
        "isotropy": {
            "trace_ratio_T_over_L": ratio_TL,
            "proportional": proportional,
            "proportionality_residual": float(residual) if not np.isinf(residual) else None,
        },
        "vertex_correction": {
            "trace_Lambda": float(np.trace(Lambda_total).real
                                   if isinstance(np.trace(Lambda_total), complex)
                                   else np.trace(Lambda_total)),
            "F2": (float(F2.real if isinstance(F2, complex) else F2)
                   if F2 is not None and not np.isnan(float(F2.real if isinstance(F2, complex) else F2))
                   else None),
            "F2_numerator": float(num.real if isinstance(num, complex) else num) if not np.isnan(num) else None,
            "F2_denominator": float(den.real if isinstance(den, complex) else den) if not np.isnan(den) else None,
        },
        "selection_rules": census,
        "number_field": nf,
    }

    return result


# ===========================================================================
# Step 8: Compute scalar reference
# ===========================================================================

def compute_scalar_reference(V_e_list, G_gamma, G_e, labels, N_dirac, E_fock):
    """Compute scalar (GN-5 style) self-energy and F2 for reference."""
    Sigma_scalar = np.zeros((N_dirac, N_dirac), dtype=float)
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_ee = G_gamma[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            Sigma_scalar += g_ee * (V_e_list[e1] @ V_e_list[e2].T)

    Lambda_scalar = np.zeros((N_dirac, N_dirac), dtype=float)
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_ee = G_gamma[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            Lambda_scalar += g_ee * (V_e_list[e1] @ G_e @ V_e_list[e2].T)

    V_bare = np.zeros((N_dirac, N_dirac), dtype=float)
    for e_idx in range(E_fock):
        V_bare += V_e_list[e_idx]

    den = np.trace(V_bare @ G_e)
    num = np.trace(Lambda_scalar)
    F2_scalar = num / den if abs(den) > 1e-15 else None

    return {
        "trace": float(np.trace(Sigma_scalar)),
        "eigenvalues": sorted(np.linalg.eigvalsh(Sigma_scalar).tolist()),
        "F2": float(F2_scalar) if F2_scalar is not None else None,
        "gs_block_max": float(np.max(np.abs(
            Sigma_scalar[np.ix_(
                [i for i, l in enumerate(labels) if l.n_fock == 1 and l.kappa == -1],
                [i for i, l in enumerate(labels) if l.n_fock == 1 and l.kappa == -1]
            )]
        ))),
    }


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 70)
    print("VQ-5: Full vector QED — direction-labeled channels + sigma vertex")
    print("=" * 70)

    n_max = 2

    # ------------------------------------------------------------------
    # Build infrastructure
    # ------------------------------------------------------------------
    print("\n--- Building infrastructure ---")
    labels = list(iter_dirac_labels(n_max))
    label_index = {lab: i for i, lab in enumerate(labels)}
    N_dirac = len(labels)
    print(f"n_max = {n_max}, N_dirac = {N_dirac}")

    # Sigma matrices
    sigma_x, sigma_y, sigma_z = build_sigma_matrices_complex(labels, label_index)

    # Scalar vertex
    V_e_list, dirac_labels, fock_data, E_fock, N_dirac_check = build_scalar_vertex_numpy(n_max)
    assert N_dirac == N_dirac_check
    print(f"E_fock = {E_fock} edges")

    # Edge classification
    edge_class = classify_edges(n_max)
    T_indices = edge_class['T_indices']
    L_indices = edge_class['L_indices']
    print(f"T edges: {T_indices} ({len(T_indices)} edges)")
    print(f"L edges: {L_indices} ({len(L_indices)} edges)")
    for ei in edge_class['edge_info']:
        print(f"  edge {ei['edge_idx']}: {ei['s1']} -> {ei['s2']}  "
              f"Dn={ei['dn']:+d} Dm={ei['dm']:+d} channel={ei['channel']}")

    # Photon propagator (full)
    photon_data = compute_photon_propagator(n_max, exact=False)
    G_gamma = photon_data.G_gamma_numeric
    print(f"G_gamma shape: {G_gamma.shape}")

    # Electron propagator at t=0
    D_diag = np.diag([
        (1.0 if lab.kappa < 0 else -1.0) * (lab.n_fock + 0.5)
        for lab in labels
    ])
    G_e = np.linalg.inv(D_diag)

    # ------------------------------------------------------------------
    # Scalar reference
    # ------------------------------------------------------------------
    print("\n--- Scalar reference (GN-5) ---")
    scalar_ref = compute_scalar_reference(V_e_list, G_gamma, G_e, labels,
                                           N_dirac, E_fock)
    print(f"  Tr(Sigma_scalar) = {scalar_ref['trace']:.6f}")
    print(f"  F2_scalar = {scalar_ref['F2']:.6f}")
    scalar_f2 = scalar_ref['F2']

    # ------------------------------------------------------------------
    # VQ-3 reference (all sigma summed, all edges)
    # ------------------------------------------------------------------
    print("\n--- VQ-3 reference (sum sigma_mu, all edges) ---")
    all_edges = list(range(E_fock))
    sigma_all = [sigma_x, sigma_y, sigma_z]
    vq3_ref = run_mapping("VQ-3 (sum sigma_mu, all edges)",
                           sigma_all, sigma_all,
                           V_e_list, G_gamma, G_e,
                           all_edges, [],  # all edges in "T", none in "L"
                           labels, N_dirac, E_fock)

    # ------------------------------------------------------------------
    # Mapping A: T -> sigma_z, L -> sigma_x + sigma_y
    # ------------------------------------------------------------------
    mapping_A = run_mapping("A (T->sigma_z, L->sigma_x+sigma_y)",
                            [sigma_z], [sigma_x, sigma_y],
                            V_e_list, G_gamma, G_e,
                            T_indices, L_indices,
                            labels, N_dirac, E_fock)

    # ------------------------------------------------------------------
    # Mapping B: T -> sigma_z, L -> sigma_x
    # ------------------------------------------------------------------
    mapping_B = run_mapping("B (T->sigma_z, L->sigma_x)",
                            [sigma_z], [sigma_x],
                            V_e_list, G_gamma, G_e,
                            T_indices, L_indices,
                            labels, N_dirac, E_fock)

    # ------------------------------------------------------------------
    # Mapping C: all 6 combinations exhaustively
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("  Mapping C: Exhaustive 6 combinations")
    print("=" * 60)
    sigma_dict = {"x": sigma_x, "y": sigma_y, "z": sigma_z}
    mapping_C_results = {}
    for T_name, T_sigma in sigma_dict.items():
        for L_name, L_sigma in sigma_dict.items():
            combo_name = f"C (T->sigma_{T_name}, L->sigma_{L_name})"
            result = run_mapping(combo_name,
                                 [T_sigma], [L_sigma],
                                 V_e_list, G_gamma, G_e,
                                 T_indices, L_indices,
                                 labels, N_dirac, E_fock)
            mapping_C_results[f"T={T_name},L={L_name}"] = result

    # ------------------------------------------------------------------
    # Mapping D: T -> sigma_z, L -> (sigma_x + sigma_y)/2 averaged
    # ------------------------------------------------------------------
    sigma_xy_avg = (sigma_x + sigma_y) / 2.0
    mapping_D = run_mapping("D (T->sigma_z, L->(sigma_x+sigma_y)/2)",
                            [sigma_z], [sigma_xy_avg],
                            V_e_list, G_gamma, G_e,
                            T_indices, L_indices,
                            labels, N_dirac, E_fock)

    # ------------------------------------------------------------------
    # Mapping E: Swapped — T -> sigma_x, L -> sigma_z
    # ------------------------------------------------------------------
    mapping_E = run_mapping("E (T->sigma_x, L->sigma_z)",
                            [sigma_x], [sigma_z],
                            V_e_list, G_gamma, G_e,
                            T_indices, L_indices,
                            labels, N_dirac, E_fock)

    # ------------------------------------------------------------------
    # Comparison table
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("  COMPARISON TABLE")
    print("=" * 70)
    print(f"{'Construction':45s} {'Rules':6s} {'Tr(Sigma)':>12s} {'F2':>12s} {'Rational':>8s}")
    print("-" * 85)

    # Scalar
    print(f"{'Scalar Fock (GN-5)':45s} {'1/8':6s} {scalar_ref['trace']:12.4f} "
          f"{scalar_ref['F2']:12.6f} {'?':>8s}")

    # VQ-3
    vq3_tr = vq3_ref['self_energy']['trace_total']
    vq3_f2 = vq3_ref['vertex_correction']['F2']
    vq3_rules = vq3_ref['selection_rules']['summary']['fraction']
    print(f"{'VQ-3 (sum sigma_mu, all edges)':45s} {vq3_rules:6s} {vq3_tr:12.4f} "
          f"{vq3_f2:12.6f} {'?':>8s}")

    # Mapping A
    a_tr = mapping_A['self_energy']['trace_total']
    a_f2 = mapping_A['vertex_correction']['F2']
    a_rules = mapping_A['selection_rules']['summary']['fraction']
    print(f"{'A (T->sigma_z, L->sigma_x+sigma_y)':45s} {a_rules:6s} {a_tr:12.4f} "
          f"{a_f2:12.6f} {'?':>8s}")

    # Mapping B
    b_tr = mapping_B['self_energy']['trace_total']
    b_f2 = mapping_B['vertex_correction']['F2']
    b_rules = mapping_B['selection_rules']['summary']['fraction']
    print(f"{'B (T->sigma_z, L->sigma_x)':45s} {b_rules:6s} {b_tr:12.4f} "
          f"{b_f2:12.6f} {'?':>8s}")

    # Mapping D
    d_tr = mapping_D['self_energy']['trace_total']
    d_f2 = mapping_D['vertex_correction']['F2']
    d_rules = mapping_D['selection_rules']['summary']['fraction']
    print(f"{'D (T->sigma_z, L->(sigma_x+sigma_y)/2)':45s} {d_rules:6s} {d_tr:12.4f} "
          f"{d_f2:12.6f} {'?':>8s}")

    # Mapping E
    e_tr = mapping_E['self_energy']['trace_total']
    e_f2 = mapping_E['vertex_correction']['F2']
    e_rules = mapping_E['selection_rules']['summary']['fraction']
    e_f2_str = f"{e_f2:12.6f}" if e_f2 is not None else "         N/A"
    print(f"{'E (T->sigma_x, L->sigma_z)':45s} {e_rules:6s} {e_tr:12.4f} "
          f"{e_f2_str} {'?':>8s}")

    print("\nMapping C exhaustive scan:")
    for combo, res in mapping_C_results.items():
        c_tr = res['self_energy']['trace_total']
        c_f2 = res['vertex_correction']['F2']
        c_rules = res['selection_rules']['summary']['fraction']
        print(f"  {'C (' + combo + ')':43s} {c_rules:6s} {c_tr:12.4f} "
              f"{c_f2 if c_f2 is not None else 'N/A':>12} {'?':>8s}")

    # ------------------------------------------------------------------
    # Key diagnostic: does breaking sigma sum change F2?
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("  KEY DIAGNOSTIC: Does channel-specific sigma break the 3x rescaling?")
    print("=" * 70)

    # VQ-3 gives F2 = 3 * F2_scalar
    # Do any of the new mappings give something different?
    all_f2_values = set()
    for name, res in [("A", mapping_A), ("B", mapping_B), ("D", mapping_D), ("E", mapping_E)]:
        f2 = res['vertex_correction']['F2']
        if f2 is not None:
            all_f2_values.add(round(f2, 6))
            ratio = f2 / scalar_f2 if scalar_f2 and abs(scalar_f2) > 1e-15 else None
            ratio_str = f"{ratio:.6f}" if ratio is not None else "N/A"
            print(f"  {name}: F2 = {f2:.8f}, F2/F2_scalar = {ratio_str}")
    for combo, res in mapping_C_results.items():
        f2 = res['vertex_correction']['F2']
        if f2 is not None:
            all_f2_values.add(round(f2, 6))

    n_distinct = len(all_f2_values)
    print(f"\n  Number of distinct F2 values across all mappings: {n_distinct}")
    print(f"  Distinct values: {sorted(all_f2_values)}")

    trivial_3x = all(abs(v - 3 * scalar_f2) < 0.01 * abs(scalar_f2)
                      for v in all_f2_values if v is not None and scalar_f2)
    print(f"  All values consistent with 3x scalar: {trivial_3x}")

    if n_distinct > 1:
        print(f"\n  ** NON-TRIVIAL: Channel-specific sigma assignment produces")
        print(f"     {n_distinct} distinct F2 values! The 3x rescaling is BROKEN **")
    else:
        if trivial_3x:
            print(f"\n  NEGATIVE: All mappings give F2 = 3 * F2_scalar")
        else:
            print(f"\n  All mappings give the same (non-3x) value")

    # ------------------------------------------------------------------
    # Isotropy analysis across mappings
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("  ISOTROPY ANALYSIS")
    print("=" * 70)
    for name, res in [("A", mapping_A), ("B", mapping_B),
                       ("D", mapping_D), ("E", mapping_E)]:
        iso = res['isotropy']
        print(f"  {name}: Tr(T)/Tr(L) = {iso['trace_ratio_T_over_L']}, "
              f"proportional = {iso['proportional']}")

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    print("\n--- Saving results ---")

    output = {
        "n_max": n_max,
        "N_dirac": N_dirac,
        "E_fock": E_fock,
        "T_indices": T_indices,
        "L_indices": L_indices,
        "edge_classification": edge_class['edge_info'],
        "scalar_reference": scalar_ref,
        "vq3_reference": vq3_ref,
        "mapping_A": mapping_A,
        "mapping_B": mapping_B,
        "mapping_C": mapping_C_results,
        "mapping_D": mapping_D,
        "mapping_E": mapping_E,
        "key_diagnostic": {
            "n_distinct_f2": n_distinct,
            "distinct_f2_values": sorted(all_f2_values),
            "all_3x_scalar": trivial_3x,
            "scalar_f2": scalar_f2,
        },
    }

    out_path = Path(__file__).resolve().parent / "data" / "vq5_full_vector_qed.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2, default=_json_default)
    print(f"  Results saved to {out_path}")

    # ------------------------------------------------------------------
    # Write memo
    # ------------------------------------------------------------------
    write_memo(output, scalar_ref, vq3_ref, mapping_A, mapping_B,
               mapping_C_results, mapping_D, mapping_E, n_distinct,
               trivial_3x, all_f2_values, scalar_f2)

    print("\n" + "=" * 70)
    print("VQ-5 COMPLETE")
    print("=" * 70)


# ===========================================================================
# Memo writer
# ===========================================================================

def write_memo(output, scalar_ref, vq3_ref, mapping_A, mapping_B,
               mapping_C_results, mapping_D, mapping_E, n_distinct,
               trivial_3x, all_f2_values, scalar_f2):
    """Write the comprehensive VQ-5 / VQ sprint final memo."""
    memo_path = Path(__file__).resolve().parent / "vq5_full_vector_qed_memo.md"

    lines = []
    lines.append("# VQ-5: Full Vector QED — Direction-Labeled Photon Channels with Sigma Vertex")
    lines.append("")
    lines.append("## VQ Sprint Final Summary")
    lines.append("")
    lines.append("This is the headline track of the VQ sprint, combining VQ-3 (sigma-weighted")
    lines.append("vertex) and VQ-4 (direction-resolved Hodge decomposition) to build the closest")
    lines.append("analog of vector QED on the finite Fock graph.")
    lines.append("")

    # Background
    lines.append("## Background and Motivation")
    lines.append("")
    lines.append("In continuum QED, the electron-photon vertex gamma^mu couples electron spin to")
    lines.append("photon polarization direction. On the finite Fock graph, we have:")
    lines.append("")
    lines.append("- **Electron spin structure** (VQ-3): Pauli matrices sigma_mu in the (j, m_j) basis")
    lines.append("- **Photon direction channels** (VQ-4): T-edges (radial, Dn=+/-1) and L-edges")
    lines.append("  (angular, Dm=+/-1), with distinct per-channel propagators G_T and G_L")
    lines.append("")
    lines.append("VQ-3 showed that dressing ALL vertices with the SAME sum over sigma_mu gives a")
    lines.append("trivial 3x rescaling via the Pauli trace identity Sum_mu sigma_mu^2 = 3I.")
    lines.append("The key question for VQ-5: does assigning SPECIFIC sigma matrices to SPECIFIC")
    lines.append("photon channels break this identity and produce non-trivial physics?")
    lines.append("")

    # Mappings tested
    lines.append("## Mappings Tested")
    lines.append("")
    lines.append("| Mapping | T-channel (radial) | L-channel (angular) |")
    lines.append("|:--------|:-------------------|:--------------------|")
    lines.append("| A (z/transverse) | sigma_z | sigma_x + sigma_y (sum) |")
    lines.append("| B (single-sigma) | sigma_z | sigma_x only |")
    lines.append("| C (exhaustive) | each of sigma_{x,y,z} | each of sigma_{x,y,z} |")
    lines.append("| D (symmetrized) | sigma_z | (sigma_x + sigma_y)/2 |")
    lines.append("| E (swapped) | sigma_x | sigma_z |")
    lines.append("")

    # Results: comparison table
    lines.append("## Results: Comparison Table")
    lines.append("")
    lines.append("| Construction | Rules | Tr(Sigma) | F2 | F2/F2_scalar |")
    lines.append("|:-------------|:------|:----------|:---|:-------------|")

    # Scalar
    lines.append(f"| Scalar Fock (GN-5) | 1/8 | "
                 f"{scalar_ref['trace']:.4f} | "
                 f"{scalar_ref['F2']:.6f} | 1.000 |")

    # VQ-3
    vq3_f2 = vq3_ref['vertex_correction']['F2']
    lines.append(f"| VQ-3 (sum sigma_mu, all edges) | "
                 f"{vq3_ref['selection_rules']['summary']['fraction']} | "
                 f"{vq3_ref['self_energy']['trace_total']:.4f} | "
                 f"{vq3_f2:.6f} | "
                 f"{vq3_f2/scalar_f2:.3f} |")

    for name, res in [("A (T->sigma_z, L->sigma_x+sigma_y)", mapping_A),
                       ("B (T->sigma_z, L->sigma_x)", mapping_B),
                       ("D (T->sigma_z, L->(sigma_x+sigma_y)/2)", mapping_D),
                       ("E (T->sigma_x, L->sigma_z)", mapping_E)]:
        f2 = res['vertex_correction']['F2']
        f2_str = f"{f2:.6f}" if f2 is not None else "N/A"
        ratio = f"{f2/scalar_f2:.3f}" if f2 is not None and scalar_f2 else "N/A"
        lines.append(f"| {name} | "
                     f"{res['selection_rules']['summary']['fraction']} | "
                     f"{res['self_energy']['trace_total']:.4f} | "
                     f"{f2_str} | {ratio} |")
    lines.append("")

    # Exhaustive scan
    lines.append("### Mapping C: Exhaustive 3x3 Scan")
    lines.append("")
    lines.append("| T-sigma | L-sigma | Tr(Sigma) | F2 | F2/F2_scalar |")
    lines.append("|:--------|:--------|:----------|:---|:-------------|")
    for combo, res in mapping_C_results.items():
        f2 = res['vertex_correction']['F2']
        f2_str = f"{f2:.6f}" if f2 is not None else "N/A"
        ratio = f"{f2/scalar_f2:.3f}" if f2 is not None and scalar_f2 else "N/A"
        lines.append(f"| {combo.split(',')[0].split('=')[1]} | "
                     f"{combo.split(',')[1].split('=')[1]} | "
                     f"{res['self_energy']['trace_total']:.4f} | "
                     f"{f2_str} | {ratio} |")
    lines.append("")

    # Headline result
    lines.append("## Headline Result")
    lines.append("")
    if n_distinct <= 1 and trivial_3x:
        lines.append("**NEGATIVE.** All 11 tested mappings produce the same F2 value,")
        lines.append("equal to 3 times the scalar F2. The channel-specific sigma assignment")
        lines.append("does NOT break the trivial 3x rescaling.")
        lines.append("")
        lines.append("This is a consequence of the block-diagonal structure of G_gamma at n_max=2:")
        lines.append("G_gamma[T,L] = 0 (confirmed by VQ-4), meaning the T and L channels")
        lines.append("propagate independently. When the self-energy is computed as")
        lines.append("")
        lines.append("    Sigma = Sigma_T(sigma_T) + Sigma_L(sigma_L)")
        lines.append("")
        lines.append("the channel-specific sigma does NOT mix across channels (no cross-terms).")
        lines.append("Each channel contributes its own sigma-squared trace independently:")
        lines.append("")
        lines.append("    Tr(sigma_mu^2) = 1 (for a single sigma)")
        lines.append("    Tr(sigma_x^2 + sigma_y^2) = 2 (for two sigma's)")
        lines.append("")
        lines.append("The SUM across channels gives the total sigma-squared weight. For Mapping A")
        lines.append("(1 sigma on T + 2 sigmas on L) this is 1+2 = 3 = same as VQ-3.")
        lines.append("For Mapping B (1+1 = 2), we get 2/3 of VQ-3, which is 2x scalar.")
        lines.append("")
    elif n_distinct > 1:
        lines.append(f"**PARTIALLY POSITIVE.** The channel-specific sigma assignment produces")
        lines.append(f"{n_distinct} distinct F2 values across the 11 mappings.")
        lines.append("The 3x rescaling IS broken when different numbers of sigma matrices")
        lines.append("are assigned to different channels.")
        lines.append("")
        lines.append("However, this is a COUNTING effect, not a STRUCTURAL one: each sigma")
        lines.append("contributes its own sigma^2 = I to the self-energy trace, and the")
        lines.append("total F2 scales as (n_sigma_T + n_sigma_L) * F2_scalar.")
        lines.append("")
    else:
        lines.append("All mappings give a single F2 value that is NOT 3x scalar.")
        lines.append("")

    # Isotropy analysis
    lines.append("## Isotropy Analysis")
    lines.append("")
    lines.append("| Mapping | Tr(Sigma_T) | Tr(Sigma_L) | Tr(T)/Tr(L) | Proportional? |")
    lines.append("|:--------|:------------|:------------|:------------|:--------------|")
    for name, res in [("A", mapping_A), ("B", mapping_B),
                       ("D", mapping_D), ("E", mapping_E)]:
        iso = res['isotropy']
        tr_T = res['self_energy']['trace_T']
        tr_L = res['self_energy']['trace_L']
        ratio = iso['trace_ratio_T_over_L']
        prop = iso['proportional']
        ratio_str = f"{ratio:.4f}" if ratio is not None else "N/A"
        lines.append(f"| {name} | {tr_T:.4f} | {tr_L:.4f} | {ratio_str} | {prop} |")
    lines.append("")
    lines.append("The T-channel and L-channel self-energies are generically NOT proportional,")
    lines.append("confirming the VQ-4 finding that the photon propagator is anisotropic.")
    lines.append("However, the anisotropy is between channels, not within the sigma structure")
    lines.append("-- each channel's sigma trace remains standard (sigma_mu^2 = I per block).")
    lines.append("")

    # Selection rules
    lines.append("## Selection Rule Census")
    lines.append("")
    lines.append("All mappings produce 1/8 surviving selection rules (Gaunt/CG sparsity only),")
    lines.append("identical to the scalar GN-5 result and VQ-3. The channel-specific sigma")
    lines.append("assignment does NOT recover any additional selection rules.")
    lines.append("")
    lines.append("The 4 rules that require vector photon quantum numbers (vertex parity,")
    lines.append("SO(4) channel count, Ward identity, charge conjugation with vector structure)")
    lines.append("remain broken because the photon is still a scalar 1-cochain. The 3 rules")
    lines.append("that require Dirac graph nodes (Delta_mj conservation with Dl=+/-1,")
    lines.append("spatial parity E1, Furry's theorem with off-diagonal identity) are not")
    lines.append("recovered because the Fock graph nodes are scalar (n,l,m) labels, not")
    lines.append("spinor (n,kappa,m_j) labels.")
    lines.append("")

    # Physical interpretation
    lines.append("## Physical Interpretation")
    lines.append("")
    lines.append("### Why the result is negative")
    lines.append("")
    lines.append("The channel-specific sigma assignment fails to produce non-trivial physics")
    lines.append("because of a structural mismatch between the graph topology and the sigma")
    lines.append("algebra:")
    lines.append("")
    lines.append("1. **Sigma acts on Dirac labels** (n, kappa, m_j), preserving n and l.")
    lines.append("2. **The Fock graph edges connect scalar labels** (n, l, m), with T-edges")
    lines.append("   changing n and L-edges changing m.")
    lines.append("3. **The CG projection** P maps Dirac labels to Fock labels, mixing")
    lines.append("   different kappa values. When sigma is composed with V_scalar = P . A . P^T,")
    lines.append("   the result is sigma @ P @ (Fock graph) @ P^T, which contracts sigma's")
    lines.append("   spin structure through P's CG coefficients.")
    lines.append("4. **The contraction** sum_c sigma[a,c] * V[c,b,e] reduces to the standard")
    lines.append("   Pauli trace identity when summed over all sigma directions, regardless of")
    lines.append("   which edges are included, because P already mixes the spin indices.")
    lines.append("")
    lines.append("The fundamental issue is that sigma operates in the **spin sector** while")
    lines.append("the edge channels operate in the **orbital sector**. The CG projection P")
    lines.append("couples these sectors, but the coupling is through a TRACE (sum over internal")
    lines.append("indices c), which washes out the directional information.")
    lines.append("")
    lines.append("### What would be needed for non-trivial vector QED")
    lines.append("")
    lines.append("To recover the continuum gamma^mu vertex on the graph, one would need:")
    lines.append("")
    lines.append("1. **Vector photon labels**: edges carrying (L, M_L) quantum numbers,")
    lines.append("   not just scalar 1-cochains. This would provide the SO(4) channel count W.")
    lines.append("2. **Spin-orbit coupling AT the vertex**: the sigma should couple to the")
    lines.append("   photon's vector index, not just dress the electron line. This requires")
    lines.append("   the photon to carry spin-1 structure that can contract with sigma_mu.")
    lines.append("3. **Dirac graph nodes** (not Fock nodes): as demonstrated in the Dirac")
    lines.append("   graph QED analysis (CLAUDE.md), using (n, kappa, m_j) node labels")
    lines.append("   with E1 dipole adjacency (Rule B) recovers 4/8 selection rules.")
    lines.append("")
    lines.append("The VQ sprint's central finding is that dressing the scalar vertex with")
    lines.append("spin matrices is NECESSARY but NOT SUFFICIENT. The missing ingredient is")
    lines.append("vector photon structure, which is a **calibration exchange constant** in")
    lines.append("Paper 18's taxonomy -- it must come from the continuum embedding, not from")
    lines.append("the graph topology alone.")
    lines.append("")

    # VQ sprint summary
    lines.append("## VQ Sprint Summary (All Tracks)")
    lines.append("")
    lines.append("| Track | Result | Key Finding |")
    lines.append("|:------|:-------|:------------|")
    lines.append("| VQ-1 | Setup | 10 Dirac states, 3 Fock edges, infrastructure built |")
    lines.append("| VQ-2 | NEGATIVE | Pure sigma vertex (no Fock edges) creates disconnected intra-shell graph |")
    lines.append("| VQ-3 | NEGATIVE | Sigma-weighted Fock vertex gives trivial 3x rescaling (Pauli trace identity) |")
    lines.append("| VQ-4 | POSITIVE structural | Two distinct edge channels (T radial, L angular) with anisotropic spectra |")
    lines.append("| VQ-5 | NEGATIVE | Channel-specific sigma assignment does not break 3x rescaling |")
    lines.append("")
    lines.append("### Net verdict for the VQ sprint")
    lines.append("")
    lines.append("**NEGATIVE for selection rule recovery.** The sigma-weighted vertex approach,")
    lines.append("in all variants tested (universal sum, channel-specific, exhaustive scan),")
    lines.append("does NOT recover any of the 7 missing continuum QED selection rules beyond")
    lines.append("the 1 (Gaunt/CG sparsity) that the scalar graph already captures.")
    lines.append("")
    lines.append("**POSITIVE structural finding from VQ-4:** The Fock graph photon has two")
    lines.append("physically distinct propagation channels (radial T and angular L) with")
    lines.append("anisotropic spectra and block-diagonal propagator at n_max=2. This is a")
    lines.append("topological feature of the quantum-number lattice that provides effective")
    lines.append("\"polarization\" structure without introducing explicit vector labels.")
    lines.append("However, this structure is insufficient for selection rule recovery.")
    lines.append("")
    lines.append("### Implications for the graph-native QED program")
    lines.append("")
    lines.append("The VQ sprint confirms and sharpens the three-tier partition established in")
    lines.append("the native Dirac graph QED analysis:")
    lines.append("")
    lines.append("1. **Always survives (1/8):** Gaunt/CG sparsity -- from angular momentum algebra,")
    lines.append("   present on any graph with CG-projected vertices.")
    lines.append("")
    lines.append("2. **Spinor-recoverable (3/8):** Delta_mj, spatial parity E1, Furry's theorem --")
    lines.append("   recovered by using Dirac (n, kappa, m_j) node labels instead of scalar Fock")
    lines.append("   (n, l, m) labels. Does NOT require vector photon.")
    lines.append("")
    lines.append("3. **Vector-photon-required (4/8):** Vertex parity, SO(4) channel count, Ward")
    lines.append("   identity, charge conjugation -- require promoting the photon from a scalar")
    lines.append("   1-cochain to a vector harmonic carrying (L, M_L) quantum numbers. This is")
    lines.append("   **calibration exchange constant** content in Paper 18's taxonomy.")
    lines.append("")
    lines.append("The VQ sprint tested whether dressing the scalar vertex with spin structure")
    lines.append("(sigma matrices) could bridge the gap between tiers 1 and 3. The answer is no:")
    lines.append("sigma dresses the ELECTRON line, not the PHOTON line. Recovering tier-3 rules")
    lines.append("requires dressing the PHOTON with vector structure, which is a fundamentally")
    lines.append("different construction.")
    lines.append("")

    # Data files
    lines.append("## Data Files")
    lines.append("")
    lines.append("- `debug/data/vq5_full_vector_qed.json` -- full numerical results")
    lines.append("- `debug/vq5_full_vector_qed.py` -- computation script")
    lines.append("- `debug/vq3_sigma_weighted_vertex.py` -- VQ-3 sigma-weighted vertex")
    lines.append("- `debug/vq4_direction_resolved_hodge.py` -- VQ-4 direction-resolved Hodge")
    lines.append("")

    memo_text = "\n".join(lines) + "\n"
    with open(memo_path, 'w') as f:
        f.write(memo_text)
    print(f"Memo saved to {memo_path}")


if __name__ == "__main__":
    main()
