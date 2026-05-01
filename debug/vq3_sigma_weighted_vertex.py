"""
VQ-3: Sigma-weighted Fock graph vertex for graph-native QED on S^3.
====================================================================

Tests whether inserting the Pauli spin matrix sigma_mu at the vertex
-- keeping the SAME Fock graph edges and SAME photon propagator from
existing graph-native QED -- recovers any of the 4 missing continuum
QED selection rules.

The construction combines sigma (spin structure, intra-shell) with
the scalar CG vertex (inter-shell Fock graph connectivity):

    V_mu[a,b,e] = sum_c  sigma_mu[a,c] * V_scalar[c,b,e]

This fixes VQ-2's problem where sigma alone created a disconnected
intra-shell-only graph: here the Fock graph edges (inter-shell T+/T-
and intra-shell L+/L-) are preserved, but each vertex gets dressed
with a spin matrix.

Self-energy:
    Sigma_vec[a,b] = sum_mu sum_{e,e'} G_gamma[e,e']
                     * (sigma_mu @ V_e)[a,:] @ (sigma_mu @ V_{e'})^dagger[:,b]

Vertex correction:
    Lambda_vec[a,b] = sum_mu sum_{e,e'} G_gamma[e,e']
                      * (sigma_mu @ V_e)[a,:] @ G_e @ (sigma_mu @ V_{e'})^dagger[:,b]

Anomalous magnetic moment:
    F_2 = Tr(Lambda_vec) / Tr(V_bare_vec @ G_e)

where V_bare_vec = sum_mu sum_e (sigma_mu @ V_e).

This script does NOT modify any production code.
"""

import numpy as np
import json
import sys
import os
from fractions import Fraction
from collections import defaultdict

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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
from sympy.physics.wigner import clebsch_gordan
from sympy import Rational


# ===========================================================================
# Helper: CG coefficient as float
# ===========================================================================

def cg_float(j1, m1, j2, m2, j, m):
    """CG coefficient <j1,m1; j2,m2 | j,m> as float."""
    return float(clebsch_gordan(j1, j2, j, m1, m2, m))


# ===========================================================================
# Step 1: Build sigma matrices in the DiracLabel basis (complex)
# ===========================================================================

def build_sigma_matrices_complex(labels, label_index):
    """Build sigma_x, sigma_y, sigma_z as complex numpy arrays.

    sigma_x: real symmetric
    sigma_y: purely imaginary Hermitian
    sigma_z: real symmetric

    Sigma matrices act on the SPIN part, preserving n and l.
    In the coupled (j, m_j) basis, they can change kappa (j changes
    while l stays fixed).
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

    # Standard convention: sigma_+ = (sigma_x + i*sigma_y)/2
    # so sigma_x = sigma_+ + sigma_-, sigma_y = -i*(sigma_+ - sigma_-)
    sigma_x = sigma_plus + sigma_minus
    sigma_y = -1j * (sigma_plus - sigma_minus)

    return sigma_x, sigma_y, sigma_z


def verify_sigma_algebra(sigma_x, sigma_y, sigma_z, labels):
    """Quick verification of Pauli algebra."""
    blocks = defaultdict(list)
    for i, lab in enumerate(labels):
        blocks[(lab.n_fock, lab.l)].append(i)

    results = {}

    # sigma^2 = 3I within each block
    sigma_sq = sigma_x @ sigma_x + sigma_y @ sigma_y + sigma_z @ sigma_z
    all_pass = True
    for key, indices in blocks.items():
        idx = np.array(indices)
        block = sigma_sq[np.ix_(idx, idx)]
        expected = 3.0 * np.eye(len(idx))
        if np.max(np.abs(block - expected)) > 1e-10:
            all_pass = False
    results["sigma_squared_3I"] = all_pass

    # [sigma_x, sigma_y] = 2i*sigma_z
    comm = sigma_x @ sigma_y - sigma_y @ sigma_x
    comm_pass = True
    for key, indices in blocks.items():
        idx = np.array(indices)
        err = np.max(np.abs(comm[np.ix_(idx, idx)] - 2j * sigma_z[np.ix_(idx, idx)]))
        if err > 1e-10:
            comm_pass = False
    results["commutation_ok"] = comm_pass

    # Hermiticity
    results["sigma_x_hermitian"] = bool(np.allclose(sigma_x, sigma_x.conj().T, atol=1e-12))
    results["sigma_y_hermitian"] = bool(np.allclose(sigma_y, sigma_y.conj().T, atol=1e-12))
    results["sigma_z_hermitian"] = bool(np.allclose(sigma_z, sigma_z.conj().T, atol=1e-12))

    # Real/imaginary structure
    results["sigma_x_real"] = bool(np.allclose(sigma_x.imag, 0, atol=1e-12))
    results["sigma_y_purely_imaginary"] = bool(np.allclose(sigma_y.real, 0, atol=1e-12))
    results["sigma_z_real"] = bool(np.allclose(sigma_z.imag, 0, atol=1e-12))

    return results


# ===========================================================================
# Step 2: Build scalar vertex matrices from production code (numpy)
# ===========================================================================

def build_scalar_vertex_numpy(n_max):
    """Build the scalar CG vertex matrices V_e from the production code,
    converted to numpy arrays.

    Returns
    -------
    V_e_list : list of np.ndarray, shape (N_dirac, N_dirac)
        One matrix per Fock edge.
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
    # vertex_tensor_to_matrices returns sympy Matrices; convert to numpy
    V_mats_sympy = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    V_e_list = []
    for V_sym in V_mats_sympy:
        V_np = np.array(V_sym.tolist(), dtype=float)
        V_e_list.append(V_np)

    return V_e_list, dirac_labels, fock_data, E_fock, N_dirac


# ===========================================================================
# Step 3: Build sigma-weighted vertex
# ===========================================================================

def build_sigma_weighted_vertex(sigma_list, V_e_list, N_dirac, E_fock):
    """Build sigma-weighted vertex matrices.

    V_mu_e[a,b] = sum_c sigma_mu[a,c] * V_scalar_e[c,b]
               = (sigma_mu @ V_e)[a,b]

    Parameters
    ----------
    sigma_list : list of 3 arrays [sigma_x, sigma_y, sigma_z]
    V_e_list : list of E numpy arrays, each (N_dirac, N_dirac)
    N_dirac : int
    E_fock : int

    Returns
    -------
    V_sigma : list of 3 lists, each containing E numpy arrays (N_dirac, N_dirac)
        V_sigma[mu][e] = sigma_mu @ V_e
    """
    V_sigma = []
    for sigma_mu in sigma_list:
        V_mu_list = []
        for V_e in V_e_list:
            V_mu_e = sigma_mu @ V_e
            V_mu_list.append(V_mu_e)
        V_sigma.append(V_mu_list)
    return V_sigma


# ===========================================================================
# Step 4: Compute self-energy and vertex correction
# ===========================================================================

def compute_sigma_weighted_self_energy(V_sigma, G_gamma, N_dirac, E_fock):
    """Compute the sigma-weighted self-energy.

    Sigma_vec[a,b] = sum_mu sum_{e,e'} G_gamma[e,e']
                     * (sigma_mu @ V_e)[a,:] @ (sigma_mu @ V_{e'})^dagger[:,b]

    Uses Hermitian contraction (V^dagger) since sigma_mu are Hermitian.
    """
    Sigma = np.zeros((N_dirac, N_dirac), dtype=complex)
    for mu_idx in range(3):
        for e1 in range(E_fock):
            for e2 in range(E_fock):
                g_ee = G_gamma[e1, e2]
                if abs(g_ee) < 1e-15:
                    continue
                V1 = V_sigma[mu_idx][e1]
                V2 = V_sigma[mu_idx][e2]
                Sigma += g_ee * (V1 @ V2.conj().T)
    return Sigma


def compute_sigma_weighted_vertex_correction(V_sigma, G_gamma, G_e, N_dirac, E_fock):
    """Compute the sigma-weighted vertex correction.

    Lambda_vec[a,b] = sum_mu sum_{e,e'} G_gamma[e,e']
                      * (sigma_mu @ V_e)[a,:] @ G_e @ (sigma_mu @ V_{e'})^dagger[:,b]
    """
    Lambda_total = np.zeros((N_dirac, N_dirac), dtype=complex)
    for mu_idx in range(3):
        for e1 in range(E_fock):
            for e2 in range(E_fock):
                g_ee = G_gamma[e1, e2]
                if abs(g_ee) < 1e-15:
                    continue
                V1 = V_sigma[mu_idx][e1]
                V2 = V_sigma[mu_idx][e2]
                Lambda_total += g_ee * (V1 @ G_e @ V2.conj().T)
    return Lambda_total


def extract_f2(Lambda_total, V_sigma, G_e, N_dirac, E_fock):
    """Extract the anomalous magnetic moment F_2.

    F_2 = Tr(Lambda_total) / Tr(V_bare_vec @ G_e)

    where V_bare_vec = sum_mu sum_e (sigma_mu @ V_e).
    """
    V_bare = np.zeros((N_dirac, N_dirac), dtype=complex)
    for mu_idx in range(3):
        for e_idx in range(E_fock):
            V_bare += V_sigma[mu_idx][e_idx]

    denominator = np.trace(V_bare @ G_e)
    numerator = np.trace(Lambda_total)

    if abs(denominator) < 1e-15:
        return None, V_bare, numerator, denominator

    F2 = numerator / denominator
    return F2, V_bare, numerator, denominator


# ===========================================================================
# Step 5: Selection rule census
# ===========================================================================

def selection_rule_census(V_sigma, V_e_scalar, sigma_list, labels, fock_data,
                          Sigma, G_e, G_gamma, N_dirac, E_fock):
    """Test 8 continuum QED selection rules on the sigma-weighted vertex."""
    N = N_dirac
    E = E_fock
    results = {}

    # Collect all nonzero sigma-weighted vertex entries for analysis
    all_vertex_entries = []
    for mu_idx in range(3):
        for e_idx in range(E):
            V = V_sigma[mu_idx][e_idx]
            for a in range(N):
                for b in range(N):
                    if abs(V[a, b]) > 1e-12:
                        all_vertex_entries.append({
                            "mu": mu_idx,
                            "edge": e_idx,
                            "a": a, "b": b,
                            "value_abs": float(abs(V[a, b])),
                            "lab_a": labels[a],
                            "lab_b": labels[b],
                        })

    # -----------------------------------------------------------------------
    # 1. Delta m_j conservation
    # -----------------------------------------------------------------------
    dm_j_values = set()
    for entry in all_vertex_entries:
        dm = entry["lab_a"].two_m_j - entry["lab_b"].two_m_j
        dm_j_values.add(dm)
    survives_dmj = dm_j_values.issubset({-2, -1, 0, 1, 2})
    results["delta_mj_conservation"] = {
        "delta_2mj_values": sorted(dm_j_values),
        "verdict": ("|Delta(m_j)| <= 1 only" if survives_dmj
                    else f"BROKEN: values {sorted(dm_j_values)}"),
        "survives": bool(survives_dmj),
    }

    # -----------------------------------------------------------------------
    # 2. Spatial parity E1 (Delta l odd)
    # -----------------------------------------------------------------------
    # The scalar vertex V_e connects Fock nodes in adjacent shells via T+/T-
    # (Delta n = +/-1, Delta l = 0) and within-shell L+/L- (Delta n = 0, Delta l = 0).
    # Sigma preserves l. So the sigma-weighted vertex inherits the scalar vertex's
    # l structure. Check l_a + l_b parity for all nonzero entries.
    parity_check = {"odd": 0, "even": 0}
    for entry in all_vertex_entries:
        l_sum = entry["lab_a"].l + entry["lab_b"].l
        if l_sum % 2 == 1:
            parity_check["odd"] += 1
        else:
            parity_check["even"] += 1

    # E1 transition requires l_a + l_b odd (Delta l = +/-1).
    # The Fock graph scalar vertex has Delta l = 0 (same l on both ends),
    # so l_a + l_b is always even. Sigma preserves l. So combined vertex
    # also has l_a + l_b even -> E1 parity NOT satisfied.
    results["spatial_parity_E1"] = {
        "l_sum_parity": parity_check,
        "verdict": ("ALL even (no E1 transitions)" if parity_check["odd"] == 0
                    else f"MIXED: {parity_check}"),
        "survives": bool(parity_check["odd"] > 0),
        "note": ("Fock graph has Dl=0; sigma preserves l; combined has Dl=0 "
                 "on both vertex legs => l_a+l_b always even"),
    }

    # -----------------------------------------------------------------------
    # 3. Gaunt/CG sparsity
    # -----------------------------------------------------------------------
    total_possible = N * N * E * 3  # 3 polarizations, all pairs, all edges
    nonzero_count = len(all_vertex_entries)
    density = nonzero_count / total_possible if total_possible > 0 else 0
    results["gaunt_cg_sparsity"] = {
        "total_possible": total_possible,
        "nonzero": nonzero_count,
        "density_percent": float(100 * density),
        "survives": True,  # CG triangle inequality always produces sparsity
    }

    # -----------------------------------------------------------------------
    # 4. Vertex parity (n_1 + n_2 + q odd in continuum)
    # -----------------------------------------------------------------------
    # Check the n_fock values on the vertex legs
    dn_values = defaultdict(int)
    n_sum_parity = {"even": 0, "odd": 0}
    for entry in all_vertex_entries:
        dn = abs(entry["lab_a"].n_fock - entry["lab_b"].n_fock)
        dn_values[dn] += 1
        s = entry["lab_a"].n_fock + entry["lab_b"].n_fock
        if s % 2 == 0:
            n_sum_parity["even"] += 1
        else:
            n_sum_parity["odd"] += 1

    # For vertex parity: the scalar vertex has Delta(n)=0 and Delta(n)=+/-1
    # (from L+/L- and T+/T- edges). Sigma preserves n. So the combined
    # vertex has the same n1+n2 distribution as the scalar vertex.
    # The scalar vertex has BOTH parities (n+n = even for same-shell,
    # n+(n+1) = odd for adjacent-shell). So vertex parity is NOT enforced.
    results["vertex_parity"] = {
        "n_sum_parity": n_sum_parity,
        "delta_n_values": dict(dn_values),
        "verdict": ("MIXED (both even and odd n-sums present)"
                    if n_sum_parity["even"] > 0 and n_sum_parity["odd"] > 0
                    else "RESTRICTED"),
        "survives": bool(n_sum_parity["even"] == 0 or n_sum_parity["odd"] == 0),
    }

    # -----------------------------------------------------------------------
    # 5. SO(4) channel count
    # -----------------------------------------------------------------------
    # The scalar vertex inherits SO(4) channel structure from the Fock graph.
    # Sigma dresses the vertex but doesn't change the photon edge structure.
    # SO(4) channel count W requires vector photon quantum numbers —
    # the sigma-weighted vertex still uses scalar photon (1-cochain).
    results["so4_channel_count"] = {
        "verdict": ("N/A -- photon is still a scalar 1-cochain; "
                    "SO(4) channel count requires vector photon harmonics"),
        "survives": False,
    }

    # -----------------------------------------------------------------------
    # 6. Charge conjugation (C)
    # -----------------------------------------------------------------------
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

    if np.allclose(Sigma.imag, 0, atol=1e-10):
        Sigma_work = Sigma.real
    else:
        Sigma_work = Sigma
    CSC = C_mat @ Sigma_work @ C_mat
    c_symmetric = bool(np.allclose(CSC, Sigma_work, atol=1e-10))
    results["charge_conjugation"] = {
        "C_symmetric": c_symmetric,
        "verdict": "C-symmetric" if c_symmetric else "C broken",
        "survives": c_symmetric,
    }

    # -----------------------------------------------------------------------
    # 7. Furry's theorem (tadpole = 0)
    # -----------------------------------------------------------------------
    tadpole = {}
    mu_names = ['x', 'y', 'z']
    for mu_idx, mu_name in enumerate(mu_names):
        tad = np.zeros((N, N), dtype=complex)
        for e_idx in range(E):
            tad += V_sigma[mu_idx][e_idx]
        tadpole_trace = np.trace(tad)
        tadpole_norm = np.max(np.abs(tad))
        tadpole[mu_name] = {
            "trace": float(abs(tadpole_trace)),
            "max_abs": float(tadpole_norm),
            "is_zero": bool(tadpole_norm < 1e-10),
        }
    all_zero = all(tadpole[m]["is_zero"] for m in mu_names)
    results["furry_theorem"] = {
        "per_polarization": tadpole,
        "all_tadpoles_zero": all_zero,
        "verdict": ("Furry satisfied" if all_zero else "Furry VIOLATED"),
        "survives": all_zero,
    }

    # -----------------------------------------------------------------------
    # 8. Ward identity
    # -----------------------------------------------------------------------
    D = np.diag([
        (1.0 if lab.kappa < 0 else -1.0) * (lab.n_fock + 0.5)
        for lab in labels
    ])
    if np.allclose(Sigma.imag, 0, atol=1e-10):
        comm = D @ Sigma.real - Sigma.real @ D
    else:
        comm = D @ Sigma - Sigma @ D
    comm_norm = np.max(np.abs(comm))
    results["ward_identity"] = {
        "commutator_D_Sigma_norm": float(comm_norm),
        "commutes": bool(comm_norm < 1e-10),
        "verdict": "[D, Sigma] = 0" if comm_norm < 1e-10 else "[D, Sigma] != 0",
        "survives": False,  # Ward identity needs full vertex correction comparison
        "note": "Ward identity is a vertex-propagator relation, not just [D,Sigma]",
    }

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    rule_names = [
        "delta_mj_conservation",
        "spatial_parity_E1",
        "gaunt_cg_sparsity",
        "vertex_parity",
        "so4_channel_count",
        "charge_conjugation",
        "furry_theorem",
        "ward_identity",
    ]
    survived = sum(1 for r in rule_names if results[r].get("survives", False))
    results["summary"] = {
        "total_rules": 8,
        "survived": survived,
        "fraction": f"{survived}/8",
    }

    return results


# ===========================================================================
# Step 6: Algebraic / number field analysis
# ===========================================================================

def analyze_number_field(Sigma, Lambda_total, F2, labels):
    """Determine the number field of the sigma-weighted QED quantities."""
    results = {}

    # Self-energy
    is_real = bool(np.allclose(Sigma.imag, 0, atol=1e-10))
    results["Sigma_is_real"] = is_real

    if is_real:
        S = Sigma.real
    else:
        S = Sigma

    # Try to identify nonzero entries as simple fractions or sqrt(rational)
    entries_info = []
    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            v = S[i, j]
            if isinstance(v, complex):
                v = v.real if abs(v.imag) < 1e-12 else v
            if abs(v) > 1e-12:
                # Check rational
                best_frac = None
                best_err = 1e-6
                for q in range(1, 300):
                    p = round(float(v.real if isinstance(v, complex) else v) * q)
                    err = abs(float(v.real if isinstance(v, complex) else v) - p / q)
                    if err < best_err:
                        best_err = err
                        best_frac = (p, q)
                entries_info.append({
                    "i": i, "j": j,
                    "value": float(v.real if isinstance(v, complex) else v),
                    "fraction": f"{best_frac[0]}/{best_frac[1]}" if best_frac else None,
                    "residual": float(best_err) if best_frac else None,
                })

    all_rational = all((e.get("residual") or 1.0) < 1e-8 for e in entries_info)
    results["Sigma_appears_rational"] = all_rational

    # Check for sqrt content: if not rational, check v^2 is rational
    if not all_rational:
        sqrt_candidates = []
        for e in entries_info:
            if (e.get("residual") or 0) > 1e-8:
                v2 = e["value"] ** 2
                best_frac2 = None
                best_err2 = 1e-6
                for q in range(1, 300):
                    p = round(v2 * q)
                    err = abs(v2 - p / q)
                    if err < best_err2:
                        best_err2 = err
                        best_frac2 = (p, q)
                if best_err2 < 1e-8:
                    sqrt_candidates.append({
                        "value": e["value"],
                        "v_squared_fraction": f"{best_frac2[0]}/{best_frac2[1]}",
                    })
        results["sqrt_rational_candidates"] = sqrt_candidates[:10]

    results["nonzero_entries_sample"] = entries_info[:15]

    # Lambda number field
    results["Lambda_is_real"] = bool(np.allclose(Lambda_total.imag, 0, atol=1e-10))

    # F2 number field
    if F2 is not None:
        f2_val = float(F2.real if isinstance(F2, complex) else F2)
        results["F2_value"] = f2_val
        # Check if rational
        best_frac_f2 = None
        best_err_f2 = 1e-6
        for q in range(1, 500):
            p = round(f2_val * q)
            err = abs(f2_val - p / q)
            if err < best_err_f2:
                best_err_f2 = err
                best_frac_f2 = (p, q)
        results["F2_fraction"] = f"{best_frac_f2[0]}/{best_frac_f2[1]}" if best_frac_f2 else None
        results["F2_fraction_residual"] = float(best_err_f2) if best_frac_f2 else None

        # Check sqrt(rational)
        f2_sq = f2_val ** 2
        best_frac_f2sq = None
        best_err_f2sq = 1e-6
        for q in range(1, 500):
            p = round(f2_sq * q)
            err = abs(f2_sq - p / q)
            if err < best_err_f2sq:
                best_err_f2sq = err
                best_frac_f2sq = (p, q)
        results["F2_squared_fraction"] = f"{best_frac_f2sq[0]}/{best_frac_f2sq[1]}" if best_frac_f2sq else None
        results["F2_squared_residual"] = float(best_err_f2sq) if best_frac_f2sq else None

    # Comparison to scalar QED number field
    results["comparison"] = {
        "scalar_fock_QED_F2": "5*sqrt(2)/3 in Q[sqrt(2),sqrt(3),sqrt(6)]",
        "dirac_graph_QED": "Q[sqrt(2),sqrt(17),sqrt(41),sqrt(881)]",
        "sigma_weighted_QED": ("Q (rationals)" if all_rational
                                else "algebraic -- see sqrt candidates"),
    }

    return results


# ===========================================================================
# Step 7: Compare to scalar and VQ-2
# ===========================================================================

def compare_to_scalar(Sigma_sigma, V_e_scalar, G_gamma, labels, N_dirac, E_fock):
    """Compare sigma-weighted self-energy to the scalar self-energy."""
    # Compute scalar self-energy for comparison
    Sigma_scalar = np.zeros((N_dirac, N_dirac), dtype=float)
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_ee = G_gamma[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            Sigma_scalar += g_ee * (V_e_scalar[e1] @ V_e_scalar[e2].T)

    results = {}
    results["scalar_trace"] = float(np.trace(Sigma_scalar))
    results["sigma_weighted_trace"] = float(np.trace(Sigma_sigma.real))
    results["trace_ratio"] = (float(np.trace(Sigma_sigma.real)) /
                               float(np.trace(Sigma_scalar))
                               if abs(np.trace(Sigma_scalar)) > 1e-15 else None)

    results["scalar_eigenvalues"] = sorted(np.linalg.eigvalsh(Sigma_scalar).tolist())
    if np.allclose(Sigma_sigma.imag, 0, atol=1e-10):
        results["sigma_eigenvalues"] = sorted(np.linalg.eigvalsh(Sigma_sigma.real).tolist())
    else:
        results["sigma_eigenvalues"] = sorted(np.linalg.eigvals(Sigma_sigma).real.tolist())

    # Ground state block comparison
    gs_idx = [i for i, lab in enumerate(labels) if lab.n_fock == 1 and lab.kappa == -1]
    if len(gs_idx) > 0:
        gs_scalar = Sigma_scalar[np.ix_(gs_idx, gs_idx)]
        gs_sigma = Sigma_sigma[np.ix_(gs_idx, gs_idx)]
        results["gs_block_scalar"] = gs_scalar.tolist()
        results["gs_block_sigma"] = [[float(gs_sigma[i, j].real)
                                       for j in range(gs_sigma.shape[1])]
                                      for i in range(gs_sigma.shape[0])]

    return results


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 70)
    print("VQ-3: Sigma-weighted Fock graph vertex for graph-native QED on S^3")
    print("=" * 70)

    n_max = 2

    # ------------------------------------------------------------------
    # Step 1: Build sigma matrices
    # ------------------------------------------------------------------
    print("\n--- Step 1: Sigma matrices ---")
    labels = list(iter_dirac_labels(n_max))
    label_index = {lab: i for i, lab in enumerate(labels)}
    N = len(labels)
    print(f"n_max = {n_max}, N_dirac = {N}")
    for i, lab in enumerate(labels):
        print(f"  [{i}] n={lab.n_fock}, kappa={lab.kappa}, l={lab.l}, "
              f"j={float(lab.j):.1f}, m_j={float(lab.m_j):.1f}")

    sigma_x, sigma_y, sigma_z = build_sigma_matrices_complex(labels, label_index)
    sigma_results = verify_sigma_algebra(sigma_x, sigma_y, sigma_z, labels)
    print(f"  sigma^2 = 3I: {sigma_results['sigma_squared_3I']}")
    print(f"  commutation: {sigma_results['commutation_ok']}")
    print(f"  sigma_x real: {sigma_results['sigma_x_real']}")
    print(f"  sigma_y imag: {sigma_results['sigma_y_purely_imaginary']}")
    print(f"  sigma_z real: {sigma_results['sigma_z_real']}")

    # ------------------------------------------------------------------
    # Step 2: Build scalar vertex from production code
    # ------------------------------------------------------------------
    print("\n--- Step 2: Scalar vertex from production code ---")
    V_e_scalar, dirac_labels, fock_data, E_fock, N_dirac = build_scalar_vertex_numpy(n_max)
    print(f"  E_fock = {E_fock} edges")
    print(f"  N_dirac = {N_dirac} states")

    # Quick check: how many nonzero entries in scalar vertex
    scalar_nonzero = sum(np.count_nonzero(np.abs(V) > 1e-12) for V in V_e_scalar)
    print(f"  Scalar vertex nonzero entries (total across all edges): {scalar_nonzero}")

    # Edge list from fock_data
    print(f"  Fock graph edges: {fock_data.edges}")

    # ------------------------------------------------------------------
    # Step 3: Build sigma-weighted vertex
    # ------------------------------------------------------------------
    print("\n--- Step 3: Sigma-weighted vertex ---")
    sigma_list = [sigma_x, sigma_y, sigma_z]
    V_sigma = build_sigma_weighted_vertex(sigma_list, V_e_scalar, N_dirac, E_fock)

    # Count nonzero entries
    sigma_nonzero = 0
    for mu_idx in range(3):
        for e_idx in range(E_fock):
            sigma_nonzero += np.count_nonzero(np.abs(V_sigma[mu_idx][e_idx]) > 1e-12)
    print(f"  Sigma-weighted vertex nonzero entries (total): {sigma_nonzero}")
    print(f"  Ratio to scalar: {sigma_nonzero / scalar_nonzero:.2f}x" if scalar_nonzero > 0 else "")

    # Check n, l transitions in the sigma-weighted vertex
    dn_dl_combos = set()
    for mu_idx in range(3):
        for e_idx in range(E_fock):
            V = V_sigma[mu_idx][e_idx]
            for a in range(N_dirac):
                for b in range(N_dirac):
                    if abs(V[a, b]) > 1e-12:
                        dn = labels[a].n_fock - labels[b].n_fock
                        dl = labels[a].l - labels[b].l
                        dk = labels[a].kappa - labels[b].kappa
                        dn_dl_combos.add((dn, dl, dk))
    print(f"  (Delta_n, Delta_l, Delta_kappa) values: {sorted(dn_dl_combos)}")

    # ------------------------------------------------------------------
    # Step 4: Get propagators
    # ------------------------------------------------------------------
    print("\n--- Step 4: Propagators ---")
    photon_data = compute_photon_propagator(n_max, exact=False)
    G_gamma = photon_data.G_gamma_numeric
    print(f"  G_gamma shape: {G_gamma.shape}")
    print(f"  G_gamma rank: {np.linalg.matrix_rank(G_gamma, tol=1e-10)}")

    # Electron propagator at t=0 (diagonal, D^{-1})
    D = np.diag([
        (1.0 if lab.kappa < 0 else -1.0) * (lab.n_fock + 0.5)
        for lab in labels
    ])
    G_e = np.linalg.inv(D)
    print(f"  G_e diagonal: {[f'{G_e[i,i]:.4f}' for i in range(N_dirac)]}")

    # ------------------------------------------------------------------
    # Step 5: Self-energy
    # ------------------------------------------------------------------
    print("\n--- Step 5: Self-energy ---")
    Sigma = compute_sigma_weighted_self_energy(V_sigma, G_gamma, N_dirac, E_fock)

    is_real = bool(np.allclose(Sigma.imag, 0, atol=1e-10))
    print(f"  Is real: {is_real}")

    if is_real:
        Sigma_work = Sigma.real
        is_hermitian = bool(np.allclose(Sigma_work, Sigma_work.T, atol=1e-10))
    else:
        Sigma_work = Sigma
        is_hermitian = bool(np.allclose(Sigma_work, Sigma_work.conj().T, atol=1e-10))
    print(f"  Is Hermitian: {is_hermitian}")

    tr_sigma = np.trace(Sigma_work)
    print(f"  Trace: {tr_sigma:.8f}")

    if is_real and is_hermitian:
        eigs = np.linalg.eigvalsh(Sigma_work)
    else:
        eigs = np.linalg.eigvals(Sigma_work)
    eigs_sorted = sorted(eigs.real)
    print(f"  Eigenvalues: {[f'{e:.6f}' for e in eigs_sorted]}")
    n_zero = sum(1 for e in eigs_sorted if abs(e) < 1e-10)
    print(f"  Zero eigenvalues: {n_zero}")
    psd = all(e > -1e-10 for e in eigs_sorted)
    print(f"  Positive semidefinite: {psd}")

    # Ground state block
    gs_idx = [i for i, lab in enumerate(labels) if lab.n_fock == 1 and lab.kappa == -1]
    gs_block = Sigma_work[np.ix_(gs_idx, gs_idx)]
    gs_max = np.max(np.abs(gs_block))
    gs_zero = bool(gs_max < 1e-10)
    print(f"  GS block (n=1,kappa=-1): max|entry| = {gs_max:.8f}")
    print(f"  GS structural zero: {gs_zero}")
    print(f"  GS block values:")
    for i in range(gs_block.shape[0]):
        vals = [f"{float(gs_block[i,j].real if isinstance(gs_block[i,j], complex) else gs_block[i,j]):.8f}"
                for j in range(gs_block.shape[1])]
        print(f"    [{', '.join(vals)}]")

    # ------------------------------------------------------------------
    # Step 6: Vertex correction and F2
    # ------------------------------------------------------------------
    print("\n--- Step 6: Vertex correction and F_2 ---")
    Lambda_total = compute_sigma_weighted_vertex_correction(
        V_sigma, G_gamma, G_e, N_dirac, E_fock)

    Lambda_real = bool(np.allclose(Lambda_total.imag, 0, atol=1e-10))
    print(f"  Lambda is real: {Lambda_real}")

    tr_lambda = np.trace(Lambda_total)
    print(f"  Tr(Lambda): {tr_lambda:.8f}")

    F2, V_bare, num, den = extract_f2(Lambda_total, V_sigma, G_e, N_dirac, E_fock)
    print(f"  Tr(Lambda) = {num:.8f}")
    print(f"  Tr(V_bare @ G_e) = {den:.8f}")
    if F2 is not None:
        print(f"  F_2 = {float(F2.real if isinstance(F2, complex) else F2):.8f}")
    else:
        print(f"  F_2 = undefined (zero denominator)")

    # Compare to scalar F2 (5*sqrt(2)/3 ~ 2.357)
    scalar_f2 = 5 * np.sqrt(2) / 3
    if F2 is not None:
        f2_val = float(F2.real if isinstance(F2, complex) else F2)
        print(f"  Scalar F_2 = {scalar_f2:.6f} (5*sqrt(2)/3)")
        print(f"  Ratio sigma/scalar F_2 = {f2_val / scalar_f2:.6f}")

    # ------------------------------------------------------------------
    # Step 7: Selection rule census
    # ------------------------------------------------------------------
    print("\n--- Step 7: Selection rule census ---")
    census = selection_rule_census(
        V_sigma, V_e_scalar, sigma_list, labels, fock_data,
        Sigma_work, G_e, G_gamma, N_dirac, E_fock)

    for rule_name in ["delta_mj_conservation", "spatial_parity_E1",
                       "gaunt_cg_sparsity", "vertex_parity",
                       "so4_channel_count", "charge_conjugation",
                       "furry_theorem", "ward_identity"]:
        info = census[rule_name]
        survives = info.get("survives", False)
        verdict = info.get("verdict", "N/A")
        mark = "PASS" if survives else "FAIL"
        print(f"  {rule_name:30s} [{mark}] {verdict}")

    print(f"\n  TOTAL: {census['summary']['fraction']} selection rules survive")

    # ------------------------------------------------------------------
    # Step 8: Number field analysis
    # ------------------------------------------------------------------
    print("\n--- Step 8: Number field analysis ---")
    nf = analyze_number_field(Sigma_work, Lambda_total, F2, labels)
    print(f"  Sigma real: {nf['Sigma_is_real']}")
    print(f"  Sigma rational: {nf['Sigma_appears_rational']}")
    print(f"  Lambda real: {nf['Lambda_is_real']}")
    if F2 is not None:
        print(f"  F2 = {nf.get('F2_value', 'N/A')}")
        print(f"  F2 fraction: {nf.get('F2_fraction', 'N/A')} "
              f"(residual {nf.get('F2_fraction_residual', 'N/A')})")
        print(f"  F2^2 fraction: {nf.get('F2_squared_fraction', 'N/A')} "
              f"(residual {nf.get('F2_squared_residual', 'N/A')})")

    # ------------------------------------------------------------------
    # Step 9: Compare to scalar QED
    # ------------------------------------------------------------------
    print("\n--- Step 9: Comparison to scalar QED ---")
    comp = compare_to_scalar(Sigma, V_e_scalar, G_gamma, labels, N_dirac, E_fock)
    print(f"  Scalar Tr(Sigma): {comp['scalar_trace']:.6f}")
    print(f"  Sigma-weighted Tr(Sigma): {comp['sigma_weighted_trace']:.6f}")
    print(f"  Trace ratio: {comp['trace_ratio']}")
    print(f"  Scalar eigenvalues: {[f'{e:.6f}' for e in comp['scalar_eigenvalues']]}")
    print(f"  Sigma eigenvalues: {[f'{e:.6f}' for e in comp['sigma_eigenvalues']]}")

    if 'gs_block_scalar' in comp:
        print(f"  Scalar GS block: {comp['gs_block_scalar']}")
        print(f"  Sigma GS block: {comp['gs_block_sigma']}")

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    print("\n--- Saving results ---")

    # Prepare JSON-serializable results
    output = {
        "n_max": n_max,
        "N_dirac": N_dirac,
        "E_fock": E_fock,
        "sigma_algebra": sigma_results,
        "sigma_weighted_vertex": {
            "scalar_nonzero_entries": int(scalar_nonzero),
            "sigma_nonzero_entries": int(sigma_nonzero),
            "delta_n_l_kappa_values": [list(x) for x in sorted(dn_dl_combos)],
        },
        "self_energy": {
            "is_real": is_real,
            "is_hermitian": is_hermitian,
            "trace": float(tr_sigma.real if isinstance(tr_sigma, complex) else tr_sigma),
            "eigenvalues": [float(e) for e in eigs_sorted],
            "n_zero_eigenvalues": n_zero,
            "positive_semidefinite": psd,
            "gs_block_zero": gs_zero,
            "gs_block_max_abs": float(gs_max),
            "gs_block_values": [[float(gs_block[i, j].real
                                       if isinstance(gs_block[i, j], complex)
                                       else gs_block[i, j])
                                 for j in range(gs_block.shape[1])]
                                for i in range(gs_block.shape[0])],
        },
        "vertex_correction": {
            "Lambda_is_real": Lambda_real,
            "trace_Lambda": float(tr_lambda.real if isinstance(tr_lambda, complex) else tr_lambda),
            "F2": float(F2.real if isinstance(F2, complex) else F2) if F2 is not None else None,
            "F2_numerator": float(num.real if isinstance(num, complex) else num),
            "F2_denominator": float(den.real if isinstance(den, complex) else den),
            "scalar_F2": scalar_f2,
            "F2_ratio_to_scalar": (float((F2 / scalar_f2).real
                                         if isinstance(F2, complex) else F2 / scalar_f2)
                                   if F2 is not None else None),
        },
        "selection_rules": census,
        "number_field": nf,
        "comparison_to_scalar": comp,
    }

    out_path = os.path.join(os.path.dirname(__file__), "data",
                            "vq3_sigma_weighted_vertex.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # Convert numpy types for JSON
    def convert(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, complex):
            return float(obj.real) if abs(obj.imag) < 1e-12 else [obj.real, obj.imag]
        elif isinstance(obj, set):
            return sorted(obj)
        elif isinstance(obj, DiracLabel):
            return f"({obj.n_fock},{obj.kappa},{obj.two_m_j})"
        return obj

    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2, default=convert)
    print(f"  Saved to {out_path}")

    print("\n" + "=" * 70)
    print("VQ-3 COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
