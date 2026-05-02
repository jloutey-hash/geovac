"""
Transverse photon self-energy from Paper 30 plaquettes.
=========================================================

THE BIG TEST: replace the scalar photon propagator G_gamma = L_1^+
(longitudinal only) with the transverse propagator G_T = (d_1 d_1^T)^+
built from Paper 30's plaquettes.

Tests whether the 4 missing QED selection rules emerge:
  1. Vertex parity (n1 + n2 + q odd)
  2. SO(4) channel count W = 0 suppression
  3. Ward identity
  4. Charge conjugation

Key hypothesis: the scalar graph QED is "topological QED" (only
longitudinal photon). The physical (transverse) photon lives in the
co-exact Hodge sector, accessible only via 2-cells (plaquettes).

Author: GeoVac research computation (2026-05-01)
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from geovac.lattice import GeometricLattice
from geovac.fock_graph_hodge import FockGraphHodge
from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels, kappa_to_l
from geovac.graph_qed_propagator import DiracGraphOperator, electron_propagator
from geovac.graph_qed_photon import build_fock_graph, compute_photon_propagator
from geovac.graph_qed_vertex import (
    build_projection_matrix,
    build_vertex_tensor,
    vertex_tensor_to_matrices,
)
from geovac.su2_wilson_gauge import enumerate_plaquettes, enumerate_oriented_edges


def build_transverse_propagator(n_max: int) -> Tuple[np.ndarray, Dict]:
    """Build the transverse photon propagator from plaquettes.

    Returns G_T = (d_1 d_1^T)^+ as an E x E matrix, where d_1 is the
    boundary operator from edges to faces (plaquettes from Paper 30).

    Returns
    -------
    G_T : np.ndarray (E x E)
        Transverse photon propagator (pseudoinverse of co-exact Laplacian).
    info : dict
        Diagnostic information.
    """
    # Build Fock graph
    lat = GeometricLattice(max_n=n_max, topological_weights=False)
    V = lat.num_states
    adj_sparse = lat.adjacency
    adj = np.array(adj_sparse.todense(), dtype=float)

    # Get edges (undirected, canonical ordering i < j)
    rows, cols = adj_sparse.nonzero()
    edge_set = set()
    for r, c in zip(rows, cols):
        if r < c:
            edge_set.add((int(r), int(c)))
    edges = sorted(edge_set)
    E = len(edges)
    edge_index = {e: i for i, e in enumerate(edges)}

    # Build incidence matrix B (V x E)
    B = np.zeros((V, E), dtype=float)
    for k, (i, j) in enumerate(edges):
        B[i, k] = 1.0
        B[j, k] = -1.0

    # L_1 = B^T B (exact / longitudinal)
    L1_exact = B.T @ B

    # Enumerate plaquettes
    plaquettes = enumerate_plaquettes(adj, max_length=8, both_orientations=False)
    F = len(plaquettes)

    if F == 0:
        return np.zeros((E, E)), {
            "n_max": n_max, "V": V, "E": E, "F": 0,
            "verdict": "No plaquettes — cannot build transverse propagator."
        }

    # Build d_1 (E x F boundary operator)
    d1 = np.zeros((E, F), dtype=float)
    for f_idx, plaq in enumerate(plaquettes):
        for oe in plaq:
            src, tgt = oe.source, oe.target
            if src < tgt:
                e_key = (src, tgt)
                sign = +1.0
            else:
                e_key = (tgt, src)
                sign = -1.0
            if e_key in edge_index:
                e_idx = edge_index[e_key]
                d1[e_idx, f_idx] += sign

    # Co-exact Laplacian: L1_coexact = d_1 d_1^T (E x E)
    L1_coexact = d1 @ d1.T

    # Full Hodge: Delta_1 = L1_exact + L1_coexact
    Delta_1 = L1_exact + L1_coexact

    # Compute pseudoinverse of L1_coexact (transverse propagator)
    # Use eigendecomposition for stable pseudoinverse
    evals, evecs = np.linalg.eigh(L1_coexact)
    tol = 1e-10
    G_T = np.zeros((E, E))
    for i in range(E):
        if abs(evals[i]) > tol:
            G_T += (1.0 / evals[i]) * np.outer(evecs[:, i], evecs[:, i])

    # Also compute pseudoinverse of full Delta_1 for comparison
    evals_full, evecs_full = np.linalg.eigh(Delta_1)
    G_full = np.zeros((E, E))
    for i in range(E):
        if abs(evals_full[i]) > tol:
            G_full += (1.0 / evals_full[i]) * np.outer(evecs_full[:, i], evecs_full[:, i])

    # Diagnostics
    rank_d1 = int(np.linalg.matrix_rank(d1, tol=1e-8))
    nz_coexact = sorted(evals[evals > tol])
    nz_full = sorted(evals_full[evals_full > tol])

    info = {
        "n_max": n_max,
        "V": V,
        "E": E,
        "F": F,
        "rank_d1": rank_d1,
        "coexact_nonzero_eigenvalues": [round(x, 4) for x in nz_coexact],
        "coexact_zero_count": int(np.sum(np.abs(evals) < tol)),
        "full_hodge_nonzero_eigenvalues": [round(x, 4) for x in nz_full[:10]],
        "G_T_norm": float(np.linalg.norm(G_T)),
        "G_full_norm": float(np.linalg.norm(G_full)),
        "B_incidence": B,
    }

    return G_T, info


def compute_transverse_self_energy(n_max: int) -> Dict:
    """Compute the self-energy using the TRANSVERSE photon propagator.

    Sigma_T[a, b] = sum_{e, e'} V_e[a, :] . G_T[e, e'] . V_{e'}[:, b]^T

    where G_T = (d_1 d_1^T)^+ is the transverse (co-exact) propagator.
    """
    print(f"\n{'='*60}")
    print(f"TRANSVERSE SELF-ENERGY at n_max={n_max}")
    print(f"{'='*60}")

    # Step 1: Build transverse propagator
    print("  Building transverse propagator from plaquettes...")
    G_T, prop_info = build_transverse_propagator(n_max)

    if prop_info["F"] == 0:
        return {"n_max": n_max, "verdict": "No plaquettes available."}

    print(f"    V={prop_info['V']}, E={prop_info['E']}, F={prop_info['F']}")
    print(f"    rank(d_1) = {prop_info['rank_d1']}")
    print(f"    Co-exact eigenvalues: {prop_info['coexact_nonzero_eigenvalues']}")

    # Step 2: Build vertex tensor (CG projection from Dirac -> Fock edges)
    print("  Building vertex coupling tensor...")
    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    fock_data = build_fock_graph(n_max)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats_sympy = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    # Convert to numpy
    V_mats = [np.array(vm.tolist(), dtype=float) for vm in V_mats_sympy]

    print(f"    N_dirac = {N_dirac}, E_fock = {E_fock}")

    # Check E consistency
    assert E_fock == prop_info['E'], (
        f"Edge count mismatch: vertex has {E_fock}, propagator has {prop_info['E']}"
    )

    # Step 3: Compute scalar self-energy (for comparison)
    print("  Computing SCALAR self-energy (L_1^+)...")
    # Build L_1 = B^T B directly in numpy (avoids sympy complex eigenvalue bug at n_max>=4)
    B_inc = prop_info['B_incidence']  # V x E incidence matrix
    L1_np = B_inc.T @ B_inc  # E x E edge Laplacian
    # Pseudoinverse
    evals_L1, evecs_L1 = np.linalg.eigh(L1_np)
    tol_pinv = 1e-10
    G_scalar = np.zeros_like(L1_np)
    for i in range(len(evals_L1)):
        if abs(evals_L1[i]) > tol_pinv:
            G_scalar += (1.0 / evals_L1[i]) * np.outer(evecs_L1[:, i], evecs_L1[:, i])


    Sigma_scalar = np.zeros((N_dirac, N_dirac))
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_ee = G_scalar[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            Sigma_scalar += g_ee * (V_mats[e1] @ V_mats[e2].T)

    # Step 4: Compute TRANSVERSE self-energy
    print("  Computing TRANSVERSE self-energy (d_1 d_1^T)^+...")
    Sigma_transverse = np.zeros((N_dirac, N_dirac))
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_ee = G_T[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            Sigma_transverse += g_ee * (V_mats[e1] @ V_mats[e2].T)

    # Step 5: Analyze both self-energies
    print("\n  ANALYSIS:")

    # Ground state indices
    gs_indices = []
    for i, lab in enumerate(dirac_labels):
        if lab.n_fock == 1 and lab.kappa == -1:
            gs_indices.append(i)

    # Ground state blocks
    gs_block_scalar = Sigma_scalar[np.ix_(gs_indices, gs_indices)]
    gs_block_transverse = Sigma_transverse[np.ix_(gs_indices, gs_indices)]

    # Diagonality check
    def diagonality(M):
        """Fraction of Frobenius norm on the diagonal."""
        diag_norm = np.linalg.norm(np.diag(np.diag(M)))
        full_norm = np.linalg.norm(M)
        if full_norm < 1e-15:
            return 1.0
        return diag_norm / full_norm

    # Eigenvalues
    evals_scalar = np.linalg.eigvalsh(Sigma_scalar)
    evals_transverse = np.linalg.eigvalsh(Sigma_transverse)

    # Zero eigenvalue count (selection rule indicator)
    tol = 1e-8
    n_zero_scalar = int(np.sum(np.abs(evals_scalar) < tol))
    n_zero_transverse = int(np.sum(np.abs(evals_transverse) < tol))

    # Block structure by n-shell
    n_shells = {}
    for i, lab in enumerate(dirac_labels):
        n = lab.n_fock
        if n not in n_shells:
            n_shells[n] = []
        n_shells[n].append(i)

    # Check if Sigma is block-diagonal by n-shell
    cross_shell_norm_scalar = 0.0
    cross_shell_norm_transverse = 0.0
    total_norm_scalar = np.linalg.norm(Sigma_scalar)
    total_norm_transverse = np.linalg.norm(Sigma_transverse)

    for n1 in n_shells:
        for n2 in n_shells:
            if n1 != n2:
                block = Sigma_scalar[np.ix_(n_shells[n1], n_shells[n2])]
                cross_shell_norm_scalar += np.linalg.norm(block)**2
                block_t = Sigma_transverse[np.ix_(n_shells[n1], n_shells[n2])]
                cross_shell_norm_transverse += np.linalg.norm(block_t)**2

    cross_shell_frac_scalar = np.sqrt(cross_shell_norm_scalar) / max(total_norm_scalar, 1e-15)
    cross_shell_frac_transverse = np.sqrt(cross_shell_norm_transverse) / max(total_norm_transverse, 1e-15)

    # Pendant-edge theorem check for scalar
    pendant_expected = 2.0 * (n_max - 1) / n_max
    gs_trace_scalar = np.trace(gs_block_scalar)

    # Print results
    print(f"\n  === SCALAR (longitudinal) self-energy ===")
    print(f"    Tr(Sigma_scalar) = {np.trace(Sigma_scalar):.6f}")
    print(f"    GS block trace = {gs_trace_scalar:.6f} (pendant: {pendant_expected:.6f})")
    print(f"    GS zero? {np.max(np.abs(gs_block_scalar)) < 1e-8}")
    print(f"    Diagonality = {diagonality(Sigma_scalar):.4f}")
    print(f"    Zero eigenvalues: {n_zero_scalar}/{N_dirac}")
    print(f"    Cross-shell fraction: {cross_shell_frac_scalar:.4f}")
    print(f"    Eigenvalues: {sorted(evals_scalar)[:6]}")

    print(f"\n  === TRANSVERSE (co-exact) self-energy ===")
    print(f"    Tr(Sigma_transverse) = {np.trace(Sigma_transverse):.6f}")
    gs_trace_transverse = np.trace(gs_block_transverse)
    print(f"    GS block trace = {gs_trace_transverse:.6f}")
    print(f"    GS zero? {np.max(np.abs(gs_block_transverse)) < 1e-8}")
    print(f"    Diagonality = {diagonality(Sigma_transverse):.4f}")
    print(f"    Zero eigenvalues: {n_zero_transverse}/{N_dirac}")
    print(f"    Cross-shell fraction: {cross_shell_frac_transverse:.4f}")
    print(f"    Eigenvalues: {sorted(evals_transverse)[:6]}")

    # THE KEY SELECTION RULE TESTS
    print(f"\n  === SELECTION RULE TESTS ===")

    # Test 1: Ground state protection (vertex parity analog)
    gs_zero_scalar = np.max(np.abs(gs_block_scalar)) < 1e-8
    gs_zero_transverse = np.max(np.abs(gs_block_transverse)) < 1e-8
    print(f"    [1] Ground state protection (Sigma(GS)=0):")
    print(f"        Scalar:     {'PASS' if gs_zero_scalar else 'FAIL'} (max |GS| = {np.max(np.abs(gs_block_scalar)):.6e})")
    print(f"        Transverse: {'PASS' if gs_zero_transverse else 'FAIL'} (max |GS| = {np.max(np.abs(gs_block_transverse)):.6e})")

    # Test 2: Diagonal structure (conservation laws -> block diagonal)
    diag_scalar = diagonality(Sigma_scalar)
    diag_transverse = diagonality(Sigma_transverse)
    print(f"    [2] Diagonal dominance (conservation -> diag):")
    print(f"        Scalar:     {diag_scalar:.4f}")
    print(f"        Transverse: {diag_transverse:.4f}")

    # Test 3: Cross-shell suppression
    print(f"    [3] Cross-shell suppression (n-shell block-diag):")
    print(f"        Scalar:     {1-cross_shell_frac_scalar:.4f} within-shell")
    print(f"        Transverse: {1-cross_shell_frac_transverse:.4f} within-shell")

    # Test 4: Positive semi-definiteness
    min_eval_scalar = min(evals_scalar)
    min_eval_transverse = min(evals_transverse)
    psd_scalar = min_eval_scalar > -1e-8
    psd_transverse = min_eval_transverse > -1e-8
    print(f"    [4] Positive semi-definite:")
    print(f"        Scalar:     {'PASS' if psd_scalar else 'FAIL'} (min eval = {min_eval_scalar:.6e})")
    print(f"        Transverse: {'PASS' if psd_transverse else 'FAIL'} (min eval = {min_eval_transverse:.6e})")

    # Test 5: Check if transverse Sigma is PROPORTIONAL to continuum structure
    # (diagonal with Sigma(n) growing with n)
    diag_vals_transverse = np.diag(Sigma_transverse)
    diag_by_shell = {}
    for i, lab in enumerate(dirac_labels):
        n = lab.n_fock
        if n not in diag_by_shell:
            diag_by_shell[n] = []
        diag_by_shell[n].append(diag_vals_transverse[i])

    print(f"\n  === SHELL-RESOLVED DIAGONAL ===")
    print(f"    {'Shell':<8} {'Mean Sigma_T':<15} {'Std':<12} {'Mean Sigma_S':<15}")
    diag_vals_scalar = np.diag(Sigma_scalar)
    diag_by_shell_scalar = {}
    for i, lab in enumerate(dirac_labels):
        n = lab.n_fock
        if n not in diag_by_shell_scalar:
            diag_by_shell_scalar[n] = []
        diag_by_shell_scalar[n].append(diag_vals_scalar[i])

    for n in sorted(diag_by_shell.keys()):
        vals_t = diag_by_shell[n]
        vals_s = diag_by_shell_scalar[n]
        print(f"    n={n:<5} {np.mean(vals_t):<15.6f} {np.std(vals_t):<12.6f} {np.mean(vals_s):<15.6f}")

    # Compile results
    results = {
        "n_max": n_max,
        "propagator_info": prop_info,
        "N_dirac": N_dirac,
        "E_fock": E_fock,
        "scalar": {
            "trace": float(np.trace(Sigma_scalar)),
            "gs_block_trace": float(gs_trace_scalar),
            "gs_zero": bool(gs_zero_scalar),
            "gs_block_max": float(np.max(np.abs(gs_block_scalar))),
            "diagonality": float(diag_scalar),
            "n_zero_eigenvalues": n_zero_scalar,
            "cross_shell_fraction": float(cross_shell_frac_scalar),
            "psd": bool(psd_scalar),
            "min_eigenvalue": float(min_eval_scalar),
            "eigenvalues": sorted([float(x) for x in evals_scalar]),
            "pendant_expected": float(pendant_expected),
        },
        "transverse": {
            "trace": float(np.trace(Sigma_transverse)),
            "gs_block_trace": float(gs_trace_transverse),
            "gs_zero": bool(gs_zero_transverse),
            "gs_block_max": float(np.max(np.abs(gs_block_transverse))),
            "diagonality": float(diag_transverse),
            "n_zero_eigenvalues": n_zero_transverse,
            "cross_shell_fraction": float(cross_shell_frac_transverse),
            "psd": bool(psd_transverse),
            "min_eigenvalue": float(min_eval_transverse),
            "eigenvalues": sorted([float(x) for x in evals_transverse]),
        },
        "selection_rules": {
            "gs_protection_scalar": bool(gs_zero_scalar),
            "gs_protection_transverse": bool(gs_zero_transverse),
            "diagonality_scalar": float(diag_scalar),
            "diagonality_transverse": float(diag_transverse),
            "cross_shell_suppression_scalar": float(1 - cross_shell_frac_scalar),
            "cross_shell_suppression_transverse": float(1 - cross_shell_frac_transverse),
            "psd_scalar": bool(psd_scalar),
            "psd_transverse": bool(psd_transverse),
        },
        "shell_diagonal": {
            str(n): {
                "transverse_mean": float(np.mean(diag_by_shell[n])),
                "transverse_std": float(np.std(diag_by_shell[n])),
                "scalar_mean": float(np.mean(diag_by_shell_scalar[n])),
            }
            for n in sorted(diag_by_shell.keys())
        },
    }

    return results


def main():
    """Run the transverse QED self-energy computation."""
    print("=" * 72)
    print("TRANSVERSE PHOTON QED: Testing selection rule recovery")
    print("from Paper 30 plaquettes")
    print("=" * 72)

    t0 = time.time()
    all_results = {}

    # n_max=3: 4 plaquettes, rank(d1)=2
    results_3 = compute_transverse_self_energy(n_max=3)
    all_results["n_max_3"] = results_3

    # n_max=4: 33 plaquettes, rank(d1)=8
    results_4 = compute_transverse_self_energy(n_max=4)
    all_results["n_max_4"] = results_4

    # Summary
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)

    for key in ["n_max_3", "n_max_4"]:
        r = all_results[key]
        if "verdict" in r:
            print(f"\n  {key}: {r['verdict']}")
            continue

        print(f"\n  {key}:")
        print(f"    GS protection recovered? {r['selection_rules']['gs_protection_transverse']}")
        print(f"    Diagonality improved? {r['selection_rules']['diagonality_transverse']:.4f} vs scalar {r['selection_rules']['diagonality_scalar']:.4f}")
        print(f"    Cross-shell suppressed? {r['selection_rules']['cross_shell_suppression_transverse']:.4f} vs scalar {r['selection_rules']['cross_shell_suppression_scalar']:.4f}")

    # Save
    output_path = Path(__file__).parent / "data" / "transverse_qed_self_energy.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def json_serializer(x):
        if isinstance(x, np.floating):
            return float(x)
        if isinstance(x, np.integer):
            return int(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
        return x

    # Strip B_incidence from propagator_info before saving
    for key in all_results:
        if isinstance(all_results[key], dict) and 'propagator_info' in all_results[key]:
            all_results[key]['propagator_info'].pop('B_incidence', None)

    with open(output_path, "w") as f:
        json.dump(all_results, f, indent=2, default=json_serializer)

    elapsed = time.time() - t0
    print(f"\n  Elapsed: {elapsed:.1f}s")
    print(f"  Results saved: {output_path}")

    return all_results


if __name__ == "__main__":
    main()
