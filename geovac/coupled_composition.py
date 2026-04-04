"""
Coupled Composition — Cross-Block ERI Replacement for PK (Track CB)
====================================================================

Builds a "coupled composed" qubit Hamiltonian that keeps the block
decomposition of the composed architecture but replaces the
Phillips-Kleinman pseudopotential with explicit cross-block two-electron
integrals.

Key idea: PK was needed because classical solvers couldn't handle
inter-group antisymmetry.  On a quantum computer, antisymmetry is
automatic in second quantization (creation/annihilation operators enforce
it).  We can replace PK with explicit cross-block ERIs and let the
quantum solver handle the physics.

Physics:
- Within-block h1 and ERIs are unchanged from composed architecture
- PK is removed entirely from h1
- Cross-block ERIs are added: <a_i b_j | c_i d_j> (Coulomb-type)
  where electron 1 density stays on center i and electron 2 on center j
- Gaunt selection rules apply to cross-block ERIs (same angular filtering)

Limitations of this scoping implementation:
- Only "Coulomb-type" cross-block ERIs are computed (electron densities
  stay on their respective centers).  "Exchange-type" cross-block ERIs
  (mixed-center charge densities) are not included.  For same-center
  sub-blocks (e.g., Li_core Z=3 and LiH_bond_Li Z=1), these could be
  nonzero but are expected to be smaller than the Coulomb terms.

Author: GeoVac Development Team
Date: April 2026
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from openfermion import jordan_wigner

from geovac.composed_qubit import (
    build_composed_hamiltonian,
    _enumerate_states,
)
from geovac.cross_block_mp2 import compute_cross_block_eri
from geovac.molecular_spec import MolecularSpec
from geovac.qubit_encoding import build_fermion_op_from_integrals


def build_coupled_hamiltonian(
    spec: MolecularSpec,
    n_grid: int = 2000,
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    Build a coupled composed Hamiltonian: block decomposition + cross-block ERIs, no PK.

    Parameters
    ----------
    spec : MolecularSpec
        Molecular specification.  PK parameters in the spec are ignored
        (PK is removed; cross-block ERIs replace it).
    n_grid : int
        Radial grid points for cross-block ERI integration.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        Keys: M, Q, N_pauli, h1, eri, nuclear_repulsion, qubit_op,
        fermion_op, wall_time_s, cross_block_eri_count,
        within_block_eri_count, cross_block_partial_1norm,
        pk_removed_1norm, blocks, comparison.
    """
    t0 = time.perf_counter()

    # ------------------------------------------------------------------
    # 1. Build standard composed Hamiltonian WITHOUT PK
    # ------------------------------------------------------------------
    result_no_pk = build_composed_hamiltonian(
        spec, pk_in_hamiltonian=False, verbose=False,
    )

    M = result_no_pk['M']
    Q = result_no_pk['Q']
    h1_no_pk = result_no_pk['h1'].copy()
    h1_pk = result_no_pk['h1_pk']
    eri_within = result_no_pk['eri'].copy()
    nuclear_repulsion = result_no_pk['nuclear_repulsion']

    if verbose:
        print(f"[coupled] {spec.name}: M={M}, Q={Q}")
        print(f"[coupled] Within-block ERIs: "
              f"{int(np.count_nonzero(np.abs(eri_within) > 1e-15))}")

    # ------------------------------------------------------------------
    # 2. Compute cross-block ERIs
    # ------------------------------------------------------------------
    n_blocks = len(spec.blocks)
    eri_coupled = eri_within.copy()
    cross_block_count = 0
    cross_block_magnitudes: List[float] = []
    cross_block_details: List[Dict[str, Any]] = []

    for i in range(n_blocks):
        for j in range(i + 1, n_blocks):
            t_cross = time.perf_counter()
            # Returns physicist notation: {(p, q, r, s): val}
            # where p,r ∈ block i sub-blocks, q,s ∈ block j sub-blocks
            cross_eris = compute_cross_block_eri(spec, i, j, n_grid=n_grid)

            pair_count = 0
            pair_magnitudes: List[float] = []

            for (p, q, r, s), val in cross_eris.items():
                if abs(val) < 1e-15:
                    continue

                pair_count += 1
                pair_magnitudes.append(abs(val))

                # Physicist <pq|rs> -> Chemist (pr|qs)
                eri_coupled[p, r, q, s] = val
                # Symmetry: (pr|qs) = (qs|pr)
                eri_coupled[q, s, p, r] = val

            cross_block_count += pair_count
            cross_block_magnitudes.extend(pair_magnitudes)

            dt = time.perf_counter() - t_cross
            detail = {
                'block_i': i,
                'block_j': j,
                'label_i': spec.blocks[i].label,
                'label_j': spec.blocks[j].label,
                'n_cross_eri': pair_count,
                'wall_time_s': dt,
            }
            if pair_magnitudes:
                detail['max_magnitude'] = max(pair_magnitudes)
                detail['mean_magnitude'] = np.mean(pair_magnitudes)
            cross_block_details.append(detail)

            if verbose:
                print(f"[coupled] Cross-block ({spec.blocks[i].label} x "
                      f"{spec.blocks[j].label}): {pair_count} ERIs, "
                      f"{dt:.2f}s")

    # Symmetrize full ERI tensor
    eri_coupled = 0.5 * (eri_coupled + eri_coupled.transpose(2, 3, 0, 1))

    within_count = int(np.count_nonzero(np.abs(eri_within) > 1e-15))

    if verbose:
        print(f"[coupled] Cross-block ERIs added: {cross_block_count}")
        print(f"[coupled] Total ERI nonzeros: "
              f"{int(np.count_nonzero(np.abs(eri_coupled) > 1e-15))}")

    # ------------------------------------------------------------------
    # 3. Build fermion operator and JW transform
    # ------------------------------------------------------------------
    if verbose:
        print(f"[coupled] Building fermion operator...")
    fermion_op = build_fermion_op_from_integrals(
        h1_no_pk, eri_coupled, nuclear_repulsion,
    )

    if verbose:
        print(f"[coupled] Jordan-Wigner transform...")
    qubit_op = jordan_wigner(fermion_op)
    N_pauli = len(qubit_op.terms)

    # ------------------------------------------------------------------
    # 4. Compute metrics
    # ------------------------------------------------------------------
    one_norm = sum(abs(c) for c in qubit_op.terms.values())

    # Cross-block partial 1-norm (just the cross-block ERI contribution)
    cross_partial_1norm = sum(cross_block_magnitudes) if cross_block_magnitudes else 0.0

    # PK 1-norm (what was removed)
    pk_entries = np.abs(h1_pk)
    pk_removed_1norm = float(np.sum(pk_entries[pk_entries > 1e-15]))

    # ------------------------------------------------------------------
    # 5. Build comparison with standard composed (with PK)
    # ------------------------------------------------------------------
    result_with_pk = build_composed_hamiltonian(
        spec, pk_in_hamiltonian=True, verbose=False,
    )
    N_pauli_composed = result_with_pk['N_pauli']
    one_norm_composed = sum(abs(c) for c in result_with_pk['qubit_op'].terms.values())

    elapsed = time.perf_counter() - t0

    if verbose:
        print(f"\n{'='*60}")
        print(f"COUPLED COMPOSITION RESULTS — {spec.name}")
        print(f"{'='*60}")
        print(f"  Qubits:          {Q}")
        print(f"  Pauli (composed): {N_pauli_composed}")
        print(f"  Pauli (coupled):  {N_pauli}")
        print(f"  Pauli ratio:      {N_pauli / max(N_pauli_composed, 1):.2f}x")
        print(f"  1-norm (composed): {one_norm_composed:.4f} Ha")
        print(f"  1-norm (coupled):  {one_norm:.4f} Ha")
        print(f"  1-norm ratio:      {one_norm / max(one_norm_composed, 1e-15):.2f}x")
        print(f"  Cross-block ERIs:  {cross_block_count}")
        print(f"  Within-block ERIs: {within_count}")
        print(f"  Cross/within ratio: {cross_block_count / max(within_count, 1):.3f}")
        print(f"  Wall time:         {elapsed:.1f}s")

    results: Dict[str, Any] = {
        'M': M,
        'Q': Q,
        'N_pauli': N_pauli,
        'N_pauli_composed': N_pauli_composed,
        'one_norm': one_norm,
        'one_norm_composed': one_norm_composed,
        'pauli_ratio': N_pauli / max(N_pauli_composed, 1),
        'one_norm_ratio': one_norm / max(one_norm_composed, 1e-15),
        'nuclear_repulsion': nuclear_repulsion,
        'wall_time_s': elapsed,
        'cross_block_eri_count': cross_block_count,
        'within_block_eri_count': within_count,
        'cross_block_partial_1norm': cross_partial_1norm,
        'pk_removed_1norm': pk_removed_1norm,
        'cross_block_details': cross_block_details,
        'h1': h1_no_pk,
        'h1_pk': h1_pk,
        'eri': eri_coupled,
        'qubit_op': qubit_op,
        'fermion_op': fermion_op,
        'blocks': result_no_pk['blocks'],
        'spec_name': spec.name,
    }

    return results


def coupled_fci_energy(
    results: Dict[str, Any],
    n_electrons: int,
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    Diagonalize the coupled Hamiltonian via sector-restricted FCI.

    Uses the h1 and ERI matrices directly in the (N_up, N_down) sector,
    avoiding the full 2^Q Hilbert space.  For M=15 spatial orbitals and
    4 electrons, the FCI dimension is C(15,2)^2 = 11,025.

    Parameters
    ----------
    results : dict
        Output from build_coupled_hamiltonian.
    n_electrons : int
        Total number of electrons.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with E_coupled and comparison metrics.
    """
    from itertools import combinations
    from scipy.sparse import lil_matrix
    from scipy.sparse.linalg import eigsh

    M = results['M']
    h1 = results['h1']
    eri = results['eri']
    nuclear_repulsion = results['nuclear_repulsion']

    n_up = n_electrons // 2
    n_down = n_electrons // 2

    # Generate determinant basis: alpha and beta occupation strings
    alpha_strings = list(combinations(range(M), n_up))
    beta_strings = list(combinations(range(M), n_down))
    n_alpha = len(alpha_strings)
    n_beta = len(beta_strings)
    n_det = n_alpha * n_beta

    if verbose:
        print(f"  Sector-restricted FCI: M={M}, N_el={n_electrons}")
        print(f"  Alpha strings: {n_alpha}, Beta strings: {n_beta}")
        print(f"  FCI dimension: {n_det}")

    # Map occupation tuples to indices
    alpha_idx = {s: i for i, s in enumerate(alpha_strings)}
    beta_idx = {s: i for i, s in enumerate(beta_strings)}

    # Build FCI Hamiltonian in the sector
    H_fci = lil_matrix((n_det, n_det))

    def det_index(a_idx: int, b_idx: int) -> int:
        return a_idx * n_beta + b_idx

    # Diagonal: one-body + Coulomb - Exchange
    for ai, alpha in enumerate(alpha_strings):
        for bi, beta in enumerate(beta_strings):
            I = det_index(ai, bi)
            E_diag = nuclear_repulsion

            # One-body
            for p in alpha:
                E_diag += h1[p, p]
            for p in beta:
                E_diag += h1[p, p]

            # Two-body: alpha-alpha
            for i_idx in range(n_up):
                for j_idx in range(i_idx + 1, n_up):
                    p, q = alpha[i_idx], alpha[j_idx]
                    E_diag += eri[p, p, q, q] - eri[p, q, q, p]

            # Two-body: beta-beta
            for i_idx in range(n_down):
                for j_idx in range(i_idx + 1, n_down):
                    p, q = beta[i_idx], beta[j_idx]
                    E_diag += eri[p, p, q, q] - eri[p, q, q, p]

            # Two-body: alpha-beta (Coulomb only, no exchange)
            for p in alpha:
                for q in beta:
                    E_diag += eri[p, p, q, q]

            H_fci[I, I] = E_diag

    # Off-diagonal: single excitations within alpha
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        for p in alpha:
            for r in range(M):
                if r in alpha_set:
                    continue
                # Single excitation p -> r in alpha
                new_alpha = tuple(sorted((alpha_set - {p}) | {r}))
                if new_alpha not in alpha_idx:
                    continue
                ai_new = alpha_idx[new_alpha]

                # Phase from reordering
                phase = _excitation_phase(alpha, p, r)

                for bi, beta in enumerate(beta_strings):
                    I = det_index(ai, bi)
                    J = det_index(ai_new, bi)

                    # h1 contribution
                    val = phase * h1[r, p]

                    # Two-body: sum over occupied alpha (except p, plus r)
                    new_alpha_set = (alpha_set - {p}) | {r}
                    for q in alpha:
                        if q == p:
                            continue
                        val += phase * (eri[r, p, q, q] - eri[r, q, q, p])
                    # Coulomb with beta
                    for q in beta:
                        val += phase * eri[r, p, q, q]

                    if abs(val) > 1e-14:
                        H_fci[I, J] += val

    # Off-diagonal: single excitations within beta
    for bi, beta in enumerate(beta_strings):
        beta_set = set(beta)
        for p in beta:
            for r in range(M):
                if r in beta_set:
                    continue
                new_beta = tuple(sorted((beta_set - {p}) | {r}))
                if new_beta not in beta_idx:
                    continue
                bi_new = beta_idx[new_beta]

                phase = _excitation_phase(beta, p, r)

                for ai, alpha in enumerate(alpha_strings):
                    I = det_index(ai, bi)
                    J = det_index(ai, bi_new)

                    val = phase * h1[r, p]

                    new_beta_set = (beta_set - {p}) | {r}
                    for q in beta:
                        if q == p:
                            continue
                        val += phase * (eri[r, p, q, q] - eri[r, q, q, p])
                    for q in alpha:
                        val += phase * eri[r, p, q, q]

                    if abs(val) > 1e-14:
                        H_fci[I, J] += val

    # Off-diagonal: double excitations within alpha
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        occ = list(alpha)
        for i1 in range(n_up):
            for i2 in range(i1 + 1, n_up):
                p, q = occ[i1], occ[i2]
                for r in range(M):
                    if r in alpha_set:
                        continue
                    for s in range(r + 1, M):
                        if s in alpha_set:
                            continue
                        new_alpha = tuple(sorted(
                            (alpha_set - {p, q}) | {r, s}))
                        if new_alpha not in alpha_idx:
                            continue
                        ai_new = alpha_idx[new_alpha]

                        phase = _double_excitation_phase(
                            alpha, p, q, r, s)
                        val = phase * (eri[r, p, s, q] - eri[r, q, s, p])

                        if abs(val) > 1e-14:
                            for bi in range(n_beta):
                                I = det_index(ai, bi)
                                J = det_index(ai_new, bi)
                                H_fci[I, J] += val

    # Off-diagonal: double excitations within beta
    for bi, beta in enumerate(beta_strings):
        beta_set = set(beta)
        occ = list(beta)
        for i1 in range(n_down):
            for i2 in range(i1 + 1, n_down):
                p, q = occ[i1], occ[i2]
                for r in range(M):
                    if r in beta_set:
                        continue
                    for s in range(r + 1, M):
                        if s in beta_set:
                            continue
                        new_beta = tuple(sorted(
                            (beta_set - {p, q}) | {r, s}))
                        if new_beta not in beta_idx:
                            continue
                        bi_new = beta_idx[new_beta]

                        phase = _double_excitation_phase(
                            beta, p, q, r, s)
                        val = phase * (eri[r, p, s, q] - eri[r, q, s, p])

                        if abs(val) > 1e-14:
                            for ai in range(n_alpha):
                                I = det_index(ai, bi)
                                J = det_index(ai, bi_new)
                                H_fci[I, J] += val

    # Off-diagonal: alpha-beta double excitations (one in each spin)
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        for p_a in alpha:
            for r_a in range(M):
                if r_a in alpha_set:
                    continue
                new_alpha = tuple(sorted((alpha_set - {p_a}) | {r_a}))
                if new_alpha not in alpha_idx:
                    continue
                ai_new = alpha_idx[new_alpha]
                phase_a = _excitation_phase(alpha, p_a, r_a)

                for bi, beta in enumerate(beta_strings):
                    beta_set = set(beta)
                    for p_b in beta:
                        for r_b in range(M):
                            if r_b in beta_set:
                                continue
                            new_beta = tuple(sorted(
                                (beta_set - {p_b}) | {r_b}))
                            if new_beta not in beta_idx:
                                continue
                            bi_new = beta_idx[new_beta]
                            phase_b = _excitation_phase(beta, p_b, r_b)

                            # Only Coulomb for opposite spin
                            val = (phase_a * phase_b
                                   * eri[r_a, p_a, r_b, p_b])

                            if abs(val) > 1e-14:
                                I = det_index(ai, bi)
                                J = det_index(ai_new, bi_new)
                                H_fci[I, J] += val

    H_fci = H_fci.tocsr()

    if verbose:
        print(f"  FCI matrix nnz: {H_fci.nnz}")
        print(f"  Diagonalizing...")

    n_states = min(6, n_det - 1)
    E_vals, _ = eigsh(H_fci, k=n_states, which='SA')
    E_vals = np.sort(E_vals)
    E_gs = float(E_vals[0])

    if verbose:
        print(f"  E_coupled (FCI ground) = {E_gs:.10f} Ha")
        if len(E_vals) > 1:
            print(f"  First excited          = {E_vals[1]:.10f} Ha")
            print(f"  Gap                    = {E_vals[1] - E_vals[0]:.6f} Ha")

    return {
        'exact_diag': True,
        'E_coupled': E_gs,
        'eigenvalues': E_vals.tolist(),
        'n_electrons': n_electrons,
        'n_det': n_det,
    }


def _excitation_phase(occ: tuple, p: int, r: int) -> int:
    """Phase factor for single excitation p -> r in an occupation string."""
    # Count occupied orbitals between p and r (exclusive)
    lo, hi = min(p, r), max(p, r)
    count = sum(1 for o in occ if lo < o < hi)
    return (-1) ** count


def _double_excitation_phase(
    occ: tuple, p: int, q: int, r: int, s: int,
) -> int:
    """Phase for double excitation {p,q} -> {r,s} (p<q, r<s)."""
    # Build intermediate strings and count transpositions
    occ_set = set(occ)
    # Remove p, then q
    temp1 = tuple(sorted(occ_set - {p}))
    phase1 = _excitation_phase(occ, p, p)  # position of p
    idx_p = list(occ).index(p)

    temp2 = tuple(sorted(occ_set - {p, q}))
    # Insert r, then s
    temp3 = tuple(sorted((occ_set - {p, q}) | {r}))
    temp4 = tuple(sorted((occ_set - {p, q}) | {r, s}))

    # Phase = (-1)^(n_transpositions)
    # Easier: compute via sequential single excitation phases
    # Remove p from occ: count occ between start and p
    phase_remove_p = sum(1 for o in occ if o < p)
    remaining_after_p = tuple(sorted(occ_set - {p}))
    phase_remove_q = sum(1 for o in remaining_after_p if o < q)

    remaining = tuple(sorted(occ_set - {p, q}))
    phase_add_r = sum(1 for o in remaining if o < r)
    remaining_plus_r = tuple(sorted(set(remaining) | {r}))
    phase_add_s = sum(1 for o in remaining_plus_r if o < s)

    total_phase = phase_remove_p + phase_remove_q + phase_add_r + phase_add_s
    return (-1) ** total_phase


def run_coupled_scoping(
    system: str = 'LiH',
    max_n: int = 2,
    n_grid: int = 2000,
    verbose: bool = True,
    save_data: bool = True,
) -> Dict[str, Any]:
    """
    Run the full coupled composition scoping investigation for a system.

    Parameters
    ----------
    system : str
        One of 'LiH', 'BeH2', 'H2O'.
    max_n : int
        Maximum principal quantum number.
    n_grid : int
        Radial grid points.
    verbose : bool
        Print progress.
    save_data : bool
        Save results to debug/data/.

    Returns
    -------
    dict with all investigation results.
    """
    from geovac.composed_qubit import lih_spec, beh2_spec, h2o_spec

    spec_funcs = {
        'LiH': lih_spec,
        'BeH2': beh2_spec,
        'H2O': h2o_spec,
    }
    # Valence electron counts (core electrons are frozen in E_core constant)
    n_electrons_map = {
        'LiH': 2,
        'BeH2': 4,
        'H2O': 8,
    }

    if system not in spec_funcs:
        raise ValueError(f"Unknown system: {system}. Choose from {list(spec_funcs)}")

    spec_func = spec_funcs[system]
    n_electrons = n_electrons_map[system]

    # Build spec with PK (the coupled builder will ignore PK)
    spec = spec_func(max_n_core=max_n, max_n_val=max_n, include_pk=True)

    if verbose:
        print(f"\n{'='*60}")
        print(f"COUPLED COMPOSITION SCOPING — {system}")
        print(f"{'='*60}")

    # Phase 3: Build coupled Hamiltonian
    results = build_coupled_hamiltonian(spec, n_grid=n_grid, verbose=verbose)

    # Phase 4: FCI via sector-restricted diagonalization
    fci_results = coupled_fci_energy(results, n_electrons, verbose=verbose)

    if fci_results.get('exact_diag'):
        # Also get composed FCI for comparison (with PK, same electron count)
        result_with_pk = build_composed_hamiltonian(
            spec, pk_in_hamiltonian=True, verbose=False,
        )
        fci_composed = coupled_fci_energy(
            result_with_pk, n_electrons, verbose=False,
        )
        if fci_composed.get('exact_diag'):
            fci_results['E_composed'] = fci_composed['E_coupled']

        # Also get composed without PK (to isolate PK vs cross-block effect)
        result_no_pk = build_composed_hamiltonian(
            spec, pk_in_hamiltonian=False, verbose=False,
        )
        fci_no_pk = coupled_fci_energy(
            result_no_pk, n_electrons, verbose=False,
        )
        if fci_no_pk.get('exact_diag'):
            fci_results['E_no_pk_no_cross'] = fci_no_pk['E_coupled']

        if verbose:
            print(f"\n  --- FCI Accuracy ({n_electrons} valence electrons) ---")
            print(f"  E_coupled (no PK + cross ERIs)  = {fci_results['E_coupled']:.6f} Ha")
            if 'E_no_pk_no_cross' in fci_results:
                print(f"  E_no_pk_no_cross                = {fci_results['E_no_pk_no_cross']:.6f} Ha")
            if 'E_composed' in fci_results:
                print(f"  E_composed (with PK, no cross)  = {fci_results['E_composed']:.6f} Ha")
            if system == 'LiH':
                E_exact = -8.0706
                print(f"  E_exact (LiH)                   = {E_exact:.6f} Ha")
                err_coupled = abs(fci_results['E_coupled'] - E_exact)
                print(f"  |E_coupled - exact|  = {err_coupled:.6f} Ha "
                      f"({100 * err_coupled / abs(E_exact):.2f}%)")
                if 'E_composed' in fci_results:
                    err_composed = abs(fci_results['E_composed'] - E_exact)
                    print(f"  |E_composed - exact| = {err_composed:.6f} Ha "
                          f"({100 * err_composed / abs(E_exact):.2f}%)")

    # QWC group analysis
    from geovac.measurement_grouping import qwc_groups
    qwc_coupled = qwc_groups(results['qubit_op'])
    n_qwc_coupled = len(qwc_coupled)

    result_with_pk_for_qwc = build_composed_hamiltonian(
        spec, pk_in_hamiltonian=True, verbose=False,
    )
    qwc_composed = qwc_groups(result_with_pk_for_qwc['qubit_op'])
    n_qwc_composed = len(qwc_composed)

    if verbose:
        print(f"\n  --- QWC Groups ---")
        print(f"  Composed: {n_qwc_composed}")
        print(f"  Coupled:  {n_qwc_coupled}")
        print(f"  Ratio:    {n_qwc_coupled / max(n_qwc_composed, 1):.2f}x")

    # Assemble full results
    full_results = {
        'system': system,
        'max_n': max_n,
        'n_grid': n_grid,
        'Q': results['Q'],
        'M': results['M'],
        'N_pauli_coupled': results['N_pauli'],
        'N_pauli_composed': results['N_pauli_composed'],
        'pauli_ratio': results['pauli_ratio'],
        'one_norm_coupled': results['one_norm'],
        'one_norm_composed': results['one_norm_composed'],
        'one_norm_ratio': results['one_norm_ratio'],
        'n_qwc_coupled': n_qwc_coupled,
        'n_qwc_composed': n_qwc_composed,
        'qwc_ratio': n_qwc_coupled / max(n_qwc_composed, 1),
        'cross_block_eri_count': results['cross_block_eri_count'],
        'within_block_eri_count': results['within_block_eri_count'],
        'cross_within_ratio': (results['cross_block_eri_count']
                               / max(results['within_block_eri_count'], 1)),
        'cross_block_partial_1norm': results['cross_block_partial_1norm'],
        'pk_removed_1norm': results['pk_removed_1norm'],
        'cross_block_details': results['cross_block_details'],
        'fci': fci_results,
        'wall_time_s': results['wall_time_s'],
    }

    if save_data:
        data_dir = Path(__file__).parent.parent / 'debug' / 'data'
        data_dir.mkdir(parents=True, exist_ok=True)
        outpath = data_dir / f'coupled_composition_{system.lower()}.json'

        # Serialize (skip non-JSON types)
        save_dict = {k: v for k, v in full_results.items()
                     if not isinstance(v, np.ndarray)}
        with open(outpath, 'w') as f:
            json.dump(save_dict, f, indent=2, default=str)
        if verbose:
            print(f"\n  Data saved to {outpath}")

    return full_results
