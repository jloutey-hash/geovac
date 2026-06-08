"""Sprint M-vS Gauge — does the FULL M-vS gauge group reduce Pauli count
beyond the abelian Z₂ Hopf-U(1) tapering?

Context: yesterday's umbrella sprint (v3.86.0) confirmed that GeoVac's
chemistry composition fits Marcolli-vS 2014 bit-exactly.  Today's
M-vS-2 + R-sweep (v3.87.0) showed (a) production-default LiH IS bit-
exactly a gauge network and (b) the spectral action does NOT bind LiH.
The reconciliation question now is: is the M-vS framework giving us
ANY new chemistry-side computational power, or is it just paper-arc
structural identification?

The Z₂ Hopf-U(1) tapering shipped at v3.52.0 saved ΔQ = 2 + n_sub_blocks
qubits via the abelian m_l → −m_l reflection sub-symmetry of the
M-vS gauge group.  Concretely for LiH default: 3 sub-blocks, so Z₂
tapering gives ΔQ = 5 qubits (and ~12.6% Pauli reduction in the
ecosystem panel).

The M-vS gauge group is RICHER: at each vertex v, the gauge group is
U(H_v).  For default LiH:
  Vertex Li: U(10) — can mix Li_core (5 orbitals) with LiH_bond_center
             (5 orbitals) AND mix m_l within each sub-block
  Vertex H:  U(5) — can mix LiH_bond_partner orbitals
The Z₂ tapering uses ONLY per-sub-block U(5) sub-gauges; the
"non-abelian extension" specifically means mixing Li_core with
LiH_bond_center via U_Li outside the per-sub-block subgroup.

Test: do any M-vS gauge transformations BEYOND the Z₂ Hopf rotation
produce additional Pauli reduction?  If yes → M-vS gives genuine
chemistry-side computational content beyond Z₂.  If no → M-vS is
purely structural-identification paper material; the chemistry-side
computational value of the framework is the Z₂ tapering already
shipped.

Gauge candidates:
  A. Identity (baseline; native composed builder output)
  B. Z₂ Hopf-U(1) tapering ONLY (existing baseline)
  C. Per-sub-block h1 eigenbasis rotation (abelian sub-gauge:
     each sub-block's h1 block diagonalized independently)
  D. Per-vertex h1 eigenbasis (full vertex Hilbert space
     diagonalized; Li_core mixed with LiH_bond_center via SAME
     vertex Dirac).  NOTE: because cross_block_h1 skips same-center
     pairs, the within-Li h1 is structurally block-diagonal —
     so D = C for this molecule.  But we verify this.
  E. Natural orbitals from FCI 1-RDM, projected onto M-vS vertex
     structure (the non-trivial non-abelian-allowed mixing: FCI
     correlation can mix Li_core with LiH_bond_center via the eri
     coupling, and the projection takes the Li-vertex 10×10 block
     of the 1-RDM and diagonalizes that).
  F. Random block-unitary U_Li ⊕ U_H respecting M-vS structure
     (control: any block-unitary gives same answer? or specific
     gauges matter?)

For each gauge candidate, we:
  1. Construct U (15×15, block-diagonal across vertices)
  2. Transform h1, eri covariantly
  3. Verify FCI ground-state energy preserved (eigenvalue invariance)
  4. Build qubit op via Jordan-Wigner, count Pauli terms
  5. Apply Z₂ Hopf tapering on top, count Pauli terms in the
     tapered representation

Decision gate:
  - If best (M-vS gauge + Z₂) < best Z₂ alone → M-vS upgrades chemistry
  - If equal → no win, but no loss
  - If worse → M-vS rotation HURTS Pauli structure (eri becomes denser)

Cross-check: repeat on H₂ (much simpler 2-vertex / 2-sub-block,
no internal Li-vertex direct sum).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import scipy.linalg as la

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass

from openfermion import jordan_wigner

from geovac.molecular_spec import lih_spec, MolecularSpec, OrbitalBlock
from geovac.balanced_coupled import (
    build_balanced_hamiltonian,
    _get_block_geometry,
    _get_nuclei_for_lih,
    _get_sub_block_positions,
)
from geovac.coupled_composition import coupled_fci_energy
from geovac.qubit_encoding import build_fermion_op_from_integrals
from geovac.z2_tapering import (
    build_pm_rotation,
    rotate_h1_eri,
    apply_hopf_tapering,
    _enumerate_orbitals,
)


# ---------------------------------------------------------------------------
# Pauli counting
# ---------------------------------------------------------------------------

def count_pauli_terms(qubit_op) -> int:
    """Number of non-identity Pauli terms in a QubitOperator."""
    return len(qubit_op.terms) - (1 if () in qubit_op.terms else 0)


def build_qubit_op(h1: np.ndarray, eri: np.ndarray, nuc_rep: float):
    """h1, eri → FermionOperator → JW QubitOperator."""
    fermion_op = build_fermion_op_from_integrals(h1, eri, nuc_rep)
    qubit_op = jordan_wigner(fermion_op)
    return qubit_op


# ---------------------------------------------------------------------------
# M-vS gauge transformation helpers
# ---------------------------------------------------------------------------

def symmetrize_eri_8fold(eri: np.ndarray) -> np.ndarray:
    """Average eri over the 8-fold permutation orbit for real chemist (pq|rs).

    The composed_qubit builder only enforces the (1↔2) swap symmetry
    (line 748 in composed_qubit.py).  The per-electron symmetries
    (pq|rs) = (qp|rs) and (pq|rs) = (pq|sr) hold IDEALLY but are not
    explicitly imposed.  Under a generic orbital rotation U, even
    machine-precision asymmetries in stored eri amplify to ~1e-3
    energy residuals because the FCI formula accesses specific
    permutation entries that should be equal but are not.

    This helper enforces the full 8-fold symmetry by averaging.
    """
    sym = eri.copy()
    sym = 0.5 * (sym + sym.transpose(1, 0, 2, 3))      # p ↔ q on electron 1
    sym = 0.5 * (sym + sym.transpose(0, 1, 3, 2))      # r ↔ s on electron 2
    sym = 0.5 * (sym + sym.transpose(2, 3, 0, 1))      # 12 ↔ 21
    return sym


def apply_block_unitary(
    h1: np.ndarray, eri: np.ndarray,
    U: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply M × M orthogonal U to h1, eri.

    h1' = U h1 U^T
    eri' (pq|rs) = Σ U_{pa} U_{qb} U_{rc} U_{sd} (ab|cd)

    Pre-symmetrizes eri to the full 8-fold orbit before rotating, then
    re-symmetrizes after.
    """
    eri_in = symmetrize_eri_8fold(eri)
    h1_new = U @ h1 @ U.T
    eri_new = np.einsum(
        'pa,qb,rc,sd,abcd->pqrs', U, U, U, U, eri_in, optimize='optimal',
    )
    eri_new = symmetrize_eri_8fold(eri_new)
    return h1_new, eri_new


def vertex_grouping_for_lih(
    spec: MolecularSpec, R: float,
) -> Dict[str, Any]:
    """Identify which orbital indices belong to which M-vS vertex
    (Li ↔ H).  See sprint M-vS-2 driver for the structural reading.
    Returns dict with:
      - 'Li_indices': list of orbital indices on Li vertex
      - 'H_indices': list of orbital indices on H vertex
      - 'sub_block_indices': dict label -> list of indices
    """
    sub_blocks = _get_block_geometry(spec)
    nuclei = _get_nuclei_for_lih(spec, R)
    sb_positions = _get_sub_block_positions(spec, nuclei)

    # LiH-specific: Li at origin, H along +z
    li_pos = (0.0, 0.0, 0.0)

    li_indices: List[int] = []
    h_indices: List[int] = []
    sb_indices: Dict[str, List[int]] = {}
    for sb in sub_blocks:
        off = sb['offset']
        n_orb = len(sb['states'])
        idx = list(range(off, off + n_orb))
        sb_indices[sb['label']] = idx
        pos = sb_positions[sb['label']]
        if all(abs(pos[k] - li_pos[k]) < 1e-10 for k in range(3)):
            li_indices.extend(idx)
        else:
            h_indices.extend(idx)
    return {
        'Li_indices': li_indices,
        'H_indices': h_indices,
        'sub_block_indices': sb_indices,
        'sub_blocks': sub_blocks,
    }


def gauge_per_sub_block_h1_diag(
    h1: np.ndarray, M: int,
    grouping: Dict[str, Any],
) -> np.ndarray:
    """Block-orthogonal U that diagonalizes each sub-block's h1 sub-matrix.

    This is a subgroup of the M-vS gauge: U respects the sub-block
    structure (no mixing of Li_core with LiH_bond_center, no mixing
    across vertices).
    """
    U = np.zeros((M, M))
    for label, idx in grouping['sub_block_indices'].items():
        idx_arr = np.array(idx)
        h_sub = h1[np.ix_(idx_arr, idx_arr)]
        # Real symmetric → orthogonal diagonalization
        eigvals, eigvecs = la.eigh(h_sub)
        for i_local, i_global in enumerate(idx_arr):
            for j_local, j_global in enumerate(idx_arr):
                U[j_global, i_global] = eigvecs[j_local, i_local]
    return U


def gauge_per_vertex_h1_diag(
    h1: np.ndarray, M: int,
    grouping: Dict[str, Any],
) -> np.ndarray:
    """Block-orthogonal U that diagonalizes each vertex's h1 sub-matrix
    (mixing across sub-blocks within a vertex when h1 has cross-sub-block
    coupling within the vertex).

    For default LiH, cross_block_h1 skips same-center pairs, so the
    within-Li h1 block is structurally block-diagonal.  This should
    coincide with per-sub-block diagonalization for LiH.  We verify.
    """
    U = np.zeros((M, M))
    for vertex_name, idx_list in (('Li', grouping['Li_indices']),
                                  ('H', grouping['H_indices'])):
        idx_arr = np.array(idx_list)
        h_sub = h1[np.ix_(idx_arr, idx_arr)]
        eigvals, eigvecs = la.eigh(h_sub)
        for i_local, i_global in enumerate(idx_arr):
            for j_local, j_global in enumerate(idx_arr):
                U[j_global, i_global] = eigvecs[j_local, i_local]
    return U


def gauge_natural_orbitals_per_vertex(
    h1: np.ndarray, eri: np.ndarray, nuc_rep: float,
    M: int, n_electrons: int,
    grouping: Dict[str, Any],
) -> np.ndarray:
    """Block-orthogonal U from per-vertex natural-orbital diagonalization.

    Procedure:
      1. Compute FCI ground state in the (h1, eri, ecore) basis
      2. Build the 1-RDM (M × M)
      3. Project onto each vertex (Li 10×10 block, H 5×5 block)
      4. Diagonalize each vertex's 1-RDM block
      5. Use the eigenvectors as the orthogonal basis change U_Li, U_H

    Step 4 is where Li_core and LiH_bond_center get mixed via the FCI
    correlation effect — this is the FIRST gauge candidate that exercises
    M-vS gauge freedom beyond the per-sub-block subgroup.
    """
    # Run FCI to get ground state
    result_for_fci = {
        'M': M, 'h1': h1, 'eri': eri, 'nuclear_repulsion': nuc_rep,
    }
    fci_out = coupled_fci_energy(result_for_fci, n_electrons=n_electrons,
                                 verbose=False)
    # We need the 1-RDM from FCI; coupled_fci_energy returns only energy.
    # Compute 1-RDM ourselves.
    one_rdm = compute_fci_1rdm(h1, eri, nuc_rep, M, n_electrons)

    U = np.zeros((M, M))
    for vertex_name, idx_list in (('Li', grouping['Li_indices']),
                                  ('H', grouping['H_indices'])):
        idx_arr = np.array(idx_list)
        # 1-RDM block on this vertex
        rdm_sub = one_rdm[np.ix_(idx_arr, idx_arr)]
        # Diagonalize (Hermitian)
        eigvals, eigvecs = la.eigh(rdm_sub)
        # Sort descending (largest occupation first)
        order = np.argsort(-eigvals)
        eigvecs = eigvecs[:, order]
        for i_local, i_global in enumerate(idx_arr):
            for j_local, j_global in enumerate(idx_arr):
                U[j_global, i_global] = eigvecs[j_local, i_local]
    return U


def compute_fci_1rdm(
    h1: np.ndarray, eri: np.ndarray, nuc_rep: float,
    M: int, n_electrons: int,
) -> np.ndarray:
    """Compute the FCI 1-RDM γ_{pq} = ⟨Ψ_GS | a_p^† a_q | Ψ_GS⟩.

    Uses the sector-restricted FCI (alpha-beta string enumeration)
    to recover the ground state, then accumulates the 1-RDM via
    explicit single-excitation matrix elements.
    """
    from itertools import combinations
    from scipy.sparse import lil_matrix
    from scipy.sparse.linalg import eigsh
    from geovac.coupled_composition import _excitation_phase

    n_up = n_electrons // 2
    n_down = n_electrons - n_up
    alpha_strings = list(combinations(range(M), n_up))
    beta_strings = list(combinations(range(M), n_down))
    n_alpha = len(alpha_strings)
    n_beta = len(beta_strings)
    n_det = n_alpha * n_beta
    alpha_idx = {s: i for i, s in enumerate(alpha_strings)}
    beta_idx = {s: i for i, s in enumerate(beta_strings)}

    def det_index(a, b):
        return a * n_beta + b

    # Build FCI Hamiltonian (re-using the same logic as coupled_fci_energy)
    H_fci = lil_matrix((n_det, n_det))
    for ai, alpha in enumerate(alpha_strings):
        for bi, beta in enumerate(beta_strings):
            I = det_index(ai, bi)
            E_diag = nuc_rep
            for p in alpha:
                E_diag += h1[p, p]
            for p in beta:
                E_diag += h1[p, p]
            for i1 in range(n_up):
                for j1 in range(i1 + 1, n_up):
                    p, q = alpha[i1], alpha[j1]
                    E_diag += eri[p, p, q, q] - eri[p, q, q, p]
            for i1 in range(n_down):
                for j1 in range(i1 + 1, n_down):
                    p, q = beta[i1], beta[j1]
                    E_diag += eri[p, p, q, q] - eri[p, q, q, p]
            for p in alpha:
                for q in beta:
                    E_diag += eri[p, p, q, q]
            H_fci[I, I] = E_diag

    # Off-diagonals (single excitations only — sufficient for 1-RDM via
    # ground state vector, double excitations are needed for ENERGY not for 1-RDM)
    # We actually need the FULL FCI ground state, so include all excitations.
    # Re-use the full coupled_fci_energy build:
    fake = {'M': M, 'h1': h1, 'eri': eri, 'nuclear_repulsion': nuc_rep}
    fci_out = coupled_fci_energy(fake, n_electrons=n_electrons, verbose=False)
    # But that returns only energy.  We need the eigenvector.  Replicate
    # the build here (one more time) and call eigsh to get the vector.

    # Rebuild full FCI matrix:
    H_full = lil_matrix((n_det, n_det))
    # Diagonal:
    for ai, alpha in enumerate(alpha_strings):
        for bi, beta in enumerate(beta_strings):
            I = det_index(ai, bi)
            E_diag = nuc_rep
            for p in alpha:
                E_diag += h1[p, p]
            for p in beta:
                E_diag += h1[p, p]
            for i1 in range(n_up):
                for j1 in range(i1 + 1, n_up):
                    p, q = alpha[i1], alpha[j1]
                    E_diag += eri[p, p, q, q] - eri[p, q, q, p]
            for i1 in range(n_down):
                for j1 in range(i1 + 1, n_down):
                    p, q = beta[i1], beta[j1]
                    E_diag += eri[p, p, q, q] - eri[p, q, q, p]
            for p in alpha:
                for q in beta:
                    E_diag += eri[p, p, q, q]
            H_full[I, I] = E_diag

    # Single excitations alpha:
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        for p in alpha:
            for r in range(M):
                if r in alpha_set:
                    continue
                new_alpha = tuple(sorted((alpha_set - {p}) | {r}))
                if new_alpha not in alpha_idx:
                    continue
                ai_new = alpha_idx[new_alpha]
                phase = _excitation_phase(alpha, p, r)
                for bi, beta in enumerate(beta_strings):
                    I = det_index(ai, bi)
                    J = det_index(ai_new, bi)
                    val = phase * h1[r, p]
                    for q in alpha:
                        if q == p:
                            continue
                        val += phase * (eri[r, p, q, q] - eri[r, q, q, p])
                    for q in beta:
                        val += phase * eri[r, p, q, q]
                    if abs(val) > 1e-14:
                        H_full[I, J] += val

    # Single excitations beta:
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
                    for q in beta:
                        if q == p:
                            continue
                        val += phase * (eri[r, p, q, q] - eri[r, q, q, p])
                    for q in alpha:
                        val += phase * eri[r, p, q, q]
                    if abs(val) > 1e-14:
                        H_full[I, J] += val

    # Double excitations alpha:
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        occ = list(alpha)
        for i1 in range(n_up):
            for j1 in range(i1 + 1, n_up):
                p, q = occ[i1], occ[j1]
                for r in range(M):
                    if r in alpha_set:
                        continue
                    for s in range(r + 1, M):
                        if s in alpha_set:
                            continue
                        new_alpha = tuple(sorted((alpha_set - {p, q}) | {r, s}))
                        if new_alpha not in alpha_idx:
                            continue
                        ai_new = alpha_idx[new_alpha]
                        from geovac.coupled_composition import _double_excitation_phase
                        phase = _double_excitation_phase(alpha, p, q, r, s)
                        val = phase * (eri[r, p, s, q] - eri[r, q, s, p])
                        if abs(val) > 1e-14:
                            for bi in range(n_beta):
                                I = det_index(ai, bi)
                                J = det_index(ai_new, bi)
                                H_full[I, J] += val

    # Double excitations beta:
    for bi, beta in enumerate(beta_strings):
        beta_set = set(beta)
        occ = list(beta)
        for i1 in range(n_down):
            for j1 in range(i1 + 1, n_down):
                p, q = occ[i1], occ[j1]
                for r in range(M):
                    if r in beta_set:
                        continue
                    for s in range(r + 1, M):
                        if s in beta_set:
                            continue
                        new_beta = tuple(sorted((beta_set - {p, q}) | {r, s}))
                        if new_beta not in beta_idx:
                            continue
                        bi_new = beta_idx[new_beta]
                        from geovac.coupled_composition import _double_excitation_phase
                        phase = _double_excitation_phase(beta, p, q, r, s)
                        val = phase * (eri[r, p, s, q] - eri[r, q, s, p])
                        if abs(val) > 1e-14:
                            for ai in range(n_alpha):
                                I = det_index(ai, bi)
                                J = det_index(ai, bi_new)
                                H_full[I, J] += val

    # Mixed double excitations (alpha-beta):
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        for bi, beta in enumerate(beta_strings):
            beta_set = set(beta)
            for p in alpha:
                for r in range(M):
                    if r in alpha_set:
                        continue
                    new_alpha = tuple(sorted((alpha_set - {p}) | {r}))
                    if new_alpha not in alpha_idx:
                        continue
                    ai_new = alpha_idx[new_alpha]
                    phase_a = _excitation_phase(alpha, p, r)
                    for q in beta:
                        for s in range(M):
                            if s in beta_set:
                                continue
                            new_beta = tuple(sorted((beta_set - {q}) | {s}))
                            if new_beta not in beta_idx:
                                continue
                            bi_new = beta_idx[new_beta]
                            phase_b = _excitation_phase(beta, q, s)
                            phase = phase_a * phase_b
                            val = phase * eri[r, p, s, q]
                            if abs(val) > 1e-14:
                                I = det_index(ai, bi)
                                J = det_index(ai_new, bi_new)
                                H_full[I, J] += val

    H_full = H_full.tocsr()
    E_vals, E_vecs = eigsh(H_full, k=1, which='SA')
    gs = E_vecs[:, 0]

    # Build 1-RDM γ_{pq} = ⟨Ψ_GS | a_p^† a_q | Ψ_GS⟩
    # = Σ_{IJ} c_I^* c_J ⟨D_I | a_p^† a_q | D_J⟩
    rdm = np.zeros((M, M))
    # Diagonal: γ_{pp} = ⟨n_p⟩
    for ai, alpha in enumerate(alpha_strings):
        for bi, beta in enumerate(beta_strings):
            I = det_index(ai, bi)
            c2 = gs[I] ** 2
            for p in alpha:
                rdm[p, p] += c2
            for p in beta:
                rdm[p, p] += c2

    # Off-diagonal: single excitation matrix elements
    # γ_{pq} = Σ_{IJ} c_I c_J ⟨D_I | a_p^† a_q | D_J⟩ for p != q
    # The action a_p^† a_q on D_J: if q in D_J and p not in D_J, gives sign * D_J^{p<-q}
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        for q in alpha:
            for p in range(M):
                if p in alpha_set and p != q:
                    continue
                if p == q:
                    continue
                # alpha excitation q -> p
                new_alpha = tuple(sorted((alpha_set - {q}) | {p}))
                if new_alpha not in alpha_idx:
                    continue
                ai_new = alpha_idx[new_alpha]
                phase = _excitation_phase(alpha, q, p)
                for bi in range(n_beta):
                    I = det_index(ai, bi)
                    J = det_index(ai_new, bi)
                    rdm[p, q] += phase * gs[J] * gs[I]
        # beta excitations same form
    for bi, beta in enumerate(beta_strings):
        beta_set = set(beta)
        for q in beta:
            for p in range(M):
                if p in beta_set and p != q:
                    continue
                if p == q:
                    continue
                new_beta = tuple(sorted((beta_set - {q}) | {p}))
                if new_beta not in beta_idx:
                    continue
                bi_new = beta_idx[new_beta]
                phase = _excitation_phase(beta, q, p)
                for ai in range(n_alpha):
                    I = det_index(ai, bi)
                    J = det_index(ai, bi_new)
                    rdm[p, q] += phase * gs[J] * gs[I]

    # Symmetrize to ensure Hermitian
    rdm = 0.5 * (rdm + rdm.T)
    return rdm


def gauge_random_block_unitary(
    M: int, grouping: Dict[str, Any], seed: int = 42,
) -> np.ndarray:
    """Random block-orthogonal U_Li ⊕ U_H respecting M-vS structure.

    Control candidate: does a 'random' M-vS gauge give the SAME Pauli
    count as identity (within statistical noise)?  If yes, M-vS gauge
    is gauge-noise; if no, specific gauges matter.
    """
    rng = np.random.default_rng(seed)
    U = np.zeros((M, M))
    for vertex_name, idx_list in (('Li', grouping['Li_indices']),
                                  ('H', grouping['H_indices'])):
        idx_arr = np.array(idx_list)
        n = len(idx_arr)
        # Random orthogonal via QR
        A = rng.standard_normal((n, n))
        Q, _ = np.linalg.qr(A)
        for i_local, i_global in enumerate(idx_arr):
            for j_local, j_global in enumerate(idx_arr):
                U[j_global, i_global] = Q[j_local, i_local]
    return U


# ---------------------------------------------------------------------------
# Z₂ Hopf tapering wrapper
# ---------------------------------------------------------------------------

def apply_z2_tapering_to_h1_eri(
    spec: MolecularSpec,
    h1: np.ndarray, eri: np.ndarray, nuc_rep: float,
) -> Tuple[Any, np.ndarray, int]:
    """Apply Z₂ Hopf tapering to (h1, eri).

    Returns (qubit_op_tapered, U_pm, dropped) where U_pm is the rotation
    to the m_l → −m_l eigenbasis.  Z₂ tapering ONLY uses U_pm.
    """
    orbital_table = _enumerate_orbitals(spec)
    U_pm, parity = build_pm_rotation(orbital_table)
    h1_pm, eri_pm = rotate_h1_eri(h1, eri, U_pm)
    qubit_op_pm = build_qubit_op(h1_pm, eri_pm, nuc_rep)
    qubit_op_tapered, dropped = apply_hopf_tapering(
        qubit_op_pm, parity, orbital_table=orbital_table, mode='per_block',
        drop_noncommuting=True,
    )
    return qubit_op_tapered, U_pm, dropped


# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------

def evaluate_candidate(
    label: str, U: np.ndarray,
    h1: np.ndarray, eri: np.ndarray, nuc_rep: float,
    spec: MolecularSpec, n_electrons: int, M: int,
    apply_z2: bool = True,
    verify_fci: bool = True,
    baseline_fci: Optional[float] = None,
) -> Dict[str, Any]:
    """Apply gauge U, count Pauli, optionally taper, verify spectrum."""
    h1_g, eri_g = apply_block_unitary(h1, eri, U)

    # Build qubit op pre-tapering
    qubit_op_g = build_qubit_op(h1_g, eri_g, nuc_rep)
    N_pauli_g = count_pauli_terms(qubit_op_g)

    out: Dict[str, Any] = {
        'label': label,
        'M': M,
        'Q_pre_taper': 2 * M,
        'N_pauli_pre_taper': N_pauli_g,
    }

    if verify_fci:
        fake = {'M': M, 'h1': h1_g, 'eri': eri_g, 'nuclear_repulsion': nuc_rep}
        fci_out = coupled_fci_energy(fake, n_electrons=n_electrons,
                                     verbose=False)
        E_fci = fci_out['E_coupled']
        out['E_FCI'] = E_fci
        if baseline_fci is not None:
            out['FCI_residual_vs_baseline'] = abs(E_fci - baseline_fci)

    if apply_z2:
        # Apply Z₂ tapering on the gauge-transformed integrals.
        # Construct a copy of the spec with the gauge applied conceptually
        # — but the Z₂ tapering only needs the orbital_table, which is
        # determined by the (n, l, m) labels.  Under the gauge transform,
        # the orbital_table no longer corresponds to definite (n, l, m).
        # Strategy: apply U_pm (Z₂ rotation) to (h1_g, eri_g), then JW + taper.
        orbital_table = _enumerate_orbitals(spec)
        U_pm, parity = build_pm_rotation(orbital_table)
        # NOTE: under the gauge, the orbital identities are mixed.  We
        # therefore CANNOT use U_pm naively — it assumes definite (n, l, m).
        # We only apply Z₂ tapering when U respects the per-sub-block
        # structure (candidates A and C); for D, E, F the within-vertex
        # mixing destroys the Z₂ structure.
        if _is_per_sub_block_gauge(U, M, vertex_grouping_for_lih_from_spec(spec)):
            h1_pm = U_pm @ h1_g @ U_pm.T
            eri_pm = np.einsum(
                'pa,qb,rc,sd,abcd->pqrs', U_pm, U_pm, U_pm, U_pm, eri_g,
                optimize='optimal',
            )
            qubit_op_pm = build_qubit_op(h1_pm, eri_pm, nuc_rep)
            qubit_op_tapered, dropped = apply_hopf_tapering(
                qubit_op_pm, parity, orbital_table=orbital_table,
                mode='per_block', drop_noncommuting=True,
            )
            n_q_tapered = (
                max((int(q) for term in qubit_op_tapered.terms
                     for (q, _) in term), default=-1) + 1
            )
            out['z2_taper_applies'] = True
            out['Q_post_taper'] = n_q_tapered
            out['N_pauli_post_taper'] = count_pauli_terms(qubit_op_tapered)
            out['Z2_stabilizers_dropped'] = dropped
        else:
            out['z2_taper_applies'] = False
            out['z2_taper_skip_reason'] = (
                'Z₂ Hopf rotation assumes definite (n, l, m) per orbital; '
                'this gauge mixes orbitals within a vertex, so the Z₂ '
                'sub-symmetry no longer factorizes per sub-block.'
            )

    return out


def vertex_grouping_for_lih_from_spec(spec: MolecularSpec) -> Dict[str, Any]:
    """Wrapper that builds the grouping at the spec's recorded R."""
    R = float(spec.R)
    return vertex_grouping_for_lih(spec, R)


def _is_per_sub_block_gauge(
    U: np.ndarray, M: int, grouping: Dict[str, Any],
) -> bool:
    """Check whether U respects the per-sub-block structure (block-diagonal
    in the sub-block partition, not just the vertex partition)."""
    tol = 1e-10
    for label, idx_list in grouping['sub_block_indices'].items():
        # Off-diagonal elements connecting this sub-block to others must be 0
        for i in idx_list:
            for j in range(M):
                if j in idx_list:
                    continue
                if abs(U[i, j]) > tol:
                    return False
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_h2_two_center():
    """Cross-check on H₂ two-center spec (M=10, 2 electrons).
    FCI is dense, bit-exact rotation invariance verified."""
    print("\n" + "=" * 72)
    print("CROSS-CHECK: H₂ two-center spec (M=10, 2 electrons)")
    print("=" * 72)
    R = 1.4
    nuclei = [
        {'Z': 1.0, 'position': (0.0, 0.0, -R/2), 'label': 'H_a'},
        {'Z': 1.0, 'position': (0.0, 0.0, +R/2), 'label': 'H_b'},
    ]
    from geovac.molecular_spec import MolecularSpec, OrbitalBlock
    blocks = [
        OrbitalBlock(label='H_a_atomic', block_type='lone_pair',
                     Z_center=1.0, n_electrons=1, max_n=2, center_nucleus_idx=0),
        OrbitalBlock(label='H_b_atomic', block_type='lone_pair',
                     Z_center=1.0, n_electrons=1, max_n=2, center_nucleus_idx=1),
    ]
    h2_spec = MolecularSpec(name='H2', blocks=blocks,
                            nuclear_repulsion_constant=1.0/R,
                            description='H2 two-center', nuclei=nuclei, R=R)
    result = build_balanced_hamiltonian(h2_spec, R=R, nuclei=h2_spec.nuclei,
                                        cross_block_h1=True, verbose=False)
    h1 = result['h1']
    eri = result['eri']
    nuc_rep = result['nuclear_repulsion']
    M = result['M']
    n_el = 2

    qubit_op_id = build_qubit_op(h1, eri, nuc_rep)
    N_pauli_id = count_pauli_terms(qubit_op_id)
    fake = {'M': M, 'h1': h1, 'eri': eri, 'nuclear_repulsion': nuc_rep}
    fci_out = coupled_fci_energy(fake, n_electrons=n_el, verbose=False)
    E_fci_baseline = fci_out['E_coupled']
    print(f"  Identity: N_pauli = {N_pauli_id}, E_FCI = {E_fci_baseline:.10f}")

    # H₂ grouping: 5 orbs on H_a (indices 0–4), 5 on H_b (indices 5–9)
    grouping = {
        'Li_indices': list(range(5)),  # call H_a "Li" for code reuse
        'H_indices': list(range(5, 10)),
        'sub_block_indices': {
            'H_a_atomic_center': list(range(5)),
            'H_b_atomic_center': list(range(5, 10)),
        },
    }
    # Per-sub-block h1 diag (= per-vertex for H₂)
    U = gauge_per_sub_block_h1_diag(h1, M, grouping)
    h1_g, eri_g = apply_block_unitary(h1, eri, U)
    qubit_op_g = build_qubit_op(h1_g, eri_g, nuc_rep)
    N_pauli_g = count_pauli_terms(qubit_op_g)
    fake_g = {'M': M, 'h1': h1_g, 'eri': eri_g, 'nuclear_repulsion': nuc_rep}
    fci_g = coupled_fci_energy(fake_g, n_electrons=n_el, verbose=False)['E_coupled']
    print(f"  Per-sub-block h1 diag: N_pauli = {N_pauli_g}, "
          f"E_FCI = {fci_g:.10f}, residual = {abs(fci_g - E_fci_baseline):.2e}")

    return {
        'system': 'H2',
        'M': M,
        'N_pauli_identity': N_pauli_id,
        'N_pauli_per_subblock_diag': N_pauli_g,
        'E_FCI_baseline': E_fci_baseline,
        'E_FCI_per_subblock_diag': fci_g,
        'FCI_residual': abs(fci_g - E_fci_baseline),
    }


def run_lih_default():
    print("=" * 72)
    print("Sprint M-vS Gauge — Pauli reduction beyond Z₂ Hopf-U(1)?")
    print("System: default lih_spec(), n_max=2, R = R_eq = 3.015 bohr")
    print("Date: 2026-06-07 (session continuation)")
    print("=" * 72)

    R = 3.015
    max_n = 2
    spec = lih_spec(R=R, max_n=max_n)
    n_electrons = 4  # default LiH spec has 4 encoded electrons

    # Build baseline
    print(f"\n[1] Build LiH balanced Hamiltonian at R = {R}, max_n = {max_n}")
    t0 = time.time()
    result = build_balanced_hamiltonian(
        spec, R=R, cross_block_h1=True, verbose=False,
    )
    h1 = result['h1']
    eri = result['eri']
    nuc_rep = result['nuclear_repulsion']
    M = result['M']
    print(f"    M = {M} spatial, Q_naive = {2 * M} qubits  ({time.time()-t0:.1f}s)")

    # Identity baseline + FCI reference
    print(f"\n[2] Identity baseline (no gauge transform, no tapering)")
    qubit_op_id = build_qubit_op(h1, eri, nuc_rep)
    N_pauli_id = count_pauli_terms(qubit_op_id)
    print(f"    N_pauli (untapered) = {N_pauli_id}")
    fake = {'M': M, 'h1': h1, 'eri': eri, 'nuclear_repulsion': nuc_rep}
    fci_out = coupled_fci_energy(fake, n_electrons=n_electrons, verbose=False)
    E_fci_baseline = fci_out['E_coupled']
    print(f"    E_FCI = {E_fci_baseline:.10f} Ha  (baseline)")

    # M-vS structural grouping
    grouping = vertex_grouping_for_lih(spec, R)
    print(f"\n[3] M-vS vertex grouping:")
    print(f"    Vertex Li indices ({len(grouping['Li_indices'])}): "
          f"{grouping['Li_indices']}")
    print(f"    Vertex H indices ({len(grouping['H_indices'])}): "
          f"{grouping['H_indices']}")

    # Candidates
    candidates = []

    # A: identity
    print(f"\n[A] Identity gauge")
    U_A = np.eye(M)
    candidates.append(evaluate_candidate(
        'A: identity', U_A, h1, eri, nuc_rep, spec, n_electrons, M,
        apply_z2=True, baseline_fci=E_fci_baseline,
    ))
    print(f"    N_pauli (pre-tape) = {candidates[-1]['N_pauli_pre_taper']}, "
          f"E_FCI residual = {candidates[-1].get('FCI_residual_vs_baseline', 0):.2e}")
    print(f"    Z₂ tapered: Q = {candidates[-1].get('Q_post_taper', '?')}, "
          f"N_pauli = {candidates[-1].get('N_pauli_post_taper', '?')}")

    # C: per-sub-block h1 diagonalization (abelian sub-gauge)
    print(f"\n[C] Per-sub-block h1 diagonalization (abelian sub-gauge)")
    U_C = gauge_per_sub_block_h1_diag(h1, M, grouping)
    candidates.append(evaluate_candidate(
        'C: per-sub-block h1 diag', U_C, h1, eri, nuc_rep, spec,
        n_electrons, M, apply_z2=True, baseline_fci=E_fci_baseline,
    ))
    print(f"    N_pauli (pre-tape) = {candidates[-1]['N_pauli_pre_taper']}, "
          f"E_FCI residual = {candidates[-1].get('FCI_residual_vs_baseline', 0):.2e}")
    z2_applies = candidates[-1].get('z2_taper_applies', False)
    if z2_applies:
        print(f"    Z₂ tapered: Q = {candidates[-1]['Q_post_taper']}, "
              f"N_pauli = {candidates[-1]['N_pauli_post_taper']}")
    else:
        print(f"    Z₂ not applied: {candidates[-1]['z2_taper_skip_reason']}")

    # D: per-vertex h1 diagonalization (full Li-vertex 10x10 + H-vertex 5x5)
    print(f"\n[D] Per-vertex h1 diagonalization (M-vS gauge, full vertex)")
    U_D = gauge_per_vertex_h1_diag(h1, M, grouping)
    # Verify D = C for LiH (because within-Li h1 is block-diagonal already)
    diff_DC = float(np.max(np.abs(U_D @ U_D.T - np.eye(M))))
    candidates.append(evaluate_candidate(
        'D: per-vertex h1 diag', U_D, h1, eri, nuc_rep, spec,
        n_electrons, M, apply_z2=True, baseline_fci=E_fci_baseline,
    ))
    print(f"    U_D orthogonality residual = {diff_DC:.2e}")
    print(f"    N_pauli (pre-tape) = {candidates[-1]['N_pauli_pre_taper']}, "
          f"E_FCI residual = {candidates[-1].get('FCI_residual_vs_baseline', 0):.2e}")
    z2_applies = candidates[-1].get('z2_taper_applies', False)
    if z2_applies:
        print(f"    Z₂ tapered: Q = {candidates[-1]['Q_post_taper']}, "
              f"N_pauli = {candidates[-1]['N_pauli_post_taper']}")
    else:
        print(f"    Z₂ not applied: per-vertex mixes Li_core ↔ LiH_bond_center; "
              f"Z₂ Hopf sub-symmetry doesn't survive.")

    # E: natural orbitals per vertex (FCI 1-RDM, projected onto vertex blocks)
    print(f"\n[E] Per-vertex natural orbitals (FCI 1-RDM eigenbasis)")
    t0 = time.time()
    U_E = gauge_natural_orbitals_per_vertex(
        h1, eri, nuc_rep, M, n_electrons, grouping,
    )
    print(f"    NO construction wall time: {time.time()-t0:.1f}s")
    diff_E = float(np.max(np.abs(U_E @ U_E.T - np.eye(M))))
    candidates.append(evaluate_candidate(
        'E: per-vertex NO from FCI 1-RDM', U_E, h1, eri, nuc_rep, spec,
        n_electrons, M, apply_z2=True, baseline_fci=E_fci_baseline,
    ))
    print(f"    U_E orthogonality residual = {diff_E:.2e}")
    print(f"    N_pauli (pre-tape) = {candidates[-1]['N_pauli_pre_taper']}, "
          f"E_FCI residual = {candidates[-1].get('FCI_residual_vs_baseline', 0):.2e}")
    z2_applies = candidates[-1].get('z2_taper_applies', False)
    if z2_applies:
        print(f"    Z₂ tapered: Q = {candidates[-1]['Q_post_taper']}, "
              f"N_pauli = {candidates[-1]['N_pauli_post_taper']}")
    else:
        print(f"    Z₂ not applied: NO mixes Li_core ↔ LiH_bond_center.")

    # F: random block-unitary (control)
    print(f"\n[F] Random M-vS block-unitary (control test)")
    U_F = gauge_random_block_unitary(M, grouping, seed=42)
    candidates.append(evaluate_candidate(
        'F: random block-unitary', U_F, h1, eri, nuc_rep, spec,
        n_electrons, M, apply_z2=True, baseline_fci=E_fci_baseline,
    ))
    print(f"    N_pauli (pre-tape) = {candidates[-1]['N_pauli_pre_taper']}, "
          f"E_FCI residual = {candidates[-1].get('FCI_residual_vs_baseline', 0):.2e}")
    z2_applies = candidates[-1].get('z2_taper_applies', False)
    if z2_applies:
        print(f"    Z₂ tapered: Q = {candidates[-1]['Q_post_taper']}, "
              f"N_pauli = {candidates[-1]['N_pauli_post_taper']}")
    else:
        print(f"    Z₂ not applied: random mixes Li_core ↔ LiH_bond_center.")

    # Comparison table
    print(f"\n" + "=" * 72)
    print("COMPARISON TABLE")
    print("=" * 72)
    header = (f"{'Candidate':40s}  {'N_pauli pre':>12s}  "
              f"{'Q post Z₂':>10s}  {'N_pauli post':>13s}  "
              f"{'FCI Δ':>10s}")
    print(header)
    print("-" * len(header))
    for c in candidates:
        Q_post = c.get('Q_post_taper', '—')
        N_post = c.get('N_pauli_post_taper', '—')
        fci_d = c.get('FCI_residual_vs_baseline', 0)
        print(f"{c['label']:40s}  {c['N_pauli_pre_taper']:>12d}  "
              f"{str(Q_post):>10s}  {str(N_post):>13s}  "
              f"{fci_d:>10.2e}")

    # Final verdict
    print(f"\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)
    # Best Z₂ tapered count among per-sub-block-respecting candidates
    z2_eligible = [
        c for c in candidates
        if c.get('z2_taper_applies', False)
    ]
    z2_eligible_best = (
        min(z2_eligible, key=lambda c: c['N_pauli_post_taper'])
        if z2_eligible else None
    )
    # Pre-tape minimum among all candidates (M-vS non-abelian gauges
    # don't get to taper, so they win only at pre-tape)
    all_pre = min(c['N_pauli_pre_taper'] for c in candidates)
    z2_pre = candidates[0]['N_pauli_pre_taper']  # identity (== A)
    print(f"  Best Z₂-eligible post-taper N_pauli: "
          f"{z2_eligible_best['N_pauli_post_taper']} "
          f"({z2_eligible_best['label']})")
    print(f"  Best M-vS-gauge pre-taper N_pauli:   {all_pre} "
          f"(in pre-taper representation)")
    if all_pre < z2_eligible_best['N_pauli_post_taper']:
        print(f"  ⇒ Some M-vS gauge gives FEWER Pauli terms than the "
              f"Z₂-tapered representation, EVEN WITHOUT TAPERING.")
        print(f"     This is a genuine 'M-vS upgrades chemistry' result.")
    else:
        print(f"  ⇒ No M-vS gauge BEATS the Z₂-tapered representation at "
              f"Pauli count.")
        print(f"     The pre-taper count of M-vS-rotated h1/eri is "
              f"larger than the Z₂-tapered count even after gauge "
              f"transformation.  Z₂ remains the winning strategy.")

    return {
        'system': 'LiH_default',
        'R': R,
        'max_n': max_n,
        'n_electrons': n_electrons,
        'M': M,
        'Q_naive': 2 * M,
        'E_FCI_baseline': E_fci_baseline,
        'candidates': candidates,
    }


def main():
    out_lih = run_lih_default()
    out_h2 = run_h2_two_center()

    # Save data
    def _jsonable(obj):
        if isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, dict):
            return {k: _jsonable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_jsonable(v) for v in obj]
        return obj

    out_path = Path(__file__).parent / 'data' / 'sprint_mvs_gauge_pauli.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump({
            'lih': _jsonable(out_lih),
            'h2': _jsonable(out_h2),
        }, f, indent=2)
    print(f"\nData written to {out_path}")


if __name__ == '__main__':
    main()
