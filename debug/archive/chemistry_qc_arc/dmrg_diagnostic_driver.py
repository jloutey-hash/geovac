"""Sprint DMRG-diagnostic — Diagnostic driver

Question
--------
Does Paper 14 §sec:mpo_bond_rank's Theorem 3.2.A.unified (chi_k^H bounded by
universal interior profile {4, 16, 16, 9, 9, 9, 6, 3, 3, 2} and chi_k^H = 2
at sub-block boundaries) imply that DMRG on the composed GeoVac Hamiltonian
is structurally tractable for systems where the W1e correlation wall dominates?

Method
------
LiH composed at n_max = 2, R = 3.014 bohr (Q = 30, three 10-qubit sub-blocks,
zero cross-block ERIs).

Step 1 — re-verify the operator Schmidt rank profile chi_k^H at every cut of
         LiH (29 cuts). Reproduces the Sprint S2-v2 unified panel result; this
         is the verification gate.

Step 2 — exact 4-electron FCI on the composed Hamiltonian. Because cross-block
         ERIs vanish (F4), the Hamiltonian decouples as H = H_Li + H_Zeff + H_H,
         each block is 10 qubits, and ground state can be computed by per-block
         FCI in fixed particle number followed by trivial coupling of block
         particle numbers (total 4 electrons distributed across blocks).
         This gives the ground-truth FCI energy of the *composed* Hamiltonian.

Step 3 — explicit MPS representation of the per-block ground state by Schmidt
         decomposition: for each block's ground state |psi> on 10 qubits, run
         left-canonical SVD sweep to get bond dimensions at every cut. Compare
         to the operator Schmidt rank profile from Step 1.
         The OPERATOR Schmidt rank is an *upper* bound on the chi_k of the
         exact MPS ground state; we test how much smaller the state-side chi
         can be than the operator-side chi.

Step 4 — truncation sweep: truncate the per-block MPS at chi_max in
         {2, 4, 8, 16, 32}, recompute the variational energy, report
         energy gap vs exact FCI. This is the actual DMRG question:
         at fixed chi_max, can we reach FCI quality?

Step 5 — W1e structural reading: explicit accounting of which energy gap
         DMRG can close vs which gap is structurally outside DMRG's scope.

Output
------
- This driver writes its JSON to debug/data/dmrg_diagnostic.json
- The companion memo lives at debug/sprint_dmrg_diagnostic_memo.md

Convention
----------
n_qubits ordering is the spin-orbital ordering from build_composed_hamiltonian:
  block 0 (Li_core_center):    qubits  0..9   (5 spatial * 2 spin)
  block 1 (LiH_bond_center):   qubits 10..19
  block 2 (LiH_bond_partner):  qubits 20..29
JW order: qubit 2*p+s for spatial p, spin s in {0=alpha, 1=beta}.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Iterable

import numpy as np

from openfermion import QubitOperator
from openfermion.transforms import jordan_wigner
from openfermion.ops import FermionOperator

from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import lih_spec

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PATH = PROJECT_ROOT / "debug" / "data" / "dmrg_diagnostic.json"


# ---------------------------------------------------------------------------
# Operator Schmidt rank (chi_k from Pauli expansion) -- reused from S2-v2.
# ---------------------------------------------------------------------------
def operator_schmidt_rank(qop: QubitOperator, cut: int, rel_thr: float = 1e-10) -> int:
    """SVD-rank of the bipartite coefficient matrix M[P_L, P_R] at cut.

    Equivalent to the MPO bond dimension after maximally compact decomposition
    of the Hamiltonian (Paper 14 sec:mpo_bond_rank).
    """
    coef: dict[tuple, complex] = {}
    for pauli_tuple, c in qop.terms.items():
        left = tuple((q, op) for (q, op) in pauli_tuple if q < cut)
        right = tuple((q, op) for (q, op) in pauli_tuple if q >= cut)
        coef[(left, right)] = coef.get((left, right), 0.0) + c
    if not coef:
        return 0
    left_keys = sorted({k[0] for k in coef})
    right_keys = sorted({k[1] for k in coef})
    if not left_keys or not right_keys:
        return 0
    li = {k: i for i, k in enumerate(left_keys)}
    ri = {k: i for i, k in enumerate(right_keys)}
    M = np.zeros((len(left_keys), len(right_keys)), dtype=float)
    for (l, r), v in coef.items():
        M[li[l], ri[r]] = float(np.real_if_close(v).real)
    sv = np.linalg.svd(M, compute_uv=False)
    if sv.size == 0:
        return 0
    return int(np.sum(sv > rel_thr * max(sv[0], 1e-30)))


# ---------------------------------------------------------------------------
# Block extraction: build a sub-block FermionOperator and its dense matrix
# in a fixed particle-number sector.
# ---------------------------------------------------------------------------
def block_fermion_op(h1: np.ndarray, eri: np.ndarray,
                     block_orb_range: tuple[int, int],
                     orb_offset: int) -> FermionOperator:
    """Build the second-quantised block Hamiltonian, with local spin-orbital
    indexing 0..2*n-1 (alpha then beta interleaved as 2*p, 2*p+1)."""
    s, e = block_orb_range
    n = e - s
    fop = FermionOperator()
    # One-body
    for p in range(n):
        for q in range(n):
            t = h1[s + p, s + q]
            if abs(t) < 1e-14:
                continue
            for spin in (0, 1):
                fop += FermionOperator(((2 * p + spin, 1), (2 * q + spin, 0)), float(t))
    # Two-body (Mulliken / chemist's notation: (pq|rs) -> 0.5 a^dag_p a^dag_r a_s a_q)
    for p in range(n):
        for q in range(n):
            for r in range(n):
                for ss in range(n):
                    v = eri[s + p, s + q, s + r, s + ss]
                    if abs(v) < 1e-14:
                        continue
                    for sp in (0, 1):
                        for sr in (0, 1):
                            ops = ((2 * p + sp, 1), (2 * r + sr, 1),
                                   (2 * ss + sr, 0), (2 * q + sp, 0))
                            fop += FermionOperator(ops, 0.5 * float(v))
    return fop


def jw_pauli_op(fop: FermionOperator) -> QubitOperator:
    qop = jordan_wigner(fop)
    qop.compress()
    return qop


# ---------------------------------------------------------------------------
# Particle-projected diagonalisation of a 10-qubit block.
# ---------------------------------------------------------------------------
def fci_block(fop: FermionOperator, n_orb: int, n_elec_sectors: Iterable[int]):
    """Diagonalise the 2^(2n) sparse matrix and return ground state per
    (n_elec, S_z) sector (using number conservation only here -- spin can be
    inferred from occupation pattern).

    Returns dict: {n_elec: (E_min, psi)} where psi is in the n_elec sector,
    stored as full 2^(2n) vector with zeros outside the sector for downstream
    use.
    """
    from openfermion.linalg import get_sparse_operator

    Hmat = get_sparse_operator(fop, n_qubits=2 * n_orb)
    Hmat = (Hmat + Hmat.conj().T) * 0.5  # symmetrise (round-off cleanup)
    Hdense = Hmat.toarray()

    dim = 2 ** (2 * n_orb)
    # Pre-build particle-number tags for each basis state.
    pop = np.array([bin(i).count('1') for i in range(dim)], dtype=int)

    out = {}
    for ne in n_elec_sectors:
        idx = np.where(pop == ne)[0]
        if idx.size == 0:
            continue
        Hsub = Hdense[np.ix_(idx, idx)]
        Hsub = 0.5 * (Hsub + Hsub.conj().T)
        eig, vec = np.linalg.eigh(Hsub)
        psi_full = np.zeros(dim, dtype=complex)
        psi_full[idx] = vec[:, 0]
        out[ne] = {
            'E_min': float(np.real_if_close(eig[0]).real),
            'psi': psi_full,
            'sector_dim': int(idx.size),
        }
    return out


# ---------------------------------------------------------------------------
# State-side MPS construction: left-canonical SVD sweep on a 1D vector
# (10 qubits = chain of 10 sites, each d=2). Returns chi profile.
# ---------------------------------------------------------------------------
def mps_chi_profile_and_truncation(psi: np.ndarray, n_qubits: int,
                                    chi_max_list: list[int]):
    """Compute the exact bond dimensions of |psi> via successive SVDs.
    Then for each chi_max in chi_max_list, build the truncated MPS and report
    the L2 fidelity 1 - |<psi|psi_trunc>|^2 and the truncated norm-squared.

    Returns:
      exact_chi: list of length n_qubits-1, chi at each cut (Schmidt rank).
      trunc: dict[chi_max] -> dict with 'fidelity_gap' and 'truncated_norm_sq'.
    """
    dim = 2 ** n_qubits
    psi = psi / np.linalg.norm(psi)

    exact_chi: list[int] = []
    schmidt_spectra: list[list[float]] = []
    # Build chi profile by iterated SVD without storing all the bond tensors:
    # for each cut k in 1..n-1, reshape psi as (2^k) x (2^(n-k)) and SVD.
    for k in range(1, n_qubits):
        M = psi.reshape(2 ** k, 2 ** (n_qubits - k))
        sv = np.linalg.svd(M, compute_uv=False)
        thr = max(sv[0], 1e-30) * 1e-12
        chi = int(np.sum(sv > thr))
        exact_chi.append(chi)
        schmidt_spectra.append([float(s) for s in sv[:chi]])

    # Truncation sweep: for each chi_max, project at every cut.
    # Method: SVD-truncate the full state at every cut sequentially using
    # left-canonical form, then reconstruct and re-normalise.
    trunc_results: dict[int, dict] = {}
    for chi_max in chi_max_list:
        # Build truncated MPS by left-canonical SVD sweep.
        A_list = []
        psi_residual = psi.reshape(2, 2 ** (n_qubits - 1))
        bond = 1
        for k in range(1, n_qubits):
            # Current shape: (bond * 2, 2^(n-k))
            mat = psi_residual.reshape(bond * 2, 2 ** (n_qubits - k))
            U, sv, Vh = np.linalg.svd(mat, full_matrices=False)
            keep = min(chi_max, np.sum(sv > 1e-14 * max(sv[0], 1e-30)))
            U = U[:, :keep]
            sv = sv[:keep]
            Vh = Vh[:keep, :]
            A_list.append(U.reshape(bond, 2, keep))
            psi_residual = (np.diag(sv) @ Vh)
            bond = keep
        # Last site
        A_list.append(psi_residual.reshape(bond, 2, 1))
        # Reconstruct
        rec = A_list[0].reshape(2, -1)
        for kk in range(1, n_qubits):
            d_left, d_p, d_right = A_list[kk].shape
            rec = rec.reshape(-1, d_left) @ A_list[kk].reshape(d_left, d_p * d_right)
        rec = rec.reshape(2 ** n_qubits)
        rec_norm = np.linalg.norm(rec)
        ovlp = np.abs(np.vdot(psi, rec))
        trunc_results[chi_max] = {
            'fidelity_gap_1_minus_overlap_sq': float(1.0 - ovlp ** 2),
            'truncated_norm_sq': float(rec_norm ** 2),
            'achieved_bonds': [int(min(chi_max, c)) for c in exact_chi],
        }

    return exact_chi, schmidt_spectra, trunc_results


# ---------------------------------------------------------------------------
# Energy from truncated MPS: rebuild the truncated state, renormalise,
# and evaluate <psi_trunc|H|psi_trunc>.
# ---------------------------------------------------------------------------
def truncate_and_energy(psi: np.ndarray, Hdense: np.ndarray, n_qubits: int,
                         chi_max: int, particle_sector_indices: np.ndarray):
    """Variational energy of the chi_max-truncated MPS approximation of |psi>.
    The state is reprojected on the particle-number sector before normalisation,
    since MPS truncation does not exactly preserve particle number on a
    non-symmetric MPS implementation.
    """
    bond = 1
    A_list = []
    psi_residual = psi.reshape(2, 2 ** (n_qubits - 1))
    for k in range(1, n_qubits):
        mat = psi_residual.reshape(bond * 2, 2 ** (n_qubits - k))
        U, sv, Vh = np.linalg.svd(mat, full_matrices=False)
        keep = min(chi_max, int(np.sum(sv > 1e-14 * max(sv[0], 1e-30))))
        U = U[:, :keep]
        sv = sv[:keep]
        Vh = Vh[:keep, :]
        A_list.append(U.reshape(bond, 2, keep))
        psi_residual = (np.diag(sv) @ Vh)
        bond = keep
    A_list.append(psi_residual.reshape(bond, 2, 1))
    rec = A_list[0].reshape(2, -1)
    for kk in range(1, n_qubits):
        d_left, d_p, d_right = A_list[kk].shape
        rec = rec.reshape(-1, d_left) @ A_list[kk].reshape(d_left, d_p * d_right)
    rec = rec.reshape(2 ** n_qubits)
    # Project to particle-number sector
    mask = np.zeros(2 ** n_qubits, dtype=complex)
    mask[particle_sector_indices] = rec[particle_sector_indices]
    norm = np.linalg.norm(mask)
    if norm < 1e-14:
        return float('nan'), 0.0
    psi_t = mask / norm
    E = float(np.real_if_close(psi_t.conj() @ Hdense @ psi_t).real)
    return E, float(norm ** 2)


# ---------------------------------------------------------------------------
# Main run
# ---------------------------------------------------------------------------
def main():
    R = 3.014
    max_n = 2
    print(f'=== Sprint DMRG diagnostic — LiH composed, R = {R}, max_n = {max_n} ===')
    t0 = time.time()
    spec = lih_spec(R=R, max_n=max_n)
    H = build_composed_hamiltonian(spec)
    qop = H['qubit_op']
    M = H['M']
    Q = H['Q']
    n_qubits_total = Q
    E_nuc = H['nuclear_repulsion']
    N_pauli = H['N_pauli']
    h1 = H['h1']
    eri = H['eri']

    # Block ranges
    running = 0
    block_ranges_orb = []
    block_labels = []
    for b in H['blocks']:
        n = b['n_orbitals']
        block_ranges_orb.append((running, running + n))
        block_labels.append(b['label'])
        running += n
    block_ranges_qubit = [(2 * s, 2 * e) for (s, e) in block_ranges_orb]
    print(f'Q = {Q}, M = {M}, N_pauli = {N_pauli}, blocks = {len(block_ranges_orb)}')
    for lbl, (s, e) in zip(block_labels, block_ranges_qubit):
        print(f'   {lbl}: qubits [{s}..{e-1}] (n_qubits = {e - s})')

    # ----- Verification gate: cross-block ERIs/h1 vanish exactly -----
    inblock_eri = np.zeros((M, M, M, M), dtype=bool)
    for (s, e) in block_ranges_orb:
        inblock_eri[s:e, s:e, s:e, s:e] = True
    cross_eri_nnz = int(np.sum(np.abs(eri[~inblock_eri]) > 1e-14))
    inblock_h1 = np.zeros((M, M), dtype=bool)
    for (s, e) in block_ranges_orb:
        inblock_h1[s:e, s:e] = True
    cross_h1_nnz = int(np.sum(np.abs(h1[~inblock_h1]) > 1e-14))
    decoupled = (cross_eri_nnz == 0) and (cross_h1_nnz == 0)
    print(f'block-decoupling check: cross_eri_nnz={cross_eri_nnz}, '
          f'cross_h1_nnz={cross_h1_nnz} ({"OK" if decoupled else "FAIL"})')

    # ----- Step 1: operator Schmidt rank at every cut (full Hamiltonian) -----
    print('\n--- Step 1: operator Schmidt rank chi_k^H at every cut ---')
    chi_full_H = []
    for cut in range(1, n_qubits_total):
        chi = operator_schmidt_rank(qop, cut)
        chi_full_H.append(chi)
    # Check the universal interior profile within block 0 (cuts 1..10 on 10-qubit block).
    profile_block0_interior = chi_full_H[0:9]  # cuts 1..9
    print(f'block 0 interior cuts 1..9 chi profile: {profile_block0_interior}')
    print(f'expected {{4, 16, 16, 9, 9, 9, 6, 3, 3, 2}} for cuts 1..10 of a single block')
    boundary_chi_10 = chi_full_H[9]  # cut 10 (between block 0 and block 1)
    boundary_chi_20 = chi_full_H[19]  # cut 20 (between block 1 and block 2)
    print(f'boundary cuts: cut 10 chi = {boundary_chi_10}, '
          f'cut 20 chi = {boundary_chi_20} (Theorem 3.2.A.E predicts = 2)')

    # ----- Step 2: per-block FCI -----
    print('\n--- Step 2: per-block FCI ground state ---')
    block_fci = []
    block_fops = []
    block_qops = []
    for i, (s_orb, e_orb) in enumerate(block_ranges_orb):
        n_orb_b = e_orb - s_orb
        n_q_b = 2 * n_orb_b
        fop = block_fermion_op(h1, eri, (s_orb, e_orb), orb_offset=s_orb)
        qop_b = jw_pauli_op(fop)
        block_fops.append(fop)
        block_qops.append(qop_b)
        # Try all particle sectors 0..2*n_orb_b
        sectors = list(range(2 * n_orb_b + 1))
        out = fci_block(fop, n_orb_b, sectors)
        # Save block-FCI results
        E_per_sector = {ne: {'E_min': out[ne]['E_min'],
                              'sector_dim': out[ne]['sector_dim']}
                         for ne in sorted(out)}
        block_fci.append({
            'label': block_labels[i],
            'n_orbitals': n_orb_b,
            'n_qubits': n_q_b,
            'E_per_sector': E_per_sector,
            'psi_dict': out,
        })
        print(f'  block {i} ({block_labels[i]}): n_qubits = {n_q_b}, '
              f'sector E_min for ne in 0..{2 * n_orb_b}:')
        for ne in sorted(out):
            print(f'    ne={ne}: E_min = {out[ne]["E_min"]:.6f}, '
                  f'dim = {out[ne]["sector_dim"]}')

    # Total FCI energy = sum over blocks of E_min in best particle distribution.
    # Total electrons for LiH = 4 (Li has 3, H has 1).
    # Composed builder Z_eff/PK absorbs core electrons; we treat the block
    # particle numbers as free per-block and find total E_min across all
    # distributions summing to N_total.
    # However, "physical" allocation is determined by Z_eff of each block:
    # Li_core_center has Z = 3 (active 2e in 1s shell), LiH_bond_center has
    # Z = ~1.2 (active 2e bonding pair), LiH_bond_partner has Z = ~1.0 (0e
    # since PK already projected). The composed builder's PK + Z_eff fixes
    # particle counts per block.
    # For diagnostic purpose, we report total FCI over the realistic
    # 2 + 2 + 0 distribution (Li-core 1s^2 + bond pair).
    distributions = [
        (2, 2, 0),  # canonical Li 1s^2 + bond pair
        (2, 0, 2),
        (0, 2, 2),
        (4, 0, 0), (0, 4, 0), (0, 0, 4),
        (3, 1, 0), (1, 3, 0), (3, 0, 1), (1, 0, 3),
        (2, 1, 1), (1, 2, 1), (1, 1, 2),
    ]
    total_fci_table = []
    for (n0, n1, n2) in distributions:
        ok = True
        for i, ne in enumerate([n0, n1, n2]):
            if ne not in block_fci[i]['psi_dict']:
                ok = False
                break
        if not ok:
            continue
        E_total = sum(block_fci[i]['psi_dict'][ne]['E_min']
                      for i, ne in enumerate([n0, n1, n2]))
        E_total += E_nuc
        total_fci_table.append({
            'distribution_nelec_per_block': [n0, n1, n2],
            'E_total': E_total,
        })
    total_fci_table.sort(key=lambda d: d['E_total'])
    E_FCI_total = total_fci_table[0]['E_total']
    best_dist = total_fci_table[0]['distribution_nelec_per_block']
    print(f'\nLowest total FCI distribution (best 5 of {len(total_fci_table)}):')
    for entry in total_fci_table[:5]:
        print(f'  {entry["distribution_nelec_per_block"]}: '
              f'E = {entry["E_total"]:.6f} Ha')
    print(f'\nBest distribution: {best_dist}, E_FCI(total) = {E_FCI_total:.6f} Ha')
    print(f'Paper 19/20 reference: -8.055 Ha at max_n=3 (down to 0.20% from -8.07 exact).')
    print(f'Reference at max_n=2: -7.882 Ha (Paper 14 sec:onenorm; STO-3G FCI = -7.881).')

    # ----- Step 3: state-side chi profile (block 0 ground state) -----
    # Pick the per-block ground state in its canonical sector and compute
    # the operator Schmidt rank profile vs the *state* Schmidt rank.
    print('\n--- Step 3: state-side chi profile inside each block ---')
    state_chi_per_block = []
    for i, blk in enumerate(block_fci):
        ne = best_dist[i]
        if ne not in blk['psi_dict']:
            state_chi_per_block.append(None)
            continue
        psi = blk['psi_dict'][ne]['psi']
        n_q_b = blk['n_qubits']
        psi_n = psi / np.linalg.norm(psi)
        exact_chi = []
        schmidt = []
        for k in range(1, n_q_b):
            Mmat = psi_n.reshape(2 ** k, 2 ** (n_q_b - k))
            sv = np.linalg.svd(Mmat, compute_uv=False)
            thr = max(sv[0], 1e-30) * 1e-12
            chi = int(np.sum(sv > thr))
            exact_chi.append(chi)
            schmidt.append([float(s) for s in sv[:chi]])
        state_chi_per_block.append(exact_chi)
        print(f'  block {i} (ne={ne}): state chi at cuts 1..{n_q_b - 1}: '
              f'{exact_chi}')
        # Operator Schmidt rank on the same cuts of the block-internal Hamiltonian
        op_chi_within = []
        for k in range(1, n_q_b):
            op_chi_within.append(operator_schmidt_rank(block_qops[i], k))
        print(f'         operator chi^H within block: {op_chi_within}')

    # ----- Step 4: truncation sweep -----
    print('\n--- Step 4: bond-dim truncation sweep ---')
    chi_max_list = [1, 2, 4, 8, 16, 32]
    truncation_results = []
    from openfermion.linalg import get_sparse_operator
    for i, blk in enumerate(block_fci):
        ne = best_dist[i]
        if ne not in blk['psi_dict']:
            continue
        psi = blk['psi_dict'][ne]['psi']
        n_q_b = blk['n_qubits']
        psi_n = psi / np.linalg.norm(psi)
        Hmat = get_sparse_operator(block_fops[i], n_qubits=n_q_b).toarray()
        Hmat = 0.5 * (Hmat + Hmat.conj().T)
        pop = np.array([bin(j).count('1') for j in range(2 ** n_q_b)], dtype=int)
        sector_idx = np.where(pop == ne)[0]
        E_FCI = blk['psi_dict'][ne]['E_min']
        per_chi = {}
        for chi_max in chi_max_list:
            E, norm_sq = truncate_and_energy(psi_n, Hmat, n_q_b, chi_max, sector_idx)
            per_chi[chi_max] = {
                'E_variational': E,
                'projected_norm_sq': norm_sq,
                'energy_gap_vs_FCI_Ha': E - E_FCI,
            }
            print(f'  block {i} (ne={ne}, FCI={E_FCI:.6f}): chi_max={chi_max}, '
                  f'E = {E:.6f}, gap = {E - E_FCI:+.2e} Ha, '
                  f'sector_norm^2 = {norm_sq:.4f}')
        truncation_results.append({
            'block_idx': i,
            'block_label': blk['label'],
            'n_elec': ne,
            'E_FCI': E_FCI,
            'per_chi': per_chi,
        })

    # ----- Save -----
    t_elapsed = time.time() - t0
    # Strip psi vectors from JSON (too large)
    block_fci_to_save = []
    for blk in block_fci:
        block_fci_to_save.append({
            'label': blk['label'],
            'n_orbitals': blk['n_orbitals'],
            'n_qubits': blk['n_qubits'],
            'E_per_sector': blk['E_per_sector'],
        })

    out = {
        'meta': {
            'sprint': 'DMRG-diagnostic',
            'date_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'molecule': 'LiH composed',
            'R_bohr': R,
            'max_n': max_n,
            'Q': Q,
            'M': M,
            'N_pauli': N_pauli,
            'E_nuclear_repulsion': E_nuc,
            'block_ranges_orbital': block_ranges_orb,
            'block_ranges_qubit': block_ranges_qubit,
            'block_labels': block_labels,
            'wall_time_s': t_elapsed,
        },
        'verification_gate': {
            'cross_block_eri_nnz': cross_eri_nnz,
            'cross_block_h1_nnz': cross_h1_nnz,
            'block_decoupled_exactly': decoupled,
        },
        'step1_operator_chi_full_H': {
            'cuts': list(range(1, n_qubits_total)),
            'chi': chi_full_H,
            'block_0_interior_profile_cuts_1_to_9': profile_block0_interior,
            'expected_universal_profile': [4, 16, 16, 9, 9, 9, 6, 3, 3, 2],
            'boundary_cut_10_chi': boundary_chi_10,
            'boundary_cut_20_chi': boundary_chi_20,
            'theorem_3_2_A_E_expected': 2,
        },
        'step2_per_block_fci': block_fci_to_save,
        'step2_total_fci': {
            'distributions_tested': [d['distribution_nelec_per_block'] for d in total_fci_table],
            'distribution_energies_Ha': [d['E_total'] for d in total_fci_table],
            'lowest_energy_Ha': E_FCI_total,
            'lowest_distribution': best_dist,
            'paper_14_reference_n2': -7.882,
            'paper_19_20_reference_n3': -8.055,
        },
        'step3_state_chi_per_block': state_chi_per_block,
        'step4_truncation_sweep': truncation_results,
    }
    DATA_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(DATA_PATH, 'w') as f:
        json.dump(out, f, indent=2, default=lambda o: float(o) if hasattr(o, 'real') else str(o))
    print(f'\nWrote {DATA_PATH} ({DATA_PATH.stat().st_size / 1024:.1f} KB)')
    print(f'Total wall time: {t_elapsed:.1f} s')
    return out


if __name__ == '__main__':
    main()
