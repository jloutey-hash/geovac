"""Sprint F6 Step 2 - NaH 2-electron FCI at max_n=4 with full F3 stack.

Runs the 2-electron FCI at NaH max_n=4 (Q=120, M=60) with the full
F3 stack at R = 3.5 bohr and R = 10.0 bohr:
  - W1c screening (screened_cross_center=True)
  - Multi-zeta basis (multi_zeta_basis=True)
  - Cross-block h1 architectural extension (cross_block_h1=True)

Reports:
  - E(R=3.5), E(R=10), D_e^max_n=4 = E(10) - E(3.5)
  - Natural orbital occupations
  - Dominant natural orbital Na/H amplitude split
  - h1 lowest eigenvalues (bonding/antibonding ordering)

Gate (per F6 sprint plan):
  - D_e in [0.0375, 0.150] Ha AND bonding signature -> WALL-CLOSES
  - D_e in (0.150, 1.0] Ha -> PARTIAL (physical magnitude class)
  - D_e in (1.0, 4.37] Ha -> PARTIAL-IMPROVEMENT (improved vs max_n=3 baseline 4.37 Ha)
  - D_e <= 0 OR D_e >= 4.37 Ha -> WALL-PERSISTS

Comparison to F3 baseline:
  - F3 max_n=2 W1c+mz+xblockh1: D_e^F3 = +4.37 Ha (well depth at R_min=2 bohr)
    Natural occs: [1.9991, 0.0007] bonding 50/50 Na/H
  - F1 max_n=3 W1c+mz alone: D_e^PES = +4.37 Ha spurious binding
"""

from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh

sys.stdout.reconfigure(line_buffering=True)

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import nah_spec


def build_fci_with_eigenvector(h1, eri, n_electrons, nuclear_repulsion, M,
                                k_states: int = 6):
    """Build 2-electron singlet FCI Hamiltonian, return lowest eigenpairs."""
    n_up = n_electrons // 2
    n_down = n_electrons // 2

    alpha_strings = list(combinations(range(M), n_up))
    beta_strings = list(combinations(range(M), n_down))
    n_alpha = len(alpha_strings)
    n_beta = len(beta_strings)
    n_det = n_alpha * n_beta

    alpha_idx = {s: i for i, s in enumerate(alpha_strings)}
    beta_idx = {s: i for i, s in enumerate(beta_strings)}

    print(f"  FCI dim: {n_det} ({n_alpha} alpha x {n_beta} beta)")

    H_fci = lil_matrix((n_det, n_det))

    def det_index(a_idx, b_idx):
        return a_idx * n_beta + b_idx

    def exc_phase(occ, p, r):
        lo, hi = min(p, r), max(p, r)
        return (-1) ** sum(1 for o in occ if lo < o < hi)

    # Diagonal
    for ai, alpha in enumerate(alpha_strings):
        for bi, beta in enumerate(beta_strings):
            I = det_index(ai, bi)
            E_diag = nuclear_repulsion
            for p in alpha:
                E_diag += h1[p, p]
            for p in beta:
                E_diag += h1[p, p]
            for i_idx in range(n_up):
                for j_idx in range(i_idx + 1, n_up):
                    p, q = alpha[i_idx], alpha[j_idx]
                    E_diag += eri[p, p, q, q] - eri[p, q, q, p]
            for i_idx in range(n_down):
                for j_idx in range(i_idx + 1, n_down):
                    p, q = beta[i_idx], beta[j_idx]
                    E_diag += eri[p, p, q, q] - eri[p, q, q, p]
            for p in alpha:
                for q in beta:
                    E_diag += eri[p, p, q, q]
            H_fci[I, I] = E_diag

    # Single excitations alpha
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
                phase = exc_phase(alpha, p, r)
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
                        H_fci[I, J] += val

    # Single excitations beta
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
                phase = exc_phase(beta, p, r)
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
                        H_fci[I, J] += val

    # Alpha-beta double exc
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
                phase_a = exc_phase(alpha, p_a, r_a)
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
                            phase_b = exc_phase(beta, p_b, r_b)
                            val = (phase_a * phase_b *
                                   eri[r_a, p_a, r_b, p_b])
                            if abs(val) > 1e-14:
                                I = det_index(ai, bi)
                                J = det_index(ai_new, bi_new)
                                H_fci[I, J] += val

    H_fci = H_fci.tocsr()
    k = min(k_states, n_det - 1)
    print(f"  Diagonalizing (eigsh, k={k}, SA)...")
    t_diag = time.perf_counter()
    E_vals, V_vecs = eigsh(H_fci, k=k, which='SA')
    elapsed_diag = time.perf_counter() - t_diag
    print(f"  Diag complete ({elapsed_diag:.1f}s)")
    order = np.argsort(E_vals)
    E_vals = E_vals[order]
    V_vecs = V_vecs[:, order]

    return E_vals, V_vecs, alpha_strings, beta_strings


def compute_1rdm(psi, alpha_strings, beta_strings, M):
    """Compute spatial 1-RDM from FCI ground state."""
    n_alpha = len(alpha_strings)
    n_beta = len(beta_strings)
    rdm = np.zeros((M, M), dtype=float)

    def exc_phase(occ, p, r):
        lo, hi = min(p, r), max(p, r)
        return (-1) ** sum(1 for o in occ if lo < o < hi)

    coef = psi.reshape(n_alpha, n_beta)

    alpha_idx = {s: i for i, s in enumerate(alpha_strings)}
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        for p in alpha:
            for bi in range(n_beta):
                rdm[p, p] += coef[ai, bi] ** 2
        for p in alpha:
            for r in range(M):
                if r in alpha_set or r == p:
                    continue
                new_alpha = tuple(sorted((alpha_set - {p}) | {r}))
                if new_alpha not in alpha_idx:
                    continue
                ai_new = alpha_idx[new_alpha]
                phase = exc_phase(alpha, p, r)
                for bi in range(n_beta):
                    rdm[r, p] += phase * coef[ai_new, bi] * coef[ai, bi]

    beta_idx = {s: i for i, s in enumerate(beta_strings)}
    for bi, beta in enumerate(beta_strings):
        beta_set = set(beta)
        for p in beta:
            for ai in range(n_alpha):
                rdm[p, p] += coef[ai, bi] ** 2
        for p in beta:
            for r in range(M):
                if r in beta_set or r == p:
                    continue
                new_beta = tuple(sorted((beta_set - {p}) | {r}))
                if new_beta not in beta_idx:
                    continue
                bi_new = beta_idx[new_beta]
                phase = exc_phase(beta, p, r)
                for ai in range(n_alpha):
                    rdm[p, r] += phase * coef[ai, bi_new] * coef[ai, bi]

    rdm = 0.5 * (rdm + rdm.T)
    return rdm


def diagnose_h1_spectrum(h1, M_total):
    """Eigendecompose h1, report lowest eigenvalues and Na/H character.

    At NaH max_n=4, M=60, split is M//2 = 30 (Na) + 30 (H).
    """
    M_per_side = M_total // 2
    eigs, vecs = np.linalg.eigh(h1)
    out = {
        'eigenvalues': eigs[:8].tolist(),
        'lowest_4_decomposition': [],
    }
    for k in range(min(4, len(eigs))):
        v = vecs[:, k]
        amp_Na = float(np.sum(v[:M_per_side] ** 2))
        amp_H = float(np.sum(v[M_per_side:] ** 2))
        # Na 3s and H 1s coefficients (first orbital in each block)
        out['lowest_4_decomposition'].append({
            'eigenvalue_Ha': float(eigs[k]),
            'Na_amplitude_squared': amp_Na,
            'H_amplitude_squared': amp_H,
            'na3s_coeff': float(v[0]),  # First Na orbital (block_n=1, l=0)
            'h1s_coeff': float(v[M_per_side]),  # First H orbital
            'character': 'bonding' if (v[0] * v[M_per_side] > 0) else 'antibonding',
        })
    return out


def run_arch(spec, R, label, **balanced_kwargs):
    """Build hamiltonian, run FCI, return diagnostic."""
    n_e = sum(b.n_electrons for b in spec.blocks)
    t0 = time.perf_counter()
    print(f"  Building at R={R:.1f}...")
    ham = build_balanced_hamiltonian(spec, R=R, verbose=False, **balanced_kwargs)
    t_build = time.perf_counter() - t0
    print(f"  Build complete ({t_build:.1f}s = {t_build/60:.2f}min)")

    M = ham['M']
    h1 = ham['h1']
    eri = ham['eri']
    nuc_rep = ham['nuclear_repulsion']

    print(f"  M = {M}, n_electrons = {n_e}, nuclear_repulsion = {nuc_rep:.4f} Ha")

    E_vals, V_vecs, alpha_strings, beta_strings = build_fci_with_eigenvector(
        h1, eri, n_e, nuc_rep, M,
    )
    psi0 = V_vecs[:, 0]
    E_gs = float(E_vals[0])
    print(f"  E_gs = {E_gs:+.6f} Ha")
    print(f"  Computing 1-RDM...")
    t_rdm = time.perf_counter()
    rdm = compute_1rdm(psi0, alpha_strings, beta_strings, M)
    elapsed_rdm = time.perf_counter() - t_rdm
    print(f"  1-RDM complete ({elapsed_rdm:.1f}s)")

    nat_occs_v, nat_vecs_v = np.linalg.eigh(rdm)
    order = np.argsort(nat_occs_v)[::-1]
    nat_occs_v = nat_occs_v[order]
    nat_vecs_v = nat_vecs_v[:, order]
    M_per_side = M // 2
    dom = nat_vecs_v[:, 0]
    amp_Na = float(np.sum(dom[:M_per_side] ** 2))
    amp_H = float(np.sum(dom[M_per_side:] ** 2))

    # Top-3 components of dominant NO by |amp|
    abs_amps = np.abs(dom)
    top3_idx = np.argsort(abs_amps)[-3:][::-1]
    top3_components = [(int(i), float(dom[i])) for i in top3_idx]

    h1_spec = diagnose_h1_spectrum(h1, M)

    elapsed = time.perf_counter() - t0
    print(f"  Total run_arch elapsed: {elapsed:.1f}s = {elapsed/60:.2f}min")
    return {
        'label': label,
        'R': R,
        'E_gs': E_gs,
        'natural_occupations': nat_occs_v[:6].tolist(),
        'dominant_NO': {
            'occupation': float(nat_occs_v[0]),
            'Na_amp_squared': amp_Na,
            'H_amp_squared': amp_H,
            'na3s_coeff': float(dom[0]),
            'h1s_coeff': float(dom[M_per_side]),
            'character': 'bonding' if (dom[0] * dom[M_per_side] > 0) else 'antibonding',
            'top3_components': top3_components,
        },
        'rdm_diag_top10': sorted(np.diag(rdm).tolist(), reverse=True)[:10],
        'h1_spectrum_diag': h1_spec,
        'cross_block_h1_info': ham.get('cross_block_h1_info', {}),
        'build_time_s': t_build,
        'elapsed_s': elapsed,
    }


def main():
    out_path = Path('debug/data/sprint_f6_step2_fci.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    spec = nah_spec(max_n=4)
    print("=" * 70)
    print(f"Sprint F6 Step 2 - NaH 2e FCI at max_n=4 with full F3 stack")
    print(f"  M=60, Q=120, fci_dim=3600, 2e singlet")
    print("=" * 70)

    # Only run the full F3 stack (the architecture under test)
    architectures = [
        ('W1c+mz+xblockh1',
            dict(screened_cross_center=True, multi_zeta_basis=True,
                 cross_block_h1=True)),
    ]

    results = {}
    for label, kwargs in architectures:
        print(f"\n--- {label} ---")
        for R in [3.5, 10.0]:
            print(f"\n  ==== R = {R} bohr ====")
            r = run_arch(spec, R, label, **kwargs)
            print(f"  R={R:.1f}: E_gs={r['E_gs']:+.6f}  "
                  f"NO_dom={r['dominant_NO']['occupation']:.4f} "
                  f"({r['dominant_NO']['character']}) "
                  f"Na/H={r['dominant_NO']['Na_amp_squared']:.3f}/"
                  f"{r['dominant_NO']['H_amp_squared']:.3f}")
            print(f"          natural occs (top 6): {[f'{o:.4f}' for o in r['natural_occupations']]}")
            print(f"          h1 lowest 4 eigvals: "
                  f"{[f'{e:+.4f}' for e in r['h1_spectrum_diag']['eigenvalues'][:4]]}")
            print(f"          lowest h1 eigvec: "
                  f"Na/H = {r['h1_spectrum_diag']['lowest_4_decomposition'][0]['Na_amplitude_squared']:.3f}/"
                  f"{r['h1_spectrum_diag']['lowest_4_decomposition'][0]['H_amplitude_squared']:.3f} "
                  f"({r['h1_spectrum_diag']['lowest_4_decomposition'][0]['character']})")
            results.setdefault(label, []).append(r)

    # D_e
    label = 'W1c+mz+xblockh1'
    Es = {r['R']: r['E_gs'] for r in results[label]}
    D_e_max_n_4 = Es[10.0] - Es[3.5]
    D_e_F3 = 4.374  # max_n=2 baseline from F3 §3
    D_e_exp = 0.075  # experimental NaH

    print(f"\n{'=' * 70}")
    print(f"F6 RESULTS (max_n=4):")
    print(f"  E(R=3.5)  = {Es[3.5]:+.6f} Ha")
    print(f"  E(R=10.0) = {Es[10.0]:+.6f} Ha")
    print(f"  D_e^max_n=4 (2-point) = {D_e_max_n_4:+.6f} Ha")
    print(f"  D_e^max_n=2 baseline (F3) = {D_e_F3:+.4f} Ha")
    print(f"  D_e experimental = {D_e_exp:.4f} Ha")
    if D_e_F3 > 0:
        wall_closure_fraction = (D_e_F3 - D_e_max_n_4) / D_e_F3 * 100
        print(f"  Wall closure fraction = {wall_closure_fraction:+.2f}%")
    print(f"{'=' * 70}")

    # F6 verdict
    naturals = results[label][0]['natural_occupations'][:2]
    no_dom_char = results[label][0]['dominant_NO']['character']
    na_amp = results[label][0]['dominant_NO']['Na_amp_squared']
    h_amp = results[label][0]['dominant_NO']['H_amp_squared']
    bonding_signature = (
        no_dom_char == 'bonding' and na_amp > 0.1 and h_amp > 0.1
    )
    magnitude_in_window = 0.5 * D_e_exp <= D_e_max_n_4 <= 2.0 * D_e_exp

    if D_e_max_n_4 > 0 and bonding_signature and magnitude_in_window:
        verdict = "WALL_CLOSES"
    elif 0.15 < D_e_max_n_4 <= 1.0:
        verdict = "PARTIAL"
    elif 1.0 < D_e_max_n_4 < D_e_F3:
        verdict = "PARTIAL_IMPROVEMENT"
    elif D_e_max_n_4 >= D_e_F3 or D_e_max_n_4 <= 0:
        verdict = "WALL_PERSISTS"
    else:
        verdict = "INTERMEDIATE"  # could be PARTIAL or borderline

    print(f"F6 VERDICT: {verdict}")
    print(f"  Bonding signature: {bonding_signature}")
    print(f"  Magnitude in window [0.0375, 0.15]: {magnitude_in_window}")
    print(f"{'=' * 70}")

    out = {
        'sprint': 'F6 Step 2 - NaH 2e FCI at max_n=4',
        'date': '2026-05-23',
        'spec': 'nah_spec(max_n=4), Q=120, 2-electron singlet FCI, M=60, fci_dim=3600',
        'architectures': {k: results[k] for k in results},
        'D_e_max_n_4_Ha': D_e_max_n_4,
        'D_e_max_n_2_baseline_F3_Ha': D_e_F3,
        'D_e_experimental_Ha': D_e_exp,
        'wall_closure_fraction_pct': (D_e_F3 - D_e_max_n_4) / D_e_F3 * 100 if D_e_F3 > 0 else None,
        'verdict': verdict,
        'bonding_signature': bonding_signature,
        'magnitude_in_window': magnitude_in_window,
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote: {out_path}")


if __name__ == "__main__":
    main()
