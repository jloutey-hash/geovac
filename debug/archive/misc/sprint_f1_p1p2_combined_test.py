"""Sprint F1 Phase 2: combined W1c x multi-zeta FCI at NaH max_n=2.

Tests whether the multi-zeta basis substitution affects the FCI ground-state
energy and Na 3s natural occupation when COMBINED with W1c screened cross-V_ne.

The alpha-PES Step 2 finding was: at NaH max_n=2 with BARE cross-V_ne,
multi-zeta is bit-zero on the FCI (because Na Z=11 over-attracts H 1s,
collapsing all FCI weight onto H-side orbitals). The hypothesis under test:
with W1c screening reducing the bare-cross-V_ne over-attraction, the FCI
might engage Na 3s, and the multi-zeta substitution would then become
load-bearing.

Phase 1 (this commit) removed the NotImplementedError that blocked the
combined test. This Phase 2 driver:
  1. Runs four-way FCI at NaH max_n=2 R=3.5 and R=10.0:
     - bare
     - bare + multi-zeta
     - W1c (screened, no multi-zeta)
     - **combined**: W1c + multi-zeta
  2. Reports binding energy D_e = E(R=10) - E(R=3.5) for each
  3. Computes FCI 1-RDM and Na 3s natural occupation for the combined case
  4. Diagnoses the h1 eigenspectrum with W1c on (where does Na 3s sit?)
  5. If combined gives binding emergence, run mini-PES at R in {3.0, 3.5, 4.0}

The Na sub-block lives at orbital indices 0..4 (block_n=1: l=0; block_n=2:
l=0, l=1 m=-1,0,+1).  Index 0 = (1,0,0) = physical Na 3s.

Driver: write all numerical results to debug/data/sprint_f1_p1p2_results.json.
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
    """FCI ground state via sector-restricted diagonalization, returning
    the eigenvector (for 1-RDM analysis).

    Reuses the same algorithm as coupled_composition.coupled_fci_energy but
    keeps the eigenvector.  Spin-restricted: n_up = n_down = N/2.
    """
    n_up = n_electrons // 2
    n_down = n_electrons // 2

    alpha_strings = list(combinations(range(M), n_up))
    beta_strings = list(combinations(range(M), n_down))
    n_alpha = len(alpha_strings)
    n_beta = len(beta_strings)
    n_det = n_alpha * n_beta

    alpha_idx = {s: i for i, s in enumerate(alpha_strings)}
    beta_idx = {s: i for i, s in enumerate(beta_strings)}

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

    # Single excitations: alpha
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

    # Single excitations: beta
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

    # Alpha-beta double excitations (only contributing for n_e>=4); for n_e=2
    # there's no alpha-beta double exc beyond singles.
    # For n_e=2 with n_up=1, no same-spin doubles possible either.
    if n_up >= 2 or n_down >= 2:
        # Same-spin doubles (need n_up >= 2 and/or n_down >= 2)
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
                            phase = _double_exc_phase(alpha, p, q, r, s)
                            val = phase * (eri[r, p, s, q] - eri[r, q, s, p])
                            if abs(val) > 1e-14:
                                for bi in range(n_beta):
                                    I = det_index(ai, bi)
                                    J = det_index(ai_new, bi)
                                    H_fci[I, J] += val
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
                            phase = _double_exc_phase(beta, p, q, r, s)
                            val = phase * (eri[r, p, s, q] - eri[r, q, s, p])
                            if abs(val) > 1e-14:
                                for ai in range(n_alpha):
                                    I = det_index(ai, bi)
                                    J = det_index(ai, bi_new)
                                    H_fci[I, J] += val
    # Alpha-beta double exc (always contributes for n_up >= 1 and n_down >= 1)
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
    E_vals, V_vecs = eigsh(H_fci, k=k, which='SA')
    order = np.argsort(E_vals)
    E_vals = E_vals[order]
    V_vecs = V_vecs[:, order]

    return E_vals, V_vecs, alpha_strings, beta_strings


def _double_exc_phase(occ, p, q, r, s):
    """Phase for double excitation {p,q} -> {r,s} (p<q, r<s)."""
    occ_set = set(occ)
    phase_remove_p = sum(1 for o in occ if o < p)
    remaining_after_p = sorted(occ_set - {p})
    phase_remove_q = sum(1 for o in remaining_after_p if o < q)
    remaining = sorted(occ_set - {p, q})
    phase_add_r = sum(1 for o in remaining if o < r)
    remaining_plus_r = sorted(set(remaining) | {r})
    phase_add_s = sum(1 for o in remaining_plus_r if o < s)
    total = phase_remove_p + phase_remove_q + phase_add_r + phase_add_s
    return (-1) ** total


def compute_1rdm(psi, alpha_strings, beta_strings, M):
    """Compute spatial 1-RDM from FCI ground state.

    1-RDM[p, q] = <psi| a^dag_p a_q |psi> summed over alpha + beta spins.

    Indices: psi is flat in (alpha_idx, beta_idx) lex order with
        I = ai * n_beta + bi.
    """
    n_alpha = len(alpha_strings)
    n_beta = len(beta_strings)
    rdm = np.zeros((M, M), dtype=float)

    def exc_phase(occ, p, r):
        lo, hi = min(p, r), max(p, r)
        return (-1) ** sum(1 for o in occ if lo < o < hi)

    # Reshape psi to (n_alpha, n_beta) coefficient matrix
    coef = psi.reshape(n_alpha, n_beta)

    # Alpha contribution
    alpha_idx = {s: i for i, s in enumerate(alpha_strings)}
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        # Diagonal (p=q=occupied)
        for p in alpha:
            for bi in range(n_beta):
                rdm[p, p] += coef[ai, bi] ** 2
        # Off-diagonal: a^dag_r a_p with p in alpha, r not in alpha
        for p in alpha:
            for r in range(M):
                if r in alpha_set or r == p:
                    continue
                new_alpha = tuple(sorted((alpha_set - {p}) | {r}))
                if new_alpha not in alpha_idx:
                    continue
                ai_new = alpha_idx[new_alpha]
                phase = exc_phase(alpha, p, r)
                # <new_alpha_ai_new, bi| a^dag_r a_p |alpha_ai, bi> = phase
                for bi in range(n_beta):
                    rdm[r, p] += phase * coef[ai_new, bi] * coef[ai, bi]

    # Beta contribution
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

    # Symmetrize and return
    rdm = 0.5 * (rdm + rdm.T)
    return rdm


def diagnose_h1_spectrum(ham, label, M):
    """Diagnose h1 eigenspectrum: where does Na 3s sit?"""
    h1 = ham['h1']
    eigvals, eigvecs = np.linalg.eigh(h1)
    out = {
        'label': label,
        'lowest_5_eigvals': eigvals[:5].tolist(),
        'eigvals_full': eigvals.tolist(),
        'na_3s_amplitude_in_low5': [
            float(eigvecs[0, i] ** 2) for i in range(5)
        ],
        'h_1s_amplitude_in_low5': [
            float(eigvecs[5, i] ** 2) for i in range(min(5, M))
        ],
    }
    return out


def run_arch(spec, R, screened, multi_zeta, label, verbose=False):
    """Build hamiltonian, run FCI, return eigenvector + Na 3s occupation."""
    n_e = sum(b.n_electrons for b in spec.blocks)
    t0 = time.perf_counter()
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=4000, L_max=4,
        screened_cross_center=screened,
        multi_zeta_basis=multi_zeta,
        verbose=False,
    )
    M = ham['M']
    h1 = ham['h1']
    eri = ham['eri']
    nuc_rep = ham['nuclear_repulsion']

    E_vals, V_vecs, alpha_strings, beta_strings = build_fci_with_eigenvector(
        h1, eri, n_e, nuc_rep, M
    )
    psi0 = V_vecs[:, 0]
    E_gs = float(E_vals[0])
    rdm = compute_1rdm(psi0, alpha_strings, beta_strings, M)
    # Natural occupation numbers (eigenvalues of 1-RDM)
    nat_occs = np.sort(np.linalg.eigvalsh(rdm))[::-1]
    elapsed = time.perf_counter() - t0

    # Na 3s diagonal occupation = rdm[0, 0]
    # Total electrons = trace
    total_e = float(np.trace(rdm))

    result = {
        'label': label,
        'R': float(R),
        'screened': bool(screened),
        'multi_zeta': bool(multi_zeta),
        'E_gs': E_gs,
        'E_excited_5': E_vals[1:6].tolist() if len(E_vals) > 1 else [],
        'na_3s_diag_occ': float(rdm[0, 0]),
        'h_1s_diag_occ': float(rdm[5, 5]) if M >= 6 else None,
        'total_electrons_from_rdm': total_e,
        'natural_occupations': nat_occs.tolist(),
        'natural_occ_above_001': int(np.sum(nat_occs > 0.01)),
        'natural_occ_above_010': int(np.sum(nat_occs > 0.10)),
        'natural_occ_above_100': int(np.sum(nat_occs > 1.0)),
        'rdm_diag': np.diag(rdm).tolist(),
        'h1_diag_full': h1.diagonal().tolist(),
        'tr_h1_cross_vne': float(np.trace(ham['h1_cross_vne'])),
        'h1_cross_vne_diag': ham['h1_cross_vne'].diagonal().tolist(),
        'elapsed_s': elapsed,
        'M': int(M),
        'n_electrons': int(n_e),
    }
    # Also diagnose h1 eigenspectrum
    diag = diagnose_h1_spectrum(ham, label, M)
    result['h1_spectrum_diag'] = diag

    if verbose:
        print(f"  {label:>25}: E = {E_gs:.6f} Ha  "
              f"Na3s_occ = {result['na_3s_diag_occ']:.4f}  "
              f"H1s_occ = {result['h_1s_diag_occ']:.4f}  "
              f"[{elapsed:.1f}s]", flush=True)
    return result


def main():
    out_path = Path('debug/data/sprint_f1_p1p2_results.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Sprint F1 Phase 2: Combined W1c x multi-zeta FCI at NaH max_n=2")
    print("=" * 70)
    print()

    # Step 1: two-point FCI at R=3.5 and R=10.0 across four architectures
    R_eq = 3.5
    R_diss = 10.0
    print(f"Step 1: Two-point FCI at R_eq={R_eq} and R_diss={R_diss}")
    print(f"  (Na 3s is at orbital index 0; H 1s is at index 5)")
    print()

    spec = nah_spec(max_n=2)

    # The four architectures
    archs = [
        ('PK+SV-class baseline (bare)',        False, False),
        ('multi-zeta alone (bare + mz)',       False, True),
        ('W1c alone (screened, no mz)',        True,  False),
        ('COMBINED (W1c + multi-zeta)',        True,  True),
    ]

    # Run all four at R_eq
    print(f"--- R = {R_eq} bohr ---")
    results_eq = {}
    for label, screened, mz in archs:
        r = run_arch(spec, R_eq, screened, mz, label, verbose=True)
        results_eq[label] = r
    print()

    # Run all four at R_diss
    print(f"--- R = {R_diss} bohr ---")
    results_diss = {}
    for label, screened, mz in archs:
        r = run_arch(spec, R_diss, screened, mz, label, verbose=True)
        results_diss[label] = r
    print()

    # Compute D_e for each
    print("--- Binding energy D_e = E(R=10) - E(R=3.5) ---")
    de_table = {}
    for label, _, _ in archs:
        E_eq = results_eq[label]['E_gs']
        E_diss = results_diss[label]['E_gs']
        D_e = E_diss - E_eq
        de_table[label] = {
            'E_eq': E_eq, 'E_diss': E_diss, 'D_e': D_e,
        }
        sign = '+' if D_e > 0 else '-'
        print(f"  {label:>30}: E_eq = {E_eq:11.6f}, E_diss = {E_diss:11.6f}, "
              f"D_e = {sign}{abs(D_e):.6f} Ha")
    print()

    # Na 3s occupation deep-dive
    print("--- Na 3s diagonal occupation across architectures (at R_eq) ---")
    print(f"  (max possible: 2.0 = fully occupied; FCI ground state expected "
          f"in [0, 2])")
    print()
    for label, _, _ in archs:
        r = results_eq[label]
        print(f"  {label:>30}: "
              f"Na 3s = {r['na_3s_diag_occ']:.5f}, "
              f"H 1s = {r['h_1s_diag_occ']:.5f}, "
              f"tr(RDM) = {r['total_electrons_from_rdm']:.4f}")
    print()

    print("--- Natural occupations (sorted, top 6) at R_eq ---")
    for label, _, _ in archs:
        r = results_eq[label]
        occs = r['natural_occupations'][:6]
        print(f"  {label:>30}: {['%.4f' % o for o in occs]}")
    print()

    print("--- h1 eigenspectrum at R_eq: where does Na 3s sit? ---")
    print(f"  (Na 3s amplitude squared in the lowest-5 h1 eigenvectors)")
    for label, _, _ in archs:
        r = results_eq[label]
        diag = r['h1_spectrum_diag']
        print(f"  {label:>30}: low5 eigvals = "
              f"{['%.3f' % e for e in diag['lowest_5_eigvals']]}")
        print(f"  {'':>30}  Na3s|^2 in low5 = "
              f"{['%.4f' % a for a in diag['na_3s_amplitude_in_low5']]}")
        print(f"  {'':>30}  H1s|^2 in low5  = "
              f"{['%.4f' % a for a in diag['h_1s_amplitude_in_low5']]}")
    print()

    # Decide verdict based on Na 3s occupation in COMBINED case at R_eq
    combined_label = 'COMBINED (W1c + multi-zeta)'
    na_3s_combined = results_eq[combined_label]['na_3s_diag_occ']
    de_combined = de_table[combined_label]['D_e']

    print("=" * 70)
    print("VERDICT")
    print("=" * 70)
    print(f"  Na 3s diagonal occupation (combined, R_eq): {na_3s_combined:.5f}")
    print(f"  D_e^combined: {de_combined:+.6f} Ha")
    print()

    if na_3s_combined > 0.1 and de_combined > 0:
        verdict = 'W1C-LAYER3-EVAPORATES'
        print(f"  VERDICT: {verdict}")
        print(f"  Combined architecture engages Na 3s significantly "
              f"({na_3s_combined:.3f} > 0.1)")
        print(f"  AND produces a binding minimum (D_e = "
              f"{de_combined:.4f} > 0).")
        print(f"  W1c Layer 3 was a mutual-exclusivity artifact.")
    elif na_3s_combined > 0.1 and de_combined <= 0:
        verdict = 'PARTIAL-CLOSURE-AT-MAX_N=2'
        print(f"  VERDICT: {verdict}")
        print(f"  Combined architecture engages Na 3s significantly "
              f"({na_3s_combined:.3f} > 0.1)")
        print(f"  BUT does NOT produce binding (D_e = {de_combined:.4f}).")
        print(f"  FCI engages but cross-V_ne kernel still has residual.")
        print(f"  Next sprint: F2 cross-V_ne kernel-shape substitution.")
    else:
        verdict = 'LAYER3-CORROBORATED'
        print(f"  VERDICT: {verdict}")
        print(f"  Combined architecture does NOT engage Na 3s "
              f"({na_3s_combined:.3f} < 0.1).")
        print(f"  W1c Layer 3 finding holds at max_n=2.")
        print(f"  Next sprint: max_n=3 basis-size test, or P2 kernel-shape.")
    print()

    # Step 4 (only if verdict is W1C-LAYER3-EVAPORATES): mini-PES
    if verdict == 'W1C-LAYER3-EVAPORATES':
        print("--- Step 4: Mini-PES at R in {3.0, 3.5, 4.0} ---")
        mini_pes = {}
        for R in [3.0, 3.5, 4.0]:
            t0 = time.perf_counter()
            r = run_arch(spec, R, True, True, f'combined R={R}', verbose=True)
            mini_pes[R] = r
        print()
        print("Internal minimum check:")
        Rs = sorted(mini_pes.keys())
        Es = [mini_pes[R]['E_gs'] for R in Rs]
        i_min = int(np.argmin(Es))
        R_min = Rs[i_min]
        E_min = Es[i_min]
        has_internal = 0 < i_min < len(Rs) - 1
        print(f"  R_min = {R_min}, E_min = {E_min:.6f}, internal_min = {has_internal}")
        if has_internal:
            print(f"  PI's hoped outcome: BINDING MINIMUM EMERGES at R={R_min}")
        else:
            print(f"  PES still monotone in this range; binding not within {Rs}")
    else:
        mini_pes = None

    # Save full data
    out = {
        'metadata': {
            'sprint': 'F1 Phase 1 + 2',
            'date': '2026-05-23',
            'system': 'NaH max_n=2',
            'R_eq': R_eq,
            'R_diss': R_diss,
            'spec': 'nah_spec(max_n=2)',
            'verdict': verdict,
        },
        'results_eq': results_eq,
        'results_diss': results_diss,
        'de_table': de_table,
        'na_3s_combined_eq': na_3s_combined,
        'de_combined': de_combined,
        'mini_pes': mini_pes,
    }
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print()
    print(f"[saved] {out_path}")


if __name__ == '__main__':
    main()
