"""Sprint F3 Step 3 — NaH 2-electron FCI with full F3 stack.

Runs the 2-electron FCI at NaH max_n=2 (Q=20, M=10) with the full
F3 stack at R = 3.5 bohr and R = 10.0 bohr:
  - W1c screening (screened_cross_center=True)
  - Multi-zeta basis (multi_zeta_basis=True)
  - Cross-block h1 architectural extension (cross_block_h1=True)

Reports:
  - E(R=3.5), E(R=10), D_e = E(10) - E(3.5)
  - Natural orbital occupations (looking for [2, 0]-like bonding signature
    vs [1, 1] separately-occupied)
  - Dominant natural orbital Na/H amplitude split
  - h1 eigenvalues at R=3.5 (bonding vs antibonding ordering check)

Comparison to F1 baselines:
  - F1-P1+P2 max_n=2 W1c + mz alone (from sprint_f1_p1p2_extended_pes.json):
    R=3.5: E=-163.115 Ha, R=10: E=-162.779 Ha, D_e = +0.336 Ha
    Naturals: [1, 1] (no bonding signature)
  - F1 max_n=3 W1c + mz alone: D_e^PES = +0.71 Ha spurious binding

Decision gate after this:
  - D_e^F3 > 0 AND naturals show bonding signature (e.g., 1.9, 0.1) AND
    D_e^F3 within 2x experimental NaH 0.075 Ha [0.0375, 0.15]: WALL-CLOSES
  - D_e^F3 > 0 but naturals still [1, 1] or magnitude outside window:
    PARTIAL-CLOSURE; PES shape needed
  - D_e^F3 <= 0: cross-block h1 didn't close; STOP
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


# ---------------------------------------------------------------------------
# 2-electron FCI machinery (singlet, n_up=n_down=1)
# ---------------------------------------------------------------------------

def build_fci_with_eigenvector(h1, eri, n_electrons, nuclear_repulsion, M,
                                k_states: int = 6):
    """Build FCI Hamiltonian and return lowest k_states eigenpairs.

    Spin-restricted: n_up = n_down = N/2.
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
    E_vals, V_vecs = eigsh(H_fci, k=k, which='SA')
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


def diagnose_h1_spectrum(h1):
    """Eigendecompose h1 and report the lowest 4 eigenvalues +
    Na/H character of the lowest 4 eigenvectors.

    For NaH max_n=2 with offsets:
       Na center (block_n=1,l=0,m=0): orbital 0 -- "Na 3s"
       Na center (block_n=2,l=0,m=0): orbital 1 -- "Na 4s"
       Na center (block_n=2,l=1,m=-1,0,+1): orbitals 2,3,4 -- "Na 4p"
       H partner (block_n=1,l=0,m=0): orbital 5 -- "H 1s"
       H partner (block_n=2,l=0,m=0): orbital 6 -- "H 2s"
       H partner (block_n=2,l=1,m=-1,0,+1): orbitals 7,8,9 -- "H 2p"

    Returns dict with eigenvalues and Na/H amplitude split per lowest eigvec.
    """
    eigs, vecs = np.linalg.eigh(h1)
    out = {
        'eigenvalues': eigs.tolist(),
        'lowest_4_decomposition': [],
    }
    for k in range(min(4, len(eigs))):
        v = vecs[:, k]
        amp_Na = float(np.sum(v[:5] ** 2))
        amp_H = float(np.sum(v[5:] ** 2))
        out['lowest_4_decomposition'].append({
            'eigenvalue_Ha': float(eigs[k]),
            'Na_amplitude_squared': amp_Na,
            'H_amplitude_squared': amp_H,
            'na3s_coeff': float(v[0]),
            'h1s_coeff': float(v[5]),
            'character': 'bonding' if (v[0] * v[5] > 0) else 'antibonding',
        })
    return out


def run_arch(spec, R, label, **balanced_kwargs):
    """Build hamiltonian, run FCI, return diagnostic."""
    n_e = sum(b.n_electrons for b in spec.blocks)
    t0 = time.perf_counter()
    ham = build_balanced_hamiltonian(spec, R=R, verbose=False, **balanced_kwargs)
    M = ham['M']
    h1 = ham['h1']
    eri = ham['eri']
    nuc_rep = ham['nuclear_repulsion']

    E_vals, V_vecs, alpha_strings, beta_strings = build_fci_with_eigenvector(
        h1, eri, n_e, nuc_rep, M,
    )
    psi0 = V_vecs[:, 0]
    E_gs = float(E_vals[0])
    rdm = compute_1rdm(psi0, alpha_strings, beta_strings, M)
    nat_occs = np.sort(np.linalg.eigvalsh(rdm))[::-1]
    nat_vecs = None
    # Diagonalize RDM for natural orbital signatures
    nat_occs_v, nat_vecs_v = np.linalg.eigh(rdm)
    order = np.argsort(nat_occs_v)[::-1]
    nat_occs_v = nat_occs_v[order]
    nat_vecs_v = nat_vecs_v[:, order]
    # Dominant NO Na/H character
    dom = nat_vecs_v[:, 0]
    amp_Na = float(np.sum(dom[:5] ** 2))
    amp_H = float(np.sum(dom[5:] ** 2))

    h1_spec = diagnose_h1_spectrum(h1)

    elapsed = time.perf_counter() - t0
    return {
        'label': label,
        'R': R,
        'E_gs': E_gs,
        'natural_occupations': nat_occs.tolist(),
        'dominant_NO': {
            'occupation': float(nat_occs_v[0]),
            'Na_amp_squared': amp_Na,
            'H_amp_squared': amp_H,
            'na3s_coeff': float(dom[0]),
            'h1s_coeff': float(dom[5]),
            'character': 'bonding' if (dom[0] * dom[5] > 0) else 'antibonding',
        },
        'rdm_diag': np.diag(rdm).tolist(),
        'h1_spectrum_diag': h1_spec,
        'cross_block_h1_info': ham.get('cross_block_h1_info', {}),
        'elapsed_s': elapsed,
    }


def main():
    out_path = Path('debug/data/sprint_f3_step3_fci.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    spec = nah_spec(max_n=2)
    Q = spec.blocks[0].max_n  # just for diag
    print("=" * 70)
    print(f"Sprint F3 Step 3 — NaH 2e FCI with full F3 stack (max_n={Q})")
    print("=" * 70)

    architectures = [
        ('bare',
            dict(screened_cross_center=False, multi_zeta_basis=False,
                 cross_block_h1=False)),
        ('W1c+mz',
            dict(screened_cross_center=True, multi_zeta_basis=True,
                 cross_block_h1=False)),
        ('W1c+mz+xblockh1',
            dict(screened_cross_center=True, multi_zeta_basis=True,
                 cross_block_h1=True)),
    ]

    results = {}
    for label, kwargs in architectures:
        print(f"\n--- {label} ---")
        for R in [3.5, 10.0]:
            r = run_arch(spec, R, label, **kwargs)
            print(f"  R={R:.1f}: E_gs={r['E_gs']:+.6f}  "
                  f"NO_dom={r['dominant_NO']['occupation']:.4f} "
                  f"({r['dominant_NO']['character']}) "
                  f"Na/H={r['dominant_NO']['Na_amp_squared']:.3f}/"
                  f"{r['dominant_NO']['H_amp_squared']:.3f}")
            print(f"          natural occs: {[f'{o:.4f}' for o in r['natural_occupations'][:4]]}")
            print(f"          h1 lowest 3 eigvals: "
                  f"{[f'{e:+.4f}' for e in r['h1_spectrum_diag']['eigenvalues'][:3]]}")
            print(f"          lowest h1 eigvec: "
                  f"Na/H = {r['h1_spectrum_diag']['lowest_4_decomposition'][0]['Na_amplitude_squared']:.3f}/"
                  f"{r['h1_spectrum_diag']['lowest_4_decomposition'][0]['H_amplitude_squared']:.3f} "
                  f"({r['h1_spectrum_diag']['lowest_4_decomposition'][0]['character']})")
            results.setdefault(label, []).append(r)

    # Compute D_e for each architecture
    summary = {}
    for label in [a[0] for a in architectures]:
        Es = {r['R']: r['E_gs'] for r in results[label]}
        D_e = Es[10.0] - Es[3.5]
        summary[label] = {
            'E_R3.5': Es[3.5],
            'E_R10.0': Es[10.0],
            'D_e_2point_Ha': D_e,
            'binding_2point': D_e > 0,
        }
        print(f"\n  {label}: D_e (2-point) = {D_e:+.6f} Ha (binding={'YES' if D_e > 0 else 'NO'})")

    # F3 verdict
    f3_D_e = summary['W1c+mz+xblockh1']['D_e_2point_Ha']
    f3_naturals = results['W1c+mz+xblockh1'][0]['natural_occupations'][:2]
    f3_no_dom_char = results['W1c+mz+xblockh1'][0]['dominant_NO']['character']
    f3_na_amp = results['W1c+mz+xblockh1'][0]['dominant_NO']['Na_amp_squared']
    f3_h_amp = results['W1c+mz+xblockh1'][0]['dominant_NO']['H_amp_squared']

    # Bonding signature: dominant NO has both Na and H amplitude > 0.1 (mixing)
    bonding_signature = (
        f3_no_dom_char == 'bonding' and f3_na_amp > 0.1 and f3_h_amp > 0.1
    )

    # Magnitude in 2x experimental window
    D_e_exp = 0.075  # Ha
    magnitude_in_window = 0.5 * D_e_exp <= f3_D_e <= 2.0 * D_e_exp

    if f3_D_e > 0 and bonding_signature and magnitude_in_window:
        verdict = "WALL_CLOSES"
    elif f3_D_e > 0 and bonding_signature:
        verdict = "PARTIAL_CLOSURE_MAGNITUDE_OFF"
    elif f3_D_e > 0:
        verdict = "PARTIAL_CLOSURE_NO_BONDING_SIG"
    else:
        verdict = "WALL_PERSISTS"

    print(f"\n{'=' * 70}")
    print(f"F3 verdict (2-point gate): {verdict}")
    print(f"  D_e^F3 = {f3_D_e:+.6f} Ha (experimental NaH: {D_e_exp} Ha)")
    print(f"  Naturals: {[f'{o:.4f}' for o in f3_naturals]}")
    print(f"  Bonding signature: {bonding_signature}")
    print(f"  Magnitude in 2x window [0.0375, 0.15]: {magnitude_in_window}")
    print(f"{'=' * 70}")

    out = {
        'sprint': 'F3 Step 3 — NaH 2e FCI with cross_block_h1',
        'date': '2026-05-23',
        'spec': 'nah_spec(max_n=2), Q=20, 2-electron singlet FCI',
        'architectures': {k: results[k] for k in results},
        'summary': summary,
        'verdict': verdict,
        'F3_de_2point_Ha': f3_D_e,
        'F3_naturals': f3_naturals,
        'F3_bonding_signature': bonding_signature,
        'F3_magnitude_in_window': magnitude_in_window,
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote: {out_path}")


if __name__ == "__main__":
    main()
