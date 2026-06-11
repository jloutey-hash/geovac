"""Sprint F1 max_n=3 Step 2: single-point FCIs at NaH max_n=3.

Runs 2-electron FCI at NaH max_n=3 (Q_tot=56, M=28 spatial) with three
architectures at R = 3.5 bohr (near experimental R_eq=3.566) and
R = 10.0 bohr (dissociation limit):

  (A) COMBINED: W1c + multi-zeta (the headline configuration)
  (B) W1c alone (no multi-zeta) — for the multi-zeta differential P3 test
  (C) Bare baseline — for sanity check / regression

Reports:
  - E(R=3.5) and E(R=10.0) for each architecture
  - D_e^A = E^A(R=10) - E^A(R=3.5) = combined binding energy
  - D_e^B = E^B(R=10) - E^B(R=3.5) = W1c-alone binding energy
  - multi-zeta differential = E^B(R=3.5) - E^A(R=3.5)
  - FCI 1-RDM natural occupations for (A) — bonding/antibonding partition?
  - Diagonal occupations Na 3s, Na 3p, H 1s, H 2s, H 2p

Decision gate after this:
  - If D_e^A > 0 AND natural occupations show bonding signature (one orbital
    > 1.5, complement < 0.5): proceed to Step 3 mini-PES
  - If D_e^A > 0 but natural occupations still (1, 1): puzzling; proceed
    to mini-PES anyway
  - If D_e^A <= 0: P1 falsified; STOP and document the clean negative
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


# Reuse the FCI machinery from sprint_f1_p1p2_combined_test.py
def _double_exc_phase(occ, p, q, r, s):
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


def build_fci_with_eigenvector(h1, eri, n_electrons, nuclear_repulsion, M,
                                k_states: int = 6):
    """FCI ground state with eigenvector (for 1-RDM analysis).
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


def get_orbital_index_map(spec, ham):
    """Returns dict of orbital label -> spatial index.

    At max_n=3, the Na center sub-block has 14 orbitals:
      block_n=1: (l=0) -> 1 orbital  -> Na 3s
      block_n=2: l=0 (1 orb) + l=1 (3 orbs)  -> 4 orbitals -> Na 4s, 4p_-1, 4p_0, 4p_+1
      block_n=3: l=0 (1 orb) + l=1 (3 orbs) + l=2 (5 orbs) -> 9 orbitals -> Na 5s, 5p_x3, 5d_x5

    H partner sub-block has 14 orbitals (identical structure but Z_orb=1
    hydrogenic).

    Layout: orbitals 0-13 are Na center, 14-27 are H partner.

    By the standard balanced enumeration (CLAUDE.md §12 for NaH layout):
      Na center indices: 0=Na 3s, 1=Na 4s, 2=Na 4p_-1, 3=Na 4p_0, 4=Na 4p_+1,
                         5=Na 5s, 6=Na 5p_-1, 7=Na 5p_0, 8=Na 5p_+1,
                         9-13=Na 5d_{m=-2..+2}
      H partner indices: 14=H 1s, 15=H 2s, 16=H 2p_-1, 17=H 2p_0, 18=H 2p_+1,
                         19=H 3s, 20=H 3p_-1, 21=H 3p_0, 22=H 3p_+1,
                         23-27=H 3d_{m=-2..+2}

    Note: this is the IF-IT-FOLLOWS-THE-CONVENTION mapping. The actual
    enumeration in molecular_spec/balanced_coupled may differ; we verify
    in the test by inspecting the h1 diagonal magnitudes.
    """
    return {
        'Na_3s':   0,
        'Na_4s':   1,
        'Na_4pm':  2, 'Na_4p0':  3, 'Na_4pp': 4,
        'Na_5s':   5,
        'Na_5pm':  6, 'Na_5p0':  7, 'Na_5pp': 8,
        'Na_5dmm': 9, 'Na_5dm': 10, 'Na_5d0': 11, 'Na_5dp': 12, 'Na_5dpp': 13,
        'H_1s':   14,
        'H_2s':   15,
        'H_2pm':  16, 'H_2p0': 17, 'H_2pp': 18,
        'H_3s':   19,
        'H_3pm':  20, 'H_3p0': 21, 'H_3pp': 22,
        'H_3dmm': 23, 'H_3dm': 24, 'H_3d0': 25, 'H_3dp': 26, 'H_3dpp': 27,
    }


def run_arch(spec, R, screened, multi_zeta, label, verbose=False):
    """Build hamiltonian, run FCI, return everything for analysis."""
    n_e = sum(b.n_electrons for b in spec.blocks)
    t0 = time.perf_counter()
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=4000, L_max=4,
        screened_cross_center=screened,
        multi_zeta_basis=multi_zeta,
        verbose=False,
    )
    t_build = time.perf_counter() - t0
    M = ham['M']
    h1 = ham['h1']
    eri = ham['eri']
    nuc_rep = ham['nuclear_repulsion']

    t1 = time.perf_counter()
    E_vals, V_vecs, alpha_strings, beta_strings = build_fci_with_eigenvector(
        h1, eri, n_e, nuc_rep, M, k_states=6,
    )
    t_fci = time.perf_counter() - t1
    psi0 = V_vecs[:, 0]
    E_gs = float(E_vals[0])
    rdm = compute_1rdm(psi0, alpha_strings, beta_strings, M)
    nat_occs = np.sort(np.linalg.eigvalsh(rdm))[::-1]
    elapsed = time.perf_counter() - t0

    # Orbital index map
    orb_map = get_orbital_index_map(spec, ham)

    # Diagonal occupations of named orbitals
    named_diag = {label: float(rdm[idx, idx]) for label, idx in orb_map.items()}

    # h1 eigenspectrum: find lowest 8 eigenvectors and their dominant components
    eigvals_h1, eigvecs_h1 = np.linalg.eigh(h1)
    low8_eigvals = eigvals_h1[:8].tolist()
    low8_dominant = []
    for i in range(8):
        v = eigvecs_h1[:, i]
        amps2 = v ** 2
        top_idx = int(np.argmax(amps2))
        # Find label
        top_label = None
        for k, idx in orb_map.items():
            if idx == top_idx:
                top_label = k
                break
        low8_dominant.append({
            'i': i,
            'eigval': float(eigvals_h1[i]),
            'top_idx': top_idx,
            'top_label': top_label,
            'top_amp2': float(amps2[top_idx]),
            'second_idx': int(np.argsort(amps2)[-2]),
            'second_amp2': float(np.sort(amps2)[-2]),
        })

    # For the dominant natural orbital (highest occupation), get amplitudes
    # on the named orbitals
    rdm_eigvals, rdm_eigvecs = np.linalg.eigh(rdm)
    order = np.argsort(rdm_eigvals)[::-1]
    rdm_eigvals = rdm_eigvals[order]
    rdm_eigvecs = rdm_eigvecs[:, order]
    dominant_no = rdm_eigvecs[:, 0]
    dominant_amplitudes = {}
    for orb_label, idx in orb_map.items():
        a = float(dominant_no[idx])
        if abs(a) > 0.05:
            dominant_amplitudes[orb_label] = a

    # Second natural orbital
    second_no = rdm_eigvecs[:, 1]
    second_amplitudes = {}
    for orb_label, idx in orb_map.items():
        a = float(second_no[idx])
        if abs(a) > 0.05:
            second_amplitudes[orb_label] = a

    # Sanity check: print h1 diagonal magnitudes to verify orbital labeling
    h1_diag = h1.diagonal()

    # Compute Na-H cross-shift content in dominant NO
    # (sum of |amp|² on Na-side orbitals vs H-side)
    na_amp2_dom = sum(dominant_no[idx]**2 for k, idx in orb_map.items() if k.startswith('Na'))
    h_amp2_dom = sum(dominant_no[idx]**2 for k, idx in orb_map.items() if k.startswith('H'))
    na_amp2_2nd = sum(second_no[idx]**2 for k, idx in orb_map.items() if k.startswith('Na'))
    h_amp2_2nd = sum(second_no[idx]**2 for k, idx in orb_map.items() if k.startswith('H'))

    result = {
        'label': label,
        'R': float(R),
        'screened': bool(screened),
        'multi_zeta': bool(multi_zeta),
        'E_gs': E_gs,
        'E_excited_5': E_vals[1:6].tolist() if len(E_vals) > 1 else [],
        'natural_occupations_top_8': nat_occs[:8].tolist(),
        'natural_occ_above_001': int(np.sum(nat_occs > 0.01)),
        'natural_occ_above_010': int(np.sum(nat_occs > 0.10)),
        'natural_occ_above_050': int(np.sum(nat_occs > 0.50)),
        'natural_occ_above_150': int(np.sum(nat_occs > 1.50)),
        'natural_occ_below_050': int(np.sum((nat_occs > 0) & (nat_occs < 0.50))),
        'tr_rdm': float(np.trace(rdm)),
        'rdm_diag_named': named_diag,
        'rdm_diag_full': np.diag(rdm).tolist(),
        'h1_diag_full': h1_diag.tolist(),
        'h1_spectrum_low8': low8_dominant,
        'dominant_no_amplitudes': dominant_amplitudes,
        'second_no_amplitudes': second_amplitudes,
        'dominant_no_na_amp2_total': float(na_amp2_dom),
        'dominant_no_h_amp2_total': float(h_amp2_dom),
        'second_no_na_amp2_total': float(na_amp2_2nd),
        'second_no_h_amp2_total': float(h_amp2_2nd),
        'tr_h1_cross_vne': float(np.trace(ham['h1_cross_vne'])),
        'build_time_s': t_build,
        'fci_time_s': t_fci,
        'total_time_s': elapsed,
        'M': int(M),
        'n_electrons': int(n_e),
    }

    if verbose:
        print(f"  {label:>20}: E = {E_gs:+.6f}, "
              f"top_no={nat_occs[0]:.4f}, "
              f"2nd_no={nat_occs[1]:.4f}, "
              f"dom_NaH_split={na_amp2_dom:.3f}/{h_amp2_dom:.3f}, "
              f"[{elapsed:.1f}s]", flush=True)

    return result


def main():
    out_path = Path('debug/data/sprint_f1_maxn3_step2_fci.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Sprint F1 max_n=3 Step 2: Single-point FCIs at NaH max_n=3")
    print("=" * 70)
    print()

    spec = nah_spec(max_n=3)
    R_eq = 3.5
    R_diss = 10.0

    archs = [
        ('bare_baseline',   False, False),  # (C) bare for regression
        ('w1c_alone',       True,  False),  # (B) W1c, no mz
        ('combined_w1c_mz', True,  True),   # (A) headline configuration
    ]

    results = {
        'sprint': 'F1 max_n=3 Step 2',
        'date': '2026-05-23',
        'system': 'NaH max_n=3',
        'M_spatial': 28,
        'n_electrons': 2,
        'R_eq_bohr': R_eq,
        'R_diss_bohr': R_diss,
        'experimental_R_eq_bohr': 3.566,
        'experimental_D_e_Ha': 0.075,
        'arch_results': {},
    }

    # Run at R_eq
    print(f"--- R = {R_eq} bohr ---")
    results_eq = {}
    for label, screened, mz in archs:
        r = run_arch(spec, R_eq, screened, mz, label, verbose=True)
        results_eq[label] = r
    print()

    # Run at R_diss
    print(f"--- R = {R_diss} bohr ---")
    results_diss = {}
    for label, screened, mz in archs:
        r = run_arch(spec, R_diss, screened, mz, label, verbose=True)
        results_diss[label] = r
    print()

    # Compute D_e for each
    print("--- D_e^2pt = E(R_diss) - E(R_eq) (two-point estimate) ---")
    de_table = {}
    for label, _, _ in archs:
        E_eq = results_eq[label]['E_gs']
        E_diss = results_diss[label]['E_gs']
        D_e = E_diss - E_eq
        de_table[label] = {'E_eq': E_eq, 'E_diss': E_diss, 'D_e': D_e}
        sign = '+' if D_e > 0 else '-'
        print(f"  {label:>20}: E_eq = {E_eq:11.6f}, E_diss = {E_diss:11.6f}, "
              f"D_e = {sign}{abs(D_e):.6f} Ha")
    print()

    # Multi-zeta differential
    print("--- Multi-zeta differential at R_eq ---")
    mz_diff_eq = (results_eq['w1c_alone']['E_gs'] -
                  results_eq['combined_w1c_mz']['E_gs'])
    mz_diff_diss = (results_diss['w1c_alone']['E_gs'] -
                    results_diss['combined_w1c_mz']['E_gs'])
    print(f"  mz_diff(R={R_eq})  = E_W1c - E_W1c+mz = {mz_diff_eq:+.6f} Ha")
    print(f"  mz_diff(R={R_diss}) = E_W1c - E_W1c+mz = {mz_diff_diss:+.6f} Ha")
    print(f"  P3 range [0.02, 0.30] Ha. Actual at R_eq: |{abs(mz_diff_eq):.4f}| Ha")
    print()

    # Natural occupations + bonding partition diagnostic
    print("--- Natural occupations (top 6) at R_eq ---")
    for label, _, _ in archs:
        r = results_eq[label]
        occs = r['natural_occupations_top_8'][:6]
        print(f"  {label:>20}: {['%.4f' % o for o in occs]}")
    print()

    # Bonding signature check on combined
    print("--- Bonding partition signature on COMBINED (W1c+mz) at R_eq ---")
    combined = results_eq['combined_w1c_mz']
    occs = combined['natural_occupations_top_8']
    has_bonding = (occs[0] > 1.5 and occs[1] < 0.5)
    has_two_singly_occ = (1.5 > occs[0] > 0.5 and 1.5 > occs[1] > 0.5)
    print(f"  Top occupation: {occs[0]:.4f}, 2nd: {occs[1]:.4f}, "
          f"3rd: {occs[2]:.4f}, 4th: {occs[3]:.4f}")
    print(f"  Bonding-dominant signature (top > 1.5 AND 2nd < 0.5)? "
          f"{'YES' if has_bonding else 'NO'}")
    print(f"  Two-singly-occupied signature (both in [0.5, 1.5])? "
          f"{'YES' if has_two_singly_occ else 'NO'}")
    print()

    # Cross-shift content in dominant NO
    print("--- Cross-shift (Na-H) content in dominant natural orbital ---")
    for label, _, _ in archs:
        r = results_eq[label]
        na = r['dominant_no_na_amp2_total']
        hh = r['dominant_no_h_amp2_total']
        print(f"  {label:>20}: dom_NO  Na_amp2 = {na:.4f}, H_amp2 = {hh:.4f}")
        amps = r['dominant_no_amplitudes']
        # Show top 4 components
        srt = sorted(amps.items(), key=lambda x: -abs(x[1]))[:4]
        comps_str = ' '.join(f'{k}={v:+.3f}' for k, v in srt)
        print(f"  {' ':>20}  top: {comps_str}")
    print()

    print("--- Cross-shift (Na-H) content in 2nd natural orbital ---")
    for label, _, _ in archs:
        r = results_eq[label]
        na = r['second_no_na_amp2_total']
        hh = r['second_no_h_amp2_total']
        print(f"  {label:>20}: 2nd_NO  Na_amp2 = {na:.4f}, H_amp2 = {hh:.4f}")
        amps = r['second_no_amplitudes']
        srt = sorted(amps.items(), key=lambda x: -abs(x[1]))[:4]
        comps_str = ' '.join(f'{k}={v:+.3f}' for k, v in srt)
        print(f"  {' ':>20}  top: {comps_str}")
    print()

    # Diagonal occupations of named orbitals (combined, R_eq)
    print("--- Diagonal occupations of named orbitals (COMBINED @ R_eq) ---")
    named = results_eq['combined_w1c_mz']['rdm_diag_named']
    for k, v in sorted(named.items(), key=lambda x: -x[1]):
        if v > 0.001:
            print(f"  {k:>10}: {v:.4f}")
    print()

    # Decision gate output
    print("=" * 70)
    print("DECISION GATE")
    print("=" * 70)
    de_combined = de_table['combined_w1c_mz']['D_e']
    top_occ_combined = combined['natural_occupations_top_8'][0]
    print(f"  D_e (combined, 2-point)  = {de_combined:+.6f} Ha")
    print(f"  Top natural occupation    = {top_occ_combined:.4f}")
    print()
    if de_combined > 0 and has_bonding:
        verdict = 'PROCEED_TO_MINI_PES (binding + bonding signature)'
    elif de_combined > 0 and has_two_singly_occ:
        verdict = 'PROCEED_TO_MINI_PES (binding, but two-singly-occupied — puzzling)'
    elif de_combined <= 0:
        verdict = 'STOP_P1_FALSIFIED (no binding)'
    else:
        verdict = 'PROCEED_TO_MINI_PES (uncertain signature)'
    print(f"  Verdict: {verdict}")
    print()

    results['arch_results'] = {
        'R_eq': {label: results_eq[label] for label, _, _ in archs},
        'R_diss': {label: results_diss[label] for label, _, _ in archs},
    }
    results['D_e_2point'] = de_table
    results['multi_zeta_differential'] = {
        'at_R_eq': mz_diff_eq,
        'at_R_diss': mz_diff_diss,
        'P3_range_Ha': [0.02, 0.30],
        'P3_pass_at_R_eq': 0.02 <= abs(mz_diff_eq) <= 0.30,
    }
    results['bonding_signature'] = {
        'top_natural_occupation': top_occ_combined,
        'has_bonding_signature': has_bonding,
        'has_two_singly_occupied_signature': has_two_singly_occ,
    }
    results['decision_gate'] = {
        'D_e_combined': de_combined,
        'verdict': verdict,
    }

    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Saved to {out_path}")


if __name__ == '__main__':
    main()
