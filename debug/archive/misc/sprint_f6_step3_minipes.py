"""Sprint F6 Step 3 - NaH max_n=4 mini-PES with full F3 stack.

5-point PES at R in {2.0, 3.0, 3.566, 5.0, 10.0} bohr.
Triggered conditionally on Step 2 PARTIAL_IMPROVEMENT verdict
(D_e improved from 4.37 Ha to 3.23 Ha at max_n=4).

Purpose: characterize PES shape across R; check for internal minimum;
confirm bonding-orbital construction is robust across PES.

Compute estimate: ~5 minutes per build x 5 R-points = ~25 minutes.
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


def build_fci(h1, eri, n_electrons, nuclear_repulsion, M, k_states=4):
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


def run_arch_one_R(spec, R, **balanced_kwargs):
    n_e = sum(b.n_electrons for b in spec.blocks)
    t0 = time.perf_counter()
    print(f"  Building at R={R}...", flush=True)
    ham = build_balanced_hamiltonian(spec, R=R, verbose=False, **balanced_kwargs)
    t_build = time.perf_counter() - t0
    print(f"  Build complete ({t_build:.1f}s = {t_build/60:.2f}min)", flush=True)

    M = ham['M']
    h1 = ham['h1']
    eri = ham['eri']
    nuc_rep = ham['nuclear_repulsion']

    print(f"  FCI...", flush=True)
    t_fci = time.perf_counter()
    E_vals, V_vecs, alpha_strings, beta_strings = build_fci(
        h1, eri, n_e, nuc_rep, M,
    )
    elapsed_fci = time.perf_counter() - t_fci
    print(f"  FCI done ({elapsed_fci:.1f}s)", flush=True)
    psi0 = V_vecs[:, 0]
    E_gs = float(E_vals[0])
    rdm = compute_1rdm(psi0, alpha_strings, beta_strings, M)
    nat_occs_v, nat_vecs_v = np.linalg.eigh(rdm)
    order = np.argsort(nat_occs_v)[::-1]
    nat_occs_v = nat_occs_v[order]
    nat_vecs_v = nat_vecs_v[:, order]
    M_per_side = M // 2
    dom = nat_vecs_v[:, 0]
    amp_Na = float(np.sum(dom[:M_per_side] ** 2))
    amp_H = float(np.sum(dom[M_per_side:] ** 2))
    elapsed = time.perf_counter() - t0
    print(f"  R={R}: E_gs={E_gs:+.6f}  NO_dom={nat_occs_v[0]:.4f}  Na/H={amp_Na:.3f}/{amp_H:.3f}", flush=True)
    return {
        'R': R,
        'E_gs': E_gs,
        'dom_NO_occ': float(nat_occs_v[0]),
        'Na_amp_squared': amp_Na,
        'H_amp_squared': amp_H,
        'natural_occupations_top6': nat_occs_v[:6].tolist(),
        'na3s_coeff': float(dom[0]),
        'h1s_coeff': float(dom[M_per_side]),
        'character': 'bonding' if (dom[0] * dom[M_per_side] > 0) else 'antibonding',
        'build_time_s': t_build,
        'elapsed_s': elapsed,
    }


def main():
    out_path = Path('debug/data/sprint_f6_step3_minipes.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    spec = nah_spec(max_n=4)
    print("=" * 70)
    print(f"Sprint F6 Step 3 - NaH max_n=4 mini-PES with full F3 stack")
    print(f"  M=60, Q=120, fci_dim=3600, 2e singlet")
    print(f"  R in (2.0, 3.0, 3.566, 5.0, 10.0)")
    print("=" * 70)

    R_grid = [2.0, 3.0, 3.566, 5.0, 10.0]
    kwargs = dict(screened_cross_center=True, multi_zeta_basis=True,
                  cross_block_h1=True)

    results = []
    t_start_total = time.perf_counter()
    for R in R_grid:
        print(f"\n  ==== R = {R} bohr ====")
        r = run_arch_one_R(spec, R, **kwargs)
        results.append(r)
    t_total = time.perf_counter() - t_start_total

    print(f"\n{'=' * 70}")
    print(f"F6 mini-PES table (max_n=4 full F3 stack):")
    print(f"{'R':>8} {'E_gs':>16} {'NO_dom':>10} {'Na/H':>14} {'char':>12}")
    print("-" * 70)
    for r in results:
        print(f"{r['R']:>8.3f} {r['E_gs']:>+16.6f} {r['dom_NO_occ']:>10.4f} "
              f"{r['Na_amp_squared']:.3f}/{r['H_amp_squared']:.3f} {r['character']:>12}")

    # Diagnose PES shape
    R_arr = np.array([r['R'] for r in results])
    E_arr = np.array([r['E_gs'] for r in results])
    R_min_idx = int(np.argmin(E_arr))
    R_min = R_arr[R_min_idx]
    E_min = float(E_arr[R_min_idx])
    is_internal_min = R_min_idx > 0 and R_min_idx < len(R_arr) - 1
    # well depth (E_min - E_dissoc)
    E_dissoc = float(E_arr[-1])  # R=10
    well_depth = E_dissoc - E_min  # positive means R_min is more bound

    # Compare against the experimental NaH R_eq=3.566 bohr
    R_eq_exp = 3.566
    R_eq_dist = abs(R_min - R_eq_exp)

    D_e_exp = 0.075
    F3_baseline_D_e = 4.374
    wall_closure_fraction = (F3_baseline_D_e - well_depth) / F3_baseline_D_e * 100

    print(f"\n{'=' * 70}")
    print(f"PES SHAPE DIAGNOSIS (max_n=4 F3 stack):")
    print(f"  R_min      = {R_min:.3f} bohr (idx {R_min_idx} of {len(R_arr)})")
    print(f"  E_min      = {E_min:+.6f} Ha")
    print(f"  E_dissoc   = {E_dissoc:+.6f} Ha")
    print(f"  well depth = {well_depth:+.6f} Ha")
    print(f"  internal_minimum = {is_internal_min}")
    print(f"  R_eq_exp   = {R_eq_exp:.3f} bohr, |R_min - R_eq_exp| = {R_eq_dist:.3f}")
    print(f"  D_e_exp    = {D_e_exp:.3f} Ha")
    print(f"  F3 baseline D_e = {F3_baseline_D_e:.3f} Ha")
    print(f"  wall closure fraction = {wall_closure_fraction:+.1f}%")
    print(f"{'=' * 70}")

    # F6 verdict
    no_dom_char_at_min = results[R_min_idx]['character']
    na_amp_at_min = results[R_min_idx]['Na_amp_squared']
    h_amp_at_min = results[R_min_idx]['H_amp_squared']
    bonding_signature = (
        no_dom_char_at_min == 'bonding' and na_amp_at_min > 0.1 and h_amp_at_min > 0.1
    )
    magnitude_in_window = 0.5 * D_e_exp <= well_depth <= 2.0 * D_e_exp

    if well_depth > 0 and bonding_signature and magnitude_in_window and is_internal_min:
        verdict = "WALL_CLOSES"
    elif 0.15 < well_depth <= 1.0:
        verdict = "PARTIAL"
    elif 1.0 < well_depth < F3_baseline_D_e:
        verdict = "PARTIAL_IMPROVEMENT"
    elif well_depth >= F3_baseline_D_e or well_depth <= 0:
        verdict = "WALL_PERSISTS"
    else:
        verdict = "INTERMEDIATE"

    print(f"\nF6 verdict (PES): {verdict}")
    print(f"  Bonding signature: {bonding_signature}")
    print(f"  Magnitude in window: {magnitude_in_window}")
    print(f"  Internal minimum: {is_internal_min}")
    print(f"  Total elapsed: {t_total:.1f}s = {t_total/60:.1f}min")
    print(f"{'=' * 70}")

    out = {
        'sprint': 'F6 Step 3 - NaH max_n=4 mini-PES',
        'date': '2026-05-23',
        'spec': 'nah_spec(max_n=4), Q=120, 2e singlet FCI, M=60, fci_dim=3600',
        'architecture': 'W1c + multi-zeta + cross-block h1',
        'R_grid': R_grid,
        'pes': results,
        'R_min_bohr': R_min,
        'R_min_idx': R_min_idx,
        'E_min_Ha': E_min,
        'E_dissoc_Ha': E_dissoc,
        'well_depth_Ha': well_depth,
        'internal_minimum': is_internal_min,
        'R_eq_exp_bohr': R_eq_exp,
        'R_eq_distance_bohr': R_eq_dist,
        'D_e_exp_Ha': D_e_exp,
        'F3_baseline_D_e_Ha': F3_baseline_D_e,
        'wall_closure_fraction_pct': wall_closure_fraction,
        'bonding_signature': bonding_signature,
        'magnitude_in_window': magnitude_in_window,
        'verdict': verdict,
        'total_elapsed_s': t_total,
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote: {out_path}")


if __name__ == "__main__":
    main()
