"""
Track CUSP-1: w/d-screened CI truncation.

EP-2 identified S_B = A * (w/d)^gamma. Per-config contribution to S
is dominated by configs with large |V_IJ|^2 / (E_J - E_I)^2 in the
H_kin eigenbasis (2nd-order RS).

Strategy: build full graph-native FCI for He at n_max=7, diagonalize
H_kin to get config-space eigenstates, compute per-config 2nd-order RS
score s_I = sum_{J != I_GS} |V_IJ|^2 / (E_J - E_I_GS)^2 in that basis,
then keep only the top-K configs by score and re-diagonalize.

Compare He energy at full vs screened CI as a function of K:
- At what K does screened CI match the full CI energy to within 0.01 mHa?
- What is the energy error reduction per retained config?
- Does the screening prefer high-|m| configs (cusp signature)?

Output: debug/data/cusp1_screened_ci.json
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.casimir_ci import (
    build_graph_native_fci, _build_orbital_basis, _build_graph_h1,
    two_electron_integral,
)


def _split_h_kin_h_vee(H_full, Z, n_max, m_total=0):
    """Recover H_kin and H_vee separately by re-building."""
    h1_mat, orbitals = _build_graph_h1(Z, n_max)
    n_spatial = len(orbitals)
    configs = []
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if orbitals[i][2] + orbitals[j][2] == m_total:
                configs.append((i, j))
    n_cfg = len(configs)

    H_kin = np.zeros((n_cfg, n_cfg))
    for I, (i, j) in enumerate(configs):
        for J in range(I, n_cfg):
            p, q = configs[J]
            bra = [(i, j, 1.0)] + ([(j, i, 1.0)] if i != j else [])
            ket = [(p, q, 1.0)] + ([(q, p, 1.0)] if p != q else [])
            N_I = np.sqrt(float(len(bra))); N_J = np.sqrt(float(len(ket)))
            me = 0.0
            for a, b, s_b in bra:
                for c, d, s_k in ket:
                    if b == d:
                        me += s_b * s_k * h1_mat[a, c]
                    if a == c:
                        me += s_b * s_k * h1_mat[b, d]
            me /= (N_I * N_J)
            H_kin[I, J] = me
            H_kin[J, I] = me

    H_vee = H_full - H_kin
    return H_kin, H_vee, configs, orbitals


def _per_config_RS_score(H_kin, H_vee):
    """Per-config 2nd-order RS perturbation score relative to GS of H_kin.

    s_J = |<I_GS | V | J>|^2 / (E_J - E_I_GS)^2
    where I_GS is the H_kin ground-state config (lowest E_kin).
    """
    e_kin, U = np.linalg.eigh(H_kin)
    V_in = U.T @ H_vee @ U  # transform V to H_kin eigenbasis
    e_GS = e_kin[0]
    scores = np.zeros(len(e_kin))
    for J in range(len(e_kin)):
        if J == 0:
            scores[J] = float('inf')  # GS always retained
            continue
        gap = e_kin[J] - e_GS
        if gap < 1e-10:
            scores[J] = float('inf')
        else:
            scores[J] = (V_in[0, J] ** 2) / gap ** 2
    return scores, e_kin, U


def _truncated_energy(H_full_kineig, K, U):
    """Pick top-K config-eigenstates by score, transform, diagonalize."""
    # We work in the H_kin eigenbasis: H_full transformed = U.T @ H_full @ U.
    # Picking top-K eigenstates corresponds to truncating the basis.
    # Need score-based selection: GS always retained, then top K-1 by RS score.
    pass


def main():
    Z = 2
    n_max = 7
    m_total = 0

    print(f"CUSP-1: w/d-screened CI for He at n_max={n_max}")
    print("=" * 60)

    # Reference: full graph-native CI
    t0 = time.time()
    H_full = build_graph_native_fci(Z=Z, n_max=n_max, m_total=m_total)
    eigs_full = np.linalg.eigvalsh(H_full)
    E_full = float(eigs_full[0])
    n_cfg = H_full.shape[0]
    print(f"Full CI: n_configs = {n_cfg}, E_full = {E_full:.6f} Ha "
          f"({time.time() - t0:.1f}s)")

    E_exact_He = -2.903724377
    err_full = 100 * abs(E_full - E_exact_He) / abs(E_exact_He)
    print(f"  Error vs exact: {err_full:.4f}%")

    # Decompose H_full -> H_kin + H_vee
    H_kin, H_vee, configs, orbitals = _split_h_kin_h_vee(H_full, Z, n_max, m_total)

    # Score each config-eigenstate of H_kin by 2nd-order RS perturbation
    scores, e_kin, U = _per_config_RS_score(H_kin, H_vee)

    # Transform full H to H_kin eigenbasis
    H_in_eig = U.T @ H_full @ U

    # Sweep K: keep top K config-eigenstates by RS score
    sorted_indices = np.argsort(-scores)  # descending
    print("\nScreened CI sweep:")
    print(f"  {'K':>4}  {'E_screened':>12}  {'err%':>10}  {'dE/d_log_K':>12}")

    K_grid = []
    K_max = n_cfg
    K_curr = 5
    while K_curr <= K_max:
        K_grid.append(K_curr)
        K_curr = max(K_curr + 5, int(K_curr * 1.4))
    if K_grid[-1] != n_cfg:
        K_grid.append(n_cfg)

    sweep_results = []
    prev_E = None
    prev_K = None
    for K in K_grid:
        keep = sorted_indices[:K]
        keep = np.sort(keep)
        H_sub = H_in_eig[np.ix_(keep, keep)]
        E_K = float(np.linalg.eigvalsh(H_sub)[0])
        err_K = 100 * abs(E_K - E_exact_He) / abs(E_exact_He)
        gain = ((prev_E - E_K) / np.log(K / prev_K)) if (prev_E is not None) else None
        gain_str = f"{gain:.4e}" if gain is not None else "  -- "
        print(f"  {K:4d}  {E_K:12.6f}  {err_K:10.4f}  {gain_str:>12}")
        sweep_results.append({'K': K, 'E': E_K, 'err_pct': err_K})
        prev_E = E_K
        prev_K = K

    # Find K such that screened error is within 0.01 percentage points of full
    target = err_full + 0.01
    K_match = None
    for r in sweep_results:
        if r['err_pct'] <= target:
            K_match = r['K']
            break
    if K_match is not None:
        compression = n_cfg / K_match
        print(f"\nAccuracy match: K_match = {K_match} configs "
              f"(of {n_cfg}, compression {compression:.1f}x) reaches err+0.01pp")

    # Inspect what's in the top scoring set: are they cusp-related?
    print("\nTop 10 most-screened-in configs (by RS score)...")
    print("Note: scores are in the H_kin eigenbasis, NOT original config basis;")
    print("the GS-eigenstate has high overlap with original (1s,1s) config.")
    e_pairs = sorted([(s, idx) for idx, s in enumerate(scores) if not np.isinf(s)],
                     reverse=True)[:10]
    for s, idx in e_pairs:
        # Decompose this H_kin eigenstate back into config basis
        coeffs = U[:, idx]
        top3 = np.argsort(-np.abs(coeffs))[:3]
        labels = []
        for c in top3:
            i, j = configs[c]
            labels.append(f"({orbitals[i]},{orbitals[j]}):{coeffs[c]:.3f}")
        print(f"  score={s:.3e}  E_kin_state[{idx}]: {' + '.join(labels)}")

    out = {
        'Z': Z, 'n_max': n_max, 'm_total': m_total,
        'n_configs_full': n_cfg,
        'E_full_CI': E_full,
        'err_full_pct': err_full,
        'sweep': sweep_results,
        'K_match_within_0p01pp': K_match,
        'compression_at_match': float(n_cfg / K_match) if K_match else None,
    }
    out_path = os.path.join(os.path.dirname(__file__), 'data',
                            'cusp1_screened_ci.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved -> {out_path}")


if __name__ == '__main__':
    main()
