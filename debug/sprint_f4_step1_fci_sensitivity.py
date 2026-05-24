"""Sprint F4 Step 1 — FCI sensitivity to bonding-orbital PK shift.

Calibrated estimate: how much does the F3 well depth reduce if we add
a rank-1 PK projector dH = c * |bond><bond| in the cross-block h1?

This BYPASSES the Step 2 production wiring by injecting the predicted
PK barrier directly as a rank-1 perturbation, then re-running FCI at
several R-points. If well depth changes by < 0.5 Ha at the predicted
barrier c = +0.19 Ha (Step 1 prediction), the expected impact is weak
and Step 2 should not be pursued.
"""

from __future__ import annotations

import json
import os
import sys
from itertools import combinations

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)

from geovac.balanced_coupled import build_balanced_hamiltonian  # noqa: E402
from geovac.molecular_spec import nah_spec  # noqa: E402


def build_fci_2e_singlet(h1, eri, M, nuclear_repulsion):
    """Build 2-electron singlet FCI Hamiltonian and return lowest eigvalue.

    Spin-restricted, n_alpha = n_beta = 1.
    """
    n_up = 1
    n_down = 1
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
            for p in alpha:
                for q in beta:
                    E_diag += eri[p, p, q, q]
            H_fci[I, I] = E_diag

    # Single alpha
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
                    for q in beta:
                        val += phase * eri[r, p, q, q]
                    if abs(val) > 1e-14:
                        H_fci[I, J] += val

    # Single beta
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
                    for q in alpha:
                        val += phase * eri[r, p, q, q]
                    if abs(val) > 1e-14:
                        H_fci[I, J] += val

    # Alpha-beta double
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
                            new_beta = tuple(sorted((beta_set - {p_b}) | {r_b}))
                            if new_beta not in beta_idx:
                                continue
                            bi_new = beta_idx[new_beta]
                            phase_b = exc_phase(beta, p_b, r_b)
                            val = phase_a * phase_b * eri[r_a, p_a, r_b, p_b]
                            if abs(val) > 1e-14:
                                I = det_index(ai, bi)
                                J = det_index(ai_new, bi_new)
                                H_fci[I, J] += val

    H_fci = H_fci.tocsr()
    H_fci_sym = 0.5 * (H_fci + H_fci.T)
    if n_det == 1:
        return float(H_fci_sym.toarray()[0, 0])
    E_vals, _ = eigsh(H_fci_sym, k=1, which='SA')
    return float(E_vals[0])


def fci_with_pk_shift(spec, R, pk_shift):
    """Run FCI with rank-1 PK shift on lowest h1 eigenvector."""
    nuclei = [
        {'Z': 11.0, 'position': (0.0, 0.0, 0.0), 'label': 'Na'},
        {'Z': 1.0, 'position': (0.0, 0.0, R), 'label': 'H'},
    ]
    result = build_balanced_hamiltonian(
        spec=spec, R=R, nuclei=nuclei,
        screened_cross_center=True, multi_zeta_basis=True, cross_block_h1=True,
        verbose=False,
    )
    h1 = result['h1'].copy()
    eri = result['eri']
    M = result['M']
    nuc_rep = result['nuclear_repulsion']

    # Diagonalize h1 and add rank-1 PK shift on lowest eigenvector
    w, v = eigh(h1)
    bond_vec = v[:, 0]
    h1_shifted = h1 + pk_shift * np.outer(bond_vec, bond_vec)

    E_gs = build_fci_2e_singlet(h1_shifted, eri, M, nuc_rep)
    return E_gs


def main():
    print("=" * 78)
    print("Sprint F4 Step 1 — FCI sensitivity to bonding-orbital PK rank-1 shift")
    print("=" * 78)
    print("This INJECTS dH = c * |bond><bond| as rank-1 perturbation on h1")
    print("at each R, where |bond> is the lowest h1 eigenvector (the F3 bonding")
    print("orbital). Tests how FCI well depth responds to a PK-like shift.")
    print("=" * 78)

    spec = nah_spec(max_n=2)
    R_values = [2.0, 3.5, 10.0]  # min, midrange, dissociation
    # Test shifts 0, 0.05, 0.10, 0.20, 0.50, 1.0, 2.0, 5.0 Ha
    pk_shifts = [0.0, 0.05, 0.10, 0.20, 0.50, 1.0, 2.0, 5.0]

    results = {}
    for R in R_values:
        results[R] = {}
        for pk in pk_shifts:
            E = fci_with_pk_shift(spec, R, pk)
            results[R][pk] = E
            print(f"  R={R:5.2f} bohr, PK_shift={pk:+.2f} Ha: E_FCI = {E:+.4f} Ha")
        print()

    print("\n" + "=" * 78)
    print("Well depths D_e^(2pt) = E(R=10) - E(R=2)")
    print("=" * 78)
    print(f"{'PK_shift':>12s}  {'D_e (Ha)':>12s}  {'reduction':>12s}  {'fraction':>12s}")
    de_baseline = results[10.0][0.0] - results[2.0][0.0]
    for pk in pk_shifts:
        de = results[10.0][pk] - results[2.0][pk]
        red = de_baseline - de
        frac = red / de_baseline * 100
        print(f"{pk:>+12.2f}  {de:>+12.4f}  {red:>+12.4f}  {frac:>10.2f}%")

    # What PK shift would be needed to bring D_e into experimental range [0.0375, 0.150]?
    print("\n" + "=" * 78)
    print("Threshold analysis")
    print("=" * 78)
    print(f"  F3 baseline D_e (PK=0):                {de_baseline:.4f} Ha")
    print(f"  Experimental NaH D_e:                  ~0.075 Ha")
    print(f"  2x experimental window:                [0.0375, 0.150] Ha")
    print(f"  Step 1 predicted PK barrier (s+2p):    +0.194 Ha")
    de_pred = results[10.0][0.20] - results[2.0][0.20]
    print(f"  D_e at predicted PK shift (0.20 Ha):   {de_pred:.4f} Ha")
    print(f"  Closes wall by:                        {(de_baseline - de_pred)/de_baseline*100:.1f}%")
    if de_pred > 0.150:
        gap_ratio = de_pred / 0.075
        print(f"  Still {gap_ratio:.1f}x experimental D_e")

    out = {
        "sprint": "F4 Step 1 — FCI sensitivity to PK rank-1 shift",
        "date": "2026-05-23",
        "R_values": R_values,
        "pk_shifts": pk_shifts,
        "results": {
            str(R): {str(pk): float(E) for pk, E in res.items()}
            for R, res in results.items()
        },
        "well_depths": {
            str(pk): float(results[10.0][pk] - results[2.0][pk])
            for pk in pk_shifts
        },
        "baseline_de": de_baseline,
        "predicted_pk_barrier_Ha": 0.194,
        "de_at_predicted_pk": float(results[10.0][0.20] - results[2.0][0.20]),
        "pct_wall_closed_by_pk": (de_baseline - (results[10.0][0.20] - results[2.0][0.20])) / de_baseline * 100,
    }
    out_path = os.path.join(_HERE, "data", "sprint_f4_step1_fci_sensitivity.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
