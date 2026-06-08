"""Diagnostic: does eri have the chemist 8-fold permutation symmetry?

Run a few sanity tests:
  1. Print sample eri[i,j,k,l] vs equivalent permutations
  2. Apply U = simple orbital swap, check FCI invariance
  3. Apply U = small rotation between two orbitals, check FCI invariance
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as la

from geovac.molecular_spec import lih_spec
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy


def main():
    R = 3.015
    spec = lih_spec(R=R, max_n=2)
    result = build_balanced_hamiltonian(spec, R=R, cross_block_h1=True)
    h1 = result['h1']
    eri = result['eri']
    nuc_rep = result['nuclear_repulsion']
    M = result['M']

    # --- Test 1: 8-fold symmetry of eri ---
    print("[Test 1] Check 8-fold symmetry of eri")
    n_violations = 0
    max_asymm = 0.0
    for p in range(M):
        for q in range(p + 1):
            for r in range(M):
                for s in range(r + 1):
                    if r > p:
                        continue
                    v0 = eri[p, q, r, s]
                    perms = [
                        eri[q, p, r, s],   # p <-> q
                        eri[p, q, s, r],   # r <-> s
                        eri[r, s, p, q],   # 12 <-> 21
                        eri[q, p, s, r],
                        eri[s, r, p, q],
                        eri[r, s, q, p],
                        eri[s, r, q, p],
                    ]
                    for vp in perms:
                        if abs(v0 - vp) > 1e-8:
                            n_violations += 1
                            max_asymm = max(max_asymm, abs(v0 - vp))
    print(f"  8-fold symmetry violations (above 1e-8): {n_violations}")
    print(f"  Max asymmetry: {max_asymm:.4e}")

    # --- Test 2: Apply orbital permutation (swap orbitals 0 and 1) ---
    print("\n[Test 2] Apply orbital swap (0 <-> 1)")
    fake = {'M': M, 'h1': h1, 'eri': eri, 'nuclear_repulsion': nuc_rep}
    fci_baseline = coupled_fci_energy(fake, n_electrons=4)['E_coupled']
    print(f"  Baseline E_FCI = {fci_baseline:.10f}")

    P = np.eye(M)
    P[0, 0] = 0
    P[1, 1] = 0
    P[0, 1] = 1
    P[1, 0] = 1
    # Apply via h1' = P h1 P^T = P h1 P
    h1_p = P @ h1 @ P.T
    eri_p = np.einsum('pa,qb,rc,sd,abcd->pqrs', P, P, P, P, eri,
                      optimize='optimal')
    fake_p = {'M': M, 'h1': h1_p, 'eri': eri_p, 'nuclear_repulsion': nuc_rep}
    fci_perm = coupled_fci_energy(fake_p, n_electrons=4)['E_coupled']
    print(f"  Swapped E_FCI = {fci_perm:.10f}, diff = {abs(fci_perm - fci_baseline):.4e}")

    # --- Test 3: Apply small rotation of orbitals 0 and 1 within Li_core block ---
    print("\n[Test 3] Apply small rotation in Li_core block (theta = 0.1)")
    theta = 0.1
    R_mat = np.eye(M)
    R_mat[0, 0] = np.cos(theta)
    R_mat[0, 1] = -np.sin(theta)
    R_mat[1, 0] = np.sin(theta)
    R_mat[1, 1] = np.cos(theta)
    h1_r = R_mat @ h1 @ R_mat.T
    eri_r = np.einsum('pa,qb,rc,sd,abcd->pqrs', R_mat, R_mat, R_mat, R_mat, eri,
                      optimize='optimal')
    fake_r = {'M': M, 'h1': h1_r, 'eri': eri_r, 'nuclear_repulsion': nuc_rep}
    fci_rot = coupled_fci_energy(fake_r, n_electrons=4)['E_coupled']
    print(f"  Rotated E_FCI = {fci_rot:.10f}, diff = {abs(fci_rot - fci_baseline):.4e}")

    # --- Test 4: Symmetrize eri to 8-fold and test rotation again ---
    print("\n[Test 4] Symmetrize eri to 8-fold first, then rotate")
    def sym8(e):
        s = e.copy()
        s = 0.5 * (s + s.transpose(1, 0, 2, 3))
        s = 0.5 * (s + s.transpose(0, 1, 3, 2))
        s = 0.5 * (s + s.transpose(2, 3, 0, 1))
        return s
    eri_sym = sym8(eri)
    fake_b = {'M': M, 'h1': h1, 'eri': eri_sym, 'nuclear_repulsion': nuc_rep}
    fci_sym_baseline = coupled_fci_energy(fake_b, n_electrons=4)['E_coupled']
    print(f"  Baseline (symmetrized eri) E_FCI = {fci_sym_baseline:.10f}, "
          f"diff vs raw = {abs(fci_sym_baseline - fci_baseline):.4e}")

    eri_r_sym = np.einsum('pa,qb,rc,sd,abcd->pqrs', R_mat, R_mat, R_mat, R_mat,
                          eri_sym, optimize='optimal')
    eri_r_sym = sym8(eri_r_sym)
    fake_rs = {'M': M, 'h1': h1_r, 'eri': eri_r_sym,
               'nuclear_repulsion': nuc_rep}
    fci_rot_sym = coupled_fci_energy(fake_rs, n_electrons=4)['E_coupled']
    print(f"  Rotated (sym-eri) E_FCI = {fci_rot_sym:.10f}, "
          f"diff vs baseline = {abs(fci_rot_sym - fci_sym_baseline):.4e}")


if __name__ == '__main__':
    main()
