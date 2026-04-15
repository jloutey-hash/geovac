"""
Track EP-2N: Be analytical degenerate-PT entanglement.

EP-2j showed the bare GeoVac S^3 graph H_1 has a degenerate ground
state for Be at n_max=2 (closed-shell ^1S, M_L=0, M_S=0):
  H_1 ground subspace = span{|1s^2 2s^2>, |1s^2 (2p_-1 2p_+1)_S=0>,
                              |1s^2 2p_0^2>}
all at energy E_kin = -5Z^2/4.

We diagonalize V_ee analytically in this 3-dim subspace
(degenerate perturbation theory, first order), recover the V_ee-
corrected GS as a specific linear combination, compute its 1-RDM
analytically, and extract the entropy.

Output: debug/data/ep2N_be_analytical.json
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.casimir_ci import two_electron_integral


# ---------------------------------------------------------------------------
# Module-scope helpers (also used by tests).
# ---------------------------------------------------------------------------

def two_e_int_singlet_pair(orb_a, orb_b, orb_c, orb_d, Z):
    """<phi_a phi_b | 1/r12 | phi_c phi_d> via Gaunt-expansion engine."""
    na, la, ma = orb_a; nb, lb, mb = orb_b
    nc, lc, mc = orb_c; nd, ld, md = orb_d
    return two_electron_integral(na, la, ma, nb, lb, mb,
                                 nc, lc, mc, nd, ld, md, k_orb=Z)


def diag_E_for_pair(b_a, b_b, Z=4):
    """Diagonal V_ee for closed-shell |1s^2, b_a, b_b; spatial-singlet>.

    If b_a == b_b: closed-shell sub-shell |1s^2 b^2>.
    Else (open n=2 pair, spatial singlet): symmetric.
    """
    s1 = (1, 0, 0)
    def J(x, y): return two_e_int_singlet_pair(x, y, x, y, Z)
    def K(x, y): return two_e_int_singlet_pair(x, y, y, x, Z)
    if b_a == b_b:
        J_ss = J(s1, s1); J_bb = J(b_a, b_a)
        J_sb = J(s1, b_a); K_sb = K(s1, b_a)
        return J_ss + J_bb + 4 * J_sb - 2 * K_sb
    return (J(s1, s1) + 2*J(s1, b_a) - K(s1, b_a)
            + 2*J(s1, b_b) - K(s1, b_b)
            + J(b_a, b_b) + K(b_a, b_b))


def main():
    Z = 4
    print(f"EP-2N: Be analytical degenerate PT (Z={Z}, n_max=2)")
    print("=" * 60)

    s2 = (2, 0, 0); p_m = (2, 1, -1); p_0 = (2, 1, 0); p_p = (2, 1, +1)

    E_kin = 2 * (-Z**2 / 2 / 1) + 2 * (-Z**2 / 2 / 4)
    print(f"E_kin (degenerate, all 3 configs): {E_kin} Ha")

    E_A = diag_E_for_pair(s2, s2, Z)        # |1s^2 2s^2>
    E_B = diag_E_for_pair(p_0, p_0, Z)      # |1s^2 2p_0^2>
    E_C = diag_E_for_pair(p_m, p_p, Z)      # |1s^2 (2p_-1 2p_+1)_S=0>

    print(f"\nDiagonal V_ee (3 closed-shell-like configs):")
    print(f"  E_A (1s^2 2s^2)              = {E_A:.6f}")
    print(f"  E_B (1s^2 2p0^2)             = {E_B:.6f}")
    print(f"  E_C (1s^2 (2p_-1 2p_+1)_S=0) = {E_C:.6f}")

    # Off-diagonal couplings.
    V_AB = two_e_int_singlet_pair(s2, s2, p_0, p_0, Z)
    V_AC = np.sqrt(2.0) * two_e_int_singlet_pair(s2, s2, p_m, p_p, Z)
    V_BC = np.sqrt(2.0) * two_e_int_singlet_pair(p_0, p_0, p_m, p_p, Z)

    print(f"\nOff-diagonal V_ee:")
    print(f"  V_AB (2s^2 <-> 2p0^2)        = {V_AB:.6f}")
    print(f"  V_AC (2s^2 <-> 2p_-1 2p_+1)  = {V_AC:.6f}")
    print(f"  V_BC (2p0^2 <-> 2p_-1 2p_+1) = {V_BC:.6f}")

    H_eff = np.array([
        [E_A, V_AB, V_AC],
        [V_AB, E_B, V_BC],
        [V_AC, V_BC, E_C],
    ])
    eigs, vecs = np.linalg.eigh(H_eff)

    print(f"\n3x3 V_ee-corrected effective H eigenvalues:")
    for i, e in enumerate(eigs):
        print(f"  {i}: {e:.6f}  + E_kin = {E_kin + e:.6f}")
    gs = vecs[:, 0]
    print(f"\nGS coefficients (basis: 2s^2, 2p0^2, (2p_-1 2p_+1)_S=0):")
    print(f"  c_A = {gs[0]:.4f}  c_B = {gs[1]:.4f}  c_C = {gs[2]:.4f}")

    occ_diag = np.array([
        2.0,
        2 * gs[0]**2,         # 2s
        gs[2]**2,             # 2p_-1
        2 * gs[1]**2,         # 2p_0
        gs[2]**2,             # 2p_+1
    ])
    print(f"\nSpatial 1-RDM diagonal occupations "
          f"(1s, 2s, 2p_-1, 2p_0, 2p_+1):")
    for orb, n in zip(['1s', '2s', '2p_-1', '2p_0', '2p_+1'], occ_diag):
        print(f"  {orb}: {n:.6f}")
    print(f"  trace = {occ_diag.sum():.6f}")

    p = occ_diag / occ_diag.sum()
    p_pos = p[p > 1e-14]
    S_full = float(-np.sum(p_pos * np.log(p_pos)))

    occ_kin_diag = np.array([2.0, 2/3, 1/3, 2/3, 1/3])
    p_kin = occ_kin_diag / occ_kin_diag.sum()
    p_pos_kin = p_kin[p_kin > 1e-14]
    S_kin = float(-np.sum(p_pos_kin * np.log(p_pos_kin)))

    print(f"\nvon Neumann entropy (normalized 1-RDM):")
    print(f"  S_full (analytical) = {S_full:.6f}")
    print(f"  S_kin  (degen avg)  = {S_kin:.6f}")
    print(f"  ratio S_full/S_kin  = {S_full/S_kin:.4f}")

    out = {
        'Z': Z,
        'E_kin': float(E_kin),
        'diagonal_V': {'E_A': float(E_A), 'E_B': float(E_B), 'E_C': float(E_C)},
        'off_diagonal_V': {'V_AB': float(V_AB),
                           'V_AC': float(V_AC),
                           'V_BC': float(V_BC)},
        'H_eff_eigenvalues': [float(e) for e in eigs],
        'GS_coefficients': {'c_A': float(gs[0]),
                            'c_B': float(gs[1]),
                            'c_C': float(gs[2])},
        'occ_diagonal': occ_diag.tolist(),
        'S_full_analytical': S_full,
        'S_kin_analytical': S_kin,
    }
    out_path = os.path.join(os.path.dirname(__file__), 'data',
                            'ep2N_be_analytical.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved -> {out_path}")


if __name__ == '__main__':
    main()
