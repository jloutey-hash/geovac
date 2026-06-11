"""G6-Full: Project onto (1,1) graviton irrep, then compute S^(2).

Key idea: instead of diagonalizing S^(2) on the full 1-form space (dim ~1000),
project onto the (J_L=1, J_R=1) angular momentum subspace first. This is the
graviton irrep. The projected subspace is much smaller, and its S^(2) eigenvalues
should give the Lichnerowicz spectrum.

A Hermitian matrix V_{ij} carries SO(4) angular momentum from both bra and ket:
  (j_L^i, j_R^i) x (j_L^j, j_R^j) -> various (J_L, J_R)

The (1,1) projection picks out the specific component via CG coefficients:
  V^{(1,M_L; 1,M_R)} = sum_{i,j} CG_L * CG_R * V_{ij}

We build a basis for the (1,1) subspace of Herm(H), compute S^(2) on it,
and diagonalize.
"""

import numpy as np
from sympy.physics.wigner import clebsch_gordan
from sympy import N as sympy_N, Rational
import json
from pathlib import Path
import time

OUT_DIR = Path(__file__).parent / "data"
OUT_DIR.mkdir(exist_ok=True)


def half_int_range(j):
    m = -j
    result = []
    while m <= j + 1e-10:
        result.append(round(2*m)/2)
        m += 1.0
    return result


def cg(j1, m1, j2, m2, J, M):
    if abs(m1 + m2 - M) > 1e-10:
        return 0.0
    if J < abs(j1 - j2) - 1e-10 or J > j1 + j2 + 1e-10:
        return 0.0
    if abs(m1) > j1 + 1e-10 or abs(m2) > j2 + 1e-10 or abs(M) > J + 1e-10:
        return 0.0
    val = clebsch_gordan(
        Rational(int(2*j1), 2), Rational(int(2*j2), 2), Rational(int(2*J), 2),
        Rational(int(2*m1), 2), Rational(int(2*m2), 2), Rational(int(2*M), 2))
    return float(sympy_N(val))


class CHBasis:
    def __init__(self, n_max):
        self.n_max = n_max
        self.states = []
        for n in range(n_max + 1):
            for chi in [+1, -1]:
                j_L = (n + 1) / 2.0 if chi == +1 else n / 2.0
                j_R = n / 2.0 if chi == +1 else (n + 1) / 2.0
                for m_L in half_int_range(j_L):
                    for m_R in half_int_range(j_R):
                        self.states.append((n, chi, j_L, j_R, m_L, m_R))
        self.dim = len(self.states)
        self.D0 = np.array([chi * (n + 1.5) for (n, chi, *_) in self.states])


def build_graviton_basis(basis, J_L_target=1.0, J_R_target=1.0):
    """Build a basis for the (J_L, J_R) subspace of Herm(H).

    For each pair (i, j) of states, the matrix element e_{ij} = |i><j| + |j><i|
    carries angular momentum from the CG decomposition:
      (j_L^i, m_L^i) x (j_L^j, m_L^j) -> J_L with M_L = m_L^i - m_L^j

    We project onto the (J_L_target, J_R_target) component using CG coefficients.

    Specifically, we build basis matrices B^{M_L, M_R} where:
      B^{M_L, M_R}_{ij} = CG(j_L^j, m_L^j, J_L, M_L | j_L^i, m_L^i)
                         * CG(j_R^j, m_R^j, J_R, M_R | j_R^i, m_R^i)

    These are the "angular momentum states" in matrix space.
    """
    dim = basis.dim
    states = basis.states
    JL = J_L_target
    JR = J_R_target

    # For each (M_L, M_R) and each PAIR of sector types, build one basis element
    # Total: (2*JL+1) * (2*JR+1) * (number of sector pairs with this coupling)
    projected_matrices = []
    labels = []

    # Group states by sector (n, chi)
    sectors = {}
    for idx, (n, chi, jL, jR, mL, mR) in enumerate(states):
        key = (n, chi, jL, jR)
        if key not in sectors:
            sectors[key] = []
        sectors[key].append(idx)

    sector_keys = sorted(sectors.keys())
    print(f"    Sectors: {len(sector_keys)}")

    # For each pair of sectors (a, b) with a <= b:
    for a_idx, key_a in enumerate(sector_keys):
        na, chia, jLa, jRa = key_a
        for b_idx, key_b in enumerate(sector_keys):
            if b_idx < a_idx:
                continue
            nb, chib, jLb, jRb = key_b

            # Check if this sector pair can produce (JL, JR)
            if JL < abs(jLa - jLb) - 1e-10 or JL > jLa + jLb + 1e-10:
                continue
            if JR < abs(jRa - jRb) - 1e-10 or JR > jRa + jRb + 1e-10:
                continue

            # For each (M_L, M_R), build the projected matrix
            for M_L in half_int_range(JL):
                for M_R in half_int_range(JR):
                    V = np.zeros((dim, dim))
                    for i_local in sectors[key_a]:
                        ni, ci, jLi, jRi, mLi, mRi = states[i_local]
                        for j_local in sectors[key_b]:
                            nj, cj, jLj, jRj, mLj, mRj = states[j_local]

                            # CG: <j_L^i, m_L^i | j_L^j, m_L^j; J_L, M_L>
                            # i.e., j_L^j x J_L -> j_L^i
                            c_L = cg(jLj, mLj, JL, M_L, jLi, mLi)
                            c_R = cg(jRj, mRj, JR, M_R, jRi, mRi)

                            if abs(c_L * c_R) < 1e-15:
                                continue

                            V[i_local, j_local] += c_L * c_R

                    # Hermitianize (for a != b sector pairs)
                    if a_idx != b_idx:
                        V_herm = (V + V.T) / 2.0
                        V_anti = (V - V.T) / 2.0
                        if np.max(np.abs(V_herm)) > 1e-14:
                            projected_matrices.append(V_herm / np.linalg.norm(V_herm))
                            labels.append((key_a, key_b, M_L, M_R, 'sym'))
                        if np.max(np.abs(V_anti)) > 1e-14:
                            projected_matrices.append(V_anti / np.linalg.norm(V_anti))
                            labels.append((key_a, key_b, M_L, M_R, 'anti'))
                    else:
                        if np.max(np.abs(V)) > 1e-14:
                            V_n = V / np.linalg.norm(V)
                            projected_matrices.append(V_n)
                            labels.append((key_a, key_b, M_L, M_R, 'diag'))

    print(f"    Raw (J_L={JL}, J_R={JR})-projected matrices: {len(projected_matrices)}")

    # Gram-Schmidt for linear independence
    flat = [v.flatten() for v in projected_matrices]
    basis_vecs = []
    basis_labels = []
    for i, v in enumerate(flat):
        w = v.copy()
        for b in basis_vecs:
            w -= np.dot(w, b) * b
        n = np.linalg.norm(w)
        if n > 1e-10:
            basis_vecs.append(w / n)
            basis_labels.append(labels[i])

    result = [v.reshape((dim, dim)) for v in basis_vecs]
    print(f"    Independent (J_L={JL}, J_R={JR}) basis: dim = {len(result)}")
    return result, basis_labels


def compute_S2(D0, V, Lambda_sq, eps=1e-4):
    """S^(2)[V] via 5-point central difference."""
    D0_mat = np.diag(D0)
    def S(e):
        D = D0_mat + e * V
        eigs = np.linalg.eigvalsh(D)
        return np.sum(np.exp(-eigs**2 / Lambda_sq))
    return (-S(2*eps) + 16*S(eps) - 30*S(0) + 16*S(-eps) - S(-2*eps)) / (12*eps**2)


def compute_projected_quadratic_form(basis, grav_basis, Lambda_sq):
    """Compute S^(2) quadratic form on the graviton-projected subspace."""
    D0 = basis.D0
    n = len(grav_basis)
    K = np.zeros((n, n))

    for i in range(n):
        K[i, i] = compute_S2(D0, grav_basis[i], Lambda_sq)
        for j in range(i+1, n):
            s2_sum = compute_S2(D0, grav_basis[i] + grav_basis[j], Lambda_sq)
            K[i, j] = 0.5 * (s2_sum - K[i, i] - compute_S2(D0, grav_basis[j], Lambda_sq))
            K[j, i] = K[i, j]
        if (i+1) % 10 == 0:
            print(f"      row {i+1}/{n}")

    return K


def main():
    print("="*72)
    print("G6-Full: (1,1) graviton irrep projection + S^(2) eigenvalues")
    print("="*72)

    Lambda_sq = 6.0
    results = {}

    for n_max in [1, 2, 3]:
        print(f"\n{'='*72}")
        print(f"n_max = {n_max}")
        print(f"{'='*72}")

        t0 = time.time()
        basis = CHBasis(n_max)
        print(f"  dim_H = {basis.dim}")

        # Build (1,1) graviton basis
        print(f"\n  Building (1,1) graviton-projected basis...")
        grav_basis, grav_labels = build_graviton_basis(basis)

        if len(grav_basis) == 0:
            print("  No (1,1) content found!")
            continue

        # Compute quadratic form on graviton subspace
        n_grav = len(grav_basis)
        print(f"\n  Computing S^(2) on {n_grav}-dim graviton subspace...")
        K = compute_projected_quadratic_form(basis, grav_basis, Lambda_sq)

        # Diagonalize
        eigs = np.linalg.eigvalsh(K)
        eigs_sorted = np.sort(eigs)[::-1]

        print(f"\n  Graviton S^(2) eigenvalues:")
        print(f"    Positive: {np.sum(eigs > 1e-10)}")
        print(f"    Zero (|e|<1e-10): {np.sum(np.abs(eigs) < 1e-10)}")
        print(f"    Negative: {np.sum(eigs < -1e-10)}")

        # Show all non-zero eigenvalues with multiplicities
        clusters = []
        sorted_e = np.sort(eigs)[::-1]
        i = 0
        while i < len(sorted_e):
            val = sorted_e[i]
            count = 1
            while i + count < len(sorted_e) and abs(sorted_e[i+count] - val) < 0.0001:
                count += 1
            if abs(val) > 1e-10:
                clusters.append((val, count))
            i += count

        print(f"\n  Non-zero eigenvalue clusters (value, multiplicity):")
        for val, mult in clusters:
            print(f"    {val:+.8f}  x{mult}")

        # Ratios of positive eigenvalues
        pos_clusters = [(v, m) for v, m in clusters if v > 0.001]
        if len(pos_clusters) >= 2:
            min_pos = min(v for v, m in pos_clusters)
            print(f"\n  POSITIVE eigenvalue ratios (vs smallest={min_pos:.6f}):")
            for v, m in sorted(pos_clusters):
                print(f"    ratio = {v/min_pos:.4f}  (value={v:.6f}, mult={m})")

            # Compare to Lichnerowicz: k(k+2)-2 for k>=2: 6, 13, 22, 33
            print(f"\n  Lichnerowicz reference ratios:")
            lich = [k*(k+2)-2 for k in range(2, 8)]
            for l in lich:
                print(f"    l(l+2)-2 = {l}, ratio = {l/lich[0]:.4f}")

        elapsed = time.time() - t0
        print(f"\n  Time: {elapsed:.1f}s")

        results[f"n_max={n_max}"] = {
            "dim_H": basis.dim,
            "dim_graviton": n_grav,
            "eigenvalues": eigs_sorted.tolist(),
            "clusters": [(float(v), int(m)) for v, m in clusters],
            "time_s": elapsed
        }

    # Save
    out_path = OUT_DIR / "g6_graviton_projected.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
