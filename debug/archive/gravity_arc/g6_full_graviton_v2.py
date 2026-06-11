"""Sprint G6-Full v2 -- Direct TT graviton test on discrete substrate.

Key insight: the Lichnerowicz eigenvalue for a TT harmonic at level l emerges as
a WEIGHTED SUM of sector-level quadratic form coefficients:

  S^(2)_l = sum_{a,b} Q(lambda_a, lambda_b) * ||delta_D^(l)_{ab}||^2

where:
- Q(lambda_a, lambda_b) is the known analytical formula from G6-Diag
- ||delta_D^(l)_{ab}||^2 is the WEIGHT of the level-l TT harmonic in sector block (a,b)
- The weights come from CG coefficients for the spinor-tensor coupling on S^3

Strategy: build delta_D for specific TT harmonics using the NCG 1-form dictionary
(A = a[D,b] for scalar harmonics a, b), then compute S^(2) via finite difference.
Compare RATIOS of S^(2) across different TT levels to Lichnerowicz ratios.

Lichnerowicz on unit S^3 (rough Laplacian + curvature):
  eigenvalues = k(k+2) + 2 for k >= 2   [where k(k+2) is the scalar Laplacian eigenvalue]
  i.e., 10, 17, 26, 37, ...

Or if using -nabla^2 on TT tensors only (without curvature term):
  eigenvalues = k(k+2) - 2 for k >= 2
  i.e., 6, 13, 22, 33, ...

We test RATIOS, so the overall shift cancels.
"""

import json
from pathlib import Path
import numpy as np
from sympy.physics.wigner import clebsch_gordan
from sympy import N as sympy_N, Rational, sqrt as sym_sqrt

OUT_DIR = Path(__file__).parent / "data"
OUT_DIR.mkdir(exist_ok=True)


def half_int_range(j):
    """m-values: -j, -j+1, ..., j"""
    m = -j
    result = []
    while m <= j + 1e-10:
        result.append(round(2*m)/2)
        m += 1.0
    return result


def cg(j1, m1, j2, m2, J, M):
    """CG coefficient as float, with selection rule shortcuts."""
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


# =============================================================================
# Hilbert space and Dirac
# =============================================================================

class CHBasis:
    """Camporesi-Higuchi basis at truncation n_max."""

    def __init__(self, n_max):
        self.n_max = n_max
        self.states = []  # (n, chi, j_L, j_R, m_L, m_R)
        self.D0 = None
        self._build()

    def _build(self):
        for n in range(self.n_max + 1):
            for chi in [+1, -1]:
                j_L = (n + 1) / 2.0 if chi == +1 else n / 2.0
                j_R = n / 2.0 if chi == +1 else (n + 1) / 2.0
                for m_L in half_int_range(j_L):
                    for m_R in half_int_range(j_R):
                        self.states.append((n, chi, j_L, j_R, m_L, m_R))
        self.dim = len(self.states)
        self.D0 = np.array([chi * (n + 1.5) for (n, chi, *_) in self.states])

    def sector_indices(self, n, chi):
        """Indices of states in sector (n, chi)."""
        return [i for i, (n_, c_, *_) in enumerate(self.states) if n_ == n and c_ == chi]


# =============================================================================
# Multiplication operators (scalar harmonics on S^3)
# =============================================================================

def build_scalar_mult(basis, l, m_L_Y, m_R_Y):
    """Multiplication by scalar harmonic Y_{l, m_L_Y, m_R_Y}.

    The scalar harmonic at level l has SO(4) quantum numbers (l/2, l/2).
    Matrix element <f|Y|i> = CG_L(j_L^i, m_L^i; l/2, m_L_Y | j_L^f, m_L^f)
                            * CG_R(j_R^i, m_R^i; l/2, m_R_Y | j_R^f, m_R^f)
                            * reduced_matrix_element
    """
    dim = basis.dim
    M = np.zeros((dim, dim))
    jY = l / 2.0

    for f_idx, (nf, cf, jLf, jRf, mLf, mRf) in enumerate(basis.states):
        for i_idx, (ni, ci, jLi, jRi, mLi, mRi) in enumerate(basis.states):
            # m-selection
            if abs(mLf - mLi - m_L_Y) > 1e-10:
                continue
            if abs(mRf - mRi - m_R_Y) > 1e-10:
                continue
            # Triangle
            if jLf < abs(jLi - jY) - 1e-10 or jLf > jLi + jY + 1e-10:
                continue
            if jRf < abs(jRi - jY) - 1e-10 or jRf > jRi + jY + 1e-10:
                continue

            c_L = cg(jLi, mLi, jY, m_L_Y, jLf, mLf)
            c_R = cg(jRi, mRi, jY, m_R_Y, jRf, mRf)
            if abs(c_L) < 1e-15 or abs(c_R) < 1e-15:
                continue

            # Reduced matrix element (Wigner-Eckart normalization)
            rme = np.sqrt((2*jLi+1)*(2*jRi+1)/(max((2*jLf+1)*(2*jRf+1), 1e-30)))
            M[f_idx, i_idx] = c_L * c_R * rme

    return M


# =============================================================================
# NCG 1-forms and spectral action
# =============================================================================

def build_commutator_D_M(D0, M):
    """[D, M]_{ij} = (D0[i] - D0[j]) * M[i,j]"""
    dim = len(D0)
    diff = D0[:, None] - D0[None, :]
    return diff * M


def spectral_action_value(D_matrix, Lambda_sq):
    """Tr(exp(-D^2/Lambda^2))"""
    eigs = np.linalg.eigvalsh(D_matrix)
    return np.sum(np.exp(-eigs**2 / Lambda_sq))


def compute_S2(D0_diag, V, Lambda_sq, eps=1e-4):
    """S^(2)[V] via 5-point central difference."""
    D0_mat = np.diag(D0_diag)

    def S(e):
        return spectral_action_value(D0_mat + e * V, Lambda_sq)

    return (-S(2*eps) + 16*S(eps) - 30*S(0) + 16*S(-eps) - S(-2*eps)) / (12*eps**2)


# =============================================================================
# Build specific 1-forms and classify by SO(4) content
# =============================================================================

def build_all_one_forms(basis, l_max_alg):
    """Build all independent Hermitian 1-forms from scalar harmonics up to l_max_alg."""
    dim = basis.dim
    D0 = basis.D0

    # Build scalar harmonic multiplication operators
    print(f"    Building scalar mult operators (l=1..{l_max_alg})...")
    M_ops = {}
    for l in range(1, l_max_alg + 1):
        jY = l / 2.0
        for mL in half_int_range(jY):
            for mR in half_int_range(jY):
                M = build_scalar_mult(basis, l, mL, mR)
                if np.max(np.abs(M)) > 1e-14:
                    M_ops[(l, mL, mR)] = M

    print(f"    {len(M_ops)} non-trivial M_Y operators")

    # Build 1-forms: A = M_a @ [D, M_b] for all pairs
    # Include identity as 'a' element
    a_ops = [(None, np.eye(dim))] + [(k, M) for k, M in M_ops.items()]
    b_ops = [(k, M) for k, M in M_ops.items()]  # b=identity gives [D,1]=0

    raw_forms = []
    for _, M_a in a_ops:
        for _, M_b in b_ops:
            comm = build_commutator_D_M(D0, M_b)
            A = M_a @ comm
            if np.max(np.abs(A)) > 1e-14:
                raw_forms.append(A)

    print(f"    {len(raw_forms)} raw 1-forms")

    # Hermitianize and extract linearly independent basis
    herm_vecs = []
    for A in raw_forms:
        H_re = (A + A.T) / 2.0  # These are real matrices
        H_im = (A - A.T) / 2.0
        if np.max(np.abs(H_re)) > 1e-14:
            herm_vecs.append(H_re)
        if np.max(np.abs(H_im)) > 1e-14:
            herm_vecs.append(H_im)

    # Gram-Schmidt for linear independence
    flat = [v.flatten() for v in herm_vecs]
    basis_vecs = []
    for v in flat:
        w = v.copy()
        for b in basis_vecs:
            w -= np.dot(w, b) * b
        n = np.linalg.norm(w)
        if n > 1e-10:
            basis_vecs.append(w / n)

    result = [v.reshape((dim, dim)) for v in basis_vecs]
    print(f"    Independent Hermitian 1-form basis: dim = {len(result)}")
    return result


# =============================================================================
# Main analysis: compute S^(2) eigenspectrum on 1-form space
# =============================================================================

def run_analysis(n_max, Lambda_sq=6.0, l_max_alg=None):
    """Full G6 analysis at n_max."""
    if l_max_alg is None:
        l_max_alg = min(2 * n_max, 4)  # balance coverage vs compute time

    print(f"\n{'='*72}")
    print(f"G6-Full v2: n_max={n_max}, Lambda^2={Lambda_sq}, l_max_alg={l_max_alg}")
    print(f"{'='*72}")

    basis = CHBasis(n_max)
    print(f"  dim_H = {basis.dim}")
    print(f"  Sectors: {sorted(set(basis.D0))}")

    # Build 1-form basis
    print(f"\n  Building 1-form space...")
    one_forms = build_all_one_forms(basis, l_max_alg)
    n_forms = len(one_forms)

    if n_forms == 0:
        print("  No 1-forms found!")
        return None

    # Compute S^(2) on each basis element
    print(f"\n  Computing S^(2) diagonal elements ({n_forms} forms)...")
    s2_diag = np.zeros(n_forms)
    for i in range(n_forms):
        s2_diag[i] = compute_S2(basis.D0, one_forms[i], Lambda_sq)
        if (i+1) % 50 == 0:
            print(f"    {i+1}/{n_forms} done")

    # For the full quadratic form, we'd need off-diagonal too.
    # But first let's look at just the diagonal spectrum.
    print(f"\n  S^(2) diagonal statistics:")
    print(f"    min = {np.min(s2_diag):.8f}")
    print(f"    max = {np.max(s2_diag):.8f}")
    print(f"    mean = {np.mean(s2_diag):.8f}")

    # Cluster the diagonal values
    vals_sorted = np.sort(s2_diag)
    print(f"\n  Sorted S^(2) values (sample):")
    # Show unique clusters
    clusters = []
    current_cluster = [vals_sorted[0]]
    for v in vals_sorted[1:]:
        if abs(v - current_cluster[-1]) < 0.001:
            current_cluster.append(v)
        else:
            clusters.append((np.mean(current_cluster), len(current_cluster)))
            current_cluster = [v]
    clusters.append((np.mean(current_cluster), len(current_cluster)))

    print(f"\n  Eigenvalue clusters (value, multiplicity):")
    for val, mult in sorted(clusters):
        print(f"    {val:+.6f}  x{mult}")

    # Now compute the FULL quadratic form matrix for more precise eigenvalues
    if n_forms <= 300:
        print(f"\n  Computing full {n_forms}x{n_forms} quadratic form...")
        K = np.zeros((n_forms, n_forms))
        for i in range(n_forms):
            K[i, i] = s2_diag[i]
            for j in range(i+1, n_forms):
                # Polarization: K[i,j] = (S2(e_i+e_j) - S2(e_i) - S2(e_j))/2
                s2_sum = compute_S2(basis.D0, one_forms[i] + one_forms[j], Lambda_sq)
                K[i, j] = 0.5 * (s2_sum - K[i, i] - s2_diag[j])
                K[j, i] = K[i, j]
            if (i+1) % 20 == 0:
                print(f"    row {i+1}/{n_forms}")

        # Diagonalize
        eigs = np.linalg.eigvalsh(K)
        eigs_sorted = np.sort(eigs)[::-1]  # descending

        print(f"\n  Full quadratic form eigenvalues:")
        print(f"    Positive: {np.sum(eigs > 1e-10)}")
        print(f"    Zero: {np.sum(np.abs(eigs) < 1e-10)}")
        print(f"    Negative: {np.sum(eigs < -1e-10)}")

        # Cluster eigenvalues
        clusters_full = []
        current = [eigs_sorted[0]]
        for v in eigs_sorted[1:]:
            if abs(v - current[-1]) < 0.0005:
                current.append(v)
            else:
                clusters_full.append((np.mean(current), len(current)))
                current = [v]
        clusters_full.append((np.mean(current), len(current)))

        print(f"\n  Full QF eigenvalue clusters:")
        for val, mult in sorted(clusters_full, reverse=True):
            if abs(val) > 1e-8:
                print(f"    {val:+.8f}  x{mult}")

        # Check ratios of positive eigenvalues
        pos_clusters = [(v, m) for v, m in clusters_full if v > 0.001]
        if len(pos_clusters) >= 2:
            print(f"\n  Positive eigenvalue RATIOS (relative to smallest positive):")
            min_pos = min(v for v, m in pos_clusters)
            for v, m in sorted(pos_clusters):
                print(f"    {v/min_pos:.6f}  (value={v:.6f}, mult={m})")

        return {"n_max": n_max, "dim_H": basis.dim, "n_forms": n_forms,
                "eigenvalues": eigs_sorted.tolist(),
                "clusters": [(float(v), int(m)) for v, m in clusters_full]}

    return {"n_max": n_max, "dim_H": basis.dim, "n_forms": n_forms,
            "diagonal_clusters": [(float(v), int(m)) for v, m in clusters]}


if __name__ == "__main__":
    # Start with n_max=1 (quick validation), then n_max=2 (the real test)
    results = {}

    r1 = run_analysis(n_max=1, l_max_alg=2)
    if r1:
        results["n_max=1"] = r1

    print("\n\n" + "="*72)
    print("Now running n_max=2 (this may take a few minutes)...")
    r2 = run_analysis(n_max=2, l_max_alg=3)
    if r2:
        results["n_max=2"] = r2

    out_path = OUT_DIR / "g6_full_graviton_v2.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")
