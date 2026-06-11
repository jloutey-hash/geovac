"""Sprint G6-Full — graviton Fierz-Pauli from discrete spectral action.

Strategy:
---------
On the GeoVac discrete substrate, the spectral action's second variation is a
finite-dimensional quadratic form. The graviton (Lichnerowicz) spectrum should
emerge algebraically from:

1. Build the CH Dirac D_0 at finite n_max with explicit SO(4) state labels
2. Build multiplication operators M_Y for scalar harmonics (CG coefficients)
3. Construct the space of NCG 1-forms Omega^1_D = span{a[D,b]}
4. Compute S^(2) on Omega^1_D via finite-difference of Tr(e^{-(D_0+eps*A)^2/L^2})
5. Diagonalize and classify eigenmodes by SO(4) content
6. Check: does the (j_L,j_R)=(1,1) sector reproduce Lichnerowicz eigenvalues?

The Lichnerowicz Laplacian on TT rank-2 tensors on unit S^3:
  Delta_L T^{(l)} = [l(l+2) - 2] T^{(l)},  l >= 2
  Eigenvalues: 6, 12, 20, 30, ... (integers!)

If these emerge from the discrete computation, we have bit-exact graviton
dynamics from the algebraic framework.
"""

import json
from pathlib import Path
from itertools import product as iterproduct

import numpy as np
from scipy.linalg import eigh

# Wigner 3j symbols
try:
    from sympy.physics.wigner import wigner_3j, clebsch_gordan
    from sympy import N as sympy_N, Rational
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False

OUT_DIR = Path(__file__).parent / "data"
OUT_DIR.mkdir(exist_ok=True)


# =============================================================================
# Step 1: CH Hilbert space with explicit SO(4) state labels
# =============================================================================

def build_ch_hilbert_space(n_max):
    """Build the Camporesi-Higuchi Hilbert space at truncation n_max.

    Each state is labeled by (n, chirality, m_L, m_R) where:
      - n = 0, 1, ..., n_max (CH level)
      - chirality = +1 or -1
      - For chirality +1: j_L = (n+1)/2, j_R = n/2
      - For chirality -1: j_L = n/2, j_R = (n+1)/2
      - m_L in {-j_L, -j_L+1, ..., j_L}
      - m_R in {-j_R, -j_R+1, ..., j_R}

    Returns:
      states: list of (n, chi, j_L, j_R, m_L, m_R) tuples
      D0: diagonal array of CH eigenvalues lambda_i = chi * (n + 3/2)
    """
    states = []
    for n in range(n_max + 1):
        for chi in [+1, -1]:
            if chi == +1:
                j_L = (n + 1) / 2.0
                j_R = n / 2.0
            else:
                j_L = n / 2.0
                j_R = (n + 1) / 2.0
            for m_L in _half_int_range(j_L):
                for m_R in _half_int_range(j_R):
                    states.append((n, chi, j_L, j_R, m_L, m_R))

    dim = len(states)
    D0 = np.zeros(dim)
    for i, (n, chi, j_L, j_R, m_L, m_R) in enumerate(states):
        D0[i] = chi * (n + 1.5)

    return states, D0


def _half_int_range(j):
    """Generate m values: -j, -j+1, ..., j."""
    m = -j
    result = []
    while m <= j + 1e-10:
        result.append(m)
        m += 1.0
    return result


# =============================================================================
# Step 2: Multiplication operators via CG coefficients
# =============================================================================

def cg_float(j1, m1, j2, m2, J, M):
    """Clebsch-Gordan coefficient <j1 m1; j2 m2 | J M> as float."""
    if not HAS_SYMPY:
        raise RuntimeError("sympy required for CG coefficients")
    if abs(m1 + m2 - M) > 1e-10:
        return 0.0
    if J < abs(j1 - j2) - 1e-10 or J > j1 + j2 + 1e-10:
        return 0.0
    if abs(m1) > j1 + 1e-10 or abs(m2) > j2 + 1e-10 or abs(M) > J + 1e-10:
        return 0.0
    # Convert to Rational for sympy
    j1r = Rational(int(2*j1), 2)
    m1r = Rational(int(2*m1), 2)
    j2r = Rational(int(2*j2), 2)
    m2r = Rational(int(2*m2), 2)
    Jr = Rational(int(2*J), 2)
    Mr = Rational(int(2*M), 2)
    val = clebsch_gordan(j1r, j2r, Jr, m1r, m2r, Mr)
    return float(sympy_N(val))


def build_scalar_multiplication_operator(states, l, m_L_Y, m_R_Y):
    """Build the matrix M_Y for scalar harmonic Y_{l, m_L_Y, m_R_Y}.

    The scalar harmonic at level l transforms as (l/2, l/2) under SO(4).

    Matrix element <psi_f | M_Y | psi_i> is the triple overlap integral
    on S^3 = SU(2), which factorizes into left and right CG coefficients:

      <psi_f | Y | psi_i> = C^{j_L^f}_{m_L^f; l/2, m_L_Y; j_L^i, m_L^i}
                           * C^{j_R^f}_{m_R^f; l/2, m_R_Y; j_R^i, m_R^i}
                           * (normalization)

    The normalization involves the dimensions and is:
      N = sqrt((2*j_L^i+1)(2*j_R^i+1)(l+1)^2 / ((2*j_L^f+1)(2*j_R^f+1)))

    (This comes from the Wigner-Eckart theorem on SU(2)xSU(2).)
    """
    dim = len(states)
    M = np.zeros((dim, dim))
    j_Y_L = l / 2.0
    j_Y_R = l / 2.0

    for f_idx, (n_f, chi_f, jLf, jRf, mLf, mRf) in enumerate(states):
        for i_idx, (n_i, chi_i, jLi, jRi, mLi, mRi) in enumerate(states):
            # CG selection: m_L^f = m_L^i + m_L_Y, m_R^f = m_R^i + m_R_Y
            if abs(mLf - mLi - m_L_Y) > 1e-10:
                continue
            if abs(mRf - mRi - m_R_Y) > 1e-10:
                continue

            # Triangle inequalities
            if jLf < abs(jLi - j_Y_L) - 1e-10 or jLf > jLi + j_Y_L + 1e-10:
                continue
            if jRf < abs(jRi - j_Y_R) - 1e-10 or jRf > jRi + j_Y_R + 1e-10:
                continue

            # CG coefficients (product of left and right sectors)
            cg_L = cg_float(jLi, mLi, j_Y_L, m_L_Y, jLf, mLf)
            cg_R = cg_float(jRi, mRi, j_Y_R, m_R_Y, jRf, mRf)

            if abs(cg_L) < 1e-15 or abs(cg_R) < 1e-15:
                continue

            # Normalization from Wigner-Eckart on SU(2)xSU(2)
            # The reduced matrix element includes dimension factors
            norm = np.sqrt((2*jLi+1) * (2*jRi+1) * (l+1)**2
                          / max((2*jLf+1) * (2*jRf+1), 1e-30))

            M[f_idx, i_idx] = cg_L * cg_R * norm

    return M


# =============================================================================
# Step 3: NCG 1-forms and the graviton subspace
# =============================================================================

def build_one_form(D0, M_a, M_b):
    """Build 1-form A = M_a [D, M_b] = M_a (D M_b - M_b D).

    Since D0 is diagonal: [D, M_b]_{ij} = (lambda_i - lambda_j) * (M_b)_{ij}
    Then A = M_a @ [D, M_b].
    """
    dim = len(D0)
    # [D, M_b]_{ij} = (D0[i] - D0[j]) * M_b[i,j]
    comm_D_b = np.outer(D0, np.ones(dim)) - np.outer(np.ones(dim), D0)
    comm_D_b = comm_D_b * M_b
    A = M_a @ comm_D_b
    return A


def build_one_form_basis(states, D0, l_max_algebra):
    """Build a basis for the space of 1-forms Omega^1_D at given n_max.

    Uses scalar harmonics Y_{l, m_L, m_R} with l = 1, ..., l_max_algebra
    as algebra elements. Constructs A = a[D,b] for all pairs.

    Returns a list of Hermitian 1-forms (A + A^dagger)/2 and (A - A^dagger)/(2i).
    """
    dim = len(states)
    one_forms = []

    # Collect scalar harmonics that give non-trivial [D, Y]
    scalar_harmonics = []
    for l in range(1, l_max_algebra + 1):
        j_Y = l / 2.0
        for m_L in _half_int_range(j_Y):
            for m_R in _half_int_range(j_Y):
                scalar_harmonics.append((l, m_L, m_R))

    print(f"  Building {len(scalar_harmonics)} scalar harmonic operators...")
    M_ops = {}
    for (l, mL, mR) in scalar_harmonics:
        M = build_scalar_multiplication_operator(states, l, mL, mR)
        if np.max(np.abs(M)) > 1e-14:
            M_ops[(l, mL, mR)] = M

    print(f"  {len(M_ops)} non-trivial operators. Building 1-forms...")

    # Build 1-forms A = a[D,b] for all pairs (a, b) from the harmonics
    # Include the identity (l=0) as an algebra element
    all_ops = [(0, 0, 0, np.eye(dim))]  # identity
    for key, M in M_ops.items():
        all_ops.append((*key, M))

    raw_forms = []
    for i, (l_a, mL_a, mR_a, M_a) in enumerate(all_ops):
        for j, (l_b, mL_b, mR_b, M_b) in enumerate(all_ops):
            if l_b == 0:
                continue  # [D, 1] = 0
            A = build_one_form(D0, M_a, M_b)
            if np.max(np.abs(A)) > 1e-14:
                raw_forms.append(A)

    print(f"  {len(raw_forms)} raw 1-forms. Extracting Hermitian basis...")

    # Extract linearly independent Hermitian basis
    # Hermitianize: take real and imaginary parts
    herm_forms = []
    for A in raw_forms:
        # Hermitian part
        H_re = (A + A.conj().T) / 2.0
        H_im = (A - A.conj().T) / (2.0j)
        if np.max(np.abs(H_re)) > 1e-14:
            herm_forms.append(H_re.real)
        if np.max(np.abs(H_im)) > 1e-14:
            herm_forms.append(H_im.real)

    # Gram-Schmidt to find linearly independent set
    basis = _gram_schmidt_real(herm_forms, tol=1e-10)
    print(f"  Linearly independent Hermitian 1-form basis: dim = {len(basis)}")

    return basis


def _gram_schmidt_real(vectors, tol=1e-10):
    """Gram-Schmidt on a list of real matrices (flattened to vectors)."""
    if not vectors:
        return []
    dim = vectors[0].shape
    flat = [v.flatten() for v in vectors]
    basis_flat = []
    basis_mat = []

    for v in flat:
        # Project out existing basis
        w = v.copy()
        for b in basis_flat:
            w -= np.dot(w, b) * b
        norm = np.linalg.norm(w)
        if norm > tol:
            w /= norm
            basis_flat.append(w)
            basis_mat.append(w.reshape(dim))

    return basis_mat


# =============================================================================
# Step 4: Spectral action second variation via finite difference
# =============================================================================

def spectral_action(D_matrix, Lambda_sq):
    """Compute Tr(e^{-D^2/Lambda^2}) for a Hermitian matrix D."""
    eigenvalues = np.linalg.eigvalsh(D_matrix)
    return np.sum(np.exp(-eigenvalues**2 / Lambda_sq))


def second_variation(D0_diag, V, Lambda_sq, eps=1e-4):
    """Compute S^(2)[V] = d^2/deps^2 Tr(e^{-(D0+eps*V)^2/L^2})|_{eps=0}.

    Uses 5-point stencil for better accuracy.
    """
    dim = len(D0_diag)
    D0_mat = np.diag(D0_diag)

    def S(epsilon):
        D = D0_mat + epsilon * V
        return spectral_action(D, Lambda_sq)

    # 5-point central difference: f''(0) = (-f(2h)+16f(h)-30f(0)+16f(-h)-f(-2h))/(12h^2)
    S_0 = S(0)
    S_p1 = S(eps)
    S_m1 = S(-eps)
    S_p2 = S(2*eps)
    S_m2 = S(-2*eps)

    d2S = (-S_p2 + 16*S_p1 - 30*S_0 + 16*S_m1 - S_m2) / (12 * eps**2)
    return d2S


def compute_quadratic_form_matrix(D0_diag, basis, Lambda_sq, eps=1e-4):
    """Compute the matrix K_{ab} = S^(2)[e_a, e_b] on the 1-form basis.

    K_{ab} = d^2/deps^2 Tr(e^{-(D0 + eps_a e_a + eps_b e_b)^2/L^2})

    By polarization: K_{ab} = (1/2)[S^(2)(e_a+e_b) - S^(2)(e_a) - S^(2)(e_b)]
    """
    n_basis = len(basis)
    K = np.zeros((n_basis, n_basis))

    # Diagonal elements
    for a in range(n_basis):
        K[a, a] = second_variation(D0_diag, basis[a], Lambda_sq, eps)

    # Off-diagonal by polarization
    for a in range(n_basis):
        for b in range(a+1, n_basis):
            S_ab = second_variation(D0_diag, basis[a] + basis[b], Lambda_sq, eps)
            K[a, b] = 0.5 * (S_ab - K[a, a] - K[b, b])
            K[b, a] = K[a, b]

    return K


# =============================================================================
# Step 5: SO(4) classification of eigenmodes
# =============================================================================

def classify_so4_content(V, states):
    """Classify a Hermitian matrix V by its SO(4) = SU(2)_L x SU(2)_R content.

    For each non-zero block V_{ij}, the coupling (j_L^i, j_R^i) x (j_L^j, j_R^j)
    produces SO(4) irreps. We identify the dominant (J_L, J_R) content.

    Returns dict of {(J_L, J_R): weight} where weight = sum of |V_{ij}|^2
    for entries contributing to that irrep.
    """
    dim = len(states)
    content = {}

    for i in range(dim):
        for j in range(dim):
            if abs(V[i, j]) < 1e-14:
                continue
            n_i, chi_i, jLi, jRi, mLi, mRi = states[i]
            n_j, chi_j, jLj, jRj, mLj, mRj = states[j]

            # The matrix element V_{ij} carries angular momentum
            # Delta_m_L = m_L^i - m_L^j, Delta_m_R = m_R^i - m_R^j
            # It contributes to irreps (J_L, J_R) where
            # |jLi - jLj| <= J_L <= jLi + jLj and similarly for R
            # The specific J is determined by the m-values

            # For classification, we track which (jLi x jLj, jRi x jRj) products
            # are present. The irrep content is the range of tensor products.
            w = abs(V[i, j])**2

            # All possible (J_L, J_R) in the tensor product
            for J_L in _half_int_range_between(abs(jLi - jLj), jLi + jLj):
                for J_R in _half_int_range_between(abs(jRi - jRj), jRi + jRj):
                    key = (J_L, J_R)
                    if key not in content:
                        content[key] = 0.0
                    # Weight is proportional to |V_{ij}|^2 divided among irreps
                    n_irreps = int((min(jLi+jLj, J_L+0.1) - max(abs(jLi-jLj), J_L-0.1)) < 0.5)
                    # Simple: just accumulate by the product range
                    content[key] += w

    # Normalize
    total = sum(content.values())
    if total > 0:
        content = {k: v/total for k, v in content.items()}

    return content


def _half_int_range_between(j_min, j_max):
    """Generate half-integer values from j_min to j_max in steps of 1."""
    j = j_min
    result = []
    while j <= j_max + 1e-10:
        result.append(round(2*j)/2)
        j += 1.0
    return result


# =============================================================================
# Step 6: Main computation
# =============================================================================

def analyze_graviton_spectrum(n_max, Lambda_sq=6.0, l_max_algebra=None):
    """Run the full G6 graviton analysis at given n_max.

    1. Build CH space
    2. Build 1-form basis from algebra elements
    3. Compute S^(2) quadratic form
    4. Diagonalize and classify
    5. Compare to Lichnerowicz eigenvalues
    """
    if l_max_algebra is None:
        l_max_algebra = 2 * n_max  # capture all relevant harmonics

    print(f"\n{'='*72}")
    print(f"G6-Full Graviton Analysis at n_max = {n_max}")
    print(f"{'='*72}")

    # Step 1: Build Hilbert space
    states, D0 = build_ch_hilbert_space(n_max)
    dim_H = len(states)
    print(f"\n[Step 1] CH Hilbert space: dim = {dim_H}")
    print(f"  Eigenvalues: {sorted(set(D0))}")

    # Step 2-3: Build 1-form basis
    print(f"\n[Step 2-3] Building 1-form basis (l_max_algebra = {l_max_algebra})...")
    basis = build_one_form_basis(states, D0, l_max_algebra)
    n_basis = len(basis)

    if n_basis == 0:
        print("  ERROR: no 1-forms found!")
        return None

    # Step 4: Compute quadratic form
    print(f"\n[Step 4] Computing S^(2) quadratic form ({n_basis} x {n_basis})...")
    K = compute_quadratic_form_matrix(D0, basis, Lambda_sq)

    # Diagonalize
    eigenvalues, eigenvectors = eigh(K)

    # Sort by absolute value (most physical modes have largest |eigenvalue|)
    idx = np.argsort(-np.abs(eigenvalues))
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    print(f"\n[Step 5] Eigenvalue spectrum of S^(2):")
    print(f"  Total modes: {n_basis}")
    print(f"  Positive eigenvalues: {np.sum(eigenvalues > 1e-10)}")
    print(f"  Zero eigenvalues (|eig| < 1e-10): {np.sum(np.abs(eigenvalues) < 1e-10)}")
    print(f"  Negative eigenvalues: {np.sum(eigenvalues < -1e-10)}")
    print(f"\n  Top 20 eigenvalues:")
    for i in range(min(20, n_basis)):
        print(f"    {i:3d}: {eigenvalues[i]:+.8f}")

    # Step 5-6: Classify eigenmodes and check Lichnerowicz
    print(f"\n[Step 6] Checking against Lichnerowicz spectrum...")
    print(f"  Expected TT eigenvalues: l(l+2)-2 for l>=2 → 6, 12, 20, 30, ...")

    # Look at ratios of positive eigenvalues
    pos_eigs = eigenvalues[eigenvalues > 1e-10]
    if len(pos_eigs) >= 2:
        # Normalize by the smallest positive eigenvalue
        eig_min = np.min(pos_eigs)
        ratios = pos_eigs / eig_min
        print(f"\n  Positive eigenvalue ratios (normalized to smallest):")
        for i, r in enumerate(ratios[:15]):
            # Check against Lichnerowicz: 6, 12, 20, 30 → ratios 1, 2, 10/3, 5
            print(f"    {i:3d}: ratio = {r:.6f}")

    # Check if eigenvalue ratios match l(l+2)-2 pattern
    lichnerowicz = [l*(l+2) - 2 for l in range(2, 10)]  # 6, 12, 20, 30, ...
    lich_ratios = [x / lichnerowicz[0] for x in lichnerowicz]  # 1, 2, 10/3, 5, ...

    print(f"\n  Lichnerowicz ratios: {lich_ratios[:6]}")

    results = {
        "n_max": n_max,
        "Lambda_sq": Lambda_sq,
        "l_max_algebra": l_max_algebra,
        "dim_H": dim_H,
        "dim_one_form_basis": n_basis,
        "eigenvalues": eigenvalues.tolist(),
        "n_positive": int(np.sum(eigenvalues > 1e-10)),
        "n_zero": int(np.sum(np.abs(eigenvalues) < 1e-10)),
        "n_negative": int(np.sum(eigenvalues < -1e-10)),
        "lichnerowicz_expected": lichnerowicz[:6],
    }

    return results


# =============================================================================
# Quick diagnostic: just the sector-level formula (from G6-Diag)
# vs full 1-form computation
# =============================================================================

def quick_sector_check(n_max=2, Lambda_sq=6.0):
    """Quick sanity check: verify the S^(2) finite-difference agrees with
    the analytical sector-level formula from G6-Diag for within-sector modes.
    """
    print(f"\n{'='*72}")
    print(f"Quick sector check: finite-difference vs analytical at n_max={n_max}")
    print(f"{'='*72}")

    states, D0 = build_ch_hilbert_space(n_max)
    dim = len(states)

    # Pick a within-sector perturbation: a random Hermitian matrix
    # supported only in sector n=1, chi=+1 (j_L=1, j_R=1/2, dim=6)
    sector_indices = [i for i, (n, chi, jL, jR, mL, mR) in enumerate(states)
                      if n == 1 and chi == 1]
    print(f"  Testing sector n=1, chi=+1: {len(sector_indices)} states")

    # Build a random within-sector Hermitian perturbation
    d_sec = len(sector_indices)
    rng = np.random.default_rng(42)
    V_block = rng.standard_normal((d_sec, d_sec))
    V_block = (V_block + V_block.T) / 2  # Hermitianize
    V_block /= np.linalg.norm(V_block)  # normalize

    V = np.zeros((dim, dim))
    for a, i in enumerate(sector_indices):
        for b, j in enumerate(sector_indices):
            V[i, j] = V_block[a, b]

    # Finite-difference S^(2)
    s2_fd = second_variation(D0, V, Lambda_sq, eps=1e-4)

    # Analytical: A_lambda = a * (4*lam^2/L^4 - 2/L^2) for within-sector
    lam = D0[sector_indices[0]]  # all same lambda in sector
    a = np.exp(-lam**2 / Lambda_sq)
    A_lam = a * (4*lam**2 / Lambda_sq**2 - 2/Lambda_sq)

    # S^(2) on V with norm 1 should be A_lambda * ||V||^2 = A_lambda
    # But we need to account for the trace normalization
    # Actually S^(2)[V] = A_lambda * Tr(V^2) for within-sector V
    tr_V2 = np.trace(V @ V)
    s2_analytical = A_lam * tr_V2

    print(f"  lambda = {lam:.4f}")
    print(f"  A_lambda (analytical) = {A_lam:.8f}")
    print(f"  Tr(V^2) = {tr_V2:.8f}")
    print(f"  S^(2) analytical = {s2_analytical:.8f}")
    print(f"  S^(2) finite-diff = {s2_fd:.8f}")
    print(f"  Relative error: {abs(s2_fd - s2_analytical)/abs(s2_analytical):.2e}")

    # Also test a cross-sector perturbation
    # Sector n=0,chi=+1 (dim 2) to sector n=1,chi=+1 (dim 6)
    sec_a = [i for i, (n, chi, jL, jR, mL, mR) in enumerate(states)
             if n == 0 and chi == 1]
    sec_b = [i for i, (n, chi, jL, jR, mL, mR) in enumerate(states)
             if n == 1 and chi == 1]

    if sec_a and sec_b:
        print(f"\n  Testing cross-sector (n=0,+) to (n=1,+):")
        da, db = len(sec_a), len(sec_b)
        V_cross_block = rng.standard_normal((da, db))
        V_cross_block /= np.linalg.norm(V_cross_block)

        V_cross = np.zeros((dim, dim))
        for a, i in enumerate(sec_a):
            for b, j in enumerate(sec_b):
                V_cross[i, j] = V_cross_block[a, b]
                V_cross[j, i] = V_cross_block[a, b]  # Hermitian

        s2_cross_fd = second_variation(D0, V_cross, Lambda_sq, eps=1e-4)

        lam_a = D0[sec_a[0]]
        lam_b = D0[sec_b[0]]
        a_a = np.exp(-lam_a**2 / Lambda_sq)
        a_b = np.exp(-lam_b**2 / Lambda_sq)
        # Cross-sector formula from G6-Diag
        B_ab = -(2/Lambda_sq) * (lam_a * a_a - lam_b * a_b) / (lam_a - lam_b)

        tr_Vcross2 = np.trace(V_cross @ V_cross)
        s2_cross_analytical = B_ab * tr_Vcross2

        print(f"  lam_a = {lam_a:.4f}, lam_b = {lam_b:.4f}")
        print(f"  B_ab (analytical) = {B_ab:.8f}")
        print(f"  S^(2) analytical = {s2_cross_analytical:.8f}")
        print(f"  S^(2) finite-diff = {s2_cross_fd:.8f}")
        print(f"  Relative error: {abs(s2_cross_fd - s2_cross_analytical)/abs(s2_cross_analytical):.2e}")

    return True


# =============================================================================

def main():
    print("=" * 72)
    print("Sprint G6-Full: Graviton Fierz-Pauli from Discrete Spectral Action")
    print("=" * 72)

    # First: validate the finite-difference S^(2) against the known sector formulas
    print("\n[VALIDATION] Checking finite-difference against G6-Diag analytical...")
    quick_sector_check(n_max=2)

    # Then: run the full graviton analysis
    # Start small (n_max=1) to verify infrastructure
    print("\n\n" + "=" * 72)
    print("[MAIN ANALYSIS]")

    results = {}
    for n_max in [1, 2]:
        r = analyze_graviton_spectrum(n_max, Lambda_sq=6.0)
        if r is not None:
            results[f"n_max={n_max}"] = r

    # Save results
    out_path = OUT_DIR / "g6_full_graviton.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
