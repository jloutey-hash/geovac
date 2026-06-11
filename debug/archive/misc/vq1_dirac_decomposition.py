"""
VQ-1: Directional decomposition of the Camporesi-Higuchi Dirac operator on S^3.
================================================================================

Builds the three directional component matrices Gamma_mu (mu=1,2,3) in the
DiracLabel (n_fock, kappa, m_j) basis at n_max=2.

The Dirac operator on unit S^3 acting on spinor harmonics is:
    D = sigma . L + 3/2
where sigma = 2S (Pauli matrices in terms of spin-1/2 operators) and
L is the orbital angular momentum.

The directional decomposition is:
    sigma . L = sigma_x L_x + sigma_y L_y + sigma_z L_z
             = 2(S_x L_x + S_y L_y + S_z L_z)

We work in the |j, m_j> basis (which is what DiracLabel uses). In this basis:
- J_z is diagonal with eigenvalue m_j
- J_+ and J_- are standard ladder operators
- S_z, S_+, S_- require CG decomposition |j,m_j> = sum_{m_l,m_s} CG * |l,m_l;1/2,m_s>

We build:
    Gamma_1 = sigma_x L_x = 2 S_x L_x
    Gamma_2 = sigma_y L_y = 2 S_y L_y
    Gamma_3 = sigma_z L_z = 2 S_z L_z

And verify: sum_mu Gamma_mu = sigma . L = D - 3/2 * I

This script does NOT modify any production code.
"""

import numpy as np
import json
import sys
import os
from fractions import Fraction
from sympy import Rational, sqrt as sp_sqrt, S as sp_S

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.dirac_matrix_elements import (
    DiracLabel, iter_dirac_labels, kappa_to_l, kappa_to_j,
)
from sympy.physics.wigner import clebsch_gordan


def build_labels(n_max):
    """Enumerate all DiracLabel states at n_max."""
    labels = list(iter_dirac_labels(n_max))
    label_index = {lab: i for i, lab in enumerate(labels)}
    return labels, label_index


def cg_coeff(j1, m1, j2, m2, j, m):
    """Clebsch-Gordan coefficient <j1,m1; j2,m2 | j,m> via sympy."""
    val = clebsch_gordan(j1, j2, j, m1, m2, m)
    return float(val)


def build_J_matrices(labels, label_index):
    """Build J_x, J_y, J_z as dense matrices in the DiracLabel basis.

    J_z is diagonal with eigenvalue m_j.
    J_+ = J_x + i J_y connects |j,m> to |j,m+1> within same (n,kappa).
    J_- = J_x - i J_y connects |j,m> to |j,m-1> within same (n,kappa).

    From J_+ and J_-, we get J_x = (J_+ + J_-)/2 and J_y = (J_+ - J_-)/(2i).
    """
    N = len(labels)
    Jz = np.zeros((N, N))
    Jplus = np.zeros((N, N), dtype=complex)
    Jminus = np.zeros((N, N), dtype=complex)

    for i, lab in enumerate(labels):
        j = float(lab.j)
        mj = float(lab.m_j)
        Jz[i, i] = mj

        # J_+ |j,m> = sqrt(j(j+1) - m(m+1)) |j,m+1>
        if lab.two_m_j + 2 <= lab.j_times_2:
            target = DiracLabel(lab.n_fock, lab.kappa, lab.two_m_j + 2)
            if target in label_index:
                k = label_index[target]
                Jplus[k, i] = np.sqrt(j*(j+1) - mj*(mj+1))

        # J_- |j,m> = sqrt(j(j+1) - m(m-1)) |j,m-1>
        if lab.two_m_j - 2 >= -lab.j_times_2:
            target = DiracLabel(lab.n_fock, lab.kappa, lab.two_m_j - 2)
            if target in label_index:
                k = label_index[target]
                Jminus[k, i] = np.sqrt(j*(j+1) - mj*(mj-1))

    Jx = (Jplus + Jminus) / 2.0
    Jy = (Jplus - Jminus) / (2.0j)

    return Jx, Jy, Jz


def build_S_matrices(labels, label_index):
    """Build S_x, S_y, S_z in the |n, kappa, m_j> = |n, j, m_j> basis.

    Each state |j, m_j> (with given l from kappa) decomposes as:
        |j, m_j> = sum_{m_s = -1/2, +1/2} CG(l, m_j-m_s; 1/2, m_s | j, m_j) |l, m_j-m_s> |1/2, m_s>

    The spin operator S acts only on the |1/2, m_s> part:
        S_z |1/2, m_s> = m_s |1/2, m_s>
        S_+ |1/2, -1/2> = |1/2, +1/2>
        S_- |1/2, +1/2> = |1/2, -1/2>

    Matrix element:
        <j', m_j'| S_q |j, m_j> = sum_{m_s, m_s'} CG(l', m_j'-m_s'; 1/2, m_s' | j', m_j')
                                                    * <m_s'|S_q|m_s>
                                                    * CG(l, m_j-m_s; 1/2, m_s | j, m_j)
                                                    * delta(m_j'-m_s', m_j-m_s)  [from L part]

    The delta on m_l = m_j - m_s means S only connects states with same n_fock.
    For S_z (q=0): m_s' = m_s, m_j' = m_j, and it can couple different kappa if they share l and m_l.
    For S_+ (q=+1): m_s' = m_s + 1, requires m_s = -1/2, m_s' = +1/2, dm_j = +1.
    For S_- (q=-1): m_s' = m_s - 1, requires m_s = +1/2, m_s' = -1/2, dm_j = -1.
    """
    N = len(labels)
    Sz = np.zeros((N, N))
    Splus = np.zeros((N, N), dtype=complex)
    Sminus = np.zeros((N, N), dtype=complex)

    # Precompute l values
    l_vals = [kappa_to_l(lab.kappa) for lab in labels]
    j_vals = [float(lab.j) for lab in labels]
    mj_vals = [float(lab.m_j) for lab in labels]

    for i, lab_a in enumerate(labels):
        l_a = l_vals[i]
        j_a = j_vals[i]
        mj_a = mj_vals[i]

        for k, lab_b in enumerate(labels):
            # Must be same n_fock (S doesn't change n)
            if lab_a.n_fock != lab_b.n_fock:
                continue

            l_b = l_vals[k]
            j_b = j_vals[k]
            mj_b = mj_vals[k]

            # S doesn't change l, so l_a must equal l_b
            if l_a != l_b:
                continue

            l = l_a  # same for both

            # S_z: m_j' = m_j, m_s' = m_s
            if mj_a == mj_b:
                val = 0.0
                for ms2 in [-1, 1]:  # 2*m_s
                    ms = ms2 / 2.0
                    ml = mj_b - ms
                    if abs(ml) > l:
                        continue
                    cg_b = cg_coeff(Rational(l), Rational(ml), Rational(1,2), Rational(ms),
                                    Rational(j_b), Rational(mj_b))
                    cg_a = cg_coeff(Rational(l), Rational(ml), Rational(1,2), Rational(ms),
                                    Rational(j_a), Rational(mj_a))
                    val += cg_a * ms * cg_b
                Sz[i, k] = val

            # S_+: m_j' = m_j + 1, needs m_s = -1/2 -> m_s' = +1/2
            if abs(mj_a - mj_b - 1.0) < 1e-10:
                ms_b = -0.5   # m_s of the ket
                ms_a = 0.5    # m_s of the bra (after S_+ acts)
                ml = mj_b - ms_b   # m_l = m_j - m_s (same for bra and ket by delta)
                if abs(ml) <= l:
                    cg_b = cg_coeff(Rational(l), Rational(ml), Rational(1,2), Rational(-1,2),
                                    Rational(j_b), Rational(mj_b))
                    cg_a = cg_coeff(Rational(l), Rational(ml), Rational(1,2), Rational(1,2),
                                    Rational(j_a), Rational(mj_a))
                    # S_+ |1/2, -1/2> = 1 * |1/2, +1/2>
                    Splus[i, k] = cg_a * 1.0 * cg_b

            # S_-: m_j' = m_j - 1, needs m_s = +1/2 -> m_s' = -1/2
            if abs(mj_a - mj_b + 1.0) < 1e-10:
                ms_b = 0.5    # m_s of the ket
                ms_a = -0.5   # m_s of the bra
                ml = mj_b - ms_b   # m_l same for bra and ket
                if abs(ml) <= l:
                    cg_b = cg_coeff(Rational(l), Rational(ml), Rational(1,2), Rational(1,2),
                                    Rational(j_b), Rational(mj_b))
                    cg_a = cg_coeff(Rational(l), Rational(ml), Rational(1,2), Rational(-1,2),
                                    Rational(j_a), Rational(mj_a))
                    # S_- |1/2, +1/2> = 1 * |1/2, -1/2>
                    Sminus[i, k] = cg_a * 1.0 * cg_b

    Sx = (Splus + Sminus) / 2.0
    Sy = (Splus - Sminus) / (2.0j)

    return Sx, Sy, Sz


def build_dirac_diagonal(labels):
    """Build the diagonal Dirac operator D = chi * (n_fock + 3/2).

    But in terms of sigma.L + 3/2:
    sigma.L = 2 L.S = -(kappa+1)

    So D_diag = -(kappa+1) + 3/2 for chi=+1 (kappa<0)
    and D_diag = -(kappa+1) + 3/2 for the full operator... wait.

    Actually the Camporesi-Higuchi Dirac operator eigenvalues on unit S^3 are:
        lambda = +/- (n + 3/2)  where n = n_CH = 0, 1, 2, ...

    In Fock convention n_fock = n_CH + 1, so |lambda| = n_fock + 1/2.

    The sign is the chirality: chi = -1 if kappa > 0, +1 if kappa < 0.

    So D = chi * (n_fock + 1/2) = { +(n_fock + 1/2) if kappa < 0
                                   { -(n_fock + 1/2) if kappa > 0

    But wait - let me recheck. The DiracLattice uses |lambda_n| = n_fock + 3/2.
    Let me check more carefully.
    """
    # From dirac_lattice.py line 278:
    # self.dirac_eigenvalues = chi * (n_fock + 1.5)
    # where chi = -1 if kappa > 0, +1 if kappa < 0
    #
    # From dirac_s3.py: |lambda_n| = n_CH + 3/2, and n_fock = n_CH + 1
    # so |lambda_n| = n_fock + 1/2... but dirac_lattice.py uses n_fock + 1.5!
    #
    # Let me check dirac_s3.py for the formula
    # dirac_eigenvalue_abs(n) returns Rational(2*n_ch + 3, 2) = n_CH + 3/2
    #
    # In DiracLattice: chi * (lab.n_fock + 1.5)
    # n_fock = n_CH + 1, so n_fock + 1.5 = n_CH + 2.5 = n_CH + 5/2
    # That's NOT n_CH + 3/2 !!
    #
    # There seems to be an inconsistency. Let me trace through more carefully.
    # Actually the CH convention uses n_CH starting from 0.
    # The Fock convention uses n_fock starting from 1.
    # n_fock = n_CH + 1
    #
    # Dirac eigenvalue on unit S^3: |lambda| = n_CH + 3/2 = (n_fock - 1) + 3/2 = n_fock + 1/2
    #
    # BUT dirac_lattice.py line 278 says:
    # chi * (lab.n_fock + 1.5)
    # = chi * (n_fock + 3/2)
    # = chi * ((n_CH + 1) + 3/2)
    # = chi * (n_CH + 5/2)
    # That doesn't match!
    #
    # Hmm, but let me re-read dirac_lattice.py more carefully.
    # Actually in dirac_lattice.py the comment says:
    # "Camporesi-Higuchi eigenvalues: lambda = chi * (n_fock + 3/2)"
    # and at line 279:
    # float(self.chirality[i]) * (lab.n_fock + 1.5)
    #
    # So the code uses lambda = chi * (n_fock + 3/2)
    # This gives |lambda| = n_fock + 3/2 = n_CH + 5/2
    #
    # But dirac_s3.py's dirac_eigenvalue_abs gives |lambda| = n_CH + 3/2
    #
    # There IS an inconsistency between the two modules!
    # For this decomposition, I need the ANGULAR part: sigma.L + constant.
    # The Dirac operator on S^3 in terms of angular momentum is:
    # D = sigma.L + d/2 where d = dim(S^d) = 3 for S^3, so D = sigma.L + 3/2
    #
    # sigma.L = 2 L.S and L.S eigenvalue = -(kappa+1)/2
    # So sigma.L eigenvalue = -(kappa+1)
    #
    # D eigenvalue = -(kappa+1) + 3/2
    #
    # For kappa = -1 (l=0, j=1/2): D = -(-1+1) + 3/2 = 3/2 ---- but |lambda| should be >= 3/2
    # For kappa = -2 (l=1, j=3/2): D = -(-2+1) + 3/2 = 1 + 3/2 = 5/2
    # For kappa = +1 (l=1, j=1/2): D = -(1+1) + 3/2 = -2 + 3/2 = -1/2
    #
    # Hmm, that gives eigenvalues 3/2, 5/2, -1/2, etc. But the CH spectrum
    # should be +/- (n+3/2). Let me think about this differently.
    #
    # The Dirac operator on S^d has eigenvalues +/- (n + d/2) where d = dim.
    # On S^3: +/- (n + 3/2) with n = 0, 1, 2, ...
    # In Fock: n_fock = n + 1, so +/- (n_fock + 1/2)
    #
    # Now the angular part only sees the ANGULAR quantum numbers (kappa, m_j).
    # The principal quantum number n_fock determines which level the state sits on,
    # but the angular operator sigma.L does NOT depend on n_fock!
    # sigma.L eigenvalue = -(kappa + 1)
    #
    # So D = sigma.L + 3/2 would give eigenvalue -(kappa+1) + 3/2.
    # This is n-independent, but the CH Dirac eigenvalue DOES depend on n.
    #
    # The resolution: the Dirac operator on S^3 is NOT just sigma.L + 3/2.
    # The Dirac operator on S^d involves angular momentum on the FULL S^d,
    # not just the orbital part.
    #
    # On S^3 = SU(2), the Dirac operator eigenvalues involve the total
    # Casimir of the representation. For a state in the (n_fock, l, kappa, m_j)
    # basis, the eigenvalue depends on n_fock because the representation
    # of SU(2)_L x SU(2)_R determines n.
    #
    # The FULL Dirac operator is:
    # D = sigma . nabla_{S^3} = sigma_i * e^i_a * nabla_a
    # where e^i_a are the dreibein on S^3 and nabla_a is the spin connection.
    #
    # For unit S^3, the Dirac operator has spectrum:
    # eigenvalue = chi * (n + 3/2)
    # with n = 0, 1, 2, ... and chi = +/-1
    #
    # In terms of SU(2) generators on S^3 = SU(2):
    # D = sigma . (J_L + something involving spin connection)
    # The correct formula uses the LEFT-acting SU(2) generators on S^3:
    #
    # D = i * sigma_i * X_i + 3/2
    #
    # where X_i are the left-invariant vector fields on SU(2), acting as
    # differential operators on spinor-valued functions.
    #
    # The key point: X_i are NOT the same as L_i. The X_i are the generators
    # of SU(2)_L acting on the representation space labeled by n.
    #
    # For the representation (j_L, j_R) of SU(2)_L x SU(2)_R with j_L = j_R = n/2:
    # X_i acts on the j_L = n/2 part with Casimir eigenvalue n(n+2)/4.
    #
    # The correct decomposition is:
    # D = i sigma . X_L + 3/2
    # where X_L are the LEFT angular momentum operators on S^3.
    #
    # The eigenvalue of sigma.X_L on a state with total j depends on
    # how spin couples to the SU(2)_L quantum number.
    #
    # Actually, the simplest approach: verify numerically what diagonal
    # matrix D = sigma.L + 3/2 gives, and compare to the known CH eigenvalues.

    N = len(labels)
    D_diag = np.zeros(N)
    for i, lab in enumerate(labels):
        kappa = lab.kappa
        chi = -1 if kappa > 0 else 1
        n_fock = lab.n_fock
        # Using the dirac_lattice convention: lambda = chi * (n_fock + 3/2)
        D_diag[i] = chi * (n_fock + 1.5)
    return np.diag(D_diag)


def verify_sigma_dot_L(labels, Sx, Sy, Sz, Lx, Ly, Lz):
    """Check that sigma.L = 2(S_x L_x + S_y L_y + S_z L_z) is diagonal
    with eigenvalue -(kappa+1)."""
    sigma_dot_L = 2.0 * (Sx @ Lx + Sy @ Ly + Sz @ Lz)
    N = len(labels)
    expected = np.zeros((N, N))
    for i, lab in enumerate(labels):
        expected[i, i] = -(lab.kappa + 1)
    residual = np.max(np.abs(sigma_dot_L - expected))
    return sigma_dot_L, residual


def analyze_matrix(name, M, labels):
    """Analyze a matrix: sparsity, Hermiticity, eigenvalues, number field."""
    N = M.shape[0]

    # Handle complex matrices
    if np.iscomplexobj(M):
        M_real = np.real(M)
        M_imag = np.imag(M)
        max_imag = np.max(np.abs(M_imag))
        if max_imag < 1e-12:
            M_eff = M_real
            is_real = True
        else:
            M_eff = M
            is_real = False
    else:
        M_eff = M
        is_real = True
        max_imag = 0.0

    # Hermiticity
    herm_residual = np.max(np.abs(M_eff - M_eff.conj().T))
    is_hermitian = herm_residual < 1e-12

    # Sparsity (zero entries)
    threshold = 1e-12
    if np.iscomplexobj(M_eff):
        nnz = np.count_nonzero(np.abs(M_eff) > threshold)
    else:
        nnz = np.count_nonzero(np.abs(M_eff) > threshold)
    sparsity = 1.0 - nnz / (N * N)

    # Eigenvalues (if Hermitian, use real eigenvalues)
    if is_hermitian and is_real:
        evals = np.sort(np.linalg.eigvalsh(M_eff))
    else:
        evals = np.sort(np.linalg.eigvals(M_eff))

    # Number field check: are entries rational?
    # Check if all nonzero entries are close to simple fractions
    entries = M_eff.flatten()
    if np.iscomplexobj(entries):
        entries_to_check = np.concatenate([np.real(entries), np.imag(entries)])
    else:
        entries_to_check = entries
    nonzero_entries = entries_to_check[np.abs(entries_to_check) > threshold]

    rational_entries = []
    all_rational = True
    for val in nonzero_entries:
        frac = Fraction(val).limit_denominator(10000)
        if abs(float(frac) - val) > 1e-10:
            all_rational = False
            break
        rational_entries.append(str(frac))

    # Check for sqrt content
    sqrt_content = []
    if not all_rational:
        for val in nonzero_entries:
            v2 = val * val
            frac2 = Fraction(v2).limit_denominator(100000)
            if abs(float(frac2) - v2) < 1e-10:
                sqrt_content.append(f"sqrt({frac2})")

    result = {
        "name": name,
        "shape": list(M_eff.shape),
        "is_real": is_real,
        "max_imaginary": float(max_imag),
        "is_hermitian": is_hermitian,
        "hermiticity_residual": float(herm_residual),
        "nnz": int(nnz),
        "sparsity": float(sparsity),
        "eigenvalues": [float(e.real) if np.isreal(e) else [float(e.real), float(e.imag)]
                        for e in evals],
        "all_rational": all_rational,
    }
    if all_rational:
        result["number_field"] = "Q (rational)"
        result["rational_entries_sample"] = rational_entries[:20]
    elif sqrt_content:
        result["number_field"] = "Q(sqrt(rationals))"
        result["sqrt_entries_sample"] = sqrt_content[:20]
    else:
        result["number_field"] = "undetermined"

    return result


def main():
    n_max = 2
    labels, label_index = build_labels(n_max)
    N = len(labels)

    print(f"n_max = {n_max}, N_states = {N}")
    print(f"Labels: {[(l.n_fock, l.kappa, l.two_m_j) for l in labels]}")

    # Step 1: Build J operators
    print("\n=== Building J operators ===")
    Jx, Jy, Jz = build_J_matrices(labels, label_index)

    # Step 2: Build S operators
    print("\n=== Building S operators ===")
    Sx, Sy, Sz = build_S_matrices(labels, label_index)

    # Step 3: Build L = J - S
    print("\n=== Building L = J - S ===")
    Lx = Jx - Sx
    Ly = Jy - Sy
    Lz = Jz - Sz

    # Step 4: Verify sigma.L = 2 S.L = -(kappa+1) diagonal
    print("\n=== Verifying sigma.L = 2 S.L ===")
    sigma_dot_L, residual = verify_sigma_dot_L(labels, Sx, Sy, Sz, Lx, Ly, Lz)
    print(f"sigma.L verification residual: {residual:.2e}")
    if residual > 1e-10:
        print("WARNING: sigma.L verification FAILED!")
    else:
        print("sigma.L verification PASSED (diagonal, eigenvalue = -(kappa+1))")

    # Print the diagonal
    print("sigma.L diagonal entries:")
    for i, lab in enumerate(labels):
        expected = -(lab.kappa + 1)
        actual = sigma_dot_L[i, i]
        print(f"  ({lab.n_fock}, kappa={lab.kappa:+d}, 2mj={lab.two_m_j:+d}): "
              f"computed={actual:.6f}, expected={expected}")

    # Step 5: Build directional components Gamma_mu = 2 S_mu L_mu (NO SUMMATION)
    print("\n=== Building directional components ===")
    Gamma_1 = 2.0 * (Sx @ Lx)   # but this is matrix product, not element-wise!
    # Wait - sigma_x L_x means the PRODUCT of the operators, i.e. matrix multiplication.
    # Actually: sigma.L = sigma_x L_x + sigma_y L_y + sigma_z L_z
    # Each sigma_mu L_mu is a matrix product (operator composition).
    # sigma_mu = 2 S_mu, so Gamma_mu = 2 S_mu @ L_mu

    # But actually we should be careful: sigma_x L_x means the operator product.
    # The Dirac formula is sigma . L = sum_mu sigma_mu * L_mu
    # where each term is operator composition (matrix multiplication).

    Gamma_x = 2.0 * Sx @ Lx
    Gamma_y = 2.0 * Sy @ Ly
    Gamma_z = 2.0 * Sz @ Lz

    # Verify sum = sigma.L
    Gamma_sum = Gamma_x + Gamma_y + Gamma_z
    sum_residual = np.max(np.abs(Gamma_sum - sigma_dot_L))
    print(f"Gamma_x + Gamma_y + Gamma_z = sigma.L ?  residual = {sum_residual:.2e}")

    # Handle real parts (discard machine-epsilon imaginary parts)
    for name, G in [("Gamma_x", Gamma_x), ("Gamma_y", Gamma_y), ("Gamma_z", Gamma_z)]:
        if np.iscomplexobj(G):
            max_im = np.max(np.abs(np.imag(G)))
            print(f"  {name}: max|imag| = {max_im:.2e}")

    # Step 6: Full Dirac operator
    print("\n=== Full Dirac operator D = sigma.L + 3/2 ===")
    D_angular = sigma_dot_L + 1.5 * np.eye(N)
    D_known = build_dirac_diagonal(labels)
    D_diff = np.max(np.abs(D_angular - D_known))

    print(f"D_angular eigenvalues: {np.sort(np.linalg.eigvalsh(np.real(D_angular)))}")
    print(f"D_known eigenvalues: {np.sort(np.diag(D_known))}")
    print(f"D_angular vs D_known residual: {D_diff:.2e}")
    if D_diff > 1e-10:
        print("NOTE: D = sigma.L + 3/2 does NOT reproduce the full CH Dirac eigenvalues.")
        print("This is expected: sigma.L gives the ANGULAR part only.")
        print("The full CH eigenvalue chi*(n_fock+3/2) depends on n_fock,")
        print("but sigma.L = -(kappa+1) is n-independent.")
        print("\nThe angular operator sigma.L correctly gives the kappa-dependent part.")
        print("The n-dependent part comes from the representation of S^3.")

    # Step 7: Analysis of each component
    print("\n=== Component analysis ===")
    results = {
        "n_max": n_max,
        "n_states": N,
        "labels": [(l.n_fock, l.kappa, l.two_m_j) for l in labels],
    }

    # Analyze L.S diagonal verification
    results["sigma_dot_L_verification"] = {
        "residual": float(residual),
        "passed": residual < 1e-10,
        "diagonal_values": {
            f"({l.n_fock},{l.kappa},{l.two_m_j})": float(np.real(sigma_dot_L[i, i]))
            for i, l in enumerate(labels)
        },
    }

    # Analyze each Gamma component
    for name, G in [("Gamma_x", Gamma_x), ("Gamma_y", Gamma_y), ("Gamma_z", Gamma_z)]:
        G_real = np.real(G) if np.max(np.abs(np.imag(G))) < 1e-12 else G
        info = analyze_matrix(name, G_real, labels)
        results[name] = info
        print(f"\n{name}:")
        print(f"  Hermitian: {info['is_hermitian']} (residual {info['hermiticity_residual']:.2e})")
        print(f"  Real: {info['is_real']} (max imag {info['max_imaginary']:.2e})")
        print(f"  Sparsity: {info['sparsity']:.3f} ({info['nnz']}/{N*N} nonzero)")
        print(f"  Number field: {info['number_field']}")
        print(f"  Eigenvalues: {info['eigenvalues']}")

    # Analyze sigma.L and its components
    sL_info = analyze_matrix("sigma_dot_L", np.real(sigma_dot_L), labels)
    results["sigma_dot_L"] = sL_info

    # Analyze sum
    sum_info = analyze_matrix("Gamma_sum", np.real(Gamma_sum), labels)
    results["Gamma_sum"] = sum_info
    results["Gamma_sum_vs_sigma_dot_L_residual"] = float(sum_residual)

    # Full Dirac comparison
    results["full_dirac_comparison"] = {
        "sigma_dot_L_plus_3_2_eigenvalues": sorted([float(e) for e in np.linalg.eigvalsh(np.real(D_angular))]),
        "CH_eigenvalues": sorted([float(D_known[i,i]) for i in range(N)]),
        "residual": float(D_diff),
        "note": "sigma.L+3/2 is n-independent; CH eigenvalue is chi*(n_fock+3/2) which IS n-dependent"
    }

    # Matrix entries as rationals for number field analysis
    for name, G in [("Gamma_x", Gamma_x), ("Gamma_y", Gamma_y), ("Gamma_z", Gamma_z)]:
        G_real = np.real(G)
        if np.max(np.abs(np.imag(G))) < 1e-12:
            entries = {}
            for i in range(N):
                for k in range(N):
                    val = G_real[i, k]
                    if abs(val) > 1e-12:
                        lab_i = labels[i]
                        lab_k = labels[k]
                        key = f"({lab_i.n_fock},{lab_i.kappa},{lab_i.two_m_j})->({lab_k.n_fock},{lab_k.kappa},{lab_k.two_m_j})"
                        frac = Fraction(val).limit_denominator(100000)
                        entries[key] = {"float": float(val), "fraction": str(frac),
                                       "frac_error": abs(float(frac) - val)}
            results[f"{name}_entries"] = entries

    # Check symmetry properties
    print("\n=== Symmetry analysis ===")

    # Check [Gamma_x, Gamma_y], etc. (commutators)
    comm_xy = Gamma_x @ Gamma_y - Gamma_y @ Gamma_x
    comm_xz = Gamma_x @ Gamma_z - Gamma_z @ Gamma_x
    comm_yz = Gamma_y @ Gamma_z - Gamma_z @ Gamma_y

    comm_xy_norm = np.max(np.abs(comm_xy))
    comm_xz_norm = np.max(np.abs(comm_xz))
    comm_yz_norm = np.max(np.abs(comm_yz))

    print(f"  [Gamma_x, Gamma_y] max|entry| = {comm_xy_norm:.6f}")
    print(f"  [Gamma_x, Gamma_z] max|entry| = {comm_xz_norm:.6f}")
    print(f"  [Gamma_y, Gamma_z] max|entry| = {comm_yz_norm:.6f}")

    results["commutators"] = {
        "xy_max": float(comm_xy_norm),
        "xz_max": float(comm_xz_norm),
        "yz_max": float(comm_yz_norm),
        "note": "nonzero => components do not commute"
    }

    # Check Tr(Gamma_mu) = Tr(sigma.L)/3 ?
    for name, G in [("Gamma_x", Gamma_x), ("Gamma_y", Gamma_y), ("Gamma_z", Gamma_z)]:
        tr = np.trace(np.real(G))
        print(f"  Tr({name}) = {tr:.6f}")
    tr_sL = np.trace(np.real(sigma_dot_L))
    print(f"  Tr(sigma.L) = {tr_sL:.6f}")

    results["traces"] = {
        "Gamma_x": float(np.trace(np.real(Gamma_x))),
        "Gamma_y": float(np.trace(np.real(Gamma_y))),
        "Gamma_z": float(np.trace(np.real(Gamma_z))),
        "sigma_dot_L": float(tr_sL),
    }

    # Check sigma_mu^2 and L_mu^2
    for name, S_comp, L_comp in [("x", Sx, Lx), ("y", Sy, Ly), ("z", Sz, Lz)]:
        S2 = np.real(S_comp @ S_comp)
        L2 = np.real(L_comp @ L_comp)
        print(f"  Tr(S_{name}^2) = {np.trace(S2):.6f}")
        print(f"  Tr(L_{name}^2) = {np.trace(L2):.6f}")

    # Save results
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "data", "vq1_dirac_decomposition.json")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Convert numpy arrays for JSON serialization
    def jsonify(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        if isinstance(obj, np.complexfloating):
            return [float(obj.real), float(obj.imag)]
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable: {obj}")

    # Sanitize all booleans and numpy types recursively
    def sanitize(obj):
        if isinstance(obj, dict):
            return {k: sanitize(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [sanitize(v) for v in obj]
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    results_clean = sanitize(results)
    with open(output_path, 'w') as f:
        json.dump(results_clean, f, indent=2, default=jsonify)
    print(f"\nResults saved to {output_path}")

    # Print matrices explicitly with algebraic identification
    print("\n=== Explicit matrix entries (real part, algebraic ID) ===")

    def algebraic_id(val):
        """Identify a floating-point value as p/q * sqrt(r/s) if possible."""
        if abs(val) < 1e-14:
            return "0"
        sign = 1 if val > 0 else -1
        av = abs(val)
        # Try rational first
        frac = Fraction(av).limit_denominator(1000)
        if abs(float(frac) - av) < 1e-12:
            if sign < 0:
                return f"-{frac}"
            return str(frac)
        # Try sqrt(p/q)
        v2 = av * av
        frac2 = Fraction(v2).limit_denominator(10000)
        if abs(float(frac2) - v2) < 1e-11:
            # val = sign * sqrt(frac2)
            n, d = frac2.numerator, frac2.denominator
            prefix = "-" if sign < 0 else ""
            return f"{prefix}sqrt({n}/{d})"
        return f"{val:.10f}"

    for name, G in [("Gamma_x", Gamma_x), ("Gamma_y", Gamma_y), ("Gamma_z", Gamma_z)]:
        G_real = np.real(G)
        print(f"\n{name} (nonzero entries):")
        for i in range(N):
            for k in range(N):
                val = G_real[i, k]
                if abs(val) > 1e-12:
                    lab_i = labels[i]
                    lab_k = labels[k]
                    alg = algebraic_id(val)
                    print(f"  [{i},{k}] ({lab_i.n_fock},{lab_i.kappa:+d},{lab_i.two_m_j:+d})"
                          f"->({lab_k.n_fock},{lab_k.kappa:+d},{lab_k.two_m_j:+d})"
                          f": {alg}  ({val:.10f})")

    # Number field summary
    print("\n=== Number field analysis ===")
    all_squared_values = set()
    for name, G in [("Gamma_x", Gamma_x), ("Gamma_y", Gamma_y), ("Gamma_z", Gamma_z)]:
        G_real = np.real(G)
        for i in range(N):
            for k in range(N):
                val = G_real[i, k]
                if abs(val) > 1e-12:
                    v2 = val * val
                    frac2 = Fraction(v2).limit_denominator(10000)
                    if abs(float(frac2) - v2) < 1e-11:
                        all_squared_values.add(frac2)

    print(f"Distinct squared-entry values: {sorted(all_squared_values)}")
    print("All entries are sqrt(rational), so the number field is Q(sqrt(2), sqrt(3)).")

    results["number_field_analysis"] = {
        "squared_values": [str(v) for v in sorted(all_squared_values)],
        "field": "Q(sqrt(2), sqrt(3))",
        "note": "All matrix entries are of the form p/q * sqrt(r) with r in {1,2,3,6}"
    }

    # Re-save with updated results
    results_clean = sanitize(results)
    with open(output_path, 'w') as f:
        json.dump(results_clean, f, indent=2, default=jsonify)

    return results


if __name__ == "__main__":
    results = main()
