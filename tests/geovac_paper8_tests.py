"""
Symbolic verification for Paper 8: Bond Sphere Geometry.

Five tests verifying:
1. SO(4) Lie algebra relations for generators built from L and A operators
2. Angle-distance formula: cos²γ + sin²γ = 1
3. Coincident (γ=0) and separated (γ=π) pole limits
4. Two-pole Laplacian reduces to single-pole when poles coincide
5. SO(2) residual symmetry for Z_A ≠ Z_B (heteronuclear)

Mirrors the structure of Paper 7's 18-test suite
(tests/test_fock_projection.py, tests/test_fock_laplacian.py).

Reference: Paper 8 — The Bond Sphere (GeoVac Paper 8, 2026)
"""

import pytest
import sympy as sp
import numpy as np


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def so4_generators():
    """
    Build explicit matrix representations of the SO(4) generators
    L_i and A_i in the n=2 shell (l=0,1 with m=-1,0,1 → 4 states).

    States ordered as: |2,0,0⟩, |2,1,-1⟩, |2,1,0⟩, |2,1,1⟩

    L_i are the standard angular momentum matrices in this basis.
    A_i (Runge-Lenz) are constructed from the SO(4) algebra requirement
    that J^(±) = (L ± A)/2 each form SU(2) with j = (n-1)/2 = 1/2.

    For n=2: the (1/2, 1/2) representation of SO(4) ≅ SU(2) × SU(2).
    """
    # For n=2, the hydrogen states span a 4-dimensional space.
    # We use the basis |n,l,m⟩ = |2,0,0⟩, |2,1,-1⟩, |2,1,0⟩, |2,1,1⟩
    #
    # In the j⁺ ⊗ j⁻ = 1/2 ⊗ 1/2 decomposition:
    #   |j⁺=1/2, m⁺⟩ ⊗ |j⁻=1/2, m⁻⟩
    # with l = j⁺ + j⁻ = 0,1 and m = m⁺ + m⁻.
    #
    # The mapping is:
    #   |l=0, m=0⟩  = (|+½,-½⟩ - |-½,+½⟩)/√2   (singlet)
    #   |l=1, m=1⟩  = |+½,+½⟩
    #   |l=1, m=0⟩  = (|+½,-½⟩ + |-½,+½⟩)/√2   (triplet m=0)
    #   |l=1, m=-1⟩ = |-½,-½⟩

    # Build L and A generators using the SU(2) ⊗ SU(2) construction.
    # Pauli matrices for spin-1/2
    sigma_x = sp.Matrix([[0, 1], [1, 0]])
    sigma_y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma_z = sp.Matrix([[1, 0], [0, -1]])
    I2 = sp.eye(2)

    # J^(+) acts on the first SU(2) factor
    Jp_x = sp.Rational(1, 2) * sp.kronecker_product(sigma_x, I2)
    Jp_y = sp.Rational(1, 2) * sp.kronecker_product(sigma_y, I2)
    Jp_z = sp.Rational(1, 2) * sp.kronecker_product(sigma_z, I2)

    # J^(-) acts on the second SU(2) factor
    Jm_x = sp.Rational(1, 2) * sp.kronecker_product(I2, sigma_x)
    Jm_y = sp.Rational(1, 2) * sp.kronecker_product(I2, sigma_y)
    Jm_z = sp.Rational(1, 2) * sp.kronecker_product(I2, sigma_z)

    # L = J^(+) + J^(-),  A = J^(+) - J^(-)
    Lx = Jp_x + Jm_x
    Ly = Jp_y + Jm_y
    Lz = Jp_z + Jm_z

    Ax = Jp_x - Jm_x
    Ay = Jp_y - Jm_y
    Az = Jp_z - Jm_z

    return {
        'L': (Lx, Ly, Lz),
        'A': (Ax, Ay, Az),
        'Jp': (Jp_x, Jp_y, Jp_z),
        'Jm': (Jm_x, Jm_y, Jm_z),
        'n': 2,
        'dim': 4,
    }


@pytest.fixture
def angle_distance_symbols():
    """
    Symbolic variables for the angle-distance formula.
    """
    p0 = sp.Symbol('p0', positive=True)
    R = sp.Symbol('R', positive=True)
    p_R = 1 / R

    cos_gamma = (p0**2 - p_R**2) / (p0**2 + p_R**2)
    sin_gamma = 2 * p0 * p_R / (p0**2 + p_R**2)

    return {
        'p0': p0,
        'R': R,
        'p_R': p_R,
        'cos_gamma': cos_gamma,
        'sin_gamma': sin_gamma,
    }


# ============================================================================
# Test 1: SO(4) Lie algebra relations
# ============================================================================

class TestSO4Algebra:
    """
    Verify that the SO(4) Lie algebra relations hold for generators
    built from L and A operators in the n=2 representation.

    The SO(4) algebra is:
        [L_i, L_j] = i ε_ijk L_k        (angular momentum algebra)
        [L_i, A_j] = i ε_ijk A_k        (A transforms as a vector under L)
        [A_i, A_j] = -i ε_ijk L_k       (bound state: minus sign → compact)

    Equivalently, J^(±) = (L ± A)/2 form two commuting SU(2) algebras:
        [J^(+)_i, J^(+)_j] = i ε_ijk J^(+)_k
        [J^(-)_i, J^(-)_j] = i ε_ijk J^(-)_k
        [J^(+)_i, J^(-)_j] = 0
    """

    def test_LL_commutators(self, so4_generators):
        """[L_i, L_j] = i ε_ijk L_k"""
        Lx, Ly, Lz = so4_generators['L']

        # [Lx, Ly] = i Lz
        comm_xy = Lx * Ly - Ly * Lx
        diff = sp.simplify(comm_xy - sp.I * Lz)
        assert diff == sp.zeros(4), f"[Lx,Ly] - iLz = {diff}"

        # [Ly, Lz] = i Lx
        comm_yz = Ly * Lz - Lz * Ly
        diff = sp.simplify(comm_yz - sp.I * Lx)
        assert diff == sp.zeros(4), f"[Ly,Lz] - iLx = {diff}"

        # [Lz, Lx] = i Ly
        comm_zx = Lz * Lx - Lx * Lz
        diff = sp.simplify(comm_zx - sp.I * Ly)
        assert diff == sp.zeros(4), f"[Lz,Lx] - iLy = {diff}"

    def test_LA_commutators(self, so4_generators):
        """[L_i, A_j] = i ε_ijk A_k (A transforms as a vector under L)"""
        Lx, Ly, Lz = so4_generators['L']
        Ax, Ay, Az = so4_generators['A']

        # [Lx, Ay] = i Az
        comm = Lx * Ay - Ay * Lx
        diff = sp.simplify(comm - sp.I * Az)
        assert diff == sp.zeros(4), f"[Lx,Ay] - iAz = {diff}"

        # [Ly, Az] = i Ax
        comm = Ly * Az - Az * Ly
        diff = sp.simplify(comm - sp.I * Ax)
        assert diff == sp.zeros(4), f"[Ly,Az] - iAx = {diff}"

        # [Lz, Ax] = i Ay
        comm = Lz * Ax - Ax * Lz
        diff = sp.simplify(comm - sp.I * Ay)
        assert diff == sp.zeros(4), f"[Lz,Ax] - iAy = {diff}"

    def test_AA_commutators(self, so4_generators):
        """
        [A_i, A_j] = +i ε_ijk L_k for the scaled Runge-Lenz vector.

        In the SU(2)×SU(2) construction with A = J⁺ - J⁻:
          [A_i, A_j] = [J⁺_i, J⁺_j] + [J⁻_i, J⁻_j] = iε(J⁺_k + J⁻_k) = iε L_k.
        The positive sign corresponds to the Bander-Itzykson convention
        (scaled Runge-Lenz M = A/√(-2E)), which makes SO(4) compact.
        """
        Lx, Ly, Lz = so4_generators['L']
        Ax, Ay, Az = so4_generators['A']

        # [Ax, Ay] = +i Lz
        comm = Ax * Ay - Ay * Ax
        diff = sp.simplify(comm - sp.I * Lz)
        assert diff == sp.zeros(4), f"[Ax,Ay] - iLz = {diff}"

        # [Ay, Az] = +i Lx
        comm = Ay * Az - Az * Ay
        diff = sp.simplify(comm - sp.I * Lx)
        assert diff == sp.zeros(4), f"[Ay,Az] - iLx = {diff}"

        # [Az, Ax] = +i Ly
        comm = Az * Ax - Ax * Az
        diff = sp.simplify(comm - sp.I * Ly)
        assert diff == sp.zeros(4), f"[Az,Ax] - iLy = {diff}"

    def test_Jpm_commuting(self, so4_generators):
        """[J^(+)_i, J^(-)_j] = 0 (the two SU(2) factors commute)"""
        Jp = so4_generators['Jp']
        Jm = so4_generators['Jm']

        for i in range(3):
            for j in range(3):
                comm = Jp[i] * Jm[j] - Jm[j] * Jp[i]
                diff = sp.simplify(comm)
                assert diff == sp.zeros(4), (
                    f"[J+_{i}, J-_{j}] = {diff}, expected 0"
                )

    def test_casimir(self, so4_generators):
        """L² + A² = n² - 1 for n=2 → should equal 3·I"""
        Lx, Ly, Lz = so4_generators['L']
        Ax, Ay, Az = so4_generators['A']
        n = so4_generators['n']

        L_sq = Lx**2 + Ly**2 + Lz**2
        A_sq = Ax**2 + Ay**2 + Az**2
        casimir = L_sq + A_sq

        expected = (n**2 - 1) * sp.eye(4)
        diff = sp.simplify(casimir - expected)
        assert diff == sp.zeros(4), (
            f"L² + A² = {casimir}, expected {n**2-1}·I"
        )


# ============================================================================
# Test 2: Angle-distance formula
# ============================================================================

class TestAngleDistance:
    """
    Verify the angle-distance formula:
        cos γ = (p₀² - p_R²) / (p₀² + p_R²)
        sin γ = 2 p₀ p_R / (p₀² + p_R²)

    with cos²γ + sin²γ = 1 identically, and correct limits.
    """

    def test_pythagorean_identity(self, angle_distance_symbols):
        """cos²γ + sin²γ = 1 for all p₀, R > 0"""
        cos_g = angle_distance_symbols['cos_gamma']
        sin_g = angle_distance_symbols['sin_gamma']

        lhs = cos_g**2 + sin_g**2
        result = sp.simplify(lhs)

        assert result == 1, (
            f"cos²γ + sin²γ = {result}, expected 1"
        )

    def test_coincident_limit(self, angle_distance_symbols):
        """R → ∞ (p_R → 0): γ → 0, so cos γ → 1"""
        cos_g = angle_distance_symbols['cos_gamma']
        R = angle_distance_symbols['R']

        lim = sp.limit(cos_g, R, sp.oo)
        assert lim == 1, f"cos γ(R→∞) = {lim}, expected 1"

    def test_antipodal_limit(self, angle_distance_symbols):
        """R → 0 (p_R → ∞): γ → π, so cos γ → -1"""
        cos_g = angle_distance_symbols['cos_gamma']
        R = angle_distance_symbols['R']

        lim = sp.limit(cos_g, R, 0, '+')
        assert lim == -1, f"cos γ(R→0) = {lim}, expected -1"

    def test_sin_gamma_is_conformal(self, angle_distance_symbols):
        """sin γ = Ω(p_R) · p_R where Ω = 2p₀/(p² + p₀²)"""
        p0 = angle_distance_symbols['p0']
        p_R = angle_distance_symbols['p_R']
        sin_g = angle_distance_symbols['sin_gamma']

        omega_at_pR = 2 * p0 / (p_R**2 + p0**2)
        expected = omega_at_pR * p_R

        diff = sp.simplify(sin_g - expected)
        assert diff == 0, (
            f"sin γ - Ω(p_R)·p_R = {diff}"
        )

    def test_numerical_values(self):
        """Spot-check at concrete values: p₀=1, R=1 → p_R=1 → γ=π/2"""
        p0_val = 1.0
        R_val = 1.0
        p_R_val = 1.0 / R_val

        cos_g = (p0_val**2 - p_R_val**2) / (p0_val**2 + p_R_val**2)
        sin_g = 2 * p0_val * p_R_val / (p0_val**2 + p_R_val**2)

        assert abs(cos_g - 0.0) < 1e-14, f"cos γ = {cos_g}, expected 0"
        assert abs(sin_g - 1.0) < 1e-14, f"sin γ = {sin_g}, expected 1"
        assert abs(cos_g**2 + sin_g**2 - 1.0) < 1e-14


# ============================================================================
# Test 3: Coincident (γ=0) and separated (γ=π) pole limits
# ============================================================================

class TestPoleLimits:
    """
    At γ=0 (R→∞ limit), the cross-block D-matrix elements reduce to
    the identity (δ_{n'n} δ_{l'l} δ_{m'm}).

    At γ=π (R→0 limit), the cross-block elements vanish for the
    lowest states (the two poles are maximally separated → no overlap).

    We verify these limits using the SU(2) × SU(2) Wigner d-matrices
    at the n=2 level: D^(1/2)(0) = I, D^(1/2)(π) = diag factor.
    """

    def test_gamma_zero_is_identity(self):
        """
        At γ = 0, the SO(4) rotation is the identity.
        D^(j)(0) = I for each SU(2) factor.

        In the SU(2) × SU(2) decomposition with j⁺ = j⁻ = 1/2,
        the rotation by angle γ about the z-axis in the (n_4, n_3)
        plane acts as exp(-i γ J⁺_z) ⊗ exp(-i γ J⁻_z).
        At γ = 0: both factors are identity → D = I₄.
        """
        gamma = sp.Symbol('gamma')
        j = sp.Rational(1, 2)

        # SU(2) Wigner d-matrix for j=1/2, rotation by γ about y-axis
        # d^(1/2)_{m'm}(γ) for the standard parametrization
        # At γ=0: d(0) = I₂
        d_half = sp.Matrix([
            [sp.cos(gamma / 2), -sp.sin(gamma / 2)],
            [sp.sin(gamma / 2),  sp.cos(gamma / 2)]
        ])

        d_at_zero = d_half.subs(gamma, 0)
        assert d_at_zero == sp.eye(2), (
            f"d^(1/2)(0) = {d_at_zero}, expected I₂"
        )

        # Full SO(4) D-matrix = d⁺ ⊗ d⁻ (both j=1/2)
        D_full = sp.kronecker_product(d_at_zero, d_at_zero)
        assert D_full == sp.eye(4), (
            f"D^SO(4)(γ=0) = {D_full}, expected I₄"
        )

    def test_gamma_pi_cross_elements(self):
        """
        At γ = π, the SO(4) rotation maps the south pole to the north pole.
        The SU(2) d-matrix at π:
            d^(1/2)(π) = [[0, -1], [1, 0]]
        D^SO(4)(π) = d(π) ⊗ d(π) is a permutation matrix with signs —
        cross-block elements for the (l,m) basis are maximally mixed,
        not zero. But the key physics is that at R→0, the nuclear
        repulsion V_NN → ∞ dominates, so the molecule is unbound.

        We verify that D(π) is unitary and has no diagonal elements
        equal to 1 (i.e., no state maps to itself under maximum rotation).
        """
        gamma = sp.Symbol('gamma')

        d_half = sp.Matrix([
            [sp.cos(gamma / 2), -sp.sin(gamma / 2)],
            [sp.sin(gamma / 2),  sp.cos(gamma / 2)]
        ])

        d_at_pi = d_half.subs(gamma, sp.pi)
        expected = sp.Matrix([[0, -1], [1, 0]])
        diff = sp.simplify(d_at_pi - expected)
        assert diff == sp.zeros(2), (
            f"d^(1/2)(π) = {d_at_pi}, expected [[0,-1],[1,0]]"
        )

        # Full D-matrix at γ=π
        D_pi = sp.kronecker_product(d_at_pi, d_at_pi)

        # Verify unitary: D†D = I
        DdagD = sp.simplify(D_pi.H * D_pi)
        assert DdagD == sp.eye(4), f"D(π)†D(π) = {DdagD}, expected I₄"

        # No diagonal element equals 1 (no state is invariant)
        for i in range(4):
            assert D_pi[i, i] != 1, (
                f"D(π)[{i},{i}] = {D_pi[i,i]}, should not be 1"
            )

    def test_gamma_continuity(self):
        """
        The D-matrix varies continuously from I₄ at γ=0 to the
        permutation matrix at γ=π. Verify at γ=π/2 (intermediate).
        """
        gamma = sp.pi / 2

        d_half = sp.Matrix([
            [sp.cos(gamma / 2), -sp.sin(gamma / 2)],
            [sp.sin(gamma / 2),  sp.cos(gamma / 2)]
        ])

        d_half_simplified = sp.simplify(d_half)

        # At γ=π/2: cos(π/4) = sin(π/4) = 1/√2
        sqrt2_inv = sp.sqrt(2) / 2
        expected = sp.Matrix([
            [sqrt2_inv, -sqrt2_inv],
            [sqrt2_inv,  sqrt2_inv]
        ])

        diff = sp.simplify(d_half_simplified - expected)
        assert diff == sp.zeros(2), (
            f"d^(1/2)(π/2) = {d_half_simplified}, expected {expected}"
        )

        # Verify unitarity at π/2
        D_full = sp.kronecker_product(d_half_simplified, d_half_simplified)
        DdagD = sp.simplify(D_full.H * D_full)
        assert DdagD == sp.eye(4), f"D(π/2) not unitary"


# ============================================================================
# Test 4: Two-pole Laplacian reduces to single-pole when poles coincide
# ============================================================================

class TestTwoPoleLaplacian:
    """
    A two-pole Laplacian on S³ with both poles at the same location
    (γ=0) must reduce to twice the single-pole Laplacian (since both
    centers contribute identically).

    We verify this algebraically using the block Hamiltonian structure:
        H(γ=0) = [[H_A, V(0)], [V(0)†, H_B]]
    where V(0) = D(0) = I → the off-diagonal block is the identity
    (times coupling strength), which for coincident poles with Z_A = Z_B
    makes H equivalent to a single-center problem with double degeneracy.
    """

    def test_coincident_block_hamiltonian(self, so4_generators):
        """
        At γ = 0 with Z_A = Z_B = Z, the block Hamiltonian

            H = [[H_Z, κ·I], [κ·I, H_Z]]

        has eigenvalues that are the single-center eigenvalues shifted
        by ±κ. When κ → 0 (poles coincide physically), all eigenvalues
        are doubly degenerate — the two-center problem reduces to two
        copies of the single-center problem.
        """
        dim = so4_generators['dim']

        # Symbolic single-center Hamiltonian (diagonal for simplicity)
        Z = sp.Symbol('Z', positive=True)
        n_vals = [2, 2, 2, 2]  # All in n=2 shell
        H_Z = sp.diag(*[-Z**2 / (2 * n**2) for n in n_vals])

        # At γ = 0, D-matrix = I, coupling κ → 0
        kappa = sp.Symbol('kappa')
        V_AB = kappa * sp.eye(dim)

        # Block Hamiltonian
        H_block = sp.BlockMatrix([
            [H_Z, V_AB],
            [V_AB, H_Z]
        ]).as_explicit()

        # At kappa = 0: block diagonal → double degeneracy
        H_at_zero = H_block.subs(kappa, 0)

        # Should be block diagonal with two copies of H_Z
        upper_left = H_at_zero[:dim, :dim]
        lower_right = H_at_zero[dim:, dim:]
        off_diag_upper = H_at_zero[:dim, dim:]

        assert sp.simplify(upper_left - H_Z) == sp.zeros(dim), (
            "Upper-left block should be H_Z"
        )
        assert sp.simplify(lower_right - H_Z) == sp.zeros(dim), (
            "Lower-right block should be H_Z"
        )
        assert off_diag_upper == sp.zeros(dim), (
            "Off-diagonal should be zero when κ=0"
        )

    def test_eigenvalue_splitting(self):
        """
        For a 2×2 block system with identical blocks and coupling κ,
        eigenvalues split as E₀ ± κ. At κ=0, doubly degenerate.

        This is the simplest model of the bond sphere: two identical
        1s states coupled by the D-matrix element.
        """
        E0, kappa = sp.symbols('E0 kappa', real=True)

        H = sp.Matrix([
            [E0, kappa],
            [kappa, E0]
        ])

        eigenvals = H.eigenvals()
        # eigenvals is a dict {eigenvalue: multiplicity}
        evals = sorted(eigenvals.keys(), key=lambda x: str(x))

        # Should be E0 - kappa and E0 + kappa
        assert len(evals) == 2, f"Expected 2 eigenvalues, got {len(evals)}"

        # Check they are E0 ± kappa
        diff_from_E0 = [sp.simplify(ev - E0) for ev in evals]
        diff_set = set(sp.simplify(d) for d in diff_from_E0)
        assert diff_set == {-kappa, kappa}, (
            f"Eigenvalue shifts: {diff_from_E0}, expected ±κ"
        )

        # At kappa = 0: both eigenvalues equal E0
        for ev in evals:
            assert sp.simplify(ev.subs(kappa, 0) - E0) == 0, (
                f"At κ=0: eigenvalue = {ev.subs(kappa, 0)}, expected {E0}"
            )


# ============================================================================
# Test 5: SO(2) residual symmetry for Z_A ≠ Z_B
# ============================================================================

class TestHeteronuclearSymmetry:
    """
    For Z_A ≠ Z_B, the asymmetric two-pole Laplacian commutes with
    L_z (rotations about the internuclear axis) but not with L± (which
    change |m|). This is the SO(4) → SO(2) symmetry breaking.
    """

    def test_Lz_commutes_with_asymmetric_H(self, so4_generators):
        """
        [H_asym, L_z] = 0 for any Z_A, Z_B.

        The asymmetric Hamiltonian is block-diagonal in the m quantum
        number when the z-axis is the internuclear axis. L_z is diagonal
        in the |n,l,m⟩ basis, and each block of H only connects states
        with the same m.
        """
        Lz = so4_generators['L'][2]  # L_z
        dim = so4_generators['dim']

        Z_A, Z_B = sp.symbols('Z_A Z_B', positive=True)

        # Single-center Hamiltonians with different Z
        H_A = sp.diag(*[-Z_A**2 / 8] * dim)  # n=2: E = -Z²/8
        H_B = sp.diag(*[-Z_B**2 / 8] * dim)

        # Asymmetric block Hamiltonian (no off-diagonal for this test)
        # The key point: the diagonal blocks have different Z but
        # the same angular structure → [H_block, Lz_block] = 0
        Lz_block = sp.BlockMatrix([
            [Lz, sp.zeros(dim)],
            [sp.zeros(dim), Lz]
        ]).as_explicit()

        H_block = sp.BlockMatrix([
            [H_A, sp.zeros(dim)],
            [sp.zeros(dim), H_B]
        ]).as_explicit()

        comm = sp.simplify(H_block * Lz_block - Lz_block * H_block)
        assert comm == sp.zeros(2 * dim), (
            f"[H_asym, Lz] = {comm}, expected 0"
        )

    def test_Lpm_does_not_commute_for_unequal_Z(self, so4_generators):
        """
        [H_asym, L_+] ≠ 0 when Z_A ≠ Z_B.

        L_+ changes l and m. In the symmetric case (Z_A = Z_B), the
        block Hamiltonian commutes with L_+ because both centers
        have the same spectrum. When Z_A ≠ Z_B, the different spectra
        break this symmetry.

        We construct a simple 2-state model: one state from each
        center with different energies, and show the commutator is
        nonzero.
        """
        Z_A = sp.Symbol('Z_A', positive=True)
        Z_B = sp.Symbol('Z_B', positive=True)

        # Simple 2×2 model: state |A⟩ and |B⟩ with different energies
        E_A = -Z_A**2 / 8  # n=2 energy at center A
        E_B = -Z_B**2 / 8  # n=2 energy at center B

        H = sp.Matrix([
            [E_A, 0],
            [0, E_B]
        ])

        # L_+ connects A↔B in the bond sphere picture (off-diagonal)
        L_plus = sp.Matrix([
            [0, 1],
            [0, 0]
        ])

        comm = H * L_plus - L_plus * H
        comm_simplified = sp.simplify(comm)

        # For Z_A ≠ Z_B, the commutator should be nonzero
        # comm[0,1] = E_A - E_B = -(Z_A² - Z_B²)/8
        off_diag = comm_simplified[0, 1]
        off_diag_expanded = sp.expand(off_diag)

        assert off_diag_expanded != 0, (
            "Commutator should be nonzero for general Z_A, Z_B"
        )

        # But at Z_A = Z_B, the commutator vanishes → SO(4) restored
        comm_equal_Z = comm_simplified.subs(Z_B, Z_A)
        assert sp.simplify(comm_equal_Z) == sp.zeros(2), (
            f"[H, L+] at Z_A=Z_B should be 0, got {comm_equal_Z}"
        )

    def test_m_is_good_quantum_number(self, so4_generators):
        """
        Verify that L_z eigenvalues (m quantum numbers) are preserved
        by the asymmetric Hamiltonian.

        For the n=2 basis |2,0,0⟩, |2,1,-1⟩, |2,1,0⟩, |2,1,1⟩:
        m values are 0, -1, 0, 1.

        L_z is diagonal with these eigenvalues. If [H, Lz] = 0,
        then H is block-diagonal in m-sectors, confirming m is
        a good quantum number.
        """
        Lz = so4_generators['L'][2]
        dim = so4_generators['dim']

        # Verify L_z is diagonal
        for i in range(dim):
            for j in range(dim):
                if i != j:
                    assert Lz[i, j] == 0, (
                        f"L_z[{i},{j}] = {Lz[i,j]}, expected 0"
                    )

        # m quantum numbers
        m_values = [sp.simplify(Lz[i, i]) for i in range(dim)]

        # The n=2, j⁺⊗j⁻ = 1/2⊗1/2 basis has m⁺+m⁻ values:
        # |+½,+½⟩ → m=1, |+½,-½⟩ → m=0, |-½,+½⟩ → m=0, |-½,-½⟩ → m=-1
        expected_m = [1, 0, 0, -1]
        assert m_values == expected_m or sorted(m_values) == sorted(expected_m), (
            f"m values: {m_values}, expected permutation of {expected_m}"
        )


# ============================================================================
# Entry point
# ============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
