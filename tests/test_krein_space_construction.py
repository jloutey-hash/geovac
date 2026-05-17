"""Tests for `geovac.krein_space_construction` (Sprint L2-B).

These tests verify:

  (a) Cl(3, 1) gamma matrices in the West-coast chiral basis satisfy
      the Clifford algebra {g^mu, g^nu} = 2 eta^{mu nu} I exactly.
  (b) gamma^5 is the chirality grading and anticommutes with all g^mu.
  (c) The Krein space K = H_GV (x) L^2(R_t)_cutoff at n_max in {1, 2, 3}
      and N_t in {1, 11, 21} satisfies the Krein axioms:
        J^2 = +I
        J* J = I   (J unitary on the Hilbert level)
        J = J*    (J Hermitian)
        Krein inner product is conjugate-symmetric on random pairs
      all bit-exactly or to machine precision.
  (d) K splits as K = K^+ (+) K^- with each block of dimension dim / 2.
  (e) THE LOAD-BEARING FALSIFIER: the Riemannian-limit reduction at
      N_t = 1 recovers H_GV^{n_max} (Paper 32 def:H_GV) bit-identically:
      dim matches, basis labels match, and the Krein-space J reduces to
      the spatial fundamental symmetry J_spatial verbatim.  Also
      verified that the I_{N_t} temporal slot of J makes the time-slice
      diagonal block coincide with J_spatial for arbitrary N_t.

If the Riemannian-limit test fails, STOP and escalate.  See module
docstring and `debug/sprint_l2a_scoping_memo.md` Section 5.7.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.full_dirac_operator_system import (
    full_dirac_basis,
    full_dirac_dim,
)
from geovac.krein_space_construction import (
    GammaMatrices,
    KreinSpace,
    gamma_chiral,
    identity_temporal,
    pauli_matrices,
    spatial_fundamental_symmetry,
    temporal_grid,
    verify_gamma_anticommutation,
)


# ===========================================================================
# 1. Pauli matrices and Cl(3, 1) gamma-matrix algebra
# ===========================================================================


def test_pauli_matrices_squared_to_identity():
    """sigma^i squared equals I for i = 1, 2, 3."""
    s1, s2, s3 = pauli_matrices()
    I2 = np.eye(2, dtype=np.complex128)
    assert np.allclose(s1 @ s1, I2, atol=1e-15)
    assert np.allclose(s2 @ s2, I2, atol=1e-15)
    assert np.allclose(s3 @ s3, I2, atol=1e-15)


def test_pauli_matrices_anticommute():
    """{sigma^i, sigma^j} = 2 delta^{ij} I."""
    s = pauli_matrices()
    I2 = np.eye(2, dtype=np.complex128)
    for i in range(3):
        for j in range(3):
            ac = s[i] @ s[j] + s[j] @ s[i]
            target = 2.0 * (1.0 if i == j else 0.0) * I2
            assert np.allclose(ac, target, atol=1e-15)


def test_pauli_matrices_commutator_identity():
    """[sigma^i, sigma^j] = 2 i epsilon^{ijk} sigma^k."""
    s1, s2, s3 = pauli_matrices()
    eps_check = [
        (s1, s2, s3),  # [s1, s2] = 2i s3
        (s2, s3, s1),  # [s2, s3] = 2i s1
        (s3, s1, s2),  # [s3, s1] = 2i s2
    ]
    for a, b, c in eps_check:
        comm = a @ b - b @ a
        assert np.allclose(comm, 2j * c, atol=1e-15)


def test_gamma_anticommutation_west_coast():
    """{g^mu, g^nu} = 2 eta^{mu nu} I_4 with eta = diag(+1, -1, -1, -1)."""
    g = gamma_chiral()
    ok, residual = verify_gamma_anticommutation(g)
    assert ok, f"Clifford anticommutation failed at residual {residual:.3e}"
    assert residual == 0.0, (
        f"expected exact (zero) anticommutation residual in fixed Pauli "
        f"basis, got {residual:.3e}"
    )


def test_gamma_zero_squared_plus_identity():
    """(g^0)^2 = +I_4 (timelike, West-coast signature +)."""
    g = gamma_chiral()
    I4 = np.eye(4, dtype=np.complex128)
    assert np.allclose(g.g0 @ g.g0, I4, atol=1e-15)


def test_gamma_spatial_squared_minus_identity():
    """(g^i)^2 = -I_4 for i = 1, 2, 3 (spacelike, West-coast signature -)."""
    g = gamma_chiral()
    I4 = np.eye(4, dtype=np.complex128)
    for gi in (g.g1, g.g2, g.g3):
        assert np.allclose(gi @ gi, -I4, atol=1e-15)


def test_gamma5_squared_plus_identity():
    """(g^5)^2 = +I_4."""
    g = gamma_chiral()
    I4 = np.eye(4, dtype=np.complex128)
    assert np.allclose(g.g5 @ g.g5, I4, atol=1e-15)


def test_gamma5_anticommutes_with_all_gamma_mu():
    """{g^5, g^mu} = 0 for mu = 0, 1, 2, 3."""
    g = gamma_chiral()
    Z4 = np.zeros((4, 4), dtype=np.complex128)
    for gm in g.list_spacetime():
        ac = g.g5 @ gm + gm @ g.g5
        assert np.allclose(ac, Z4, atol=1e-15)


def test_gamma5_is_diagonal_chirality_grading():
    """g^5 = diag(-I_2, +I_2) (Peskin-Schroeder chiral basis convention)."""
    g = gamma_chiral()
    expected = np.diag([-1.0, -1.0, +1.0, +1.0]).astype(np.complex128)
    assert np.allclose(g.g5, expected, atol=1e-15)


def test_gamma_zero_is_chirality_swap_off_diagonal():
    """g^0 = [[0, I_2], [I_2, 0]] (off-diagonal swap in chiral basis)."""
    g = gamma_chiral()
    I2 = np.eye(2)
    Z2 = np.zeros((2, 2))
    expected = np.block([[Z2, I2], [I2, Z2]]).astype(np.complex128)
    assert np.allclose(g.g0, expected, atol=1e-15)


def test_gamma_matrices_hermitian_spacelike():
    """g^0 is Hermitian; g^i are anti-Hermitian in West-coast convention.

    In West-coast (chiral basis): (g^0)^dagger = g^0; (g^i)^dagger = -g^i.
    Equivalently, g^mu dagger = g^0 g^mu g^0 (standard identity).
    """
    g = gamma_chiral()
    assert np.allclose(g.g0, g.g0.conj().T, atol=1e-15)
    for gi in (g.g1, g.g2, g.g3):
        assert np.allclose(gi.conj().T, -gi, atol=1e-15)


# ===========================================================================
# 2. Temporal grid and identity_temporal
# ===========================================================================


def test_temporal_grid_symmetric():
    """Grid on [-T, T] with N_t = 21 is symmetric around 0 and includes 0."""
    grid = temporal_grid(1.0, 21)
    assert grid.shape == (21,)
    assert grid[0] == -1.0
    assert grid[-1] == 1.0
    assert grid[10] == 0.0
    assert np.allclose(grid + grid[::-1], 0.0, atol=1e-15)


def test_temporal_grid_N_t_1_is_zero():
    """N_t = 1 returns the singleton [0.0] (Riemannian-limit collapse)."""
    grid = temporal_grid(1.0, 1)
    assert grid.shape == (1,)
    assert grid[0] == 0.0


def test_temporal_grid_rejects_invalid_inputs():
    """T_max <= 0 or N_t < 1 raise ValueError."""
    with pytest.raises(ValueError):
        temporal_grid(0.0, 11)
    with pytest.raises(ValueError):
        temporal_grid(-1.0, 11)
    with pytest.raises(ValueError):
        temporal_grid(1.0, 0)


def test_identity_temporal_is_identity():
    """identity_temporal(N_t) is I_{N_t}."""
    for N_t in (1, 5, 21):
        I = identity_temporal(N_t)
        assert I.shape == (N_t, N_t)
        assert np.allclose(I, np.eye(N_t), atol=1e-15)


# ===========================================================================
# 3. Spatial fundamental symmetry
# ===========================================================================


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_spatial_fundamental_symmetry_swaps_chirality(n_max):
    """J_spatial is a permutation matrix swapping chirality labels at fixed
    (n_fock, l, two_m_j).
    """
    basis = full_dirac_basis(n_max)
    J = spatial_fundamental_symmetry(basis)

    # Permutation-matrix check: each column has exactly one 1, rest 0.
    col_sums = np.abs(J).sum(axis=0)
    assert np.allclose(col_sums, 1.0)

    # Explicit chirality-swap check: J |chi=+1, l, m_j> = |chi=-1, l, m_j>.
    label_to_idx = {b: i for i, b in enumerate(basis)}
    for j_idx, label in enumerate(basis):
        from geovac.full_dirac_operator_system import FullDiracLabel
        flipped = FullDiracLabel(
            n_fock=label.n_fock,
            l=label.l,
            two_m_j=label.two_m_j,
            chirality=-label.chirality,
        )
        i_idx = label_to_idx[flipped]
        assert J[i_idx, j_idx] == 1.0
        # All other entries in column j_idx are zero
        col = J[:, j_idx].copy()
        col[i_idx] = 0.0
        assert np.allclose(col, 0.0)


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_spatial_fundamental_symmetry_squared_to_identity(n_max):
    """J_spatial^2 = I (involution from chirality swap)."""
    basis = full_dirac_basis(n_max)
    J = spatial_fundamental_symmetry(basis)
    I = np.eye(len(basis), dtype=np.complex128)
    assert np.allclose(J @ J, I, atol=1e-15)


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_spatial_fundamental_symmetry_hermitian(n_max):
    """J_spatial is Hermitian (real symmetric permutation matrix)."""
    basis = full_dirac_basis(n_max)
    J = spatial_fundamental_symmetry(basis)
    assert np.allclose(J, J.conj().T, atol=1e-15)


# ===========================================================================
# 4. KreinSpace Krein axioms at n_max in {1, 2, 3} and N_t in {1, 11, 21}
# ===========================================================================


PARAM_GRID = [
    (n_max, N_t) for n_max in (1, 2, 3) for N_t in (1, 11, 21)
]


@pytest.mark.parametrize("n_max,N_t", PARAM_GRID)
def test_krein_J_squared_identity_bit_exact(n_max, N_t):
    """J^2 = +I bit-exact (Kronecker product of permutation matrices)."""
    K = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    ok, residual = K.verify_J_squared_identity(tol=1e-14)
    assert ok, f"J^2 != I at (n_max={n_max}, N_t={N_t}): residual={residual:.3e}"
    assert residual == 0.0, (
        f"expected bit-exact J^2 = I, got residual {residual:.3e}"
    )


@pytest.mark.parametrize("n_max,N_t", PARAM_GRID)
def test_krein_J_unitary_bit_exact(n_max, N_t):
    """J* J = I bit-exact (J is unitary as a Hilbert-space operator)."""
    K = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    ok, residual = K.verify_J_unitary(tol=1e-14)
    assert ok, f"J* J != I at (n_max={n_max}, N_t={N_t}): {residual:.3e}"
    assert residual == 0.0


@pytest.mark.parametrize("n_max,N_t", PARAM_GRID)
def test_krein_J_hermitian_bit_exact(n_max, N_t):
    """J = J* bit-exact."""
    K = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    ok, residual = K.verify_J_hermitian(tol=1e-14)
    assert ok, f"J != J* at (n_max={n_max}, N_t={N_t}): {residual:.3e}"
    assert residual == 0.0


@pytest.mark.parametrize("n_max,N_t", PARAM_GRID)
def test_krein_inner_product_hermitian(n_max, N_t):
    """<phi, psi>_K = conj(<psi, phi>_K) on randomly sampled vectors."""
    K = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    ok, residual = K.verify_krein_inner_product_hermitian(
        n_samples=20, seed=42, tol=1e-9
    )
    assert ok, (
        f"Krein form not Hermitian at (n_max={n_max}, N_t={N_t}): "
        f"max residual {residual:.3e}"
    )


# ===========================================================================
# 5. K^+ / K^- positive/negative split
# ===========================================================================


@pytest.mark.parametrize("n_max,N_t", PARAM_GRID)
def test_krein_split_dimensions_equal(n_max, N_t):
    """K^+ and K^- each have dimension dim / 2 (gamma^0 chirality-swap)."""
    K = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    d_plus, d_minus = K.split_dimensions()
    assert d_plus == K.dim // 2
    assert d_minus == K.dim // 2
    assert d_plus + d_minus == K.dim


@pytest.mark.parametrize("n_max,N_t", PARAM_GRID)
def test_krein_split_completeness_bit_exact(n_max, N_t):
    """P_+ + P_- = I, P_+ P_- = 0, P_+^2 = P_+, P_-^2 = P_- all bit-exact."""
    K = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    ok, residual = K.verify_split_completeness(tol=1e-14)
    assert ok, (
        f"K split completeness failed at (n_max={n_max}, N_t={N_t}): "
        f"residual {residual:.3e}"
    )
    assert residual == 0.0


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_krein_split_positive_definite_on_K_plus(n_max):
    """Krein form restricted to K^+ is positive-semidefinite: <v, v>_K >= 0."""
    K = KreinSpace(n_max=n_max, N_t=5, T_max=1.0)
    P_plus, _ = K.positive_negative_split()
    # P_plus has rank dim/2; take its columns as a (non-orthonormal)
    # spanning set of K^+, then evaluate <v, v>_K for each.
    # Equivalently, eigenvalues of P_plus * J * P_plus restricted to K^+
    # should all be +1 (since J|K^+ = +I).
    M = P_plus @ K.J @ P_plus  # restriction of J to K^+
    # eigenvalues should be either 0 or +1
    eigvals = np.linalg.eigvalsh(0.5 * (M + M.conj().T))
    # Eigenvalues split: dim/2 at 0 (orthogonal complement), dim/2 at +1.
    assert np.all(eigvals >= -1e-12), (
        f"Krein form on K^+ has negative eigenvalue: min = {eigvals.min():.3e}"
    )


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_krein_split_negative_definite_on_K_minus(n_max):
    """Krein form restricted to K^- is negative-semidefinite: <v, v>_K <= 0."""
    K = KreinSpace(n_max=n_max, N_t=5, T_max=1.0)
    _, P_minus = K.positive_negative_split()
    M = P_minus @ K.J @ P_minus  # = -P_minus (since J|K^- = -I)
    eigvals = np.linalg.eigvalsh(0.5 * (M + M.conj().T))
    # Eigenvalues split: dim/2 at 0 (orthogonal complement), dim/2 at -1.
    assert np.all(eigvals <= 1e-12), (
        f"Krein form on K^- has positive eigenvalue: max = {eigvals.max():.3e}"
    )


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_krein_J_eigenvalues_are_plus_minus_one(n_max):
    """J spectrum at (n_max, N_t = 5) is exactly the multiset {+1, -1}^dim/2."""
    K = KreinSpace(n_max=n_max, N_t=5, T_max=1.0)
    # J is Hermitian unitary, eigenvalues are real +/- 1.
    eigvals = np.sort(np.linalg.eigvalsh(0.5 * (K.J + K.J.conj().T)))
    half = K.dim // 2
    expected = np.concatenate([-np.ones(half), np.ones(half)])
    assert np.allclose(eigvals, expected, atol=1e-12)


# ===========================================================================
# 6. THE LOAD-BEARING FALSIFIER:
#    Riemannian limit at N_t = 1 recovers H_GV bit-identically
# ===========================================================================


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_riemannian_limit_check_passes_bit_identically(n_max):
    """LOAD-BEARING: Krein space at N_t = 1 equals H_GV^{n_max} bit-identically.

    This is the falsifier from Sprint L2-A Section 5.7 / L2-B brief.  If
    this fails, the Camporesi-Higuchi spatial spinor bundle is
    incompatible with the Cl(3, 1) gamma matrix embedding -- escalate.
    """
    K = KreinSpace(n_max=n_max, N_t=1, T_max=1.0)
    ok, details = K.riemannian_limit_check(tol=1e-14)
    assert ok, (
        f"RIEMANNIAN LIMIT CHECK FAILED at n_max={n_max}: details={details}. "
        f"This is a MAJOR STRUCTURAL FINDING; escalate per L2-B brief."
    )
    assert details["dim_match"]
    assert details["basis_match"]
    assert details["J_match_residual"] == 0.0
    assert details["time_slice_check_residual"] == 0.0


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_riemannian_limit_dim_matches_paper32_formula(n_max):
    """At N_t = 1: dim(K) = N_Dirac(n_max) = (2/3) n_max (n_max+1)(n_max+2).

    This is the Paper 32 def:H_GV formula for the Camporesi-Higuchi
    spinor space dimension.
    """
    K = KreinSpace(n_max=n_max, N_t=1, T_max=1.0)
    N_Dirac = (2 * n_max * (n_max + 1) * (n_max + 2)) // 3
    assert K.dim == N_Dirac
    assert K.dim == full_dirac_dim(n_max)


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_riemannian_limit_basis_labels_match(n_max):
    """At N_t = 1: K.basis_spatial == full_dirac_basis(n_max) label-by-label."""
    K = KreinSpace(n_max=n_max, N_t=1, T_max=1.0)
    ref = full_dirac_basis(n_max)
    assert len(K.basis_spatial) == len(ref)
    for i, (kb, rb) in enumerate(zip(K.basis_spatial, ref)):
        assert kb == rb, (
            f"basis label mismatch at index {i}: krein={kb}, ref={rb}"
        )


@pytest.mark.parametrize("n_max,N_t", PARAM_GRID)
def test_time_slice_block_equals_J_spatial(n_max, N_t):
    """For arbitrary N_t, the time-slice diagonal block of K.J equals J_spatial.

    This verifies the Kronecker structure J = J_spatial (x) I_{N_t}
    decouples cleanly into independent temporal slices at the level of
    the fundamental symmetry.
    """
    K = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    J_spatial_ref = spatial_fundamental_symmetry(full_dirac_basis(n_max))
    # In our Kronecker convention J = kron(J_spatial, I_{N_t}), the time-
    # slice block at temporal index k=0 is J[0::N_t, 0::N_t] = J_spatial.
    slice_k0 = K.J[0 :: N_t, 0 :: N_t]
    assert np.allclose(slice_k0, J_spatial_ref, atol=1e-15)
    # Bit-identical residual:
    assert float(np.linalg.norm(slice_k0 - J_spatial_ref)) == 0.0


# ===========================================================================
# 7. Constructor / API smoke tests
# ===========================================================================


def test_KreinSpace_constructor_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        KreinSpace(n_max=0, N_t=11)
    with pytest.raises(ValueError):
        KreinSpace(n_max=2, N_t=0)
    with pytest.raises(ValueError):
        KreinSpace(n_max=2, N_t=11, T_max=0.0)
    with pytest.raises(ValueError):
        KreinSpace(n_max=2, N_t=11, T_max=-1.0)


def test_KreinSpace_repr_contains_key_fields():
    K = KreinSpace(n_max=2, N_t=11, T_max=1.5)
    r = repr(K)
    assert "n_max=2" in r
    assert "N_t=11" in r
    assert "T_max=1.5" in r
    assert "dim_spatial=16" in r
    assert "dim=176" in r


def test_krein_inner_product_shape_check():
    K = KreinSpace(n_max=1, N_t=5, T_max=1.0)
    psi = np.ones(K.dim, dtype=np.complex128)
    phi = np.ones(K.dim, dtype=np.complex128)
    val = K.krein_inner_product(psi, phi)
    # <1, J 1> = sum of all J entries; J permutation -> sum = dim
    assert np.isclose(val.real, K.dim)
    # Wrong-shape inputs raise
    with pytest.raises(ValueError):
        K.krein_inner_product(psi[:-1], phi)
    with pytest.raises(ValueError):
        K.krein_inner_product(psi, phi[:-1])


def test_positive_negative_split_raises_if_J_not_involution():
    """The split builder sanity-checks J^2 = I."""
    K = KreinSpace(n_max=1, N_t=5, T_max=1.0)
    # Deliberately corrupt J (we should not be able to call split)
    K.J = K.J + 1e-3 * np.eye(K.dim, dtype=np.complex128)
    with pytest.raises(ValueError):
        K.positive_negative_split(tol=1e-12)
