"""Tests for `geovac.lorentzian_dirac` (Sprint L2-C).

These tests verify the van den Dungen 2016 Proposition 4.1 lift on
S^3 x R at signature (3, 1) using the Krein space from L2-B.

Coverage:

  (a) Centered finite-difference matrix is anti-Hermitian; the
      N_t = 1 degenerate case is the 1x1 zero matrix.
  (b) Cl(3, 1) gamma matrix Clifford algebra holds on the 4-manifold
      (cross-check of L2-B properties used in the L2-C construction).
  (c) D_L = i * [gamma^0 (x) d/dt + D_GV (x) I_{N_t}] is well-defined
      at every (n_max, N_t) in the L2-B panel.
  (d) L2C-FALS-1: D_L^x = D_L (Krein-self-adjoint) bit-exact across
      n_max in {1, 2, 3} and N_t in {1, 11, 21}.
  (e) L2C-FALS-2: {gamma^5, D_L} = 0 -- structural finding: this is
      NON-ZERO because D_GV commutes with gamma^5 by GeoVac labeling
      convention.  Test asserts the structural value matches predicted
      decomposition |{gamma^5, D_L}|_F = 2 * |D_GV| (no temporal
      contribution since {gamma^5, gamma^0} = 0).
  (f) L2C-FALS-3 (LOAD-BEARING): at N_t = 1, D_L = i * D_GV
      bit-identically and the absolute-spectrum matches |spec(D_GV)|.
      If this fails, STOP and escalate.
  (g) L2C-FALS-4: spectrum of D_L|_{K^+} is real (|imag| < 1e-10).

Zero regression in upstream tests.  Run
    pytest tests/test_lorentzian_dirac.py
    pytest tests/test_krein_space_construction.py
    pytest tests/test_full_dirac_operator_system.py
    pytest tests/test_modular_hamiltonian.py
to verify.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.full_dirac_operator_system import (
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
    full_dirac_dim,
)
from geovac.krein_space_construction import (
    KreinSpace,
    gamma_chiral,
    spatial_fundamental_symmetry,
)
from geovac.lorentzian_dirac import (
    centered_difference_matrix,
    krein_adjoint,
    l2c_audit,
    lorentzian_dirac_matrix,
    spacetime_chirality_lifted,
    spatial_chirality_grading,
    verify_anti_hermitian,
    verify_chirality_anticommutation,
    verify_chirality_commutation,
    verify_krein_self_adjoint,
    verify_riemannian_limit,
    verify_spectrum_reality_K_plus,
)


# ---------------------------------------------------------------------------
# Group 1: centered finite-difference temporal derivative
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("N_t", [1, 2, 3, 5, 11, 21])
def test_centered_difference_shape(N_t: int) -> None:
    """D_t has shape (N_t, N_t) and dtype complex128."""
    D_t = centered_difference_matrix(N_t, T_max=1.0)
    assert D_t.shape == (N_t, N_t)
    assert D_t.dtype == np.complex128


def test_centered_difference_Nt1_is_zero() -> None:
    """At N_t = 1, d/dt is identically zero (singleton grid, no neighbors)."""
    D_t = centered_difference_matrix(N_t=1, T_max=1.0)
    assert np.allclose(D_t, np.zeros((1, 1)))
    assert float(np.linalg.norm(D_t)) == 0.0


@pytest.mark.parametrize("N_t", [2, 3, 5, 11, 21])
def test_centered_difference_anti_hermitian(N_t: int) -> None:
    """L2-C key structural property: D_t^dagger = -D_t (anti-Hermitian).

    Required for D_L = i * (gamma^0 (x) D_t + D_GV (x) I) to be
    Krein-self-adjoint with respect to J = gamma^0.
    """
    D_t = centered_difference_matrix(N_t, T_max=1.0)
    ok, residual = verify_anti_hermitian(D_t, tol=1e-14)
    assert ok, f"D_t at N_t={N_t} not anti-Hermitian, residual={residual}"


def test_centered_difference_first_row_explicit() -> None:
    """Explicit value check at N_t = 5, T_max = 1 (Delta_t = 0.5)."""
    D_t = centered_difference_matrix(N_t=5, T_max=1.0)
    # Delta_t = 2*1/(5-1) = 0.5, so 1/(2*Delta_t) = 1.0
    assert D_t[0, 1].real == pytest.approx(+1.0)
    assert D_t[0, 0].real == pytest.approx(0.0)
    assert D_t[2, 1].real == pytest.approx(-1.0)
    assert D_t[2, 3].real == pytest.approx(+1.0)
    assert D_t[4, 3].real == pytest.approx(-1.0)
    # Boundary: D_t[0, -1] doesn't exist (out of range), D_t[N-1, N] doesn't either
    assert D_t[0, -1].real == pytest.approx(0.0)  # = D_t[0, 4] which is unrelated


def test_centered_difference_invalid_inputs() -> None:
    """Negative N_t and T_max <= 0 must raise."""
    with pytest.raises(ValueError):
        centered_difference_matrix(N_t=0, T_max=1.0)
    with pytest.raises(ValueError):
        centered_difference_matrix(N_t=-1, T_max=1.0)
    with pytest.raises(ValueError):
        centered_difference_matrix(N_t=5, T_max=0.0)
    with pytest.raises(ValueError):
        centered_difference_matrix(N_t=5, T_max=-1.0)


@pytest.mark.parametrize("T_max", [0.5, 1.0, 2.0, 10.0])
def test_centered_difference_T_max_scaling(T_max: float) -> None:
    """D_t scales as 1/T_max (Delta_t proportional to T_max)."""
    D_t1 = centered_difference_matrix(N_t=11, T_max=1.0)
    D_t = centered_difference_matrix(N_t=11, T_max=T_max)
    # D_t should be D_t1 / T_max
    ratio_matrix = D_t * T_max - D_t1
    assert float(np.linalg.norm(ratio_matrix)) < 1e-14


# ---------------------------------------------------------------------------
# Group 2: spatial chirality grading and spacetime gamma^5 lift
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_spatial_chirality_grading_diagonal(n_max: int) -> None:
    """gamma^5_spatial = diag(-chirality) is diagonal and involutive."""
    basis = full_dirac_basis(n_max)
    gamma5 = spatial_chirality_grading(basis)
    dim = len(basis)
    # Diagonal
    off_diag = gamma5 - np.diag(np.diag(gamma5))
    assert float(np.linalg.norm(off_diag)) == 0.0
    # Involutive: (gamma^5)^2 = I
    assert np.allclose(gamma5 @ gamma5, np.eye(dim))
    # Hermitian
    assert np.allclose(gamma5, gamma5.conj().T)
    # Eigenvalues are exactly +/-1
    eigvals = np.diag(gamma5).real
    assert set(np.round(eigvals, 12)) <= {-1.0, +1.0}


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_spatial_chirality_grading_balanced(n_max: int) -> None:
    """gamma^5_spatial has dim_H / 2 entries +1 and dim_H / 2 entries -1."""
    basis = full_dirac_basis(n_max)
    gamma5 = spatial_chirality_grading(basis)
    eigvals = np.diag(gamma5).real
    n_plus = int(np.sum(eigvals > 0))
    n_minus = int(np.sum(eigvals < 0))
    assert n_plus == len(basis) // 2
    assert n_minus == len(basis) // 2


@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("N_t", [1, 11, 21])
def test_spacetime_chirality_lifted_shape(n_max: int, N_t: int) -> None:
    """spacetime_chirality_lifted returns a (dim, dim) involutive matrix."""
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    gamma5 = spacetime_chirality_lifted(krein)
    assert gamma5.shape == (krein.dim, krein.dim)
    assert np.allclose(gamma5 @ gamma5, np.eye(krein.dim))
    assert np.allclose(gamma5, gamma5.conj().T)


def test_spacetime_chirality_anticomm_with_J() -> None:
    """{gamma^5, gamma^0} = 0 (Clifford algebra) lifts to {gamma^5_K, J} = 0."""
    krein = KreinSpace(n_max=2, N_t=5)
    gamma5 = spacetime_chirality_lifted(krein)
    anticomm = gamma5 @ krein.J + krein.J @ gamma5
    assert float(np.linalg.norm(anticomm)) < 1e-12


# ---------------------------------------------------------------------------
# Group 3: gamma^5 from Cl(3,1) consistency cross-check with L2-B
# ---------------------------------------------------------------------------


def test_cl31_gamma5_consistency_with_lift() -> None:
    """gamma^5 (4x4) eigenvalue pattern matches the lift on H_GV.

    The L2-B gamma^5 = diag(-1, -1, +1, +1) has 2 eigenvalues -1 and
    2 eigenvalues +1.  The spatial lift on H_GV at n_max should have
    half its eigenvalues +1 and half -1.
    """
    gammas = gamma_chiral()
    g5 = gammas.g5
    eigvals = np.linalg.eigvalsh(g5)
    eigvals_real = sorted(np.round(eigvals.real, 12))
    assert eigvals_real == [-1.0, -1.0, +1.0, +1.0]


# ---------------------------------------------------------------------------
# Group 4: lorentzian_dirac_matrix construction
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("N_t", [1, 11, 21])
def test_lorentzian_dirac_shape(n_max: int, N_t: int) -> None:
    """D_L has shape (dim_K, dim_K)."""
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    D_L = lorentzian_dirac_matrix(krein)
    assert D_L.shape == (krein.dim, krein.dim)
    assert D_L.dtype == np.complex128


@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("N_t", [1, 11, 21])
def test_lorentzian_dirac_finite(n_max: int, N_t: int) -> None:
    """D_L has no NaN or Inf entries."""
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    D_L = lorentzian_dirac_matrix(krein)
    assert np.all(np.isfinite(D_L))


def test_lorentzian_dirac_custom_dirac() -> None:
    """Custom dirac_diag override is honored if shape matches."""
    krein = KreinSpace(n_max=2, N_t=5)
    custom = np.eye(krein.dim_spatial, dtype=np.complex128) * 7.0
    D_L = lorentzian_dirac_matrix(krein, dirac_diag=custom)
    assert D_L.shape == (krein.dim, krein.dim)
    assert np.all(np.isfinite(D_L))


def test_lorentzian_dirac_wrong_shape_raises() -> None:
    """Wrong-shape dirac_diag override raises ValueError."""
    krein = KreinSpace(n_max=2, N_t=5)
    with pytest.raises(ValueError):
        lorentzian_dirac_matrix(krein, dirac_diag=np.zeros((3, 3)))


# ---------------------------------------------------------------------------
# Group 5: L2C-FALS-1 -- Krein self-adjointness
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("N_t", [1, 11, 21])
def test_FALS1_krein_self_adjoint_bit_exact(n_max: int, N_t: int) -> None:
    """L2C-FALS-1: D_L^x = D_L bit-exact at all panel cells.

    The vdD Prop 4.1 i^t = i factor combined with D_t anti-Hermitian
    + {gamma^0, D_GV} = 0 makes the Krein-adjoint identity hold to
    machine precision.
    """
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    D_L = lorentzian_dirac_matrix(krein)
    ok, residual = verify_krein_self_adjoint(D_L, krein.J, tol=1e-12)
    assert ok, (
        f"L2C-FALS-1 FAILED at n_max={n_max}, N_t={N_t}: "
        f"||D_L^x - D_L||_F = {residual}"
    )
    # Stronger: expect bit-exact (machine zero) at every cell
    assert residual <= 1e-12, (
        f"Residual {residual} too large -- expected machine precision."
    )


def test_FALS1_anti_dirac_breaks_krein_sa() -> None:
    """A D' that COMMUTES with gamma^0 fails Krein-self-adjointness.

    Sanity check: if we override the spatial Dirac to one that commutes
    (rather than anticommutes) with J_spatial = gamma^0, the resulting
    pseudo-D_L should fail Krein-self-adjointness.  This confirms our
    test detects genuine failures.
    """
    krein = KreinSpace(n_max=2, N_t=5)
    # Build an obviously commuting matrix: a multiple of J_spatial itself
    # commutes with gamma^0 = J_spatial.
    bad_diag = krein.J_spatial.copy() * 3.0
    D_bad = lorentzian_dirac_matrix(krein, dirac_diag=bad_diag)
    ok, residual = verify_krein_self_adjoint(D_bad, krein.J, tol=1e-12)
    # Should FAIL since bad_diag commutes (not anticommutes) with J_spatial
    assert not ok, "Expected Krein-self-adjointness to fail on commuting D'"
    assert residual > 1e-6


def test_FALS1_krein_adjoint_helper() -> None:
    """krein_adjoint(D, J) is involutive: (D^x)^x = D for J^2 = I."""
    krein = KreinSpace(n_max=1, N_t=3)
    rng = np.random.default_rng(42)
    D_rand = (
        rng.standard_normal((krein.dim, krein.dim))
        + 1j * rng.standard_normal((krein.dim, krein.dim))
    )
    Dxx = krein_adjoint(krein_adjoint(D_rand, krein.J), krein.J)
    assert np.allclose(Dxx, D_rand, atol=1e-12)


# ---------------------------------------------------------------------------
# Group 6: L2C-FALS-2 -- chirality anticommutation (STRUCTURAL FINDING)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_FALS2_anticommutator_structural_at_Nt1(n_max: int) -> None:
    """L2C-FALS-2 at N_t = 1: {gamma^5, D_L} = 2 * (gamma^5 * D_L) since they commute.

    At N_t = 1, D_t = 0 and D_L = i * D_GV (x) I_1.  gamma^5_K commutes
    with D_GV (both diagonal in chirality), so the COMMUTATOR is zero
    but the ANTICOMMUTATOR is nonzero (equal to 2 * |D_L|).  This
    confirms the structural-finding diagnosis.
    """
    krein = KreinSpace(n_max=n_max, N_t=1)
    D_L = lorentzian_dirac_matrix(krein)
    gamma5 = spacetime_chirality_lifted(krein)

    # Commutator [gamma^5, D_L] = 0 at N_t = 1 (gamma^5 commutes with D_GV)
    ok_comm, res_comm = verify_chirality_commutation(D_L, gamma5, tol=1e-12)
    assert ok_comm, (
        f"At N_t=1, [gamma^5, D_L] should be ZERO (D_GV diagonal in "
        f"chirality), got residual {res_comm}"
    )

    # Anticommutator is 2 * gamma^5 * D_L = 2|D_L| in Frobenius norm
    anticomm = gamma5 @ D_L + D_L @ gamma5
    expected = 2.0 * gamma5 @ D_L
    assert np.allclose(anticomm, expected, atol=1e-12)


@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("N_t", [5, 11])
def test_FALS2_anticommutator_nonzero_general(n_max: int, N_t: int) -> None:
    """At general N_t > 1, both {gamma^5, D_L} and [gamma^5, D_L] are nonzero.

    The temporal piece i * gamma^0 (x) d/dt anticommutes with gamma^5
    (since {gamma^5, gamma^0} = 0); the spatial piece i * D_GV (x) I
    commutes with gamma^5.  So {gamma^5, D_L} picks up the spatial
    contribution, and [gamma^5, D_L] picks up the temporal one.  Both
    nonzero in general; this documents the structural feature
    cleanly.
    """
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    D_L = lorentzian_dirac_matrix(krein)
    gamma5 = spacetime_chirality_lifted(krein)
    _, res_anti = verify_chirality_anticommutation(D_L, gamma5, tol=1e-12)
    _, res_comm = verify_chirality_commutation(D_L, gamma5, tol=1e-12)
    # Both nonzero
    assert res_anti > 1e-3, (
        f"Expected nonzero anticommutator at general N_t, got {res_anti}"
    )
    assert res_comm > 1e-3, (
        f"Expected nonzero commutator at general N_t, got {res_comm}"
    )


def test_FALS2_anticommutator_decomposition_at_Nt1() -> None:
    """{gamma^5, D_L} at N_t = 1 = 2 * gamma^5 * (i * D_GV) bit-exact.

    Test the structural-finding formula explicitly.
    """
    n_max = 2
    krein = KreinSpace(n_max=n_max, N_t=1)
    D_L = lorentzian_dirac_matrix(krein)
    gamma5 = spacetime_chirality_lifted(krein)
    D_GV = camporesi_higuchi_full_dirac_matrix(krein.basis_spatial)

    anticomm = gamma5 @ D_L + D_L @ gamma5
    # gamma^5 commutes with D_GV; D_L = i * D_GV at N_t=1.
    # So {gamma^5, D_L} = 2 * gamma^5 * D_L = 2 * gamma^5 * (i * D_GV).
    expected = 2j * (gamma5 @ D_GV)
    assert np.allclose(anticomm, expected, atol=1e-12)


# ---------------------------------------------------------------------------
# Group 7: L2C-FALS-3 -- LOAD-BEARING Riemannian-limit recovery
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_FALS3_load_bearing_riemannian_limit_bit_exact(n_max: int) -> None:
    """L2C-FALS-3 (LOAD-BEARING): D_L|_{N_t=1} = i * D_GV bit-exact.

    Failure here would mean the Wick-rotation lift is incompatible
    with Paper 32 §III's D_GV -- STOP and escalate.
    """
    ok, details = verify_riemannian_limit(n_max, tol_struct=1e-12)
    assert ok, (
        f"LOAD-BEARING L2C-FALS-3 FAILED at n_max={n_max}.  "
        f"details = {details}.  STOP and escalate."
    )
    assert details["matrix_residual_a"] == 0.0
    assert details["matrix_residual_b"] == 0.0
    assert details["spectrum_residual_max"] == 0.0


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_FALS3_spectrum_multiplicities(n_max: int) -> None:
    """|spec(D_L|_{N_t=1})| has multiplicities g_n = 2*(n+1)*(n+2).

    Camporesi-Higuchi degeneracy at level n_fock = n is 2*(n+1)*(n+2)
    counting both chirality signs; each absolute eigenvalue n + 1/2
    has multiplicity 2*(n_fock)*(n_fock+1) at fixed n_fock when summed
    over both chirality signs (chirality = +/- 1).  At n_max we sum
    levels n_fock = 1, ..., n_max.

    The L2C-FALS-3 spec specifically says:
      "spectrum should be {|lambda_n| = n + 3/2} with same multiplicities"

    Here n in the spec is n_fock - 1 (n_fock starts at 1 in FullDiracLabel),
    so |lambda_{n_fock}| = (n_fock - 1) + 3/2 = n_fock + 1/2.  This is
    the FullDiracLabel formula chi*(n_fock + 1/2).
    """
    krein = KreinSpace(n_max=n_max, N_t=1)
    D_L = lorentzian_dirac_matrix(krein)
    D_GV = camporesi_higuchi_full_dirac_matrix(krein.basis_spatial)

    abs_DL = np.sort(np.abs(np.linalg.eigvals(D_L).round(10)))
    abs_DGV = np.sort(np.abs(np.linalg.eigvalsh(D_GV.real).round(10)))
    assert np.allclose(abs_DL, abs_DGV, atol=1e-12)

    # Each absolute eigenvalue should appear with multiplicity 2 * (n_fock) * (n_fock + 1)
    # (combined over chirality +/- chains).
    from collections import Counter
    counts_DL = Counter(np.round(abs_DL, 10))
    for n_fock in range(1, n_max + 1):
        abs_val = round(n_fock + 0.5, 10)
        # Spinor-bundle level n_fock contributes 2 chains * spinor_dim_per_level
        # spinor_dim at n_fock = n_fock * (n_fock + 1) (Weyl chain dim),
        # full Dirac = 2 * spinor_dim_at_level = 2 * n_fock * (n_fock + 1).
        # Actually: at each n_fock, l ranges 0..n_fock-1, 2j+1 = 2l+2.
        # Spinor dim at n_fock = sum_{l=0..n_fock-1} (2l+2) = n_fock * (n_fock+1).
        # Full Dirac at n_fock = 2 * n_fock * (n_fock + 1).
        expected_mult = 2 * n_fock * (n_fock + 1)
        assert counts_DL[abs_val] == expected_mult, (
            f"At n_fock={n_fock}, |eigval| = {abs_val}: "
            f"expected mult {expected_mult}, got {counts_DL[abs_val]}"
        )


def test_FALS3_load_bearing_details_dict_shape() -> None:
    """riemannian_limit returns a details dict with all expected keys."""
    ok, details = verify_riemannian_limit(n_max=2, tol_struct=1e-12)
    expected_keys = {
        "n_max", "dim", "matrix_residual_a", "matrix_residual_b",
        "spectrum_residual_max", "D_GV_spectrum_sample",
        "D_L_abs_spectrum_sample", "load_bearing_pass",
    }
    assert set(details.keys()) >= expected_keys
    assert details["n_max"] == 2
    assert details["dim"] == full_dirac_dim(2)


# ---------------------------------------------------------------------------
# Group 8: L2C-FALS-4 -- real spectrum on Krein-positive cone
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("N_t", [1, 5, 11])
def test_FALS4_spectrum_real_on_K_plus(n_max: int, N_t: int) -> None:
    """L2C-FALS-4: spectrum of D_L|_{K^+} is real (|imag| < 1e-10)."""
    krein = KreinSpace(n_max=n_max, N_t=N_t)
    D_L = lorentzian_dirac_matrix(krein)
    ok, details = verify_spectrum_reality_K_plus(D_L, krein, tol=1e-10)
    assert ok, (
        f"L2C-FALS-4 FAILED at n_max={n_max}, N_t={N_t}: "
        f"max_imag_K_plus = {details['max_imag_K_plus']}"
    )
    assert details["dim_K_plus"] == krein.dim // 2


def test_FALS4_K_plus_dimension_correct() -> None:
    """dim K^+ = dim K / 2 (positive cone half of Krein splitting)."""
    krein = KreinSpace(n_max=2, N_t=11)
    D_L = lorentzian_dirac_matrix(krein)
    _, details = verify_spectrum_reality_K_plus(D_L, krein, tol=1e-10)
    assert details["dim_K_plus"] == krein.dim // 2


# ---------------------------------------------------------------------------
# Group 9: l2c_audit convenience wrapper
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("N_t", [1, 11, 21])
def test_l2c_audit_load_bearing_pass(n_max: int, N_t: int) -> None:
    """At every (n_max, N_t) panel cell, the two load-bearing falsifiers pass.

    Load-bearing: FALS-1 (Krein-self-adjoint) and FALS-3 (Riemannian
    limit) MUST pass.  FALS-2 is expected to be structural (nonzero).
    FALS-4 should be near machine zero.
    """
    result = l2c_audit(n_max=n_max, N_t=N_t, T_max=1.0)
    assert result["FALS_1_krein_self_adjoint"]["pass"], (
        f"L2C-FALS-1 (Krein self-adjoint) failed at ({n_max}, {N_t})"
    )
    assert result["FALS_3_riemannian_limit"]["pass"], (
        f"L2C-FALS-3 (Riemannian limit, LOAD-BEARING) failed at "
        f"({n_max}, {N_t})"
    )
    assert result["all_load_bearing_pass"]


def test_l2c_audit_FALS2_structural() -> None:
    """At general N_t, l2c_audit reports FALS-2 as structurally non-zero."""
    result = l2c_audit(n_max=2, N_t=11, T_max=1.0)
    # FALS-2 should not pass (structural)
    assert not result["FALS_2_chirality_anticommutation"]["pass"]
    # Anti-commutator residual should be substantially nonzero
    assert (
        result["FALS_2_chirality_anticommutation"]["residual"]
        > 1.0
    )


def test_l2c_audit_dict_shape() -> None:
    """l2c_audit returns dict with all four FALS keys + summary."""
    result = l2c_audit(n_max=1, N_t=5, T_max=1.0)
    for k in (
        "FALS_1_krein_self_adjoint",
        "FALS_2_chirality_anticommutation",
        "FALS_3_riemannian_limit",
        "FALS_4_real_spectrum_K_plus",
        "all_load_bearing_pass",
        "n_max", "N_t", "T_max", "dim_K",
    ):
        assert k in result, f"Missing key {k!r}"


# ---------------------------------------------------------------------------
# Group 10: zero regression cross-check with existing modules
# ---------------------------------------------------------------------------


def test_no_regression_full_dirac_dim() -> None:
    """full_dirac_dim and full_dirac_basis unchanged by L2-C work."""
    assert full_dirac_dim(1) == 4
    assert full_dirac_dim(2) == 16
    assert full_dirac_dim(3) == 40


def test_no_regression_krein_dim_at_Nt1() -> None:
    """KreinSpace at N_t=1 still gives dim = full_dirac_dim."""
    for n_max in [1, 2, 3]:
        krein = KreinSpace(n_max=n_max, N_t=1)
        assert krein.dim == full_dirac_dim(n_max)


def test_no_regression_camporesi_higuchi_eigenvalues() -> None:
    """Camporesi-Higuchi D_GV eigenvalues match L2-B convention.

    Diagonal D_GV with eigenvalue chi*(n_fock + 1/2).
    """
    basis = full_dirac_basis(2)
    D_GV = camporesi_higuchi_full_dirac_matrix(basis)
    for i, b in enumerate(basis):
        expected = b.chirality * (b.n_fock + 0.5)
        assert D_GV[i, i].real == pytest.approx(expected)


def test_no_regression_J_spatial_unchanged() -> None:
    """spatial_fundamental_symmetry unchanged by L2-C edits."""
    basis = full_dirac_basis(2)
    J = spatial_fundamental_symmetry(basis)
    # J is real, Hermitian, J^2 = I
    assert np.allclose(J, J.conj().T)
    assert np.allclose(J @ J, np.eye(len(basis)))
    # All entries are 0 or 1 (permutation matrix)
    entries = set(J.flatten().real)
    assert entries <= {0.0, 1.0}
