"""Tests for geovac/chirality_grading.py.

Sprint G3-A: Z_2 chirality grading γ_GV on the full-Dirac truncated
operator system.

Test panel
==========

(A) Constructor / shape sanity:
    - γ has the right dimension at n_max ∈ {1, 2, 3}.
    - γ matches sigma_x ⊗ I (or sigma_z ⊗ I) explicitly at n_max = 1.

(B) Z_2 grading axioms:
    - γ² = I exactly.
    - γ is Hermitian.
    - γ is unitary.

(C) Anticommutation with the truthful CH Dirac:
    - {γ_sigma_x, D_truthful} = 0 exactly.
    - {γ_sigma_z, D_truthful} ≠ 0 (negative control: sigma_z commutes with D).

(D) Anticommutation with the offdiag CH Dirac:
    - {γ_sigma_x, D_offdiag} ≠ 0 (the within-chirality ladder is
      symmetric across blocks → {γ,X} = 2X ≠ 0 generically).

(E) Operator-system preservation γ M γ ⊆ O:
    - For sigma_x convention, γ M γ = M exactly for every generator
      (since M = I_2 ⊗ M_Weyl block-diagonal commutes with sigma_x ⊗ I).
    - All generators pass the contains() check at every n_max.

(F) Commutation with multipliers [γ, M] = 0:
    - For sigma_x convention, [γ, M] = 0 exactly for every M.

(G) Audit driver returns structured dict.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.chirality_grading import (
    ChiralityGrading,
    audit_chirality_grading,
    build_gamma_GV,
)
from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
    camporesi_higuchi_offdiag_dirac_matrix,
    full_dirac_basis,
    full_dirac_dim,
)
from geovac.spinor_operator_system import spinor_dim


# ---------------------------------------------------------------------------
# (A) Constructor / shape sanity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_dim_matches_full_dirac(n_max):
    """γ_GV has the right shape: (full_dirac_dim, full_dirac_dim)."""
    g = build_gamma_GV(n_max)
    expected_dim = 2 * spinor_dim(n_max)
    assert g.matrix.shape == (expected_dim, expected_dim)
    assert g.dim == expected_dim


def test_sigma_x_matches_swap_at_nmax_1():
    """At n_max = 1, dim_weyl = 2, so γ_sigma_x is a 4x4 swap of two 2x2 blocks."""
    g = build_gamma_GV(1, convention="sigma_x")
    # Expected: kron(sigma_x, I_2) = [[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]]
    expected = np.array(
        [
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [1, 0, 0, 0],
            [0, 1, 0, 0],
        ],
        dtype=np.complex128,
    )
    assert np.allclose(g.matrix, expected)


def test_sigma_z_matches_diag_at_nmax_1():
    """At n_max = 1, γ_sigma_z = diag(+1, +1, -1, -1)."""
    g = build_gamma_GV(1, convention="sigma_z")
    expected = np.diag([1, 1, -1, -1]).astype(np.complex128)
    assert np.allclose(g.matrix, expected)


def test_invalid_n_max_raises():
    with pytest.raises(ValueError):
        build_gamma_GV(0)


def test_invalid_convention_raises():
    with pytest.raises(ValueError):
        build_gamma_GV(2, convention="bogus")


# ---------------------------------------------------------------------------
# (B) Z_2 grading axioms
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("convention", ["sigma_x", "sigma_z"])
def test_gamma_squared_is_identity(n_max, convention):
    """γ² = I exactly (machine zero) for both conventions."""
    g = build_gamma_GV(n_max, convention=convention)
    assert g.verify_grading_squared() < 1e-14


@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("convention", ["sigma_x", "sigma_z"])
def test_gamma_hermitian(n_max, convention):
    """γ = γ† exactly."""
    g = build_gamma_GV(n_max, convention=convention)
    assert g.verify_hermitian() < 1e-14


@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("convention", ["sigma_x", "sigma_z"])
def test_gamma_unitary(n_max, convention):
    """γγ† = I exactly."""
    g = build_gamma_GV(n_max, convention=convention)
    assert g.verify_unitary() < 1e-14


# ---------------------------------------------------------------------------
# (C) Anticommutation with the truthful CH Dirac
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_sigma_x_anticommutes_truthful(n_max):
    """{γ_sigma_x, D_truthful} = 0 exactly. This is the load-bearing axiom."""
    g = build_gamma_GV(n_max, convention="sigma_x")
    res = g.verify_anticommutes_dirac("truthful")
    assert res < 1e-14, (
        f"sigma_x convention must give {{γ, D_truthful}} = 0 at n_max={n_max}, "
        f"got {res:.3e}"
    )


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_sigma_z_does_not_anticommute_truthful(n_max):
    """Negative control: {γ_sigma_z, D_truthful} ≠ 0.

    γ_sigma_z commutes with D_truthful (both are diagonal in chirality);
    so {γ_sigma_z, D_truthful} = 2 D_truthful ≠ 0 generically.
    """
    g = build_gamma_GV(n_max, convention="sigma_z")
    res = g.verify_anticommutes_dirac("truthful")
    # Expected: 2 * max|D_truthful| = 2 * (n_max + 0.5) (top Dirac eigenvalue magnitude).
    expected_lower_bound = 2.0 * (n_max + 0.5) - 1e-9
    assert res > expected_lower_bound, (
        f"sigma_z convention should NOT anticommute with D_truthful at "
        f"n_max={n_max}; got {res:.3e}"
    )


# ---------------------------------------------------------------------------
# (D) Anticommutation with the offdiag CH Dirac
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_sigma_x_does_not_anticommute_offdiag(n_max):
    """For the offdiag CH Dirac, the within-chirality offdiag part is
    block-diagonal symmetric, so it commutes with sigma_x ⊗ I. The
    anticommutator is therefore nonzero. Honest disclosure of the
    expected R3.5/R3.2 mismatch.

    n_max = 1 has no offdiag entries (only one shell), so we skip it.
    """
    g = build_gamma_GV(n_max, convention="sigma_x")
    res = g.verify_anticommutes_dirac("offdiag")
    assert res > 1e-3, (
        f"sigma_x convention should NOT anticommute with D_offdiag at "
        f"n_max={n_max} (within-chirality offdiag part is symmetric); "
        f"got {res:.3e}"
    )


# ---------------------------------------------------------------------------
# (E) Operator-system preservation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_sigma_x_preserves_O(n_max):
    """For sigma_x convention, γ M γ = M exactly for every generator."""
    g = build_gamma_GV(n_max, convention="sigma_x")
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    ok, worst, failures = g.verify_preserves_operator_system(op_sys, tol=1e-10)
    assert ok, (
        f"sigma_x must preserve O at n_max={n_max}; failures: {failures}"
    )
    # Stronger: γ M γ = M exactly, so worst residual should be 0.
    assert worst < 1e-12, (
        f"sigma_x should give γMγ = M exactly; worst residual {worst:.3e}"
    )


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_sigma_z_preserves_O(n_max):
    """sigma_z also preserves O (same reason: γMγ = M for block-diagonal M).
    Diagnostic: the operator-system preservation does NOT distinguish
    the two conventions; only the anticommutator with D does.
    """
    g = build_gamma_GV(n_max, convention="sigma_z")
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    ok, worst, failures = g.verify_preserves_operator_system(op_sys, tol=1e-10)
    assert ok, (
        f"sigma_z preserves O at n_max={n_max}; failures: {failures}"
    )
    assert worst < 1e-12


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_sigma_x_gMg_equals_M_exactly(n_max):
    """Strong statement: γ M γ = M exactly (matrix equality, not just span)."""
    g = build_gamma_GV(n_max, convention="sigma_x")
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    for M in op_sys.multiplier_matrices:
        diff = g.matrix @ M @ g.matrix - M
        assert np.max(np.abs(diff)) < 1e-12


# ---------------------------------------------------------------------------
# (F) Commutation with multipliers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_sigma_x_commutes_with_all_multipliers(n_max):
    """[γ_sigma_x, M] = 0 exactly for every scalar generator M."""
    g = build_gamma_GV(n_max, convention="sigma_x")
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    for M in op_sys.multiplier_matrices:
        res = g.verify_commutes_with_multiplier(M)
        assert res < 1e-12


# ---------------------------------------------------------------------------
# (G) Audit driver
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_audit_returns_dict(n_max):
    """audit_chirality_grading returns the expected keys."""
    result = audit_chirality_grading(n_max, convention="sigma_x")
    expected_keys = {
        "n_max", "convention", "dim_H", "dim_O", "n_generators",
        "gamma_squared_minus_I", "hermitian_residual", "unitary_residual",
        "anticommutator_truthful", "anticommutator_offdiag",
        "preserves_O", "preserves_O_worst_residual", "preserves_O_failures",
        "max_commutator_with_M",
    }
    assert expected_keys.issubset(set(result.keys()))


def test_audit_sigma_x_truthful_zero():
    """sigma_x convention's truthful anticommutator is machine zero at all
    tested n_max."""
    for n_max in (1, 2, 3):
        result = audit_chirality_grading(n_max, convention="sigma_x")
        assert result["anticommutator_truthful"] < 1e-14


def test_audit_sigma_x_dim_matches():
    for n_max in (1, 2, 3):
        result = audit_chirality_grading(n_max, convention="sigma_x")
        assert result["dim_H"] == 2 * spinor_dim(n_max)
