"""Tests for geovac.operator_system_lorentzian.

Sprint L3a-1 (2026-05-17): the Lorentzian truncated operator system on
the BBB Krein spectral triple at signature (3, 1), built on the L2-B
Krein space and the FullDiracTruncatedOperatorSystem chirality-doubled
spinor lift.

Verifies, at n_max in {1, 2, 3} and N_t in {1, 3}:

  - Construction: instantiation, basic dimensions, multiplier count.
  - *-closure: every generator's conjugate transpose lies in O^L.
  - LOAD-BEARING Riemannian limit at N_t = 1: O^L reduces bit-exactly
    to FullDiracTruncatedOperatorSystem(n_max), with max_residual = 0.0
    in float64.
  - Propagation number: prop = 2 under the "achievable envelope"
    convention (Weyl-block-diagonal subspace of dim dim_Weyl^2 * N_t),
    matching Paper 32 §III prop=2.  prop = infinity under the "full"
    envelope convention (B(K) of dim dim_K^2), reflecting the
    structural fact that scalar multipliers can never reach chirality-
    flipping operators on FullDirac.
  - Krein-positive restriction: every multiplier preserves K^+ in this
    construction, so the restriction is trivial.
  - Witness pair: a = M^{2,1,0}^spat (x) g_0, b = a^*, a*b NOT in O^L.
  - Identity in O^L (constant function f = 1 acts as identity).
  - Zero regression in upstream Sprint L2 tests (verified externally).
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    full_dirac_dim,
)
from geovac.operator_system_lorentzian import (
    LorentzianTruncatedOperatorSystem,
    restrict_to_lorentzian_wedge,
    temporal_grid_symmetric,
    temporal_multiplier_matrices,
    witness_pair_lorentzian,
)


# ---------------------------------------------------------------------------
# Temporal multiplier basics
# ---------------------------------------------------------------------------


def test_temporal_grid_symmetric_basic():
    """Symmetric grid construction matches L2-B convention."""
    assert temporal_grid_symmetric(1)[0] == 0.0
    grid = temporal_grid_symmetric(3, T_max=1.0)
    assert grid.shape == (3,)
    assert grid[0] == -1.0
    assert grid[1] == 0.0
    assert grid[2] == +1.0


def test_temporal_grid_symmetric_validation():
    """Rejects N_t < 1 or T_max <= 0."""
    with pytest.raises(ValueError):
        temporal_grid_symmetric(0)
    with pytest.raises(ValueError):
        temporal_grid_symmetric(3, T_max=0.0)
    with pytest.raises(ValueError):
        temporal_grid_symmetric(3, T_max=-1.0)


def test_temporal_multiplier_matrices_basic():
    """Polynomial basis t^p, p = 0..N_t-1.  At N_t=1: identity only."""
    mats = temporal_multiplier_matrices(1)
    assert len(mats) == 1
    assert mats[0].shape == (1, 1)
    assert np.allclose(mats[0], np.eye(1))

    mats = temporal_multiplier_matrices(3, T_max=1.0)
    assert len(mats) == 3
    # All should be diagonal
    for m in mats:
        assert np.allclose(m - np.diag(np.diag(m)), 0)
    # g_0 = 1 (identity)
    assert np.allclose(mats[0], np.eye(3))
    # g_1 = t = [-1, 0, 1]
    assert np.allclose(np.diag(mats[1]), np.array([-1.0, 0.0, 1.0]))
    # g_2 = t^2 = [1, 0, 1]
    assert np.allclose(np.diag(mats[2]), np.array([1.0, 0.0, 1.0]))


# ---------------------------------------------------------------------------
# Construction basics
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (2, 1), (1, 3), (2, 3)])
def test_construction_basic(n_max, N_t):
    """Instantiation succeeds; basic dims and multiplier count consistent."""
    O = LorentzianTruncatedOperatorSystem(n_max=n_max, N_t=N_t)
    assert O.n_max == n_max
    assert O.N_t == N_t
    assert O.dim_K == full_dirac_dim(n_max) * N_t
    assert O.envelope_dim == O.dim_K ** 2
    # multiplier_labels = (NLM, p) tuples
    assert len(O.multiplier_labels) == len(O.multiplier_matrices)
    # At N_t = 1, multiplier count matches FullDirac (no temporal expansion).
    if N_t == 1:
        ref = FullDiracTruncatedOperatorSystem(n_max)
        assert len(O.multiplier_labels) == len(ref.multiplier_labels)
    # Achievable envelope dim
    dim_Weyl = O.dim_spatial // 2
    assert O.achievable_envelope_dim == dim_Weyl * dim_Weyl * N_t


def test_construction_validation():
    """Rejects n_max < 1 or N_t < 1 or T_max <= 0."""
    with pytest.raises(ValueError):
        LorentzianTruncatedOperatorSystem(n_max=0, N_t=1)
    with pytest.raises(ValueError):
        LorentzianTruncatedOperatorSystem(n_max=2, N_t=0)
    with pytest.raises(ValueError):
        LorentzianTruncatedOperatorSystem(n_max=2, N_t=1, T_max=0.0)


# ---------------------------------------------------------------------------
# *-closure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (2, 1), (1, 3), (2, 3)])
def test_star_closed(n_max, N_t):
    """For every M_i, M_i^dagger lies in O^L."""
    O = LorentzianTruncatedOperatorSystem(n_max=n_max, N_t=N_t)
    star_ok, failures = O.is_star_closed()
    assert star_ok, (
        f"*-closure violated at n_max={n_max}, N_t={N_t}: "
        f"{len(failures)} failures"
    )


# ---------------------------------------------------------------------------
# Identity in O^L
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(1, 1), (2, 1), (1, 3), (2, 3)])
def test_identity_in_O(n_max, N_t):
    """The identity matrix lies in O^L (constant multiplier f = 1)."""
    O = LorentzianTruncatedOperatorSystem(n_max=n_max, N_t=N_t)
    in_O, residual = O.identity_in_O()
    assert in_O, f"identity not in O at n_max={n_max}, N_t={N_t}, residual={residual:.3e}"


# ---------------------------------------------------------------------------
# LOAD-BEARING: Riemannian limit at N_t = 1 (bit-exact recovery)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_riemannian_limit_bit_exact(n_max):
    """LOAD-BEARING: at N_t = 1, O^L reduces bit-exactly to
    FullDiracTruncatedOperatorSystem(n_max).

    Verifies:
      (a) dim_K == full_dirac_dim(n_max)
      (b) multiplier count match
      (c) max residual ||M_self - M_ref||_F = 0.0 exactly in float64
      (d) dim(O) match

    Failure would mean the scalar-multiplier spinor lift is structurally
    incompatible with the L2-B Krein-space spinor basis at N_t = 1, a
    major structural finding that would break Paper 43's foundation.
    """
    O = LorentzianTruncatedOperatorSystem(n_max=n_max, N_t=1)
    ok, details = O.verify_riemannian_limit()
    assert ok, (
        f"Riemannian-limit FAIL at n_max={n_max}: details={details}"
    )
    # Bit-exact match
    assert details["max_residual"] == 0.0, (
        f"max_residual != 0.0 at n_max={n_max}: "
        f"{details['max_residual']:.3e} (load-bearing falsifier)"
    )
    assert details["dim_match"]
    assert details["multiplier_count_match"]
    assert details["dim_O_match"]
    assert details["n_unmatched"] == 0


def test_riemannian_limit_from_higher_N_t():
    """verify_riemannian_limit can be called on an N_t > 1 system; it
    builds a fresh N_t = 1 reference internally."""
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=3)
    ok, details = O.verify_riemannian_limit()
    assert ok
    assert details["N_t"] == 1  # recursively built N_t = 1 reference
    assert details["max_residual"] == 0.0


# ---------------------------------------------------------------------------
# Propagation number (under achievable envelope)
# ---------------------------------------------------------------------------


def test_propagation_number_nmax_2_Nt_1_achievable_envelope():
    """prop = 2 under achievable envelope at n_max = 2, N_t = 1.

    Matches Paper 32 §III prop:propagation_2 verbatim, lifted to
    the chirality-doubled FullDirac construction.

    Achievable envelope = dim_Weyl^2 * N_t = 8^2 * 1 = 64.
    dim sequence should be [14, 64].
    """
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=1)
    assert O.achievable_envelope_dim == 64
    prop, dims = O.compute_propagation_number(
        envelope="achievable", max_k=4
    )
    assert prop == 2, (
        f"prop != 2 under achievable envelope; got prop={prop}, dims={dims}"
    )
    assert dims == [14, 64]


def test_propagation_number_nmax_2_Nt_1_full_envelope_is_infinity():
    """prop = infinity under FULL envelope at n_max = 2, N_t = 1.

    Structural finding: scalar multipliers act on H_GV chirality-
    doubled as M (+) M, so products O^L^k can never reach the full
    M_{dim_K}(C) = M_{16}(C) envelope of dim 256.  The achievable
    subspace is dim_Weyl^2 = 64, a factor of 4 smaller than the full
    envelope.  Reported by compute_propagation_number(envelope='full')
    as prop = -1 (saturated below target).
    """
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=1)
    prop, dims = O.compute_propagation_number(
        envelope="full", max_k=4
    )
    assert prop == -1, (
        f"prop should be infinity (-1) under full envelope; got {prop}"
    )
    # Saturated dim
    assert dims[-1] < O.envelope_dim
    assert dims[-1] == dims[-2]  # saturation


def test_propagation_number_nmax_2_Nt_3_full_envelope_infinity():
    """prop = infinity at N_t > 1 (full envelope).

    Achievable envelope at N_t > 1: under the polynomial-basis
    temporal multipliers, dim(O^L^k) saturates at dim_Weyl^2 * N_t,
    short of the full achievable dim_Weyl^2 * N_t^2.  This is the
    structural finding that the COMMUTATIVE temporal subalgebra
    blocks propagation under the strict 'full' envelope.
    """
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=3)
    prop, dims = O.compute_propagation_number(
        envelope="full", max_k=3
    )
    assert prop == -1, (
        f"prop should be infinity under full envelope at N_t > 1; got {prop}"
    )
    assert dims[-1] == dims[-2]  # saturation
    assert dims[-1] < O.envelope_dim


def test_propagation_number_nmax_2_Nt_3_achievable_envelope():
    """prop = 2 at N_t = 3 (achievable envelope) -- the paper's
    prop_ach = 2 for ALL N_t >= 1 (not just N_t = 1).  The achievable
    envelope is dim_Weyl^2 * N_t and O^L^2 fills it (dim seq 42, 192)."""
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=3)
    prop, dims = O.compute_propagation_number(envelope="achievable", max_k=3)
    assert prop == 2, f"prop_ach should be 2 at N_t=3; got {prop}"
    assert dims[-1] == 192  # O^2 fills the achievable envelope dim_Weyl^2 * N_t = 64*3


def test_propagation_number_nmax_3_achievable_envelope():
    """prop = 2 at n_max = 3, N_t = 1 under achievable envelope.

    Achievable envelope = dim_Weyl^2 = 20^2 = 400.
    dim sequence: [55, 400] (matches Paper 32 §III dim(O_3) = 55,
    dim(O_3^2) = N_scalar^2 = 196 -> 196 * 2 = ... wait.

    Actually the achievable on FullDirac is the Weyl-block-diagonal
    subspace, not 2 x the scalar.  At n_max = 3, dim_Weyl = 20 (cumulative
    spinor_dim), so achievable envelope = 20^2 = 400.  The scalar
    operator system reaches 14^2 = 196.  The discrepancy is a structural
    feature of the spinor lift: the operator system on the spinor bundle
    has more degrees of freedom than the scalar one because the spinor
    basis carries (l, m_j) instead of just (l, m).

    The test verifies that the prop number stays at 2 (matching the
    Paper 32 §III result), but the dim sequence reflects the spinor-
    bundle structure.
    """
    O = LorentzianTruncatedOperatorSystem(n_max=3, N_t=1)
    assert O.achievable_envelope_dim == 400
    prop, dims = O.compute_propagation_number(
        envelope="achievable", max_k=3
    )
    assert prop == 2, (
        f"prop != 2 at n_max=3, N_t=1 under achievable envelope; "
        f"got prop={prop}, dims={dims}"
    )
    assert dims[-1] == 400


# ---------------------------------------------------------------------------
# Krein-positive restriction (verified trivial in this construction)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max,N_t", [(2, 1), (2, 3), (3, 1)])
def test_krein_positive_restriction_trivial(n_max, N_t):
    """All scalar multipliers in O^L preserve the Krein-positive cone K^+.

    Structural reason: scalar multipliers M^spat = M (+) M (chirality-
    doubled) commute with J_spatial (chirality-swap) because J swaps the
    two equal copies of M.  Therefore [J, multiplier] = 0 for every
    multiplier in O^L, and the Krein-positive restriction is trivial.

    Test verifies that ALL multipliers are in the preserver set.
    """
    O = LorentzianTruncatedOperatorSystem(n_max=n_max, N_t=N_t)
    preserving_indices, preserving_labels = O.krein_positive_preservers()
    assert len(preserving_indices) == len(O.multiplier_matrices), (
        f"Some multipliers do NOT preserve K^+ at n_max={n_max}, N_t={N_t}: "
        f"{len(preserving_indices)} / {len(O.multiplier_matrices)} preserve"
    )


def test_restrict_to_krein_positive_returns_equivalent_system():
    """Since restriction is trivial, the sub-system has the same dim."""
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=3)
    sub = O.restrict_to_krein_positive()
    assert sub.dim == O.dim
    assert len(sub.multiplier_matrices) == len(O.multiplier_matrices)


# ---------------------------------------------------------------------------
# Witness pair: a, b in O^L with ab NOT in O^L
# ---------------------------------------------------------------------------


def test_witness_pair_nmax_2_Nt_1():
    """Lorentzian lift of Paper 32 §III witness pair: a*b not in O^L.

    Take a = M^{N=2, L=1, M=0}^spat (x) g_0, b = a^*.  The product a*b
    factorizes as (M^spat M^spat *) (x) g_0^2; the spatial factor is the
    Riemannian witness (Paper 32 §III) which is NOT in O^spat.  Therefore
    a*b is NOT in O^L.
    """
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=1)
    a, b, ab, test = witness_pair_lorentzian(
        O, N_target=2, L_target=1, M_target=0, p_target=0
    )
    assert a is not None, "M^{2,1,0,0} generator not in O^L"
    in_O_L, residual = test
    assert not in_O_L, (
        f"Witness pair unexpectedly in O^L "
        f"(residual {residual:.3e})"
    )
    # Significant residual, pinned to the paper's ~38% non-closure (Thm 8.1: 0.3812)
    assert 0.30 < residual < 0.45, (
        f"ab residual {residual:.4f} should be ~0.38 (the paper's ~38% non-closure)"
    )


def test_witness_pair_components_in_O_L():
    """The witness components a and b are individually in O^L (closed
    under *) even though a*b is not (not closed under multiplication).
    """
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=1)
    a, b, ab, _ = witness_pair_lorentzian(O)
    a_in, _ = O.contains(a)
    b_in, _ = O.contains(b)
    assert a_in
    assert b_in


def test_witness_pair_nmax_2_Nt_3():
    """Same witness pair holds at N_t = 3 with p = 0 (constant temporal)."""
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=3)
    a, b, ab, test = witness_pair_lorentzian(
        O, N_target=2, L_target=1, M_target=0, p_target=0
    )
    assert a is not None
    in_O_L, residual = test
    assert not in_O_L
    assert 0.30 < residual < 0.45  # paper Thm 8.1: ~0.3589 at N_t=3


# ---------------------------------------------------------------------------
# Connection to L2-E hemispheric wedge
# ---------------------------------------------------------------------------


def test_restrict_to_lorentzian_wedge_basic():
    """Wedge-restriction connects O^L to the Sprint L2-E hemispheric wedge.

    Verifies basic well-formedness of the wedge restriction:
      - wedge dimension > 0
      - each wedge-block matrix has shape (dim_W_L, dim_W_L)
      - wedge linear-span dim is bounded by both dim(O) and dim_W_L^2.
    """
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=1)
    info = restrict_to_lorentzian_wedge(O)
    dim_W_L = info["dim_W_L"]
    assert dim_W_L > 0
    assert dim_W_L <= O.dim_K
    for M in info["wedge_block_matrices"]:
        assert M.shape == (dim_W_L, dim_W_L)
    span_dim = info["wedge_linear_span_dim"]
    assert span_dim <= O.dim
    assert span_dim <= dim_W_L * dim_W_L


@pytest.mark.slow
def test_restrict_to_lorentzian_wedge_Nt_3():
    """Wedge restriction at N_t > 1 builds the joint spatial * temporal wedge."""
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=3)
    info = restrict_to_lorentzian_wedge(O)
    dim_W_L = info["dim_W_L"]
    assert dim_W_L > 0
    # At N_t = 3 symmetric grid, P_t_positive selects t >= 0 indices: t = 0, +1
    # so 2 of 3 grid points, dim_W_L = 8 (spatial wedge) * 2 = 16
    # Check: dim_W_L = (dim_spatial / 2) * 2 = 16
    # (spatial wedge dim = 8 at n_max=2; positive-t dim = 2 at N_t=3)
    assert dim_W_L == 16


# ---------------------------------------------------------------------------
# Identity is the (1, 0, 0, 0) multiplier
# ---------------------------------------------------------------------------


def test_identity_multiplier_label():
    """The (N=1, L=0, M=0, p=0) generator is proportional to the identity."""
    O = LorentzianTruncatedOperatorSystem(n_max=2, N_t=1)
    for label, mat in O.basis_matrices:
        if label == (1, 0, 0, 0):
            # mat should be a scalar multiple of identity
            assert np.allclose(mat - np.diag(np.diag(mat)), 0)
            diag_vals = np.diag(mat)
            assert np.allclose(diag_vals, diag_vals[0])
            assert abs(diag_vals[0]) > 1e-10
            return
    pytest.fail("(1, 0, 0, 0) generator not found in O^L")
