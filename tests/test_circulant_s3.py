"""Tests for geovac.circulant_s3: WH1 round-3 falsification comparator.

Verifies that the commutative C*-subalgebra A_circ_N of diagonal matrices has
propagation number 1 at every tested N, contrasting with prop(O_{n_max}) = 2
from WH1 round 2.

This is a *plausibility / falsification* test: confirming that prop = 1 here
demonstrates that the round-2 prop = 2 result for the Connes-vS truncated
operator system on Fock-projected S^3 is structurally specific to the
spectral-truncation construction, not a generic property of any finite
truncation of S^3.

See `debug/wh1_round3_falsification_memo.md` for the full discussion.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.circulant_s3 import (
    CirculantS3Truncation,
    circulant_for_geovac,
    compare_to_geovac,
)


# ---------------------------------------------------------------------------
# Basic algebraic properties
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_points", [1, 2, 5, 14, 30])
def test_dim_equals_n_points(n_points):
    """A_circ_N has complex dimension N: the N elementary diagonals are linearly
    independent."""
    A = CirculantS3Truncation(n_points=n_points)
    assert A.dim == n_points
    assert A.envelope_dim == n_points
    assert A.ambient_envelope_dim == n_points * n_points


@pytest.mark.parametrize("n_points", [1, 2, 5, 14, 30])
def test_identity_in_algebra(n_points):
    """The identity 1 = sum_i E_{i, i} lies in A_circ_N at machine precision."""
    A = CirculantS3Truncation(n_points=n_points)
    in_algebra, residual = A.identity_in_algebra()
    assert in_algebra, f"identity not in A_circ_N at N={n_points}, residual={residual}"
    assert residual < 1e-12


@pytest.mark.parametrize("n_points", [1, 2, 5, 14])
def test_multiplicative_closure(n_points):
    """A_circ_N is closed under matrix multiplication (E_i E_j = delta_{ij} E_i)."""
    A = CirculantS3Truncation(n_points=n_points)
    closed, failures = A.verify_multiplicative_closure()
    assert closed, f"A_circ_N not multiplicatively closed at N={n_points}: failures={failures}"


@pytest.mark.parametrize("n_points", [1, 2, 5, 14, 30])
def test_star_closure(n_points):
    """A_circ_N is closed under adjoints (each E_i is real-diagonal so E_i^* = E_i)."""
    A = CirculantS3Truncation(n_points=n_points)
    closed, failures = A.verify_star_closure()
    assert closed, f"A_circ_N not *-closed at N={n_points}: failures={failures}"


# ---------------------------------------------------------------------------
# Propagation number — the headline test
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_points", [1, 2, 5, 14, 30])
def test_propagation_intrinsic_envelope_is_one(n_points):
    """prop(A_circ_N) = 1 with the intrinsic envelope (= A_circ_N itself).

    This is the apples-to-apples Connes-vS Definition 2.39 reading: A_circ_N is
    its own C*-envelope as an abstract C*-algebra, so (A_circ_N)^1 = A_circ_N
    is trivially the C*-algebra, and prop = 1.
    """
    A = CirculantS3Truncation(n_points=n_points)
    result = A.compute_propagation_number(envelope="intrinsic")
    assert result.prop == 1, (
        f"Expected prop = 1 for the commutative C*-algebra A_circ_N (intrinsic "
        f"envelope reading) at N={n_points}, got prop = {result.prop}, "
        f"dim_sequence = {result.dim_sequence}"
    )
    # Also verify dim(A^1) = N (= envelope_dim).
    assert result.dim_sequence[0] == n_points


@pytest.mark.parametrize("n_points", [2, 5, 14])
def test_propagation_ambient_envelope_is_infinite(n_points):
    """In the ambient-envelope reading (target = M_N(C), N^2), A_circ_N is
    commutative and never fills the non-commutative M_N(C) for N > 1, so prop
    is infinity (returned as -1).

    This is consistent: the commutative algebra cannot generate M_N(C) by
    *any* number of products of its own elements.  This is the "structural
    obstruction" reading that complements the round-2 finding.
    """
    A = CirculantS3Truncation(n_points=n_points)
    result = A.compute_propagation_number(envelope="ambient", max_k=4)
    assert result.prop == -1, (
        f"Expected prop = infinity (-1) for the commutative A_circ_N -> M_N(C) "
        f"ambient envelope at N={n_points}, got prop = {result.prop}, "
        f"dim_sequence = {result.dim_sequence}"
    )
    # Sanity: dim(A^k) saturates at N for all k.
    for k, d in enumerate(result.dim_sequence, start=1):
        assert d == n_points, (
            f"Expected dim(A^{k}) = N = {n_points}, got {d}; A_circ_N is not "
            f"saturating at N (multiplicative-closure check probably failed)."
        )


def test_n_points_one_trivial_case():
    """N = 1 is the trivial case: A_circ_1 = M_1(C) = C, both intrinsic and
    ambient envelopes coincide, prop = 1."""
    A = CirculantS3Truncation(n_points=1)
    intrinsic = A.compute_propagation_number(envelope="intrinsic")
    ambient = A.compute_propagation_number(envelope="ambient")
    assert intrinsic.prop == 1
    assert ambient.prop == 1
    assert intrinsic.dim_sequence[0] == 1
    assert ambient.dim_sequence[0] == 1


# ---------------------------------------------------------------------------
# Validation: prop(GeoVac) is NOT 1 (cross-check round 2 from the other side)
# ---------------------------------------------------------------------------


def test_geovac_propagation_remains_two_at_nmax2():
    """Reproduce round 2: prop(GeoVac O_{n_max=2}) = 2.

    This is the apples-to-apples comparison anchor for the round-3 falsification
    table.
    """
    result = compare_to_geovac(n_max=2, match="envelope")
    assert result["geovac_prop"] == 2, (
        f"Round 2 reproduction failed: expected prop(GeoVac O_{{n_max=2}}) = 2, "
        f"got {result['geovac_prop']}"
    )


# ---------------------------------------------------------------------------
# Cross-comparison: structural difference confirmed
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_structural_difference_envelope_matched(n_max):
    """At dim-matched N (= GeoVac envelope dim), GeoVac has prop = 2 while
    circulant has intrinsic prop = 1.  This is the headline round-3 result.
    """
    result = compare_to_geovac(n_max=n_max, match="envelope")
    # GeoVac side
    assert result["geovac_prop"] == 2, (
        f"GeoVac prop != 2 at n_max={n_max}; round-2 reproduction failed."
    )
    # Circulant side (intrinsic envelope)
    assert result["circulant_intrinsic_prop"] == 1, (
        f"Circulant intrinsic prop != 1 at n_max={n_max}; "
        f"the commutative-C*-algebra triviality is broken."
    )
    # Conclusive
    assert result["structural_difference"] is True, (
        "structural_difference should be True (GeoVac prop = 2, circulant prop = 1)"
    )
    # Dimension match
    assert result["circulant_N"] == result["geovac_N"]


def test_structural_difference_operator_system_matched():
    """At dim-matched dim(O), GeoVac O has prop = 2 while a circulant of the
    same intrinsic *-algebra dimension has prop = 1.
    """
    result = compare_to_geovac(n_max=2, match="operator_system")
    assert result["geovac_prop"] == 2
    assert result["circulant_intrinsic_prop"] == 1
    # Dimension matched on the operator-system side
    assert result["circulant_N"] == result["geovac_dim_O"]
    assert result["structural_difference"] is True


# ---------------------------------------------------------------------------
# Algorithmic robustness
# ---------------------------------------------------------------------------


def test_propagation_returns_at_k1_for_circulant():
    """The propagation algorithm should converge at k = 1 for a C*-algebra
    (intrinsic envelope), not get stuck or saturate at k > 1.
    """
    A = CirculantS3Truncation(n_points=14)
    result = A.compute_propagation_number(envelope="intrinsic")
    assert result.prop == 1
    assert len(result.dim_sequence) == 1, (
        f"Expected dim_sequence to terminate at k = 1 for C*-algebra, "
        f"got dim_sequence = {result.dim_sequence}"
    )


def test_circulant_for_geovac_envelope_match():
    """circulant_for_geovac(n_max, match='envelope') gives N = N(GeoVac)."""
    from geovac.operator_system import hilbert_dim
    for n_max in (1, 2, 3, 4):
        circ = circulant_for_geovac(n_max, match="envelope")
        assert circ.n_points == hilbert_dim(n_max)


def test_circulant_for_geovac_operator_system_match():
    """circulant_for_geovac(n_max, match='operator_system') gives N = dim(O_GV)."""
    from geovac.operator_system import TruncatedOperatorSystem
    for n_max in (1, 2):  # keep small for test speed
        gv = TruncatedOperatorSystem(n_max)
        circ = circulant_for_geovac(n_max, match="operator_system")
        assert circ.n_points == gv.dim
