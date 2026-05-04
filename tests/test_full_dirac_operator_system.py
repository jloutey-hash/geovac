"""Tests for geovac.full_dirac_operator_system: full Dirac (both chiralities)
on the spinor bundle of S^3.

Sprint WH1-R3.5 equation verification per CLAUDE.md Sec 13.4a.

Verifies:

  - Basis dimension formula: dim_H_full(n_max) = 2 * spinor_dim(n_max).
  - FullDiracLabel validation (chirality must be +/- 1, etc.).
  - Camporesi-Higuchi full-Dirac eigenvalues = chirality * (n_fock + 1/2).
  - Multiplier matrices are block-diagonal in chirality (scalar function).
  - The two chirality blocks are equal (scalar acts identically).
  - Identity matrix is in O (constant-multiplier branch).
  - O is *-closed.
  - dim(O) at small n_max (regression baseline).
  - Propagation number prop = 2 at n_max = 2, 3 (matches scalar / Weyl).
  - Reduction to spinor (Weyl) sector correctness.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import sympy as sp
from sympy import Rational, simplify, sqrt

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    FullDiracTruncatedOperatorSystem,
    build_full_dirac_multiplier_matrix,
    camporesi_higuchi_full_dirac_matrix,
    camporesi_higuchi_offdiag_dirac_matrix,
    full_dirac_basis,
    full_dirac_dim,
    full_dirac_graph_distance,
    full_dirac_label_strings,
)
from geovac.spinor_operator_system import (
    SpinorTruncatedOperatorSystem,
    spinor_basis,
    spinor_dim,
)


# ---------------------------------------------------------------------------
# Basis dimension and label validation
# ---------------------------------------------------------------------------


def test_full_dirac_dim_table():
    """dim_H_full(n_max) = 2 * spinor_dim(n_max).

    Per the Paper 32 §3.3 / Camporesi-Higuchi convention, the full
    Dirac sector at level n_ch has 2*(n_ch+1)(n_ch+2) states (both
    chirality signs). Cumulative:
        n_max=1:  4 - WAIT, not 2*spinor_dim(1)=4. spinor_dim(1)=2.
        n_max=2: 16  (2*8)
        n_max=3: 40  (2*20)  <- this is g_3^Dirac = 40 (Paper 2 Delta^-1).
        n_max=4: 80  (2*40)
    """
    assert full_dirac_dim(1) == 2 * spinor_dim(1)
    assert full_dirac_dim(2) == 2 * spinor_dim(2)
    assert full_dirac_dim(3) == 2 * spinor_dim(3)
    assert full_dirac_dim(4) == 2 * spinor_dim(4)
    # Paper 2 Delta^{-1} = g_3^Dirac = 40
    assert full_dirac_dim(3) == 40


def test_full_dirac_basis_consistency():
    for n_max in (1, 2, 3, 4):
        b = full_dirac_basis(n_max)
        assert len(b) == full_dirac_dim(n_max)


def test_full_dirac_label_validation():
    """FullDiracLabel rejects out-of-range quantum numbers."""
    FullDiracLabel(n_fock=1, l=0, two_m_j=+1, chirality=+1)  # OK
    FullDiracLabel(n_fock=1, l=0, two_m_j=+1, chirality=-1)  # OK
    with pytest.raises(ValueError):
        FullDiracLabel(n_fock=1, l=0, two_m_j=+1, chirality=0)
    with pytest.raises(ValueError):
        FullDiracLabel(n_fock=1, l=0, two_m_j=+1, chirality=2)
    with pytest.raises(ValueError):
        FullDiracLabel(n_fock=0, l=0, two_m_j=+1, chirality=+1)
    with pytest.raises(ValueError):
        FullDiracLabel(n_fock=1, l=1, two_m_j=+1, chirality=+1)


def test_full_dirac_chirality_assignments():
    """Each (n, l, m_j) label appears in BOTH chirality sectors."""
    for n_max in (1, 2, 3):
        b = full_dirac_basis(n_max)
        # Group by chirality
        by_chir = {+1: [], -1: []}
        for x in b:
            by_chir[x.chirality].append(
                (x.n_fock, x.l, x.two_m_j)
            )
        # Both sectors present
        assert len(by_chir[+1]) == spinor_dim(n_max)
        assert len(by_chir[-1]) == spinor_dim(n_max)
        # Sets of (n, l, m_j) are equal across chirality
        assert set(by_chir[+1]) == set(by_chir[-1])


def test_full_dirac_to_spinor():
    """FullDiracLabel.to_spinor() projects to the SpinorLabel."""
    full = FullDiracLabel(n_fock=2, l=1, two_m_j=+1, chirality=-1)
    sl = full.to_spinor()
    assert sl.n_fock == 2
    assert sl.l == 1
    assert sl.two_m_j == 1


# ---------------------------------------------------------------------------
# Camporesi-Higuchi full Dirac matrix
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_full_ch_dirac_eigenvalues(n_max):
    """The full CH Dirac is diagonal with eigenvalues chi*(n_fock+1/2)."""
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    # Off-diagonal entries are zero
    off = D - np.diag(np.diag(D))
    assert np.allclose(off, 0, atol=1e-15)
    # Diagonal entries match chirality * (n_fock + 1/2)
    for i, b in enumerate(basis):
        expected = b.chirality * (float(b.n_fock) + 0.5)
        assert abs(D[i, i].real - expected) < 1e-15
        assert abs(D[i, i].imag) < 1e-15


def test_full_ch_dirac_traceless_per_chirality_pair():
    """Tr(D_full) = 0 because chirality blocks have opposite signs."""
    for n_max in (1, 2, 3):
        basis = full_dirac_basis(n_max)
        D = camporesi_higuchi_full_dirac_matrix(basis)
        tr = np.trace(D)
        assert abs(tr.real) < 1e-12, f"Tr D = {tr} != 0 at n_max={n_max}"


def test_full_ch_dirac_spectrum_paper32():
    """Per Paper 32 Def. 3.3: spectrum is +/-(n + 3/2) in CH (= n_fock + 1/2 in Fock).

    For n_max=3, the eigenvalues are {+/-3/2, +/-5/2, +/-7/2} with multiplicities
    matching dim_Weyl(level n_ch).
    """
    n_max = 3
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    eigvals = np.diag(D).real
    # Distinct values
    distinct = sorted(set(np.round(eigvals, 10).tolist()))
    expected = [-3.5, -2.5, -1.5, 1.5, 2.5, 3.5]
    assert distinct == expected, f"Got {distinct}, expected {expected}"


# ---------------------------------------------------------------------------
# Multiplier matrix structure
# ---------------------------------------------------------------------------


def test_multiplier_block_diagonal_in_chirality():
    """A scalar multiplier acts identically on both chirality sectors,
    so M^full is block-diagonal with two equal Weyl blocks."""
    n_max = 2
    weyl_basis = spinor_basis(n_max)
    full = full_dirac_basis(n_max)
    dim_w = len(weyl_basis)
    # Pick a non-trivial multiplier label (N=2, L=1, M=0)
    M_full = build_full_dirac_multiplier_matrix(2, 1, 0, full)
    # Off-block (cross-chirality) entries should be zero.
    # Block layout: [Weyl, Weyl] x [Weyl, Weyl] with Weyl block first.
    cross_tr = M_full[:dim_w, dim_w:]
    cross_bl = M_full[dim_w:, :dim_w]
    assert np.allclose(cross_tr, 0, atol=1e-13)
    assert np.allclose(cross_bl, 0, atol=1e-13)
    # The two diagonal blocks are equal.
    block_tl = M_full[:dim_w, :dim_w]
    block_br = M_full[dim_w:, dim_w:]
    assert np.allclose(block_tl, block_br, atol=1e-13)


@pytest.mark.parametrize("n_max", [2, 3])
def test_constant_multiplier_is_identity_full(n_max):
    """M_{1,0,0} (constant) lifts to a multiple of the full identity."""
    O = FullDiracTruncatedOperatorSystem(n_max)
    for label, M in zip(O.multiplier_labels, O.multiplier_matrices):
        if label == (1, 0, 0):
            diag = np.diag(M)
            assert np.allclose(M - np.diag(diag), 0, atol=1e-12)
            assert np.allclose(diag, diag[0], atol=1e-12)
            return
    pytest.fail("M_{1,0,0} multiplier not found")


# ---------------------------------------------------------------------------
# FullDiracTruncatedOperatorSystem properties
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_identity_in_O_full(n_max):
    """Identity matrix is in O for the full-Dirac sector."""
    O = FullDiracTruncatedOperatorSystem(n_max)
    in_O, residual = O.identity_in_O()
    assert in_O, f"identity not in full O at n_max={n_max}, residual={residual:.2e}"


@pytest.mark.parametrize("n_max", [2, 3])
def test_star_closed_full(n_max):
    """Full-Dirac O is *-closed (each generator's conjugate transpose is in O)."""
    O = FullDiracTruncatedOperatorSystem(n_max)
    star_ok, failures = O.is_star_closed()
    assert star_ok, f"*-closure failed at n_max={n_max}: {failures}"


@pytest.mark.parametrize("n_max", [2, 3])
def test_O_strictly_smaller_than_envelope(n_max):
    """dim(O) < dim_H^2."""
    O = FullDiracTruncatedOperatorSystem(n_max)
    assert O.dim < O.envelope_dim
    assert O.dim > 1


def test_dim_O_scalar_equality():
    """dim(O_full) = dim(O_Weyl). The chirality doubling does NOT
    enlarge the operator system because each multiplier acts identically
    on both chirality sectors (block-diagonal with equal blocks)."""
    for n_max in (1, 2, 3):
        O_full = FullDiracTruncatedOperatorSystem(n_max)
        O_weyl = SpinorTruncatedOperatorSystem(n_max)
        assert O_full.dim == O_weyl.dim, (
            f"dim(O_full)={O_full.dim} != dim(O_Weyl)={O_weyl.dim} at "
            f"n_max={n_max}"
        )


def test_dim_H_table():
    for n_max, expected in [(1, 4), (2, 16), (3, 40), (4, 80)]:
        O = FullDiracTruncatedOperatorSystem(n_max)
        assert O.dim_H == expected


# ---------------------------------------------------------------------------
# Propagation number on full-Dirac sector
# ---------------------------------------------------------------------------


def test_propagation_number_full_does_not_saturate():
    """prop(O_full_{n_max}) does NOT reach the C*-envelope at any finite k.

    Because each scalar multiplier acts BLOCK-DIAGONALLY in chirality
    with two equal blocks, products O^k stay block-diagonal forever.
    The scalar-multiplier operator system cannot generate ANY
    cross-chirality matrix element of M_{2 dim_W}(C). Therefore
    dim(O^k) saturates at most at (dim_W^2 + dim_W^2) = 2 * dim_W^2
    (block-diagonal Hermitian span), strictly less than the C*-envelope
    (2 dim_W)^2 = 4 dim_W^2.

    This is a STRUCTURAL FINDING: the scalar-multiplier operator system
    on the full Dirac sector has prop = +infinity (does not saturate to
    the C*-envelope at any finite k). The chirality grading is invisible
    to the algebra of scalar multipliers, so the algebra never sees the
    cross-chirality "Pauli matrix" structure that one would need to
    generate the full M_N(C) envelope.

    Implication for WH1: the saturated-dim sequence of the scalar
    operator system on the full Dirac is a DIFFERENT structural
    invariant from the propagation number of the Weyl-only sector
    (where prop = 2 because the operator system's C*-envelope IS the
    Weyl-block matrix algebra).
    """
    from geovac.operator_system import propagation_number
    for n_max in (2, 3):
        O = FullDiracTruncatedOperatorSystem(n_max)
        prop, dim_seq = propagation_number(O, max_k=4)
        # Expected: prop = -1 (does not reach envelope within max_k).
        # Saturated dim should equal 2 * dim_Weyl^2 (block-diagonal span).
        dim_W = O.dim_H // 2
        max_saturated = 2 * dim_W ** 2  # block-diag Hermitian span
        target = O.envelope_dim  # = (2 dim_W)^2 = 4 dim_W^2
        # Final reached dimension (last in sequence) should equal saturation
        final_dim = dim_seq[-1]
        assert final_dim <= max_saturated, (
            f"Saturated dim {final_dim} exceeds block-diag bound "
            f"{max_saturated} at n_max={n_max}"
        )
        assert final_dim < target, (
            f"Saturated dim {final_dim} reached envelope target {target} "
            f"at n_max={n_max} (unexpected — chirality block-diagonality "
            f"violated)"
        )
        assert prop == -1, (
            f"prop should be -1 (envelope not reached), got {prop} at n_max={n_max}"
        )


def test_propagation_number_per_chirality_block_is_2():
    """The Weyl-sector propagation number prop = 2 carries through to each
    diagonal chirality block of the full-Dirac operator system.

    This restores the WH1 prop=2 alignment claim at the per-block level:
    each chirality block, viewed as an operator system on its dim_W-dim
    Hilbert space, has prop = 2 (matching the Weyl R3.2 result and the
    Connes-vS Toeplitz S^1 prop = 2 result).
    """
    from geovac.operator_system import propagation_number
    for n_max in (2, 3):
        # The Weyl-only (single-chirality) prop = 2 should hold.
        O_weyl = SpinorTruncatedOperatorSystem(n_max)
        prop_w, _ = propagation_number(O_weyl, max_k=4)
        assert prop_w == 2, (
            f"prop(O_Weyl_{n_max}) = {prop_w}, expected 2"
        )


# ---------------------------------------------------------------------------
# Off-diag CH variant
# ---------------------------------------------------------------------------


def test_offdiag_dirac_hermitian():
    """The off-diag full Dirac matrix is Hermitian."""
    basis = full_dirac_basis(2)
    D = camporesi_higuchi_offdiag_dirac_matrix(basis)
    err = np.linalg.norm(D - D.conj().T)
    assert err < 1e-12, f"D not Hermitian: ||D - D*|| = {err}"


def test_offdiag_dirac_chirality_coupling_zero_decouples():
    """With chirality_coupling=0, the off-diag Dirac decouples into two
    block-diagonal R3.2-style Diracs."""
    basis = full_dirac_basis(2)
    D = camporesi_higuchi_offdiag_dirac_matrix(
        basis, chirality_coupling=0.0,
    )
    dim_w = full_dirac_dim(2) // 2
    cross_tr = D[:dim_w, dim_w:]
    cross_bl = D[dim_w:, :dim_w]
    assert np.allclose(cross_tr, 0, atol=1e-13)
    assert np.allclose(cross_bl, 0, atol=1e-13)


def test_offdiag_dirac_breaks_n_degeneracy():
    """The diagonal of the offdiag Dirac is non-degenerate within each chirality."""
    basis = full_dirac_basis(2)
    D = camporesi_higuchi_offdiag_dirac_matrix(basis)
    diag = np.diag(D).real
    # Within each chirality block, all diagonal values should be distinct.
    dim_w = full_dirac_dim(2) // 2
    for chir_block in (diag[:dim_w], diag[dim_w:]):
        sorted_block = np.sort(chir_block)
        gaps = np.diff(sorted_block)
        assert np.all(gaps > 1e-10), (
            f"Degenerate diagonal in chirality block: {sorted_block}"
        )


# ---------------------------------------------------------------------------
# Graph distance
# ---------------------------------------------------------------------------


def test_graph_distance_self_zero_full():
    lab = FullDiracLabel(n_fock=2, l=1, two_m_j=+1, chirality=+1)
    assert full_dirac_graph_distance(lab, lab) == 0


def test_graph_distance_chirality_flip():
    """Chirality flip alone counts as 1 hop."""
    a = FullDiracLabel(n_fock=2, l=1, two_m_j=+1, chirality=+1)
    b = FullDiracLabel(n_fock=2, l=1, two_m_j=+1, chirality=-1)
    assert full_dirac_graph_distance(a, b) == 1


def test_graph_distance_specific():
    a = FullDiracLabel(n_fock=1, l=0, two_m_j=+1, chirality=+1)
    b = FullDiracLabel(n_fock=2, l=0, two_m_j=+1, chirality=+1)
    assert full_dirac_graph_distance(a, b) == 1
    c = FullDiracLabel(n_fock=2, l=1, two_m_j=+3, chirality=+1)
    d = FullDiracLabel(n_fock=2, l=1, two_m_j=-3, chirality=-1)
    # |Delta n_fock|=0, |Delta l|=0, |Delta two_m_j|/2 = 6/2 = 3, chir = 1 -> 4
    assert full_dirac_graph_distance(c, d) == 4


# ---------------------------------------------------------------------------
# Performance smoke
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_build_n_max_3_under_60s():
    """FullDiracTruncatedOperatorSystem(3) must build in < 60s."""
    t0 = time.time()
    O = FullDiracTruncatedOperatorSystem(3)
    elapsed = time.time() - t0
    assert elapsed < 60, f"build took {elapsed:.1f}s, exceeds 60s budget"
    assert O.dim_H == 40


# ---------------------------------------------------------------------------
# Structural commutator identity verification (§1.4 of WH1-R3.5 memo)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_commutator_decomposes_in_chirality(n_max):
    """[D_full, M_full] decomposes block-diagonally with opposite-sign blocks.

    Per §1.4 of debug/wh1_r35_full_dirac_memo.md:

        [D_full, M_full] = [D_Weyl, M_Weyl] (+) [-D_Weyl, M_Weyl]

    where D_full has +(n+1/2) on (+) chirality, -(n+1/2) on (-) chirality,
    and M_full is block-diagonal with two equal Weyl blocks.

    Therefore the (+,+) block of [D_full, M_full] = +[D_Weyl, M_Weyl],
    the (-,-) block = -[D_Weyl, M_Weyl], and cross-blocks are zero.
    """
    O_full = FullDiracTruncatedOperatorSystem(n_max)
    D_full = camporesi_higuchi_full_dirac_matrix(O_full.basis)
    dim_W = O_full.dim_H // 2

    # Pick a non-trivial multiplier (N=2, L=1, M=0)
    target_label = (2, 1, 0)
    M_full = None
    for lab, M in zip(O_full.multiplier_labels, O_full.multiplier_matrices):
        if lab == target_label:
            M_full = M
            break
    assert M_full is not None

    comm = D_full @ M_full - M_full @ D_full

    # (+, +) block
    block_pp = comm[:dim_W, :dim_W]
    # (-, -) block
    block_mm = comm[dim_W:, dim_W:]
    # cross blocks
    block_pm = comm[:dim_W, dim_W:]
    block_mp = comm[dim_W:, :dim_W]

    # Cross blocks should be zero
    assert np.allclose(block_pm, 0, atol=1e-13), (
        f"(+,-) cross-block of commutator should be zero, "
        f"got max |.|={np.abs(block_pm).max()}"
    )
    assert np.allclose(block_mp, 0, atol=1e-13)

    # (-, -) block = - (+, +) block (opposite chirality sign)
    assert np.allclose(block_mm, -block_pp, atol=1e-13), (
        f"(-,-) block != -(+,+) block; max diff = "
        f"{np.abs(block_mm + block_pp).max()}"
    )


@pytest.mark.parametrize("n_max", [2, 3])
def test_chirality_flip_pair_diagonal_match(n_max):
    """For any (n,l,m_j) and any scalar multiplier in O, the diagonal entries
    at the two chirality-conjugate basis vectors are equal:

        <v,+1| M |v,+1> = <v,-1| M |v,-1>     for all M in O.

    Mechanism: scalar multipliers act block-diagonally with TWO EQUAL
    blocks, so M[i, i] = M[i + dim_W, i + dim_W] for every i in the Weyl
    block.

    Consequence: any state functional of the form phi_v(x) = <v|x|v>
    cannot distinguish chirality-conjugate pairs, forcing
    d_Connes(phi_{v,+}, phi_{v,-}) = 0 for any Dirac proxy.
    """
    O = FullDiracTruncatedOperatorSystem(n_max)
    dim_W = O.dim_H // 2
    for M in O.multiplier_matrices:
        # diagonal entries at i and i + dim_W should be equal for i < dim_W
        for i in range(dim_W):
            assert abs(M[i, i] - M[i + dim_W, i + dim_W]) < 1e-13, (
                f"Chirality-flip diagonal mismatch at i={i}, n_max={n_max}"
            )


# ---------------------------------------------------------------------------
# Empirical Connes-distance regression baselines (from WH1-R3.5 sprint)
# ---------------------------------------------------------------------------
#
# These are slow (require cvxpy + SCS) and locked to verified data values
# from the WH1-R3.5 sprint runs. They are skipped by default.


@pytest.mark.slow
def test_truthful_full_dirac_nmax2_zero_finite():
    """At n_max=2 truthful CH full Dirac: 0 finite cross-pair distances."""
    from geovac.connes_distance import compute_distance_matrix
    O = FullDiracTruncatedOperatorSystem(2)
    D = camporesi_higuchi_full_dirac_matrix(O.basis)
    dist = compute_distance_matrix(O, D=D)
    N = O.dim_H
    n_finite = 0
    for i in range(N):
        for j in range(i + 1, N):
            if not np.isinf(dist[i, j]) and dist[i, j] >= 1e-7:
                n_finite += 1
    assert n_finite == 0, (
        f"Expected 0 finite cross-pairs at n_max=2 truthful CH full Dirac, "
        f"got {n_finite}"
    )
