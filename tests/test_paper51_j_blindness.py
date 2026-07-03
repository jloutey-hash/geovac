"""Verification of Paper 51 Theorem thm:j_blindness (J-blindness of S^(2)).

Paper 51 §6.x (sec:fierz_pauli) states: the Gaussian spectral-action second
variation S^(2) restricted to the (1, 1) subspace of Herm(H_{n_max}) has
**identical eigenvalue spectra** across J = 0, 1, 2 (Schur's lemma from the
SO(4) -> diagonal SU(2) invariance of D_0). Multiplicities follow

    within-sector modes (||[D_0, V]||^2 = 0): 1 : 1 : 3
    cross-sector modes (||[D_0, V]||^2 > 0):  1 : 3 : 5

Numerical check from Paper 51 §6:
    n_max = 2: analytical S^(2) weight vs 5-point stencil, max rel err 1.24e-6.

This test reuses the existing decomposition machinery in
tests/gravity_support/g6_fierz_pauli.py (build_fp_basis, s2_weight_matrix,
CHBasis, compute_s2_stencil) — migrated from debug/ on 2026-07-03 per the
durability policy (debug/ is prune-by-design; see
tests/gravity_support/README.md) — and verifies the load-bearing claims
directly:

    (A) Per-sector eigenvalue SET equality across J = 0, 1, 2.
    (B) Within-sector multiplicity ratio 1 : 1 : 3.
    (C) Cross-sector multiplicity ratio 1 : 3 : 5.
    (D) Analytical S^(2) weight vs 5-point stencil sanity check.

Per CLAUDE.md §13.4a, this is a symbolic-identity / numerical cross-check
verification of the theorem statement.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Import the decomposition machinery from its permanent home
# tests/gravity_support/ (the wh7_support import pattern).
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent / "gravity_support"))

import g6_fierz_pauli as g6  # noqa: E402

CHBasis = g6.CHBasis
build_fp_basis = g6.build_fp_basis
s2_weight_matrix = g6.s2_weight_matrix
compute_s2_stencil = g6.compute_s2_stencil


# ---------------------------------------------------------------------------
# Helpers.
# ---------------------------------------------------------------------------

def _s2_eigenvalues(fp_basis_J, W):
    """Diagonalise S^(2) restricted to a J-sector and return sorted eigenvalues."""
    if not fp_basis_J:
        return np.array([])
    V_stack = np.array(fp_basis_J)
    K_s2 = np.einsum('aij,ij,bij->ab', V_stack, W, V_stack)
    eigs = np.linalg.eigvalsh(K_s2)
    return np.sort(eigs)[::-1]


def _cluster(values, tol_abs=1e-6, tol_rel=1e-3):
    """Group a sorted-descending numeric array into (value, multiplicity)."""
    out = []
    if len(values) == 0:
        return out
    i = 0
    arr = np.asarray(values)
    while i < len(arr):
        v = arr[i]
        j = i + 1
        while j < len(arr) and abs(arr[j] - v) < max(tol_abs, tol_rel * abs(v)):
            j += 1
        out.append((float(v), j - i))
        i = j
    return out


def _split_within_cross(fp_basis_J, D0, tol_lap=0.5):
    """Partition the FP basis vectors into within-sector and cross-sector."""
    diff = D0[:, None] - D0[None, :]
    within = []
    cross = []
    for V in fp_basis_J:
        lap_val = float(np.sum(diff ** 2 * V ** 2))
        if lap_val < tol_lap:
            within.append(V)
        else:
            cross.append(V)
    return within, cross


def _eigvals_set(eigs, tol_abs=1e-6, tol_rel=1e-3):
    """Return the SET of distinct eigenvalues (as sorted tuple) at given tolerance."""
    clusters = _cluster(np.sort(eigs)[::-1], tol_abs=tol_abs, tol_rel=tol_rel)
    return tuple(round(v, 6) for v, _ in clusters)


# ---------------------------------------------------------------------------
# Fixtures.
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def n_max_1_data():
    """Build (1,1) decomposition at n_max = 1 (fast; dim_H = 16)."""
    basis = CHBasis(n_max=1)
    fp = build_fp_basis(basis)
    Lambda_sq = 6.0
    W = s2_weight_matrix(basis.D0, Lambda_sq)
    return {"basis": basis, "fp": fp, "W": W, "Lambda_sq": Lambda_sq}


# ---------------------------------------------------------------------------
# Tests.
# ---------------------------------------------------------------------------

def test_11_subspace_dimensions_n_max_1(n_max_1_data):
    """Sanity: (1,1) sector dimensions at n_max=1 match Paper 51 Table (§6.1):
    J=0: 6, J=1: 14, J=2: 26, total 46.
    """
    fp = n_max_1_data["fp"]
    assert len(fp[0]) == 6
    assert len(fp[1]) == 14
    assert len(fp[2]) == 26
    assert len(fp[0]) + len(fp[1]) + len(fp[2]) == 46


def test_j_blindness_total_eigenvalue_set_equal_n_max_1(n_max_1_data):
    """Theorem thm:j_blindness, primary content: the SET of distinct S^(2)
    eigenvalues is the SAME across J = 0, 1, 2 at n_max = 1.
    """
    fp = n_max_1_data["fp"]
    W = n_max_1_data["W"]
    eigs_J0 = _s2_eigenvalues(fp[0], W)
    eigs_J1 = _s2_eigenvalues(fp[1], W)
    eigs_J2 = _s2_eigenvalues(fp[2], W)

    set_J0 = _eigvals_set(eigs_J0)
    set_J1 = _eigvals_set(eigs_J1)
    set_J2 = _eigvals_set(eigs_J2)

    assert set_J0 == set_J1, (
        f"J=0 and J=1 distinct eigenvalue sets differ:\n"
        f"  J=0: {set_J0}\n  J=1: {set_J1}"
    )
    assert set_J1 == set_J2, (
        f"J=1 and J=2 distinct eigenvalue sets differ:\n"
        f"  J=1: {set_J1}\n  J=2: {set_J2}"
    )


def test_within_sector_multiplicity_ratio_1_1_3_n_max_1(n_max_1_data):
    """Theorem thm:j_blindness, within-sector ratio: multiplicities 1:1:3
    across J=0,1,2 (Paper 51 §6.1 — Hermitianisation eliminates 4 of 9 modes
    per within-sector block, leaving 1+1+3 = 5 instead of 1+3+5 = 9).
    """
    fp = n_max_1_data["fp"]
    W = n_max_1_data["W"]
    D0 = n_max_1_data["basis"].D0

    mults_J = {}
    for J in (0, 1, 2):
        within, _ = _split_within_cross(fp[J], D0)
        eigs = _s2_eigenvalues(within, W)
        # All within-sector modes at n_max=1 land at one cluster (+0.1274...).
        clusters = _cluster(eigs)
        # Sum of multiplicities across distinct within-sector clusters:
        mults_J[J] = sum(m for _, m in clusters)

    # Paper claim: 1:1:3 cross-J for any within-sector eigenvalue.
    # At n_max=1, all within-sector modes share a single cluster
    # (only |lambda|=3/2 within-pair), so the totals already exhibit 1:1:3.
    base = mults_J[0]
    assert mults_J[0] == base
    assert mults_J[1] == base, f"Within-sector 1:1 violated: J=0={mults_J[0]}, J=1={mults_J[1]}"
    assert mults_J[2] == 3 * base, (
        f"Within-sector 1:3 violated: J=0={mults_J[0]}, J=2={mults_J[2]}, "
        f"expected J=2 = 3 * J=0 = {3 * base}"
    )


def test_cross_sector_multiplicity_ratio_1_3_5_n_max_1(n_max_1_data):
    """Theorem thm:j_blindness, cross-sector ratio: multiplicities 1:3:5
    across J=0,1,2 (Paper 51 §6.1 — full irrep dimensions dim(J) = 2J+1).
    """
    fp = n_max_1_data["fp"]
    W = n_max_1_data["W"]
    D0 = n_max_1_data["basis"].D0

    mults_J = {}
    for J in (0, 1, 2):
        _, cross = _split_within_cross(fp[J], D0)
        eigs = _s2_eigenvalues(cross, W)
        clusters = _cluster(eigs)
        mults_J[J] = sum(m for _, m in clusters)

    base = mults_J[0]
    assert mults_J[1] == 3 * base, (
        f"Cross-sector 1:3 violated: J=0={mults_J[0]}, J=1={mults_J[1]}, "
        f"expected J=1 = 3 * J=0 = {3 * base}"
    )
    assert mults_J[2] == 5 * base, (
        f"Cross-sector 1:5 violated: J=0={mults_J[0]}, J=2={mults_J[2]}, "
        f"expected J=2 = 5 * J=0 = {5 * base}"
    )


def test_analytical_weight_vs_stencil_n_max_1(n_max_1_data):
    """Paper 51 §6.1 claims analytical eq:s2_weight matches the 5-point central
    -difference stencil to relative error 1.24e-6 at n_max=2. We verify at
    n_max=1 on a sample of J=2 basis vectors at relative error < 1e-4
    (stencil noise sets the floor; the paper's tighter rel-err is at n_max=2,
    not the goal here — we just confirm the analytical S^(2) is self-consistent
    with first-principles directional second derivative).
    """
    fp = n_max_1_data["fp"]
    W = n_max_1_data["W"]
    D0 = n_max_1_data["basis"].D0
    Lambda_sq = n_max_1_data["Lambda_sq"]

    n_test = min(5, len(fp[2]))
    max_rel_err = 0.0
    for idx in range(n_test):
        V = fp[2][idx]
        s2_analytical = float(np.sum(W * V * V))
        s2_stencil = compute_s2_stencil(D0, V, Lambda_sq, eps=1e-4)
        if abs(s2_stencil) > 1e-10:
            rel_err = abs(s2_analytical - s2_stencil) / abs(s2_stencil)
            max_rel_err = max(max_rel_err, rel_err)
    assert max_rel_err < 1e-4, (
        f"Analytical vs stencil max relative error {max_rel_err:.2e} > 1e-4"
    )


@pytest.mark.slow
def test_j_blindness_at_n_max_2():
    """Stronger test at n_max=2 (dim_H = 40, takes ~30-60s).

    Sanity: (1,1) sector dimensions at n_max=2 per Paper 51 Table §6.1:
    J=0: 16, J=1: 40, J=2: 72, total 128.

    Then the same J-blindness check on the full eigenvalue set.
    """
    basis = CHBasis(n_max=2)
    fp = build_fp_basis(basis)
    Lambda_sq = 6.0
    W = s2_weight_matrix(basis.D0, Lambda_sq)

    # Dimensions per Paper 51 Table (n_max=2 row).
    assert len(fp[0]) == 16
    assert len(fp[1]) == 40
    assert len(fp[2]) == 72

    eigs_J0 = _s2_eigenvalues(fp[0], W)
    eigs_J1 = _s2_eigenvalues(fp[1], W)
    eigs_J2 = _s2_eigenvalues(fp[2], W)

    set_J0 = _eigvals_set(eigs_J0)
    set_J1 = _eigvals_set(eigs_J1)
    set_J2 = _eigvals_set(eigs_J2)

    assert set_J0 == set_J1, f"n_max=2 J=0/J=1 sets differ: {set_J0} vs {set_J1}"
    assert set_J1 == set_J2, f"n_max=2 J=1/J=2 sets differ: {set_J1} vs {set_J2}"
