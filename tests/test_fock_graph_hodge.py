"""Tests for geovac/fock_graph_hodge.py (Hodge decomposition of the S³ Fock graph).

All structural tests use exact sympy arithmetic.  Numerical
comparisons use numpy float64 with tol 1e-9.

Coverage:
  - L₀ = D − A identity (exact)
  - SVD identity: nonzero spectra of L₀ and L₁ match
  - Betti formula: β₁ = E − V + c
  - π-free certificate: all eigenvalues are rational integers
  - Incidence matrix shape, ±1 entries, column sum = 0
  - Node count V = n_max(n_max+1)(2n_max+1)/6
  - Edge count (regression values at n_max=2,3)
  - Spectrum comparison with continuum hodge1_s3
  - Bochner-Weitzenböck gap = 2n+1
"""

from __future__ import annotations

import pytest
import numpy as np
import sympy as sp
from sympy import Integer

from geovac.fock_graph_hodge import FockGraphHodge
from geovac.hodge1_s3 import hodge1_eigenvalue


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def hodge2():
    """FockGraphHodge at n_max=2 (sympy)."""
    return FockGraphHodge(2, use_sympy=True)


@pytest.fixture(scope="module")
def hodge3():
    """FockGraphHodge at n_max=3 (sympy)."""
    return FockGraphHodge(3, use_sympy=True)


# ---------------------------------------------------------------------------
# 1. Node / edge counts
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_max,expected_V", [
    (1, 1),
    (2, 5),
    (3, 14),
    (4, 30),
])
def test_node_count(n_max, expected_V):
    """V = Σₙ₌₁ⁿᵐᵃˣ n² = n_max(n_max+1)(2n_max+1)/6."""
    h = FockGraphHodge(n_max)
    assert h.n_nodes == expected_V


@pytest.mark.parametrize("n_max,expected_E", [
    # Compute by hand / from lattice:
    # n_max=2: (1,0,0)↔(2,0,0), (2,1,-1)↔(2,1,0)↔(2,1,1),
    #          (1,0,0)↔ no l transitions; radial (1,0,0)-(2,0,0)=1,
    #          angular within n=2,l=1: (m=-1)-(m=0), (m=0)-(m=1) = 2
    #          radial (2,l,m)↔(1,l,m) when l<1: only (2,0,0)-(1,0,0) already counted
    #          radial (2,1,-1)-(3,..) doesn't exist. So total = 1+2 = 3.
    #          Wait, also (1,0,0)-(2,0,0)=1, (2,1,-1)-(2,1,0)=1, (2,1,0)-(2,1,1)=1
    #          = 3 edges.
    (2, 3),
    # n_max=3: extend. For n=3: states (3,0,0),(3,1,-1),(3,1,0),(3,1,1),
    #          (3,2,-2),(3,2,-1),(3,2,0),(3,2,1),(3,2,2) = 9 states.
    # Angular at n=3,l=1: 2 edges; l=2: 4 edges = 6 angular
    # Angular at n=2,l=1: 2 edges (already counted in n_max=2)
    # Radial n=2->n=3 (same l,m): (2,0,0)-(3,0,0), (2,1,-1)-(3,1,-1),
    #   (2,1,0)-(3,1,0),(2,1,1)-(3,1,1) = 4 radial n=2->3
    # Radial n=1->n=2: 1 (already counted)
    # Total = 3 (from n_max=2) + 6 (new angular at n=3) + 4 (radial 2->3) = 13
    (3, 13),
])
def test_edge_count(n_max, expected_E):
    """Edge count regression at n_max=2,3."""
    h = FockGraphHodge(n_max)
    assert h.n_edges == expected_E


# ---------------------------------------------------------------------------
# 2. Incidence matrix structure
# ---------------------------------------------------------------------------

def test_incidence_shape(hodge2):
    """B must be V×E."""
    B = hodge2.incidence
    assert B.shape == (hodge2.n_nodes, hodge2.n_edges)


def test_incidence_entries_pm1_or_0(hodge2):
    """Every entry of B is in {-1, 0, +1}."""
    B = hodge2.incidence
    V, E = B.shape
    for i in range(V):
        for j in range(E):
            assert B[i, j] in (Integer(-1), Integer(0), Integer(1))


def test_incidence_column_sum_zero(hodge2):
    """Each column of B sums to 0 (one +1 and one -1 per edge)."""
    B = hodge2.incidence
    E = B.cols
    for j in range(E):
        col_sum = sum(B[i, j] for i in range(B.rows))
        assert col_sum == Integer(0), f"Column {j} sums to {col_sum}"


def test_incidence_column_nnz_two(hodge2):
    """Each column of B has exactly 2 non-zero entries (one +1, one -1)."""
    B = hodge2.incidence
    V, E = B.shape
    for j in range(E):
        nonzero = [i for i in range(V) if B[i, j] != 0]
        assert len(nonzero) == 2, f"Column {j} has {len(nonzero)} nonzero entries"


# ---------------------------------------------------------------------------
# 3. L₀ = D − A (exact)
# ---------------------------------------------------------------------------

def test_L0_equals_D_minus_A_nmax2(hodge2):
    """L₀ = B·Bᵀ equals D − A exactly at n_max=2."""
    assert hodge2.verify_L0_equals_D_minus_A()


def test_L0_equals_D_minus_A_nmax3(hodge3):
    """L₀ = B·Bᵀ equals D − A exactly at n_max=3."""
    assert hodge3.verify_L0_equals_D_minus_A()


def test_L0_symmetric(hodge2):
    """L₀ must be symmetric."""
    L0 = hodge2.node_laplacian
    V = L0.rows
    for i in range(V):
        for j in range(i + 1, V):
            assert L0[i, j] == L0[j, i], f"L0[{i},{j}] != L0[{j},{i}]"


def test_L0_row_sum_zero(hodge2):
    """Each row of L₀ sums to 0 (Laplacian property: L·1 = 0)."""
    L0 = hodge2.node_laplacian
    V = L0.rows
    for i in range(V):
        row_sum = sum(L0[i, j] for j in range(V))
        assert row_sum == Integer(0), f"Row {i} sums to {row_sum}"


def test_L0_diagonal_nonneg(hodge2):
    """Diagonal entries of L₀ must be non-negative (= degree)."""
    L0 = hodge2.node_laplacian
    for i in range(L0.rows):
        assert L0[i, i] >= 0


def test_L1_symmetric(hodge2):
    """L₁ = BᵀB must be symmetric."""
    L1 = hodge2.edge_laplacian
    E = L1.rows
    for i in range(E):
        for j in range(i + 1, E):
            assert L1[i, j] == L1[j, i], f"L1[{i},{j}] != L1[{j},{i}]"


# ---------------------------------------------------------------------------
# 4. SVD identity
# ---------------------------------------------------------------------------

def test_svd_identity_nmax2(hodge2):
    """Nonzero spectra of L₀ and L₁ agree at n_max=2."""
    assert hodge2.verify_svd_identity()


def test_svd_identity_nmax3(hodge3):
    """Nonzero spectra of L₀ and L₁ agree at n_max=3."""
    assert hodge3.verify_svd_identity()


def test_svd_identity_exact_sympy_nmax2(hodge2):
    """SVD identity holds in exact sympy arithmetic at n_max=2.

    The nonzero eigenvalues of L₀ and L₁ must be identical
    as sets (with multiplicity) in exact integer arithmetic.
    """
    ev0 = hodge2.node_laplacian_spectrum()
    ev1 = hodge2.edge_laplacian_spectrum()

    # Filter out zeros
    nz0 = [e for e in ev0 if e != Integer(0)]
    nz1 = [e for e in ev1 if e != Integer(0)]

    assert len(nz0) == len(nz1), (
        f"Nonzero multiplicity mismatch: |nz(L₀)|={len(nz0)}, |nz(L₁)|={len(nz1)}"
    )
    # Sort by float value for comparison
    nz0_sorted = sorted(nz0, key=float)
    nz1_sorted = sorted(nz1, key=float)
    for a, b in zip(nz0_sorted, nz1_sorted):
        assert a == b, f"Eigenvalue mismatch: L₀ has {a}, L₁ has {b}"


# ---------------------------------------------------------------------------
# 5. Betti numbers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_max,expected_c", [
    # The Fock graph decomposes into n_max angular-momentum sectors,
    # one per l = 0, 1, ..., n_max-1.  T± connects (n,l,m)↔(n±1,l,m)
    # (same l) and L± connects (n,l,m)↔(n,l,m±1) (same n,l).  States
    # with different l are never directly connected, so the graph has
    # exactly n_max connected components.
    (1, 1),
    (2, 2),
    (3, 3),
    (4, 4),
])
def test_betti_0_equals_n_max(n_max, expected_c):
    """β₀ = n_max (one connected component per angular momentum sector l)."""
    h = FockGraphHodge(n_max)
    assert h.betti_0 == expected_c


@pytest.mark.parametrize("n_max", [1, 2, 3, 4])
def test_betti_1_formula(n_max):
    """β₁ = E − V + c verifies for n_max = 1..4."""
    h = FockGraphHodge(n_max)
    assert h.verify_betti_formula()


@pytest.mark.parametrize("n_max,expected_beta1", [
    # n_max=1: V=1, E=0, c=1, β₁ = 0-1+1 = 0
    (1, 0),
    # n_max=2: V=5, E=3, c=2, β₁ = 3-5+2 = 0  (both sectors are trees)
    (2, 0),
    # n_max=3: V=14, E=13, c=3, β₁ = 13-14+3 = 2
    (3, 2),
    # n_max=4: V=30, E=34, c=4, β₁ = 34-30+4 = 8
    (4, 8),
])
def test_betti_1_values(n_max, expected_beta1):
    """Known β₁ values at small n_max (= E − V + n_max)."""
    h = FockGraphHodge(n_max)
    assert h.betti_1 == expected_beta1


# ---------------------------------------------------------------------------
# 6. π-free certificate and algebraic-integer property
# ---------------------------------------------------------------------------

def test_pi_free_nmax2(hodge2):
    """No eigenvalue of L₀ or L₁ contains π at n_max=2.

    At n_max=2 all eigenvalues happen to be rational integers because
    the l=1 sector is a simple path P₃ (m=-1,0,+1) whose eigenvalues
    are {0, 1, 3} and the l=0 sector is a single edge whose Laplacian
    has eigenvalues {0, 2}.
    """
    assert hodge2.verify_pi_free()


def test_pi_free_nmax3(hodge3):
    """No eigenvalue of L₀ or L₁ contains π at n_max=3.

    At n_max=3 the l=2 sector is a path P₅ (m=-2,...,2) whose
    eigenvalues involve √5 (algebraic irrationals), but no π.
    The π-free certificate passes because algebraic irrationals like
    √5 are not transcendentals.
    """
    assert hodge3.verify_pi_free()


def test_algebraic_integer_nmax2(hodge2):
    """All L₀ entries are integers → eigenvalues are algebraic integers."""
    from geovac.fock_graph_hodge import FockGraphHodge
    assert hodge2.verify_algebraic_integer()


def test_algebraic_integer_nmax3(hodge3):
    """All L₀ entries are integers → eigenvalues are algebraic integers."""
    assert hodge3.verify_algebraic_integer()


def test_l0_entries_are_integers(hodge3):
    """All entries of L₀ must be sympy.Integer (integer matrix)."""
    L0 = hodge3.node_laplacian
    for i in range(L0.rows):
        for j in range(L0.cols):
            assert isinstance(L0[i, j], (sp.Integer, int)), (
                f"L0[{i},{j}] = {L0[i,j]!r} is not an integer"
            )


def test_nmax3_has_irrational_eigenvalues(hodge3):
    """At n_max=3 some eigenvalues are irrational algebraic (l=2 path P₅ sector).

    This is NOT a failure — it documents that the Fock graph spectrum is
    algebraic over ℚ but not purely rational.  The l=2 sector at n_max=3
    is a path graph P₅ (m=-2,-1,0,1,2) whose eigenvalues are
    2 − 2cos(kπ/5) for k=0,...,4, which are quadratic irrationals
    involving √5 (specifically (3±√5)/2 and (5±√5)/2).

    These are algebraic integers — free of π, e, and all transcendentals —
    consistent with the π-free certificate.

    This test guards against regression if the algebraic structure changes.
    """
    evals = hodge3.node_laplacian_spectrum()
    # At least one eigenvalue must be irrational (not is_rational)
    irrational_evals = [
        ev for ev in evals
        if hasattr(ev, 'is_rational') and ev.is_rational is False
    ]
    assert len(irrational_evals) > 0, (
        "Expected irrational algebraic eigenvalues from l=2 path sector "
        "at n_max=3, but none found"
    )
    # Those irrational eigenvalues must not contain π
    for ev in irrational_evals:
        ev_expr = sp.sympify(ev)
        assert not ev_expr.has(sp.pi), (
            f"Irrational eigenvalue {ev} contains π — violates π-free certificate"
        )


def test_spectrum_nonneg_nmax2(hodge2):
    """All eigenvalues must be non-negative (L₀ and L₁ are PSD)."""
    for ev in hodge2.node_laplacian_spectrum():
        assert float(ev) >= -1e-12, f"Negative eigenvalue {ev}"
    for ev in hodge2.edge_laplacian_spectrum():
        assert float(ev) >= -1e-12, f"Negative eigenvalue {ev}"


def test_spectrum_nonneg_nmax3(hodge3):
    """All eigenvalues must be non-negative at n_max=3."""
    for ev in hodge3.node_laplacian_spectrum():
        assert float(ev) >= -1e-12, f"Negative eigenvalue {ev}"
    for ev in hodge3.edge_laplacian_spectrum():
        assert float(ev) >= -1e-12, f"Negative eigenvalue {ev}"


# ---------------------------------------------------------------------------
# 7. Known spectra (regression)
# ---------------------------------------------------------------------------

def test_L0_spectrum_nmax2_zero_count(hodge2):
    """L₀ at n_max=2 has exactly β₀ = 2 zero eigenvalues.

    The graph has 2 connected components (l=0 and l=1 sectors),
    so there are exactly 2 zero eigenvalues of L₀.
    """
    evals = hodge2.node_laplacian_spectrum()
    evals_float = sorted(float(e) for e in evals)
    n_zeros = sum(1 for e in evals_float if abs(e) < 1e-9)
    assert n_zeros == 2, f"Expected 2 zero eigenvalues (β₀=n_max=2), got {n_zeros}"


def test_L0_spectrum_nmax3_zero_count(hodge3):
    """L₀ at n_max=3 has exactly β₀ = 3 zero eigenvalues."""
    evals = hodge3.node_laplacian_spectrum()
    n_zeros = sum(1 for e in evals if float(e) < 1e-9)
    assert n_zeros == 3, f"Expected 3 zero eigenvalues (β₀=n_max=3), got {n_zeros}"


def test_L0_spectrum_zero_count_equals_betti0(hodge2, hodge3):
    """Number of zero eigenvalues of L₀ equals β₀ for both n_max=2,3."""
    for h in [hodge2, hodge3]:
        evals = h.node_laplacian_spectrum()
        n_zeros = sum(1 for e in evals if float(e) < 1e-9)
        assert n_zeros == h.betti_0, (
            f"n_max={h.n_max}: zero-count {n_zeros} != betti_0 {h.betti_0}"
        )


def test_L0_spectrum_algebraic_not_all_integer(hodge3):
    """At n_max=3 the L₀ spectrum contains irrational algebraic entries.

    Specifically, the l=2 sector (path P₅) contributes eigenvalues
    (3±√5)/2 and (5±√5)/2, confirming the spectrum is algebraic
    over ℚ but not purely rational.  This is the key structural finding:
    the graph spectrum is π-free but not integer-valued beyond l=1.
    """
    evals = hodge3.node_laplacian_spectrum()
    irrational_count = sum(
        1 for ev in evals
        if hasattr(ev, 'is_rational') and ev.is_rational is False
    )
    assert irrational_count == 4, (
        f"Expected 4 irrational eigenvalues from l=2 P₅ sector, got {irrational_count}"
    )


# ---------------------------------------------------------------------------
# 8. Continuum Hodge-1 comparison (Bochner-Weitzenböck gap)
# ---------------------------------------------------------------------------

def test_hodge1_comparison_gap_formula(hodge3):
    """For each scalar level n, the gap μₙ − λₙ = 2n+1.

    λₙ = n²−1 (scalar Laplacian / graph edge eigenvalue)
    μₙ = n(n+2) (continuum Hodge-1 eigenvalue)
    gap = μₙ − λₙ = n(n+2) − (n²−1) = 2n + 1
    """
    cmp = hodge3.compare_with_hodge1_spectrum()
    assert len(cmp) > 0, "No comparison entries returned"
    for entry in cmp:
        n = entry["n"]
        gap = entry["gap_ricci"]
        expected_gap = 2 * n + 1
        assert gap == expected_gap, (
            f"At n={n}: gap={gap}, expected 2n+1={expected_gap}"
        )


def test_hodge1_comparison_scalar_eigenvalue(hodge3):
    """Graph scalar eigenvalues must match λₙ = n²−1."""
    cmp = hodge3.compare_with_hodge1_spectrum()
    for entry in cmp:
        n = entry["n"]
        assert entry["lambda_scalar"] == n * n - 1


def test_hodge1_comparison_mu_value(hodge3):
    """Continuum Hodge-1 eigenvalues must match μₙ = n(n+2)."""
    cmp = hodge3.compare_with_hodge1_spectrum()
    for entry in cmp:
        n = entry["n"]
        assert entry["mu_hodge1"] == n * (n + 2)


def test_ricci_shift_positive(hodge3):
    """The continuum Hodge-1 eigenvalue μₙ > scalar eigenvalue λₙ for all n."""
    cmp = hodge3.compare_with_hodge1_spectrum()
    for entry in cmp:
        assert entry["gap_ricci"] > 0, (
            f"Ricci shift not positive at n={entry['n']}"
        )


# ---------------------------------------------------------------------------
# 9. Summary dict
# ---------------------------------------------------------------------------

def test_summary_dict_keys(hodge2):
    """Summary dict must contain all expected keys."""
    s = hodge2.summary()
    required_keys = [
        "n_max", "n_nodes", "n_edges", "betti_0", "betti_1",
        "betti_formula_check", "svd_identity_check",
        "L0_equals_D_minus_A", "pi_free_certificate",
        "node_spectrum_exact", "edge_spectrum_exact",
        "hodge1_comparison", "bochner_weitzenbock_note",
    ]
    for k in required_keys:
        assert k in s, f"Key '{k}' missing from summary"


def test_summary_all_checks_pass_nmax2(hodge2):
    """All verification checks must pass at n_max=2."""
    s = hodge2.summary()
    assert s["betti_formula_check"] is True
    assert s["svd_identity_check"] is True
    assert s["L0_equals_D_minus_A"] is True
    assert s["pi_free_certificate"] is True


def test_summary_all_checks_pass_nmax3(hodge3):
    """All verification checks must pass at n_max=3."""
    s = hodge3.summary()
    assert s["betti_formula_check"] is True
    assert s["svd_identity_check"] is True
    assert s["L0_equals_D_minus_A"] is True
    assert s["pi_free_certificate"] is True


# ---------------------------------------------------------------------------
# 10. Numerical consistency (numpy vs sympy)
# ---------------------------------------------------------------------------

def test_numpy_sympy_spectrum_agree_nmax2(hodge2):
    """Numpy and sympy spectra agree for L₀ at n_max=2."""
    ev_sym = sorted(float(e) for e in hodge2.node_laplacian_spectrum())
    ev_np = sorted(hodge2.node_laplacian_spectrum_numpy().tolist())
    assert len(ev_sym) == len(ev_np)
    for a, b in zip(ev_sym, ev_np):
        assert abs(a - b) < 1e-9, f"Discrepancy: sympy={a}, numpy={b}"


def test_numpy_sympy_spectrum_agree_nmax3(hodge3):
    """Numpy and sympy spectra agree for L₀ at n_max=3."""
    ev_sym = sorted(float(e) for e in hodge3.node_laplacian_spectrum())
    ev_np = sorted(hodge3.node_laplacian_spectrum_numpy().tolist())
    assert len(ev_sym) == len(ev_np)
    for a, b in zip(ev_sym, ev_np):
        assert abs(a - b) < 1e-9, f"Discrepancy: sympy={a}, numpy={b}"
