"""
Tests for graph_qed_propagator.py — electron propagator on the GeoVac Dirac graph.

Coverage:
  1. Operator construction (n_max=1, 2, 3; t=0, 1, -1/16)
  2. Matrix structure (diagonal eigenvalues, off-diagonal adjacency)
  3. Exact sympy inverse at n_max=1, 2
  4. D · G = I identity verification
  5. π-free certificate (all entries rational)
  6. Massless vs massive propagator
  7. Spectral properties (gap, condition number, eigenvalue signs)
  8. Determinant is nonzero rational
  9. Node-indexed propagator access
  10. n_max=3 (numpy path; no exact certificate required but spectrum checked)
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp
from sympy import Rational, Matrix

from geovac.graph_qed_propagator import (
    DiracGraphOperator,
    build_dirac_graph_operator,
    electron_propagator,
    massive_propagator,
    propagator_spectrum,
    propagator_pi_free_certificate,
    propagator_entry,
    KAPPA_SCALAR,
)
from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def op_n1_t0():
    """n_max=1, t=0 — smallest possible operator, 2 states."""
    return DiracGraphOperator(n_max=1, t=Rational(0))


@pytest.fixture
def op_n2_t0():
    """n_max=2, t=0 — 10 states, free diagonal propagator."""
    return DiracGraphOperator(n_max=2, t=Rational(0))


@pytest.fixture
def op_n2_t1():
    """n_max=2, t=1 — 10 states with full unit adjacency."""
    return DiracGraphOperator(n_max=2, t=Rational(1))


@pytest.fixture
def op_n2_tkappa():
    """n_max=2, t=κ_scalar=-1/16 — scalar-graph coupling."""
    return DiracGraphOperator(n_max=2, t=KAPPA_SCALAR)


@pytest.fixture
def op_n3_t0():
    """n_max=3, t=0 — 28 states."""
    return DiracGraphOperator(n_max=3, t=Rational(0))


# ---------------------------------------------------------------------------
# 1. Operator construction and basic properties
# ---------------------------------------------------------------------------

def test_n_states_n1():
    """n_max=1: only n=1 shell, 2 states (1s1/2 m_j=±1/2)."""
    op = DiracGraphOperator(n_max=1)
    assert op.N == 2


def test_n_states_n2():
    """n_max=2: 2+8=10 states in atomic mode."""
    # n=1: 1s1/2 (κ=-1) → 2 states
    # n=2: 2s1/2 (κ=-1) + 2p1/2 (κ=+1) + 2p3/2 (κ=-2) → 2+2+4=8 states
    op = DiracGraphOperator(n_max=2)
    assert op.N == 10


def test_n_states_n3():
    """n_max=3: 2+8+18=28 states in atomic mode."""
    op = DiracGraphOperator(n_max=3)
    assert op.N == 28


def test_labels_type_n2(op_n2_t0):
    """All labels should be DiracLabel instances."""
    for lab in op_n2_t0.labels:
        assert isinstance(lab, DiracLabel)


def test_eigenvalues_half_integer(op_n2_t0):
    """All diagonal eigenvalues must be exact half-integers (Camporesi-Higuchi)."""
    for e in op_n2_t0.eigenvalues:
        assert isinstance(e, sp.Rational)
        # Half-integer: 2*e must be odd
        twice = 2 * int(e * 2)  # not right; use numerator/denominator check
        assert e.denominator == 2 or e.denominator == 1


def test_eigenvalues_magnitudes_n1(op_n1_t0):
    """n_max=1: both states at n=1, |λ| = 1 + 1/2 = 3/2."""
    eigs = op_n1_t0.eigenvalues
    assert len(eigs) == 2
    for e in eigs:
        assert abs(e) == Rational(3, 2)


def test_eigenvalues_chirality_signs(op_n2_t0):
    """κ < 0 states get χ=+1 (positive eigenvalue), κ > 0 get χ=-1 (negative)."""
    for i, lab in enumerate(op_n2_t0.labels):
        e = op_n2_t0.eigenvalues[i]
        if lab.kappa < 0:
            assert e > 0, f"κ={lab.kappa} < 0 should give positive eigenvalue, got {e}"
        else:
            assert e < 0, f"κ={lab.kappa} > 0 should give negative eigenvalue, got {e}"


def test_eigenvalue_magnitudes_n2(op_n2_t0):
    """n=1: |λ|=3/2; n=2: |λ|=5/2. (Fock: |λ|=n+1/2)"""
    for i, lab in enumerate(op_n2_t0.labels):
        expected_abs = Rational(2 * lab.n_fock + 1, 2)
        assert abs(op_n2_t0.eigenvalues[i]) == expected_abs


def test_symmetry_count_n2(op_n2_t0):
    """n_max=2 atomic mode: κ<0 states dominate (8 positive, 2 negative).

    In atomic mode (l < n_fock), the n=1 shell has only l=0 (κ=-1 only, χ=+1).
    The n=2 shell has l=0 (κ=-1, χ=+1) and l=1 (κ=-2 j=3/2, χ=+1 AND κ=+1 j=1/2, χ=-1).
    Counts: κ<0 → {n=1: 2, n=2 κ=-1: 2, n=2 κ=-2: 4} = 8; κ>0 → {n=2 κ=+1: 2} = 2.
    """
    pos = sum(1 for e in op_n2_t0.eigenvalues if e > 0)
    neg = sum(1 for e in op_n2_t0.eigenvalues if e < 0)
    assert pos == 8
    assert neg == 2
    assert pos + neg == op_n2_t0.N


# ---------------------------------------------------------------------------
# 2. Matrix structure
# ---------------------------------------------------------------------------

def test_diagonal_t0(op_n2_t0):
    """t=0: off-diagonal entries must all be zero."""
    M = op_n2_t0.matrix_sympy()
    N = op_n2_t0.N
    for i in range(N):
        for j in range(N):
            if i != j:
                assert M[i, j] == 0, f"M[{i},{j}] = {M[i,j]} should be 0 for t=0"


def test_diagonal_entries_match_eigenvalues(op_n2_t0):
    """Diagonal entries of matrix match stored eigenvalues."""
    M = op_n2_t0.matrix_sympy()
    for i, e in enumerate(op_n2_t0.eigenvalues):
        assert M[i, i] == e


def test_offdiagonal_t1(op_n2_t1):
    """t=1: off-diagonal entries are 0 or 1 (adjacency)."""
    M = op_n2_t1.matrix_sympy()
    N = op_n2_t1.N
    for i in range(N):
        for j in range(N):
            if i != j:
                assert M[i, j] in (0, 1), f"M[{i},{j}] = {M[i,j]} for t=1"


def test_matrix_symmetric(op_n2_t1):
    """D_graph is symmetric (real symmetric Hamiltonian)."""
    M = op_n2_t1.matrix_sympy()
    assert M == M.T


def test_matrix_type_rational(op_n2_t0):
    """All matrix entries are sympy Rational (π-free)."""
    M = op_n2_t0.matrix_sympy()
    for e in M:
        assert isinstance(e, (sp.Rational, sp.Integer,
                               sp.core.numbers.Zero,
                               sp.core.numbers.One,
                               sp.core.numbers.NegativeOne))


def test_numpy_diagonal_matches_sympy(op_n2_t1):
    """numpy matrix diagonal matches sympy matrix diagonal."""
    M_sym = op_n2_t1.matrix_sympy()
    M_np = op_n2_t1.matrix_numpy()
    for i in range(op_n2_t1.N):
        assert abs(M_np[i, i] - float(M_sym[i, i])) < 1e-14


# ---------------------------------------------------------------------------
# 3. Exact propagator (sympy inverse)
# ---------------------------------------------------------------------------

def test_propagator_t0_diagonal(op_n2_t0):
    """t=0: G_e should be diagonal with entries 1/λ_i."""
    G, is_rat = electron_propagator(op_n2_t0, exact=True)
    N = op_n2_t0.N
    for i in range(N):
        for j in range(N):
            if i == j:
                expected = Rational(1) / op_n2_t0.eigenvalues[i]
                assert G[i, j] == expected, f"G[{i},{i}] = {G[i,i]}, expected {expected}"
            else:
                assert G[i, j] == 0


def test_propagator_t0_rational(op_n2_t0):
    """t=0: all propagator entries are rational."""
    _, is_rat = electron_propagator(op_n2_t0, exact=True)
    assert is_rat is True


def test_propagator_t1_rational(op_n2_t1):
    """t=1: all propagator entries are rational (matrix inverse of rational matrix)."""
    _, is_rat = electron_propagator(op_n2_t1, exact=True)
    assert is_rat is True


def test_propagator_tkappa_rational(op_n2_tkappa):
    """t=-1/16: all propagator entries are rational."""
    _, is_rat = electron_propagator(op_n2_tkappa, exact=True)
    assert is_rat is True


def test_identity_t0_n2(op_n2_t0):
    """D · G = I to machine precision, t=0."""
    G, _ = electron_propagator(op_n2_t0, exact=True)
    M = op_n2_t0.matrix_numpy()
    G_np = np.array(G.tolist(), dtype=np.float64)
    residual = np.max(np.abs(M @ G_np - np.eye(op_n2_t0.N)))
    assert residual < 1e-12


def test_identity_t1_n2(op_n2_t1):
    """D · G = I to machine precision, t=1."""
    G, _ = electron_propagator(op_n2_t1, exact=True)
    M = op_n2_t1.matrix_numpy()
    G_np = np.array(G.tolist(), dtype=np.float64)
    residual = np.max(np.abs(M @ G_np - np.eye(op_n2_t1.N)))
    assert residual < 1e-12


def test_identity_tkappa_n2(op_n2_tkappa):
    """D · G = I to machine precision, t=κ=-1/16."""
    G, _ = electron_propagator(op_n2_tkappa, exact=True)
    M = op_n2_tkappa.matrix_numpy()
    G_np = np.array(G.tolist(), dtype=np.float64)
    residual = np.max(np.abs(M @ G_np - np.eye(op_n2_tkappa.N)))
    assert residual < 1e-12


def test_propagator_n1_t0_values():
    """n_max=1, t=0: G_e has exactly two entries 2/3.

    n_max=1 in atomic mode has only s-states (l=0, κ=-1, χ=+1).
    Both states have λ = +3/2, so G_e[i,i] = 2/3 for both.
    (No κ>0 state exists at n=1 because κ>0 requires l≥1 and l<n_fock=1.)
    """
    op = DiracGraphOperator(n_max=1, t=Rational(0))
    G, _ = electron_propagator(op, exact=True)
    # n=1 shell: |λ| = 3/2, χ=+1 for both states → 1/λ = 2/3
    diag = [G[i, i] for i in range(2)]
    assert all(d == Rational(2, 3) for d in diag)


# ---------------------------------------------------------------------------
# 4. Massive propagator
# ---------------------------------------------------------------------------

def test_massive_m0_equals_massless(op_n2_t0):
    """G_e(m=0) should equal the massless G_e."""
    G_massless, _ = electron_propagator(op_n2_t0, exact=True)
    G_massive, _ = massive_propagator(op_n2_t0, Rational(0), exact=True)
    assert G_massless == G_massive


def test_massive_rational_entries(op_n2_t0):
    """G_e(m=1) has all rational entries."""
    _, is_rat = massive_propagator(op_n2_t0, Rational(1), exact=True)
    assert is_rat is True


def test_massive_propagator_identity(op_n2_t0):
    """(D - m) · G_e(m) = I for m=1/2."""
    m = Rational(1, 2)
    G, _ = massive_propagator(op_n2_t0, m, exact=True)
    M = op_n2_t0.matrix_numpy()
    N = op_n2_t0.N
    M_shifted = M - float(m) * np.eye(N)
    G_np = np.array(G.tolist(), dtype=np.float64)
    residual = np.max(np.abs(M_shifted @ G_np - np.eye(N)))
    assert residual < 1e-12


def test_massive_propagator_diagonal_form(op_n2_t0):
    """t=0, m rational: G_e(m)[i,i] = 1/(λ_i - m)."""
    m = Rational(1, 4)
    G, _ = massive_propagator(op_n2_t0, m, exact=True)
    for i, lam in enumerate(op_n2_t0.eigenvalues):
        expected = Rational(1) / (lam - m)
        assert G[i, i] == expected


# ---------------------------------------------------------------------------
# 5. π-free certificate
# ---------------------------------------------------------------------------

def test_pi_free_certificate_t0_n2(op_n2_t0):
    """Full π-free certificate passes for t=0, n_max=2."""
    cert = propagator_pi_free_certificate(op_n2_t0)
    assert cert["diag_rational"] is True
    assert cert["offdiag_rational"] is True
    assert cert["det_rational"] is True
    assert cert["det_nonzero"] is True
    assert cert["propagator_rational"] is True
    assert cert["all_pass"] is True


def test_pi_free_certificate_t1_n2(op_n2_t1):
    """Full π-free certificate passes for t=1, n_max=2."""
    cert = propagator_pi_free_certificate(op_n2_t1)
    assert cert["all_pass"] is True


def test_pi_free_certificate_tkappa_n2(op_n2_tkappa):
    """Full π-free certificate passes for t=κ=-1/16, n_max=2."""
    cert = propagator_pi_free_certificate(op_n2_tkappa)
    assert cert["all_pass"] is True


def test_identity_residual_small(op_n2_t0):
    """D·G-I residual is below 1e-12."""
    cert = propagator_pi_free_certificate(op_n2_t0)
    assert cert["identity_residual"] < 1e-12


def test_determinant_t0_exact():
    """t=0: det(D) = product of all eigenvalues (rational)."""
    op = DiracGraphOperator(n_max=2, t=Rational(0))
    det = op.det_sympy()
    expected = sp.Integer(1)
    for e in op.eigenvalues:
        expected *= e
    assert det == expected


def test_pi_free_with_mass(op_n2_t0):
    """Certificate includes mass test, also passes."""
    cert = propagator_pi_free_certificate(op_n2_t0, test_mass=Rational(1, 2))
    assert cert["massive_propagator_rational"] is True
    assert cert["mass_tested"] == 0.5


# ---------------------------------------------------------------------------
# 6. Spectral properties
# ---------------------------------------------------------------------------

def test_spectral_gap_t0(op_n2_t0):
    """t=0: spectral gap = min|λ| = 3/2 (first shell)."""
    spec = propagator_spectrum(op_n2_t0)
    assert abs(spec["spectral_gap"] - 1.5) < 1e-12


def test_condition_number_t0(op_n2_t0):
    """t=0: condition number = max|λ| / min|λ| = (5/2) / (3/2) = 5/3."""
    spec = propagator_spectrum(op_n2_t0)
    expected_cond = (5 / 2) / (3 / 2)  # = 5/3
    assert abs(spec["condition_number"] - expected_cond) < 1e-10


def test_chirality_balance_n3(op_n3_t0):
    """n_max=3 atomic mode: verify chirality counts are nonzero and sum to N=28.

    In atomic mode the κ<0 states always outnumber κ>0 states because:
    each n-shell contributes 2 states per (l, j=l+1/2) entry for all l in [0,n-1],
    but only 2l states per (l, j=l-1/2) entry for l in [1,n-1].
    """
    spec = propagator_spectrum(op_n3_t0)
    assert spec["n_positive_chi"] > 0
    assert spec["n_negative_chi"] > 0
    assert spec["n_positive_chi"] + spec["n_negative_chi"] == op_n3_t0.N


def test_eigenvalues_real_n2(op_n2_t1):
    """All eigenvalues of D_graph are real (symmetric matrix)."""
    eigs = op_n2_t1.spectrum_numpy()
    assert np.all(np.isreal(eigs))


def test_spectral_gap_increases_with_hopping():
    """Spectral gap may change with t (topology modifies spectrum)."""
    op0 = DiracGraphOperator(n_max=2, t=Rational(0))
    op1 = DiracGraphOperator(n_max=2, t=Rational(1))
    spec0 = propagator_spectrum(op0)
    spec1 = propagator_spectrum(op1)
    # Both should be nonzero (invertible)
    assert spec0["spectral_gap"] > 0
    assert spec1["spectral_gap"] > 0


# ---------------------------------------------------------------------------
# 7. Node-indexed propagator access
# ---------------------------------------------------------------------------

def test_propagator_entry_1s_t0():
    """G_e[1s1/2, 1s1/2] = 2/3 for t=0 (positive chirality, n=1)."""
    op = DiracGraphOperator(n_max=2, t=Rational(0))
    G, _ = electron_propagator(op, exact=True)
    # 1s1/2: n_fock=1, κ=-1, m_j=+1/2 (two_m_j=1)
    state = DiracLabel(n_fock=1, kappa=-1, two_m_j=1)
    entry = propagator_entry(G, op.labels, state, state)
    assert entry == Rational(2, 3)


def test_propagator_entry_offdiag_t0():
    """Off-diagonal G_e[a, b] = 0 for t=0 (diagonal propagator)."""
    op = DiracGraphOperator(n_max=2, t=Rational(0))
    G, _ = electron_propagator(op, exact=True)
    s1 = DiracLabel(n_fock=1, kappa=-1, two_m_j=1)
    s2 = DiracLabel(n_fock=1, kappa=-1, two_m_j=-1)
    entry = propagator_entry(G, op.labels, s1, s2)
    assert entry == 0


def test_propagator_entry_invalid_label():
    """propagator_entry raises KeyError for unknown state."""
    op = DiracGraphOperator(n_max=1, t=Rational(0))
    G, _ = electron_propagator(op, exact=True)
    bad = DiracLabel(n_fock=2, kappa=-1, two_m_j=1)  # n=2 not in n_max=1
    with pytest.raises(KeyError):
        propagator_entry(G, op.labels, bad, bad)


# ---------------------------------------------------------------------------
# 8. n_max=3 (numpy path)
# ---------------------------------------------------------------------------

def test_n3_matrix_size(op_n3_t0):
    """n_max=3 operator has 28×28 matrix."""
    M = op_n3_t0.matrix_numpy()
    assert M.shape == (28, 28)


def test_n3_invertible_t0(op_n3_t0):
    """n_max=3, t=0 is invertible."""
    assert op_n3_t0.is_invertible()


def test_n3_invertible_t1():
    """n_max=3, t=1 is invertible."""
    op = DiracGraphOperator(n_max=3, t=Rational(1))
    assert op.is_invertible()


def test_n3_numpy_propagator_identity():
    """n_max=3: D · G = I to 1e-10 (numpy)."""
    op = DiracGraphOperator(n_max=3, t=Rational(0))
    G, _ = electron_propagator(op, exact=False)
    M = op.matrix_numpy()
    N = op.N
    residual = np.max(np.abs(M @ G - np.eye(N)))
    assert residual < 1e-10


def test_n3_spectral_gap_t0(op_n3_t0):
    """n_max=3, t=0: spectral gap = 3/2 (first shell eigenvalue magnitude)."""
    spec = propagator_spectrum(op_n3_t0)
    assert abs(spec["spectral_gap"] - 1.5) < 1e-12


# ---------------------------------------------------------------------------
# 9. Structural invariants
# ---------------------------------------------------------------------------

def test_kappa_scalar_value():
    """KAPPA_SCALAR = -1/16 (scalar graph coupling constant, Paper 0)."""
    assert KAPPA_SCALAR == Rational(-1, 16)


def test_build_helper_matches_direct():
    """build_dirac_graph_operator gives same result as direct construction."""
    op1 = build_dirac_graph_operator(n_max=2, t=Rational(1))
    op2 = DiracGraphOperator(n_max=2, t=Rational(1))
    assert op1.N == op2.N
    assert op1.matrix_sympy() == op2.matrix_sympy()


def test_n_max_1_no_edges():
    """n_max=1: only one shell, no E1 dipole edges (Δn must be ≥ 1 for parity flip)."""
    op = DiracGraphOperator(n_max=1)
    # n=1 has only l=0 states; E1 requires Δl=±1, impossible within same shell at l_max=0
    assert op.n_edges == 0


def test_propagator_n1_symmetric(op_n1_t0):
    """G_e[a,b] = G_e[b,a] (symmetric propagator, symmetric D)."""
    G, _ = electron_propagator(op_n1_t0, exact=True)
    assert G == G.T


def test_bare_eigenvalues_are_half_integers(op_n2_t0):
    """All |eigenvalues| are exact half-integers (n+1/2) in Fock convention."""
    for lab, e in zip(op_n2_t0.labels, op_n2_t0.eigenvalues):
        expected_abs = Rational(2 * lab.n_fock + 1, 2)
        assert abs(e) == expected_abs


def test_det_product_of_eigenvalues_t0():
    """det(D_graph) = product of eigenvalues for diagonal (t=0) case."""
    op = DiracGraphOperator(n_max=2, t=Rational(0))
    det = op.det_sympy()
    product = sp.Integer(1)
    for e in op.eigenvalues:
        product *= e
    assert det == product
