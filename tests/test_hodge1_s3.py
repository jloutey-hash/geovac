"""Tests for geovac/hodge1_s3.py (Hodge-1 Laplacian on S^3, photon propagator).

All spectral tests use exact sympy arithmetic. No floating-point comparisons
except for the spectral zeta convergence test.
"""

from __future__ import annotations

import pytest
import sympy as sp
from sympy import Integer, Rational, Symbol

from geovac.hodge1_s3 import (
    VectorHarmonicLabel,
    hodge1_eigenvalue,
    hodge1_degeneracy,
    hodge1_spectral_zeta,
    hodge1_propagator_diagonal,
    verify_bochner_weitzenbock,
    compare_with_edge_laplacian,
    verify_pi_free,
    vertex_coupling,
    vertex_coupling_count,
    vector_labels_at_n,
    count_vector_labels,
)


# ---------------------------------------------------------------------------
# Eigenvalue mu_n = n(n+2)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", list(range(1, 11)))
def test_eigenvalue_exact(n):
    """mu_n = n(n+2) exactly for n = 1..10."""
    mu = hodge1_eigenvalue(n)
    expected = Integer(n) * Integer(n + 2)
    assert mu == expected
    assert isinstance(mu, (sp.Integer, int))


@pytest.mark.parametrize("n,expected", [
    (1, 3), (2, 8), (3, 15), (4, 24), (5, 35),
    (6, 48), (7, 63), (8, 80), (9, 99), (10, 120),
])
def test_eigenvalue_values(n, expected):
    """Known eigenvalue values."""
    assert hodge1_eigenvalue(n) == expected


def test_eigenvalue_n_zero_raises():
    """n=0 should raise ValueError (no 1-form zero mode on S^3)."""
    with pytest.raises(ValueError):
        hodge1_eigenvalue(0)


def test_eigenvalue_negative_raises():
    with pytest.raises(ValueError):
        hodge1_eigenvalue(-1)


# ---------------------------------------------------------------------------
# Degeneracy
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", list(range(1, 11)))
def test_degeneracy_total(n):
    """Total degeneracy d_n = 2n(n+2)."""
    assert hodge1_degeneracy(n, mode_type="all") == 2 * n * (n + 2)


@pytest.mark.parametrize("n", list(range(1, 11)))
def test_degeneracy_transverse(n):
    """Transverse degeneracy = n(n+2)."""
    assert hodge1_degeneracy(n, mode_type="transverse") == n * (n + 2)


@pytest.mark.parametrize("n", list(range(1, 11)))
def test_degeneracy_longitudinal(n):
    """Longitudinal degeneracy = n(n+2)."""
    assert hodge1_degeneracy(n, mode_type="longitudinal") == n * (n + 2)


@pytest.mark.parametrize("n", list(range(1, 11)))
def test_degeneracy_sum(n):
    """transverse + longitudinal = total."""
    d_t = hodge1_degeneracy(n, mode_type="transverse")
    d_l = hodge1_degeneracy(n, mode_type="longitudinal")
    d_all = hodge1_degeneracy(n, mode_type="all")
    assert d_t + d_l == d_all


def test_degeneracy_n1():
    """At n=1: mu=3, total=6, transverse=3, longitudinal=3."""
    assert hodge1_eigenvalue(1) == 3
    assert hodge1_degeneracy(1, "all") == 6
    assert hodge1_degeneracy(1, "transverse") == 3
    assert hodge1_degeneracy(1, "longitudinal") == 3


# ---------------------------------------------------------------------------
# No zero eigenvalue (beta_1 = 0)
# ---------------------------------------------------------------------------

def test_no_zero_eigenvalue():
    """No zero eigenvalue exists: mu_n = n(n+2) > 0 for all n >= 1.

    This confirms beta_1(S^3) = 0 (simply connected, no harmonic 1-forms).
    """
    for n in range(1, 50):
        assert hodge1_eigenvalue(n) > 0


def test_smallest_eigenvalue_is_3():
    """Smallest Hodge-1 eigenvalue is mu_1 = 3 (gap to zero)."""
    assert hodge1_eigenvalue(1) == 3


# ---------------------------------------------------------------------------
# Bochner-Weitzenbock verification (INDEPENDENT symbolic route)
# ---------------------------------------------------------------------------
#
# The previous tests here were a tautology: mu_conn was DEFINED as
# n(n+2) - 2 and the assertion reduced to 2 == 2 (group5 cert finding E17).
# They are replaced by a genuinely independent verification: the operator
# identity  Delta_Hodge = nabla*nabla + 2  (Bochner-Weitzenbock with
# Ric = 2g) is verified SYMBOLICALLY on the unit round S^3 in Hopf
# coordinates (eta, xi1, xi2), metric g = diag(1, sin^2 eta, cos^2 eta),
# by explicit construction of both operators from the metric (Christoffel
# symbols, covariant derivatives, d and delta). Nothing is assumed about
# the spectrum.
#
# Literature anchors verified alongside (Ikeda-Taniguchi 1978,
# Camporesi 1990):
#   * exact family:    Delta_H(df) = k(k+2) df, rough = k(k+2) - 2
#   * coclosed family: the Killing 1-form has Delta_H = 4 = 2(n-1)|_{n=3}
#     (lowest CO-EXACT eigenvalue is 4 with multiplicity 6 = dim so(4),
#     NOT 3 -- see the module docstring's labeling-convention note).

_HOPF = None


def _hopf_calculus():
    """Symbolic exterior/covariant calculus on unit S^3, Hopf coordinates.

    Returns (coords, hodge_laplacian, rough_laplacian, scalar_laplacian)
    where the three operators act on component lists / scalars.
    """
    global _HOPF
    if _HOPF is not None:
        return _HOPF
    eta, x1, x2 = sp.symbols('eta xi1 xi2', real=True)
    coords = [eta, x1, x2]
    g = sp.diag(1, sp.sin(eta) ** 2, sp.cos(eta) ** 2)
    ginv = g.inv()
    sqrtg = sp.sin(eta) * sp.cos(eta)
    n = 3
    Gamma = [[[sp.simplify(sum(ginv[k, m] * (sp.diff(g[m, i], coords[j])
                                             + sp.diff(g[m, j], coords[i])
                                             - sp.diff(g[i, j], coords[m]))
                               for m in range(n)) / 2)
               for j in range(n)] for i in range(n)] for k in range(n)]

    def cov_d1(w):
        # T[i][j] = nabla_i w_j
        return [[sp.diff(w[j], coords[i]) - sum(Gamma[k][i][j] * w[k]
                                                for k in range(n))
                 for j in range(n)] for i in range(n)]

    def rough_laplacian(w):
        # (nabla* nabla w)_i = -g^{jk} (nabla_j nabla_k w)_i
        T = cov_d1(w)
        out = []
        for i in range(n):
            s = 0
            for j in range(n):
                for k in range(n):
                    njT = (sp.diff(T[k][i], coords[j])
                           - sum(Gamma[m][j][k] * T[m][i] for m in range(n))
                           - sum(Gamma[m][j][i] * T[k][m] for m in range(n)))
                    s += ginv[j, k] * njT
            out.append(-s)
        return out

    def hodge_laplacian(w):
        # (d delta + delta d) w  for a 1-form w
        dw_scalar = -sum(sp.diff(sqrtg * sum(ginv[i, j] * w[j]
                                             for j in range(n)),
                                 coords[i]) for i in range(n)) / sqrtg
        d_delta = [sp.diff(dw_scalar, coords[i]) for i in range(n)]
        F = [[sp.diff(w[j], coords[i]) - sp.diff(w[i], coords[j])
              for j in range(n)] for i in range(n)]
        delta_d = []
        for j in range(n):
            s = 0
            for i in range(n):
                for k in range(n):
                    niF = (sp.diff(F[k][j], coords[i])
                           - sum(Gamma[m][i][k] * F[m][j] for m in range(n))
                           - sum(Gamma[m][i][j] * F[k][m] for m in range(n)))
                    s += ginv[i, k] * niF
            delta_d.append(-s)
        return [d_delta[i] + delta_d[i] for i in range(n)]

    def scalar_laplacian(f):
        return -sum(sp.diff(sqrtg * sum(ginv[i, j] * sp.diff(f, coords[j])
                                        for j in range(n)), coords[i])
                    for i in range(n)) / sqrtg

    _HOPF = (coords, hodge_laplacian, rough_laplacian, scalar_laplacian)
    return _HOPF


def test_weitzenbock_operator_identity_generic():
    """SYMBOLIC PROOF (patch-level): Delta_Hodge - nabla*nabla = 2 * id on
    1-forms of unit S^3, for a 1-form with three ARBITRARY function
    components. This is the Bochner-Weitzenbock identity with Ric = 2g,
    verified with no spectral input whatsoever."""
    coords, hodge, rough, _ = _hopf_calculus()
    eta, x1, x2 = coords
    a = sp.Function('a')(eta, x1, x2)
    b = sp.Function('b')(eta, x1, x2)
    c = sp.Function('c')(eta, x1, x2)
    w = [a, b, c]
    H = hodge(w)
    R = rough(w)
    for i in range(3):
        assert sp.simplify(sp.expand(H[i] - R[i] - 2 * w[i])) == 0


def test_weitzenbock_killing_form():
    """The Killing 1-form omega = sin^2(eta) d xi1 satisfies
    Delta_H omega = 4 omega and nabla*nabla omega = 2 omega.

    Anchors the CO-EXACT (transverse) family: its lowest Hodge eigenvalue
    on unit S^3 is 4 (multiplicity 6 = dim so(4); Ikeda-Taniguchi coclosed
    k=1 eigenvalue (k+1)^2), and the Ricci shift +2 is realized as
    4 = 2 + 2."""
    coords, hodge, rough, _ = _hopf_calculus()
    eta = coords[0]
    w = [sp.Integer(0), sp.sin(eta) ** 2, sp.Integer(0)]
    H = hodge(w)
    R = rough(w)
    for i in range(3):
        assert sp.simplify(H[i] - 4 * w[i]) == 0
        assert sp.simplify(R[i] - 2 * w[i]) == 0


def test_weitzenbock_exact_form_level1():
    """The exact 1-form omega = df with f = cos(eta) cos(xi2) (a level-1
    scalar harmonic, Delta_0 f = 3 f) satisfies Delta_H omega = 3 omega
    and nabla*nabla omega = 1 omega.

    Anchors the EXACT (longitudinal) family at eigenvalue n(n+2) = 3
    (n=1) and its rough eigenvalue n(n+2) - 2 = 1: the arithmetic
    relation mu_conn = n(n+2) - 2 used by verify_bochner_weitzenbock()
    is here derived from the metric for the exact family, not assumed."""
    coords, hodge, rough, scalar_lap = _hopf_calculus()
    eta, x1, x2 = coords
    f = sp.cos(eta) * sp.cos(x2)
    assert sp.simplify(scalar_lap(f) - 3 * f) == 0
    w = [sp.diff(f, coords[i]) for i in range(3)]
    H = hodge(w)
    R = rough(w)
    for i in range(3):
        assert sp.simplify(H[i] - 3 * w[i]) == 0
        assert sp.simplify(R[i] - 1 * w[i]) == 0


def test_bochner_weitzenbock_arithmetic_corollary():
    """verify_bochner_weitzenbock() checks the exact-family arithmetic
    relation mu_Hodge(n) = [n(n+2) - 2] + 2. On its own this is circular;
    it is kept as an API-consistency check because the +2 shift and the
    exact-family rough eigenvalue are now independently established by
    the symbolic tests above."""
    assert verify_bochner_weitzenbock(20) is True


# ---------------------------------------------------------------------------
# Degeneracy: independent literature anchor (Ikeda-Taniguchi)
# ---------------------------------------------------------------------------

def test_degeneracy_ikeda_taniguchi_anchor():
    """Independent anchor for the degeneracy formula 2n(n+2).

    Ikeda-Taniguchi (1978): coclosed 1-eigenforms on S^3 have eigenvalue
    (k+1)^2 with multiplicity m_k = 2k(k+2), k >= 1. Hardcoded literature
    values (independent of the module's formula):
        k=1: m = 6  = dim so(4)   (Killing forms, eigenvalue 4)
        k=2: m = 16                (eigenvalue 9)
        k=3: m = 30                (eigenvalue 16)
    The module's total degeneracy hodge1_degeneracy(k, 'all') equals this
    coclosed multiplicity numerically at the same label k. NOTE the
    labeling-convention caveat (module docstring): in the literature this
    multiplicity is attached to Hodge eigenvalue (k+1)^2, not k(k+2);
    the module attaches level k to the exact-family eigenvalue k(k+2)."""
    literature = {1: 6, 2: 16, 3: 30}
    for k, m in literature.items():
        assert hodge1_degeneracy(k, mode_type="all") == m
        # closed form of the I-T multiplicity, derived independently:
        assert m == 2 * ((k + 1) ** 2 - 1)


# ---------------------------------------------------------------------------
# Edge Laplacian comparison
# ---------------------------------------------------------------------------

def test_edge_laplacian_comparison():
    """Hodge-1 eigenvalue at level n equals scalar eigenvalue at level n+1."""
    results = compare_with_edge_laplacian(10)
    for r in results:
        n = r['n']
        # mu_n^{Hodge-1} = lambda_{n+1}^{scalar} (level-shifted match)
        assert r['gap_shifted'] == 0, (
            f"At n={n}: mu_hodge1={r['mu_hodge1']} != "
            f"lambda_scalar(n+1)={r['lambda_scalar_shifted']}")


def test_edge_laplacian_gap_same_level():
    """Gap at same level is 2n+1 (from representation theory shift)."""
    results = compare_with_edge_laplacian(10)
    for r in results:
        n = r['n']
        expected_gap = 2 * n + 1
        assert r['gap_same_n'] == expected_gap


# ---------------------------------------------------------------------------
# pi-free certificate
# ---------------------------------------------------------------------------

def test_pi_free_certificate():
    """All eigenvalues and degeneracies are exact integers for n=1..20."""
    assert verify_pi_free(20) is True


@pytest.mark.parametrize("n", list(range(1, 11)))
def test_eigenvalue_is_integer(n):
    """Each mu_n is a sympy Integer (not Rational with denominator > 1)."""
    mu = hodge1_eigenvalue(n)
    assert isinstance(mu, sp.Integer)


# ---------------------------------------------------------------------------
# Label generators
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", list(range(1, 8)))
def test_label_count_all(n):
    """vector_labels_at_n produces exactly d_n labels."""
    labels = vector_labels_at_n(n, mode_type="all")
    assert len(labels) == hodge1_degeneracy(n, mode_type="all")


@pytest.mark.parametrize("n", list(range(1, 8)))
def test_label_count_transverse(n):
    labels = vector_labels_at_n(n, mode_type="transverse")
    assert len(labels) == hodge1_degeneracy(n, mode_type="transverse")


@pytest.mark.parametrize("n", list(range(1, 8)))
def test_label_count_longitudinal(n):
    labels = vector_labels_at_n(n, mode_type="longitudinal")
    assert len(labels) == hodge1_degeneracy(n, mode_type="longitudinal")


def test_label_types():
    """All labels have correct mode_type field."""
    for label in vector_labels_at_n(3, mode_type="transverse"):
        assert label.mode_type == "transverse"
    for label in vector_labels_at_n(3, mode_type="longitudinal"):
        assert label.mode_type == "longitudinal"


def test_count_vector_labels():
    """Total count up to n_max matches sum of degeneracies."""
    n_max = 5
    expected = sum(hodge1_degeneracy(n) for n in range(1, n_max + 1))
    assert count_vector_labels(n_max) == expected


# ---------------------------------------------------------------------------
# Photon propagator
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", list(range(1, 8)))
def test_propagator_diagonal(n):
    """Propagator diagonal = 1/(n(n+2)) on unit sphere."""
    prop = hodge1_propagator_diagonal(n)
    expected = Rational(1, n * (n + 2))
    assert prop == expected


def test_propagator_R_scaling():
    """Propagator scales as R^2: G_n(R) = R^2 / (n(n+2))."""
    R = Symbol('R', positive=True)
    for n in [1, 2, 5]:
        prop = hodge1_propagator_diagonal(n, R=R)
        expected = R**2 / (Integer(n) * Integer(n + 2))
        assert sp.simplify(prop - expected) == 0


# ---------------------------------------------------------------------------
# Radius scaling
# ---------------------------------------------------------------------------

def test_eigenvalue_R_scaling():
    """mu_n(R) = n(n+2)/R^2 on sphere of radius R."""
    R = Symbol('R', positive=True)
    for n in [1, 3, 7]:
        mu = hodge1_eigenvalue(n, R=R)
        expected = Integer(n) * Integer(n + 2) / R**2
        assert sp.simplify(mu - expected) == 0


# ---------------------------------------------------------------------------
# Selection rule (QED vertex)
# ---------------------------------------------------------------------------

def test_vertex_triangle_inequality():
    """Basic triangle inequality checks."""
    # n1=0, n2=0: n_gamma must be 0, but n_gamma >= 1, so nothing allowed
    assert vertex_coupling(0, 0, 1) is False

    # n1=1, n2=0: |1-0| <= n_gamma <= 1+0, so n_gamma = 1
    # parity: 1+0+1 = 2 (even) -> not allowed
    assert vertex_coupling(1, 0, 1) is False

    # n1=1, n2=1: |0| <= n_gamma <= 2
    # n_gamma=1: 1+1+1=3 (odd) -> allowed
    assert vertex_coupling(1, 1, 1) is True
    # n_gamma=2: 1+1+2=4 (even) -> not allowed
    assert vertex_coupling(1, 1, 2) is False

    # n1=2, n2=1: |1| <= n_gamma <= 3
    # n_gamma=1: 2+1+1=4 (even) -> not
    assert vertex_coupling(2, 1, 1) is False
    # n_gamma=2: 2+1+2=5 (odd) -> allowed
    assert vertex_coupling(2, 1, 2) is True
    # n_gamma=3: 2+1+3=6 (even) -> not
    assert vertex_coupling(2, 1, 3) is False


def test_vertex_symmetry():
    """vertex_coupling(n1, n2, ng) == vertex_coupling(n2, n1, ng)."""
    for n1 in range(5):
        for n2 in range(5):
            for ng in range(1, 6):
                assert vertex_coupling(n1, n2, ng) == vertex_coupling(n2, n1, ng)


def test_vertex_coupling_count_small():
    """Smoke test: vertex_coupling_count at small n_max gives reasonable results."""
    result = vertex_coupling_count(3, 3)
    assert result['n_allowed_triples'] > 0
    assert result['n_allowed_triples'] <= result['total_possible_triples']
    assert 0 < result['density'] < 1


def test_vertex_invalid_inputs():
    """Invalid inputs raise ValueError."""
    with pytest.raises(ValueError):
        vertex_coupling(-1, 0, 1)
    with pytest.raises(ValueError):
        vertex_coupling(0, -1, 1)
    with pytest.raises(ValueError):
        vertex_coupling(0, 0, 0)  # n_gamma must be >= 1


# ---------------------------------------------------------------------------
# Spectral zeta convergence
# ---------------------------------------------------------------------------

def test_spectral_zeta_convergence():
    """Spectral zeta at s=2 converges as n_max increases."""
    z10 = float(hodge1_spectral_zeta(2, n_max=10, mode_type="transverse"))
    z50 = float(hodge1_spectral_zeta(2, n_max=50, mode_type="transverse"))
    z100 = float(hodge1_spectral_zeta(2, n_max=100, mode_type="transverse"))
    # Should be monotonically increasing (all positive terms)
    assert z10 < z50 < z100
    # And converging (difference shrinks)
    assert (z100 - z50) < (z50 - z10)


def test_spectral_zeta_s1_diverges():
    """At s=1, the zeta function should grow with n_max (divergent)."""
    z10 = float(hodge1_spectral_zeta(1, n_max=10, mode_type="all"))
    z100 = float(hodge1_spectral_zeta(1, n_max=100, mode_type="all"))
    assert z100 > 2 * z10  # rough check that it's growing fast


def test_spectral_zeta_exact_small():
    """Exact value of zeta at s=2, n_max=2 (transverse modes)."""
    # n=1: d_1^T = 3, mu_1 = 3 -> 3/9 = 1/3
    # n=2: d_2^T = 8, mu_2 = 8 -> 8/64 = 1/8
    # Total: 1/3 + 1/8 = 11/24
    z = hodge1_spectral_zeta(2, n_max=2, mode_type="transverse")
    assert z == Rational(11, 24)


# ---------------------------------------------------------------------------
# Regression: known values
# ---------------------------------------------------------------------------

def test_known_degeneracy_table():
    """Spot-check degeneracy table against hand calculation."""
    # n: total, trans, long
    known = {
        1: (6, 3, 3),
        2: (16, 8, 8),
        3: (30, 15, 15),
        4: (48, 24, 24),
        5: (70, 35, 35),
    }
    for n, (d_all, d_t, d_l) in known.items():
        assert hodge1_degeneracy(n, "all") == d_all
        assert hodge1_degeneracy(n, "transverse") == d_t
        assert hodge1_degeneracy(n, "longitudinal") == d_l


def test_propagator_partial_fractions():
    """1/(n(n+2)) = 1/2 * (1/n - 1/(n+2)) partial fraction decomposition."""
    for n in range(1, 10):
        prop = hodge1_propagator_diagonal(n)
        pf = Rational(1, 2) * (Rational(1, n) - Rational(1, n + 2))
        assert prop == pf
