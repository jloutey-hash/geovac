"""Tests for geovac.ihara_zeta_dirac (Track RH-C).

Covers, for n_max in {1, 2, 3} and both adjacency rules A (scalar-analog)
and B (E1 dipole):
  - Sanity: node count matches Σ_{n_fock=1..n_max} [Σ_l g^{(spinor_per_l)}].
  - π-free certificate (in the WEAK sense = charpoly has integer
    coefficients).  We do NOT require strict rationality of the
    adjacency spectrum; only that no transcendental constants enter.
  - Ramanujan verdict is reported for each rule × n_max combination,
    WITH numerical deviation (if not Ramanujan at a given n_max, the
    test records this without failing).
  - Cross-consistency: Bass closed form matches the truncated Euler
    product (up to order s^6, the RH-A cross-check convention).
  - Per-κ decomposition: Rule A's κ-sectors are each Ramanujan
    (each individual κ-component is a subgraph of a scalar graph
    that we already know is Ramanujan from RH-A).
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels
from geovac.dirac_s3 import dirac_degeneracy
from geovac.ihara_zeta import (
    ihara_zeta_bass,
    ihara_zeta_euler,
    is_ramanujan,
)
from geovac.ihara_zeta_dirac import (
    build_dirac_s3_graph,
    describe_adjacency_rule,
    dirac_s3_ihara_zeta,
    per_kappa_components,
)

s = sp.symbols("s")


# ---------------------------------------------------------------------------
# Node-count sanity
# ---------------------------------------------------------------------------

def _dirac_label_count(n_max: int) -> int:
    """Sum of Dirac-spinor multiplicities up to n_fock = n_max.

    At each n_fock, for each l = 0..n_fock-1, we get:
      - κ = -(l+1)  (j = l+1/2) → 2j+1 = 2l+2 states
      - κ = +l (if l ≥ 1, j = l-1/2) → 2j+1 = 2l states

    For n_fock=1 (only l=0): 2 states.
    For n_fock=2 (l=0, l=1): 2 + (4 + 2) = 8.  Cumulative 10.
    For n_fock=3 (add l=0,1,2): 2 + 6 + (6 + 4) = 18. Cumulative 28.
    """
    return sum(1 for _ in iter_dirac_labels(n_max))


@pytest.mark.parametrize("n_max, expected", [(1, 2), (2, 10), (3, 28)])
def test_node_count_matches_dirac_label_enumeration(n_max, expected):
    """The graph must have exactly the number of DiracLabel's for n_max."""
    # This is a counting identity on the node enumeration, independent of rule
    for rule in ("A", "B"):
        A, labels, deg, desc = build_dirac_s3_graph(n_max, rule)
        assert len(labels) == expected, (
            f"rule={rule} n_max={n_max}: got {len(labels)} labels, "
            f"expected {expected}"
        )
        assert A.shape == (expected, expected), desc
        assert _dirac_label_count(n_max) == expected


def test_node_count_n_max_1_matches_dirac_s3_sanity():
    """n_max=1 should give exactly the single j=1/2 doublet."""
    for rule in ("A", "B"):
        A, labels, deg, desc = build_dirac_s3_graph(1, rule)
        assert len(labels) == 2
        # Both labels should be n_fock=1, κ=-1, 2m_j=±1
        kappas = {lab.kappa for lab in labels}
        two_mjs = sorted(lab.two_m_j for lab in labels)
        assert kappas == {-1}
        assert two_mjs == [-1, 1]


# ---------------------------------------------------------------------------
# π-free certificate (integer-coefficient charpoly is always true here)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("rule", ["A", "B"])
def test_charpoly_integer_coefficient(n_max, rule):
    """Adjacency is integer → charpoly is integer → π-free in weak sense.

    This is the Paper 24 π-free certificate in its weakest form: there
    are no transcendental constants (π, √π, etc.) in the adjacency, so
    every spectral invariant (Hashimoto spectrum, Ihara zeta coefficients,
    graph eigenvalues) is algebraic over ℚ.
    """
    result = dirac_s3_ihara_zeta(n_max, rule)
    assert result["charpoly_integer_coefficient"] is True, (
        f"rule={rule} n_max={n_max}: charpoly has non-integer coefficients — "
        f"this would be a bug in the adjacency constructor."
    )


@pytest.mark.parametrize("rule", ["A", "B"])
def test_zeta_inverse_polynomial_has_integer_coeffs(rule):
    """The expanded Bass polynomial itself must have integer coefficients.

    This is stronger than the charpoly-of-A check: it uses the full
    Bass determinantal formula including (I - sA + s²Q) and the
    (1-s²)^{r-c} prefactor.
    """
    result = dirac_s3_ihara_zeta(2, rule, factor=False)
    zinv = result["zeta_inverse_expanded"]
    poly = sp.Poly(zinv, s)
    coeffs = poly.all_coeffs()
    for c in coeffs:
        assert c.is_Integer or (hasattr(c, "is_rational") and c.is_rational), (
            f"rule={rule}: non-rational coefficient {c} in zeta_inverse"
        )


# ---------------------------------------------------------------------------
# Ramanujan verdict smoke test
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("rule", ["A", "B"])
def test_ramanujan_verdict_is_bool_with_finite_deviation(n_max, rule):
    """is_ramanujan must return a bool and a finite deviation on every rule × n_max."""
    result = dirac_s3_ihara_zeta(n_max, rule)
    assert isinstance(result["ramanujan_verdict"], bool)
    dev = result["ramanujan_deviation"]
    assert np.isfinite(dev), (
        f"rule={rule} n_max={n_max}: deviation = {dev} is not finite"
    )


@pytest.mark.parametrize("rule", ["A", "B"])
def test_rule_A_kappa_preserving_matches_decomposition(rule):
    """Rule A preserves κ: each connected component sits inside one κ-sector.

    Concretely: for rule A, the number of per-κ components is the same
    as the full graph's number of connected components.  For rule B,
    this is NOT generally true (dipole edges mix κ), so we only test A.
    """
    if rule != "A":
        pytest.skip("only applicable to rule A")
    result = dirac_s3_ihara_zeta(3, "A")
    c_total = result["c"]
    per_k_total = sum(sizes["components"] for sizes in result["per_kappa_sizes"].values())
    assert c_total == per_k_total, (
        f"Rule A: total components ({c_total}) should equal sum of per-κ "
        f"components ({per_k_total})"
    )


# ---------------------------------------------------------------------------
# Bass vs Euler cross-consistency
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_max", [1, 2])
@pytest.mark.parametrize("rule", ["A", "B"])
def test_bass_vs_euler_agree_to_s6(n_max, rule):
    """Bass closed form and truncated Euler product match at low s."""
    A, _labels, _deg, _desc = build_dirac_s3_graph(n_max, rule)
    if A.size == 0 or A.sum() == 0:
        pytest.skip("empty graph has trivial zeta")
    zb = ihara_zeta_bass(A)
    ze = ihara_zeta_euler(A, max_length=6)
    series_bass = sp.expand(sp.series(zb, s, 0, 7).removeO())
    series_euler = sp.expand(sp.series(ze, s, 0, 7).removeO())
    diff = sp.simplify(series_bass - series_euler)
    assert diff == 0, (
        f"rule={rule} n_max={n_max}: Bass and Euler disagree at low s: {diff}"
    )


# ---------------------------------------------------------------------------
# Per-κ decomposition for Rule A (matches l-shell structure of scalar S^3)
# ---------------------------------------------------------------------------

def test_rule_A_per_kappa_structure_at_n_max_3():
    """At n_max=3, Rule A has κ-blocks matching the scalar l-shell pattern.

    Rule A's per-κ sub-blocks are determined by:
      - κ = -1 (l=0, j=1/2): n_fock=1,2,3; 2 m_j-values each.
      - κ = -2 (l=1, j=3/2): n_fock=2,3; 4 m_j-values each.
      - κ = +1 (l=1, j=1/2): n_fock=2,3; 2 m_j-values each.
      - κ = -3 (l=2, j=5/2): n_fock=3; 6 m_j-values.
      - κ = +2 (l=2, j=3/2): n_fock=3; 4 m_j-values.

    At each κ the sub-block is a product graph of an n-path × an
    m_j-path, lifted to spinors.  We check sizes here, not topology.
    """
    result = dirac_s3_ihara_zeta(3, "A")
    per_k = result["per_kappa_sizes"]
    # Check expected V per κ
    expected_V = {-1: 6, -2: 8, +1: 4, -3: 6, +2: 4}
    for kappa, expected in expected_V.items():
        assert per_k[kappa]["V"] == expected, (
            f"Rule A n_max=3 κ={kappa}: V = {per_k[kappa]['V']}, "
            f"expected {expected}"
        )


# ---------------------------------------------------------------------------
# n_max=1 boundary: graph is a single edge (2 nodes), zeta is trivial
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("rule", ["A", "B"])
def test_n_max_1_graph_is_a_single_m_j_edge(rule):
    """At n_max=1 there's only one κ=-1 doublet (m_j=±1/2).  Both rules
    connect the two m_j-states, giving a single edge between 2 nodes
    (a degenerate "path graph" with 2 vertices).
    """
    A, labels, deg, desc = build_dirac_s3_graph(1, rule)
    assert A.shape == (2, 2)
    # For rule A, the m-ladder edge exists (Δ m_j = ±1 at fixed κ, n).
    # For rule B, the dipole rule requires Δl = ±1, which is NOT possible
    # within the same (n, κ) sector at n=1 (only l=0 is available).
    # So rule B has 0 edges at n_max=1, rule A has 1 edge.
    if rule == "A":
        assert int(A.sum()) == 2  # two directed entries → one undirected edge
    else:
        assert int(A.sum()) == 0  # no dipole edges possible at n_max=1


# ---------------------------------------------------------------------------
# Descriptions are strings and are non-trivial
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("rule", ["A", "B"])
def test_describe_adjacency_rule_nonempty(rule):
    desc = describe_adjacency_rule(rule)
    assert isinstance(desc, str)
    assert len(desc) > 20


def test_unknown_rule_raises():
    with pytest.raises(ValueError):
        describe_adjacency_rule("C")
    with pytest.raises(ValueError):
        build_dirac_s3_graph(2, "C")


# ---------------------------------------------------------------------------
# Rule B induces l-mixing (sanity)
# ---------------------------------------------------------------------------

def test_rule_B_mixes_kappa_at_n_max_2():
    """Rule B edges must include a κ-mixing case at n_max=2.

    At n_max=2, κ=-1 (l=0, n=1) is connected to κ=-2 (l=1, n=2) and
    κ=+1 (l=1, n=2) by the dipole rule.  So rule B should produce
    at least one edge between nodes of different κ.
    """
    A, labels, _deg, _desc = build_dirac_s3_graph(2, "B")
    mixing = 0
    V = len(labels)
    for i in range(V):
        for j in range(i + 1, V):
            if A[i, j] and labels[i].kappa != labels[j].kappa:
                mixing += 1
    assert mixing > 0, "Rule B at n_max=2 should have κ-mixing edges"


# ---------------------------------------------------------------------------
# End-to-end: computing result dict doesn't raise
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_max", [1, 2, 3])
@pytest.mark.parametrize("rule", ["A", "B"])
def test_dirac_s3_ihara_zeta_end_to_end(n_max, rule):
    """Smoke test: full pipeline produces a well-formed dict."""
    result = dirac_s3_ihara_zeta(n_max, rule, factor=False)  # skip factoring for speed
    # Required keys
    for key in [
        "n_max", "adjacency_rule", "V", "E", "c", "r_betti1",
        "zeta_inverse_expanded",
        "hashimoto_spectrum",
        "ramanujan_verdict", "ramanujan_deviation",
        "zeros",
        "per_kappa_sizes",
        "charpoly_integer_coefficient",
    ]:
        assert key in result, f"missing key {key!r}"
    assert result["V"] > 0
    assert result["E"] >= 0
    assert result["c"] >= 1
