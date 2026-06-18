"""
Tests for geovac.ihara_zeta — Ihara zeta, Bass formula, Hashimoto
edge operator, and graph-Ramanujan property.

Sanity witnesses:
  - K_4 has the known closed form
        zeta_{K_4}(s)^{-1} = (1-s)^3 (1+s)^2 (1-2s) (2s^2 + s + 1)^3
    (Terras 2011; also Stark-Terras 1996).
  - K_{3,3} is 3-regular Ramanujan; non-trivial Hashimoto eigenvalues
    lie exactly on |mu| = sqrt(2).
  - Bass determinantal formula and the truncated Euler product
    agree to order s^6.

GeoVac graphs:
  - S^3 Coulomb graph (GeometricLattice) at max_n = 2 and 3.
  - S^5 Bargmann-Segal graph (build_bargmann_graph) at N_max = 2 and 3.
"""

import numpy as np
import pytest
import sympy as sp

from geovac.ihara_zeta import (
    functional_equation_report,
    hashimoto_matrix,
    ihara_zeta_bass,
    ihara_zeta_euler,
    is_ramanujan,
    zeta_zeros_from_hashimoto,
)

s = sp.symbols("s")


# ---------------------------------------------------------------------------
# K_4 sanity witness
# ---------------------------------------------------------------------------

def _K_n(n: int) -> np.ndarray:
    A = np.ones((n, n), dtype=int) - np.eye(n, dtype=int)
    return A


def test_K4_bass_closed_form():
    """K_4 matches the Terras/Stark closed form."""
    A = _K_n(4)
    zeta = ihara_zeta_bass(A)
    zeta_inv = sp.expand(1 / zeta)
    expected = (1 - s) ** 3 * (1 + s) ** 2 * (1 - 2 * s) * (2 * s ** 2 + s + 1) ** 3
    # sign convention: Terras gives (1 - s)^3 etc.; our sign matches up to factor (-1)^(nothing).
    diff = sp.simplify(zeta_inv - expected)
    assert diff == 0, f"K_4 Bass differs from Terras form: {diff}"


def test_K4_bass_vs_euler_agree_to_s6():
    """Bass and truncated Euler agree up to s^6 on K_4."""
    A = _K_n(4)
    zb = ihara_zeta_bass(A)
    ze = ihara_zeta_euler(A, max_length=6)
    series_bass = sp.expand(sp.series(zb, s, 0, 7).removeO())
    series_euler = sp.expand(sp.series(ze, s, 0, 7).removeO())
    assert sp.simplify(series_bass - series_euler) == 0


def test_K4_is_ramanujan_zero_deviation():
    """K_4 is 3-regular Ramanujan: nontriv |mu| = sqrt(2)."""
    A = _K_n(4)
    ir, dev, expl = is_ramanujan(A)
    assert ir
    assert abs(dev) < 1e-9, expl


# ---------------------------------------------------------------------------
# K_{3,3} second witness (3-regular bipartite)
# ---------------------------------------------------------------------------

def _K33() -> np.ndarray:
    A = np.zeros((6, 6), dtype=int)
    for i in range(3):
        for j in range(3, 6):
            A[i, j] = 1
            A[j, i] = 1
    return A


def test_K33_bass_vs_euler_agree_to_s6():
    A = _K33()
    zb = ihara_zeta_bass(A)
    ze = ihara_zeta_euler(A, max_length=6)
    series_bass = sp.expand(sp.series(zb, s, 0, 7).removeO())
    series_euler = sp.expand(sp.series(ze, s, 0, 7).removeO())
    assert sp.simplify(series_bass - series_euler) == 0


def test_K33_is_ramanujan_zero_deviation():
    A = _K33()
    ir, dev, expl = is_ramanujan(A)
    assert ir
    assert abs(dev) < 1e-9, expl


# ---------------------------------------------------------------------------
# Hashimoto trace = closed non-backtracking walks (sanity)
# ---------------------------------------------------------------------------

def test_hashimoto_trace_on_K4():
    """tr(T^2) counts oriented closed NB walks of length 2.  On K_n,
    there are none (cannot return without backtracking in 2 edges), so
    tr(T^2) = 0.  tr(T^3) counts closed NB walks of length 3, i.e. the
    oriented triangles of K_n: 2 * C(n, 3) * 3! / 3 = 6 * 4 = 24 for K_4
    (each of 4 triangles traversed clockwise and counterclockwise from
    each of 3 starting edges).
    """
    A = _K_n(4)
    T = hashimoto_matrix(A).astype(float)
    assert abs(T.trace()) < 1e-9
    T2 = T @ T
    assert abs(T2.trace()) < 1e-9  # no 2-cycles (simple graph)
    T3 = T2 @ T
    # K_4 has C(4,3)=4 triangles, each giving 2 orientations * 3 starting
    # oriented edges = 24.
    assert abs(T3.trace() - 24.0) < 1e-9


# ---------------------------------------------------------------------------
# GeoVac graphs
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def s3_graphs():
    """S^3 Coulomb graph at max_n = 2 and 3."""
    from geovac.lattice import GeometricLattice
    out = {}
    for mn in (2, 3):
        L = GeometricLattice(max_n=mn)
        out[mn] = L.adjacency.toarray()
    return out


@pytest.fixture(scope="module")
def s5_graphs():
    """S^5 Bargmann-Segal graph at N_max = 2 and 3."""
    from geovac.nuclear.bargmann_graph import build_bargmann_graph
    out = {}
    for nm in (2, 3):
        g = build_bargmann_graph(nm)
        out[nm] = g.adjacency_dense()
    return out


def test_s3_max_n_2_is_disconnected_forest_ihara_is_trivial(s3_graphs):
    """S^3 max_n=2 is V=5, E=3, disconnected (3 components per l).  For
    a forest (r = 0), zeta_G(s) is the trivial rational function (no
    non-trivial zeros — the graph has no cycles).
    """
    A = (s3_graphs[2] != 0).astype(int)
    assert A.shape == (5, 5)
    from geovac.ihara_zeta import _count_components
    # Components coincide with l-shells.
    c = _count_components(A)
    assert c == 2  # l=0 path of 2 nodes, l=1 star of 3 leaves on center
    # Bass should compute without error.
    rep = functional_equation_report(A, num_components=c)
    assert rep["r_betti1"] == 0
    zeros = zeta_zeros_from_hashimoto(A)
    # A forest has no non-backtracking closed walks and T is nilpotent,
    # so zeta_G has no zeros.
    assert len(zeros) == 0


def test_s3_max_n_3_ramanujan(s3_graphs):
    """S^3 max_n=3 is Ramanujan (deviation is negative)."""
    A = (s3_graphs[3] != 0).astype(int)
    ir, dev, expl = is_ramanujan(A)
    assert ir, expl


def test_s5_N_max_2_ramanujan(s5_graphs):
    A = (s5_graphs[2] != 0).astype(int)
    ir, dev, expl = is_ramanujan(A)
    assert ir, expl


def test_s5_N_max_3_ramanujan(s5_graphs):
    A = (s5_graphs[3] != 0).astype(int)
    ir, dev, expl = is_ramanujan(A)
    assert ir, expl


def test_s3_max_n_3_bass_vs_euler_small_s(s3_graphs):
    """Bass and truncated Euler agree on the first few s^L coefficients
    of the S^3 max_n=3 graph."""
    A = (s3_graphs[3] != 0).astype(int)
    zb = ihara_zeta_bass(A)
    ze = ihara_zeta_euler(A, max_length=6)
    series_bass = sp.expand(sp.series(zb, s, 0, 7).removeO())
    series_euler = sp.expand(sp.series(ze, s, 0, 7).removeO())
    assert sp.simplify(series_bass - series_euler) == 0


def test_s5_N_max_2_bass_vs_euler_small_s(s5_graphs):
    A = (s5_graphs[2] != 0).astype(int)
    zb = ihara_zeta_bass(A)
    ze = ihara_zeta_euler(A, max_length=6)
    series_bass = sp.expand(sp.series(zb, s, 0, 7).removeO())
    series_euler = sp.expand(sp.series(ze, s, 0, 7).removeO())
    assert sp.simplify(series_bass - series_euler) == 0


def test_zeta_evaluation_at_sample_points(s3_graphs, s5_graphs):
    """Ihara zeta must evaluate to a finite non-zero value at the
    critical circle |s| = 1/sqrt(q_max) for each graph.  This is a
    smoke test that the rational function is well-formed.
    """
    for label, graphs in [("s3", s3_graphs), ("s5", s5_graphs)]:
        for key, A in graphs.items():
            A01 = (A != 0).astype(int)
            d = A01.sum(axis=1)
            q_max = max(int(d.max()) - 1, 1)
            rc = 1.0 / np.sqrt(q_max)
            zeta = ihara_zeta_bass(A01)
            val = complex(zeta.subs(s, rc / 2))  # safely inside convergence disk
            assert np.isfinite(val.real) and np.isfinite(val.imag), (
                f"{label} {key} zeta evaluation failed at s={rc/2:.4f}"
            )


def test_functional_equation_regular_identifies_q():
    """On K_4 (3-regular), the functional-equation report picks out q=2
    with the correct critical radius 1/sqrt(2)."""
    rep = functional_equation_report(_K_n(4))
    assert rep["regular"] is True
    assert rep["q_effective"] == 2
    assert abs(rep["critical_radius"] - 1 / np.sqrt(2)) < 1e-12


def test_functional_equation_irregular_flags_graph_rh_regime(s5_graphs):
    """Irregular graph: report must flag q_max vs q_min asymmetry."""
    A = (s5_graphs[3] != 0).astype(int)
    rep = functional_equation_report(A)
    assert rep["regular"] is False
    assert "q_max" in rep and "q_min" in rep
    assert rep["q_max"] > rep["q_min"]


@pytest.mark.slow
def test_s5_N3_zeta_matches_paper_closed_form(s5_graphs):
    """End-to-end (Paper 29, sec:closed_forms): the S^5 N_max=3 Ihara zeta
    inverse factors EXACTLY as the published (s^2-1)^23 * P12(s) * P22(s)
    (total degree 80), not merely to series order s^6.  Closes the closed-form
    coverage gap (a typo in the published P12/P22 would now be caught) and pins
    the corrected degree bookkeeping: the Hashimoto matrix has size 2|E|=84;
    the 4 missing degrees are zeros at s=infinity from the degree-1 leaves."""
    s = sp.symbols("s")
    zeta = ihara_zeta_bass(s5_graphs[3], symbolic=True)
    _, zinv = sp.fraction(sp.together(zeta))          # zeta = 1 / zinv
    zinv = sp.expand(zinv)

    P12 = 432*s**12 + 666*s**10 + 374*s**8 + 135*s**6 + 47*s**4 + 11*s**2 + 1
    P22 = (829440*s**22 + 2453184*s**20 + 3308104*s**18 + 2696682*s**16
           + 1470640*s**14 + 557227*s**12 + 146654*s**10 + 25709*s**8
           + 2630*s**6 + 79*s**4 - 12*s**2 - 1)

    # Published P12, P22 are EXACT factors of the graph-derived inverse zeta:
    assert sp.rem(zinv, P12, s) == 0, "published P12 is not an exact factor"
    assert sp.rem(zinv, P22, s) == 0, "published P22 is not an exact factor"
    # Total degree is 80 (NOT the Hashimoto matrix size 84):
    assert sp.Poly(zinv, s).degree() == 80
    # Cofactor is (s^2-1)^23 up to overall sign (the 46 trivial zeros s=+-1):
    cof = sp.cancel(zinv / (P12 * P22))
    assert (sp.expand(cof - (s**2 - 1)**23) == 0
            or sp.expand(cof + (s**2 - 1)**23) == 0), "cofactor != (s^2-1)^23"
