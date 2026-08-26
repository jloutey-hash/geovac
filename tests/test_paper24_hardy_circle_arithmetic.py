"""Paper 24 -- arithmetic of the Bargmann/S^5 Hardy phase circle (Layer 7 +
companion chi_{-4} NON-layer; CHANGELOG v4.103.0, driver
debug/aha_t4b_bargmann_circle.py).

Pins: (i) the Hardy generator N + 3/2 is strictly positive on the production
graph (no +/- pair, hence no traceless 2-plane, hence no Hodge circle);
(ii) the graph is bipartite (all edges |dN| = 1), so (-1)^N is the natural Z2
grading; (iii) the quarter point acts by i^N: Q-rational exactly where it
squares to +1; (iv) ungraded partition function = 1/(1-q)^3 (Hilbert series of
C[z1,z2,z3]: pure Tate, conductor 1); (v) the Z2-graded spectral zeta is
Z_-(s) = 2^(s-3) (beta(s) - beta(s-2)) -- exactly 1/4 of Paper 28 thm:chi4's
Dirac value (same index + d/2 eigenvalues at d = 3; degeneracies differ by 4)
-- via eta_H(s) = 2^s (1 - beta(s)); (vi) an integer zero-point shift gives a
pure-zeta (Tate) value instead: the chi_{-4} needs the half-integer
(metaplectic) shift.
"""
from fractions import Fraction

import mpmath as mp
import sympy as sp

from geovac.nuclear.bargmann_graph import build_bargmann_graph, shell_degeneracy


def _beta(s):
    """Dirichlet beta(s) = sum_k (-1)^k (2k+1)^(-s)."""
    return mp.nsum(lambda k: (-1) ** int(k) * (2 * k + 1) ** (-s), [0, mp.inf],
                   method="a")


def test_generator_positive_and_bipartite():
    g = build_bargmann_graph(5)
    assert all(d == Fraction(2 * N + 3, 2) for d, (N, _, _) in zip(g.diagonal, g.nodes))
    assert all(d > 0 for d in g.diagonal)                      # strictly positive
    assert len(g.adjacency) > 0
    for (i, j) in g.adjacency:
        assert abs(g.nodes[i][0] - g.nodes[j][0]) == 1         # bipartite in N


def test_quarter_point_rationality():
    """i^N is Q-rational exactly on even shells, where it squares to +1."""
    for N in range(8):
        z = 1j ** N
        rational = z.imag == 0.0
        squares_to_minus_one = (z * z) == -1
        assert rational == (N % 2 == 0)
        assert squares_to_minus_one == (N % 2 == 1)


def test_partition_function_pure_tate():
    """sum_N g_N q^N = 1/(1-q)^3: Hilbert series of C[z1,z2,z3], conductor 1."""
    q = sp.Symbol("q")
    lhs = sum(shell_degeneracy(N) * q ** N for N in range(21))
    rhs = sp.series(1 / (1 - q) ** 3, q, 0, 21).removeO()
    assert sp.expand(lhs - rhs) == 0


def test_eta_h_identity():
    """eta_H(s) = sum_N (-1)^N (N + 3/2)^(-s) = 2^s (1 - beta(s))."""
    mp.mp.dps = 40
    for s in (4, 6):
        eta_h = mp.nsum(lambda N: (-1) ** int(N) * (N + mp.mpf(3) / 2) ** (-s),
                        [0, mp.inf], method="a")
        assert abs(eta_h - 2 ** s * (1 - _beta(s))) < mp.mpf("1e-30")


def test_graded_zeta_quarter_of_dirac():
    """Z_-(s) = 2^(s-3)(beta(s) - beta(s-2)); Paper 28 thm:chi4 is
    2^(s-1)(beta(s) - beta(s-2)) -- ratio exactly 4."""
    mp.mp.dps = 40
    for s in (4, 6):
        direct = mp.nsum(lambda N: (-1) ** int(N) * shell_degeneracy(int(N))
                         * (N + mp.mpf(3) / 2) ** (-s), [0, mp.inf], method="a")
        closed = 2 ** (s - 3) * (_beta(s) - _beta(s - 2))
        assert abs(direct - closed) < mp.mpf("1e-28")
        dirac = 2 ** (s - 1) * (_beta(s) - _beta(s - 2))
        assert abs(dirac / closed - 4) < mp.mpf("1e-28")
    # headline value: Z_-(4) = 2 (beta(4) - G), G = Catalan's constant
    z4 = 2 * (_beta(4) - mp.catalan)
    assert abs(z4 - mp.mpf("0.14595791512777264")) < mp.mpf("1e-15")


def test_integer_shift_is_pure_zeta():
    """Replacing the metaplectic 3/2 by an integer shift gives eta(s) =
    (1 - 2^(1-s)) zeta(s): pure Tate.  The chi_{-4} needs the half-integer
    shift."""
    mp.mp.dps = 30
    for s in (4, 6):
        eta = mp.nsum(lambda N: (-1) ** int(N) * (N + 1) ** (-s), [0, mp.inf],
                      method="a")
        assert abs(eta - (1 - 2 ** (1 - s)) * mp.zeta(s)) < mp.mpf("1e-25")


def test_graded_zeta_dps_stable():
    """Cross-precision stability guard (the corpus's accelerated-sum lesson):
    the alternating sums agree across working precisions."""
    vals = {}
    for dps in (40, 60):
        mp.mp.dps = dps
        vals[dps] = mp.nsum(lambda N: (-1) ** int(N) * shell_degeneracy(int(N))
                            * (N + mp.mpf(3) / 2) ** (-4), [0, mp.inf], method="a")
    mp.mp.dps = 60
    assert abs(vals[40] - vals[60]) < mp.mpf("1e-30")
