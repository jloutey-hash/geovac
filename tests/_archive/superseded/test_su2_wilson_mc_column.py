"""ARCHIVED 2026-11-07 (superseded): the Monte Carlo Wilson-loop column of old Paper 30.

The values 0.126 / 0.158 / 0.332 / 0.681 came from an unconverged chain and were replaced by the
exact <W> = (I2/I1)^A (Paper 30 prop:factorization); no live paper carries them, and
tests/test_su2_wilson_gauge.py::test_single_plaquette_expectation_is_I2_over_I1_exactly rejects
them. Kept as institutional memory (CLAUDE.md sec. 14). Also the only coverage of
geovac.su2_wilson_gauge.monte_carlo_wilson_expectation, which no live claim uses.
"""
import numpy as np
import pytest

from geovac.lattice import GeometricLattice
from geovac.su2_wilson_gauge import *  # noqa: F401,F403


class TestMonteCarloWilsonLoop:
    """Pins Table 'Wilson-loop expectation' of Paper 30 (Result 3):
    <W>_MC = 0.126 / 0.158 / 0.332 / 0.681 at beta = 0.5 / 1 / 2 / 5,
    n_max = 3, loop = first p-block 4-cycle, 2000 samples after 500
    thermalization steps, seed = 12.

    The Metropolis chain is fully deterministic for a fixed seed
    (verified: two consecutive runs agree bit-exactly), so the pins use
    atol = 1e-6 (guard against libm variation across platforms; the
    paper's numbers are quoted to 3 decimals).
    """

    def _setup(self):
        from geovac.su2_wilson_gauge import monte_carlo_wilson_expectation
        lat = GeometricLattice(max_n=3)
        A = lat.adjacency.toarray()
        oriented, _ = enumerate_oriented_edges(A)
        forward = [(e.source, e.target)
                   for e in oriented if e.source < e.target]
        edge_list = [OrientedEdge(u, v) for (u, v) in forward]
        plaqs = enumerate_plaquettes(A, max_length=4, both_orientations=False)
        assert len(forward) == 13   # "13 forward links" (paper)
        assert len(plaqs) == 2
        # The loop is the first primitive 4-cycle in the p-block:
        # (2,1,-1) -> (2,1,0) -> (3,1,0) -> (3,1,-1) -> close.
        loop_states = [lat.states[e.source] for e in plaqs[0]]
        assert loop_states == [(2, 1, -1), (2, 1, 0), (3, 1, 0), (3, 1, -1)]
        return monte_carlo_wilson_expectation, plaqs, edge_list

    @pytest.mark.slow
    def test_mc_wilson_table_seed12(self):
        mc, plaqs, edge_list = self._setup()
        # (beta, mean, stderr) recomputed live 2026-07-02; means round to
        # the paper table 0.126 / 0.158 / 0.332 / 0.681.
        expected = {
            0.5: (0.12576712406708293, 0.011194572527873982),
            1.0: (0.1575325588927505, 0.010766717015609136),
            2.0: (0.3324166614690038, 0.009133638708602723),
            5.0: (0.6812474363347273, 0.00603985232651571),
        }
        paper_table = {0.5: 0.126, 1.0: 0.158, 2.0: 0.332, 5.0: 0.681}
        means = {}
        for beta, (m_exp, s_exp) in expected.items():
            m, s = mc(plaqs[0], plaqs, edge_list, beta,
                      n_samples=2000, n_thermalize=500, seed=12)
            means[beta] = m
            assert abs(m - m_exp) < 1e-6, (beta, m, m_exp)
            assert abs(s - s_exp) < 1e-6, (beta, s, s_exp)
            assert abs(round(m, 3) - paper_table[beta]) < 1e-12
        # Monotonic in beta (strong -> weak coupling)
        bs = sorted(means)
        assert all(means[bs[i]] < means[bs[i + 1]] for i in range(len(bs) - 1))

    @pytest.mark.slow
    def test_mc_seed_determinism(self):
        """Two runs with the same seed agree bit-exactly."""
        mc, plaqs, edge_list = self._setup()
        m1, s1 = mc(plaqs[0], plaqs, edge_list, 0.5,
                    n_samples=2000, n_thermalize=500, seed=12)
        m2, s2 = mc(plaqs[0], plaqs, edge_list, 0.5,
                    n_samples=2000, n_thermalize=500, seed=12)
        assert m1 == m2
        assert s1 == s2
