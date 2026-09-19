r"""Guard: the per-rate-pair V_ee X-table (Paper 12, two-block radial exponents).

WHY THIS FILE EXISTS.  The X-table generalization was validated in
``debug/multiexp_vee_xtable.py`` at five levels, but ``debug/`` is prunable by the
§9 clean-room policy, so that validation had no permanent home -- logged as a
COVERAGE GAP in ``docs/claim_test_matrix.md``.  This closes it.

WHAT WRONG ANSWERS THIS REJECTS (the §9 question for a guard is not "does it pass"):

1. **Reverting to a single rate.**  If the two electrons' rates are collapsed back
   to a shared ``2*alpha`` -- the pre-2026-09-18 behaviour, and the most likely
   regression since it is what every other call site still does -- then a table
   built at c1 != c2 would silently equal the c1 == c2 one.
   ``test_two_rate_table_differs_from_either_single_rate`` forbids that.
2. **Swapping which side carries which rate.**  The assignment (A/P side <- c1,
   B/Q side <- c2, IBP divisor = the P-side rate, tail table at c1+c2) is forced by
   the production code's structure, and getting it backwards still produces a
   plausible-looking symmetric-ish table.  ``test_two_rate_entry_matches_quadrature``
   rejects it against a reference built from the DEFINITION.
3. **Losing the transpose identity.**  ``test_mirror_rate_pair_is_the_transpose``
   pins X^{(c1,c2)}[P1][P2] == X^{(c2,c1)}[P2][P1], which is what makes 9 ordered
   rate pairs cost 6 builds; a "simplification" that rebuilds all 9 independently
   would pass every other test here while silently costing 50% more.
   **Its limitation, from the fire test rather than from speculation:** this test
   PASSES under the collapsed-rate mutation of (1), because a collapsed table is
   symmetric and the identity then holds trivially.  It guards the build-count
   optimisation, NOT the rate split.  Do not read it as collapse protection.
4. **Breaking the degenerate reduction.**  ``test_degenerate_reduces_to_production``
   and ``test_degenerate_energy_is_unchanged`` pin that a single-rate basis still
   reproduces ``pr._build_Xtab_mp`` and ``recondition_energy``'s energy.  The
   energy test is NOT redundant with the matrix tests: in the v5.14.0 sprint four
   matrix-level falsifiers passed while a fixture solved the wrong Hamiltonian.

WHY THE QUADRATURE REFERENCE IS THE LOAD-BEARING ONE.  ``debug/direct_vee_corr.py``
validated the ordered integral to 1e-27, but its ``_corr_gen`` takes a SINGLE ``c``
for both sides, so it cannot reach c1 != c2.  Test 2 therefore builds its reference
by direct 2D mpmath quadrature of the Neumann kernel -- no monomial moments, no
B-table, no IBP -- which is independent of everything the table uses.

A TRAP THIS FILE ENCODES so it is not re-discovered: ``pr.vee_mp`` does NOT
substitute a default at ``l_neumann <= 0`` (it only does
``min(l_neumann, q_max + 2*max(s))``, so 0 stays 0), while ``recondition_energy``
substitutes ``2*l_max + 4*mu_max + 10``.  Comparing the two with 0 compares
different physics: it produced a 0.155-relative "failure" that was pure harness
error.  Every test here passes a positive truncation explicitly.
"""
from __future__ import annotations

import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from geovac import neumann_vee_general_m as ngm
from geovac import prolate_recondition as pr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "debug"))
multiexp_vee_xtable = pytest.importorskip("multiexp_vee_xtable")

DPS = 40
ALPHA = 1.40
# a genuinely split pair of rates; c1 != c2 is the whole point
C1, C2 = 2.80, 2.00


def _fns(j_max: int, l_max: int, mu_max: int, alpha: float = ALPHA):
    idx = pr._product_index(j_max, l_max, mu_max)
    return [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]


class TestRateAssignment:
    """The two-rate table against a reference with no shared machinery."""

    @pytest.mark.parametrize("l,m,s,P1,P2", [
        (0, 0, 0, 0, 0),
        (0, 0, 0, 1, 0),      # asymmetric powers: exercises I1 != I2
        (0, 0, 0, 1, 2),
        (1, 0, 0, 0, 1),
        (1, 1, 1, 0, 0),      # associated-Legendre sector, not just sigma
    ])
    def test_two_rate_entry_matches_quadrature(self, l, m, s, P1, P2):
        """X at c1 != c2 == direct 2D quadrature of the Neumann kernel.

        Rejects a swapped rate assignment.  The reference integrates the
        definition (P_l on the smaller argument, Q_l on the larger) with no
        monomial moments, no B-table and no IBP.
        """
        mp.mp.dps = DPS
        p_max = max(P1, P2)
        X = multiexp_vee_xtable.build_Xtab_rated(
            [(m, s)], l, p_max, [C1, C2], {(m, s): l})
        got = X[(C1, C2, l, m, s)][P1][P2]
        ref = multiexp_vee_xtable._x_entry_quad(l, m, s, P1, P2, C1, C2)
        assert ref != 0, "degenerate reference; pick a different entry"
        assert abs(got - ref) / abs(ref) < mp.mpf('1e-20')

    def test_two_rate_table_differs_from_either_single_rate(self):
        """c1 != c2 must NOT reproduce the c1 == c2 table.

        This is the anti-regression: collapsing back to one shared rate is what
        every other call site still does, so it is the likeliest reversion, and
        it would leave every *degenerate* test here passing.
        """
        mp.mp.dps = DPS
        m, s, l, p_max = 0, 0, 1, 2
        caps = {(m, s): l}
        mixed = multiexp_vee_xtable.build_Xtab_rated(
            [(m, s)], l, p_max, [C1, C2], caps)[(C1, C2, l, m, s)]
        for c in (C1, C2):
            same = multiexp_vee_xtable.build_Xtab_rated(
                [(m, s)], l, p_max, [c], caps)[(c, c, l, m, s)]
            diffs = [abs(mixed[a][b] - same[a][b])
                     for a in range(p_max + 1) for b in range(p_max + 1)]
            scale = max(abs(same[a][b])
                        for a in range(p_max + 1) for b in range(p_max + 1))
            assert max(diffs) > mp.mpf('1e-6') * scale, (
                f"two-rate table collapsed onto the single-rate table at c={c}")

    def test_mirror_rate_pair_is_the_transpose(self):
        """X^{(c1,c2)}[P1][P2] == X^{(c2,c1)}[P2][P1].

        Pins the identity that makes 9 ordered rate pairs cost 6 builds.  Uses
        OFF-DIAGONAL powers only: at P1 == P2 the identity is trivially true and
        would not detect its loss.
        """
        mp.mp.dps = DPS
        m, s, l, p_max = 0, 0, 1, 2
        X = multiexp_vee_xtable.build_Xtab_rated(
            [(m, s)], l, p_max, [C1, C2], {(m, s): l})
        fwd, mir = X[(C1, C2, l, m, s)], X[(C2, C1, l, m, s)]
        checked = 0
        for P1 in range(p_max + 1):
            for P2 in range(p_max + 1):
                if P1 == P2:
                    continue
                den = abs(fwd[P1][P2]) or mp.mpf(1)
                assert abs(fwd[P1][P2] - mir[P2][P1]) / den < mp.mpf('1e-30')
                checked += 1
        assert checked > 0, "no off-diagonal entries were compared"


class TestDegenerateReduction:
    """A single-rate basis must reproduce production exactly."""

    def test_degenerate_reduces_to_production(self):
        """build_Xtab_rated at one rate == pr._build_Xtab_mp, elementwise."""
        mp.mp.dps = DPS
        fns = _fns(2, 2, 1)
        ms_pairs, p_max, _q, l_neu, l_caps = multiexp_vee_xtable._vee_shapes(fns, 0)
        ref = pr._build_Xtab_mp(ms_pairs, l_neu, p_max, ALPHA, l_caps)
        got = multiexp_vee_xtable.build_Xtab_rated(
            ms_pairs, l_neu, p_max, [2.0 * ALPHA], l_caps)
        assert ref, "no blocks built; the shape helper drifted"
        c = 2.0 * ALPHA
        worst = max(abs(ref[(l, m, s)][a][b] - got[(c, c, l, m, s)][a][b])
                    for (l, m, s), mat in ref.items()
                    for a in range(len(mat)) for b in range(len(mat)))
        scale = max(abs(v) for mat in ref.values() for row in mat for v in row)
        assert worst / scale < mp.mpf('1e-30')

    def test_degenerate_vee_reduces_to_production(self):
        """vee_multi at one rate == pr.vee_mp, elementwise.

        Note the EXPLICIT positive l_neumann: passing 0 compares different
        Neumann truncations (see the module docstring's trap).
        """
        mp.mp.dps = DPS
        fns = _fns(2, 2, 0)
        l_neu = 8
        ref = pr.vee_mp(fns, ALPHA, pr.R_DEFAULT, l_neu)
        got = multiexp_vee_xtable.vee_multi(fns, ALPHA, pr.R_DEFAULT, l_neu)
        n = len(fns)
        assert all(got[i, j] is not None for i in range(n) for j in range(n)), \
            "gather left entries unfilled"
        worst = max(abs(ref[i, j] - got[i, j])
                    for i in range(n) for j in range(n))
        scale = max(abs(ref[i, j]) for i in range(n) for j in range(n))
        assert worst / scale < mp.mpf('1e-30')

    def test_multirate_gather_fills_every_entry(self):
        """A genuine split must leave no entry unrouted.

        The degenerate tests cannot catch this: with one rate pair the gather
        covers the matrix however badly it selects.
        """
        mp.mp.dps = DPS
        poc = pytest.importorskip("multiexp_overlap_poc")
        alpha_of = poc.alpha_of_split(1, 1.60, 1.00)
        fns = _fns(2, 2, 0, alpha=1.60)
        rates = multiexp_vee_xtable.pair_rates(fns, alpha_of)
        assert len(rates) == 3, f"expected 3 pair rates at two blocks, got {rates}"
        V = multiexp_vee_xtable.vee_multi(
            fns, 1.60, pr.R_DEFAULT, 8, alpha_of=alpha_of)
        n = len(fns)
        assert all(V[i, j] is not None for i in range(n) for j in range(n))
        ref1 = multiexp_vee_xtable.vee_multi(fns, 1.60, pr.R_DEFAULT, 8)
        same = sum(1 for i in range(n) for j in range(n) if V[i, j] == ref1[i, j])
        assert same < n * n, "split alpha_of produced the single-rate matrix"


@pytest.mark.slow
def test_degenerate_energy_is_unchanged():
    """END-TO-END: vee_multi degenerate gives recondition_energy's own energy.

    NOT redundant with the matrix tests above.  In the v5.14.0 sprint four
    matrix-level falsifiers passed while a fixture solved the wrong Hamiltonian:
    a matrix check proves the matrix, only the energy proves the wiring.
    """
    j_max, l_max, mu_max = 3, 3, 1
    R = pr.R_DEFAULT
    with mp.workdps(DPS):
        fns = _fns(j_max, l_max, mu_max)
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1)
                  if (l + m) % 2 == 0])
        Nmu = mu_max + 1
        l_neu = 2 * l_max + 4 * mu_max + 10
        Tr, Ta = pr._transforms_per_mu("laguerre_legendre", j_max, l_max,
                                       mu_max, ALPHA)
        S_o, H1_o = pr.build_one_body_direct(j_max, l_max, mu_max, ALPHA, R,
                                             "laguerre_legendre", DPS)
        energies = []
        for fn in (pr.vee_mp, multiexp_vee_xtable.vee_multi):
            V = fn(fns, ALPHA, R, l_neu)
            V_o = pr._to_f64(pr._factored_cob(V, Nmu, Nr, Na, Tr, Ta))
            E, _cn, _nk, _sw = pr._normalized_solve(S_o, H1_o + V_o + (1.0 / R) * S_o)
            energies.append(E)
    assert abs(energies[0] - energies[1]) < 1e-12, (
        f"degenerate vee_multi changed the energy: {energies}")
    # and it must be the KNOWN energy, not merely self-consistent
    assert abs(energies[1] - (-1.1726082231)) < 1e-8
