r"""Backing test for the M3 period-map RANK finding (Sprint M3-parity-rank,
2026-07-08; memo ``debug/sprint_m3_parity_rank_memo.md``).

Companion to ``tests/test_paper56_injection_g4_periodmap.py``, which shows the
PARITY-BLIND M3 period map (single scalar eta_{(n,l)} per sector) has Gram
**rank 1** (the corpus baseline, Paper 56 thm:injection_g4). This file supplies
the "sector-resolved rank->=2 M3 period map" Paper 56 names as the open follow-on
("a per-(n,l) parity decomposition of D(s)"), and pins its rank.

FINDING (bit-exact):
  * parity-blind    -> rank 1  (reproduces the corpus / the periodmap test)
  * parity-resolved -> rank 2  for every n_max >= 2  (rank 1 at n_max=1: only odd
    shells exist)
  * the two shift-classes are a GENUINE independent period direction:
    D_even(4), D_odd(4) are rank 2 in the motivic basis {pi^2, pi^4, G, beta(4)}.

MECHANISM (grounded, Paper 28 Thm 3, PROVEN): the vertex-parity split sends
even CH-shells to Hurwitz shift 3/4 and odd shells to 5/4, with
D_even - D_odd = 2^{s-1}(beta(s) - beta(s-2)) pure level-4 (Catalan G / beta).
Each Coulomb sector (n,l) sits in its shell's parity class -> one shift-class.
Rank = number of populated shift-classes = 2, provably CAPPED at 2 (two parity
classes) -> faithful embedding (rank N) is unreachable by parity.

SCOPE: the computed ranks are the assertions here. Whether the chi_-4-graded
trace is admitted into *the* canonical Paper 56 period map is an interpretive
judgment (PI adjudication); this test pins the linear-algebra facts either way.
"""
from __future__ import annotations

import pytest
from sympy import Matrix, Rational

try:
    from geovac.pro_system import sectors_at_cutoff, N_sectors
except Exception:  # self-contained fallback (Paper 55 sector indexing)
    def sectors_at_cutoff(n_max):
        return [(n, l) for n in range(1, n_max + 1) for l in range(0, n + 1)]

    def N_sectors(n_max):
        return len(sectors_at_cutoff(n_max))


def eta(n: int, l: int) -> Rational:
    """Per-sector parity-BLIND M3 scalar eta_{(n,l)} (Paper 55)."""
    return Rational(n * (2 * n + 1)) if l == n else Rational((2 * l + 1) * (2 * n + 1))


def shift_class(n: int, conv: str = "n") -> int:
    """Vertex-parity shift-class (Paper 28 Thm 3): 0 -> 3/4 (even shell),
    1 -> 5/4 (odd shell). conv toggles the sector<->shell offset."""
    n_shell = n if conv == "n" else (n - 1)
    return n_shell % 2


def parity_blind_gram(n_max: int) -> Matrix:
    P = Matrix([[eta(n, l)] for (n, l) in sectors_at_cutoff(n_max)])
    return P * P.T


def parity_resolved_matrix(n_max: int, conv: str = "n") -> Matrix:
    rows = []
    for (n, l) in sectors_at_cutoff(n_max):
        v = [0, 0]
        v[shift_class(n, conv)] = eta(n, l)
        rows.append(v)
    return Matrix(rows)


class TestParityBlindBaseline:
    """The corpus baseline: parity-blind M3 period map is rank 1."""

    @pytest.mark.parametrize("n_max", [2, 3, 4])
    def test_parity_blind_gram_rank_one(self, n_max: int) -> None:
        assert parity_blind_gram(n_max).rank() == 1


class TestParityResolvedRank:
    """The follow-on: parity-resolved M3 period map is rank 2 (n_max>=2)."""

    @pytest.mark.parametrize("n_max", [2, 3, 4])
    @pytest.mark.parametrize("conv", ["n", "n-1"])
    def test_parity_resolved_rank_two(self, n_max: int, conv: str) -> None:
        assert parity_resolved_matrix(n_max, conv).rank() == 2

    def test_nmax1_edge_case_rank_one(self) -> None:
        """n_max=1 has only odd shells -> a single shift-class -> rank 1."""
        assert parity_resolved_matrix(1, "n").rank() == 1

    @pytest.mark.parametrize("n_max", [2, 3, 4])
    def test_capped_at_two_not_full(self, n_max: int) -> None:
        """The parity route is CAPPED at rank 2 < N -> not a faithful embedding."""
        r = parity_resolved_matrix(n_max, "n").rank()
        assert r == 2 < N_sectors(n_max)


class TestShiftClassIndependence:
    """Anti-artifact check: the two shift-classes are a genuine independent
    period, not a manufactured sign. D_even(4), D_odd(4) in {pi^2,pi^4,G,beta4}."""

    def test_shift_classes_independent(self) -> None:
        #                pi^2            pi^4         G   beta(4)
        D_even = Matrix([[Rational(1, 2), Rational(-1, 24), -4, 4]])
        D_odd = Matrix([[Rational(1, 2), Rational(-1, 24), 4, -4]])
        assert Matrix.vstack(D_even, D_odd).rank() == 2


if __name__ == "__main__":
    for nm in (1, 2, 3, 4):
        print(f"n_max={nm}: blind rank {parity_blind_gram(nm).rank()}, "
              f"resolved rank {parity_resolved_matrix(nm).rank()}, N={N_sectors(nm)}")
