"""Paper 2 ingredient tests: B = 42, the selection principle, the Z3 circulant.

Backs paper_2_alpha.tex Section II-III (spectral invariants of the Hopf
bundle) and Section VI (the circulant structure), per the group5 Tier-2
track E backing build (2026-07-02).  Companion to
tests/test_paper2_corrections.py, which pins the cubic's physical root and
the GENERALIZED S^d identity B_formal(m)/N(m) = d (test_paper2_viiid_identity_k0);
neither is duplicated here.

Coverage:
  (a) B = 42: the degeneracy-weighted Casimir trace of the S^2 base at
      n_max = 3 (Eq. eq:casimir / boxed eq:base), recomputed exactly from
      the paper's own definition, including the per-(n,l) term breakdown.
  (b) The selection principle (Eq. eq:selection): B(n_max)/N(n_max) =
      dim(S^3) = 3 holds at n_max = 3 and FAILS at n_max = 1, 2, 4, 5
      (Table tab:selection uniqueness), plus the closed forms
      eq:innersum / eq:Bclosed / eq:ratioclosed and the algebraic
      uniqueness proof (m+4)(m-3) = 0.
  (c) The Z3 circulant (Eq. eq:circulant): traceless, det M = b^3 + c^3,
      characteristic polynomial det(M - sI) = -s^3 + 3bc s + (b^3 + c^3);
      under the paper's conditions bc = K/3 and b^3 + c^3 = -1 it becomes
      s^3 - K s + 1 -- exactly the paper's cubic.  Plus the SU(2)
      structure-constant analog b = c = 1 with char poly s^3 - 3s - 2.

Every assertion is exact (sympy / Fraction); no floats, no tolerances.

Register tier: these are INTERNAL-THEOREM-level ingredient identities.
They do NOT bear on the combination rule K = pi(B + F - Delta), which
remains an Observation (CLAUDE.md 13.5 hard prohibition).
"""

from fractions import Fraction

import sympy as sp


# ---------------------------------------------------------------------------
# Paper definitions, computed exactly
# ---------------------------------------------------------------------------

def B_dw_casimir(n_max: int) -> Fraction:
    """Degeneracy-weighted Casimir trace, Eq. (eq:casimir):
    B(n_max) = sum_{n=1}^{n_max} sum_{l=0}^{n-1} (2l+1) l(l+1)."""
    return Fraction(sum((2 * l + 1) * l * (l + 1)
                        for n in range(1, n_max + 1) for l in range(n)))


def N_states(n_max: int) -> Fraction:
    """Total state count through shell n_max: N = sum_{n<=n_max} n^2."""
    return Fraction(sum(n ** 2 for n in range(1, n_max + 1)))


# ---------------------------------------------------------------------------
# (a) B = 42
# ---------------------------------------------------------------------------

class TestB42:
    """The base-manifold invariant B = 42 (boxed Eq. eq:base)."""

    def test_paper2_b42_exact(self) -> None:
        """B(3) = 42 exactly from the paper's own definition."""
        assert B_dw_casimir(3) == 42

    def test_paper2_b42_term_breakdown(self) -> None:
        """The explicit per-(n,l) decomposition printed in the paper:
        0_(1,0) + 0_(2,0) + 6_(2,1) + 0_(3,0) + 6_(3,1) + 30_(3,2) = 42."""
        term = {(n, l): (2 * l + 1) * l * (l + 1)
                for n in range(1, 4) for l in range(n)}
        assert term[(1, 0)] == 0
        assert term[(2, 0)] == 0
        assert term[(2, 1)] == 6
        assert term[(3, 0)] == 0
        assert term[(3, 1)] == 6
        assert term[(3, 2)] == 30
        assert sum(term.values()) == 42

    def test_paper2_inner_sum_closed_form(self) -> None:
        """Eq. (eq:innersum): sum_{l=0}^{n-1} (2l+1) l(l+1) = n^2(n^2-1)/2,
        proven symbolically for general integer n."""
        n, l = sp.symbols("n l", integer=True, positive=True)
        inner = sp.summation((2 * l + 1) * l * (l + 1), (l, 0, n - 1))
        assert sp.simplify(inner - n ** 2 * (n ** 2 - 1) / 2) == 0

    def test_paper2_per_shell_casimir_average(self) -> None:
        """Eq. (eq:casimir_identity): the per-shell degeneracy-weighted
        SO(3) Casimir average equals the SO(4) Casimir eigenvalue,
        (1/n^2) sum_l (2l+1) l(l+1) = (n^2 - 1)/2."""
        n, l = sp.symbols("n l", integer=True, positive=True)
        avg = sp.summation((2 * l + 1) * l * (l + 1), (l, 0, n - 1)) / n ** 2
        assert sp.simplify(avg - (n ** 2 - 1) / 2) == 0


# ---------------------------------------------------------------------------
# (b) Selection principle: B/N = 3 uniquely at n_max = 3
# ---------------------------------------------------------------------------

class TestSelectionPrinciple:
    """B(n_max)/N(n_max) = dim(S^3) = 3 selects n_max = 3 uniquely.

    The generalized S^d identity B_formal/N = d is pinned in
    test_paper2_corrections.py::test_paper2_viiid_identity_k0 and is
    deliberately NOT re-tested here; this class pins the S^3 SELECTION
    (holds at 3, fails at 1, 2, 4, 5) and the paper's closed forms.
    """

    def test_paper2_selection_holds_at_3(self) -> None:
        assert B_dw_casimir(3) / N_states(3) == Fraction(42, 14) == 3

    def test_paper2_selection_fails_off_3(self) -> None:
        """Uniqueness rows of Table tab:selection: exact ratios at
        n_max = 1, 2, 4, 5 (0, 6/5, 27/5, 42/5) -- none equals 3."""
        expected = {
            1: Fraction(0, 1),
            2: Fraction(6, 5),
            4: Fraction(27, 5),
            5: Fraction(42, 5),
        }
        for m, ratio in expected.items():
            got = B_dw_casimir(m) / N_states(m)
            assert got == ratio, (m, got)
            assert got != 3, m

    def test_paper2_selection_table_values(self) -> None:
        """Table tab:selection columns B and N, exactly."""
        assert [B_dw_casimir(m) for m in range(1, 6)] == [0, 6, 42, 162, 462]
        assert [N_states(m) for m in range(1, 6)] == [1, 5, 14, 30, 55]

    def test_paper2_B_closed_form(self) -> None:
        """Eq. (eq:Bclosed): B(m) = m(m+1)(2m+1)(m+2)(m-1)/20, symbolically."""
        m, n = sp.symbols("m n", integer=True, positive=True)
        B_sum = sp.summation(n ** 2 * (n ** 2 - 1) / 2, (n, 1, m))
        closed = m * (m + 1) * (2 * m + 1) * (m + 2) * (m - 1) / 20
        assert sp.simplify(B_sum - closed) == 0

    def test_paper2_ratio_closed_form_and_uniqueness(self) -> None:
        """Eq. (eq:ratioclosed): B(m)/N(m) = 3(m+2)(m-1)/10, and the
        algebraic uniqueness proof: B/N = 3 iff m^2 + m - 12 =
        (m+4)(m-3) = 0, whose unique positive root is m = 3."""
        m = sp.symbols("m", positive=True)
        B_closed = m * (m + 1) * (2 * m + 1) * (m + 2) * (m - 1) / 20
        N_closed = m * (m + 1) * (2 * m + 1) / 6
        ratio = sp.simplify(B_closed / N_closed)
        assert sp.simplify(ratio - sp.Rational(3, 10) * (m + 2) * (m - 1)) == 0
        # (m+2)(m-1) = 10  <=>  m^2 + m - 12 = 0  =  (m+4)(m-3)
        quadratic = sp.expand((m + 2) * (m - 1) - 10)
        assert quadratic == sp.expand(m ** 2 + m - 12)
        assert sp.expand((m + 4) * (m - 3)) == quadratic
        roots = sp.solve(sp.Eq(ratio, 3), m)
        assert roots == [3]


# ---------------------------------------------------------------------------
# (c) The Z3 circulant and the cubic's characteristic-polynomial origin
# ---------------------------------------------------------------------------

class TestCirculant:
    """Section VI: the cubic alpha^3 - K alpha + 1 = 0 as the
    characteristic polynomial of a traceless Z3 circulant."""

    @staticmethod
    def _circulant(b, c):
        return sp.Matrix([[0, b, c], [c, 0, b], [b, c, 0]])

    def test_paper2_circulant_traceless_and_det(self) -> None:
        """Conditions 1 and 3 of the paper: tr(M) = 0 (kills the s^2 term)
        and det(M) = b^3 + c^3 (so det M = -1 <=> b^3 + c^3 = -1)."""
        b, c = sp.symbols("b c")
        M = self._circulant(b, c)
        assert M.trace() == 0
        assert sp.expand(M.det() - (b ** 3 + c ** 3)) == 0

    def test_paper2_circulant_charpoly_shape(self) -> None:
        """det(M - sI) = -s^3 + 3bc s + (b^3 + c^3), symbolically -- the
        paper's Eq. after (eq:circulant).  Equivalently charpoly(s) =
        det(sI - M) = s^3 - 3bc s - (b^3 + c^3)."""
        b, c, s = sp.symbols("b c s")
        M = self._circulant(b, c)
        det_form = sp.expand((M - s * sp.eye(3)).det())
        assert sp.expand(det_form - (-s ** 3 + 3 * b * c * s
                                     + b ** 3 + c ** 3)) == 0
        charpoly = M.charpoly(s).as_expr()
        assert sp.expand(charpoly - (s ** 3 - 3 * b * c * s
                                     - (b ** 3 + c ** 3))) == 0

    def test_paper2_circulant_gives_papers_cubic(self) -> None:
        """Substituting the paper's conditions bc = K/3 and
        b^3 + c^3 = -1 into the proven char-poly shape yields exactly
        s^3 - K s + 1 -- the paper's cubic with s = alpha."""
        s, K, e2, e3 = sp.symbols("s K e2 e3")
        # proven shape in the two symmetric functions e2 = bc, e3 = b^3+c^3
        shape = s ** 3 - 3 * e2 * s - e3
        cubic = shape.subs({e2: K / 3, e3: -1})
        assert sp.expand(cubic - (s ** 3 - K * s + 1)) == 0

    def test_paper2_su2_structure_constant_analog(self) -> None:
        """The SU(2) epsilon_ijk analog b = c = 1: char poly s^3 - 3s - 2
        = (s-2)(s+1)^2, i.e. eigenvalues {2, -1, -1}."""
        s = sp.symbols("s")
        M = self._circulant(1, 1)
        charpoly = M.charpoly(s).as_expr()
        assert sp.expand(charpoly - (s ** 3 - 3 * s - 2)) == 0
        assert sp.factor(charpoly) == (s - 2) * (s + 1) ** 2
        assert M.eigenvals() == {sp.Integer(2): 1, sp.Integer(-1): 2}
