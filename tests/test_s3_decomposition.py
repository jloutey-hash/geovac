"""Frozen falsifier: S^(3) three-loop chain decomposition and closure.

Closure sprint 2026-06-11 (debug/sprint_s3_closure_memo.md; stage-1
decomposition: debug/sprint_s3_decomposition_memo.md).

Five test groups:
  (i)   Chain-collapse identity in exact Fraction arithmetic: the
        channel-count closed form I(n1,n2) = 2*min-1-[n1=n2] factorizes
        the three-loop chain cube sum exactly, cross-checked against the
        real SO(4) CG machinery at small N.
  (ii)  Trailing-1 Abel evaluator vs rigorous direct-sum brackets — the
        precision-contract regression guard (a frozen-precision term
        cache fed to nsum/Levin reproduces the stage-1 mis-convergence;
        the per-precision cache must land inside the bracket).
  (iii) The four trailing-1 t3 identifications (closure-sprint PSLQ at
        220/340 dps against complete weight-homogeneous level-2 bases).
  (iv)  Stuffle witnesses R2/R6 at 45 dps (slow) — the stage-1 residuals
        6.4e-12 / 9.8e-20 decomposed exactly into the trailing-1
        evaluator errors; with correct values both must pass.
  (v)   Global anchor bracket: rigorous lower bound (positive-term
        partial sum) + integral tail bound must contain the canonical
        S^(3) value and EXCLUDE both wrong stage-1 figures (30.6154
        decomposition-with-poisoned-entries, 30.2197 log-blind Levin).

Runtime: fast groups < 60 s total; group (iv) is marked slow.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import mpmath
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from geovac.qed_vertex import so4_channel_count, _vertex_allowed

# ---------------------------------------------------------------------------
# Pinned constants (closure sprint, debug/data/s3_closure_*.json)
# ---------------------------------------------------------------------------

# Canonical S^(3) at physical exponents (a=4, p=1), corrected assembly
# (s3_closure_recompute.py, 200 dps; first 60 digits pinned).
S3_CANONICAL_50 = (
    "31.5725612075120227547619295450689833719575844695771941287005")

# The two refuted stage-1 figures (sprint_s3_decomposition_memo.md):
S3_STALE_DECOMPOSITION = "30.61538052881950"
S3_STALE_LEVIN_ANCHOR = "30.2197"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def phi_frac(n: int) -> Fraction:
    """phi(n) = g_n/|lambda_n|^4 = 8*(o^2-1)/o^4 at o = 2n+3."""
    o = 2 * n + 3
    return Fraction(8 * (o * o - 1), o**4)


def I_closed_form(n1: int, n2: int) -> int:
    m = min(n1, n2)
    if m == 0:
        return 0
    return 2 * m - 1 - (1 if n1 == n2 else 0)


def I_channel_geovac(n1: int, n2: int) -> int:
    tot = 0
    for q in range(1, n1 + n2 + 1):
        if _vertex_allowed(n1, n2, q):
            tot += so4_channel_count(n1, n2, q)
    return tot


def lam(s: int) -> mpmath.mpf:
    return (1 - mpmath.power(2, -s)) * mpmath.zeta(s)


def t2val(b: int, c: int) -> mpmath.mpf:
    """t(b,c) = sum_{o2>o1>=1 odd} o2^{-b} o1^{-c} (descending convention)."""
    def term(j):
        j = int(j)
        return (mpmath.mpf(2 * j + 1) ** (-c) * mpmath.power(2, -b)
                * mpmath.zeta(b, mpmath.mpf(j) + mpmath.mpf(3) / 2))
    return mpmath.nsum(term, [0, mpmath.inf], method="levin")


def t3val_b3ge2(b1: int, b2: int, b3: int) -> mpmath.mpf:
    """t3(b1,b2,b3) for b3 >= 2 (log-free summand, Levin-safe)."""
    lam3 = lam(b3)

    def term(j):
        j = int(j)
        tau = mpmath.power(2, -b1) * mpmath.zeta(
            b1, mpmath.mpf(j) + mpmath.mpf(3) / 2)
        pL = lam3 - mpmath.power(2, -b3) * mpmath.zeta(
            b3, mpmath.mpf(j) + mpmath.mpf(1) / 2)
        return mpmath.mpf(2 * j + 1) ** -b2 * tau * pL
    return mpmath.nsum(term, [1, mpmath.inf], method="levin")


class Trailing1:
    """t3(b1,b2,1) via Abel summation with per-precision term caches.

    Mirrors debug/s3_closure_trailing1.py (v3).  The per-precision cache
    honors mpmath nsum's contract: the Levin transform re-requests terms
    at elevated internal precision; serving frozen-precision terms lets
    the transform's cancellations amplify rounding error (the stage-1
    failure mode, observed at 1e-4..1e-11).
    """

    NTERMS0 = 2048

    def __init__(self, b1: int, b2: int) -> None:
        self.b1, self.b2 = b1, b2
        self._cache: dict = {}

    def _a_direct(self, j):
        j = int(j)
        return (mpmath.mpf(2 * j + 1) ** -self.b2 * mpmath.power(2, -self.b1)
                * mpmath.zeta(self.b1, mpmath.mpf(j) + mpmath.mpf(3) / 2))

    def _new_state(self) -> dict:
        return {
            "A1": mpmath.nsum(self._a_direct, [1, mpmath.inf],
                              method="levin"),
            "z": mpmath.zeta(self.b1, mpmath.mpf(5) / 2),
            "j_done": 1,
            "cum": mpmath.mpf(0),
            "terms": [],
        }

    def _fill(self, st: dict, k_target: int) -> None:
        terms = st["terms"]
        while len(terms) < k_target:
            k = len(terms) + 1
            terms.append((st["A1"] - st["cum"]) / (2 * k - 1))
            j = st["j_done"]
            a_j = (mpmath.mpf(2 * j + 1) ** -self.b2
                   * mpmath.power(2, -self.b1) * st["z"])
            st["cum"] += a_j
            st["z"] -= (mpmath.mpf(j) + mpmath.mpf(3) / 2) ** -self.b1
            st["j_done"] += 1

    def term(self, k):
        k = int(k)
        prec = mpmath.mp.prec
        st = self._cache.get(prec)
        if st is None:
            st = self._new_state()
            self._fill(st, self.NTERMS0)
            self._cache[prec] = st
        if k > len(st["terms"]):
            self._fill(st, 2 * k)
        return st["terms"][k - 1]

    def value(self):
        return mpmath.nsum(self.term, [1, mpmath.inf], method="levin")


def trailing1_bracket(b1: int, b2: int, J: int = 20000):
    """Rigorous bracket for t3(b1,b2,1): direct partial sum (strict lower
    bound, positive terms) + integral tail bound (driver docstring has
    the derivation)."""
    z = mpmath.zeta(b1, mpmath.mpf(5) / 2)
    H = mpmath.mpf(1)
    pref = mpmath.power(2, -b1)
    partial = mpmath.mpf(0)
    for j in range(1, J + 1):
        partial += mpmath.mpf(2 * j + 1) ** -b2 * pref * z * H
        z -= (mpmath.mpf(j) + mpmath.mpf(3) / 2) ** -b1
        H += mpmath.mpf(1) / (2 * j + 1)
    m_exp = b1 + b2 - 1
    Cz = 1 + mpmath.mpf(b1 - 1) / (J + mpmath.mpf(1) / 2)
    pref_tail = mpmath.power(2, -b1) * mpmath.power(2, -b2) * Cz / (b1 - 1)
    Jm = mpmath.mpf(J)
    tail = (pref_tail * Jm ** (1 - m_exp)
            * ((mpmath.log(Jm) / 2 + mpmath.mpf(3) / 2) / (m_exp - 1)
               + mpmath.mpf(1) / (2 * (m_exp - 1) ** 2)))
    return partial, partial + tail


# ===========================================================================
# Group (i): chain-collapse identity, exact Fractions
# ===========================================================================

class TestChainCollapse:
    """The channel-count closed form factorizes the chain cube sum exactly:

    sum_{n1,n2,n3<=N} I(n1,n2) I(n2,n3) phi(n1) phi(n2) phi(n3)
        == sum_{m<=N} phi(m) * L_N(m)^2,
    L_N(m) = sum_{n<=N} I(n,m) phi(n) = 2 A_N(m) - P_N - phi(m).
    """

    def test_factorization_exact_fraction(self) -> None:
        N = 25
        phis = {n: phi_frac(n) for n in range(1, N + 1)}
        direct = Fraction(0)
        for n2 in range(1, N + 1):
            inner = sum(I_closed_form(n, n2) * phis[n]
                        for n in range(1, N + 1))
            direct += phis[n2] * inner * inner
        P_N = sum(phis.values())
        factorized = Fraction(0)
        for m in range(1, N + 1):
            A_N = sum(min(n, m) * phis[n] for n in range(1, N + 1))
            L_N = 2 * A_N - P_N - phis[m]
            factorized += phis[m] * L_N * L_N
        assert direct == factorized

    def test_chain_vs_cg_machinery_small_N(self) -> None:
        """Closed-form chain == real SO(4) CG chain at N=6 (exact)."""
        N = 6
        phis = {n: phi_frac(n) for n in range(1, N + 1)}
        cg = Fraction(0)
        cf = Fraction(0)
        for n2 in range(1, N + 1):
            inner_cg = sum(I_channel_geovac(n, n2) * phis[n]
                           for n in range(1, N + 1))
            inner_cf = sum(I_closed_form(n, n2) * phis[n]
                           for n in range(1, N + 1))
            cg += phis[n2] * inner_cg * inner_cg
            cf += phis[n2] * inner_cf * inner_cf
        assert cg == cf


# ===========================================================================
# Group (ii): trailing-1 evaluator vs rigorous brackets
# ===========================================================================

class TestTrailing1Evaluator:
    """The Abel evaluator must land inside the rigorous direct-sum
    bracket — the gate that caught the stage-1 mis-convergence (old
    t3(2,1,1) was below its bracket by 2.9e-4) and the v2 frozen-
    precision bug (all six below)."""

    def test_t3_211_in_bracket(self) -> None:
        mpmath.mp.dps = 40
        lo, hi = trailing1_bracket(2, 1)
        v = Trailing1(2, 1).value()
        assert lo <= v <= hi, (
            "t3(2,1,1) = %s outside rigorous bracket [%s, %s]"
            % (mpmath.nstr(v, 20), mpmath.nstr(lo, 20), mpmath.nstr(hi, 20)))

    def test_t3_431_in_bracket(self) -> None:
        mpmath.mp.dps = 40
        lo, hi = trailing1_bracket(4, 3)
        v = Trailing1(4, 3).value()
        assert lo <= v <= hi

    def test_stage1_value_excluded(self) -> None:
        """The stage-1 cached t3(2,1,1) (error -2.9e-4) must fail the
        bracket — regression guard against silently re-trusting Levin
        on log-modulated summands."""
        mpmath.mp.dps = 40
        lo, _ = trailing1_bracket(2, 1)
        stale = mpmath.mpf("0.13952268")  # stage-A cache value (rounded)
        assert stale < lo


# ===========================================================================
# Group (iii): the four trailing-1 identifications
# ===========================================================================

class TestTrailing1Identifications:
    """PSLQ identifications from the closure sprint (G4, residuals
    2.6e-216 .. 7.3e-339; debug/data/s3_closure_trailing1.json):

      t3(2,1,1) = (1/2) Li4(1/2) + (1/48) ln^4 2 + (1/24) pi^2 ln^2 2
                  - (19/5760) pi^4                            [w4, depth<=1]
      t3(4,1,1) = (1/192) pi^4 ln^2 2 - (1/64) pi^2 ln2 z3
                  - (11/322560) pi^6 - (1/2) t2(5,1) - (7/128) z3^2   [w6]
      t3(3,2,1) = (3/64) pi^2 ln2 z3 - (1/9216) pi^6
                  - (1/2) t2(5,1) - (49/256) z3^2                     [w6]
      t3(3,4,1) = (1/768) pi^4 ln2 z3 + (5/128) pi^2 ln2 z5
                  + (61/2903040) pi^8 - (1/256) pi^2 z3^2
                  + (1/2) t2(5,3) - (1/2) t2(7,1) - (217/512) z3 z5   [w8]

    The k=3 headline: every component reduces to CATALOGUED level-2
    objects; the genuine depth-2 generators t2(5,1), t2(5,3), t2(7,1)
    appear with coefficients +-1/2 — one rung up from S_min's depth<=1,
    still classical.
    """

    DPS = 45
    TOL = mpmath.mpf(10) ** -30

    def test_t3_211_closed_form(self) -> None:
        """w4: pure depth-1 monomials (fast — the headline poisoned value)."""
        mpmath.mp.dps = self.DPS
        pi, ln2 = mpmath.pi, mpmath.log(2)
        formula = (
            mpmath.polylog(4, mpmath.mpf(1) / 2) / 2
            + ln2 ** 4 / 48
            + pi ** 2 * ln2 ** 2 / 24
            - mpmath.mpf(19) / 5760 * pi ** 4
        )
        v = Trailing1(2, 1).value()
        assert abs(v - formula) < self.TOL, (
            "t3(2,1,1) closed form: residual %s"
            % mpmath.nstr(abs(v - formula), 5))

    @pytest.mark.slow
    def test_t3_411_closed_form(self) -> None:
        mpmath.mp.dps = self.DPS
        pi, ln2, z3 = mpmath.pi, mpmath.log(2), mpmath.zeta(3)
        formula = (
            pi ** 4 * ln2 ** 2 / 192
            - pi ** 2 * ln2 * z3 / 64
            - mpmath.mpf(11) / 322560 * pi ** 6
            - t2val(5, 1) / 2
            - mpmath.mpf(7) / 128 * z3 ** 2
        )
        v = Trailing1(4, 1).value()
        assert abs(v - formula) < self.TOL

    @pytest.mark.slow
    def test_t3_321_closed_form(self) -> None:
        mpmath.mp.dps = self.DPS
        pi, ln2, z3 = mpmath.pi, mpmath.log(2), mpmath.zeta(3)
        formula = (
            mpmath.mpf(3) / 64 * pi ** 2 * ln2 * z3
            - pi ** 6 / 9216
            - t2val(5, 1) / 2
            - mpmath.mpf(49) / 256 * z3 ** 2
        )
        v = Trailing1(3, 2).value()
        assert abs(v - formula) < self.TOL

    @pytest.mark.slow
    def test_t3_341_closed_form(self) -> None:
        mpmath.mp.dps = self.DPS
        pi, ln2 = mpmath.pi, mpmath.log(2)
        z3, z5 = mpmath.zeta(3), mpmath.zeta(5)
        formula = (
            pi ** 4 * ln2 * z3 / 768
            + mpmath.mpf(5) / 128 * pi ** 2 * ln2 * z5
            + mpmath.mpf(61) / 2903040 * pi ** 8
            - pi ** 2 * z3 ** 2 / 256
            + t2val(5, 3) / 2
            - t2val(7, 1) / 2
            - mpmath.mpf(217) / 512 * z3 * z5
        )
        v = Trailing1(3, 4).value()
        assert abs(v - formula) < self.TOL


# ===========================================================================
# Group (iv): stuffle witnesses R2, R6  (slow)
# ===========================================================================

@pytest.mark.slow
class TestStuffleWitnesses:
    """R2: lam(3) t2(2,1) = t3(3,2,1)+t3(2,3,1)+t3(2,1,3)+t2(5,1)+t2(2,4)
    R6: lam(3) t2(4,1) = t3(3,4,1)+t3(4,3,1)+t3(4,1,3)+t2(7,1)+t2(4,4)

    Stage-1 residuals 6.4e-12 / 9.8e-20 decomposed exactly into the
    trailing-1 evaluator errors (closure memo table); with the fixed
    evaluator both must pass at 1e-30."""

    DPS = 45
    TOL = mpmath.mpf(10) ** -30

    def test_r2(self) -> None:
        mpmath.mp.dps = self.DPS
        lhs = lam(3) * t2val(2, 1)
        rhs = (Trailing1(3, 2).value() + Trailing1(2, 3).value()
               + t3val_b3ge2(2, 1, 3) + t2val(5, 1) + t2val(2, 4))
        assert abs(lhs - rhs) < self.TOL, (
            "R2 residual %s" % mpmath.nstr(abs(lhs - rhs), 5))

    def test_r6(self) -> None:
        mpmath.mp.dps = self.DPS
        lhs = lam(3) * t2val(4, 1)
        rhs = (Trailing1(3, 4).value() + Trailing1(4, 3).value()
               + t3val_b3ge2(4, 1, 3) + t2val(7, 1) + t2val(4, 4))
        assert abs(lhs - rhs) < self.TOL, (
            "R6 residual %s" % mpmath.nstr(abs(lhs - rhs), 5))


# ===========================================================================
# Group (v): global anchor bracket
# ===========================================================================

class TestAnchorBracket:
    """Rigorous bracket for S^(3) = sum_m phi(m) L(m)^2 (factorized chain,
    closed-form full-range L): monotone partial sum (lower) + integral
    tail bound 32[(ln N+1)^2+2(ln N+1)+2]/N (upper).  Derivations in
    debug/s3_closure_anchor.py.  Must contain the canonical value and
    exclude BOTH wrong stage-1 figures."""

    N = 20000

    def _bracket(self):
        mpmath.mp.dps = 30
        pi = mpmath.pi
        P = pi ** 2 - pi ** 4 / 12 - mpmath.mpf(64) / 81
        head = mpmath.mpf(0)
        cum_phi = mpmath.mpf(0)
        cum_nphi = mpmath.mpf(0)
        for m in range(1, self.N + 1):
            x = mpmath.mpf(2 * m + 3)
            x2 = x * x
            pm = 8 * (x2 - 1) / (x2 * x2)
            A = cum_nphi + m * (P - cum_phi)
            L = 2 * A - P - pm
            head += pm * L * L
            cum_phi += pm
            cum_nphi += m * pm
        u = mpmath.log(self.N) + 1
        tail = 32 * (u * u + 2 * u + 2) / self.N
        return head, head + tail

    def test_bracket_excludes_stale_figures(self) -> None:
        lo, hi = self._bracket()
        for stale in (S3_STALE_DECOMPOSITION, S3_STALE_LEVIN_ANCHOR):
            v = mpmath.mpf(stale)
            assert not (lo <= v <= hi), (
                "stale figure %s inside bracket [%s, %s]"
                % (stale, mpmath.nstr(lo, 15), mpmath.nstr(hi, 15)))

    def test_bracket_contains_canonical(self) -> None:
        lo, hi = self._bracket()
        v = mpmath.mpf(S3_CANONICAL_50)
        assert lo <= v <= hi, (
            "canonical S^(3) %s outside bracket [%s, %s]"
            % (S3_CANONICAL_50[:18], mpmath.nstr(lo, 15), mpmath.nstr(hi, 15)))
