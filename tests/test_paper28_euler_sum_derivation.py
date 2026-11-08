"""Paper 28 sec:euler_cancel -- derivation-level guards (2026-11-07 guard pass).

Backs claims that were printed with no test, or pinned only to frozen
decimals of the printed closed form:

  * the Sommerfeld coefficients c_2(n) .. c_5(n) (sec. "Closed forms for
    c_p(n)", "D_4 at weight 7", "D_5 at weight 9");
  * D_2, D_3 (single-zeta collapse);
  * Theorem thm:zeta3_complement (eq:shell_complement,
    eq:dirichlet_complement, the s = 4 evaluation);
  * the S_{r,s} Euler-sum tables at weights 7, 9, 11;
  * eq:D4_decomp, eq:D5_decomp, eq:D6_decomp and eq:D4_three_tier;
  * obs:product_rule through p = 6 (zeta(2) zeta(2p-3) cancels; the number
    of surviving products).

Independent routes (none shares a code path with the printed formula):

  R1  exact series.  sum_j d(n,j) E_nj/mc^2 expanded in x = (Z alpha)^2 by
      truncated power-series arithmetic in Fractions, straight from the exact
      Dirac formula (no sympy, no geovac).  This gives c_p(n) exactly at each
      integer n.
  R2  exact ansatz fit.  c_p(n) = a n^{-(2p-2)} + b n^{-(2p-1)}
      + sum_{r=1}^{2p-4} beta_r H_{n-1}^{(r)} / n^{2p-1-r}; the coefficients
      are solved from R1 at n = 1..dim by exact Gaussian elimination and the
      ansatz is then CHECKED at n = dim+1 .. 20 (so the basis is tested, not
      assumed).  For p = 6 this is the only route to c_6, which the paper
      does not print.
  R3  numeric MZVs.  zeta(s, r) = sum_{m>k>=1} m^-s k^-r is evaluated as
      sum_{k>=1} k^-r zeta(s, k+1) (Hurwitz), a log-free summand, by mpmath
      nsum at 60 dps.  S_{r,s} = zeta(s, r) + zeta(r+s).  D_p then follows
      from R2's coefficients: sum_n H_{n-1}^{(r)} n^{-s} = zeta(s, r).
  R4  exact rational bookkeeping.  D_p re-assembled from R2's coefficients
      and the PRINTED S_{r,s} table rows (plus Euler's classical
      S_{1,m} formula for S_{1,8}, which the weight-9 table omits), giving
      the printed D_p decomposition as an exact rational vector.

Runtime: about 10 s in total.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, List, Tuple

import pytest
from mpmath import mp, mpf, nsum, inf, zeta, pi

# ---------------------------------------------------------------------------
# R1: exact series of sum_j d(n,j) E_nj / mc^2 in x = (Z alpha)^2
# ---------------------------------------------------------------------------

ORDER = 7  # keep x^0 .. x^7


def _mul(a: List[Fr], b: List[Fr]) -> List[Fr]:
    c = [Fr(0)] * (ORDER + 1)
    for i, ai in enumerate(a):
        if ai:
            for j in range(ORDER + 1 - i):
                c[i + j] += ai * b[j]
    return c


def _inv(a: List[Fr]) -> List[Fr]:
    b = [Fr(0)] * (ORDER + 1)
    b[0] = 1 / a[0]
    for k in range(1, ORDER + 1):
        b[k] = -sum(a[i] * b[k - i] for i in range(1, k + 1)) / a[0]
    return b


def _binom(alpha: Fr, i: int) -> Fr:
    r = Fr(1)
    for t in range(i):
        r = r * (alpha - t) / (t + 1)
    return r


def _one_plus_pow(y: List[Fr], alpha: Fr) -> List[Fr]:
    """(1 + y)^alpha for a series y with y[0] = 0."""
    out = [Fr(0)] * (ORDER + 1)
    yp = [Fr(1)] + [Fr(0)] * ORDER
    for i in range(ORDER + 1):
        coef = _binom(alpha, i)
        for k in range(ORDER + 1):
            out[k] += coef * yp[k]
        yp = _mul(yp, y)
    return out


def shell_series(n: int) -> List[Fr]:
    """Coefficients of sum_j d(n,j) E_nj/mc^2 in powers of x.

    E/mc^2 = (1 + x / (n_r + sqrt(kappa^2 - x))^2)^(-1/2) - 1, kappa = j+1/2,
    n_r = n - kappa, d = 2(2j+1) = 4 kappa for kappa < n and 2 kappa at
    kappa = n."""
    tot = [Fr(0)] * (ORDER + 1)
    x = [Fr(0), Fr(1)] + [Fr(0)] * (ORDER - 1)
    for kap in range(1, n + 1):
        d = 4 * kap if kap < n else 2 * kap
        # sqrt(kap^2 - x) = kap * sum_i binom(1/2, i) (-x/kap^2)^i
        u = [_binom(Fr(1, 2), i) * (-1) ** i * Fr(kap) ** (1 - 2 * i)
             for i in range(ORDER + 1)]
        u[0] += n - kap
        iu = _inv(u)
        E = _one_plus_pow(_mul(x, _mul(iu, iu)), Fr(-1, 2))
        E[0] -= 1
        for k in range(ORDER + 1):
            tot[k] += d * E[k]
    return tot


N_SERIES = 20
_SER: Dict[int, List[Fr]] = {}


def ser(n: int) -> List[Fr]:
    if n not in _SER:
        _SER[n] = shell_series(n)
    return _SER[n]


def H(n: int, r: int) -> Fr:
    """H_{n-1}^{(r)} = sum_{k=1}^{n-1} k^-r."""
    return sum((Fr(1, k ** r) for k in range(1, n)), Fr(0))


# ---------------------------------------------------------------------------
# Printed closed forms (Paper 28 sec:euler_cancel)
# ---------------------------------------------------------------------------

def c2_printed(n: int) -> Fr:
    return -Fr(5 * n - 4, 4 * n ** 3)


def c3_printed(n: int) -> Fr:
    return (Fr(19 * n - 20, 8 * n ** 5) - Fr(3, 2) * H(n, 1) / n ** 4
            - Fr(1, 2) * H(n, 2) / n ** 3)


def c4_printed(n: int) -> Fr:
    return (-Fr(5 * (41 * n - 40), 64 * n ** 7) + Fr(15, 4) * H(n, 1) / n ** 6
            - Fr(1, 4) * H(n, 2) / n ** 5 - Fr(3, 4) * H(n, 3) / n ** 4
            - Fr(1, 4) * H(n, 4) / n ** 3)


def c5_printed(n: int) -> Fr:
    return (Fr(7 * (71 * n - 72), 128 * n ** 9)
            - Fr(105, 16) * H(n, 1) / n ** 8 + Fr(45, 16) * H(n, 2) / n ** 7
            + Fr(5, 4) * H(n, 3) / n ** 6 - Fr(3, 8) * H(n, 4) / n ** 5
            - Fr(15, 32) * H(n, 5) / n ** 4 - Fr(5, 32) * H(n, 6) / n ** 3)


PRINTED_C = {2: c2_printed, 3: c3_printed, 4: c4_printed, 5: c5_printed}

# ---------------------------------------------------------------------------
# R2: exact ansatz fit of c_p(n)
# ---------------------------------------------------------------------------


def _basis(p: int, n: int) -> List[Fr]:
    w = 2 * p - 1
    return ([Fr(1, n ** (w - 1)), Fr(1, n ** w)]
            + [H(n, r) / n ** (w - r) for r in range(1, 2 * p - 3)])


def _solve(A: List[List[Fr]], b: List[Fr]) -> List[Fr]:
    m = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(m):
        piv = next(r for r in range(col, m) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        for r in range(m):
            if r != col and M[r][col] != 0:
                f = M[r][col] / M[col][col]
                M[r] = [x - f * y for x, y in zip(M[r], M[col])]
    return [M[i][m] / M[i][i] for i in range(m)]


_FIT: Dict[int, List[Fr]] = {}


def fit(p: int) -> List[Fr]:
    """[a_{2p-2}, a_{2p-1}, beta_1, ..., beta_{2p-4}] from R1, checked on
    extra n."""
    if p not in _FIT:
        dim = 2 + (2 * p - 4)
        A = [_basis(p, n) for n in range(1, dim + 1)]
        b = [ser(n)[p] for n in range(1, dim + 1)]
        coef = _solve(A, b)
        for n in range(dim + 1, N_SERIES + 1):
            val = sum(c * f for c, f in zip(coef, _basis(p, n)))
            assert val == ser(n)[p], ("ansatz fails", p, n)
        _FIT[p] = coef
    return _FIT[p]


# ---------------------------------------------------------------------------
# R3: numeric MZVs
# ---------------------------------------------------------------------------

DPS = 60


@pytest.fixture(autouse=True)
def _mp_hermetic():
    saved = mp.dps
    yield
    mp.dps = saved


_MZV: Dict[Tuple[int, int], object] = {}


def mzv2(s: int, r: int):
    """zeta(s, r) = sum_{m > k >= 1} m^-s k^-r (s >= 2), at DPS."""
    key = (s, r)
    if key not in _MZV:
        with mp.workdps(DPS + 10):
            v = nsum(lambda k: k ** (-r) * zeta(s, k + 1), [1, inf])
        _MZV[key] = v
    return _MZV[key]


def D_numeric(p: int):
    """D_p = sum_n c_p(n) from R2's coefficients and R3's MZVs."""
    coef = fit(p)
    w = 2 * p - 1
    v = mpf(coef[0].numerator) / coef[0].denominator * zeta(w - 1) \
        + mpf(coef[1].numerator) / coef[1].denominator * zeta(w)
    for r in range(1, 2 * p - 3):
        c = coef[r + 1]
        if c:
            v += mpf(c.numerator) / c.denominator * mzv2(w - r, r)
    return v


def _q(f: Fr):
    return mpf(f.numerator) / f.denominator


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestSommerfeldCoefficients:
    """sec:euler_cancel "Closed forms for c_p(n)" and the c_4, c_5 displays."""

    @pytest.mark.parametrize("p", [2, 3, 4, 5])
    def test_printed_cp_equals_exact_series(self, p: int) -> None:
        """Printed c_p(n) == exact series coefficient (R1) for n = 1..20.

        Rejects: any wrong printed coefficient (e.g. c_3 with +1/2 H^(2),
        c_4 with 41n-40 -> 41n-41, c_5 with -15/32 -> +15/32), and a
        degeneracy convention d(n,j) = 2j+1 throughout."""
        c = PRINTED_C[p]
        for n in range(1, N_SERIES + 1):
            assert c(n) == ser(n)[p], (p, n, c(n), ser(n)[p])

    def test_c1_is_minus_one(self) -> None:
        """c_1(n) = -1 at every n (sum_j d(n,j) / (2 n^2) = 2n^2/(2n^2)): the
        series' normalization anchor.  Rejects: a degeneracy-weight error
        (any d(n,j) that does not sum to 2n^2)."""
        for n in range(1, N_SERIES + 1):
            assert ser(n)[1] == -1

    @pytest.mark.parametrize("p", [2, 3, 4, 5])
    def test_fit_reproduces_printed_coefficients(self, p: int) -> None:
        """R2's exact fit equals the printed coefficients term by term.
        Rejects: a printed closed form that agrees with the series only at
        the sampled n by accident (the fit basis is checked on extra n)."""
        coef = fit(p)
        w = 2 * p - 1
        printed = {
            2: [-Fr(5, 4), Fr(1)],
            3: [Fr(19, 8), -Fr(20, 8), -Fr(3, 2), -Fr(1, 2)],
            4: [-Fr(205, 64), Fr(200, 64), Fr(15, 4), -Fr(1, 4), -Fr(3, 4),
                -Fr(1, 4)],
            5: [Fr(497, 128), -Fr(504, 128), -Fr(105, 16), Fr(45, 16),
                Fr(5, 4), -Fr(3, 8), -Fr(15, 32), -Fr(5, 32)],
        }[p]
        assert len(coef) == 2 + (2 * p - 4) == len(printed)
        assert coef == printed, (p, w, coef)

    def test_c6_has_no_H3_term(self) -> None:
        """sec. "D_6 at weight 11": "seven nontrivial sums (r = 1..8,
        r != 3) contribute through the beta_r coefficients of c_6(n)".
        Rejects: a nonzero beta_3 in c_6 (which would put S_{3,8} into D_6
        and change the zeta(3) zeta(8) column)."""
        coef = fit(6)
        betas = coef[2:]
        assert len(betas) == 8
        assert betas[2] == 0
        assert all(b != 0 for i, b in enumerate(betas) if i != 2)


class TestLowOrderDirichletSums:
    """D_2 = -(5/4) zeta(2) + zeta(3); D_3 = (19/8) zeta(4) - (11/4) zeta(5)."""

    def test_D2(self) -> None:
        """Rejects: D_2 with +zeta(3) -> -zeta(3) or -5/4 -> -5/2."""
        mp.dps = DPS
        assert abs(D_numeric(2) - (-mpf(5) / 4 * zeta(2) + zeta(3))) < mpf(10) ** -50

    def test_D3(self) -> None:
        """Rejects: D_3 = 19/8 zeta(4) - 5/2 zeta(5) (the rational part alone,
        dropping the Euler-sum contributions) and any surviving
        zeta(2) zeta(3) product."""
        mp.dps = DPS
        d3 = D_numeric(3)
        assert abs(d3 - (mpf(19) / 8 * zeta(4) - mpf(11) / 4 * zeta(5))) < mpf(10) ** -50
        assert abs(d3 - (mpf(19) / 8 * zeta(4) - mpf(5) / 2 * zeta(5))) > mpf(10) ** -3


class TestZeta3Complementarity:
    """Theorem thm:zeta3_complement (eq:shell_complement,
    eq:dirichlet_complement)."""

    @staticmethod
    def shell_sum(n: int) -> Fr:
        """S(n) = sum_j d(n,j) (n/(j+1/2) - 3/4), direct over j."""
        tot = Fr(0)
        for kap in range(1, n + 1):
            d = 4 * kap if kap < n else 2 * kap
            tot += d * (Fr(n, kap) - Fr(3, 4))
        return tot

    def test_shell_sum_closed_form_and_c2(self) -> None:
        """S(n) = n(5n-4)/2 and c_2(n) = -S(n)/(2n^4) (c_2 from R1).
        Rejects: S(n) = n(5n-4) (factor 2) or c_2 = -S/n^4."""
        for n in range(1, N_SERIES + 1):
            S = self.shell_sum(n)
            assert S == Fr(n * (5 * n - 4), 2)
            assert ser(n)[2] == -S / (2 * n ** 4)

    def test_shell_complement(self) -> None:
        """eq:shell_complement: S(n) + 2n(n+1) = 9n^2/2 = |lambda_0|^2 * 2n^2.
        Rejects: the single-chirality degeneracy 2(n+1)(n+2) (the CH
        convention g_n, which does not cancel the linear term) and a
        constant other than 9/2."""
        for n in range(1, N_SERIES + 1):
            S = self.shell_sum(n)
            assert S + 2 * n * (n + 1) == Fr(9 * n * n, 2)
            assert Fr(9 * n * n, 2) == Fr(3, 2) ** 2 * (2 * n * n)
            assert S + 2 * (n + 1) * (n + 2) != Fr(9 * n * n, 2)

    @pytest.mark.parametrize("s", [4, 5, 6, 7])
    def test_dirichlet_complement(self, s: int) -> None:
        """eq:dirichlet_complement: D_S(s) + D_g(s) = (9/2) zeta(s-2), with
        each side summed directly (nsum over n of the j-summed S(n)).
        Rejects: a surviving zeta(s-1) term."""
        mp.dps = 40
        DS = nsum(lambda n: (mpf(5) / 2 * n * n - 2 * n) / n ** s, [1, inf])
        Dg = nsum(lambda n: 2 * n * (n + 1) / n ** s, [1, inf])
        assert abs(DS + Dg - mpf(9) / 2 * zeta(s - 2)) < mpf(10) ** -30
        assert abs(DS - (mpf(5) / 2 * zeta(s - 2) - 2 * zeta(s - 1))) < mpf(10) ** -30
        # the closed form of S(n) used in the nsum is the direct j-sum
        for n in range(1, 8):
            assert self.shell_sum(n) == Fr(5, 2) * n * n - 2 * n

    def test_s4_evaluation(self) -> None:
        """D_S(4) = (5/2) pi^2/6 - 2 zeta(3); D_g(4) = 2 zeta(2) + 2 zeta(3);
        sum = 3 pi^2 / 4.  Rejects: D_S(4) = (5/2) pi^2/6 + 2 zeta(3)."""
        mp.dps = 40
        DS = nsum(lambda n: (mpf(5) / 2 * n * n - 2 * n) / n ** 4, [1, inf])
        Dg = nsum(lambda n: 2 * n * (n + 1) / n ** 4, [1, inf])
        assert abs(DS - (mpf(5) / 2 * pi ** 2 / 6 - 2 * zeta(3))) < mpf(10) ** -30
        assert abs(Dg - (2 * zeta(2) + 2 * zeta(3))) < mpf(10) ** -30
        assert abs(DS + Dg - 3 * pi ** 2 / 4) < mpf(10) ** -30


# Printed S_{r,s} tables: row -> coefficients in the weight-w basis
# w = 7: [z7, z2 z5, z3 z4]
# w = 9: [z9, z2 z7, z3 z6, z4 z5]
# w = 11: [z11, z2 z9, z3 z8, z4 z7, z5 z6]
S_TABLE: Dict[int, Dict[Tuple[int, int], List[Fr]]] = {
    7: {
        (1, 6): [Fr(4), Fr(-1), Fr(-1)],
        (2, 5): [Fr(-10), Fr(5), Fr(2)],
        (3, 4): [Fr(18), Fr(-10), Fr(0)],
        (4, 3): [Fr(-17), Fr(10), Fr(1)],
        (5, 2): [Fr(11), Fr(-4), Fr(-2)],
    },
    9: {
        (2, 7): [Fr(-35, 2), Fr(7), Fr(2), Fr(4)],
        (3, 6): [Fr(85, 2), Fr(-21), Fr(0), Fr(-6)],
        (4, 5): [Fr(-125, 2), Fr(35), Fr(0), Fr(5)],
        (5, 4): [Fr(127, 2), Fr(-35), Fr(0), Fr(-4)],
        (6, 3): [Fr(-83, 2), Fr(21), Fr(1), Fr(6)],
        (7, 2): [Fr(37, 2), Fr(-6), Fr(-2), Fr(-4)],
    },
    11: {
        (1, 10): [Fr(6), Fr(-1), Fr(-1), Fr(-1), Fr(-1)],
        (2, 9): [Fr(-27), Fr(9), Fr(2), Fr(6), Fr(4)],
        (4, 7): [Fr(-329, 2), Fr(84), Fr(0), Fr(21), Fr(4)],
        (5, 6): [Fr(463, 2), Fr(-126), Fr(0), Fr(-21), Fr(0)],
        (6, 5): [Fr(-461, 2), Fr(126), Fr(0), Fr(21), Fr(1)],
        (7, 4): [Fr(331, 2), Fr(-84), Fr(0), Fr(-20), Fr(-4)],
        (8, 3): [Fr(-82), Fr(36), Fr(1), Fr(15), Fr(6)],
    },
}


def _basis_values(w: int):
    out = [zeta(w)]
    for a in range(2, (w - 1) // 2 + 1):
        out.append(zeta(a) * zeta(w - a))
    return out


def _euler_S1(m: int) -> List[Fr]:
    """Euler: 2 S_{1,m} = (m+2) zeta(m+1) - sum_{k=1}^{m-2} zeta(m-k) zeta(k+1),
    as a vector in the weight-(m+1) basis [z_w, z2 z_{w-2}, z3 z_{w-3}, ...]."""
    w = m + 1
    vec = [Fr(0)] * ((w - 1) // 2)
    vec[0] = Fr(m + 2, 2)
    for k in range(1, m - 1):
        a, b = sorted((m - k, k + 1))
        vec[a - 1] -= Fr(1, 2)
    return vec


class TestEulerSumTables:
    """The S_{r,s} = sum_n H_n^(r) / n^s tables at weights 7, 9, 11."""

    @pytest.mark.parametrize("w,rs", [(w, rs) for w in S_TABLE for rs in S_TABLE[w]])
    def test_table_row_numeric(self, w: int, rs: Tuple[int, int]) -> None:
        """Each printed row == zeta(s, r) + zeta(r+s) by R3, to 1e-45.
        Rejects: any single wrong table entry (e.g. S_{3,4} zeta(3)zeta(4)
        coefficient 0 -> 1, S_{8,3} zeta(3)zeta(8) coefficient 1 -> 0)."""
        mp.dps = DPS
        r, s = rs
        val = mzv2(s, r) + zeta(r + s)
        vec = S_TABLE[w][rs]
        pred = sum(_q(c) * b for c, b in zip(vec, _basis_values(w)))
        assert abs(val - pred) < mpf(10) ** -45, (rs, val, pred)

    def test_euler_S1_formula_matches_table(self) -> None:
        """Euler's S_{1,m} formula (route for S_{1,8}) reproduces the
        printed S_{1,6} and S_{1,10} rows exactly.  Rejects: a mis-stated
        Euler formula used in R4."""
        assert _euler_S1(6) == S_TABLE[7][(1, 6)]
        assert _euler_S1(10) == S_TABLE[11][(1, 10)]


# Printed D_p decompositions: [z_{w-1} coefficient, then the weight-w basis]
D_PRINTED = {
    4: [-Fr(205, 64), Fr(71, 8), Fr(0), -Fr(9, 2)],
    5: [Fr(497, 128), -Fr(467, 16), Fr(0), Fr(385, 32), Fr(75, 8)],
    6: [-Fr(2289, 512), Fr(1589, 16), Fr(0), -Fr(1617, 64), -Fr(1785, 64),
        -Fr(2065, 64)],
}


def _D_exact_from_tables(p: int) -> List[Fr]:
    """R4: [coefficient of zeta(w-1), weight-w basis vector] of D_p, built
    from R2's coefficients and the printed S tables (zeta(s, r) =
    S_{r,s} - zeta(w))."""
    w = 2 * p - 1
    coef = fit(p)
    nb = (w - 1) // 2
    vec = [Fr(0)] * nb
    vec[0] += coef[1]
    for r in range(1, 2 * p - 3):
        b = coef[r + 1]
        if b == 0:
            continue
        s = w - r
        row = S_TABLE[w].get((r, s))
        if row is None:
            assert r == 1, ("table row missing", r, s)
            row = _euler_S1(s)
        row = list(row)
        row[0] -= 1  # zeta(s, r) = S_{r,s} - zeta(w)
        vec = [v + b * x for v, x in zip(vec, row)]
    return [coef[0]] + vec


class TestDirichletDecompositions:
    """eq:D4_decomp, eq:D5_decomp, eq:D6_decomp, eq:D4_three_tier,
    obs:product_rule."""

    @pytest.mark.parametrize("p", [4, 5, 6])
    def test_printed_decomposition_numeric(self, p: int) -> None:
        """Printed D_p == sum_n c_p(n) by R1 -> R2 -> R3, to 1e-45 (the old
        pin only evaluated the printed closed form against its own frozen
        decimal).  Rejects: any wrong printed coefficient, e.g. D_4 with
        -9/2 -> -7/2, D_6 with -1617/64 -> -1615/64."""
        mp.dps = DPS
        w = 2 * p - 1
        vec = D_PRINTED[p]
        pred = _q(vec[0]) * zeta(w - 1) + sum(
            _q(c) * b for c, b in zip(vec[1:], _basis_values(w)))
        assert abs(D_numeric(p) - pred) < mpf(10) ** -45, (p, D_numeric(p), pred)

    @pytest.mark.parametrize("p", [4, 5, 6])
    def test_printed_decomposition_exact_from_tables(self, p: int) -> None:
        """R4: the printed S tables + the exact c_p coefficients give the
        printed D_p as an exact rational vector, so the zeta(2) zeta(2p-3)
        cancellation is exact ("cancels algebraically").  Rejects: a
        D_p that agrees numerically only to the PSLQ residual."""
        assert _D_exact_from_tables(p) == D_PRINTED[p]

    @pytest.mark.parametrize("p", [2, 3, 4, 5, 6])
    def test_product_survival_rule(self, p: int) -> None:
        """obs:product_rule through p = 6: zeta(2) zeta(2p-3) has coefficient
        0, and exactly max(0, floor((2p-5)/2)) products survive, each
        zeta(a) zeta(b) with a odd >= 3, b even >= 4.  Rejects: a surviving
        zeta(2) product, i.e. every product surviving, floor((2p-3)/2)
        (4 at p = 6, not 3; 3 at p = 5, not 2)."""
        if p <= 3:
            survivors = []  # D_2, D_3 are single-zeta (TestLowOrder)
        else:
            vec = _D_exact_from_tables(p)[1:]
            w = 2 * p - 1
            assert vec[1] == 0  # zeta(2) zeta(w-2)
            survivors = [a for a, c in zip(range(2, (w - 1) // 2 + 1), vec[1:])
                         if c != 0]
            for a in survivors:
                b = w - a
                odd, even = (a, b) if a % 2 else (b, a)
                assert odd >= 3 and even >= 4
        assert len(survivors) == max(0, (2 * p - 5) // 2)

    def test_D4_three_tier(self) -> None:
        """eq:D4_three_tier: -(205/64) zeta(6) = -41 pi^6/12096 and
        -(9/2) zeta(3) zeta(4) = -(pi^4/20) zeta(3), exactly (zeta(6) =
        pi^6/945, zeta(4) = pi^4/90).  Rejects: 41/12096 -> 41/12095."""
        assert Fr(205, 64) * Fr(1, 945) == Fr(41, 12096)
        assert Fr(9, 2) * Fr(1, 90) == Fr(1, 20)
        mp.dps = 40
        assert abs(zeta(6) - pi ** 6 / 945) < mpf(10) ** -35
        assert abs(zeta(4) - pi ** 4 / 90) < mpf(10) ** -35
