"""Helper for tests/test_paper58_pachucki_exchange.py -- Pachucki's ODE at the exchange point, exact.

Source (prior art): K. Pachucki, Phys. Rev. A 80, 032520 (2009), arXiv:0905.1635.  Every formula below is
transcribed from the arXiv TeX source; the label noted at each is the arXiv equation label (the journal numbering is
shifted by one or more).  Nothing here reads or imports debug/.

  master integral (label 01)
      f(r) = int d3r1/4pi int d3r2/4pi  e^{-u3 r1A}/r1A  e^{-u2 r1B}/r1B  e^{-w2 r2A}/r2A  e^{-w3 r2B}/r2B  r/r12
  ODE (label 28)                r f'' + f' - p^2 r f + F = 0
  p^2 (label 20), a+-, b+- (labels 22, 23), F (label 30), F_i^- (labels 31-34)
  exchange point (Sec. "Special case: exchanged atomic orbitals"):  u3 = w2 = a,  u2 = w3 = b
  f0 = f(exchange point) = int_r^oo F^(-1)(s)/s ds (label 59),  F^(-1) (label 60),  dF^(-1)/dr = F.

Derivative equations.  Apply d_S (S a set of exponents) to label 28 and set the exponents to the exchange point,
where p^2 = 0:

      (r f_S')' = -F_S + r * sum_{0 != U subset S} (d_U p^2) f_{S \\ U}                              (*)

Function algebra.  A function is a dict  key -> {power k: Fraction}.  Keys

      ('E', lam, j)     r^k  LC^j  e^{lam r}          LC = ln((a+b)/(b-a)), kept as an INDETERMINATE (j = 0, 1)
      ('Ei', lam, mu)   r^k  e^{lam r}  Ei(mu r)
      ('f',)            r^k  f0(r)                    with the exact rule  f0' = -F^(-1)/r  (label 59)

All arithmetic is exact (fractions.Fraction; python-flint for the linear algebra).
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from functools import lru_cache
from math import gcd
from typing import Dict, List, Sequence, Tuple

import flint
import sympy as sp

r = sp.Symbol("r", positive=True)
u2, u3, w2, w3 = sp.symbols("u2 u3 w2 w3", positive=True)
THETA = (u3, u2, w2, w3)            # order: exponent of r1A, r1B, r2A, r2B
LCS = sp.Symbol("LC")

Func = Dict[tuple, Dict[int, Fraction]]
SUBSETS = [al for al in itertools.product((0, 1), repeat=4)]
PREFACTOR = 16                      # K = PREFACTOR a^3 b^3 / R * d_u3 d_u2 d_w2 d_w3 f
PMIN, PMAX = -6, 6                  # Laurent window of the ansatz


# ------------------------------------------------------------------ Pachucki's formulas (transcription)
def p2_expr() -> sp.Expr:
    """p^2, arXiv label 20."""
    return ((u2**2 - u3**2 + w2**2 - w3**2) * (u2**2 * w2**2 - u3**2 * w3**2)
            / ((u2**2 - u3**2) * (w2**2 - w3**2)))


def ln_argument() -> sp.Expr:
    """The argument of ln|.| in arXiv labels 33 and 34 (the same in both)."""
    return ((u2 + w3 - u3 - w2) * (u2 + u3) * (w2 + w3)
            / ((u2 + u3 + w2 + w3) * (u2 - u3) * (w2 - w3)))


def source_F() -> sp.Expr:
    """F(r), arXiv label 30, with a+-, b+- of labels 22, 23 and F_i^- of labels 31-34.

    ln|x| is entered as ln(-x): x < 0 in a neighbourhood of the exchange point with a < b (checked by the test).
    """
    Ei, E = sp.Ei, sp.exp
    S = u2 + u3 + w2 + w3
    lnabs = sp.log(-ln_argument())
    F1m = E(r * (u3 - w2)) * Ei(-r * (u2 + u3)) - E(r * (w2 - u3)) * Ei(-r * (w2 + w3))             # label 31
    F2m = E(r * (w3 - u2)) * Ei(-r * (w2 + w3)) - E(r * (u2 - w3)) * Ei(-r * (u2 + u3))             # label 32
    F3m = (E(-r * (u3 + w2)) * (lnabs - Ei(r * (u3 + w2 - u2 - w3)) + Ei(r * (w2 - w3)) + Ei(r * (u3 - u2)))
           - E(r * (u3 + w2)) * Ei(-r * S))                                                         # label 33
    F4m = (E(-r * (u2 + w3)) * (lnabs - Ei(r * (u2 + w3 - u3 - w2)) + Ei(r * (w3 - w2)) + Ei(r * (u2 - u3)))
           - E(r * (u2 + w3)) * Ei(-r * S))                                                         # label 34
    ap = w2 / (w2**2 - w3**2) + u3 / (u2**2 - u3**2)                                                # label 22
    am = w2 / (w2**2 - w3**2) - u3 / (u2**2 - u3**2)
    bp = u2 / (u2**2 - u3**2) + w3 / (w2**2 - w3**2)                                                # label 23
    bm = u2 / (u2**2 - u3**2) - w3 / (w2**2 - w3**2)
    return sp.Rational(1, 2) * (ap * F1m + bp * F2m + am * F3m + bm * F4m)                          # label 30


def exchange_point(a, b) -> dict:
    return {u3: a, w2: a, u2: b, w3: b}


@lru_cache(maxsize=None)
def _dsym(al: tuple) -> Tuple[sp.Expr, sp.Expr]:
    """(d_al F, d_al p^2), symbolic in the four exponents; one variable at a time from a cached parent."""
    if sum(al) == 0:
        return source_F(), p2_expr()
    i = max(j for j in range(4) if al[j])
    par = tuple(x - (1 if j == i else 0) for j, x in enumerate(al))
    e, q = _dsym(par)
    return sp.diff(e, THETA[i]), sp.diff(q, THETA[i])


# ------------------------------------------------------------------ function algebra
def _frac(x) -> Fraction:
    x = sp.Rational(x)
    return Fraction(int(x.p), int(x.q))


def clean(A: Func) -> Func:
    out = {}
    for k, d in A.items():
        dd = {p: c for p, c in d.items() if c != 0}
        if dd:
            out[k] = dd
    return out


def fadd(A: Func, B: Func, s=1) -> Func:
    out = {k: dict(v) for k, v in A.items()}
    for k, d in B.items():
        o = out.setdefault(k, {})
        for p, c in d.items():
            o[p] = o.get(p, 0) + s * c
    return clean(out)


def fscale(A: Func, c, shift: int = 0) -> Func:
    return clean({k: {p + shift: c * v for p, v in d.items()} for k, d in A.items()})


class Algebra:
    """The function algebra at one pair (a, b), a < b rational."""

    def __init__(self, a, b):
        self.a, self.b = Fraction(a), Fraction(b)
        assert 0 < self.a < self.b
        self.FM1 = self.fm1()

    # -- arXiv label 60 (u = a, w = b); ln|(w-u)/(w+u)| = -LC
    def fm1(self) -> Func:
        a, b = self.a, self.b
        c = 1 / (2 * (a * a - b * b))
        return clean({
            ("Ei", 2 * b, -2 * (a + b)): {0: c},
            ("Ei", 2 * a, -2 * (a + b)): {0: -c},
            ("Ei", -2 * b, 2 * (b - a)): {0: -c},
            ("Ei", -2 * b, (b - a)): {0: 2 * c},
            ("E", -2 * b, 1): {0: c},
            ("Ei", -2 * a, 2 * (a - b)): {0: c},
            ("Ei", -2 * a, (a - b)): {0: -2 * c},
            ("E", -2 * a, 1): {0: -c},
        })

    def keys11(self) -> List[tuple]:
        """The 11 functions (LC-free keys)."""
        a, b = self.a, self.b
        return [("f",), ("Ei", Fraction(0), -(a + b)), ("E", -(a + b), 0), ("E", -2 * a, 0), ("E", -2 * b, 0),
                ("Ei", 2 * a, -2 * (a + b)), ("Ei", 2 * b, -2 * (a + b)),
                ("Ei", -2 * a, (a - b)), ("Ei", -2 * a, 2 * (a - b)),
                ("Ei", -2 * b, (b - a)), ("Ei", -2 * b, 2 * (b - a))]

    def keys_lc(self) -> List[tuple]:
        """The LC-multiples that occur (LC is a constant: these are the same functions over Q(LC))."""
        return [("E", -2 * self.a, 1), ("E", -2 * self.b, 1)]

    def diff(self, A: Func) -> Func:
        out: Func = {}

        def add(k, p, c):
            o = out.setdefault(k, {})
            o[p] = o.get(p, 0) + c

        for k, d in A.items():
            for p, c in d.items():
                if p != 0:
                    add(k, p - 1, p * c)
                if k[0] == "E":
                    if k[1] != 0:
                        add(k, p, k[1] * c)
                elif k[0] == "Ei":
                    lam, mu = k[1], k[2]
                    if lam != 0:
                        add(k, p, lam * c)
                    add(("E", lam + mu, 0), p - 1, c)          # d Ei(mu r) = e^{mu r} / r
                elif k[0] == "f":
                    for kk, dd in self.FM1.items():            # f0' = -F^(-1)/r, label 59
                        for pp, cc in dd.items():
                            add(kk, p + pp - 1, -c * cc)
                else:
                    raise KeyError(k)
        return clean(out)

    def L0(self, A: Func) -> Func:
        """(r A')'"""
        return self.diff(fscale(self.diff(A), 1, 1))

    # -- sympy expression at the exchange point -> Func
    def parse(self, expr: sp.Expr) -> Func:
        a, b = self.a, self.b
        lcval = sp.Rational(str((a + b) / (b - a)))
        expr = expr.xreplace({sp.log(lcval): LCS})
        assert not expr.has(sp.log), "unexpected logarithm"
        expr = sp.expand(sp.powsimp(sp.expand(expr), combine="exp"))
        out: Func = {}
        for t in sp.Add.make_args(expr):
            coeff, k, lam, mu, j = sp.Integer(1), 0, sp.Integer(0), None, 0
            for fac in sp.Mul.make_args(t):
                if fac == LCS:
                    j += 1
                elif not fac.has(r):
                    assert not fac.has(LCS)
                    coeff *= fac
                elif fac == r:
                    k += 1
                elif fac.is_Pow and fac.base == r:
                    k += int(fac.exp)
                elif isinstance(fac, sp.exp):
                    lam += sp.expand(fac.args[0] / r)
                elif isinstance(fac, sp.Ei):
                    assert mu is None
                    mu = sp.expand(fac.args[0] / r)
                else:
                    raise ValueError(f"unparsed factor {fac} in {t}")
            assert not sp.sympify(lam).has(r)
            if mu is None:
                key = ("E", _frac(lam), j)
            else:
                assert j == 0 and not mu.has(r)
                key = ("Ei", _frac(lam), _frac(mu))
            o = out.setdefault(key, {})
            o[k] = o.get(k, 0) + _frac(coeff)
        return clean(out)

    def sources(self) -> Tuple[Dict[tuple, Func], Dict[tuple, Fraction]]:
        """(F_S as Func, d_S p^2 as Fraction) at the exchange point, for the 16 subsets S."""
        sub = exchange_point(sp.Rational(str(self.a)), sp.Rational(str(self.b)))
        Fs, P2 = {}, {}
        for al in SUBSETS:
            e, q = _dsym(al)
            Fs[al] = self.parse(e.subs(sub))
            P2[al] = _frac(sp.nsimplify(q.subs(sub)))
        return Fs, P2

    def rhs(self, al: tuple, Fs, P2, fS) -> Func:
        """Right-hand side of (*) for the subset al, given the lower f_T."""
        src = fscale(Fs[al], -1)
        for be in SUBSETS:
            if any(y > x for x, y in zip(al, be)) or be == al:
                continue
            U = tuple(x - y for x, y in zip(al, be))
            if P2[U] != 0:
                src = fadd(src, fscale(fS[be], P2[U], 1))      # the p^2 coupling term
        return src

    def solve_L0(self, src: Func, keys: Sequence[tuple]) -> Tuple[Func, int]:
        """u = sum_{key in keys, PMIN <= k <= PMAX} x r^k B_key with (r u')' = src, exactly.

        Returns (u, number of free ansatz columns).  Raises RuntimeError if the system is inconsistent.
        """
        cols = [(k, p) for k in keys for p in range(PMIN, PMAX + 1)]
        imgs = [self.L0({k: {p: Fraction(1)}}) for (k, p) in cols]
        rows = sorted({(k, p) for im in imgs + [src] for k, d in im.items() for p in d}, key=str)
        ridx = {rw: i for i, rw in enumerate(rows)}
        n = len(cols)
        M = [[Fraction(0)] * (n + 1) for _ in rows]
        for j, im in enumerate(imgs):
            for k, d in im.items():
                for p, c in d.items():
                    M[ridx[(k, p)]][j] = c
        for k, d in src.items():
            for p, c in d.items():
                M[ridx[(k, p)]][n] = c
        Rm, piv = rref_q(M)
        if n in piv:
            raise RuntimeError("no solution in the ansatz")
        sol: Func = {}
        for i, pc in enumerate(piv):
            val = Rm[i][n]
            if val != 0:
                k, p = cols[pc]
                sol.setdefault(k, {})[p] = val
        return clean(sol), n - len(piv)

    def solve_all(self, keys: Sequence[tuple] = None):
        """Solve (*) for the 15 nonempty subsets, lowest first.  Returns (fS, free-column counts, Fs, P2)."""
        if keys is None:
            keys = self.keys11() + self.keys_lc()
        Fs, P2 = self.sources()
        fS = {(0, 0, 0, 0): {("f",): {0: Fraction(1)}}}
        free = {}
        for al in sorted(SUBSETS, key=lambda t: (sum(t), t)):
            if sum(al) == 0:
                continue
            fS[al], free[al] = self.solve_L0(self.rhs(al, Fs, P2, fS), keys)
        return fS, free, Fs, P2

    def K(self, fS) -> Func:
        """K(R) = (1s_A^a 1s_B^b | 1s_A^a 1s_B^b) = PREFACTOR a^3 b^3 / R * f_{1111}(R)."""
        return fscale(fS[(1, 1, 1, 1)], PREFACTOR * self.a**3 * self.b**3, -1)

    def derivatives(self, A: Func, n: int) -> List[Func]:
        out = [A]
        for _ in range(n):
            out.append(self.diff(out[-1]))
        return out


# ------------------------------------------------------------------ exact linear algebra (python-flint)
def _to_fmpz_rows(M: List[List[Fraction]]) -> List[List[int]]:
    out = []
    for row in M:
        den = 1
        for c in row:
            d = Fraction(c).denominator
            den = den * d // gcd(den, d)
        out.append([int(Fraction(c) * den) for c in row])
    return out


def rref_q(M: List[List[Fraction]]):
    """Reduced row echelon form over Q.  Returns (rows as Fractions, pivot columns)."""
    Z = _to_fmpz_rows(M)
    A = flint.fmpz_mat(len(Z), len(Z[0]), [x for row in Z for x in row])
    Rz, den, rank = A.rref()
    den = int(den)
    out, piv = [], []
    for i in range(rank):
        row = [int(Rz[i, j]) for j in range(len(Z[0]))]
        pc = next(j for j, x in enumerate(row) if x != 0)
        piv.append(pc)
        out.append([Fraction(x, den) for x in row])
    return out, piv


def rank_q(M: List[List[Fraction]]) -> int:
    Z = _to_fmpz_rows(M)
    return int(flint.fmpz_mat(len(Z), len(Z[0]), [x for row in Z for x in row]).rank())


def nullspace_q(M: List[List[Fraction]]) -> List[List[int]]:
    """Basis of the right null space, as integer vectors."""
    Z = _to_fmpz_rows(M)
    A = flint.fmpz_mat(len(Z), len(Z[0]), [x for row in Z for x in row])
    N, nullity = A.nullspace()
    return [[int(N[i, j]) for i in range(len(Z[0]))] for j in range(nullity)]


# ------------------------------------------------------------------ operators  L = sum_i c_i(R) D^i
def all_keys(ders: Sequence[Func]) -> List[tuple]:
    return sorted({k for A in ders for k in A}, key=str)


def coordinate_matrix_at(ders: Sequence[Func], keys: Sequence[tuple], x: Fraction) -> List[List[Fraction]]:
    """Rows = derivatives, columns = keys, entries = the Laurent coordinate evaluated at R = x."""
    return [[sum((c * x**p for p, c in A.get(k, {}).items()), Fraction(0)) for k in keys] for A in ders]


def operator_system(ders: Sequence[Func], order: int, degree: int) -> List[List[Fraction]]:
    """Linear system for the unknowns c_ij (i <= order, j <= degree; column i*(degree+1)+j):
    for every key and every power of R, the coordinate of sum_i c_i(R) D^i K vanishes."""
    rows: Dict[tuple, List[Fraction]] = {}
    ncol = (order + 1) * (degree + 1)
    for i in range(order + 1):
        for k, d in ders[i].items():
            for p, c in d.items():
                for j in range(degree + 1):
                    rows.setdefault((k, p + j), [Fraction(0)] * ncol)[i * (degree + 1) + j] += c
    return [rows[kk] for kk in sorted(rows, key=str)]


def operators_annihilating(ders: Sequence[Func], order: int, degree: int) -> List[List[List[int]]]:
    """All operators of order <= `order` and coefficient degree <= `degree` that annihilate ders[0] (a basis),
    each as [c_0, ..., c_order], c_i an ascending integer coefficient list, content removed."""
    out = []
    for v in nullspace_q(operator_system(ders, order, degree)):
        g = 0
        for x in v:
            g = gcd(g, x)
        v = [x // g for x in v]
        out.append([v[i * (degree + 1):(i + 1) * (degree + 1)] for i in range(order + 1)])
    return out


def apply_operator_mod(coeffs: Sequence[Sequence[int]], ders: Sequence[Func], prime: int) -> Dict[tuple, int]:
    """Coordinates of sum_i c_i(R) D^i K modulo `prime` (only the nonzero ones are returned)."""
    acc: Dict[tuple, int] = {}
    for i, ci in enumerate(coeffs):
        for k, d in ders[i].items():
            for p, c in d.items():
                c = Fraction(c)
                assert c.denominator % prime != 0
                cm = c.numerator * pow(c.denominator, -1, prime) % prime
                for j, x in enumerate(ci):
                    if x:
                        acc[(k, p + j)] = (acc.get((k, p + j), 0) + cm * int(x)) % prime
    return {kk: v for kk, v in acc.items() if v}


def trim(c: Sequence[int]) -> List[int]:
    c = list(c)
    while c and c[-1] == 0:
        c.pop()
    return c


def indicial_roots_at_zero(op: Sequence[Sequence[int]]) -> dict:
    """Roots (with multiplicity) of the indicial polynomial of L at R = 0."""
    k = sp.Symbol("k")
    val = [next((j for j, x in enumerate(c) if x), None) for c in op]
    wmin = min(v - i for i, v in enumerate(val) if v is not None)
    ind = sum((op[i][wmin + i] if 0 <= wmin + i < len(op[i]) else 0) * sp.ff(k, i) for i in range(len(op)))
    P = sp.Poly(sp.expand(ind), k)
    return P.degree(), sp.roots(P)


def rates_at_infinity(op: Sequence[Sequence[int]]) -> dict:
    """Roots (with multiplicity) of the characteristic polynomial at infinity, sum_i [R^dmax] c_i x^i."""
    x = sp.Symbol("x")
    dmax = max(len(trim(c)) - 1 for c in op)
    ch = sum((op[i][dmax] if dmax < len(op[i]) else 0) * x**i for i in range(len(op)))
    P = sp.Poly(ch, x)
    return P.degree(), sp.roots(P)


# ------------------------------------------------------------------ numerics (mpmath)
def evaluate(alg: Algebra, A: Func, x, f0=None):
    """A(x) in mpmath at the ambient precision; LC = ln((a+b)/(b-a)); f0 must be supplied if A has an f key."""
    import mpmath as mp

    def q(v):
        v = Fraction(v)
        return mp.mpf(v.numerator) / v.denominator

    lc = mp.log(q((alg.a + alg.b) / (alg.b - alg.a)))
    s = mp.mpf(0)
    for k, d in A.items():
        if k[0] == "f":
            base = f0
        elif k[0] == "E":
            base = mp.exp(q(k[1]) * x) * lc**k[2]
        else:
            base = mp.exp(q(k[1]) * x) * mp.ei(q(k[2]) * x)
        for p, c in d.items():
            s += q(c) * x**p * base
    return s


def f0_numeric(alg: Algebra, x):
    """f0(x) = int_x^oo F^(-1)(s)/s ds (label 59), the only quadrature."""
    import mpmath as mp
    return mp.quad(lambda s: evaluate(alg, alg.FM1, s) / s, [x, x + 1, x + 2, x + 4, x + 10, x + 30, mp.inf])
