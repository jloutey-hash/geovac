"""Helper for tests/test_paper58_minimal_order.py -- minimal order in R of Paper 58's two-centre closed forms.

Own code path, independent of tests/test_paper58_r_operators.py (`decompose`, `d_atoms`) and of debug/:

  * `Fn(expr)`: the closed form as a finite sum over ATOMS e^{cR} y with coefficients polynomial in the numerical
    constants (Euler gamma, logs of PRIMES, surds) over Q(R).  Atom labels are normalised the way Lemma GI states
    them: y = 1, y = ln R (kappa = 0), y = Ei(kappa R) with E1(kR) = -Ei(-kR).  Two spellings of one function
    (E1(2R) and -Ei(-2R)) therefore land on ONE label and their coefficients merge; logs of composite numbers are
    reduced to logs of primes before anything is split, so ln 6 - ln 2 - ln 3 is 0, not a free symbol.
  * the derivative calculus on atoms: (e^{cR} ln R)' brings e^{cR}/R, (e^{cR} Ei(kR))' brings e^{(c+k)R}/R
    (cross-checked numerically against sympy.diff in the test file);
  * LOWER bound "no ODE of order < r": an r x r minor of the coordinate matrix of f, ..., f^(r-1) at a rational R0,
    computed EXACTLY as a polynomial in the constant symbols and evaluated at the TRUE constants in interval
    arithmetic.  An interval that excludes 0 proves rank r over C(R) (given the atoms are independent over C(R):
    Lemma GI, a hand proof).  No independence of gamma, ln 2, ... is assumed;
  * UPPER bound "an ODE of order r exists": the GENERIC rank over Q(R) of the split vectors of f, ..., f^(r)
    (rows = (atom, monomial in the constants)).  A Q(R)-relation among split vectors is a relation among the true
    functions.  Never a rank at a point: a point rank only bounds the generic rank from below.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import mpmath as mp
import sympy as sp
from sympy.polys.matrices import DomainMatrix

from geovac.two_center_eri import R_s

R = R_s
FR = sp.QQ.frac_field(R)
RG = FR.gens[0]
R0_DEFAULT = sp.Rational(37, 59)

Label = Tuple[sp.Rational, str, sp.Rational]      # (c, '1' | 'y', kappa)


def _reduce_logs(expr: sp.Expr) -> sp.Expr:
    """Expand everything and write every log of a number as a sum of logs of primes."""
    e = sp.expand(sp.expand_func(expr))
    return sp.expand(sp.expand_log(e, force=True, factor=True))


def _is_number_constant(a: sp.Expr) -> bool:
    if a.free_symbols:
        return False
    if a is sp.EulerGamma or isinstance(a, sp.NumberSymbol):
        return True
    if isinstance(a, sp.log):
        return True
    if isinstance(a, sp.Pow) and not a.exp.is_integer:
        return True
    return False


class Fn:
    """f = sum over (label, monomial in constant symbols) of FR coefficient * monomial * atom."""

    def __init__(self, expr: sp.Expr, reduce_logs: bool = True):
        e = _reduce_logs(expr) if reduce_logs else sp.expand(sp.expand_log(sp.expand(sp.expand_func(expr)),
                                                                          force=True))
        consts = sorted({a for a in e.atoms(sp.log, sp.Pow, sp.NumberSymbol) if _is_number_constant(a)}
                        | ({sp.EulerGamma} if e.has(sp.EulerGamma) else set()), key=str)
        assert not e.has(sp.pi), "pi in a two-centre closed form: not expected"
        self.syms: List[sp.Symbol] = [sp.Symbol(f"k{i}") for i in range(len(consts))]
        self.true: Dict[sp.Symbol, sp.Expr] = dict(zip(self.syms, consts))
        e = sp.expand(e.xreplace(dict(zip(consts, self.syms))))
        acc: Dict[Label, sp.Expr] = {}
        for term in sp.Add.make_args(e):
            c, kind, k, coeff = sp.Integer(0), "1", sp.Integer(0), sp.Integer(1)
            for fct in sp.Mul.make_args(term):
                base, ex = fct.as_base_exp()
                if isinstance(fct, sp.exp):
                    arg = fct.args[0]
                    assert sp.diff(arg, R, 2) == 0 and arg.subs(R, 0) == 0, fct
                    c += sp.diff(arg, R)
                elif isinstance(base, sp.exp):
                    arg = base.args[0]
                    assert sp.diff(arg, R, 2) == 0 and arg.subs(R, 0) == 0, fct
                    c += sp.diff(arg, R) * ex
                elif isinstance(fct, sp.expint):
                    assert fct.args[0] == 1 and kind == "1", term
                    kk = sp.diff(fct.args[1], R)
                    assert kk > 0 and fct.args[1] == kk * R, fct
                    kind, k = "y", -kk                    # E1(kR) = -Ei(-kR)
                    coeff *= -1
                elif isinstance(fct, sp.Ei):
                    assert kind == "1", term
                    kk = sp.diff(fct.args[0], R)
                    assert kk != 0 and fct.args[0] == kk * R, fct
                    kind, k = "y", kk
                elif fct == sp.log(R):
                    assert kind == "1", term
                    kind, k = "y", sp.Integer(0)
                else:
                    assert not fct.has(sp.exp, sp.expint, sp.Ei, sp.log), fct
                    coeff *= fct
            lab = (sp.Rational(c), kind, sp.Rational(k))
            acc[lab] = acc.get(lab, 0) + coeff
        self.coords: Dict[Tuple[Label, tuple], object] = {}
        for lab, v in acc.items():
            num, den = sp.fraction(sp.together(v))
            assert not den.has(*self.syms) if self.syms else True, (lab, den)
            P = sp.Poly(sp.expand(num), *self.syms) if self.syms else None
            parts = {m: cf for m, cf in P.terms()} if P is not None else {(): sp.expand(num)}
            for m, cf in parts.items():
                el = FR.from_sympy(sp.together(cf / den))      # raises unless rational in R
                if el != 0:
                    key = (lab, tuple(m))
                    self.coords[key] = self.coords.get(key, FR.zero) + el
        self.coords = {k: v for k, v in self.coords.items() if v != 0}

    # ------------------------------------------------------------------ structure
    def labels(self) -> List[Label]:
        """Atoms carrying a coefficient that is nonzero AT THE TRUE CONSTANTS (not merely as symbols)."""
        out = []
        for lab in sorted({lab for (lab, _m) in self.coords}, key=str):
            if not coefficient_vanishes(self, lab):
                out.append(lab)
        return out

    def rates(self) -> List[sp.Rational]:
        """Distinct exponential rates of f's own atoms (the E-atoms' second rate c + kappa is NOT included)."""
        return sorted({lab[0] for lab in self.labels()})


def derivative(coords: Dict) -> Dict:
    out: Dict = {}

    def add(key, v):
        out[key] = out.get(key, FR.zero) + v

    for ((c, kind, k), m), a in coords.items():
        add(((c, kind, k), m), a.diff(RG) + FR.from_sympy(c) * a)
        if kind == "y":
            add(((c + k, "1", sp.Integer(0)), m), a / RG)
    return {key: v for key, v in out.items() if v != 0}


def derivatives(f: Fn, upto: int) -> List[Dict]:
    ds = [f.coords]
    for _ in range(upto):
        ds.append(derivative(ds[-1]))
    return ds


def closure(labels) -> List[Label]:
    """Smallest d/dR-closed atom set containing `labels`."""
    out = set(labels)
    for (c, kind, k) in list(labels):
        if kind == "y":
            out.add((c + k, "1", sp.Integer(0)))
    return sorted(out, key=str)


# ------------------------------------------------------------------ true constants, intervals
def _intervals(f: Fn, dps: int = 70):
    mp.iv.dps = dps
    out = {}
    for s, expr in f.true.items():
        v = sp.N(expr, dps + 15)
        with mp.workdps(dps + 15):                 # parse at full precision (mp.mpf(str) rounds to the context)
            x = mp.mpf(str(v))
            w = abs(x) * mp.mpf(10) ** (-(dps - 5))
            lo, hi = x - w, x + w
        out[s] = mp.iv.mpf([lo, hi])
    return out


def _iv_poly(poly: sp.Poly, ivals) -> "mp.iv.mpf":
    tot = mp.iv.mpf(0)
    for mono, cf in poly.terms():
        cf = sp.Rational(cf)
        term = mp.iv.mpf(int(cf.p)) / mp.iv.mpf(int(cf.q))
        for g, e in zip(poly.gens, mono):
            if e:
                term = term * ivals[g] ** int(e)
        tot = tot + term
    return tot


def _mono_expr(f: Fn, m: tuple) -> sp.Expr:
    out = sp.Integer(1)
    for s, e in zip(f.syms, m):
        out *= s ** e
    return out


def coefficient_vanishes(f: Fn, lab: Label, R0=sp.Rational(7, 5)) -> bool:
    """True iff the coefficient of atom `lab` evaluates to an interval containing 0 at two rational R.

    (A rational function of R that is nonzero vanishes at finitely many points; two points that both give an
    interval containing 0 are treated as 'vanishes'.  Used only to count rates; the order itself is certified by
    the minors below.)"""
    iv = _intervals(f)
    for x0 in (R0, sp.Rational(29, 13)):
        expr = sum(sp.Rational(FR.to_sympy(v).subs(R, x0)) * _mono_expr(f, m)
                   for (lb, m), v in f.coords.items() if lb == lab)
        expr = sp.expand(expr)
        if f.syms:
            val = _iv_poly(sp.Poly(expr, *f.syms), iv)
            if val.a > 0 or val.b < 0:
                return False
        elif expr != 0:
            return False
    return True


def true_matrix(f: Fn, ds: List[Dict], rows: List[Label], ncols: int, x0) -> List[List[sp.Expr]]:
    """Entries: polynomials over Q in f.syms (the coordinate of f^(j) on atom `row` at R = x0)."""
    T = []
    for lab in rows:
        row = []
        for j in range(ncols):
            tot = sp.Integer(0)
            for (lb, m), w in ds[j].items():
                if lb == lab:
                    tot += sp.Rational(FR.to_sympy(w).subs(R, x0)) * _mono_expr(f, m)
            row.append(sp.expand(tot))
        T.append(row)
    return T


class Minor:
    def __init__(self, rows, poly, interval, value):
        self.rows, self.poly, self.interval, self.value = rows, poly, interval, value

    @property
    def nonzero(self) -> bool:
        return bool(self.interval.a > 0 or self.interval.b < 0)

    def depends_on(self, f: Fn) -> List[str]:
        """Names of gamma / log constants the minor's polynomial actually contains."""
        out = []
        if self.poly is None:
            return out
        for g in self.poly.free_symbols:
            t = f.true.get(g)
            if t is not None and (t is sp.EulerGamma or isinstance(t, sp.log)):
                out.append(str(t))
        return sorted(out)


def lower_minor(f: Fn, r: int, x0=R0_DEFAULT, rows: Optional[List[Label]] = None) -> Minor:
    """An r x r minor of (f, ..., f^(r-1)) on the closure atoms at R = x0; exact, then at the true constants."""
    ds = derivatives(f, max(r - 1, 0))
    clos = closure(f.labels())
    seen = {lb for D in ds for (lb, _m) in D}
    assert seen <= set(clos) | {lb for (lb, _m) in f.coords}, (seen, clos)
    allrows = sorted(set(clos) | seen, key=str)
    T = true_matrix(f, ds, allrows, r, x0)
    iv = _intervals(f)
    if rows is None:
        # greedy numerical pivoting at the true constants (only a choice of rows; the minor is exact)
        with mp.workdps(60):
            sub = {s: sp.Float(sp.N(t, 70), 70) for s, t in f.true.items()}
            Tn = [[mp.mpf(str(sp.N(e.subs(sub), 60))) if f.syms else mp.mpf(sp.Rational(e).p) / sp.Rational(e).q
                   for e in row] for row in T]
            left = list(range(len(allrows)))
            work = {i: Tn[i][:] for i in left}
            chosen = []
            for col in range(r):
                if not left:                       # fewer atoms than r: rank < r
                    return Minor(None, None, mp.iv.mpf(0), sp.Integer(0))
                piv = max(left, key=lambda i: abs(work[i][col]))
                if abs(work[piv][col]) < mp.mpf(10) ** -30:
                    return Minor(None, None, mp.iv.mpf(0), sp.Integer(0))
                chosen.append(piv)
                left.remove(piv)
                for i in left:
                    fct = work[i][col] / work[piv][col]
                    work[i] = [a - fct * b for a, b in zip(work[i], work[piv])]
        rows = [allrows[i] for i in chosen]
    idx = [allrows.index(lb) for lb in rows]
    if f.syms:
        dom = sp.QQ[tuple(f.syms)]
        M = DomainMatrix([[dom.from_sympy(T[i][j]) for j in range(r)] for i in idx], (r, r), dom)
        poly = sp.Poly(dom.to_sympy(M.det()), *f.syms)
        ivd = _iv_poly(poly, iv)
        value = sp.expand(poly.as_expr().xreplace(f.true))
        return Minor(rows, poly, ivd, value)
    M = DomainMatrix([[sp.QQ.from_sympy(T[i][j]) for j in range(r)] for i in idx], (r, r), sp.QQ)
    d = sp.Rational(sp.QQ.to_sympy(M.det()))
    mp.iv.dps = 70
    return Minor(rows, None, mp.iv.mpf(int(d.p)) / mp.iv.mpf(int(d.q)), d)


def generic_rank(f: Fn, ncols: int) -> int:
    """Rank over Q(R) of the split vectors of f, ..., f^(ncols-1); rows (atom, constant monomial)."""
    ds = derivatives(f, ncols - 1)
    keys = sorted(set().union(*[set(D) for D in ds]), key=str)
    M = DomainMatrix([[D.get(k, FR.zero) for D in ds] for k in keys], (len(keys), ncols), FR)
    return int(M.rank())


def point_rank(f: Fn, ncols: int, x0) -> int:
    """Rank over Q of the split vectors at ONE rational point (a lower bound on generic_rank only)."""
    ds = derivatives(f, ncols - 1)
    keys = sorted(set().union(*[set(D) for D in ds]), key=str)
    rows = [[sp.QQ.from_sympy(sp.Rational(FR.to_sympy(D.get(k, FR.zero)).subs(R, x0))) for D in ds]
            for k in keys]
    return int(DomainMatrix(rows, (len(keys), ncols), sp.QQ).rank())


def order_at_most(f: Fn, r: int) -> bool:
    """An order-r operator over Q(R) exists: closure size <= r, or generic rank of f..f^(r) is <= r."""
    if len(closure(f.labels())) <= r:
        return True
    return generic_rank(f, r + 1) <= r


def order_at_least(f: Fn, r: int, points=(R0_DEFAULT, sp.Rational(11, 7))) -> bool:
    """No operator of order < r: an r x r minor nonzero at the true constants at SOME rational point.

    A nonzero value at any one point proves the minor is not identically zero; a zero at a point proves nothing,
    so a second point is tried (the minor of a true-order-3 function can vanish exactly at 37/59)."""
    return any(lower_minor(f, r, x0).nonzero for x0 in points)


def order_is(f: Fn, r: int) -> bool:
    return order_at_least(f, r) and order_at_most(f, r)


def numeric(f: Fn, ds_index: int, x0, dps: int = 40) -> mp.mpf:
    """Value of f^(ds_index) at x0 from the atom calculus (for the cross-check against sympy.diff)."""
    ds = derivatives(f, ds_index)
    with mp.workdps(dps):
        tot = mp.mpf(0)
        x = mp.mpf(sp.Rational(x0).p) / sp.Rational(x0).q
        for ((c, kind, k), m), w in ds[ds_index].items():
            coef = mp.mpf(str(sp.N((FR.to_sympy(w).subs(R, x0)) * _mono_expr(f, m).xreplace(f.true), dps + 10)))
            if kind == "1":
                atom = mp.exp(c * x)
            elif k == 0:
                atom = mp.exp(c * x) * mp.log(x)
            else:
                atom = mp.exp(c * x) * mp.ei(k * x)
            tot += coef * atom
        return tot


__all__ = ["R", "FR", "Fn", "derivative", "derivatives", "closure", "lower_minor", "generic_rank", "point_rank",
           "order_at_most", "order_at_least", "order_is", "numeric", "coefficient_vanishes"]
