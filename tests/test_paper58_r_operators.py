"""Paper 58 sec:census -- minimal linear ODEs in R for Slater-1s two-centre integrals.

Guard for the R-operators found by the derivation pass (CHANGELOG, chemistry
period inventory rows 2, 3a, 3b, 3c). Written as a SEPARATE pass from the
derivation (CLAUDE.md Sec. 9); the operators below are frozen literals, and every
closed form is rebuilt live from geovac (`qfd_core`, `two_center_eri`).

What each test rejects:

* ``test_*_annihilates``: a stored operator that does not kill the live closed
  form. Annihilation is checked by sympy's OWN derivatives of exp / E1 / Ei / log
  (not by this file's atom calculus), with gamma_E and the logarithms of
  numbers kept as independent symbols, so every coefficient must cancel
  identically -- an exact identity, no tolerance.
* ``test_*_minimal`` / ``no_lower_order``: a claimed order r for which an order r-1 operator exists.
  f, ..., f^(r-1) are vectors over Q(R) indexed by (atom, constant monomial);
  full rank r at a rational point R0 proves full generic rank (rank can only
  drop under specialisation), hence no order-(r-1) operator over Q(R).
  Assumptions (standard, stated in the paper): the atoms e^{cR}, e^{cR}E1(kR),
  e^{cR}Ei(kR), e^{cR}ln R at distinct (c, k) are linearly independent over
  Q(R); gamma_E and the ln of the integers involved are algebraically
  independent over Q. (Since 2026-11-18 the second assumption is no longer
  needed for the conclusion: the true-constant minor below certifies rank r
  by itself.  The first is Lemma GI of the minimal-order derivation, a hand
  proof, cold-reviewed.) The atom-calculus derivative is itself cross-checked
  against sympy's derivative (``test_atom_derivative_matches_sympy``).
* Since 2026-11-18 ``no_lower_order`` ALSO requires an r x r minor that is nonzero at the TRUE constants
  (interval arithmetic), and ``decompose`` writes Ei(-kR) as -E1(kR) so one function has one label; the
  one-electron upper bound is a generic rank over Q(R) (``order_at_most``), not a rank at a point.  These close
  three holes found by the cold review of the two-centre minimal-order memo; section (6) holds a control for each:

    W1  constants as independent symbols: f = (ln 6 - ln 2 - ln 3) R e^{-R} + e^{-2R} (true order 1; also with
        e^{-R}E1(R)): BEFORE no_lower_order(f, 2) was True; AFTER False (the minor's interval contains 0).
    W2  redundant spelling: f = e^{-R}E1(2R) + e^{-R}Ei(-2R) + e^{-3R} = e^{-3R}: BEFORE decompose kept two atoms
        and no_lower_order(f, 2) was True; AFTER one atom e^{-3R}, False.
    W3  rank at a point as an upper bound: f = e^{-R} + e^{-2R} + (2R/3 + 103/177) e^{-3R} (true order 3) has rank
        2 at 37/59: BEFORE the one-electron test's idiom accepted "order 2 exists"; AFTER order_at_most(f, 2) is
        False (generic rank 3).

  The split/spec ranks are kept; the true-constant minor is an extra condition, so no assertion was weakened.
  For the minimal-order results themselves (orders proved with true constants, Lemma GI) see
  tests/test_paper58_minimal_order.py.
* R = 0: exponents from the indicial polynomial; log-freeness of the LOCAL
  SOLUTION SPACE from an exact triangular Frobenius count (dimension of Laurent
  solutions == order  <=>  no log solution, when the exponents are integers);
  log-freeness of f ITSELF from the exact ln R coefficient of its expansion.
* R = infinity: Newton-polygon slope exactly 1 (irregular rank 1), and the
  characteristic roots equal the atom decay rates.
* Finite singular points != 0: each irreducible factor of the leading
  coefficient is simple and APPARENT, certified by an exact power-series count
  over Q(rho) = Q[x]/(factor): r independent analytic solutions at rho.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from typing import Dict, List, Tuple

from itertools import combinations

import mpmath as mp
import pytest
import sympy as sp
from sympy.polys.matrices import DomainMatrix

from geovac import qfd_core as Q
from geovac.two_center_eri import R_s, aabb_closed_form, hybrid_closed_form

R = R_s
FR = sp.QQ.frac_field(R)
RG = FR.gens[0]
S1 = (1, 0, 0)
P0 = (2, 1, 0)
TH = sp.Symbol("theta")

# ---------------------------------------------------------------- frozen operators
# Coefficient lists [p_0, p_1, ..., p_r] of  sum_j p_j(R) D^j.

OP_S_11 = ["R*(R + 1)", "R**2 + 3*R + 3"]

OP_AABB_11 = [
    "4*R*(4*R**2 + 6*R + 3)",
    "2*(8*R**4 + 20*R**3 + 30*R**2 + 30*R + 15)",
    "R*(8*R**3 + 24*R**2 + 30*R + 15)",
]
OP_AABB_12 = [
    "32*(3*R + 1)**2",
    "2*(144*R**3 + 312*R**2 + 250*R + 81)",
    "6*(36*R**3 + 57*R**2 + 39*R + 11)",
    "R*(36*R**2 + 51*R + 22)",
]

OP_HYB_33 = [
    "-162*R**3*(28697814*R**8 + 78121827*R**7 + 114791256*R**6 + 111202389*R**5 + 66460743*R**4 + 26651268*R**3 + 3800196*R**2 - 2003040*R - 357344)",
    "-2324522934*R**11 - 5553027009*R**10 - 7246198035*R**9 - 8734764276*R**8 - 6184201770*R**7 - 3043779120*R**6 - 929335032*R**5 + 39914208*R**4 + 63006336*R**3 + 19438848*R**2 + 3534336*R + 392704",
    "R*(258280326*R**10 + 1133563653*R**9 + 1922753538*R**8 + 1515787830*R**7 + 1013517036*R**6 + 375516648*R**5 - 9599472*R**4 - 23245056*R**3 - 15020928*R**2 - 3534336*R - 392704)",
    "R**2*(258280326*R**9 + 789189885*R**8 + 1063413441*R**7 + 1008793116*R**6 + 659577330*R**5 + 261495216*R**4 + 49319280*R**3 - 2829168*R**2 - 3534336*R - 392704)",
    "R**3*(28697814*R**8 + 87687765*R**7 + 108413964*R**6 + 87996132*R**5 + 43954326*R**4 + 9519768*R**3 + 191160*R**2 - 441792*R - 49088)",
]
OP_HYB_34 = [
    "-24*R**2*(112116629520000*R**11 + 347493062448000*R**10 + 892418798048400*R**9 + 1425184703453280*R**8 + 1072372649913908*R**7 - 148231711332984*R**6 - 1024386508785430*R**5 - 1049042840949548*R**4 - 690428081460057*R**3 - 303954657842113*R**2 - 79907662567664*R - 10053748776552)",
    "-8*(252262416420000*R**13 + 721129549518000*R**12 + 1872517974881400*R**11 + 2707359071898180*R**10 + 1919363923532888*R**9 - 143465336858336*R**8 - 1609977584758080*R**7 - 1813592913088928*R**6 - 1274755201578035*R**5 - 574915098828592*R**4 - 164400351588224*R**3 - 29202578539812*R**2 - 1850948522824*R - 30257195952)",
    "-4*R*(79415945910000*R**12 + 194754130704000*R**11 + 514752065875950*R**10 + 569524263471590*R**9 + 336158897434654*R**8 + 22737051778904*R**7 - 185268404751441*R**6 - 246027053581895*R**5 - 147390655215910*R**4 - 44998222705486*R**3 - 9275500974936*R**2 - 1771967230432*R + 60514391904)",
    "4*R*(22189749592500*R**12 + 65854964715750*R**11 + 167408593413600*R**10 + 258537201986770*R**9 + 179267284906367*R**8 - 25883621692172*R**7 - 168983617351983*R**6 - 182307214458518*R**5 - 126015073061574*R**4 - 54901126581672*R**3 - 14152194779960*R**2 - 2354109845376*R - 190395279168)",
    "R**2*(30364920495000*R**11 + 81266007280500*R**10 + 205053019051650*R**9 + 263483574052580*R**8 + 149107348374733*R**7 - 37173575124335*R**6 - 145140973601324*R**5 - 140822005753452*R**4 - 78816820371556*R**3 - 25571381265796*R**2 - 4235182193984*R - 349058011808)",
    "R**3*(2335763115000*R**10 + 6071557243500*R**9 + 14839223428050*R**8 + 16843531059110*R**7 + 6287171087561*R**6 - 5277648504780*R**5 - 8728852505588*R**4 - 6132310631872*R**3 - 2289102716660*R**2 - 382265908912*R - 31732546528)",
]

OP_EXCH_1 = [
    "24*R*(32*R**6 + 192*R**5 + 576*R**4 + 1380*R**3 + 2070*R**2 + 1890*R + 945)",
    "-16*(16*R**8 + 80*R**7 + 72*R**6 - 246*R**5 - 1035*R**4 - 2430*R**3 - 3915*R**2 - 4050*R - 2025)",
    "-2*R*(64*R**7 + 768*R**6 + 2976*R**5 + 5400*R**4 + 6480*R**3 + 5940*R**2 + 4050*R + 2025)",
    "2*R**2*(32*R**6 + 64*R**5 - 480*R**4 - 1980*R**3 - 3690*R**2 - 4050*R - 2025)",
    "R**3*(32*R**5 + 240*R**4 + 720*R**3 + 1260*R**2 + 1350*R + 675)",
]


def _ops(strs: List[str]) -> List[sp.Expr]:
    return [sp.expand(sp.sympify(s, locals={"R": R})) for s in strs]


# ---------------------------------------------------------------- live closed forms


def _A(z) -> tuple:
    return ("A", Fraction(z), 1)


def _B(z) -> tuple:
    return ("B", Fraction(z), 1)


@lru_cache(maxsize=None)
def overlap(za, zb) -> sp.Expr:
    return Q.overlap(_A(za), _B(zb), R)


@lru_cache(maxsize=None)
def aabb(za, zb) -> sp.Expr:
    return aabb_closed_form(Fraction(za), S1, S1, Fraction(zb), S1, S1)


@lru_cache(maxsize=None)
def hyb_p(za, zb) -> sp.Expr:
    return hybrid_closed_form(Fraction(za), P0, P0, S1, Fraction(zb), S1)


@lru_cache(maxsize=None)
def exch(z) -> sp.Expr:
    return Q.exchange_closed_form(Fraction(z), S1, S1, Fraction(z), S1, S1, R, tau_max=6)


# ---------------------------------------------------------------- atom calculus


def _const_symbols(expr: sp.Expr) -> Dict[sp.Expr, sp.Symbol]:
    reps: Dict[sp.Expr, sp.Symbol] = {}
    k = 0
    for a in sorted(expr.atoms(sp.log, sp.Pow), key=str):
        if a.free_symbols:
            continue
        if isinstance(a, sp.log) or (isinstance(a, sp.Pow) and not a.exp.is_integer):
            reps[a] = sp.Symbol(f"k{k}")
            k += 1
    if expr.has(sp.EulerGamma):
        reps[sp.EulerGamma] = sp.Symbol("gE")
    return reps


def decompose(expr: sp.Expr):
    """f -> {(c, kind, k): coeff}; f = sum coeff * e^{cR} * {1, ln R, E1(kR), Ei(kR)}."""
    expr = sp.expand(sp.expand_log(sp.expand(sp.expand_func(expr)), force=True))
    reps = _const_symbols(expr)
    expr = sp.expand(expr.xreplace(reps))
    atoms: dict = {}
    for term in sp.Add.make_args(expr):
        c, kind, k, coeff = sp.Integer(0), "1", sp.Integer(0), sp.Integer(1)
        for fct in sp.Mul.make_args(term):
            if isinstance(fct, sp.exp):
                assert sp.diff(fct.args[0], R, 2) == 0 and fct.args[0].subs(R, 0) == 0
                c += sp.diff(fct.args[0], R)
            elif isinstance(fct, sp.Pow) and isinstance(fct.base, sp.exp):
                c += sp.diff(fct.base.args[0] * fct.exp, R)
            elif isinstance(fct, sp.expint):
                assert fct.args[0] == 1 and kind == "1", term
                kind, k = "E1", sp.diff(fct.args[1], R)
                assert k > 0, fct
            elif isinstance(fct, sp.Ei):
                assert kind == "1", term
                kind, k = "Ei", sp.diff(fct.args[0], R)
                if k < 0:                       # W2: Ei(-kR) = -E1(kR); one label per function
                    kind, k = "E1", -k
                    coeff *= -1
            elif fct == sp.log(R):
                assert kind == "1"
                kind = "log"
            else:
                assert not fct.has(sp.exp, sp.expint, sp.Ei, sp.log), fct
                coeff *= fct
        key = (sp.nsimplify(c), kind, sp.nsimplify(k))
        atoms[key] = atoms.get(key, 0) + coeff
    atoms = {a: sp.cancel(sp.together(v)) for a, v in atoms.items()}
    return {a: v for a, v in atoms.items() if v != 0}, reps


def to_field(atoms, reps, mode: str = "split"):
    """{(atom, constant monomial): Q(R)}; mode 'spec' fixes constants at rationals."""
    csyms = sorted(set(reps.values()), key=str)
    cvals = {s: sp.Rational(7 + 3 * i, 13 + 5 * i) for i, s in enumerate(csyms)}
    out: dict = {}
    for a, v in atoms.items():
        if mode == "spec" or not csyms:
            parts = {(): sp.together(v.subs(cvals)) if csyms else sp.together(v)}
        else:
            num, den = sp.fraction(sp.together(v))
            assert not den.has(*csyms)
            P = sp.Poly(sp.expand(num), *csyms)
            parts = {m: cf / den for m, cf in P.terms()}
        for m, w in parts.items():
            el = FR.from_sympy(sp.together(w))
            if el != 0:
                out[(a, m)] = out.get((a, m), FR.zero) + el
    return out


def d_atoms(atoms):
    out: dict = {}

    def add(key, v):
        out[key] = out.get(key, FR.zero) + v

    for ((c, kind, k), m), a in atoms.items():
        add(((c, kind, k), m), a.diff(RG) + FR.from_sympy(c) * a)
        if kind == "log":
            add(((c, "1", 0), m), a / RG)
        elif kind == "E1":
            add(((c - k, "1", 0), m), -a / RG)
        elif kind == "Ei":
            add(((c + k, "1", 0), m), a / RG)
    return {key: v for key, v in out.items() if v != 0}


def deriv_vectors(expr: sp.Expr, upto: int, mode: str = "split"):
    atoms, reps = decompose(expr)
    ds = [to_field(atoms, reps, mode)]
    for _ in range(upto):
        ds.append(d_atoms(ds[-1]))
    return ds


def rank_at(ds, x0) -> int:
    keys = sorted(set().union(*[set(D) for D in ds]), key=str)
    rows = [[sp.QQ.from_sympy(sp.Rational(FR.to_sympy(D.get(kk, FR.zero)).subs(R, x0)))
             for D in ds] for kk in keys]
    return DomainMatrix(rows, (len(keys), len(ds)), sp.QQ).rank()


def _true_constant_interval(expr: sp.Expr, dps: int = 70):
    """mpmath interval around the true value of a constant (gamma, log of a number, surd)."""
    v = sp.N(expr, dps + 15)
    with mp.workdps(dps + 15):
        x = mp.mpf(str(v))
        w = abs(x) * mp.mpf(10) ** (-(dps - 5))
        lo, hi = x - w, x + w
    mp.iv.dps = dps
    return mp.iv.mpf([lo, hi])


def true_minor_nonzero(expr: sp.Expr, r: int, points=(sp.Rational(37, 59), sp.Rational(11, 7))) -> bool:
    """W1: some r x r minor of the coordinates of f, ..., f^(r-1) is nonzero AT THE TRUE CONSTANTS.

    The minor is computed exactly as a polynomial in the constant symbols and then evaluated at the true values
    of gamma, the logs and the surds in interval arithmetic; it counts only if the interval excludes 0.  Unlike the
    split/spec modes this assumes nothing about the arithmetic of gamma, ln 2, ln 3, ... (a coefficient such as
    ln 6 - ln 2 - ln 3, independent as symbols and zero in truth, is refused).  A nonzero value at any one point
    proves the minor is not identically zero; all row choices are tried.
    """
    atoms, reps = decompose(expr)
    csyms = sorted(set(reps.values()), key=str)
    true = {v: k for k, v in reps.items()}
    ds = [to_field(atoms, reps, "split")]
    for _ in range(r - 1):
        ds.append(d_atoms(ds[-1]))
    keys = sorted({a for D in ds for (a, _m) in D}, key=str)
    if len(keys) < r:
        return False
    ivals = {s: _true_constant_interval(true[s]) for s in csyms}
    for x0 in points:
        T = []
        for a in keys:
            row = []
            for j in range(r):
                tot = sp.Integer(0)
                for (aa, m), w in ds[j].items():
                    if aa == a:
                        mono = sp.Integer(1)
                        for s_, e in zip(csyms, m):
                            mono *= s_ ** e
                        tot += sp.Rational(FR.to_sympy(w).subs(R, x0)) * mono
                row.append(sp.expand(tot))
            T.append(row)
        for rows in combinations(range(len(keys)), r):
            if csyms:
                dom = sp.QQ[tuple(csyms)]
                det = sp.Poly(dom.to_sympy(DomainMatrix([[dom.from_sympy(T[i][j]) for j in range(r)] for i in rows],
                                                        (r, r), dom).det()), *csyms)
                mp.iv.dps = 70
                tot = mp.iv.mpf(0)
                for mono, cf in det.terms():
                    cf = sp.Rational(cf)
                    term = mp.iv.mpf(int(cf.p)) / mp.iv.mpf(int(cf.q))
                    for g, e in zip(det.gens, mono):
                        if e:
                            term = term * ivals[g] ** int(e)
                    tot = tot + term
                if tot.a > 0 or tot.b < 0:
                    return True
            else:
                M_ = DomainMatrix([[sp.QQ.from_sympy(T[i][j]) for j in range(r)] for i in rows], (r, r), sp.QQ)
                if M_.det() != 0:
                    return True
    return False


def no_lower_order(expr: sp.Expr, r: int) -> bool:
    """True iff f, ..., f^(r-1) are Q(R)-independent: both constant modes AND the true constants (W1)."""
    x0 = sp.Rational(37, 59)
    return (all(rank_at(deriv_vectors(expr, r - 1, mode), x0) == r for mode in ("split", "spec"))
            and true_minor_nonzero(expr, r))


def generic_rank(ds) -> int:
    """W3: rank over Q(R) (not at a point) of the split vectors in ds."""
    keys = sorted(set().union(*[set(D) for D in ds]), key=str)
    return int(DomainMatrix([[D.get(kk, FR.zero) for D in ds] for kk in keys], (len(keys), len(ds)), FR).rank())


def order_at_most(expr: sp.Expr, r: int) -> bool:
    """W3: an order-r operator over Q(R) exists iff the split vectors of f..f^(r) have GENERIC rank <= r (a Q(R)
    relation among split vectors is a relation among the true functions).  A rank at one point only bounds the
    generic rank from below, so it cannot prove this."""
    return generic_rank(deriv_vectors(expr, r)) <= r


def annihilates(polys: List[sp.Expr], expr: sp.Expr) -> bool:
    """Route 1: sympy's own derivatives of the closed form, then exact atom collection."""
    res = sum(p * sp.diff(expr, R, j) for j, p in enumerate(polys))
    atoms, _ = decompose(res)
    return len(atoms) == 0


def annihilates_atoms(polys: List[sp.Expr], expr: sp.Expr) -> bool:
    """Route 2 (fast): the atom calculus d_atoms (itself checked against sympy below)."""
    ds = deriv_vectors(expr, len(polys) - 1)
    Pf = [FR.from_sympy(p) for p in polys]
    keys = set().union(*[set(D) for D in ds])
    return all(sum((Pf[j] * D.get(kk, FR.zero) for j, D in enumerate(ds)), FR.zero) == 0
               for kk in keys)


def atom_rates(expr: sp.Expr) -> set:
    atoms, _ = decompose(expr)
    rs = set()
    for (c, kind, k) in atoms:
        rs.add(c)
        if kind == "E1":
            rs.add(c - k)
        if kind == "Ei":
            rs.add(c + k)
    return rs


def lnR_coefficient(expr: sp.Expr) -> sp.Expr:
    """Exact coefficient B(R) of ln R in f's expansion at R = 0 (constants re-inserted).

    E1(kR) = Ein(kR) - gamma - ln k - ln R ;  Ei(kR) = gamma + ln(kR) - Ein(-kR).
    """
    atoms, reps = decompose(expr)
    back = {v: k for k, v in reps.items()}
    B = sp.Integer(0)
    for (c, kind, k), a in atoms.items():
        sgn = {"log": 1, "E1": -1, "Ei": 1}.get(kind, 0)
        B += sgn * a.xreplace(back) * sp.exp(c * R)
    return sp.simplify(B)


# ---------------------------------------------------------------- singular points


def _ord0(p: sp.Expr) -> int:
    return min(m[0] for m in sp.Poly(p, R).monoms())


def indicial_at_zero(polys) -> Tuple[bool, sp.Poly, int]:
    r = len(polys) - 1
    vals = [(_ord0(p) - k) if p != 0 else None for k, p in enumerate(polys)]
    v = min(x for x in vals if x is not None)
    ind = sum(sp.Poly(p, R).coeff_monomial(R ** _ord0(p)) * sp.ff(TH, k)
              for k, p in enumerate(polys) if vals[k] == v)
    return vals[r] == v, sp.Poly(sp.expand(ind), TH), v


def exponents_at_zero(polys) -> Dict[int, int]:
    fuchs, ind, _ = indicial_at_zero(polys)
    assert fuchs, "R = 0 not regular singular"
    return {int(k): int(m) for k, m in sp.roots(ind).items()}


def laurent_dim_at_zero(polys, normalise: Tuple[int, sp.Rational] | None = None):
    """Dimension of the space of Laurent solutions sum_{m>=emin} a_m R^m.

    Exact: the recursion is triangular (equation for R^{m+v} involves a_k, k <= m,
    with diagonal ind(m)), and ind(m) != 0 for m > emax, so truncating at emax
    loses nothing. Returns (dim, solution-or-None).
    """
    exps = exponents_at_zero(polys)
    emin, emax = min(exps), max(exps)
    _, _, v = indicial_at_zero(polys)
    N = emax + (12 if normalise else 0)
    a = {m: sp.Symbol(f"a{m - emin}") for m in range(emin, N + 1)}
    y = sum(a[m] * R ** m for m in a)
    res = sp.expand(sum(p * sp.diff(y, R, j) for j, p in enumerate(polys)) * R ** (-v))
    eqs = [res.coeff(R, m) for m in range(emin, N + 1)]
    unknowns = list(a.values())
    M, _ = sp.linear_eq_to_matrix(eqs, unknowns)
    dim = len(unknowns) - M.rank()
    if normalise is None:
        return dim, None
    m0, val = normalise
    sol = sp.linsolve(eqs + [a[m] for m in a if m < 0] + [a[m0] - val], unknowns)
    sol = list(sol)
    assert len(sol) == 1 and not any(s.free_symbols for s in sol[0])
    return dim, {m: sol[0][i] for i, m in enumerate(a)}


def infinity_rank_and_rates(polys):
    r = len(polys) - 1
    w = [sp.Poly(p, R).degree() - j for j, p in enumerate(polys) if p != 0]
    js = [j for j, p in enumerate(polys) if p != 0]
    top = w[js.index(r)]
    rank = max(sp.Rational(wj - top, r - j) for wj, j in zip(w, js) if j < r)
    D = max(sp.Poly(p, R).degree() for p in polys if p != 0)
    lam = sp.Symbol("lam")
    chi = sum(sp.Poly(p, R).LC() * lam ** k for k, p in enumerate(polys)
              if p != 0 and sp.Poly(p, R).degree() == D)
    return rank, set(sp.roots(sp.Poly(chi, lam)).keys()), sp.Poly(chi, lam).degree()


def _mod(p: sp.Poly, fac: sp.Poly) -> sp.Poly:
    return p.rem(fac)


def _rank_mod(rows: List[List[sp.Poly]], fac: sp.Poly) -> int:
    """Rank over the field Q[x]/(fac) (fac irreducible) by Gaussian elimination."""
    M = [[_mod(e, fac) for e in row] for row in rows]
    rank, ncol = 0, len(M[0]) if M else 0
    for col in range(ncol):
        piv = next((i for i in range(rank, len(M)) if not M[i][col].is_zero), None)
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        inv = sp.Poly(sp.invert(M[rank][col].as_expr(), fac.as_expr(), fac.gen),
                      fac.gen, domain="QQ")
        for i in range(len(M)):
            if i != rank and not M[i][col].is_zero:
                f = _mod(M[i][col] * inv, fac)
                M[i] = [_mod(M[i][j] - f * M[rank][j], fac) for j in range(ncol)]
        rank += 1
    return rank


def finite_points_apparent(polys) -> List[Tuple[str, int]]:
    """Every root rho != 0 of p_r: simple, and r analytic solutions at rho (exact).

    Returns [(factor, e_last)] where {0..r-2, e_last} are the exponents.
    """
    r = len(polys) - 1
    x = sp.Symbol("x")
    out = []
    for fac_e, mult in sp.factor_list(polys[-1])[1]:
        fac = sp.Poly(fac_e.subs(R, x), x, domain="QQ")
        if fac.degree() == 0 or fac.as_expr() == x:
            continue
        assert mult == 1, f"multiple root factor {fac_e}"
        # Taylor coefficients of p_j(x + t) in t, as polys in x mod fac
        tay = []
        for p in polys:
            P = sp.Poly(p.subs(R, x), x, domain="QQ")
            cs, Pi, fact = [], P, 1
            for i in range(r + 2):
                cs.append(_mod(Pi * sp.Rational(1, fact), fac))
                Pi = Pi.diff(x)
                fact *= i + 1
            tay.append(cs)
        assert tay[r][0].is_zero and not tay[r][1].is_zero
        # exponent e from the indicial factor c (m - r + 1) + p_{r-1}(rho)
        e_poly = _mod(tay[r - 1][0] * sp.Poly(sp.invert(tay[r][1].as_expr(), fac.as_expr(), x),
                                              x, domain="QQ"), fac)
        assert e_poly.degree() <= 0, "exponent not rational"
        e = r - 1 - e_poly.as_expr()
        assert e.is_Integer and e >= r - 1, f"non-apparent exponent {e}"
        e = int(e)
        rows = []
        for s in range(1 - r, e + 2 - r):              # coefficient of t^s
            row = []
            for k in range(e + 1):
                acc = sp.Poly(0, x, domain="QQ")
                for j in range(r + 1):
                    i = s - k + j
                    if 0 <= i < len(tay[j]) and k >= j:
                        acc += tay[j][i] * sp.ff(k, j)
                row.append(acc)
            rows.append(row)
        dim = (e + 1) - _rank_mod(rows, fac)
        assert dim == r, f"only {dim} analytic solutions at a root of {fac_e}: not apparent"
        out.append((str(fac_e), e))
    return out


# ================================================================ (1) one-electron


def test_overlap_operator_annihilates_live_closed_form():
    S = overlap(1, 1)
    assert sp.simplify(S - (1 + R + R ** 2 / 3) * sp.exp(-R)) == 0
    assert annihilates(_ops(OP_S_11), S)


def test_overlap_operator_minimal_order_one():
    assert no_lower_order(overlap(1, 1), 1)


ONE_E = [
    ("S", lambda a, b: Q.overlap(_A(a), _B(b), R)),
    ("T", lambda a, b: Q.kinetic(_A(a), _B(b), R)),
    ("V_A ab", lambda a, b: Q._inv_r(_A(a), _B(b), "A", R)),
    ("V_B aa", lambda a, b: Q._inv_r(_A(a), _A(a), "B", R)),
]


@pytest.mark.parametrize("za,zb", [(1, 1), (1, 2), (Fraction(3, 2), 1)])
@pytest.mark.parametrize("name,fn", ONE_E, ids=[n for n, _ in ONE_E])
def test_one_electron_order_equals_number_of_rates(name, fn, za, zb):
    """Order (exactly: f..f^(r) dependent, f..f^(r-1) not) == # distinct rates.

    Note the rate 0 (the 1/R tail of <a|1/r_B|a>) counts: that element has order
    2 even at equal exponents.
    """
    f = fn(za, zb)
    r = len(atom_rates(f))
    assert no_lower_order(f, r)
    assert order_at_most(f, r)                                         # order-r exists (W3: generic rank)
    assert rank_at(deriv_vectors(f, r), sp.Rational(37, 59)) == r     # (implied by the two lines above)
    if name != "V_B aa":
        assert r == (1 if za == zb else 2)


# ================================================================ (2) (AA|BB)


@pytest.mark.parametrize("ops,za,zb,r,exps", [
    (OP_AABB_11, 1, 1, 2, {-1: 1, 0: 1}),
    (OP_AABB_12, 1, 2, 3, {-1: 1, 0: 1, 1: 1}),
])
def test_aabb_operator(ops, za, zb, r, exps):
    f, P = aabb(za, zb), _ops(ops)
    assert len(P) - 1 == r
    assert annihilates(P, f)
    assert no_lower_order(f, r)
    assert exponents_at_zero(P) == exps
    assert laurent_dim_at_zero(P)[0] == r            # no log solution at R = 0
    assert lnR_coefficient(f) == 0


def test_aabb_odd_taylor_family_from_operator_alone():
    """Analytic solution at R = 0 with J(0) = 5/8: coefficients from the ODE only."""
    P = _ops(OP_AABB_11)
    _, a = laurent_dim_at_zero(P, normalise=(0, sp.Rational(5, 8)))
    assert a[1] == a[3] == a[5] == 0
    assert a[7] == sp.Rational(2, 315)
    assert a[9] == sp.Rational(16, 14175)
    # and the closed form agrees (independent route)
    ser = sp.series(aabb(1, 1), R, 0, 11).removeO()
    assert all(ser.coeff(R, m) == a[m] for m in range(11))


# ================================================================ (3) hybrid


@pytest.mark.parametrize("ops,za,zb,r", [
    (OP_HYB_33, 3, 3, 4),          # rate coincidence a_c = a_d
    (OP_HYB_34, 3, 4, 5),          # a_d > a_c (Ei branch)
])
def test_hybrid_p_operator(ops, za, zb, r):
    f, P = hyb_p(za, zb), _ops(ops)
    assert len(P) - 1 == r
    assert annihilates(P, f)
    assert no_lower_order(f, r)
    assert lnR_coefficient(f) == 0                   # gamma + ln R pair cancels in Ein


def test_hybrid_coincidence_carries_the_log_pair():
    """At a_c = a_d the closed form carries a bare ln R atom that the E1 cancels."""
    atoms, _ = decompose(hyb_p(3, 3))
    assert any(kind == "log" for (_c, kind, _k) in atoms)
    assert any(kind == "E1" and c == -3 for (c, kind, _k) in atoms)


# ================================================================ (4) exchange


def test_exchange_operator_order_four():
    f, P = exch(1), _ops(OP_EXCH_1)
    assert len(P) - 1 == 4
    assert annihilates_atoms(P, f)
    assert no_lower_order(f, 4)


@pytest.mark.slow
def test_exchange_operator_annihilates_by_sympy_derivatives():
    """Same identity by the independent route (~40 s)."""
    assert annihilates(_ops(OP_EXCH_1), exch(1))


def test_exchange_log_solution_at_zero():
    P = _ops(OP_EXCH_1)
    assert exponents_at_zero(P) == {-1: 1, 0: 1, 4: 1, 9: 1}
    assert laurent_dim_at_zero(P)[0] == 3            # one solution carries ln R


@pytest.mark.parametrize("z", [1, 2])
def test_exchange_R9_lnR_leading(z):
    B = lnR_coefficient(exch(z))
    ser = sp.series(B, R, 0, 11).removeO()
    for m in range(9):
        assert sp.simplify(ser.coeff(R, m)) == 0
    assert sp.nsimplify(ser.coeff(R, 9)) == sp.Rational(8, 3375) * sp.Integer(z) ** 10


# ================================================================ (5) global singularities

ALL = [
    ("S11", OP_S_11, lambda: overlap(1, 1)),
    ("AABB11", OP_AABB_11, lambda: aabb(1, 1)),
    ("AABB12", OP_AABB_12, lambda: aabb(1, 2)),
    ("HYB33", OP_HYB_33, lambda: hyb_p(3, 3)),
    ("HYB34", OP_HYB_34, lambda: hyb_p(3, 4)),
    ("EXCH1", OP_EXCH_1, lambda: exch(1)),
]


@pytest.mark.parametrize("name,ops,fn", ALL, ids=[a[0] for a in ALL])
def test_infinity_irregular_rank_one_rates_match(name, ops, fn):
    P = _ops(ops)
    rank, rates, deg = infinity_rank_and_rates(P)
    assert rank == 1
    assert deg == len(P) - 1                         # all r exponential behaviours
    assert rates == atom_rates(fn())


@pytest.mark.parametrize("name,ops,fn", ALL, ids=[a[0] for a in ALL])
def test_finite_singular_points_apparent(name, ops, fn):
    P = _ops(ops)
    pts = finite_points_apparent(P)
    assert pts, "expected at least one finite non-zero singular point"


def test_atom_derivative_matches_sympy():
    """The atom calculus used by the minimality rank tests IS d/dR."""
    for f in (aabb(1, 2), hyb_p(3, 3), hyb_p(3, 4), exch(1)):
        a1 = d_atoms(deriv_vectors(f, 0)[0])
        atoms, reps = decompose(sp.diff(f, R))
        a2 = to_field(atoms, reps, "split")
        assert set(a1) == set(a2) and all(a1[k] == a2[k] for k in a1)


# ================================================================ (6) controls for the minimality machinery (W1-W3)
# Added 2026-11-18 after the cold review of the two-centre minimal-order memo found three ways the routines above
# could accept a wrong order.  Each control is a function whose true order is known by construction.


@pytest.mark.parametrize("variant", ["R_exp", "E1"])
def test_no_lower_order_refuses_dependent_constants(variant):
    """W1 rejects: gamma / logs of numbers read as independent symbols.  The coefficient ln 6 - ln 2 - ln 3 is
    ZERO, so f has order 1; the split mode alone reported order >= 2 (review G1, memo K2)."""
    c = sp.log(6) - sp.log(2) - sp.log(3)
    f = (c * R * sp.exp(-R) if variant == "R_exp" else c * sp.exp(-R) * sp.expint(1, R)) + sp.exp(-2 * R)
    assert rank_at(deriv_vectors(f, 1, "split"), sp.Rational(37, 59)) == 2      # the trap: split mode says 2
    assert not true_minor_nonzero(f, 2)
    assert not no_lower_order(f, 2)
    assert no_lower_order(f, 1)


def test_no_lower_order_refuses_a_redundant_spelling():
    """W2 rejects: E1(kR) next to Ei(-kR) kept as two atoms.  f = e^{-R}E1(2R) + e^{-R}Ei(-2R) + e^{-3R} IS
    e^{-3R}, order 1 (review G3)."""
    f = sp.exp(-R) * sp.expint(1, 2 * R) + sp.exp(-R) * sp.Ei(-2 * R) + sp.exp(-3 * R)
    atoms, _ = decompose(f)
    assert set(atoms) == {(-3, "1", 0)}
    assert not no_lower_order(f, 2)
    assert no_lower_order(f, 1)


def test_order_upper_bound_is_not_read_at_a_point():
    """W3 rejects: 'an order-r operator exists' read from a rank at one point.  The plant
    e^{-R} + e^{-2R} + (2R/3 + 103/177) e^{-3R} has true order 3, but f, f', f'' have rank 2 at R = 37/59, so the
    old idiom rank_at(..., r) == r accepted order 2 (review G4)."""
    f = sp.exp(-R) + sp.exp(-2 * R) + (2 * R / 3 + sp.Rational(103, 177)) * sp.exp(-3 * R)
    assert rank_at(deriv_vectors(f, 2), sp.Rational(37, 59)) == 2            # the trap: point rank says 2
    assert no_lower_order(f, 2)
    assert not order_at_most(f, 2)
    assert order_at_most(f, 3) and true_minor_nonzero(f, 3)    # order exactly 3 (minor nonzero at 11/7)
    # (no_lower_order(f, 3) itself is False here: its split/spec ranks are read at 37/59 only, where this
    #  plant is singular -- a false NEGATIVE, never a false positive; the panel cases are regular at 37/59.)


# ================================================================ FIRE TESTS (2026-11-18, section (6) and the W fixes)
# In-memory plants (a pytest plugin in the session scratchpad patched this module after collection; no committed
# file edited) unless marked [text], which went through debug/qa/fire_test.py (in place, restored from a byte copy).
# "panel" = test_one_electron_order_equals_number_of_rates with every parametrisation's closed form replaced by the
# plant and the rate count forced to 2 (the wrong claim "order 2").
#
#   plant                                                        -k selector                          result
#   W1 fix reverted (true_minor_nonzero -> True)                 refuses_dependent_constants          FIRED
#   W1 fix reverted                                              one_electron/aabb/hybrid_p/exchange  NOT FIRED (expected:
#                                                                  _operator_order (22 panel tests)     panel constants are
#                                                                                                       Q-independent)
#   W1 plant (ln6 - ln2 - ln3) R e^{-R} + e^{-2R} in panel        one_electron_order                   FIRED (no_lower_order)
#   same, W1 fix reverted                                        one_electron_order                   NOT FIRED (the hole)
#   W2 plant e^{-R}E1(2R) + e^{-R}Ei(-2R) + e^{-3R} in panel      one_electron_order                   FIRED (no_lower_order)
#   [text] W2 fix reverted (`if k < 0:` -> `if False:`)          redundant_spelling                   FIRED
#   W3 plant e^{-R} + e^{-2R} + (2R/3 + 103/177)e^{-3R} in panel  one_electron_order                   FIRED (order_at_most)
#   same, W3 fix reverted (order_at_most -> point rank)          one_electron_order                   NOT FIRED (the hole)
#   W3 fix reverted                                              not_read_at_a_point                  FIRED
# Clean runs after the plants: fast 41 passed, 1 skipped (18 s); --slow 42 passed (49 s).
