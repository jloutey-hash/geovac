"""AHA Track-3 helper: exact term census + trans-series decomposition of the
HYBRID two-centre ERI class closed form (geovac.two_center_eri.hybrid_closed_form).

Every closed-form term has the shape

    coeff * R^k * exp(-c R) * X ,      X in { 1 , E_1(lam R) , log(q) }

with coeff algebraic over Q(Z_A, Z_B) and c, lam, q rational.  This module parses
that shape exactly (sympy, no floats) and regroups it into the large-R trans-series

    F(R) = sum_A e^{-A R} phi_A(1/R) ,   A = c + lam  (E_1 terms) or A = c (elementary)

which is what the resurgence analysis in debug/aha_t3_object2_hybrid.py consumes.
"""
from __future__ import annotations

from collections import defaultdict

import sympy as sp

from geovac.two_center_eri import R_s


class Term:
    __slots__ = ("coeff", "k", "c", "lam", "logarg")

    def __init__(self, coeff, k, c, lam, logarg):
        self.coeff = sp.sympify(coeff)     # algebraic prefactor (already exact)
        self.k = sp.Integer(k)             # power of R
        self.c = sp.sympify(c)             # exp(-c R)
        self.lam = None if lam is None else sp.sympify(lam)     # E_1(lam R)
        self.logarg = None if logarg is None else sp.sympify(logarg)

    @property
    def kind(self):
        if self.lam is not None:
            return "E1"
        if self.logarg is not None:
            return "log"
        return "elem"

    @property
    def action(self):
        """Exponential action A of the trans-series sector this term belongs to."""
        return self.c + (self.lam if self.lam is not None else 0)

    def __repr__(self):
        x = (f"E1({self.lam}R)" if self.lam is not None
             else (f"log({self.logarg})" if self.logarg is not None else "1"))
        return f"Term({self.coeff}, R^{self.k}, exp(-{self.c}R), {x})"


def parse_closed_form(expr, R=R_s):
    """Exact parse of the expanded closed form into Term objects.

    Raises if any term does not match the expected {exp, E_1, log} shape -- that
    is the point: an unmatched term would mean the seed set is bigger than claimed.
    """
    expr = sp.expand(expr)
    terms = expr.args if expr.is_Add else (expr,)
    out = []
    for t in terms:
        rest = t
        lam = None
        logarg = None
        c = sp.Integer(0)

        e1s = [f for f in t.atoms(sp.Function) if isinstance(f, sp.expint)]
        if len(e1s) > 1:
            raise ValueError(f"more than one E_1 in a single term: {t}")
        if e1s:
            f = e1s[0]
            if f.args[0] != 1:
                raise ValueError(f"expint of order != 1: {f}")
            lam = sp.simplify(f.args[1] / R)
            if lam.has(R):
                raise ValueError(f"E_1 argument not linear in R: {f}")
            rest = rest / f

        lgs = [f for f in t.atoms(sp.log)]
        if len(lgs) > 1:
            raise ValueError(f"more than one log in a single term: {t}")
        if lgs:
            f = lgs[0]
            if f.args[0].has(R):
                raise ValueError(f"R-DEPENDENT LOG (breaks the census): {f}")
            logarg = f.args[0]
            rest = rest / f

        exs = [f for f in t.atoms(sp.exp)]
        if len(exs) > 1:
            raise ValueError(f"more than one exp in a single term: {t}")
        if exs:
            f = exs[0]
            c = sp.simplify(-f.args[0] / R)
            if c.has(R):
                raise ValueError(f"exp argument not linear in R: {f}")
            rest = rest / f

        rest = sp.powsimp(sp.simplify(rest))
        coeff, k = rest.as_coeff_exponent(R)
        if coeff.has(R):
            raise ValueError(f"unparsed R dependence in {t} -> {rest}")
        out.append(Term(coeff, k, c, lam, logarg))
    return out


def group_prefactors(terms):
    """dict (c, lam) -> {k: coeff} for the E_1 terms (the divergence carriers)."""
    g = defaultdict(dict)
    for t in terms:
        if t.kind == "E1":
            key = (t.c, t.lam)
            g[key][t.k] = g[key].get(t.k, sp.Integer(0)) + t.coeff
    return {k: v for k, v in g.items()}


def group_logs(terms):
    """dict (c, logarg) -> {k: coeff} for the log terms."""
    g = defaultdict(dict)
    for t in terms:
        if t.kind == "log":
            key = (t.c, t.logarg)
            g[key][t.k] = g[key].get(t.k, sp.Integer(0)) + t.coeff
    return {k: v for k, v in g.items()}


def group_elementary(terms):
    g = defaultdict(dict)
    for t in terms:
        if t.kind == "elem":
            g[t.c][t.k] = g[t.c].get(t.k, sp.Integer(0)) + t.coeff
    return {k: v for k, v in g.items()}


def sector_actions(terms):
    return sorted({t.action for t in terms}, key=lambda z: float(z))


def field_of_coefficients(terms, gens=(sp.sqrt(3),)):
    """Return (all_in_field, offenders). Tests coeff in Q + sum Q*gen."""
    offenders = []
    for t in terms:
        c = sp.sympify(t.coeff)
        if c.is_Rational:
            continue
        ok = False
        for g in gens:
            q = sp.radsimp(sp.cancel(c / g))
            if q.is_Rational:
                ok = True
                break
        if not ok:
            offenders.append(t)
    return (len(offenders) == 0), offenders


def phi_coefficients(terms, A, N):
    """Coefficients a_1..a_N of  phi_A(1/R) = sum_{m} a_m R^{-m}  (exact, sympy).

    Elementary and log terms of action A contribute their (finite) Laurent
    coefficients; each E_1 term contributes the full Gevrey-1 tail
        R^k e^{-cR} E_1(lam R) = e^{-(c+lam)R} R^k sum_n (-1)^n n! /(lam R)^{n+1}.
    """
    a = [sp.Integer(0)] * (N + 1)          # index m = power R^{-m}
    for t in terms:
        if t.action != A:
            continue
        base = t.coeff * (sp.log(t.logarg) if t.logarg is not None else 1)
        if t.kind in ("elem", "log"):
            m = -int(t.k)
            if 0 <= m <= N:
                a[m] += base
        else:
            lam = t.lam
            for n in range(0, N + 1):
                m = n + 1 - int(t.k)
                if 0 <= m <= N:
                    a[m] += base * (-1) ** n * sp.factorial(n) / lam ** (n + 1)
    return a
