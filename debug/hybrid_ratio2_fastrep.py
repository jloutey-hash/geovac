"""Speed patches for the second pipeline, applied IN-PROCESS only (no file of geovac or of hybrid_ratio_* is edited).

1. Collector field.  `hybrid_ratio_rep` keeps every coefficient in sympy's pure-Python sparse field Q(u, R); profiling
   (debug/data/hybrid_ratio2_profile_rep_H6.txt) put ~90% of the collection time in its heuristic gcd.  `patch_field`
   rebinds the collector's module globals K / ONE to FK below: the same field Q(u, R), with numerator and denominator
   held as python-flint fmpq_mpoly (gcds in C).  Only the arithmetic backend changes; the walk, the keys and the
   constants are the collector's own code.  `to_sympy_K` converts a finished coefficient back to the sympy field
   (once per key) so serialisation and `hybrid_ratio_analyse` are unchanged.

2. Integral accumulation.  `T._integrate_with_e1` accumulates `total += ...` over its canonical terms, and the
   repeated flatten of a growing sympy Add dominated its time (436 of 584 s on the profiled piece).  It is linear
   over the terms of its (expanded) input, so `integrate_collect` calls it ONE TERM AT A TIME and collects each
   result with the collector, summing the collected dicts.  Same function, same calls per term.
"""
from __future__ import annotations

import flint
import sympy as sp

import hybrid_ratio_rep as HR
from hybrid_ratio_engine import u as U, R as RS

CTX = flint.fmpq_mpoly_ctx.get(("u", "R"), "lex")
_PU, _PR = CTX.gens()
_ONEP = CTX.from_dict({(0, 0): flint.fmpq(1)})
_ZEROP = CTX.from_dict({})


def _poly_from_sympy(e):
    if e == 0:
        return _ZEROP
    P = sp.Poly(e, U, RS, domain=sp.QQ)
    return CTX.from_dict({m: flint.fmpq(int(c.p), int(c.q)) for m, c in P.terms()})


class FE:
    """Element of Q(u, R): num/den, gcd-reduced, den with leading coefficient 1."""
    __slots__ = ("n", "d")

    def __init__(self, n, d=None, reduce=True):
        if d is None:
            d = _ONEP
        if reduce and not d == _ONEP:
            if n == _ZEROP:
                d = _ONEP
            else:
                g = n.gcd(d)
                if not g == _ONEP:
                    n = n / g
                    d = d / g
                lc = d.leading_coefficient()
                if lc != 1:
                    n = n / lc
                    d = d / lc
        self.n, self.d = n, d

    def __add__(self, o):
        if not isinstance(o, FE):
            o = FK(o)
        if self.d == o.d:
            return FE(self.n + o.n, self.d)
        return FE(self.n * o.d + o.n * self.d, self.d * o.d)

    __radd__ = __add__

    def __neg__(self):
        return FE(-self.n, self.d, reduce=False)

    def __sub__(self, o):
        return self + (-(o if isinstance(o, FE) else FK(o)))

    def __mul__(self, o):
        if not isinstance(o, FE):
            o = FK(o)
        g1 = self.n.gcd(o.d) if not o.d == _ONEP else _ONEP
        g2 = o.n.gcd(self.d) if not self.d == _ONEP else _ONEP
        n1 = self.n / g1 if not g1 == _ONEP else self.n
        d2 = o.d / g1 if not g1 == _ONEP else o.d
        n2 = o.n / g2 if not g2 == _ONEP else o.n
        d1 = self.d / g2 if not g2 == _ONEP else self.d
        return FE(n1 * n2, d1 * d2)

    __rmul__ = __mul__

    def __pow__(self, k: int):
        k = int(k)
        if k >= 0:
            return FE(self.n ** k, self.d ** k, reduce=False)
        assert not self.n == _ZEROP
        return FE(self.d ** (-k), self.n ** (-k))

    def __eq__(self, o):
        if isinstance(o, FE):
            return self.n * o.d == o.n * self.d
        if o == 0:
            return self.n == _ZEROP
        return self == FK(o)

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash((str(self.n), str(self.d)))

    def to_sympy_K(self):
        K = HR_K_ORIG
        num = K.ring.from_dict({tuple(int(x) for x in m): sp.QQ(int(c.p), int(c.q)) for m, c in self.n.to_dict().items()})
        den = K.ring.from_dict({tuple(int(x) for x in m): sp.QQ(int(c.p), int(c.q)) for m, c in self.d.to_dict().items()})
        return K(num) / K(den)


class _FK:
    """Stand-in for the collector's field object K: K.one, K.zero, K(x), K.from_expr(e)."""

    @property
    def one(self):
        return FE(_ONEP, reduce=False)

    @property
    def zero(self):
        return FE(_ZEROP, reduce=False)

    def __call__(self, x):
        if isinstance(x, FE):
            return x
        x = sp.sympify(x)
        assert x.is_Rational, x
        return FE(CTX.from_dict({(0, 0): flint.fmpq(int(x.p), int(x.q))}) if x != 0 else _ZEROP, reduce=False)

    def from_expr(self, e):
        e = sp.sympify(e)
        if e.is_Rational:
            return self(e)
        n, d = sp.fraction(sp.together(e))
        return FE(_poly_from_sympy(sp.expand(n)), _poly_from_sympy(sp.expand(d)))


FK = _FK()
HR_K_ORIG = HR.K


def patch_field():
    HR.K = FK
    HR.ONE = FK.one
    HR._CACHE.clear()


def integrate_collect(T, expr, v, X, log=None):
    """sum over the Add-args t of expand(expr) of rep(T._integrate_with_e1(t, v, X)) -- the collected form of
    T._integrate_with_e1(expr, v, X), term by term."""
    import time
    t0 = time.time()
    args = sp.Add.make_args(sp.expand(expr))
    tot: dict = {}
    for i, t in enumerate(args):
        r = HR.rep(T._integrate_with_e1(t, v, X))
        tot = HR.add(tot, r)
        if log and (i + 1) % 500 == 0:
            log(f"    {i + 1}/{len(args)} terms {time.time() - t0:.1f}s, {len(tot)} keys")
    return tot
