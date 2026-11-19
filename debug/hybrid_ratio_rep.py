"""Fast exact collector for the parametric hybrid closed forms (replaces sympy expand + cancel, which time out).

rep(expr) walks the UNEXPANDED expression tree and returns {key: coefficient}, key = (c, kind, kappa, consts):

    atom      e^{cR} * y,   y = 1 (kind '1') | ln R ('log', kappa 0) | Ei(kappa R) ('Ei';  E1(kR) = -Ei(-kR))
    c, kappa  sympy expressions in u (expanded, linear)
    consts    frozenset of (name, exponent): Euler gamma 'GAM', pi 'PI', logs of irreducible u-polynomials or primes
              'LOG:<factor>', fractional powers of primes 'P:<p>' (exponent kept in [0, 1))
    coeff     element of Q(u, R) (sympy sparse rational-function field: fast exact arithmetic)

Products of two y-atoms are refused (assert): the closed forms are linear in E1 / Ei / ln R.
"""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp
from sympy.polys.fields import field

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "debug"))
from hybrid_ratio_engine import u, R  # noqa: E402

K, KU, KR = field([u, R], sp.QQ)
ONE = K.one
ZERO_KEY = (sp.Integer(0), "1", sp.Integer(0), frozenset())


def to_K(e):
    return K.from_expr(e) if not isinstance(e, type(ONE)) else e


def _norm_consts(d: dict):
    """Normalise a const-exponent dict; returns (frozenset, rational multiplier)."""
    mult = sp.Integer(1)
    out = {}
    for k, v in d.items():
        v = sp.Rational(v)
        if k.startswith("P:"):
            p = sp.Integer(k[2:])
            fl = sp.floor(v)
            mult *= p ** fl
            v = v - fl
        if v != 0:
            out[k] = v
    return frozenset(out.items()), mult


def _mul_keys(k1, k2):
    c1, kind1, kap1, s1 = k1
    c2, kind2, kap2, s2 = k2
    assert kind1 == "1" or kind2 == "1", "product of two y-atoms"
    kind, kap = (kind1, kap1) if kind1 != "1" else (kind2, kap2)
    d = dict(s1)
    for n, e in s2:
        d[n] = d.get(n, 0) + e
    s, mult = _norm_consts(d)
    return (sp.expand(c1 + c2), kind, kap, s), mult


def add(a, b):
    o = dict(a)
    for k, v in b.items():
        o[k] = o.get(k, K.zero) + v
    return {k: v for k, v in o.items() if v != 0}


def mul(a, b):
    o: dict = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            k, m = _mul_keys(k1, k2)
            o[k] = o.get(k, K.zero) + v1 * v2 * K(m)
    return {k: v for k, v in o.items() if v != 0}


def scalar(e):
    return {ZERO_KEY: to_K(e)}


def const_key(name, exp=1):
    s, m = _norm_consts({name: exp})
    return {(sp.Integer(0), "1", sp.Integer(0), s): K(m)}


LOGNAMES: dict = {}


def log_const(arg):
    """log of a positive constant (u-expression or number) -> sum of LOG atoms of irreducible factors."""
    arg = sp.factor(sp.together(arg))
    c, facs = sp.factor_list(arg)
    c = sp.Rational(c)
    assert c > 0, arg
    out: dict = {}
    pieces = [(sp.Integer(p), e) for p, e in sp.factorint(c.p).items()]
    pieces += [(sp.Integer(p), -e) for p, e in sp.factorint(c.q).items()]
    pieces += [(sp.expand(f), e) for f, e in facs]
    for f, e in pieces:
        name = f"LOG:{f}"
        LOGNAMES[name] = f
        out = add(out, {k: v * K(e) for k, v in const_key(name).items()})
    return out


def is_scalar(r):
    return all(k[1] == "1" and k[0] == 0 and not k[3] for k in r)


def lin(arg):
    arg = sp.expand(arg)
    c = sp.expand(sp.diff(arg, R))
    rest = sp.expand(arg - c * R)
    return c, rest


_CACHE: dict = {}


def rep(e):
    if e in _CACHE:
        return _CACHE[e]
    r = _rep(e)
    _CACHE[e] = r
    return r


def _rep(e):
    if e.is_Number:
        return scalar(e) if e != 0 else {}
    if e == u or e == R:
        return scalar(e)
    if e is sp.pi:
        return const_key("PI")
    if e is sp.EulerGamma:
        return const_key("GAM")
    if isinstance(e, sp.Add):
        out: dict = {}
        for a in e.args:
            out = add(out, rep(a))
        return out
    if isinstance(e, sp.Mul):
        out = scalar(1)
        for a in e.args:
            out = mul(out, rep(a))
        return out
    if isinstance(e, sp.exp):
        c, rest = lin(e.args[0])
        assert rest == 0, f"exp with constant part {e}"
        return {(c, "1", sp.Integer(0), frozenset()): K.one}
    if e.func is sp.expint:
        assert e.args[0] == 1
        k, rest = lin(e.args[1])
        assert rest == 0 and k != 0, e
        return {(sp.Integer(0), "Ei", sp.expand(-k), frozenset()): -K.one}
    if isinstance(e, sp.Ei):
        k, rest = lin(e.args[0])
        assert rest == 0 and k != 0, e
        return {(sp.Integer(0), "Ei", k, frozenset()): K.one}
    if isinstance(e, sp.log):
        a = sp.factor(e.args[0])
        if not a.has(R):
            return log_const(a)
        pw = sp.degree(sp.numer(sp.together(a)), R) - sp.degree(sp.denom(sp.together(a)), R)
        rest = sp.cancel(a / R ** pw)
        assert not rest.has(R), e
        out = {(sp.Integer(0), "log", sp.Integer(0), frozenset()): K(pw)}
        return add(out, log_const(rest))
    if isinstance(e, sp.Pow):
        b, x = e.args
        if b.is_Number and not x.is_Integer:
            assert b.is_Rational and b > 0 and x.is_Rational, e
            b = sp.Rational(b)
            out = scalar(1)
            for p, m in sp.factorint(b.p).items():
                out = mul(out, const_key(f"P:{p}", m * x))
            for p, m in sp.factorint(b.q).items():
                out = mul(out, const_key(f"P:{p}", -m * x))
            return out
        if b is sp.pi:
            return const_key("PI", x)
        assert x.is_Integer, e
        rb = rep(b)
        if x >= 0:
            out = scalar(1)
            for _ in range(int(x)):
                out = mul(out, rb)
            return out
        if isinstance(b, sp.exp):
            c, rest = lin(b.args[0] * x)
            assert rest == 0
            return {(c, "1", sp.Integer(0), frozenset()): K.one}
        assert len(rb) == 1 and ZERO_KEY in rb, f"cannot invert {b}"
        return {ZERO_KEY: rb[ZERO_KEY] ** int(x)}
    if e.free_symbols <= {u, R} and not e.has(sp.exp, sp.log, sp.Ei, sp.expint):
        return scalar(e)
    raise ValueError(f"unhandled {e!r}")


def rep_parts(parts, log=print):
    import time
    t0 = time.time()
    tot: dict = {}
    for i, p in enumerate(parts):
        tot = add(tot, rep(p))
        log(f"  rep part {i}: {time.time() - t0:.1f}s, {len(tot)} keys")
    return tot
