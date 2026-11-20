"""S^(4) corrected multiple-t evaluator: single-nsum + closed factors + Abel.

Replaces the retired nested-cascade evaluator (memo S1.6b: nested nsum is
intractable past one acceleration layer).  Architecture mirrors the proven
k=3 t3val: ONE accelerated sum, with the largest variable a CLOSED Hurwitz
tail and the smallest a CLOSED partial; trailing-1s peeled by Abel so the
single remaining accelerated sum is a clean decaying tail (no log
modulation -> Levin safe, the §3 trap avoided).

Descending Hoffman convention:
    t_d(b1..bd) = sum_{o_d > .. > o_1 >= 1 odd} o_d^-b1 .. o_1^-bd
Closed blocks (o odd):
    tau_b(x)  = sum_{o>x odd}  o^-b = 2^-b zeta(b, x/2+1)            (b>=2)
    pH_b(x)   = sum_{o>=x odd}  o^-b = 2^-b zeta(b, x/2)             (b>=2)
    pL_b(x)   = sum_{o<x odd}  o^-b = lam(b) - pH_b(x)              (b>=2)
    pL_1(x)   = sum_{o<x odd}  o^-1 = (digamma(x/2)-digamma(1/2))/2

Strategy by trailing structure (k = #trailing 1s at the small end):
  k=0: middle-variable single nsum, closed tau & pL on the ends.
       t2(b1,b2)        = sum_{o} o^-b1 pL_b2(o)
       t3(b1,b2,b3)     = sum_{o} o^-b2 tau_b1(o) pL_b3(o)
       t4(b1,b2,b3,b4)  = sum_{o3} tau_b1(o3) o3^-b2 W(o3),
           W(o3) = sum_{o2<o3} o2^-b3 pL_b4(o2)  [cumprefix of CLOSED w]
  k>=1: Abel — make the smallest variable's o1^-1 explicit, inner = the
       (d-1)-fold TAIL above o1 (a clean decaying tail, computed by the
       k=0 machinery restricted to vars > o1).  Recurse on remaining
       trailing 1s.

Every accelerated sum is a tail with pure algebraic decay; the cumprefix
W is a cumulative sum of CLOSED terms (NOT a nested nsum).
"""
from __future__ import annotations

from functools import lru_cache

from mpmath import mp, mpf, nsum, inf, zeta, digamma

GUARD = 15


class _Pref:
    """Precision-aware cumulative prefix of a CLOSED term g(o) over odds.
    cum(idx) = sum of g over the first idx odds (o = 1,3,..,2*idx-1).
    Rebuilds from scratch when mp.prec changes — Levin re-requests terms at
    elevated internal precision and a cache frozen at the calling precision
    silently serves stale values (the k=3 v2->v3 precision-contract lesson,
    memo S1.6b)."""

    __slots__ = ("g", "prec", "arr")

    def __init__(self, g):
        self.g = g
        self.prec = None
        self.arr = [mpf(0)]

    def cum(self, idx):
        if self.prec != mp.prec:
            self.prec = mp.prec
            self.arr = [mpf(0)]
        while len(self.arr) <= idx:
            o = 2 * len(self.arr) - 1
            self.arr.append(self.arr[-1] + self.g(o))
        return self.arr[idx]


def _levin(term, start=0):
    """nsum(term,[s,inf],levin) with leading zeros skipped (Levin rejects a
    zero leading weight).  term(k) for integer k; advances start past any
    initial exact zeros (capped)."""
    s = start
    while s < start + 8 and term(s) == 0:
        s += 1
    return nsum(term, [s, inf], method='levin')


def _lam(b):
    return (1 - mpf(2) ** (-b)) * zeta(b)


def tau(b, x):
    """sum_{o > x, odd} o^-b,  b >= 2,  x a positive odd int (or 0)."""
    return mpf(2) ** (-b) * zeta(b, mpf(x) / 2 + 1)


def pL(b, x):
    """sum_{o < x, odd} o^-b,  x positive odd int."""
    if b == 1:
        return (digamma(mpf(x) / 2) - digamma(mpf(1) / 2)) / 2
    return _lam(b) - mpf(2) ** (-b) * zeta(b, mpf(x) / 2)


# --------------------------------------------------------------------------
# k=0 (no trailing 1) single-nsum forms
# --------------------------------------------------------------------------

def _t2_k0(b1, b2):
    return _levin(lambda k: mpf(2 * int(k) + 1) ** (-b1)
                  * pL(b2, 2 * int(k) + 1), 1)


def _t3_k0(b1, b2, b3):
    return _levin(lambda k: mpf(2 * int(k) + 1) ** (-b2)
                  * tau(b1, 2 * int(k) + 1) * pL(b3, 2 * int(k) + 1), 1)


class _Const:
    """One-shot constant: evaluated once (at the elevated working precision set
    by eval_t) and cached forever.  Levin's internal precision bump is small
    (<~15 digits, k=3 evidence), so a constant computed at dps+trailing-guard
    stays accurate beyond the target — recomputing it on every Levin precision
    step (the first design) was correct but catastrophically slow for nested
    Abel (double/triple-trailing t4).  Correctness rests on the PREFIX being
    precision-aware (_Pref); the constant only needs generous fixed headroom."""
    __slots__ = ("thunk", "v")

    def __init__(self, thunk):
        self.thunk = thunk
        self.v = None

    def __call__(self):
        if self.v is None:
            self.v = self.thunk()
        return self.v


def _t4_k0(b1, b2, b3, b4):
    """outer over o3; inner W(o3) = sum_{o2<o3} o2^-b3 pL_b4(o2) (cumprefix)."""
    W = _Pref(lambda o: mpf(o) ** (-b3) * pL(b4, o))
    return _levin(lambda k: tau(b1, 2 * int(k) + 1)
                  * mpf(2 * int(k) + 1) ** (-b2) * W.cum((2 * int(k) + 1 - 1) // 2),
                  1)


# --------------------------------------------------------------------------
# Abel layer: peel trailing 1s.  t_d(..,1) = sum_{o1 odd} o1^-1 TAIL(o1),
# TAIL(o1) = (full) - (cumprefix over o<=o1).  full and prefix both
# precision-aware (Levin re-requests at elevated prec; memo S1.6b).
# --------------------------------------------------------------------------

def _t3_trailing1(b1, b2):
    """t3(b1,b2,1) = sum_{o1} o1^-1 B(o1), B(o1)=sum_{o2>o1} o2^-b2 tau_b1(o2)."""
    full = _Const(lambda: _t2_k0(b1, b2))
    bpref = _Pref(lambda o: mpf(o) ** (-b2) * tau(b1, o))
    return _levin(lambda k: mpf(2 * int(k) - 1) ** (-1)
                  * (full() - bpref.cum((2 * int(k) - 1 + 1) // 2)), 1)


def _t3_trailing11(b1):
    """t3(b1,1,1) = sum_{o1} o1^-1 F(o1), F(o1)=sum_{o2>o1} o2^-1 tau_b1(o2);
       double Abel.  Ffull = sum_{o2>=1} o2^-1 tau_b1(o2) = t2(b1,1)."""
    Ffull = _Const(lambda: _t2_trailing1(b1))
    fpref = _Pref(lambda o: mpf(o) ** (-1) * tau(b1, o))
    return _levin(lambda k: mpf(2 * int(k) - 1) ** (-1)
                  * (Ffull() - fpref.cum((2 * int(k) - 1 + 1) // 2)), 1)


def _t4_trailing1(b1, b2, b3):
    """t4(b1,b2,b3,1) = sum_{o1} o1^-1 C(o1),
       C(o1)=sum_{o2>o1} o2^-b3 D(o2), D(o2)=sum_{o3>o2} o3^-b2 tau_b1(o3)."""
    Dtot = _Const(lambda: _t2_k0(b1, b2))          # sum_{o3>=1} o3^-b2 tau_b1(o3)
    d3pref = _Pref(lambda o: mpf(o) ** (-b2) * tau(b1, o))

    def D(o2):
        return Dtot() - d3pref.cum((o2 + 1) // 2)

    Cfull = _Const(lambda: _t3_k0(b1, b2, b3))     # = sum_{o2>=1} o2^-b3 D(o2)
    cpref = _Pref(lambda o: mpf(o) ** (-b3) * D(o))
    return _levin(lambda k: mpf(2 * int(k) - 1) ** (-1)
                  * (Cfull() - cpref.cum((2 * int(k) - 1 + 1) // 2)), 1)


def _t4_trailing2(b1, b2):
    """t4(b1,b2,1,1) = sum_{o1} o1^-1 sum_{o2>o1} o2^-1 E(o2),
       E(o2)=sum_{o3>o2} o3^-b2 tau_b1(o3).  Double Abel."""
    Etot = _Const(lambda: _t2_k0(b1, b2))
    e3pref = _Pref(lambda o: mpf(o) ** (-b2) * tau(b1, o))

    def E(o2):
        return Etot() - e3pref.cum((o2 + 1) // 2)

    Ffull = _Const(lambda: _t3_trailing1(b1, b2))  # = sum_{o2>=1} o2^-1 E(o2)
    fpref = _Pref(lambda o: mpf(o) ** (-1) * E(o))
    return _levin(lambda k: mpf(2 * int(k) - 1) ** (-1)
                  * (Ffull() - fpref.cum((2 * int(k) - 1 + 1) // 2)), 1)


def _t4_trailing3(b1):
    """t4(b1,1,1,1) = triple Abel over o1<o2<o3 with o3^-1 tau_b1(o3) core."""
    Gfull = _Const(lambda: _levin(
        lambda k: mpf(2 * int(k) + 1) ** (-1) * tau(b1, 2 * int(k) + 1), 0))
    g3pref = _Pref(lambda o: mpf(o) ** (-1) * tau(b1, o))

    def G(o2):
        return Gfull() - g3pref.cum((o2 + 1) // 2)

    Hfull = _Const(lambda: _levin(
        lambda k: mpf(2 * int(k) + 1) ** (-1) * G(2 * int(k) + 1), 0))
    hpref = _Pref(lambda o: mpf(o) ** (-1) * G(o))
    return _levin(lambda k: mpf(2 * int(k) - 1) ** (-1)
                  * (Hfull() - hpref.cum((2 * int(k) - 1 + 1) // 2)), 1)


# --------------------------------------------------------------------------

def _t2_trailing1(b1):
    """t2(b1,1) = sum_{o1} o1^-1 tau_b1(o1) -- Abel form, log-free
    (the direct sum_o o^-b1 pL_1(o) is log-modulated -> Levin trap)."""
    return _levin(lambda k: mpf(2 * int(k) - 1) ** (-1)
                  * tau(b1, 2 * int(k) - 1), 1)


def eval_t(args, dps):
    args = tuple(int(a) for a in args)
    old = mp.dps
    d = len(args)
    ntr = 0
    for b in reversed(args):
        if b == 1:
            ntr += 1
        else:
            break
    # trailing-Abel tail subtractions need cancellation headroom; nested Abel
    # (double/triple-trailing t4) stacks one-shot constants, so scale guard
    # with the trailing depth
    mp.dps = dps + (40 + 15 * ntr if ntr else GUARD)
    try:
        if d == 2:
            if ntr == 1:
                v = _t2_trailing1(args[0])
            else:
                v = _t2_k0(*args)
        elif d == 3:
            if ntr == 0:
                v = _t3_k0(*args)
            elif ntr == 1:
                v = _t3_trailing1(args[0], args[1])
            else:  # t3(b1,1,1) double-trailing
                v = _t3_trailing11(args[0])
        elif d == 4:
            if ntr == 0:
                v = _t4_k0(*args)
            elif ntr == 1:
                v = _t4_trailing1(args[0], args[1], args[2])
            elif ntr == 2:
                v = _t4_trailing2(args[0], args[1])
            else:
                v = _t4_trailing3(args[0])
        else:
            raise ValueError("depth %d unsupported" % d)
        return mp.nstr(v, dps, strip_zeros=False)
    finally:
        mp.dps = old


if __name__ == "__main__":
    import json
    import time
    from mpmath import mpmathify
    k3 = json.load(open("debug/data/s3_pslq_cache.json"))
    DPS = 50
    print("=== core validation vs k=3 cache @ %d dps ===" % DPS)
    worst = 0.0
    n = 0
    for key, sval in sorted(k3.items()):
        if "@220" not in key:
            continue
        name = key.split("@")[0]
        tup = tuple(int(x) for x in name[name.index("(") + 1:-1].split(","))
        t0 = time.time()
        try:
            v = mpmathify(eval_t(tup, DPS))
        except NotImplementedError:
            continue
        r = abs(v - mpmathify(sval))
        worst = max(worst, float(r))
        n += 1
        flag = "" if float(r) < 1e-45 else "  <-- FAIL"
        print("  t%-14s res %.2e  (%.1fs)%s" % (str(tup), float(r),
                                                time.time() - t0, flag))
    print("validated %d symbols; worst residual %.2e" % (n, worst))
