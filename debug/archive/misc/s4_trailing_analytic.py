"""Fast arbitrary-precision evaluator for log-modulated trailing-1 multiple-t-values.

Descending Hoffman convention (odd o):
    t_d(b1..bd) = sum_{o_d > .. > o_1 >= 1 odd} o_d^-b1 .. o_1^-bd

The five constants handled here are the trailing-1 forms whose direct single
sums carry  ln(o)  factors (via  pL(1,o)  inside  e2 / e3), which silently
poison  mpmath nsum method='levin'  (the project's §3 log-modulation trap).

Closed single-sum forms (verified vs the brute definition):
    e2(o) = (p1^2 - p2)/2 ,  e3(o) = (p1^3 - 3 p1 p2 + 2 p3)/6 ,  p_k = pL(k,o)
    t2(b1,1)      = sum_{o}  o^-b1            * pL(1,o)
    t3(b1,b2,1)   = sum_{o}  tau(b1,o) o^-b2  * pL(1,o)
    t3(b1,1,1)    = sum_{o}  o^-b1            * e2(o)
    t4(b1,b2,1,1) = sum_{o}  tau(b1,o) o^-b2  * e2(o)
    t4(b1,1,1,1)  = sum_{o}  o^-b1            * e3(o)

Method (analytic tail; NO mpmath.quad / sumem / Levin):
    sum_o = head[o <= 2N-1, summed term-by-term in exact mpf]
          + tail[o >= 2N+1].
For the tail every building block is replaced by its asymptotic series in 1/o:
    pL(1,o)  = (1/2) ln o  +  Lconst  +  (decaying series, NO log)
    pL(k>=2,o) = lambda(k) - 2^-k zeta(k,o/2)   (Hurwitz asymptotic, decays)
    tau(b,o)  = 2^-b zeta(b,o/2+1)              (Hurwitz asymptotic, decays)
The product summand collapses to a FINITE sum  sum_{p,j} C_{p,j} o^-p (ln o)^j
(j <= #trailing-ones, p over a finite asymptotic window).  Each piece is summed
in closed form over odd o >= 2N+1:
    S(p,j) = Lsum(p,j) - head_partial(p,j),
    Lsum(p,j) = sum_{o odd>=1} o^-p (ln o)^j = (-1)^j d^j/ds^j[(1-2^-s)zeta(s)]|_{s=p}.

The summand is expanded numerically as a *log-power series in 1/o* (a small
truncated-series algebra below), so no hand bookkeeping of the products is
needed -- the multiplications are done by the series ring, which is exact up
to the chosen truncation order.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Dict, Tuple

from mpmath import mp, mpf, zeta, digamma, log, diff, bernoulli

# --------------------------------------------------------------------------
# Closed building blocks (match debug/s4_mt_eval.py exactly)
# --------------------------------------------------------------------------


def _lam(b):
    return (1 - mpf(2) ** (-b)) * zeta(b)


def tau(b, x):
    """sum_{o > x, odd} o^-b,  b >= 2,  x a non-negative odd int (or 0)."""
    return mpf(2) ** (-b) * zeta(b, mpf(x) / 2 + 1)


def pL(b, x):
    """sum_{o < x, odd} o^-b,  x positive odd int."""
    if b == 1:
        return (digamma(mpf(x) / 2) - digamma(mpf(1) / 2)) / 2
    return _lam(b) - mpf(2) ** (-b) * zeta(b, mpf(x) / 2)


# --------------------------------------------------------------------------
# Lsum(p,j) = sum_{o odd>=1} o^-p (ln o)^j  (p >= 2 integer, j >= 0)
# --------------------------------------------------------------------------


def Lsum(p, j):
    if j == 0:
        return (1 - mpf(2) ** (-p)) * zeta(p)
    f = lambda s: (1 - mpf(2) ** (-s)) * zeta(s)
    return (-1) ** j * diff(f, mpf(p), j)


# --------------------------------------------------------------------------
# Tiny log-power series algebra:  f(o) = sum_{p, j}  c[(p,j)] * o^-p (ln o)^j
# p ranges over a *real* exponent set (here always integers >= some p0),
# j is a non-negative integer (the ln-power, bounded by #trailing ones).
# Stored as dict {(p, j): coeff}.  We carry p only up to PMAX.
# --------------------------------------------------------------------------


class LPSeries:
    __slots__ = ("c",)

    def __init__(self, c=None):
        self.c: Dict[Tuple[int, int], object] = {} if c is None else c

    def add_term(self, p, j, coeff):
        key = (p, j)
        self.c[key] = self.c.get(key, mpf(0)) + coeff

    def __add__(self, other):
        r = LPSeries(dict(self.c))
        for k, v in other.c.items():
            r.c[k] = r.c.get(k, mpf(0)) + v
        return r

    def scale(self, s):
        return LPSeries({k: v * s for k, v in self.c.items()})

    def mul(self, other, pmax):
        """Multiply two series; drop terms with exponent p > pmax."""
        r: Dict[Tuple[int, int], object] = {}
        for (p1, j1), v1 in self.c.items():
            for (p2, j2), v2 in other.c.items():
                p = p1 + p2
                if p > pmax:
                    continue
                key = (p, j1 + j2)
                r[key] = r.get(key, mpf(0)) + v1 * v2
        return LPSeries(r)


def _hurwitz_series(s, half_shift, pmax, order=24):
    """Asymptotic series of  zeta(s, o/2 + half_shift)  in powers of  o^-1.

    Using a = o/2 + half_shift, and the Hurwitz asymptotic
        zeta(s,a) ~ a^(1-s)/(s-1) + a^-s/2
                    + sum_{k>=1} B_{2k}/(2k) (s)_{2k-1} a^(-s-2k+1).
    Each a^(-s-m) is re-expanded:  a = (o/2)(1 + 2*half_shift/o), so
        a^(-w) = (2/o)^w (1 + 2c/o)^(-w),  c = half_shift,
    and (1+2c/o)^(-w) = sum_n C(-w,n) (2c)^n o^-n.  All exponents land on
    integers + (s itself is integer here, b/k), giving integer p = w + n.
    Returns an LPSeries (j = 0, pure decay, NO log).
    """
    s = mpf(s)
    c = mpf(half_shift)
    res = LPSeries()

    def binom(top, n):
        # generalized binomial C(top, n) for real top, integer n>=0
        out = mpf(1)
        for i in range(n):
            out *= (top - i)
        from math import factorial
        return out / factorial(n)

    def add_power(w, pref):
        # add  pref * a^(-w)  expanded in o^-1, a = (o/2)(1+2c/o)
        # a^(-w) = 2^w o^-w (1+2c/o)^-w
        base = pref * mpf(2) ** w
        # how many n-terms until exponent exceeds pmax
        n = 0
        while True:
            p = int(w) + n
            if p > pmax:
                break
            coeff = base * binom(-w, n) * (2 * c) ** n
            res.add_term(p, 0, coeff)
            n += 1

    # leading a^(1-s)/(s-1): exponent w = s-1
    add_power(s - 1, mpf(1) / (s - 1))
    # a^-s/2: w = s
    add_power(s, mpf(1) / 2)
    # Bernoulli tail
    for k in range(1, order + 1):
        w = s + 2 * k - 1  # a^(-s-2k+1) => a^-w
        if int(w) > pmax:
            break
        rf = mpf(1)
        for jj in range(2 * k - 1):
            rf *= (s + jj)
        pref = bernoulli(2 * k) / (2 * k) * rf
        add_power(w, pref)
    return res


def tau_series(b, pmax, order=24):
    """tau(b,o) = 2^-b zeta(b, o/2 + 1)  as LPSeries in o^-1 (no log)."""
    s = _hurwitz_series(b, 1, pmax, order)  # zeta(b, o/2+1): half_shift=1
    return s.scale(mpf(2) ** (-b))


def pLk_series(b, pmax, order=24):
    """pL(b>=2, o) = lambda(b) - 2^-b zeta(b, o/2)  as LPSeries (no log).
    The constant lambda(b) is the p=0, j=0 term."""
    res = LPSeries()
    res.add_term(0, 0, _lam(b))
    h = _hurwitz_series(b, 0, pmax, order)  # zeta(b, o/2): half_shift=0
    res = res + h.scale(-(mpf(2) ** (-b)))
    return res


def pL1_series(pmax, order=24):
    """pL(1,o) = (psi(o/2) - psi(1/2))/2  as LPSeries: (1/2)ln o + const +
    decaying series (the ONLY block carrying a ln o).

    psi(z) ~ ln z - 1/(2z) - sum_{k>=1} B_2k/(2k) z^-2k,  z = o/2.
        ln(o/2) = ln o - ln 2
    => pL1 = (1/2)[ln o - ln2 - 1/o*? ...] ... let's expand exactly:
        psi(o/2) = ln(o/2) - 1/(2*(o/2)) - sum_k B2k/(2k) (o/2)^-2k
                 = ln o - ln2 - 1/o - sum_k B2k/(2k) 2^2k o^-2k
        pL1 = (psi(o/2) - psi(1/2))/2,  psi(1/2) = -gamma - 2 ln 2
            = (1/2) ln o + (1/2)(gamma + ln2) - 1/(2o)
              - (1/2) sum_k B2k/(2k) 2^2k o^-2k
    """
    from mpmath import euler
    res = LPSeries()
    # (1/2) ln o  -> term (p=0, j=1) coeff 1/2  (o^0 (ln o)^1)
    res.add_term(0, 1, mpf(1) / 2)
    # constant (p=0,j=0): (1/2)(gamma + ln2)
    res.add_term(0, 0, (euler + log(mpf(2))) / 2)
    # -1/(2 o): p=1
    res.add_term(1, 0, -mpf(1) / 2)
    # -(1/2) sum_k B2k/(2k) 2^2k o^-2k
    for k in range(1, order + 1):
        p = 2 * k
        if p > pmax:
            break
        res.add_term(p, 0, -(bernoulli(2 * k) / (2 * k)) * mpf(2) ** (2 * k) / 2)
    return res


# --------------------------------------------------------------------------
# e2, e3 as LPSeries (Newton's identities on p_k = pL(k,o))
# --------------------------------------------------------------------------


def e2_series(pmax, order=24):
    p1 = pL1_series(pmax, order)
    p2 = pLk_series(2, pmax, order)
    return (p1.mul(p1, pmax) + p2.scale(mpf(-1))).scale(mpf(1) / 2)


def e3_series(pmax, order=24):
    p1 = pL1_series(pmax, order)
    p2 = pLk_series(2, pmax, order)
    p3 = pLk_series(3, pmax, order)
    p1p1 = p1.mul(p1, pmax)
    p1cub = p1p1.mul(p1, pmax)
    p1p2 = p1.mul(p2, pmax)
    return (p1cub + p1p2.scale(mpf(-3)) + p3.scale(mpf(2))).scale(mpf(1) / 6)


# --------------------------------------------------------------------------
# Tail evaluator: given the summand as an LPSeries (in o^-1, ln o) times an
# overall o^-q0 prefactor, sum over odd o >= o_lo = 2N+1.
# --------------------------------------------------------------------------


@lru_cache(maxsize=None)
def _head_logpow(o_lo_idx, p, j, prec, jflag):
    """sum_{o odd, 1<=o<=2*o_lo_idx-1} o^-p (ln o)^j  -- partial of Lsum.
    Cached by (#odds, p, j, prec).  o_lo_idx = number of odd terms included."""
    tot = mpf(0)
    for i in range(o_lo_idx):
        o = 2 * i + 1
        if o == 1:
            # ln 1 = 0 contributes 0 for j>=1; for j=0 it is 1
            if j == 0:
                tot += mpf(1)
            continue
        tot += mpf(o) ** (-p) * (log(o) ** j if j else mpf(1))
    return tot


def _sum_series_tail(series: LPSeries, q0, N):
    """sum_{o odd >= 2N+1}  o^-q0 * series(o).
    series has terms c[(p,j)] o^-p (ln o)^j; multiply by o^-q0 -> exponent p+q0.
    Each (p+q0, j) summed via Lsum minus the head partial over o<=2N-1.
    Returns the tail value (mpf)."""
    o_lo_idx = N  # odds 1,3,..,2N-1 are the head (o < 2N+1); count = N
    prec = mp.prec
    total = mpf(0)
    for (p, j), coeff in series.c.items():
        P = p + q0
        if coeff == 0:
            continue
        # need P >= 2 for convergence of Lsum; with N>=1 the actual exponents
        # are always >= 2 for these constants (b1>=2 or trailing structure).
        full = Lsum(P, j)
        head = _head_logpow(o_lo_idx, P, j, prec, 1 if j else 0)
        total += coeff * (full - head)
    return total


def _head_direct(term_fn, N):
    """sum over odd o = 1,3,..,2N-1 of term_fn(o) (exact mpf, no asymptotics)."""
    tot = mpf(0)
    for i in range(N):
        o = 2 * i + 1
        tot += term_fn(o)
    return tot


# --------------------------------------------------------------------------
# The five constants
# --------------------------------------------------------------------------


def _params(dps):
    """Return (N, order, pmax) tuned so tail truncation < 10^-(dps+25)."""
    # tail term magnitude ~ (1/(2N))^pmax-ish; with N modest and a deep
    # asymptotic window pmax, the residual is set by the first DROPPED term.
    # Empirically pmax ~ dps/log10(2N) + slack;  order >= pmax/2 + slack.
    import math
    N = 256
    pmax = int(dps / math.log10(2 * N)) + 14
    order = pmax // 2 + 6
    return N, order, pmax


def _t2_trailing1(b1, dps):
    N, order, pmax = _params(dps)
    # head: sum_{o<2N+1} o^-b1 pL(1,o)
    head = _head_direct(lambda o: mpf(o) ** (-b1) * pL(1, o), N)
    # tail: sum_{o>=2N+1} o^-b1 pL1_series(o)
    s = pL1_series(pmax, order)
    tail = _sum_series_tail(s, b1, N)
    return head + tail


def _t3_trailing1(b1, b2, dps):
    N, order, pmax = _params(dps)
    head = _head_direct(lambda o: tau(b1, o) * mpf(o) ** (-b2) * pL(1, o), N)
    # tail: tau(b1,o) o^-b2 pL1(o); tau is a series in o^-1 (q-shift built into
    # tau_series exponents), times o^-b2, times pL1 series.
    s_tau = tau_series(b1, pmax, order)
    s_pl1 = pL1_series(pmax, order)
    summand = s_tau.mul(s_pl1, pmax)
    tail = _sum_series_tail(summand, b2, N)
    return head + tail


def _t3_double_trailing(b1, dps):
    """t3(b1,1,1) = sum_o o^-b1 e2(o)."""
    N, order, pmax = _params(dps)
    head = _head_direct(
        lambda o: mpf(o) ** (-b1) * ((pL(1, o) ** 2 - pL(2, o)) / 2), N)
    s = e2_series(pmax, order)
    tail = _sum_series_tail(s, b1, N)
    return head + tail


def _t4_double_trailing(b1, b2, dps):
    """t4(b1,b2,1,1) = sum_o tau(b1,o) o^-b2 e2(o)."""
    N, order, pmax = _params(dps)
    head = _head_direct(
        lambda o: tau(b1, o) * mpf(o) ** (-b2)
        * ((pL(1, o) ** 2 - pL(2, o)) / 2), N)
    s_tau = tau_series(b1, pmax, order)
    s_e2 = e2_series(pmax, order)
    summand = s_tau.mul(s_e2, pmax)
    tail = _sum_series_tail(summand, b2, N)
    return head + tail


def _t4_triple_trailing(b1, dps):
    """t4(b1,1,1,1) = sum_o o^-b1 e3(o)."""
    N, order, pmax = _params(dps)
    head = _head_direct(
        lambda o: mpf(o) ** (-b1)
        * ((pL(1, o) ** 3 - 3 * pL(1, o) * pL(2, o) + 2 * pL(3, o)) / 6), N)
    s = e3_series(pmax, order)
    tail = _sum_series_tail(s, b1, N)
    return head + tail


# --------------------------------------------------------------------------
# Public dispatcher
# --------------------------------------------------------------------------


def eval_trailing(args, dps):
    """Evaluate one of the five trailing-1 forms; return mpf decimal string.

    args is the Hoffman tuple (b1, .., bd).  Recognised:
        (b1, 1)            -> t2(b1,1)
        (b1, b2, 1)        -> t3(b1,b2,1)   (b2 != 1)
        (b1, 1, 1)         -> t3(b1,1,1)
        (b1, b2, 1, 1)     -> t4(b1,b2,1,1) (b2 != 1)
        (b1, 1, 1, 1)      -> t4(b1,1,1,1)
    """
    args = tuple(int(a) for a in args)
    old = mp.dps
    mp.dps = dps + 30
    try:
        d = len(args)
        if d == 2 and args[1] == 1:
            v = _t2_trailing1(args[0], dps)
        elif d == 3 and args[1] != 1 and args[2] == 1:
            v = _t3_trailing1(args[0], args[1], dps)
        elif d == 3 and args[1] == 1 and args[2] == 1:
            v = _t3_double_trailing(args[0], dps)
        elif d == 4 and args[1] != 1 and args[2] == 1 and args[3] == 1:
            v = _t4_double_trailing(args[0], args[1], dps)
        elif d == 4 and args[1] == 1 and args[2] == 1 and args[3] == 1:
            v = _t4_triple_trailing(args[0], dps)
        else:
            raise ValueError("unsupported trailing form: %s" % (args,))
        return mp.nstr(v, dps, strip_zeros=False)
    finally:
        mp.dps = old


# --------------------------------------------------------------------------
# Gates
# --------------------------------------------------------------------------

if __name__ == "__main__":
    import json
    import time
    from mpmath import mpmathify

    cache = json.load(open("debug/data/s3_pslq_cache.json"))

    results = {"gate1": [], "gate2": None, "gate3": [], "timings": {}}
    print("=" * 78)
    print("GATE 1: BIT-EXACT vs k=3 cache @ 220 dps  (threshold < 1e-200)")
    print("=" * 78)
    DPS = 220
    g1_pass = True
    # all t2(b,1), t3(b1,b2,1), and double-trailing t3(b,1,1) present in cache
    targets = []
    for key in cache:
        if "@220" not in key:
            continue
        name = key.split("@")[0]
        if not (name.startswith("t2(") or name.startswith("t3(")):
            continue
        tup = tuple(int(x) for x in name[name.index("(") + 1:-1].split(","))
        # only trailing-1 forms this module owns
        if tup[-1] != 1:
            continue
        # we own: (b,1); (b1,b2,1) b2!=1; (b,1,1)
        if len(tup) == 2:
            targets.append((name, tup, cache[key]))
        elif len(tup) == 3 and tup[2] == 1:
            targets.append((name, tup, cache[key]))
    for name, tup, sval in sorted(targets):
        t0 = time.time()
        v = mpmathify(eval_trailing(tup, DPS))
        dt = time.time() - t0
        r = abs(v - mpmathify(sval))
        ok = float(r) < 1e-200
        g1_pass = g1_pass and ok
        results["gate1"].append({"name": name, "res": mp.nstr(r, 5),
                                 "time": round(dt, 2), "pass": ok})
        results["timings"][name] = round(dt, 2)
        print("  %-14s res %s  (%5.1fs)  %s"
              % (name, mp.nstr(r, 3), dt, "PASS" if ok else "FAIL <<<"))
    print("  GATE 1: %s" % ("PASS" if g1_pass else "FAIL"))

    print("=" * 78)
    print("GATE 2: t4(4,1,1,1) vs rigorous ground truth 0.000118762689570762...")
    print("=" * 78)
    REF = "0.000118762689570762"
    t0 = time.time()
    v411 = eval_trailing((4, 1, 1, 1), 220)
    dt = time.time() - t0
    results["timings"]["t4(4,1,1,1)"] = round(dt, 2)
    v411m = mpmathify(v411)
    # match the 15 given digits
    g2_res = abs(v411m - mpmathify(REF))
    g2_pass = float(g2_res) < 5e-19
    results["gate2"] = {"value40": mp.nstr(v411m, 40), "res_vs_ref15": mp.nstr(g2_res, 5),
                        "time": round(dt, 2), "pass": g2_pass}
    print("  t4(4,1,1,1) = %s" % mp.nstr(v411m, 40))
    print("  ref(15dig)  = %s" % REF)
    print("  |diff|      = %s   (%5.1fs)  %s"
          % (mp.nstr(g2_res, 3), dt, "PASS" if g2_pass else "FAIL <<<"))

    print("=" * 78)
    print("GATE 3: two-precision self-consistency (220 vs 150)  threshold < 1e-140")
    print("=" * 78)
    hard = [(2, 1, 1, 1), (2, 3, 1, 1), (4, 3, 1, 1), (2, 1, 1)]
    g3_pass = True
    for tup in hard:
        t0 = time.time()
        v220 = mpmathify(eval_trailing(tup, 220))
        v150 = mpmathify(eval_trailing(tup, 150))
        dt = time.time() - t0
        r = abs(v220 - v150)
        ok = float(r) < 1e-140
        g3_pass = g3_pass and ok
        results["gate3"].append({"name": "t%s" % (tup,), "res": mp.nstr(r, 5),
                                 "time": round(dt, 2), "pass": ok})
        results["timings"]["t%s" % (tup,)] = round(dt, 2)
        print("  t%-12s res %s  (%5.1fs)  %s"
              % (str(tup), mp.nstr(r, 3), dt, "PASS" if ok else "FAIL <<<"))
    print("  GATE 3: %s" % ("PASS" if g3_pass else "FAIL"))

    print("=" * 78)
    print("GATE 4: speed (each constant < ~30s @ 220 dps)")
    print("=" * 78)
    slow = {k: t for k, t in results["timings"].items() if t > 30}
    g4_pass = len(slow) == 0
    if slow:
        for k, t in slow.items():
            print("  %-14s  %5.1fs  SLOW <<<" % (k, t))
    print("  GATE 4: %s  (max %.1fs)" % ("PASS" if g4_pass else "FAIL",
          max(results["timings"].values()) if results["timings"] else 0))

    print("=" * 78)
    allpass = g1_pass and g2_pass and g3_pass and g4_pass
    print("OVERALL: %s" % ("ALL GATES PASS" if allpass else "SOME GATES FAILED"))
    results["overall"] = allpass
    results["gate_summary"] = {"gate1": g1_pass, "gate2": g2_pass,
                               "gate3": g3_pass, "gate4": g4_pass}
    results["t4_4111_40"] = mp.nstr(v411m, 40)
    results["t4_2111_40"] = mp.nstr(mpmathify(eval_trailing((2, 1, 1, 1), 220)), 40)
    json.dump(results, open("debug/data/s4_trailing_analytic_validation.json", "w"),
              indent=2)
    print("wrote debug/data/s4_trailing_analytic_validation.json")
