"""Manual Euler-Maclaurin evaluator for log-modulated trailing multiple-t.

The Levin/Abel route under-converges on the slowest log-modulated trailing
constants (b1=2 family); mpmath.sumem silently mis-handles the boundary.
This module sums sum_{o odd>=1} g(o) for a CLOSED smooth g by a manual,
controlled Euler-Maclaurin tail (explicit head + integral + Bernoulli
corrections), which is exact to working precision for smooth g with
log-polynomial modulation.

Trailing constants reduce to single sums with CLOSED summands:
    t2(b1,1)       = sum_o  o^-b1 pL_1(o)
    t3(b1,b2,1)    = sum_o  tau_b1(o) o^-b2 pL_1(o)
    t3(b1,1,1)     = sum_o  o^-b1 e2(o)
    t4(b1,b2,1,1)  = sum_o  tau_b1(o) o^-b2 e2(o)
    t4(b1,1,1,1)   = sum_o  o^-b1 e3(o)
with e2(o)=(p1^2-p2)/2, e3(o)=(p1^3-3 p1 p2+2 p3)/6, p_k=pL_k(o) (closed).
(Single-trailing t4(b1,b2,b3,1) has a tail factor, not closed -> stays Abel.)
"""
from __future__ import annotations

from mpmath import mp, mpf, quad, diff, bernoulli, inf, fac

from s4_mt_eval import tau, pL


def sum_em_odd(g, M=1500, J=6):
    """sum_{o=1,3,5,...} g(o), g closed smooth.  o=2m+1; manual E-M from m=M.
    Large explicit head (cheap) + few low-order Bernoulli corrections keeps
    `diff` order low (2J-1) and stable; E-M tail error ~ B_{2J+2}/(2J+2)!
    h^{(2J+1)}(M) is negligible for M this large (the convergence-radius
    factor (2*pi*M)^{-2J} dominates)."""
    head = mpf(0)
    for m in range(M):
        head += g(2 * m + 1)

    def h(x):                      # continuous reparametrization
        return g(2 * x + 1)

    integral = quad(h, [M, inf])
    s = integral + h(M) / 2
    for j in range(1, J + 1):
        s -= bernoulli(2 * j) / fac(2 * j) * diff(h, M, 2 * j - 1)
    return head + s


def _e2(o):
    p1, p2 = pL(1, o), pL(2, o)
    return (p1 ** 2 - p2) / 2


def _e3(o):
    p1, p2, p3 = pL(1, o), pL(2, o), pL(3, o)
    return (p1 ** 3 - 3 * p1 * p2 + 2 * p3) / 6


def t2_t1(b1, M=64, J=18):
    return sum_em_odd(lambda o: mpf(o) ** (-b1) * pL(1, o), M, J)


def t3_t1(b1, b2, M=64, J=18):
    return sum_em_odd(lambda o: tau(b1, o) * mpf(o) ** (-b2) * pL(1, o), M, J)


def t3_t11(b1, M=64, J=18):
    return sum_em_odd(lambda o: mpf(o) ** (-b1) * _e2(o), M, J)


def t4_t11(b1, b2, M=64, J=18):
    return sum_em_odd(lambda o: tau(b1, o) * mpf(o) ** (-b2) * _e2(o), M, J)


def t4_t111(b1, M=64, J=18):
    return sum_em_odd(lambda o: mpf(o) ** (-b1) * _e3(o), M, J)


if __name__ == "__main__":
    import json
    from mpmath import mpmathify
    mp.dps = 60
    k3 = json.load(open("debug/data/s3_pslq_cache.json"))

    def ref(name):
        return mpmathify(k3[name + "@220"])

    print("=== E-M vs k=3 cache (bit-exact target) ===")
    # t2(b,1)
    for b1 in (2, 4, 6):
        v = t2_t1(b1)
        print("  t2(%d,1)  res %.2e" % (b1, float(abs(v - ref("t2(%d,1)" % b1)))))
    # t3(b1,b2,1) trailing  (cache has e.g. t3(2,3,1),(3,2,1),(3,4,1),(4,3,1))
    for (b1, b2) in [(2, 3), (3, 2), (3, 4), (4, 3)]:
        nm = "t3(%d,%d,1)" % (b1, b2)
        if nm + "@220" in k3:
            v = t3_t1(b1, b2)
            print("  %s  res %.2e" % (nm, float(abs(v - ref(nm)))))
    print("=== E-M vs rigorous e3 ground-truth ===")
    print("  t4(4,1,1,1) E-M = %s  vs 0.000118762689570762: %.2e"
          % (mp.nstr(t4_t111(4), 18),
             float(abs(t4_t111(4) - mpf("0.000118762689570762")))))
    # two-precision self-consistency on the hard one
    a = t4_t111(2)
    mp.dps = 90
    b = t4_t111(2)
    print("  t4(2,1,1,1) E-M two-prec |v90-v60| = %.2e ; v = %s"
          % (float(abs(a - b)), mp.nstr(b, 40)))
