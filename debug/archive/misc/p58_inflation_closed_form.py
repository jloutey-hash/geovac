r"""Closed-form for the Paper 58 permitted-density inflation factor R(n_max).

Target: replace the COUNTED 13.8x (n_max=2) / 15.0x (n_max=3) two-body
inflation with a closed form and settle whether R plateaus or grows.

Predicate is IDENTICAL to tests/test_paper58_census.py (validated here against
the two KNOWN test anchors: n_max=2 -> genuine/builder = 2944/214 = 13.8x;
n_max=3 -> 114280/7600 = 15.0x). Two centers A,B; orbitals (n,l,m),
1<=n<=n_max, 0<=l<n, -l<=m<=l. ERI quartet (p,q,r,s), charge-distribution
sides (p,q),(r,s):

  m-rule (global M_L):     m_p + m_r == m_q + m_s
  same-center side (x,y):  multipole L in side_L_range(l_x,l_y, m_y-m_x)
  builder: all four on ONE center AND L1 cap L2 != empty
  genuine: m-rule; each same-center side nonempty L; cross sides free
           (all-four-one-center still requires the intersect)

Hand-proved reductions (checked here):
  c(mu)      = (n-|mu|)(n-|mu|+1)/2                      [orbitals with that m]
  d(t)       = sum_m c(m) c(t-m)  = (c*c)(t)
  one_center = 2 * sum_t d(t)^2
  m_rule     = 16 * sum_t d(t)^2  = 8 * one_center
  genuine    = m_rule - one_center + builder = 7*one_center + builder
  R          = genuine/builder    = 1 + 7*(one_center/builder)
Nonempty-L removal never fires (|M| <= l_x+l_y always), so builder is the only
piece needing the (l,l',M) L-range intersection count.

Run: python debug/p58_inflation_closed_form.py
"""
from __future__ import annotations

from collections import defaultdict
from itertools import product

import sympy as sp

P = lambda *a: print(*a, flush=True)


# ---------------------------------------------------------------- brute (test-identical)
def basis(n_max):
    return [(n, l, m)
            for n in range(1, n_max + 1)
            for l in range(n)
            for m in range(-l, l + 1)]


def side_L_range(l1, l2, M):
    lo = max(abs(l1 - l2), abs(M))
    return [L for L in range(lo, l1 + l2 + 1) if (l1 + l2 + L) % 2 == 0]


def brute_census(n_max):
    orbs = ([("A",) + b for b in basis(n_max)]
            + [("B",) + b for b in basis(n_max)])
    counts = dict(dense=len(orbs) ** 4, m_rule=0, one_center=0,
                  genuine=0, builder=0)
    for p, q, r, s in product(orbs, repeat=4):
        if p[3] + r[3] != q[3] + s[3]:
            continue
        counts["m_rule"] += 1
        p_same, r_same = p[0] == q[0], r[0] == s[0]
        L1 = side_L_range(p[2], q[2], q[3] - p[3]) if p_same else None
        L2 = side_L_range(r[2], s[2], s[3] - r[3]) if r_same else None
        if p_same and r_same and p[0] == r[0]:
            counts["one_center"] += 1
            if set(L1) & set(L2):
                counts["builder"] += 1
                counts["genuine"] += 1
        else:
            if (L1 is None or L1) and (L2 is None or L2):
                counts["genuine"] += 1
    return counts


# ---------------------------------------------------------------- fast exact (factorized)
def c_count(mu, n_max):
    mu = abs(mu)
    if mu > n_max - 1:
        return 0
    k = n_max - mu
    return k * (k + 1) // 2


def sum_d_squared(n_max):
    ms = range(-(n_max - 1), n_max)
    tot = 0
    for t in range(-2 * (n_max - 1), 2 * (n_max - 1) + 1):
        d = sum(c_count(m, n_max) * c_count(t - m, n_max) for m in ms)
        tot += d * d
    return tot


def _pairs_by_transfer(n_max):
    """One center: transfer M -> {parity: list of (lo', top)}. L-set never empty."""
    orbs = basis(n_max)
    D = defaultdict(lambda: {0: [], 1: []})
    for (_nx, lx, mx) in orbs:
        for (_ny, ly, my) in orbs:
            Mt = my - mx
            top = lx + ly
            lo = max(abs(lx - ly), abs(Mt))
            loP = lo if (top - lo) % 2 == 0 else lo + 1
            D[Mt][top % 2].append((loP, top))
    return D


class _BIT:
    def __init__(self, n):
        self.n = n
        self.t = [0] * (n + 1)

    def add(self, i, v=1):
        i += 1
        while i <= self.n:
            self.t[i] += v
            i += i & (-i)

    def pref(self, i):  # sum over [0, i]
        i += 1
        s = 0
        while i > 0:
            s += self.t[i]
            i -= i & (-i)
        return s


def _count_pairs(A, B, topmax):
    """# ordered (a in A, b in B) with a.lo <= b.top AND b.lo <= a.top."""
    if not A or not B:
        return 0
    A = sorted(A, key=lambda x: x[1])          # by a.top asc
    B = sorted(B, key=lambda x: x[0])          # by b.lo asc
    bit = _BIT(topmax + 1)
    total = 0
    j = 0
    for lo_a, top_a in A:
        while j < len(B) and B[j][0] <= top_a:  # add all b with b.lo <= a.top
            bit.add(B[j][1])                     # index by b.top
            j += 1
        # among those, count b.top >= a.lo  ==  added - (# b.top < a.lo)
        below = bit.pref(lo_a - 1) if lo_a - 1 >= 0 else 0
        total += j - below
    return total


def builder_percenter(n_max):
    D = _pairs_by_transfer(n_max)
    topmax = 2 * (n_max - 1)
    tot = 0
    for Mt, byp in D.items():
        other = D.get(-Mt, {0: [], 1: []})
        for par in (0, 1):
            tot += _count_pairs(byp[par], other[par], topmax)
    return tot


def fast_counts(n_max):
    M = sum(n * n for n in range(1, n_max + 1))
    S = sum_d_squared(n_max)
    one_center = 2 * S
    m_rule = 16 * S
    builder = 2 * builder_percenter(n_max)
    genuine = m_rule - one_center + builder
    return dict(dense=(2 * M) ** 4, m_rule=m_rule, one_center=one_center,
                genuine=genuine, builder=builder)


# ---------------------------------------------------------------- identification
def finite_diff_degree(vals, maxdeg=16):
    seq = list(vals)
    for d in range(maxdeg + 1):
        if len(set(seq)) == 1:
            return d, seq[0]
        seq = [seq[i + 1] - seq[i] for i in range(len(seq) - 1)]
    return None, None


def fit_and_verify(name, ns, vals):
    """Try a single polynomial; if held-out fails, split by parity of n."""
    n = sp.Symbol('n')
    deg, _ = finite_diff_degree(vals)
    P(f"  [{name}] finite-difference degree = {deg}")
    if deg is not None and deg + 2 <= len(ns):
        k = deg + 1
        poly = sp.interpolate(list(zip(ns[:k], vals[:k])), n)
        poly = sp.expand(poly)
        ok = all(int(poly.subs(n, x)) == v for x, v in zip(ns, vals))
        if ok:
            P(f"  [{name}] single polynomial fits ALL {len(ns)} points:")
            P(f"          {poly}")
            return ("poly", poly)
    # parity split
    P(f"  [{name}] no single polynomial -> trying parity split")
    res = {}
    for par, tag in ((0, "even"), (1, "odd")):
        sub = [(x, v) for x, v in zip(ns, vals) if x % 2 == par]
        dp, _ = finite_diff_degree([v for _, v in sub])
        if dp is None or dp + 2 > len(sub):
            P(f"  [{name}/{tag}] insufficient points (deg={dp}, npts={len(sub)})")
            continue
        k = dp + 1
        poly = sp.expand(sp.interpolate(sub[:k], n))
        ok = all(int(poly.subs(n, x)) == v for x, v in sub)
        P(f"  [{name}/{tag}] deg={dp} fits all {len(sub)}: {ok}")
        P(f"          {poly}")
        res[tag] = poly
    return ("quasi", res)


def main():
    P("=== validation: fast == brute at the two KNOWN anchors ===")
    for nm in (2, 3):
        b, f = brute_census(nm), fast_counts(nm)
        P(f"  n_max={nm}: fast==brute? {b == f}")
        if b != f:
            P("   brute:" + str(b))
            P("   fast :" + str(f))
            raise SystemExit("FAST COUNTER DISAGREES -- stop.")
    assert fast_counts(2)["genuine"] == 2944 and fast_counts(2)["builder"] == 214
    assert fast_counts(3)["genuine"] == 114280 and fast_counts(3)["builder"] == 7600
    P("  anchors reproduced: (2944/214)=13.8x, (114280/7600)=15.0x")

    NMAX = 26
    ns = list(range(2, NMAX + 1))
    data = {k: [] for k in ("one_center", "builder", "genuine", "m_rule")}
    Rs, sd2 = [], []
    P(f"\n=== exact counts n_max = 2..{NMAX} ===")
    for nm in ns:
        f = fast_counts(nm)
        for k in data:
            data[k].append(f[k])
        sd2.append(sum_d_squared(nm))
        R = f["genuine"] / f["builder"]
        Rs.append(R)
        P(f"  n={nm:2d}  one_center={f['one_center']:>14d}  builder={f['builder']:>13d}"
          f"  R={R:9.5f}  oc/b={f['one_center']/f['builder']:9.6f}")

    P("\n=== closed-form identification ===")
    fit_and_verify("sum_d2", ns, sd2)
    kind_oc, oc = fit_and_verify("one_center", ns, data["one_center"])
    kind_b, bld = fit_and_verify("builder", ns, data["builder"])

    P("\n=== asymptotics of R = 1 + 7*(one_center/builder) ===")
    n = sp.Symbol('n')
    if kind_oc == "poly" and kind_b == "poly":
        ratio = sp.simplify(oc / bld)
        R_expr = sp.simplify(1 + 7 * ratio)
        P("  R(n) = " + str(R_expr))
        lim = sp.limit(R_expr, n, sp.oo)
        P("  lim_{n->oo} R = " + str(lim) + f"  (float {float(lim):.6f})")
        P(f"  deg(one_center)={sp.degree(sp.Poly(oc, n))}, "
          f"deg(builder)={sp.degree(sp.Poly(bld, n))}")
        for nm in (2, 3, 5, 10, 25, 100, 1000):
            P(f"     n={nm:4d}  R={float(R_expr.subs(n, nm)):.5f}")
    else:
        P("  (quasi-polynomial branch; using even-n and odd-n leading behavior)")
        for tag in ("even", "odd"):
            o = oc[tag] if kind_oc == "quasi" else oc
            b = bld[tag] if kind_b == "quasi" else bld
            R_expr = sp.simplify(1 + 7 * o / b)
            lim = sp.limit(R_expr, n, sp.oo)
            P(f"  [{tag}] R(n)={R_expr};  lim={lim} (float {float(lim):.6f})")


if __name__ == "__main__":
    main()
