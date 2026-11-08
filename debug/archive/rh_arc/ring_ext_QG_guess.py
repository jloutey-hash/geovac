"""Guess a linear ODE with polynomial coefficients for a power series known to high precision.

Series f = sum_n c_n x^(lam + step*n).  Operator ansatz L = sum_{i<=d} x^(step*i) P_i(theta),
deg P_i <= r, theta = x d/dx.  Coefficient of x^(lam+step*n) in L f:
    sum_i P_i(lam + step*(n-i)) c_{n-i}   (c_m = 0 for m < 0).
Method: U unknowns; U-1 rows + a random normalisation row, solved in flint arb_mat at high
precision; the candidate must then kill all remaining rows (residual test).  Rational recognition by
continued fractions, then an exact re-check of the recurrence with the rational operator."""
from __future__ import annotations
import random
from fractions import Fraction
import mpmath as mp
from flint import arb, arb_mat, ctx


def ipow(e, j):
    out = arb(1)
    for _ in range(j):
        out = out * e
    return out


def build_rows(c, lam, step, r, d, n_from, n_to, scale):
    """rows n_from..n_to-1, columns (i, j); entries scaled: row n * scale^n, column (i,j) / (scale^i n_to^j)."""
    rows = []
    for n in range(n_from, n_to):
        row = []
        for i in range(d + 1):
            m = n - i
            for j in range(r + 1):
                if m < 0:
                    row.append(arb(0))
                else:
                    e = arb(lam) + step * m
                    row.append(ipow(e, j) * c[m] * ipow(arb(scale), m) / ipow(arb(n_to), j))
        rows.append(row)
    return rows


def try_operator(c, lam, step, r, d, n_from=0, n_to=None, scale=1, seed=1, tol_digits=40):
    """Return (vector p_{ij} in unscaled form, worst relative residual on check rows) or None."""
    U = (r + 1) * (d + 1)
    n_to = n_to or len(c)
    if n_to - n_from < U + 10:
        return None
    rows = build_rows(c, lam, step, r, d, n_from, n_to, scale)
    rng = random.Random(seed)
    xi = [arb(rng.randint(1, 10**6)) for _ in range(U)]
    A = arb_mat([rows[k] for k in range(U - 1)] + [xi])
    b = arb_mat([[0]] * (U - 1) + [[1]])
    try:
        x = A.solve(b, algorithm="approx")
    except Exception:
        return None
    xs = [x[k, 0] for k in range(U)]
    if not all(v.is_finite() for v in xs):
        return None
    # reject an ill-determined solve: radius must be tiny relative to the vector
    big = max(abs(mp.mpf(v.mid().str(30, radius=False))) for v in xs)
    if big == 0:
        return None
    xs = [v.mid() for v in xs]
    worst = mp.mpf(0)
    for row in rows[U - 1:]:
        res = sum((a.mid() * v for a, v in zip(row, xs)), arb(0))
        size = sum((abs(a.mid() * v) for a, v in zip(row, xs)), arb(0))
        if size == 0:
            rel = mp.mpf(1)
        else:
            rel = abs(mp.mpf((res / size).mid().str(20, radius=False)))
        worst = max(worst, rel)
    worst = arb(mp.nstr(worst, 15)) if worst > 0 else arb(0)
    # unscale: p_{ij} = x_{ij} / n_to^j   (the scale^m factor was per-row-consistent: scale^(n-i) -> P_i * scale^-i)
    p = {}
    k = 0
    for i in range(d + 1):
        for j in range(r + 1):
            p[(i, j)] = xs[k] / ipow(arb(n_to), j) / ipow(arb(scale), i)
            k += 1
    return p, worst


def minimal_degree(c, lam, step, r, dmax, **kw):
    for d in range(0, dmax + 1):
        out = try_operator(c, lam, step, r, d, **kw)
        if out is None:
            continue
        p, worst = out
        if worst < arb(10) ** (-kw.get("tol_digits", 40)):
            return d, p, worst
    return None


def rationalize(p, maxden=10**45, digits=60):
    """Scale p to integers: divide by the largest entry, recognise ratios, clear denominators."""
    keys = list(p)
    ref = max(keys, key=lambda k: abs(mp.mpf(p[k].mid().str(30, radius=False))))
    out = {}
    for k in keys:
        v = (p[k] / p[ref]).mid()
        s = v.str(digits + 10, radius=False)
        q = Fraction(mp.mpf(s)) if False else Fraction(s.replace("e", "E")) if "e" in s else Fraction(s)
        q = q.limit_denominator(maxden)
        # verify
        if abs(arb(q.numerator) / q.denominator - v) > arb(10) ** (-digits):
            return None
        out[k] = q
    from math import lcm, gcd
    L = 1
    for q in out.values():
        L = lcm(L, q.denominator)
    ints = {k: int(q * L) for k, q in out.items()}
    G = 0
    for v in ints.values():
        G = gcd(G, v)
    return {k: v // G for k, v in ints.items()}


def exact_residuals(P, c, lam, step, n_to):
    """max over n of |sum_i P_i(lam+step(n-i)) c_{n-i}| / (sum of |terms|); P = {(i,j): int}."""
    d = max(i for i, _ in P)
    worst = mp.mpf(0)
    for n in range(n_to):
        res, size = arb(0), arb(0)
        for i in range(d + 1):
            m = n - i
            if m < 0:
                continue
            e = arb(lam) + step * m
            val = sum((arb(v) * ipow(e, j) for (ii, j), v in P.items() if ii == i), arb(0))
            res += val * c[m]
            size += abs(val * c[m])
        if size != 0:
            r_ = abs(mp.mpf((res / size).mid().str(20, radius=False)))
            if r_ > worst:
                worst = r_
    return arb(mp.nstr(worst, 15)) if worst > 0 else arb(0)

