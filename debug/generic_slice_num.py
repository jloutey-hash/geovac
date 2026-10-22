"""GS numerics (record docs/predictions/2026-10-21_generic_slice.md): G_b(t) = int x K1(t x) J(x) J(b x) dx
and theta_t^j G_b on ONE double-exponential grid x = exp(s - e^-s), trapezoid in s, with
f(x) = x J(x) J(b x) cached.  Two step sizes (h, h/2) give the quadrature error.
apply(Lc, t): sum_k Lc[k](t) theta^k G_b(t) for a theta-form coefficient list (low -> high), returns
(value, largest single term)."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
from three_electron_p3e_check import J as Jfun  # noqa: E402

HERE = Path(__file__).parent
y = sp.Symbol("y")


def theta_k1(n):
    """theta^j K1(y) = a_j(y) K1 + b_j(y) theta K1, theta^2 K1 = (1 + y^2) K1."""
    al, be = [sp.Integer(1)], [sp.Integer(0)]
    for _ in range(n):
        a, b = al[-1], be[-1]
        al.append(sp.expand(y * sp.diff(a, y) + b * (1 + y**2)))
        be.append(sp.expand(a + y * sp.diff(b, y)))
    return al, be


class Grid:
    def __init__(self, b, dps, hinv, smin=-5, smax=5.2):
        mp.mp.dps = dps
        self.b = mp.mpf(b)
        h = mp.mpf(1) / hinv
        self.nodes = []
        k = int(smin * hinv)
        while k * h <= smax:
            s = k * h
            x = mp.exp(s - mp.exp(-s))
            w = h * x * (1 + mp.exp(-s))
            self.nodes.append((x, w, k % 2 == 0))
            k += 1
        t0 = time.time()
        self.f = [x * Jfun(x) * Jfun(self.b * x) for x, _, _ in self.nodes]
        self.build_time = time.time() - t0

    def moments(self, t, n):
        """I_j = int f(x) theta_t^j K1(t x) dx for j <= n, fine and coarse (2h) rule."""
        t = mp.mpf(t)
        al, be = theta_k1(n)
        fa = [sp.lambdify(y, a, "mpmath") for a in al]
        fb = [sp.lambdify(y, b, "mpmath") for b in be]
        fine = [mp.mpf(0)] * (n + 1)
        coarse = [mp.mpf(0)] * (n + 1)
        for (x, w, even), fx in zip(self.nodes, self.f):
            yy = t * x
            k0, k1 = mp.besselk(0, yy), mp.besselk(1, yy)
            th1 = -yy * k0 - k1
            for j in range(n + 1):
                v = fx * (fa[j](yy) * k1 + fb[j](yy) * th1) * w
                fine[j] += v
                if even:
                    coarse[j] += 2 * v
        return fine, coarse


def apply(grid, Lc, t, cache={}):
    key = (id(grid), str(t), len(Lc))
    if key not in cache:
        cache[key] = grid.moments(t, len(Lc) - 1)
    fine, coarse = cache[key]
    tv = mp.mpf(t)
    tt = sp.Symbol("t")
    cs = [mp.mpf(sp.N(sp.sympify(c).subs(tt, sp.Rational(str(t))), mp.mp.dps + 10)) for c in Lc]
    terms = [c * I for c, I in zip(cs, fine)]
    val = sum(terms)
    valc = sum(c * I for c, I in zip(cs, coarse))
    return val, max(abs(x) for x in terms), abs(val - valc)
