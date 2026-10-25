"""IB record control: the kite with the Coulomb line removed factorizes into two one-loop bubbles,
X(z) = B_rr(z) B_ff(z), indices (3/2, 3/2) masses (a, a) and (2, 2) masses (g, g), D = 3.
Feynman parameters: B_{nu1 nu2}(z; m) = pi^(3/2) Gamma(nu1+nu2-3/2)/(Gamma(nu1)Gamma(nu2))
  int_0^1 x^(nu1-1)(1-x)^(nu2-1) (m^2 + x(1-x) z)^(3/2-nu1-nu2) dx,
so its Taylor coefficients are binomials times Beta functions: exact, B_rr in pi Q, B_ff in pi^(3/2) Q
up to the constant, hence X_n / X_0 in Q.  The reducer's operator for X (shift (1, 1, 1, 1, -1)) must
annihilate that exact series modulo P; a wrong operator (gamma changed) must not.
usage: python ib_sector_control.py [g] [gwrong]"""
from __future__ import annotations

import random
import sys
import time
from fractions import Fraction as F
from math import comb
from multiprocessing import Pool

from flint import nmod_poly

from ib_reduce import P, Reducer, shifts_box, dz_terms, fr
from ib_kite import kite_family
from ib_kite_de import fit_rational, RF, cyclic_operator

XT = (1, 1, 1, 1, -1)


def gbinom(alpha: F, n: int) -> F:
    out = F(1)
    for k in range(n):
        out = out * (alpha - k) / (k + 1)
    return out


def beta_half(n):
    """B(n + 3/2, n + 3/2) / pi as a rational: Gamma(n+3/2)^2 / Gamma(2n+3) / pi."""
    g = F(1, 2)  # Gamma(3/2)/sqrt(pi)
    for k in range(n):
        g = g * (F(3, 2) + k)
    from math import factorial
    return g * g / factorial(2 * n + 2)


def beta_int(n):
    from math import factorial
    return F(factorial(n + 1) ** 2, factorial(2 * n + 3))


def series_X(g: F, N: int):
    """X_n up to a common constant: B_rr coefficients binom(-3/2, n) B(n+3/2, n+3/2) (m = 1),
    B_ff coefficients binom(-5/2, n) g^(-5-2n) B(n+2, n+2)."""
    brr = [gbinom(F(-3, 2), n) * beta_half(n) for n in range(N)]
    bff = [gbinom(F(-5, 2), n) * beta_int(n) / g ** (5 + 2 * n) for n in range(N)]
    return [sum(brr[i] * bff[n - i] for i in range(n + 1)) for n in range(N)]


def sample(args):
    zi, g = args
    fam = kite_family(F(1), F(1), g, g)
    fam.setup({("w", "w"): zi}, 3)
    seeds = [n for n in shifts_box(5, 5, 1, fam.integer) if not fam.is_zero(n)]
    Rd = Reducer(fam)
    Rd.seed(seeds, sym=lambda n: (n[1], n[0], n[3], n[2], n[4]))
    c = Rd.reduce({XT: 1})
    masters, rows, todo = set(c), {}, list(c)
    while todo:
        m = todo.pop()
        r = Rd.reduce(dz_terms(fam, m, "w", zi))
        rows[m] = r
        for n in r:
            if n not in masters:
                masters.add(n)
                todo.append(n)
    ms = sorted(masters, key=Rd.key)
    return zi, ms, {m: c.get(m, 0) for m in ms}, {m: rows[m] for m in ms}


def operator_for(g: F, nsamp=60):
    rng = random.Random(99)
    zs = [rng.randrange(10**6, P - 1) for _ in range(nsamp)]
    with Pool(6, maxtasksperchild=1) as pool:
        res = pool.map(sample, [(z, g) for z in zs], chunksize=1)
    assert len({tuple(r[1]) for r in res}) == 1
    ms = res[0][1]
    fitp, testp = res[:-8], res[-8:]
    xs = [r[0] for r in fitp]

    def rebuild(get):
        ys = [get(r) for r in fitp]
        if not any(ys):
            return RF(nmod_poly([0], P))
        N, D = fit_rational(xs, ys, len(xs) - 3)
        for r in testp:
            assert (N(r[0]) - D(r[0]) * get(r)) == 0
        return RF(N, D)

    c = [rebuild(lambda r, m=m: r[2][m]) for m in ms]
    A = [[rebuild(lambda r, m=m, n=n: r[3][m].get(n, 0)) for n in ms] for m in ms]
    a = cyclic_operator(c, A)
    den = nmod_poly([1], P)
    for x in a:
        den = den * x.d // den.gcd(x.d)
    return len(ms), [x.n * (den // x.d) for x in a]


def residual(polys, xs_series):
    """max over N of |sum_k sum_i p_{k,i} (N-i+k)!/(N-i)! x_{N-i+k}| mod P (0 = annihilated)."""
    Nmax = len(xs_series)
    x = [fr(v) for v in xs_series]
    order = len(polys) - 1
    bad = 0
    for N in range(0, Nmax - order - 1):
        acc = 0
        for k, p in enumerate(polys):
            co = [int(c) for c in p.coeffs()]
            for i, c in enumerate(co):
                m = N - i
                if m < 0 or not c or m + k >= Nmax:
                    continue
                ff = 1
                for j in range(k):
                    ff *= (m + k - j)
                acc = (acc + c * ff % P * x[m + k]) % P
        if acc:
            bad += 1
    return bad


def main():
    g = F(sys.argv[1]) if len(sys.argv) > 1 else F(3, 2)
    gw = F(sys.argv[2]) if len(sys.argv) > 2 else g + F(1, 1000)
    t0 = time.time()
    xs = series_X(g, 120)
    nm, L = operator_for(g)
    print(f"masters {nm}, operator order {len(L)-1}, degrees {[p.degree() for p in L]} ({time.time()-t0:.0f}s)")
    print("residual count on the exact series X(g):", residual(L, xs), "of", 120 - len(L))
    xw = series_X(gw, 120)
    print("residual count on X at gamma + 1/1000 (must be nonzero):", residual(L, xw))
    Lw = operator_for(gw)[1]
    print("operator built at gamma + 1/1000 applied to X(g) (must be nonzero):", residual(Lw, xs))


if __name__ == "__main__":
    main()
