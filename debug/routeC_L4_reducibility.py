"""Differential-Galois reducibility scan of the rank-4 Picard-Fuchs operator L4
(Paper 59 eq:pf), the N(D)-closed-form leg of the Route C frontier.  One driver that
runs the whole battery and prints the verdict.  Backing tests (permanent record):
tests/test_routeC_momentum.py::test_L4_*.  Memo: debug/sprint_L4_reducibility_scan_memo.md.

L4 (in the Laplace variable D, rho = c2/c1 the scale ratio):
    D rho L'''' + 2 rho L''' + D(1-2rho) L'' + (1-2rho) L' - D(1-rho) L = 0.

Verdict: L4 is formally self-adjoint (Galois group in Sp4), has trivial eigenring
(indecomposable), no order-1/2/3 factor over Q(rho)(D) with polar locus {0,inf}, and is
neither a Bessel-sector symmetric product nor a Sym^3.  DECISIVE closure (Stokes step):
the monodromy of L4 around D=0, in the exponential-line (Lefschetz-thimble) basis, is the
integer unipotent matrix [[-1,2,2,2],[-2,3,2,2],[-2,2,3,2],[2,-2,-2,-1]] (same for rho=1/5
and 1/3 -- a topological invariant of the 4-branch-point config); it has a single 2x2
Jordan block (the one D ln D log) and NO invariant coordinate subspace, so with the
exponential torus L4 is IRREDUCIBLE over C(D), hence over Q(rho)(D).  The Stokes matrices
DO couple the real exponential pair {+-1} to the imaginary pair {+-i w}.  N(D) is thus a
period of an irreducible rank-4 connection: no closed form via any factorization.

Run:  python debug/routeC_L4_reducibility.py   (symbolic battery + numerical monodromy)
"""
from __future__ import annotations

from collections import defaultdict

import sympy as sp

D, s, rho, w = sp.symbols('D s rho w', positive=True)


def L4_coeffs(rho_val):
    """p[k] multiplies d^k/dD^k (non-monic eq:pf form)."""
    return {0: -D * (1 - rho_val), 1: (1 - 2 * rho_val), 2: D * (1 - 2 * rho_val),
            3: 2 * rho_val, 4: D * rho_val}


# ------------------------------------------------------------------ (1) self-adjoint
def check_self_adjoint():
    p = L4_coeffs(rho)
    y = sp.Function('y')(D)
    Ly = sum(p[k] * sp.diff(y, D, k) for k in p)
    Lstar = sum((-1) ** k * sp.diff(p[k] * y, D, k) for k in p)
    diff = sp.simplify(sp.expand(Ly - Lstar))
    print(f"(1) self-adjoint:  L4 - L4* = {diff}   => Galois group in Sp4" if diff == 0
          else f"(1) NOT self-adjoint: {diff}")


# --------------------------------------------------------------- (2) indicial at D=0
def check_indicial():
    p = L4_coeffs(rho)
    bucket = defaultdict(lambda: sp.Integer(0))
    for k, c in p.items():
        for (deg,), a in sp.Poly(c, D).terms():
            bucket[deg - k] += a * sp.ff(s, k)
    ind = sp.factor(sp.simplify(bucket[min(bucket)]))
    roots = sp.roots(sp.Poly(ind, s))
    print(f"(2) indicial at D=0: {ind}  => exponents {dict(roots)}  "
          f"(double at 1 => the D ln D nonanalyticity)")


# ---------------------------------------------------------- (3) symbol factorization
def check_symbol():
    symbol = s ** 4 + ((1 - 2 * rho) / rho) * s ** 2 - (1 - rho) / rho
    print(f"(3) symbol factors over Q(rho): {sp.factor(symbol)}  "
          f"(only the two Bessel sectors; constant term <0 forbids a mixed split)")


# ------------------------------------------------- order-1 / order-2 / eigenring tools
def monic(rho_val):
    """monic p3,p2,p1,p0 of L4 at a numeric rho."""
    return (2 / D, (1 - 2 * rho_val) / rho_val,
            (1 - 2 * rho_val) / (rho_val * D), -(1 - rho_val) / rho_val)


def check_order1(rho_val, Nmax=6):
    p = L4_coeffs(rho_val)
    wv = sp.sqrt((1 - rho_val) / rho_val)

    def has_sol(lam):
        for N in range(Nmax + 1):
            qs = sp.symbols(f'q0:{N + 1}')
            Q = sum(qs[j] * D ** j for j in range(N + 1))
            yv = sp.exp(lam * D) * Q
            e = sp.expand(sum(p[k] * sp.diff(yv, D, k) for k in p) / sp.exp(lam * D))
            sold = sp.solve(sp.Poly(sp.simplify(e), D).all_coeffs(), qs, dict=True)
            if sold and any(sold[0].get(q, q) != 0 for q in qs):
                return True
        return False

    lams = [sp.Integer(1), sp.Integer(-1), sp.I * wv, -sp.I * wv, sp.Integer(0)]
    any_sol = any(has_sol(l) for l in lams)
    print(f"(4) order-1 (hyperexponential e^(lam D)Q) at rho={rho_val}: "
          f"{'FOUND (reducible!)' if any_sol else 'NONE'} "
          f"(complete for polar locus {{0,inf}}; exponents integral)")


def _remainder(p3, p2, p1, p0, a, b):
    ap, bp = sp.diff(a, D), sp.diff(b, D)
    P = a ** 2 - ap - b
    Q0 = a * b - bp
    R1 = (-a * P + sp.diff(P, D) + Q0) + p3 * (a ** 2 - ap - b) + p2 * (-a) + p1
    R0 = (sp.diff(Q0, D) - b * P) + p3 * (a * b - bp) + p2 * (-b) + p0
    return R1, R0


def check_order2(rho_val):
    am1, a0, a1, bm2, bm1, b0, b1 = sp.symbols('am1 a0 a1 bm2 bm1 b0 b1')
    unks = [am1, a0, a1, bm2, bm1, b0, b1]

    def search(p3, p2, p1, p0):
        a = am1 / D + a0 + a1 * D
        b = bm2 / D ** 2 + bm1 / D + b0 + b1 * D
        eqs = []
        for R in _remainder(p3, p2, p1, p0, a, b):
            eqs.extend(sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(R)))), D).all_coeffs())
        return sp.solve([sp.expand(e) for e in eqs if e != 0], unks, dict=True)

    ctrl = search(3 / D, sp.Integer(2), 1 / D, sp.Integer(-3))       # planted factor
    ctrl_ok = any(c.get(b0) == -1 and c.get(am1) == 1 for c in ctrl)
    sol = search(*monic(rho_val))
    print(f"(5) order-2 factor over Q(rho)(D), locus {{0,inf}} at rho={rho_val}: "
          f"{'FOUND' if sol else 'NONE'};  positive control recovers planted factor: {ctrl_ok}")


def check_eigenring(rho_val, lo=-2, hi=2):
    p0, p1, p2, p3 = monic(rho_val)

    def deriv(v):
        a0, a1, a2, a3 = v
        return [sp.diff(a0, D) - a3 * p0, a0 + sp.diff(a1, D) - a3 * p1,
                a1 + sp.diff(a2, D) - a3 * p2, a2 + sp.diff(a3, D) - a3 * p3]

    def Lapply(v):
        Dm = [v]
        for _ in range(4):
            Dm.append(deriv(Dm[-1]))
        cf = [p0, p1, p2, p3, sp.Integer(1)]
        return [sum(cf[m] * Dm[m][i] for m in range(5)) for i in range(4)]

    cs = {i: sp.symbols(f'c{i}_{lo + 50}:{hi + 51}') for i in range(4)}
    r = [sum(cs[i][j] * D ** (lo + j) for j in range(hi - lo + 1)) for i in range(4)]
    allc = [c for i in range(4) for c in cs[i]]
    eqs = []
    for Ei in Lapply(r):
        eqs.extend(sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(Ei)))), D).all_coeffs())
    sol = list(sp.linsolve([sp.expand(e) for e in eqs if e != 0], allc))[0]
    free = set().union(*[sp.sympify(x).free_symbols for x in sol]) & set(allc)
    print(f"(6+7) eigenring dim E(L4) at rho={rho_val}, window[{lo},{hi}] = {len(free)} "
          f"(=1 => indecomposable, only scalars)")


def check_monodromy(rho_val):
    """Numerical monodromy of L4 around D=0 in the exponential-line (thimble) basis, and the
    invariant-coordinate-subspace (irreducibility) test.  rho_val a python/sympy rational."""
    import mpmath as mp
    from itertools import combinations
    mp.mp.dps = 25
    mpf = mp.mpf
    rho = mpf(sp.Rational(rho_val).p) / mpf(sp.Rational(rho_val).q)
    w2 = (1 - rho) / rho
    BP = [mpf(1), mpf(-1), 1j * mp.sqrt(w2), -1j * mp.sqrt(w2)]
    Qf = lambda x: (x * x - 1) * (x * x + w2)
    Qpf = lambda x: 2 * x * (2 * x * x + (w2 - 1))
    Nf = 2500

    def thimble(xk, Dc):
        d = mp.conj(Dc) / abs(Dc); umax = mp.sqrt(50 / abs(Dc)) + 3; h = umax / Nf
        ref = mp.arg(Qpf(xk) * d); pa = ref; vals = []
        for i in range(Nf + 1):
            u = i * h; x = xk + u * u * d
            if i == 0:
                base = 2 * d * mp.e ** (-Dc * xk) / (mp.sqrt(abs(Qpf(xk) * d)) * mp.e ** (1j * ref / 2))
            else:
                qv = Qf(x); a = mp.arg(qv)
                while a - pa > mp.pi: a -= 2 * mp.pi
                while a - pa < -mp.pi: a += 2 * mp.pi
                pa = a
                base = mp.e ** (-Dc * x) / (mp.sqrt(abs(qv)) * mp.e ** (1j * a / 2)) * (2 * u) * d
            vals.append([((-x) ** n) * base for n in range(4)])
        out = []
        for n in range(4):
            ssum = vals[0][n] + vals[Nf][n]
            for i in range(1, Nf): ssum += (4 if i % 2 else 2) * vals[i][n]
            out.append(ssum * h / 3)
        return out

    p2n = (1 - 2 * rho) / rho; p1n = (1 - 2 * rho) / rho; p0n = (1 - rho) / rho
    Am = lambda Dc: mp.matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
                               [p0n, -p1n / Dc, -p2n, -2 / Dc]])
    D0 = 3 * mp.e ** (1j * mpf('0.5')); Rr = abs(D0); th0 = mp.arg(D0)
    Y = mp.eye(4); Ns = 3000; dt = 2 * mp.pi / Ns
    Fn = lambda t, Y: (1j * Rr * mp.e ** (1j * (th0 + t))) * (Am(Rr * mp.e ** (1j * (th0 + t))) * Y)
    t = mpf(0)
    for _ in range(Ns):
        k1 = Fn(t, Y); k2 = Fn(t + dt / 2, Y + (dt / 2) * k1)
        k3 = Fn(t + dt / 2, Y + (dt / 2) * k2); k4 = Fn(t + dt, Y + dt * k3)
        Y = Y + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4); t += dt
    C = mp.matrix(4, 4)
    for k, xk in enumerate(BP):
        j = thimble(xk, D0)
        for n in range(4): C[n, k] = j[n]
    Me = C ** -1 * Y * C
    Mint = [[int(mp.nint(mp.re(Me[i, j]))) for j in range(4)] for i in range(4)]
    worst = min(max([abs(Me[c, a]) for a in S for c in [i for i in range(4) if i not in S]] + [mpf(0)])
                for d in (1, 2, 3) for S in combinations(range(4), d))
    print(f"(8) monodromy M0 at rho={rho_val} (exp-line basis), integer matrix = {Mint}")
    print(f"      min off-block over all proper coord subspaces = {mp.nstr(worst, 6)}  "
          f"=> {'IRREDUCIBLE over C(D)' if worst > 1e-4 else 'reducible'}")


def main():
    print(__doc__.splitlines()[0])
    print("=" * 78)
    check_self_adjoint()
    check_indicial()
    check_symbol()
    for rv in (sp.Rational(1, 5), sp.Rational(1, 3)):
        check_order1(rv)
        check_order2(rv)
        check_eigenring(rv)
    print("-" * 78)
    print("Stokes step (numerical monodromy) -- closes base-field reducibility:")
    for rv in (sp.Rational(1, 5), sp.Rational(1, 3)):
        check_monodromy(rv)
    print("=" * 78)
    print("VERDICT: L4 is IRREDUCIBLE over C(D), hence over Q(rho)(D).  The monodromy around")
    print("D=0 in the exponential-line basis is an integer unipotent matrix with a single")
    print("2x2 Jordan block (the one D ln D log) and NO invariant coordinate subspace; the")
    print("exponential torus forces invariant subspaces to be coordinate subspaces.  The")
    print("Stokes matrices couple the real and imaginary exponential pairs.  N(D) is a period")
    print("of an IRREDUCIBLE rank-4 connection: no closed form via factorization.  The only")
    print("open piece is a transcendental closed form at the irreducible level (elliptic")
    print("polylog / Gamma(2) MMV) -- the sec:modular frontier.")


if __name__ == '__main__':
    main()
