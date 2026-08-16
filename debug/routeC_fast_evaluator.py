"""Fast high-precision evaluator for the collinear integrated T2 (Rung-3 frontier,
obstacle ii).  Replaces the slow adaptive nested quadrature (times out >260s at dps>=18)
with a FIXED tensor Gauss-Legendre grid + sin^2 substitution (kills the sqrt(s) endpoint
non-analyticity so GL converges fast) + s<->t symmetry (halves the fiber evaluations).

Collinear geometry: X=0, Y=(0,0,1), Z=(0,0,-1), 1s zeta=1 => D1=D2=1, |W|=s+t.
    T2 = (8/pi) int_0^1 ds int_0^1 dt  J(s,t),   J symmetric in s,t
    J(s,t) = int_0^inf dk  j0(k(s+t)) P(s,k) P(t,k)
    P(x,k) = c e^{-Delta}(1/Delta^3 + 3/Delta^4 + 3/Delta^5),  c=x(1-x), Delta=sqrt(c k^2+1)
Substitution s=sin^2(phi), ds=sin(2phi) dphi, phi in [0, pi/2] (smooths BOTH endpoints).
"""
from __future__ import annotations
import sys
import time
import mpmath as mp


def P(x, k):
    c = x * (1 - x)
    Del = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-Del) * (1 / Del ** 3 + 3 / Del ** 4 + 3 / Del ** 5)


def J(s, t):
    b = s + t
    def f(k):
        j0 = mp.sin(k * b) / (k * b) if k * b > mp.mpf('1e-30') else mp.mpf(1)
        return j0 * P(s, k) * P(t, k)
    return mp.quad(f, [0, 1, 3, 8, 20, mp.inf])


def gl(N):
    # mpmath has no direct gauss-legendre; build from Legendre poly roots + weights.
    from mpmath import legendre, diff
    # roots of P_N via Newton from cos((k-0.25)pi/(N+0.5)) seeds
    roots = []
    for k in range(1, N + 1):
        x = mp.cos(mp.pi * (k - mp.mpf('0.25')) / (N + mp.mpf('0.5')))
        for _ in range(80):
            f = legendre(N, x)
            fp = N * (x * legendre(N, x) - legendre(N - 1, x)) / (x * x - 1)   # P_N'(x)
            dx = f / fp
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        roots.append(x)
    ws = []
    for x in roots:
        fp = N * (x * legendre(N, x) - legendre(N - 1, x)) / (x * x - 1)       # P_N'(x)
        ws.append(2 / ((1 - x * x) * fp * fp))
    return roots, ws


def T2(N):
    xs, ws = gl(N)
    half = mp.pi / 4
    phi = [half * (x + 1) for x in xs]
    wphi = [half * w for w in ws]
    s = [mp.sin(p) ** 2 for p in phi]
    jac = [mp.sin(2 * p) for p in phi]
    tot = mp.mpf(0)
    for i in range(N):
        wi = wphi[i] * jac[i]
        tot += wi * wi * J(s[i], s[i])
        for j in range(i):
            tot += 2 * wi * (wphi[j] * jac[j]) * J(s[i], s[j])
    return (8 / mp.pi) * tot


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    Ns = [int(a) for a in sys.argv[2:]] or [40, 56]
    mp.mp.dps = dps
    print(f"Fast evaluator: dps={dps}, GL tensor + sin^2 sub + symmetry")
    prev = None
    for N in Ns:
        t0 = time.time()
        v = T2(N)
        dt = time.time() - t0
        agree = "" if prev is None else f"  |d vs prev N| {mp.nstr(abs(v-prev),3)}"
        print(f"  N={N:3d}: {mp.nstr(v, dps-3)}  ({dt:.1f}s){agree}")
        prev = v


if __name__ == '__main__':
    main()
