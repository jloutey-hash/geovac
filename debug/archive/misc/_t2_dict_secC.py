"""Section C of the T2 Euclidean dictionary, with fixed Gauss-Legendre grids so the
nested (position-space) integrals are affordable.

C0  j0(k|W|) = (1/2|W|) int_{-|W|}^{|W|} cos(kv) dv                    [exact kernel identity]
C1  f_i(b) = int_0^inf cos(kb)P(x_i,k)dk  ==  sqrt(c) int_1^inf w(u) u K_1(R)/R du
            (the position-space / propagator form; = -d/dp K_0 of the D1 propagator, u-averaged)
C2  UNSMEARED physical overlap (the |W|->0 fibre), both factors in PROPAGATOR form:
        int_0^inf P(s,k)P(t,k) dk = (2/pi) int_0^inf f_1(b) f_2(b) db
C3  FULL physical fibre with the j0 box smearing:
        J(s,t) = (2/pi) int_0^inf db f_1(b) <f_2>_{|W|}(b)
"""
from __future__ import annotations
import sys
import time
import mpmath as mp

sys.path.insert(0, 'debug')
from t2_euclidean_dictionary import P, w_u, f_direct, J_direct, g_of  # noqa: E402

_GL = {}


def gl(N):
    key = (N, mp.mp.dps)
    if key in _GL:
        return _GL[key]
    roots, ws = [], []
    for k in range(1, N + 1):
        x = mp.cos(mp.pi*(k - mp.mpf('0.25'))/(N + mp.mpf('0.5')))
        for _ in range(60):
            f = mp.legendre(N, x)
            fp = N*(x*mp.legendre(N, x) - mp.legendre(N - 1, x))/(x*x - 1)
            dx = f/fp
            x -= dx
            if abs(dx) < mp.mpf(10)**(-mp.mp.dps - 5):
                break
        fp = N*(x*mp.legendre(N, x) - mp.legendre(N - 1, x))/(x*x - 1)
        roots.append(x)
        ws.append(2/((1 - x*x)*fp*fp))
    _GL[key] = (roots, ws)
    return roots, ws


def glint(f, a, b, N):
    xs, ws = gl(N)
    h, m = (b - a)/2, (b + a)/2
    return h*sum(w*f(m + h*x) for x, w in zip(xs, ws))


def f_prop_fast(x, b, N=44):
    """sqrt(c) int_1^inf w(u) u K_1(R)/R du, R = sqrt(u^2 + b^2/c)  -- fixed GL, 3 panels."""
    c = x*(1 - x)
    b2c = b*b/c

    def h(u):
        R = mp.sqrt(u*u + b2c)
        return w_u(u)*u*mp.besselk(1, R)/R
    return mp.sqrt(c)*(glint(h, 1, 4, N) + glint(h, 4, 14, N) + glint(h, 14, 55, N))


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    mp.mp.dps = dps
    print(f"Section C -- position-space (propagator) form of the fibre, dps={dps}")

    # ---- C0
    k, bW = mp.mpf('1.7'), mp.mpf('0.9')
    d = abs(mp.sin(k*bW)/(k*bW) - mp.quad(lambda v: mp.cos(k*v), [-bW, bW])/(2*bW))
    print(f"\nC0  kernel identity j0 = box average:  |d| = {mp.nstr(d, 3)}", flush=True)

    # ---- C1
    print("\nC1  f_i(b): cosine transform (momentum) vs sqrt(c) int w(u) u K_1(R)/R du (position)")
    worst = mp.mpf(0)
    for xs, bs in [('0.3', '0.5'), ('0.2', '1.2'), ('0.45', '0.05'), ('0.1', '2.0'), ('0.5', '0.8')]:
        x, b = mp.mpf(xs), mp.mpf(bs)
        v1, v2 = f_direct(x, b), f_prop_fast(x, b)
        rd = abs(v1 - v2)/abs(v1)
        worst = max(worst, rd)
        print(f"    x={xs:>5} b={bs:>5}  f={mp.nstr(v1, 18)}   rel.d={mp.nstr(rd, 3)}", flush=True)
    print(f"    worst: {mp.nstr(worst, 3)}")

    # ---- C2 : unsmeared overlap, both factors in propagator form
    print("\nC2  int_0^inf P(s,k)P(t,k)dk  ==  (2/pi) int_0^inf f_1(b) f_2(b) db   [both f in K_1 form]")
    worst2 = mp.mpf(0)
    for ss, ts in [('0.3', '0.2'), ('0.4', '0.4'), ('0.25', '0.6'), ('0.15', '0.45')]:
        s, t = mp.mpf(ss), mp.mpf(ts)
        c1, c2 = s*(1 - s), t*(1 - t)
        A = mp.sqrt(c1) + mp.sqrt(c2)
        lhs = mp.quad(lambda k: P(s, k)*P(t, k), [0, 1/A, 2/A, 4/A, 8/A, 20/A, mp.inf])
        a1 = 1/mp.sqrt(c1)
        Bmax = 32/a1
        t0 = time.time()
        rhs = (2/mp.pi)*(glint(lambda b: f_prop_fast(s, b)*f_prop_fast(t, b), 0, Bmax/4, 24)
                         + glint(lambda b: f_prop_fast(s, b)*f_prop_fast(t, b), Bmax/4, Bmax, 24))
        rd = abs(lhs - rhs)/abs(lhs)
        worst2 = max(worst2, rd)
        print(f"    s={ss} t={ts}  momentum={mp.nstr(lhs, 18)}  position={mp.nstr(rhs, 18)}"
              f"  rel.d={mp.nstr(rd, 3)}  [{time.time()-t0:.0f}s]", flush=True)
    print(f"    worst: {mp.nstr(worst2, 3)}")

    # ---- C3 : full fibre with the box smearing
    print("\nC3  J(s,t) = (2/pi) int_0^inf db f_1(b) <f_2>_{|W|}(b)   [f_1 propagator, f_2 propagator]")
    worst3 = mp.mpf(0)
    for ss, ts in [('0.3', '0.2'), ('0.4', '0.4')]:
        s, t = mp.mpf(ss), mp.mpf(ts)
        bW = s + t
        c1 = s*(1 - s)
        a1 = 1/mp.sqrt(c1)
        Jd = J_direct(s, t)
        Bmax = 32/a1 + 2*bW
        t0 = time.time()

        def smeared(b):
            return glint(lambda v: f_prop_fast(t, abs(b + v)), -bW, bW, 18)/(2*bW)
        Jp = (2/mp.pi)*(glint(lambda b: f_prop_fast(s, b)*smeared(b), 0, Bmax/4, 20)
                        + glint(lambda b: f_prop_fast(s, b)*smeared(b), Bmax/4, Bmax, 20))
        rd = abs(Jd - Jp)/abs(Jd)
        worst3 = max(worst3, rd)
        print(f"    s={ss} t={ts} |W|={mp.nstr(bW, 6)}  J_momentum={mp.nstr(Jd, 18)}"
              f"  J_position={mp.nstr(Jp, 18)}  rel.d={mp.nstr(rd, 3)}  [{time.time()-t0:.0f}s]",
              flush=True)
    print(f"    worst: {mp.nstr(worst3, 3)}")


if __name__ == '__main__':
    main()
