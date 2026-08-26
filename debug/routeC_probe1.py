"""Probe: fixed tanh-sinh k-grid vs adaptive mp.quad for a sample J(s,t)."""
import time
import mpmath as mp


def P(x, k):
    c = x * (1 - x)
    Del = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-Del) * (1 / Del ** 3 + 3 / Del ** 4 + 3 / Del ** 5)


def J_adaptive(s, t):
    b = s + t
    def f(k):
        j0 = mp.sin(k * b) / (k * b) if k * b > mp.mpf('1e-30') else mp.mpf(1)
        return j0 * P(s, k) * P(t, k)
    return mp.quad(f, [0, 1, 3, 8, 20, mp.inf])


def J_fixed(s, t, nodes):
    b = s + t
    tot = mp.mpf(0)
    for k, w in nodes:
        kb = k * b
        j0 = mp.sin(kb) / kb if kb > mp.mpf('1e-40') else mp.mpf(1) - kb * kb / 6
        tot += w * j0 * P(s, k) * P(t, k)
    return tot


def main():
    dps = 30
    mp.mp.dps = dps
    s, t = mp.mpf('0.3'), mp.mpf('0.5')

    t0 = time.time()
    jad = J_adaptive(s, t)
    print(f"adaptive: {mp.nstr(jad, dps-3)}  ({time.time()-t0:.2f}s)")

    rule = mp.calculus.quadrature.TanhSinh(mp.mp)
    for degree in [5, 6, 7, 8]:
        nodes = rule.get_nodes(0, mp.inf, degree, mp.mp.prec)
        t0 = time.time()
        jf = J_fixed(s, t, nodes)
        dt = time.time() - t0
        print(f"degree={degree} npts={len(nodes):5d}  {mp.nstr(jf, dps-3)}  diff={mp.nstr(abs(jf-jad),3)}  ({dt:.3f}s)")


if __name__ == '__main__':
    main()
