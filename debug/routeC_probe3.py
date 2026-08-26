"""Probe 3: Mobius-map GL for k-integral at a SLOW-DECAY corner (small s), using
mpmath's fast built-in GaussLegendre.calc_nodes for the standard [-1,1] rule."""
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
    return mp.quad(f, [0, 1, 3, 8, 20, 200, mp.inf])


def std_gl_nodes(degree, prec):
    rule = mp.calculus.quadrature.GaussLegendre(mp.mp)
    return rule.calc_nodes(degree, prec)


def k_grid_mobius(std_nodes, K0):
    nodes = []
    kmax = mp.mpf(0)
    for x, w in std_nodes:
        onepx = 1 + x
        k = K0 * (1 - x) / onepx
        jac = 2 * K0 / (onepx * onepx)
        nodes.append((k, w * jac))
        if k > kmax:
            kmax = k
    return nodes, kmax


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
    s, t = mp.mpf('1e-4'), mp.mpf('0.5')
    print(f"s={s}, t={t}, c_s={s*(1-s)}, sqrt(c_s)={mp.sqrt(s*(1-s))}")

    t0 = time.time()
    jad = J_adaptive(s, t)
    print(f"adaptive: {mp.nstr(jad, dps-3)}  ({time.time()-t0:.2f}s)")

    for degree in [7, 8, 9, 10]:
        std_nodes = std_gl_nodes(degree, mp.mp.prec)
        n = len(std_nodes)
        for K0 in [10, 30, 60]:
            nodes, kmax = k_grid_mobius(std_nodes, mp.mpf(K0))
            t0 = time.time()
            jf = J_fixed(s, t, nodes)
            dt = time.time() - t0
            print(f"degree={degree} n={n:5d} K0={K0:3d} kmax={mp.nstr(kmax,6):>12s}  {mp.nstr(jf, dps-3)}  diff={mp.nstr(abs(jf-jad),3)}  ({dt:.3f}s)")


if __name__ == '__main__':
    main()
