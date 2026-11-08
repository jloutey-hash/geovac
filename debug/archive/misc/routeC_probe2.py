"""Probe 2: rational (Mobius) map k = K0*(1-x)/(1+x) + plain Gauss-Legendre for the
semi-infinite k-integral. Compare against adaptive mp.quad at a sample (s,t)."""
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


def gl_std(N):
    from mpmath import legendre
    roots = []
    for k in range(1, N + 1):
        x = mp.cos(mp.pi * (k - mp.mpf('0.25')) / (N + mp.mpf('0.5')))
        for _ in range(100):
            f = legendre(N, x)
            fp = N * (x * legendre(N, x) - legendre(N - 1, x)) / (x * x - 1)
            dx = f / fp
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        roots.append(x)
    ws = []
    for x in roots:
        fp = N * (x * legendre(N, x) - legendre(N - 1, x)) / (x * x - 1)
        ws.append(2 / ((1 - x * x) * fp * fp))
    return roots, ws


def k_grid_mobius(M, K0):
    """k = K0*(1-x)/(1+x), x in (-1,1); dk = K0 * (-2)/(1+x)^2 dx (take abs for weight)."""
    xs, ws = gl_std(M)
    nodes = []
    for x, w in zip(xs, ws):
        onepx = 1 + x
        k = K0 * (1 - x) / onepx
        jac = 2 * K0 / (onepx * onepx)
        nodes.append((k, w * jac))
    return nodes


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

    for K0 in [1, 2, 3]:
        for M in [40, 80, 160, 240]:
            nodes = k_grid_mobius(M, mp.mpf(K0))
            t0 = time.time()
            jf = J_fixed(s, t, nodes)
            dt = time.time() - t0
            print(f"K0={K0} M={M:4d}  {mp.nstr(jf, dps-3)}  diff={mp.nstr(abs(jf-jad),3)}  ({dt:.3f}s)")


if __name__ == '__main__':
    main()
