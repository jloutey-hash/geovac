"""Probe 5: leaner, timed pipeline. Plain-loop accumulation (no fsum generator
overhead). Compare plain-GL vs sin2-GL outer substitution, converge in N and M."""
import time
import mpmath as mp


def P(x, k):
    c = x * (1 - x)
    Del = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-Del) * (1 / Del ** 3 + 3 / Del ** 4 + 3 / Del ** 5)


_gl_cache = {}


def std_gl_nodes(degree, prec):
    key = (degree, prec)
    if key in _gl_cache:
        return _gl_cache[key]
    rule = mp.calculus.quadrature.GaussLegendre(mp.mp)
    nodes = rule.calc_nodes(degree, prec)
    _gl_cache[key] = nodes
    return nodes


def k_grid_mobius(std_nodes, K0):
    nodes = []
    for x, w in std_nodes:
        onepx = 1 + x
        k = K0 * (1 - x) / onepx
        jac = 2 * K0 / (onepx * onepx)
        nodes.append((k, w * jac))
    return nodes


def outer_grid_plain(degree):
    std_nodes = std_gl_nodes(degree, mp.mp.prec)
    s_list = [(x + 1) / 2 for x, w in std_nodes]
    w_list = [w / 2 for x, w in std_nodes]
    return s_list, w_list


def outer_grid_sin2(degree):
    std_nodes = std_gl_nodes(degree, mp.mp.prec)
    half = mp.pi / 4
    s_list = []
    w_list = []
    for x, w in std_nodes:
        phi = half * (x + 1)
        s_list.append(mp.sin(phi) ** 2)
        w_list.append(half * w * mp.sin(2 * phi))
    return s_list, w_list


def compute_T2(s_list, w_list, knodes, verbose=False):
    N = len(s_list)
    M = len(knodes)
    ks = [kn[0] for kn in knodes]
    wp = [kn[1] / kn[0] for kn in knodes]

    t0 = time.time()
    U = [[mp.mpf(0)] * M for _ in range(N)]
    V = [[mp.mpf(0)] * M for _ in range(N)]
    for i, si in enumerate(s_list):
        Ui = U[i]
        Vi = V[i]
        for m in range(M):
            k = ks[m]
            Pv = P(si, k)
            Ui[m] = mp.sin(k * si) * Pv
            Vi[m] = mp.cos(k * si) * Pv
    t_build = time.time() - t0
    if verbose:
        print(f"    build U,V: {t_build:.3f}s")

    t0 = time.time()
    total = mp.mpf(0)
    for i in range(N):
        Ui, Vi, wi = U[i], V[i], w_list[i]
        Aii = mp.mpf(0)
        for m in range(M):
            Aii += wp[m] * Ui[m] * Vi[m]
        Jii = Aii / s_list[i]
        total += wi * wi * Jii
        for j in range(i):
            Uj, Vj, wj = U[j], V[j], w_list[j]
            Aij = mp.mpf(0)
            Aji = mp.mpf(0)
            for m in range(M):
                Aij += wp[m] * Ui[m] * Vj[m]
                Aji += wp[m] * Uj[m] * Vi[m]
            Jij = (Aij + Aji) / (s_list[i] + s_list[j])
            total += 2 * wi * wj * Jij
    t_sum = time.time() - t0
    if verbose:
        print(f"    bilinear sum: {t_sum:.3f}s")
    return (8 / mp.pi) * total


def main():
    import sys
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    mp.mp.dps = dps
    kdeg = int(sys.argv[2]) if len(sys.argv) > 2 else 7  # n=192
    K0 = mp.mpf(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf(60)
    degs = [int(a) for a in sys.argv[4:]] or [5, 6]

    std_k = std_gl_nodes(kdeg, mp.mp.prec)
    knodes = k_grid_mobius(std_k, K0)
    kmax = max(kn[0] for kn in knodes)
    print(f"dps={dps} K0={K0} M={len(knodes)} kmax={mp.nstr(kmax,6)}", flush=True)
    for label, outer_fn in [("sin2", outer_grid_sin2), ("plain", outer_grid_plain)]:
        for deg in degs:
            s_list, w_list = outer_fn(deg)
            N = len(s_list)
            t0 = time.time()
            v = compute_T2(s_list, w_list, knodes, verbose=False)
            dt = time.time() - t0
            print(f"  {label:5s} N={N:4d}  T2={mp.nstr(v, dps-3)}  ({dt:.2f}s)", flush=True)


if __name__ == '__main__':
    main()
