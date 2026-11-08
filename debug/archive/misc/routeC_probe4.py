"""Probe 4: full 2D T2 via separable sin/cos trick, comparing plain-GL-on-[0,1]
vs sin^2-substituted GL for the outer (s,t) integral, with a Mobius-mapped k grid.
"""
import time
import mpmath as mp


def P(x, k):
    c = x * (1 - x)
    Del = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-Del) * (1 / Del ** 3 + 3 / Del ** 4 + 3 / Del ** 5)


def std_gl_nodes(degree, prec):
    rule = mp.calculus.quadrature.GaussLegendre(mp.mp)
    return rule.calc_nodes(degree, prec)


def k_grid_mobius(std_nodes, K0):
    nodes = []
    for x, w in std_nodes:
        onepx = 1 + x
        k = K0 * (1 - x) / onepx
        jac = 2 * K0 / (onepx * onepx)
        nodes.append((k, w * jac))
    return nodes


def outer_grid_plain(degree):
    """Plain GL on [0,1], no substitution."""
    std_nodes = std_gl_nodes(degree, mp.mp.prec)
    s_list = []
    w_list = []
    for x, w in std_nodes:
        s_list.append((x + 1) / 2)
        w_list.append(w / 2)
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


def compute_T2(s_list, w_list, knodes):
    N = len(s_list)
    M = len(knodes)
    ks = [kn[0] for kn in knodes]
    wp = [kn[1] / kn[0] for kn in knodes]  # w'_m = w_m / k_m

    # Build U,V matrices: U[i][m]=sin(k_m s_i) P(s_i,k_m); V[i][m]=cos(k_m s_i) P(s_i,k_m)
    U = []
    V = []
    for si in s_list:
        Ui = []
        Vi = []
        for k in ks:
            Pv = P(si, k)
            Ui.append(mp.sin(k * si) * Pv)
            Vi.append(mp.cos(k * si) * Pv)
        U.append(Ui)
        V.append(Vi)

    # A[i][j] = sum_m wp_m U[i][m] V[j][m]
    total = mp.mpf(0)
    for i in range(N):
        Ui = U[i]
        Vi = V[i]
        wi = w_list[i]
        # diagonal j=i: Asym[i][i] = 2*A[i][i]; J = Asym/(2 s_i)
        Aii = mp.fsum(wp[m] * Ui[m] * Vi[m] for m in range(M))
        Jii = (2 * Aii) / (2 * s_list[i])
        total += wi * wi * Jii
        for j in range(i):
            Uj = U[j]
            Vj = V[j]
            Aij = mp.fsum(wp[m] * Ui[m] * Vj[m] for m in range(M))
            Aji = mp.fsum(wp[m] * Uj[m] * Vi[m] for m in range(M))
            Jij = (Aij + Aji) / (s_list[i] + s_list[j])
            total += 2 * wi * w_list[j] * Jij
    return (8 / mp.pi) * total


def main():
    dps = 30
    mp.mp.dps = dps

    # k-grid
    kdeg = 9  # n = 3*2^8 = 768
    std_k = std_gl_nodes(kdeg, mp.mp.prec)
    print(f"k std nodes: {len(std_k)}")

    for K0 in [80]:
        knodes = k_grid_mobius(std_k, mp.mpf(K0))
        kmax = max(kn[0] for kn in knodes)
        print(f"K0={K0} kmax={mp.nstr(kmax,8)}")

        for label, outer_fn in [("plain", outer_grid_plain), ("sin2", outer_grid_sin2)]:
            for deg in [5, 6]:
                s_list, w_list = outer_fn(deg)
                N = len(s_list)
                smin = min(s_list)
                t0 = time.time()
                v = compute_T2(s_list, w_list, knodes)
                dt = time.time() - t0
                print(f"{label:5s} deg={deg} N={N:4d} smin={mp.nstr(smin,4)}  T2={mp.nstr(v, dps-3)}  ({dt:.2f}s)")


if __name__ == '__main__':
    main()
