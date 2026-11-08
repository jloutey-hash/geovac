"""Probe 11: tanh-sinh (double-exponential) quadrature for the OUTER (s,t)
integral instead of GL, on the finite domain [0,1] directly (no sin^2 needed --
TS's own doubly-exponential endpoint clustering should handle the accumulating
complex-branch-point pathology near s=0,1 better than polynomial GL)."""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe5 import P
from routeC_probe8 import k_grid_sinh_paneled


def ts_outer_nodes(degree, prec):
    rule = mp.calculus.quadrature.TanhSinh(mp.mp)
    return rule.get_nodes(0, 1, degree, prec)


def compute_T2_ts(s_list, w_list, knodes):
    N = len(s_list)
    M = len(knodes)
    ks = [kn[0] for kn in knodes]
    wp = [kn[1] / kn[0] for kn in knodes]

    U = [[mp.mpf(0)] * M for _ in range(N)]
    V = [[mp.mpf(0)] * M for _ in range(N)]
    for i, si in enumerate(s_list):
        Ui, Vi = U[i], V[i]
        for m in range(M):
            k = ks[m]
            Pv = P(si, k)
            Ui[m] = mp.sin(k * si) * Pv
            Vi[m] = mp.cos(k * si) * Pv

    total = mp.mpf(0)
    for i in range(N):
        Ui, Vi, wi = U[i], V[i], w_list[i]
        Aii = mp.mpf(0)
        for m in range(M):
            Aii += wp[m] * Ui[m] * Vi[m]
        total += wi * wi * (Aii / s_list[i])
        for j in range(i):
            Uj, Vj, wj = U[j], V[j], w_list[j]
            Aij = mp.mpf(0)
            Aji = mp.mpf(0)
            for m in range(M):
                Aij += wp[m] * Ui[m] * Vj[m]
                Aji += wp[m] * Uj[m] * Vi[m]
            total += 2 * wi * wj * (Aij + Aji) / (s_list[i] + s_list[j])
    return (8 / mp.pi) * total


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    mp.mp.dps = dps
    pdeg = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    tsdegs = [int(a) for a in sys.argv[3:]] or [4, 5, 6, 7]

    knodes = k_grid_sinh_paneled(mp.mpf(14), mp.mpf(2), pdeg)
    print(f"dps={dps} k: pdeg={pdeg} M={len(knodes)}", flush=True)

    prev = None
    for tsdeg in tsdegs:
        nodes = ts_outer_nodes(tsdeg, mp.mp.prec)
        s_list = [n[0] for n in nodes]
        w_list = [n[1] for n in nodes]
        N = len(s_list)
        t0 = time.time()
        v = compute_T2_ts(s_list, w_list, knodes)
        dt = time.time() - t0
        d = "" if prev is None else mp.nstr(abs(v - prev), 4)
        print(f"  tsdeg={tsdeg} N={N:5d}  T2={mp.nstr(v, dps-3)}  diff_prev={d}  ({dt:.2f}s)", flush=True)
        prev = v


if __name__ == '__main__':
    main()
