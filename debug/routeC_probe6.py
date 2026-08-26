"""Probe 6: systematic (K0, kdeg, N) convergence sweep with explicit numeric diffs."""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe5 import (P, std_gl_nodes, k_grid_mobius, outer_grid_plain,
                            outer_grid_sin2, compute_T2)


def main():
    dps = 60
    mp.mp.dps = dps
    N_deg = 5  # N=48 sin2
    s_list, w_list = outer_grid_sin2(N_deg)
    print(f"N={len(s_list)} (sin2)")

    results = {}
    for kdeg in [7, 8, 9, 10]:
        std_k = std_gl_nodes(kdeg, mp.mp.prec)
        for K0 in [20, 40, 60, 100]:
            knodes = k_grid_mobius(std_k, mp.mpf(K0))
            t0 = time.time()
            v = compute_T2(s_list, w_list, knodes)
            dt = time.time() - t0
            results[(kdeg, K0)] = v
            print(f"  kdeg={kdeg} M={len(knodes):5d} K0={K0:4d}  T2={mp.nstr(v, dps-3)}  ({dt:.2f}s)", flush=True)

    # diffs vs the largest (kdeg,K0)
    best = results[(10, 100)]
    print("\nDiffs vs kdeg=10,K0=100:")
    for key, v in results.items():
        print(f"  {key}: {mp.nstr(abs(v-best), 4)}")


if __name__ == '__main__':
    main()
