"""Is the ~0.05 mHa plateau an alpha artifact or the r12-power (p<=1) floor?
Variational alpha scan at a fixed (j,l) basis, mpf-orthogonalized solve.
"""
import os
import sys
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geovac import prolate_recondition as pr
import prolate_r12_mpf as m
from r12ci_first_energy import solve_canonical, NUC, E_EXACT, DE_EXACT

R = 1.4011


def build_basis_alpha(j_max, l_max, alpha, p_set=(0, 1)):
    rlist = [(j, k) for j in range(j_max + 1) for k in range(j_max + 1)]
    alist = [(l, mm) for l in range(l_max + 1) for mm in range(l_max + 1)
             if (l + mm) % 2 == 0]
    return [(pr.ProductFn(j, l, k, mm, 0, alpha), p)
            for p in p_set for (j, k) in rlist for (l, mm) in alist]


def scan(j_max, l_max, alphas):
    print(f"alpha scan at (j={j_max}, l={l_max}):", flush=True)
    best = None
    for a in alphas:
        basis = build_basis_alpha(j_max, l_max, a)
        t = time.time()
        S, H = m.assemble_mixed(basis, R, a)   # float64 (fast; enough to locate opt)
        bf = solve_canonical(S, H)
        et = bf[0] + NUC
        err = (E_EXACT - et) * 1000
        de = (-1.0 - et) / DE_EXACT * 100
        print(f"  alpha={a:.2f}  E_tot={et:.7f}  D_e%={de:.4f}  err={err:+.4f} mHa"
              f"  kept={bf[2]}/{len(basis)}  [{time.time()-t:.0f}s]", flush=True)
        if best is None or et < best[1]:
            best = (a, et, err)
    print(f"  BEST: alpha={best[0]:.2f}  err={best[2]:+.4f} mHa", flush=True)


if __name__ == "__main__":
    scan(3, 3, [1.0, 1.2, 1.4, 1.6, 1.8])
