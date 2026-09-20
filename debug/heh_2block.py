"""(a) Close the HeH+ single-alpha diagnosis: does a TWO-BLOCK (two-exponent)
basis drop the heteronuclear error sharply?  p=0 (no r12) isolates the radial-
scale question.  If two-block << single-alpha, the bottleneck is the single-alpha
basis (He Z=2 contracted vs H Z=1 diffuse), NOT the r12 machinery.
Run from root:  python debug/heh_2block.py
"""
import os
import sys
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geovac import prolate_recondition as pr
import prolate_r12_mpf as m
from r12ci_first_energy import solve_canonical

ZA, ZB = 2, 1
R = 1.4632
E_REF = -2.97869
NUC = ZA * ZB / R


def index(j_max, l_max):
    return [(j, l, k, mm) for j in range(j_max + 1) for k in range(j_max + 1)
            for l in range(l_max + 1) for mm in range(l_max + 1)]


def funcs(idx, alpha):
    return [pr.ProductFn(j, l, k, mm, 0, alpha) for (j, l, k, mm) in idx]


def single_alpha(j_max, l_max, alpha, use_mpf=True):
    """p=0 HeH+ at one exponent (baseline)."""
    basis = [(f, 0) for f in funcs(index(j_max, l_max), alpha)]
    if use_mpf:
        S, H = m.assemble_hetero(basis, R, alpha, ZA, ZB, mpf_out=True)
        e = m.solve_canonical_mpf(S, H)[0]
    else:
        S, H = m.assemble_hetero(basis, R, alpha, ZA, ZB)
        e = solve_canonical(S, H)[0]
    return e + NUC


def two_block(j_max, l_max, aa, ab, use_mpf=True):
    idx = index(j_max, l_max)
    A_funcs, B_funcs = funcs(idx, aa), funcs(idx, ab)
    if use_mpf:
        S, H = m.assemble_hetero_2block(A_funcs, B_funcs, R, ZA, ZB, mpf_out=True)
        e = m.solve_canonical_mpf(S, H)[0]
    else:
        S, H = m.assemble_hetero_2block(A_funcs, B_funcs, R, ZA, ZB)
        e = solve_canonical(S, H)[0]
    return e + NUC


if __name__ == "__main__":
    jm, lm = 2, 2   # more converged than (2,1); confirm two-block still barely helps
    print(f"HeH+ p=0, (j={jm},l={lm}); E_ref={E_REF}; single-block n="
          f"{len(index(jm,lm))}, two-block n={2*len(index(jm,lm))}\n")

    print("single-alpha baseline:")
    best_s = None
    for a in [1.4, 1.7, 2.0]:
        e = single_alpha(jm, lm, a)
        print(f"  a={a:.2f}  E_tot={e:.6f}  err={(E_REF-e)*1000:+.3f} mHa", flush=True)
        if best_s is None or e < best_s[1]:
            best_s = (a, e)
    print(f"  BEST single: a={best_s[0]:.2f}  err={(E_REF-best_s[1])*1000:+.3f} mHa\n")

    print("two-block (He-contracted a_a, H-diffuse a_b):")
    best_t = None
    for (aa, ab) in [(2.6, 1.2), (3.0, 1.0), (2.2, 0.9)]:
        e = two_block(jm, lm, aa, ab)
        print(f"  (a_a={aa:.1f}, a_b={ab:.1f})  E_tot={e:.6f}  err={(E_REF-e)*1000:+.3f} mHa", flush=True)
        if best_t is None or e < best_t[1]:
            best_t = ((aa, ab), e)
    print(f"  BEST two-block: {best_t[0]}  err={(E_REF-best_t[1])*1000:+.3f} mHa\n")

    imp = (E_REF - best_s[1]) - (E_REF - best_t[1])
    print(f"VERDICT: two-block improves single-alpha by {imp*1000:+.3f} mHa "
          f"({(E_REF-best_s[1])*1000:+.3f} -> {(E_REF-best_t[1])*1000:+.3f})")
