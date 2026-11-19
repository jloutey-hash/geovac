"""Profile the steps of one piece job (inner double, weight expand, x -> R + w substitution, integral, collect).

Usage: python debug/hybrid_ratio2_profile.py <a,b,c,d> <region> <ti> <lt|gt>
"""
from __future__ import annotations

import cProfile
import pstats
import sys
import time
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "debug"))
import hybrid_ratio_engine as E  # noqa: E402

T = E.T
case = tuple(sys.argv[1].split(","))
region, ti, branch = sys.argv[2], int(sys.argv[3]), sys.argv[4]
a, b, c, d = (E.ORB[x] for x in case)
ZA, ZB, *_ = E.region_Z(region, a, b, c, d)
coeff, L, Lp, Mp, ld, md = list(T.hybrid_terms(ZA, a, b, c, ZB, d))[ti]
x_s = T.x_s
for LL, _M, _g, rad, bb in T.multipole_decomposition(ZA, *a, ZA, *b):
    if LL == L:
        rad1, b1 = rad, bb
        break
weight = (4 * sp.pi / (2 * L + 1) * sum(cc * x_s ** (k + 2) for k, cc in rad1.items()) * sp.exp(-b1 * x_s))
t0 = time.time()
pr = cProfile.Profile()
pr.enable()
e_br = T.shell_inner_double(ZA, a, b, c, ZB, d, L, Lp, Mp, ld, md, x_s, E.R, branch=branch)
pr.disable()
print("inner", round(time.time() - t0, 1), len(sp.Add.make_args(e_br)), flush=True)
pstats.Stats(pr).sort_stats("tottime").print_stats(12)
t0 = time.time()
pr = cProfile.Profile()
pr.enable()
ex = E.expand_fixed(weight * e_br)
pr.disable()
print("weight expand", round(time.time() - t0, 1), len(sp.Add.make_args(ex)), flush=True)
pstats.Stats(pr).sort_stats("tottime").print_stats(12)
if branch == "gt":
    w_s = sp.Symbol("w", positive=True)
    t0 = time.time()
    pr = cProfile.Profile()
    pr.enable()
    ex = E.expand_fixed(ex.subs(x_s, E.R + w_s))
    pr.disable()
    print("subs expand", round(time.time() - t0, 1), len(sp.Add.make_args(ex)), flush=True)
    pstats.Stats(pr).sort_stats("tottime").print_stats(15)
