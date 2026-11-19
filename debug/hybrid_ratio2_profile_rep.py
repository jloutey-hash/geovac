"""Profile the integral and the exact collection of one v2 piece (same steps as hybrid_ratio2_job.piece).
Pickles the integral expression to the scratch path given, so collector variants can be timed on it.

Usage: python debug/hybrid_ratio2_profile_rep.py <a,b,c,d> <region> <ti> <lt|gt> <chunk> <nchunks> <pickle_out>
"""
from __future__ import annotations

import cProfile
import pickle
import pstats
import sys
import time
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "debug"))
import hybrid_ratio_engine as E  # noqa: E402
import hybrid_ratio_rep as HR  # noqa: E402

T = E.T
case = tuple(sys.argv[1].split(","))
region, ti, branch, chunk, nch = sys.argv[2], int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
a, b, c, d = (E.ORB[x] for x in case)
ZA, ZB, *_ = E.region_Z(region, a, b, c, d)
coeff, L, Lp, Mp, ld, md = list(T.hybrid_terms(ZA, a, b, c, ZB, d))[ti]
x_s = T.x_s
w_s = sp.Symbol("w", positive=True)
for LL, _M, _g, rad, bb in T.multipole_decomposition(ZA, *a, ZA, *b):
    if LL == L:
        rad1, b1 = rad, bb
        break
weight = (4 * sp.pi / (2 * L + 1) * sum(cc * x_s ** (k + 2) for k, cc in rad1.items()) * sp.exp(-b1 * x_s))
t0 = time.time()
e_br = T.shell_inner_double(ZA, a, b, c, ZB, d, L, Lp, Mp, ld, md, x_s, E.R, branch=branch)
ea = sp.Add.make_args(e_br)
e_br = sp.Add(*[t for i, t in enumerate(ea) if i % nch == chunk])
ex = E.expand_fixed(weight * e_br)
if branch == "gt":
    ex = E.expand_fixed(ex.subs(x_s, E.R + w_s))
    var, X = w_s, sp.oo
else:
    var, X = x_s, E.R
print("pre", round(time.time() - t0, 1), flush=True)
t0 = time.time()
pr = cProfile.Profile()
pr.enable()
p = T._integrate_with_e1(ex, var, X)
pr.disable()
print("integral", round(time.time() - t0, 1), len(sp.Add.make_args(p)), flush=True)
pstats.Stats(pr).sort_stats("cumtime").print_stats(25)
full = 2 * sp.pi / E.R * coeff * p
Path(sys.argv[7]).write_bytes(pickle.dumps(full))
t0 = time.time()
pr = cProfile.Profile()
pr.enable()
r = HR.rep(full)
pr.disable()
print("rep", round(time.time() - t0, 1), len(r), flush=True)
pstats.Stats(pr).sort_stats("cumtime").print_stats(30)
