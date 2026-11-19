"""Timing probe for the parametric engine (one shell_inner_double call)."""
import cProfile
import pstats
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sympy as sp  # noqa: E402

import hybrid_ratio_engine as E  # noqa: E402
import geovac.two_center_eri as T  # noqa: E402

case = tuple((sys.argv[1] if len(sys.argv) > 1 else "2p0,2p0,1s,1s").split(","))
a, b, c, d = (E.ORB[x] for x in case)
ZA, ZB, ap, ac, ad = E.region_Z("mid", a, b, c, d)
t0 = time.time()
terms = list(T.hybrid_terms(ZA, a, b, c, ZB, d))
print(terms, round(time.time() - t0, 1), flush=True)
co, L, Lp, Mp, ld, md = terms[-1]
pr = cProfile.Profile()
pr.enable()
e = T.shell_inner_double(ZA, a, b, c, ZB, d, L, Lp, Mp, ld, md, T.x_s, T.R_s, branch="lt")
pr.disable()
print("lt", round(time.time() - t0, 1), len(sp.Add.make_args(e)), flush=True)
pstats.Stats(pr).sort_stats("cumtime").print_stats(15)

if len(sys.argv) > 2:
    rad1, b1 = None, None
    for LL, _M, _g, rad, bb in T.multipole_decomposition(ZA, *a, ZA, *b):
        if LL == L:
            rad1, b1 = rad, bb
    weight = sum(cc * T.x_s ** (k + 2) for k, cc in rad1.items()) * sp.exp(-b1 * T.x_s)
    lt = sp.expand(weight * e)
    print("lt terms", len(sp.Add.make_args(lt)), flush=True)
    pr = cProfile.Profile()
    pr.enable()
    part = T._integrate_with_e1(lt, T.x_s, T.R_s)
    pr.disable()
    print("int", round(time.time() - t0, 1), flush=True)
    pstats.Stats(pr).sort_stats("tottime").print_stats(15)
