"""LT instrumented run: wall-clock share of the three-body triangle class in the LiH pi-marriage energy.

Mirrors tests/data/lih_golden/capture.py::_pi_marriage("exp", trunc, 1.5) (the golden items
pi_marriage_exp1.5_1e-2 and pi_marriage_exp1.5_1e-3), and wraps timers around existing functions by
monkeypatching module attributes.  geovac/ is not edited.  Wrappers only call the original function
and read perf_counter, so the numerics and the RNG streams are untouched.

NOT RUN in phase 1 (pre-registration record docs/predictions/2026-10-26_lih_triangle_share.md).

Usage (fresh process, single BLAS thread, like the golden gate):
    set OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
    python debug/lih_triangle_share.py --trunc 1e-2
    python debug/lih_triangle_share.py --trunc 1e-3        # about 3 h, about 8.4 GB
    python debug/lih_triangle_share.py --trunc 1e-2 --cprofile   # cross-check, slower

Output: debug/data/lih_triangle_share_<trunc>.json and a printed table.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEBUG = os.path.join(ROOT, "debug")
for _p in (ROOT, DEBUG):
    if _p not in sys.path:
        sys.path.insert(0, _p)

T = defaultdict(float)      # name -> accumulated seconds
N = defaultdict(int)        # name -> call count
EXTRA = defaultdict(int)    # other counters


def timed(mod, name, label=None):
    """Replace mod.name by a timing wrapper (nested labels are NOT subtracted; see report)."""
    orig = getattr(mod, name)
    lab = label or name

    def wrapper(*a, **k):
        t0 = time.perf_counter()
        try:
            return orig(*a, **k)
        finally:
            T[lab] += time.perf_counter() - t0
            N[lab] += 1
    wrapper.__name__ = name
    setattr(mod, name, wrapper)
    return orig


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trunc", type=float, default=1e-2)
    ap.add_argument("--gam", type=float, default=1.5)
    ap.add_argument("--cprofile", action="store_true")
    args = ap.parse_args()
    sys.argv = [sys.argv[0]]

    from geovac.lih_r12ci import configure_geminal
    configure_geminal(args.gam)                       # before the consumer imports (golden order)
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.pi_m16 as M
    import geovac.lih_r12ci.gvee_sparse as S
    import geovac.lih_r12ci.triangle as TRI
    from geovac.lih_r12ci import gVee as GV
    import geovac.lih_r12ci.enumerators as P2
    assert P2.GAM == M.GAM == GV.GAM == KG.GAM == args.gam

    # ---- timers: whole-piece level ----
    timed(M, "build_primitives"); timed(M, "build_ref_from_prims")
    timed(M, "fbar"); timed(M, "sigma2")
    timed(M, "hg_Vne"); timed(M, "hg_T"); timed(M, "assemble_c")
    timed(S, "hg_Vee_sparse")
    # g_Vee split: non-triangle (every vee_reduce_sparse call) vs triangle phase
    timed(S, "vee_reduce_sparse", "gVee_nontri_vee_reduce_sparse")
    timed(S, "eval_triangles", "TRIANGLE_eval_triangles")
    # (tri_entries is timed below, once, around a counting wrapper)
    # inside the triangle phase
    timed(M, "coul_mode_stack_any", "tri_coulomb_mode_solves")
    # K_mu builders; count cache misses of the symmetric-kernel path separately
    orig_kmu = M.Kmu_modes

    def kmu_wrapper(Fm, mu, A, key, *a, **k):
        miss = (key, A, k.get("rank", a[0] if a else 30)) not in TRI._MODE_CACHE
        t0 = time.perf_counter()
        try:
            return orig_kmu(Fm, mu, A, key, *a, **k)
        finally:
            T["tri_Kmu_modes(eigh)"] += time.perf_counter() - t0
            N["tri_Kmu_modes(eigh)_calls"] += 1
            if miss:
                EXTRA["Kmu_eigh_builds"] += 1
    M.Kmu_modes = kmu_wrapper
    timed(M, "_rand_svd_middle", "tri_rand_svd_middle")

    # ---- count distinct middle densities used (tri cm) ----
    orig_tri_entries = S.tri_entries

    def tri_entries_counting(ctx, cm, il, ir):
        EXTRA["tri_cm_seen"] += 1
        EXTRA["tri_entries_requested"] += len(il)
        return orig_tri_entries(ctx, cm, il, ir)
    S.tri_entries = tri_entries_counting
    timed(S, "tri_entries", "TRIANGLE_tri_entries")

    prof = None
    if args.cprofile:
        import cProfile
        prof = cProfile.Profile()

    t_start = time.perf_counter()
    KG.USE_EXACT_NEUMANN = True
    if prof:
        prof.enable()
    T0 = time.time()
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag="pi",
                                        trunc=args.trunc, verbose=False)
    gem = M.GeminalM("exp")
    p = S.all_pieces_sparse(gem, refm, T0, mem_cap_gb=6.0, verbose=True)
    E_R12, c, h = M.assemble_c(p, p["E0_grid"])
    if prof:
        prof.disable()
    wall = time.perf_counter() - t_start

    # ---- middle densities and mode kernels (what the record's LT2 predicts) ----
    import numpy as np
    mA = np.asarray(refm.mA, int)
    mm = {int(v): int((mA == v).sum()) for v in np.unique(mA)}
    out = dict(trunc=args.trunc, gam=args.gam, nact=int(refm.nact), E_R12=float(E_R12), c=float(c),
               h=float(h), E_trunc=float(info["E_trunc"]), wall_s=wall,
               mA_counts=mm, timers_s=dict(T), call_counts=dict(N), counters=dict(EXTRA),
               pieces={k: float(v) for k, v in p.items() if isinstance(v, (int, float, np.floating))})
    path = os.path.join(DEBUG, "data", f"lih_triangle_share_{args.trunc:g}.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nE_R12 = {E_R12!r}  nact={refm.nact}  wall={wall:.1f} s  mA counts {mm}")
    for k in sorted(T):
        print(f"  {k:42s} {T[k]:10.1f} s  {100*T[k]/wall:5.1f}%  calls {N.get(k, '')}")
    print("  counters:", dict(EXTRA))
    if prof:
        import pstats
        pstats.Stats(prof).sort_stats("cumulative").print_stats(40)
    print("saved", path)


if __name__ == "__main__":
    main()
