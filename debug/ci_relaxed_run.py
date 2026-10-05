r"""Driver: analytic ansatz (A) -- CI coefficients relaxed alongside the e^{-1.5 r12} geminal on the
truncated M=16 pi reference (geovac/lih_r12ci/ci_relaxed.py).  Gate G1 (frozen restriction == the
banked 2x2) + the relaxed energy (gate G3 numbers).

Frozen 2x2 pieces (Fbar, sigma2, g, h, E0_grid) are taken from the golden oracle
tests/data/lih_golden/golden.json (bit-identical to a fresh all_pieces_sparse run; the 90-det g_Vee
is ~3 h / 8.4 GB), or recomputed with --recompute-pieces.

Usage:  python debug/ci_relaxed_run.py 1e-2 [--recompute-pieces] [--save debug/data/ci_relaxed_1e-2.npz]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

GOLDEN = os.path.join(ROOT, "tests", "data", "lih_golden", "golden.json")
ITEM = {1e-2: "pi_marriage_exp1.5_1e-2", 1e-3: "pi_marriage_exp1.5_1e-3"}


def golden_pieces(trunc: float) -> dict:
    g = json.load(open(GOLDEN, encoding="utf-8"))
    items = g.get("items", g)
    vals = items[ITEM[trunc]]["values"]
    pieces = {k: float(v["repr"]) for k, v in vals["pieces"]["value"].items()}
    return dict(pieces=pieces, E_R12=float(vals["E_R12"]["repr"]), c=float(vals["c"]["repr"]),
                h=float(vals["h"]["repr"]))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("trunc", type=float)
    ap.add_argument("--recompute-pieces", action="store_true")
    ap.add_argument("--save", default=None)
    a = ap.parse_args()
    from geovac.lih_r12ci import configure_geminal
    configure_geminal(1.5)                                    # before the consumer imports (golden order)
    import numpy as np
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.pi_m16 as M
    import geovac.lih_r12ci.ci_relaxed as CR
    T0 = time.time()
    KG.USE_EXACT_NEUMANN = True
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=a.trunc,
                                        verbose=False)
    gem = M.GeminalM('exp')
    gold = golden_pieces(a.trunc)
    if a.recompute_pieces:
        import geovac.lih_r12ci.gvee_sparse as GS
        pieces = GS.all_pieces_sparse(gem, refm, T0, mem_cap_gb=6.0, verbose=True)
    else:
        pieces = gold["pieces"]
    print(f"trunc={a.trunc} N_det={len(info['kept'])} E_trunc={info['E_trunc']:.9f}  [{time.time()-T0:.1f}s]",
          flush=True)
    r = CR.ci_relaxed_energy(gem, refm, info, pieces, verbose=True)
    print("\nG1 identity diagnostics (new-code contraction vs frozen 2x2 pieces):")
    for k, v in r["identity"].items():
        print(f"  {k:18s} {v:+.3e}")
    print(f"  E_frozen (new code)   = {r['E_frozen']:.15f}")
    print(f"  E_R12 golden (2x2)    = {gold['E_R12']:.15f}   diff {r['E_frozen'] - gold['E_R12']:+.3e}")
    print(f"\n  E0_grid (d.H_CI.d)    = {r['E0_grid']:.9f}")
    print(f"  CI-only relaxed (grid)= {r['E_CI_relaxed_only']:.9f}")
    print(f"  E(A) relaxed          = {r['E_A']:.9f}   E(A)-E(frozen) = {(r['E_A']-r['E_frozen'])*1e3:+.4f} mHa")
    print(f"  |<d_A|d_ref>| = {r['overlap']:.6f}   geminal c: frozen {r['c_frozen']:+.5f}  relaxed {r['c_geminal']:+.5f}")
    print(f"  grid S_CI - 1: max|.| = {np.max(np.abs(r['S_CI_grid'] - np.eye(r['N']))):.2e}")
    print(f"  H_CI asym max = {np.max(np.abs(r['H_CI'] - r['H_CI'].T)):.2e}")
    for I in range(min(r['N'], 12)):
        print(f"   I={I:3d} d={r['d_ref'][I]:+.5f} relaxed {r['d_relaxed'][I]:+.5f}  s={r['s'][I]:+.6f} "
              f"t={r['t'][I]:+.6f} (T {r['tT'][I]:+.6f} Vne {r['tVne'][I]:+.6f} Vee {r['tVee'][I]:+.6f})")
    print(f"wall {time.time()-T0:.1f}s  columns {r['wall_s']:.1f}s  peak RSS {M.__dict__.get('_x', 0) or 0}")
    try:
        import psutil
        mi = psutil.Process().memory_info()
        print(f"peak RSS {getattr(mi, 'peak_wset', mi.rss)/1e9:.2f} GB")
    except Exception:
        pass
    if a.save:
        np.savez(a.save, **{k: v for k, v in r.items() if k != 'identity'},
                 identity=json.dumps(r['identity']), golden_E_R12=gold['E_R12'])
    return 0


if __name__ == "__main__":
    sys.exit(main())
