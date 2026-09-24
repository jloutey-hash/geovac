r"""Phase-2 74-det PREVIEW (UNGATED) + sparse-triangle regression on the 2-MO g_Vee.

(1) verify the sparse-triangle vee_reduce reproduces the 2-MO exact-path g_Vee (0.535896 exp /
    0.463641 linexp) -- a regression on the g_Vee generalisation before the expensive 74-det run;
(2) feed the real 74-det tensor T (Phase 0b) through the validated enumerators on the EXACT kernel
    path and report E_R12 as a PREVIEW awaiting the Phase-3 G2 VMC cross-check.  Expected window
    -8.045..-8.055 (near VMC -8.047; the pairwise ansatz cannot reach -8.062).

Run:  python debug/lih_marriage_phase2_preview.py > debug/data/lih_marriage_phase2_preview.log 2>&1; echo $? > ..._preview.exit
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import debug.lih_marriage_phase2 as M

G_VEE_2MO_EXACT = {'exp': 0.535896, 'linexp': 0.463641}   # from the fixed full G1 run (live oracle)


def main():
    T0 = time.time()
    M.KG.USE_EXACT_NEUMANN = True
    M.GV._NEU_CACHE.clear(); M._MODE_CACHE.clear()
    print("=" * 96)
    print("PHASE-2 PREVIEW (exact kernel path)")
    print("=" * 96, flush=True)

    ref2 = M.build_ref_2mo()
    print("\n(1) sparse-triangle regression on the 2-MO g_Vee (exact path):", flush=True)
    for gname in ('exp', 'linexp'):
        gem = M.Geminal(gname)
        ref2._cache.clear()
        h_Vee, g_Vee, Vee, FVee, F2Vee = M.hg_Vee(gem, ref2, M.fbar(gem, ref2))
        tgt = G_VEE_2MO_EXACT[gname]
        print(f"    {gname:7s}: g_Vee(sparse) = {g_Vee:.6f}  vs G1(dense) {tgt:.6f}  "
              f"d = {(g_Vee - tgt) * 1e3:+.4f} mHa   {'OK' if abs(g_Vee - tgt) < 1e-4 else 'MISMATCH'}"
              f"   [{time.time() - T0:.0f}s]", flush=True)

    print("\n(2) 74-det tensor through the enumerators (exact path, UNGATED):", flush=True)
    M.GV._NEU_CACHE.clear(); M._MODE_CACHE.clear()
    refA = M.build_ref_artifact()
    save = {}
    for gname in ('exp', 'linexp'):
        print(f"\n  --- PREVIEW geminal = {gname} ---", flush=True)
        gem = M.Geminal(gname)
        p = M.all_pieces(gem, refA, T0)
        E_R12, _ = M.assemble(p, p['E0_grid'])
        inwin = -8.055 <= E_R12 <= -8.045
        print(f"  >>> PREVIEW E_R12 ({gname}) = {E_R12:.6f} Ha  E0_grid={p['E0_grid']:.6f}  "
              f"h=({p['h_T']:+.4f},{p['h_Vne']:+.4f},{p['h_Vee']:+.4f})  "
              f"g=({p['g_T']:+.4f},{p['g_Vne']:+.4f},{p['g_Vee']:+.4f})  sigma2={p['sigma2']:.5f}",
              flush=True)
        print(f"      window -8.045..-8.055: {'IN' if inwin else 'OUTSIDE'}   [{time.time() - T0:.0f}s]",
              flush=True)
        save[f"preview_{gname}_E_R12"] = E_R12
        for k, v in p.items():
            if isinstance(v, (int, float)):
                save[f"preview_{gname}_{k}"] = v
    np.savez(os.path.join(HERE, "data", "lih_marriage_phase2_preview.npz"), **save)
    print(f"\nsaved preview npz;  wall {time.time() - T0:.0f} s")


if __name__ == "__main__":
    main()
