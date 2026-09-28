r"""PHASE C -- the full M=16 pi marriage energy E_R12 (deterministic, RI-free, explicit-r12 CI).

Standalone driver: launch with NO further edits as
    python debug/lih_pi_m16_energy.py
It builds the pi-containing M=16 natural-orbital reference (FCI cached from Phase A), runs the
seven m-extended enumerators (debug/lih_pi_m16.all_pieces_m, EXACT prolate-Neumann path) with the
general-m V_ee reducer, assembles the 2x2 -> E_R12, and writes the pieces + the linear-Jastrow
mixing coefficient c (for the G2-pi VMC cross-check, debug/lih_pi_m16_vmc.py).

  ** g_Vee (the <F^2 V_ee> 216-triple sum with the general-m triangle) is the wall-clock wall. **
  Use --trunc to trade accuracy for wall time (larger cut -> fewer dets -> fewer active pairs).

Options:
  --geminal linexp|exp|both   (default linexp; linexp is the cusp-correct geminal)
  --trunc  FLOAT              det truncation |c| cut (default 1e-3; try 2e-3 / 3e-3 if too slow)
Writes debug/data/lih_pi_m16_energy.npz + debug/data/lih_pi_m16_energy.log
"""
from __future__ import annotations

import os
import sys
import time
import argparse

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, "data")

import geovac.lih_r12ci.kernels as KG
import lih_pi_m16 as M
import lih_marriage_phase2 as P2
OUT = os.path.join(DATA, "lih_pi_m16_energy.npz")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--geminal', default='linexp', choices=['linexp', 'exp', 'both'])
    ap.add_argument('--trunc', type=float, default=1e-3)
    args = ap.parse_args()
    gems = ['linexp', 'exp'] if args.geminal == 'both' else [args.geminal]

    T0 = time.time()
    np.set_printoptions(linewidth=140, precision=6, suppress=True)
    KG.USE_EXACT_NEUMANN = True
    print("=" * 96)
    print(f"PHASE C -- M=16 pi marriage energy  (geminals={gems}, trunc={args.trunc:.0e}, EXACT path)")
    print("=" * 96, flush=True)

    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI)
    print(f"  M={ref['M']} Mk={ref['Mk']}  m_prim={ref['m_prim'].tolist()}", flush=True)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi',
                                        trunc=args.trunc)
    npi_pairs = int(np.sum(refm.mA != 0))
    print(f"  reference: E_trunc={info['E_trunc']:.6f}  dets kept={len(info['kept'])}  "
          f"active ordered pairs={refm.nact} ({npi_pairs} carry m!=0)  T nnz={len(refm.T_val)}",
          flush=True)
    print(f"  anchors: beat M=16 ceiling -8.032, land near VMC-Pade -8.047+-0.001, stay > -8.0705\n",
          flush=True)

    out = dict(Mk=info['Mk'], R=info['R'], Vnn=info['Vnn'], E_trunc=info['E_trunc'],
               E_full=info['E_full'], trunc=args.trunc, nact=refm.nact, npi_pairs=npi_pairs,
               T_no=info['T_no'], ct=info['ct'], m_no=info['m_no'], GAM=M.GAM,
               kept=info['kept'])
    for gname in gems:
        print(f"\n--- geminal = {gname} ---", flush=True)
        gem = M.GeminalM(gname)
        p = M.all_pieces_m(gem, refm, T0, force_general_vee=True)
        E_R12, c, h = M.assemble_c(p, p['E0_grid'])
        corr = (E_R12 - info['E_trunc']) * 1e3
        variational = E_R12 > -8.0705
        near047 = abs(E_R12 - (-8.047))
        print(f"\n  >>> E_R12 ({gname}) = {E_R12:.6f} Ha   (corr {corr:+.2f} mHa vs E_trunc, "
              f"c = {c:+.5f})", flush=True)
        print(f"      variational (> -8.0705): {variational}   |E_R12 - VMC-Pade(-8.047)| = "
              f"{near047*1e3:.1f} mHa", flush=True)
        for k, v in p.items():
            if isinstance(v, (int, float)):
                out[f"{gname}_{k}"] = v
        out[f"{gname}_E_R12"] = E_R12; out[f"{gname}_c"] = c; out[f"{gname}_h"] = h
        out[f"{gname}_Fbar"] = p['Fbar']
    os.makedirs(DATA, exist_ok=True)
    np.savez(OUT, **out)
    print(f"\n  saved {OUT}")
    print(f"  wall time {time.time() - T0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
