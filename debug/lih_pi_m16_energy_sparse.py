r"""PHASE C (O(nnz) g_Vee) -- the M=16 pi marriage energy E_R12 via debug/lih_gvee_sparse.

Copy of debug/lih_pi_m16_energy.py with g_Vee routed through the O(nnz) reducer
(lih_gvee_sparse.all_pieces_sparse; exact twin of the dense all_pieces_m, gates R1/R2).
Memory is bounded by --mem-cap-gb (chunked field / potential / triangle solves).

    python debug/lih_pi_m16_energy_sparse.py --trunc 1e-3 --geminal linexp --mem-cap-gb 6

Writes debug/data/lih_pi_m16_energy_sparse_<trunc>.npz (per-trunc, never overwritten -- a
timestamp suffix is added if the file exists) + prints E_R12, c, Fbar, E_trunc, nact, peak RSS, wall.
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
import lih_gvee_sparse as S


def _outpath(trunc: float) -> str:
    base = os.path.join(DATA, f"lih_pi_m16_energy_sparse_{trunc:.0e}")
    p = base + ".npz"
    if os.path.exists(p):
        p = base + time.strftime("_%Y%m%d_%H%M%S") + ".npz"
    return p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--geminal', default='linexp', choices=['linexp', 'exp', 'both'])
    ap.add_argument('--trunc', type=float, default=1e-3)
    ap.add_argument('--mem-cap-gb', type=float, default=6.0)
    args = ap.parse_args()
    gems = ['linexp', 'exp'] if args.geminal == 'both' else [args.geminal]
    out_path = _outpath(args.trunc)

    T0 = time.time()
    np.set_printoptions(linewidth=140, precision=6, suppress=True)
    KG.USE_EXACT_NEUMANN = True
    print("=" * 96)
    print(f"PHASE C (O(nnz) g_Vee) -- M=16 pi marriage energy  (geminals={gems}, trunc={args.trunc:.0e}, "
          f"mem cap {args.mem_cap_gb} GB, EXACT path)")
    print("=" * 96, flush=True)

    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI)
    refm, info = M.build_ref_from_prims(ref, "M=16 pi natural-orbital", cfg_tag='pi', trunc=args.trunc)
    npi_pairs = int(np.sum(refm.mA != 0))
    print(f"  reference: E_trunc={info['E_trunc']:.6f}  dets kept={len(info['kept'])}  "
          f"active ordered pairs={refm.nact} ({npi_pairs} carry m!=0)  T nnz={len(refm.T_val)}", flush=True)
    print("  anchors: beat M=16 ceiling -8.032, land near VMC-Pade -8.047+-0.001, stay > -8.0705\n", flush=True)

    out = dict(Mk=info['Mk'], R=info['R'], Vnn=info['Vnn'], E_trunc=info['E_trunc'],
               E_full=info['E_full'], trunc=args.trunc, nact=refm.nact, npi_pairs=npi_pairs,
               nnz=len(refm.T_val), T_no=info['T_no'], ct=info['ct'], m_no=info['m_no'], GAM=M.GAM,
               kept=info['kept'], mem_cap_gb=args.mem_cap_gb)
    for gname in gems:
        print(f"\n--- geminal = {gname} ---", flush=True)
        gem = M.GeminalM(gname)
        p = S.all_pieces_sparse(gem, refm, T0, mem_cap_gb=args.mem_cap_gb)
        E_R12, c, h = M.assemble_c(p, p['E0_grid'])
        corr = (E_R12 - info['E_trunc']) * 1e3
        print(f"\n  >>> E_R12 ({gname}) = {E_R12:.6f} Ha   (corr {corr:+.2f} mHa vs E_trunc "
              f"{info['E_trunc']:.6f}, c = {c:+.5f}, Fbar = {p['Fbar']:.6f})", flush=True)
        print(f"      variational (> -8.0705): {E_R12 > -8.0705}   |E_R12 - VMC-Pade(-8.047)| = "
              f"{abs(E_R12 + 8.047)*1e3:.1f} mHa   nact={refm.nact}", flush=True)
        for k, v in p.items():
            if isinstance(v, (int, float)):
                out[f"{gname}_{k}"] = v
        out[f"{gname}_E_R12"] = E_R12; out[f"{gname}_c"] = c; out[f"{gname}_h"] = h
        out[f"{gname}_Fbar"] = p['Fbar']
        np.savez(out_path, **out)                     # checkpoint after each geminal
    wall = time.time() - T0
    out['wall_s'] = wall; out['peak_rss_gb'] = S.peak_rss_gb()
    os.makedirs(DATA, exist_ok=True)
    np.savez(out_path, **out)
    print(f"\n  saved {out_path}")
    print(f"  peak RSS {S.peak_rss_gb():.2f} GB   wall time {wall:.0f} s", flush=True)


if __name__ == "__main__":
    main()
