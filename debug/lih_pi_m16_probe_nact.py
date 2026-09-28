r"""Cheap probe: active-pair count nact(trunc) for the M=16 pi reference, and the peak-memory
estimate of the g_Vee kept-pair triple field (~1.4 * nact^3 * NG * 8 bytes).  NO g_Vee is run;
this only builds the reference (fast, cached primitives) so we can pick the tightest truncation
that fits RAM before launching the long energy compute.

Run:  python debug/lih_pi_m16_probe_nact.py
"""
from __future__ import annotations

import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)

import geovac.lih_r12ci.kernels as KG
import lih_pi_m16 as M

NG = KG.NXI * KG.NETA
OVERHEAD = 1.4          # observed peak / single-array (agent: ~9 GB at nact=64 -> 6.6 GB array)


def est_peak_gb(nact: int) -> float:
    return OVERHEAD * (nact ** 3) * NG * 8 / 1e9


def main() -> None:
    print(f"NG = NXI*NETA = {KG.NXI}*{KG.NETA} = {NG};  peak ~ {OVERHEAD} * nact^3 * NG * 8 bytes")
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI)
    print(f"M={ref['M']} Mk={ref['Mk']}")
    print(f"{'trunc':>8} {'ndet':>7} {'nact':>6} {'npi(m!=0)':>10} {'est_peak_GB':>12}")
    for trunc in (1e-2, 8e-3, 6e-3, 5e-3, 4e-3, 3e-3, 2e-3, 1e-3):
        refm, info = M.build_ref_from_prims(ref, "probe", cfg_tag='pi', trunc=trunc)
        npi = int(np.sum(refm.mA != 0))
        print(f"{trunc:>8.0e} {len(info['kept']):>7} {refm.nact:>6} {npi:>10} "
              f"{est_peak_gb(refm.nact):>12.1f}", flush=True)


if __name__ == "__main__":
    main()
