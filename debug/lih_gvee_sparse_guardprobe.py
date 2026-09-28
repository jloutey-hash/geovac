"""Probe: guard-reference size vs trunc + the per-stage cost of the sparse==dense guard."""
from __future__ import annotations

import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
t0 = time.time()
import geovac.lih_r12ci.kernels as KG
import geovac.lih_r12ci.triangle as TRI
import lih_pi_m16 as M
import lih_gvee_sparse as S
print(f"import {time.time()-t0:.1f}s", flush=True)
KG.USE_EXACT_NEUMANN = True
ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
for tr in (3e-2, 5e-2, 1e-1, 2e-1):
    refm, info = M.build_ref_from_prims(ref, f"probe{tr}", cfg_tag='pi', trunc=tr, verbose=False)
    print(f"trunc {tr:.0e}: dets {len(info['kept'])} nact {refm.nact} nnz {len(refm.T_val)} "
          f"m!=0 {int(np.sum(refm.mA != 0))}", flush=True)
tr = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-1
refm, info = M.build_ref_from_prims(ref, "guard", cfg_tag='pi', trunc=tr, verbose=False)
t = time.time(); gem = M.GeminalM('linexp'); print(f"GeminalM {time.time()-t:.1f}s", flush=True)
t = time.time(); Fbar = M.fbar(gem, refm); print(f"fbar {time.time()-t:.1f}s", flush=True)
for nm in ('dense', 'sparse'):
    TRI.rng = np.random.default_rng(20260921); M._RSVD_RNG = np.random.default_rng(20260923)
    TRI._MODE_CACHE.clear(); M._TRIM_CACHE.clear()
    t = time.time()
    r = (M.hg_Vee_m(gem, refm, Fbar, force_general=True) if nm == 'dense'
         else S.hg_Vee_sparse(gem, refm, Fbar, 2.0))
    print(f"{nm} {time.time()-t:.1f}s  g={r[1]:.15f} F2={r[4]:.15f}", flush=True)
