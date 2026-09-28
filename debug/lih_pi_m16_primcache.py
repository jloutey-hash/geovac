r"""Build + cache the M=16 (pi-containing) primitive engine matrices (S, h1, T1, eri) for the
LiH pi/M=16 marriage (plan debug/lih_pi_m16_build_plan.md, Phase A).  Mirrors
lih_marriage_phase0.build_primitives but with npi=1 and its OWN cache file so the sigma
(npi=0) cache debug/data/lih_marriage_phase0_prim.npz is NOT clobbered.

Run:  python debug/lih_pi_m16_primcache.py > debug/data/lih_pi_m16_primcache.log 2>&1
Writes debug/data/lih_pi_m16_prim.npz  (the ~324 s float64 ERI + one-body build).
"""
from __future__ import annotations

import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, "data")
CACHE = os.path.join(DATA, "lih_pi_m16_prim.npz")

CFG = dict(Jb=2, Lb=1, npi=1, Jpi=1, Lpi=0, alpha=1.0, core2=(4.5, 1.6), tol=1e-11)
Z_A, Z_B = 3.0, 1.0


def main() -> None:
    import lih_core2exp_probe as P
    from prolate_allelectron_c4 import one_body_general
    from prolate_float_eri import build_eri_tensor_m_f

    t0 = time.time()
    R = float(P.R)
    orbs, tags = P.lih_orbs(CFG['Jb'], CFG['Lb'], CFG['npi'], CFG['Jpi'], CFG['Lpi'],
                            CFG['alpha'], list(CFG['core2']))
    M = len(orbs)
    print(f"M={M} orbitals; tags:")
    for i, (key, j, l) in enumerate(tags):
        print(f"  {i:2d}: key={key} j={j} l={l}")
    S, h1 = one_body_general(orbs, R, Z_A, Z_B)
    S0, T1 = one_body_general(orbs, R, 0.0, 0.0)
    assert np.allclose(S, S0, atol=1e-13), np.max(np.abs(S - S0))
    print(f"one_body_general done [{time.time()-t0:.0f}s]", flush=True)
    eri = build_eri_tensor_m_f(orbs, R, verbose=True)
    print(f"ERI engine done [{time.time()-t0:.0f}s]", flush=True)
    # dtype note: the pi ERI engine may return complex; keep as-is.
    os.makedirs(DATA, exist_ok=True)
    np.savez(CACHE, cfg=repr(CFG), M=M, R=R, S=S, h1=h1, T1=T1, eri=eri,
             tags=np.array(tags, dtype=object))
    print(f"saved {CACHE}  (S {S.dtype}, eri {np.asarray(eri).dtype}, shape {np.asarray(eri).shape})")
    print(f"wall {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
