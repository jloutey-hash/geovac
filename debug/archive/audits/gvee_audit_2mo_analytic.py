r"""g_Vee defect audit, step A: analytic per-triple <f_p f_q / r_r> on the 2-MO sigma reference.

Uses the phase-2 sigma path (P2.vee_reduce, exact Neumann) at a chosen linexp/exp gamma via the
lih_exp_geminal_run source-override (kernels.py GAM).  Writes the 216 triple values + Vee/FVee
per-pair values to debug/data/gvee_audit_2mo_<gem>_g<gam>.npz.

    python debug/gvee_audit_2mo_analytic.py --geminal linexp --gam 0.25
"""
from __future__ import annotations

import os
import sys
import argparse
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
DATA = os.path.join(HERE, "data")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--geminal', default='linexp')
    ap.add_argument('--gam', type=float, default=0.25)
    ap.add_argument('--grid', type=int, nargs=3, default=None)
    a = ap.parse_args()
    import lih_exp_geminal_run as RUN
    KG = RUN.install_kernels(a.gam, a.grid)
    import lih_marriage_phase2 as P2
    assert P2.GAM == a.gam
    KG.USE_EXACT_NEUMANN = True
    t0 = time.time()
    ref = P2.build_ref_2mo()
    gem = P2.Geminal(a.geminal)
    PAIRS6 = P2.PAIRS6
    Fbar = P2.fbar(gem, ref)
    Vee_r = np.array([P2.vee_reduce([], r, 'coul', None, gem, ref) for r in PAIRS6])
    FV = np.array([[P2._vee_pair_fc(p, r, gem, ref) for r in PAIRS6] for p in PAIRS6])
    F2V = np.zeros((6, 6, 6))
    for i, p in enumerate(PAIRS6):
        for j, q in enumerate(PAIRS6):
            for k, r in enumerate(PAIRS6):
                F2V[i, j, k] = P2._vee_triple(p, q, r, gem, ref)
    Vee = Vee_r.sum(); FVee = FV.sum(); F2Vee = F2V.sum()
    g = F2Vee - 2 * Fbar * FVee + Fbar ** 2 * Vee
    print(f"{a.geminal} g={a.gam}: Fbar={Fbar:.10f} Vee={Vee:.10f} FVee={FVee:.10f} "
          f"F2Vee={F2Vee:.10f} g_Vee={g:.10f} h_Vee={FVee-Fbar*Vee:.10f}  [{time.time()-t0:.0f}s]")
    tag = f"{a.geminal}_g{a.gam:g}" + ("" if a.grid is None else "_%dx%dx%d" % tuple(a.grid))
    np.savez(os.path.join(DATA, f"gvee_audit_2mo_{tag}.npz"), Fbar=Fbar, Vee_r=Vee_r, FV=FV, F2V=F2V,
             g_Vee=g, pairs=np.array(PAIRS6))


if __name__ == "__main__":
    main()
