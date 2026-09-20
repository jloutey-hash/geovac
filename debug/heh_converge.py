"""(ii) Confirm HeH+ converges toward the reference as the basis grows (p={0,1},
r12, full parity, mpf solve). If the error descends monotonically from above,
the exact-algebraic r12 machinery is confirmed to generalize to a heteronuclear
center quantitatively (not just structurally).
Run from root:  python debug/heh_converge.py
"""
import os
import sys
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geovac import prolate_recondition as pr
import prolate_r12_mpf as m
from heh_probe import build_basis_full, ZA, ZB, R, E_REF

ALPHA = 1.6
NUC = ZA * ZB / R


def energy(j_max, l_max):
    basis = build_basis_full(j_max, l_max, ALPHA, p_set=(0, 1))
    t = time.time()
    Smp, Hmp = m.assemble_hetero(basis, R, ALPHA, ZA, ZB, mpf_out=True)
    cond = np.linalg.cond(m._tofloat(Smp))
    e = m.solve_canonical_mpf(Smp, Hmp)[0]
    et = e + NUC
    print(f"(j={j_max},l={l_max}) n={len(basis):3d} cond={cond:.1e}  "
          f"E_tot={et:.6f}  err={(E_REF-et)*1000:+.3f} mHa  [{time.time()-t:.0f}s]",
          flush=True)
    return et


if __name__ == "__main__":
    print(f"HeH+ convergence (p={{0,1}}, a={ALPHA}, E_ref={E_REF}):", flush=True)
    for (jm, lm) in [(2, 2), (3, 2), (2, 3)]:
        energy(jm, lm)
