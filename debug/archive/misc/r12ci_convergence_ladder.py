"""Verification ladder for the r12 engine H2 energy.
Checks (i) monotone convergence from ABOVE toward E_EXACT as the basis grows,
and (ii) mpf-solve == float64-solve (conditioning not yet corrupting float64).
Run from root:  python debug/r12ci_convergence_ladder.py
"""
import os
import sys
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geovac import prolate_recondition as pr
import prolate_r12_mpf as m
from r12ci_first_energy import build_basis, solve_canonical, NUC, E_EXACT, DE_EXACT

R, ALPHA = 1.4011, 1.0


def one(jm, lm, do_mpf):
    basis = build_basis(jm, lm)
    n = len(basis)
    t0 = time.time()
    Smp, Hmp = m.assemble_mixed(basis, R, ALPHA, mpf_out=True)
    Sf, Hf = m._tofloat(Smp), m._tofloat(Hmp)
    cond = np.linalg.cond(Sf)
    ef = solve_canonical(Sf, Hf)
    etf = ef[0] + NUC
    line = (f"(j={jm},l={lm}) n={n:3d} cond={cond:.1e}  "
            f"float: E_tot={etf:.7f} D_e%={(-1.0-etf)/DE_EXACT*100:.3f} "
            f"err={(E_EXACT-etf)*1000:+.4f}mHa kept={ef[2]}/{n}")
    if do_mpf:
        em = m.solve_canonical_mpf(Smp, Hmp)
        etm = em[0] + NUC
        line += (f"  || mpf: E_tot={etm:.7f} "
                 f"err={(E_EXACT-etm)*1000:+.4f}mHa kept={em[2]}/{n} "
                 f"(dfloat={abs(etm-etf)*1e6:.2f}uHa)")
    line += f"  [{time.time()-t0:.0f}s]"
    print(line, flush=True)


if __name__ == "__main__":
    # bounded scaling ladder (optimized odd routing + mpf-orthogonalized solve):
    # grow radial x angular toward uHa; mpf cross-check throughout.  Capped where
    # the mpf S-orthogonalization (O(n^3)) stays tractable.
    for (jm, lm) in [(3, 2), (2, 4), (3, 3), (4, 3), (3, 4)]:
        one(jm, lm, True)
