"""LiH R12 F12-SYNTHESIS CEILING PROBE (diagnostic-before-engineering, 2026-09-22).

Question the diagnostic answers: how much of the Route-C residual can an r12 geminal
recover before the expensive two-center 4e mpf build?  Route C plateaus at E=-8.012
(58 mHa above exact) at the FCI determinant wall; the dominant residual is the Li 1s^2
core-core correlation, which the single analytic core orbital leaves at ZERO (so the
geminal's core contribution is DOUBLE-COUNTING-FREE).  The isolated Li+ (He-like Z=3)
core is the clean, transferable proxy for that pair.

Model: {Phi0, b_k Phi0} generalized eigenproblem, Phi0 = 1s(zeta)^2 with zeta=2.6875
(= Route C's frozen analytic core exponent -- NOT re-optimized, the realistic scenario),
in Hylleraas coordinates (r1, r2, r12); volume element 8 pi^2 r1 r2 r12; kinetic energy
via the symmetric gradient form (first derivatives only; the corpus 'bounded gradient').
Hylleraas primitives b(r1,r2,r12) in {1, u=r12, u2, t2=(r1-r2)^2, s=r1+r2, ut2}.

Finding (see debug/sprint_lih_r12_ceiling_diagnostic_memo.md):
  - r12-only {u,u2}      -> ~52% of the core deficit (PLATEAUS; more u-terms do not help)
  - radial-only {t2,s}   -> ~45-51%  (an ORBITAL degree of freedom, NOT a geminal one)
  - both {u,t2}          -> ~74-79%  (near-additive: the two halves are complementary)
  - rich {u,t2,s,u2,ut2} -> ~95-98%  (full core correlation with a small correlated basis)
So a geminal ALONE caps at ~half the core correlation; the other half needs a 2nd core
ORBITAL (radial flexibility), which Route C's ladder can supply with existing code.

Validated against known He (Z=2) Hylleraas: {u,t2}=-2.892 (classic 3-term -2.9024),
5-term -> 98% of the 57 mHa He correlation.  Float64, ~seconds.
"""
# The grids, integ, basisfn and energy moved to geovac/prolate_allelectron/helike_cusp.py
# (LiH productionization step 2a, 2026-10-04); re-exported here.
import os as _os
import sys as _sys

_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _ROOT not in _sys.path:
    _sys.path.insert(0, _ROOT)

from geovac.prolate_allelectron.helike_cusp import *        # noqa: E402,F401,F403
from geovac.prolate_allelectron.helike_cusp import R1, R2, integ, energy   # noqa: E402,F401
import numpy as np                                          # noqa: E402


# ---- quadrature self-check (norm=1, <1/r1>=zeta, T00=zeta^2) --------------- #
def _validate(zeta=2.6875, Z=3.0):
    N = zeta**3 / np.pi
    Phi2 = (N * np.exp(-zeta * (R1 + R2)))**2
    norm = integ(Phi2)
    inv_r1 = integ(Phi2 / R1) / norm
    _, _, _ = energy(Z, zeta, ['u'])
    print("quadrature self-check (zeta=2.6875):")
    print(f"  norm={norm:.8f} (1)   <1/r1>={inv_r1:.6f} (zeta={zeta})   "
          f"E0={energy(Z, zeta, ['u'])[0]:.6f} (zeta^2-2Z zeta+5zeta/8={zeta**2-2*Z*zeta+5*zeta/8:.6f})\n")


CASES = {'He (Z=2)': (2.0, 27/16, -2.90372), 'Li+ (Z=3)': (3.0, 2.6875, -7.27991)}
SETS = [
    ("r12 only        {u}",          ['u']),
    ("r12 only        {u,u2}",       ['u', 'u2']),
    ("radial only     {t2}",         ['t2']),
    ("radial only     {t2,s}",       ['t2', 's']),
    ("BOTH  {u,t2}",                 ['u', 't2']),
    ("BOTH  {u,t2,s,u2,ut2}",        ['u', 't2', 's', 'u2', 'ut2']),
]
if __name__ == "__main__":
    _validate()
    for name, (Z, zeta, exact) in CASES.items():
        E0, _, _ = energy(Z, zeta, ['u'])
        corr = E0 - exact
        print("=" * 76)
        print(f"{name}: ref E0={E0:.5f} exact={exact} corr={corr*1e3:.1f} mHa")
        print("=" * 76)
        for label, ks in SETS:
            _, ER, cond = energy(Z, zeta, ks)
            print(f"  {label:26s} E={ER:.5f}  rec={ (E0-ER)*1e3:6.1f} mHa  {100*(E0-ER)/corr:5.0f}%  cond={cond:.0e}")
        print()
