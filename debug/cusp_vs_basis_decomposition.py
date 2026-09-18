"""Is H2's residual 0.41 mHa the electron-electron CUSP, or basis incompleteness?

WHY THIS RUNS BEFORE ANY EXPLICIT-CORRELATION WORK.  Paper 12 attributes the
residual at (5,5)+delta to "radial completeness and the true electron-electron
cusp" -- but it never decomposes the two, and the lever implied by each is
different.  If the residual is mostly basis incompleteness, explicit correlation
(geminals / Hy-CI / r12) is solving the wrong problem and a larger basis is the
answer, which is now cheap; if it is mostly the slow L^-3 partial-wave tail, then
explicit correlation is the right lever.  The corpus has ~15 ledgered negative
results in the explicit-correlation family (CLAUDE.md Sec. 3), so the standing
diagnostic-before-engineering rule says measure first.

The external anchor that makes the question sharp: Tao, McCurdy & Rescigno,
Phys. Rev. A 82, 023423 (2010) reach -1.17442 Ha (0.05 mHa) on H2 in THESE SAME
prolate spheroidal coordinates with a polynomial angular basis and NO r12
dependence anywhere.  So 0.05 mHa is demonstrably reachable here without explicit
correlation, and our 0.41 mHa is 8x above that -- which is evidence, before any
new calculation, that the cusp is not the binding constraint yet.

THE MEASUREMENT.  Raising j_max and l_max together conflates the two axes, so the
ladder moves ONE at a time from the (5,5)+delta anchor:

    (5,5,2)  anchor          -- reproduces the headline in this family
    (6,5,2)  radial  +1      -- dE attributable to radial incompleteness
    (5,6,2)  angular +1      -- dE attributable to angular truncation

Run in the `gegenbauer` family throughout: it spans the identical space but
conditions the downcast solve ~1000x tighter (9.1e4 vs 9.4e10 at (5,5)+delta),
which is what keeps a float64 solve trustworthy at larger truncation.  Both
families agree on the energy to <1 uHa, so the anchor is comparable with the
published laguerre_legendre headline.

READING THE RESULT.  Let dE_rad and dE_ang be the gains.  If dE_rad + dE_ang
accounts for most of the 0.41 mHa, the residual is basis-limited and the honest
next step is more basis (or a second exponent), NOT explicit correlation.  If both
are small and the energy has plateaued well above -1.174475, the remainder is the
partial-wave tail plus the cusp, and explicit correlation is on the critical path.
A plateau is only evidence of a cusp limit if the basis axes are BOTH saturated --
one axis still moving means the basis is still the limit.

Usage:  python debug/cusp_vs_basis_decomposition.py
        python debug/cusp_vs_basis_decomposition.py 6 5 2
"""
from __future__ import annotations

import sys
import time
from typing import List, Tuple

from geovac import prolate_recondition as pr

BASIS = "gegenbauer"
LADDER: List[Tuple[str, Tuple[int, int, int]]] = [
    ("anchor      ", (5, 5, 2)),
    ("radial  +1  ", (6, 5, 2)),
    ("angular +1  ", (5, 6, 2)),
]


def run(label: str, j: int, l: int, mu: int) -> None:
    t0 = time.time()
    r = pr.recondition_energy(j, l, mu, alpha=1.0, basis=BASIS)
    dt = time.time() - t0
    print(f"  {label} ({j},{l}) mu<={mu}  N={r.n_basis:5d} kept={r.n_kept:5d}  "
          f"E={r.energy:.7f}  D_e={r.de_pct:7.3f}%  |err|={abs(r.err_mha):6.3f} mHa  "
          f"cond={r.cond_norm:.2e}  var={r.variational}  [{dt:.0f}s]", flush=True)


if __name__ == "__main__":
    print("=== cusp vs basis decomposition (basis=%s, alpha=1.0) ===" % BASIS,
          flush=True)
    print("exact = %.6f Ha; D_e_exact = %.6f Ha; TMR 2010 reached 0.05 mHa here "
          "with NO r12" % (pr.E_EXACT, pr.DE_EXACT), flush=True)
    print(flush=True)
    if len(sys.argv) >= 4:
        run("single      ", int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
    else:
        for label, (j, l, mu) in LADDER:
            run(label, j, l, mu)
        print(flush=True)
        print("reading: a residual is only CUSP-limited if BOTH basis axes have "
              "saturated.", flush=True)
        print("one axis still moving means the basis is still the limit, and "
              "explicit", flush=True)
        print("correlation would be the wrong lever.", flush=True)
