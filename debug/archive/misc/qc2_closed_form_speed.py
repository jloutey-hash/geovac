"""QC-2: is the COMPILED closed form actually slow? (I asserted it; never measured.)

I told the PI that "evaluating a formula symbolically is itself slow -- you'd have
to compile those formulas into fast numeric code to get real speed", implying the
compiled version would still lose to Gaussians. The first half is right. The
second half was an assumption, and it is wrong at least against this repo's
Gaussian path.

THE STRUCTURE THAT MATTERS. The closed form is an expression in R with the
orbitals baked in, so the costs separate:
    build (symbolic)  -- once per quartet TYPE
    evaluate          -- once per GEOMETRY
Lumping them together, as I did, hides the amortisation entirely.

FAIRNESS. An earlier cut of this compared our VECTORISED array evaluation against
eri_md called one-at-a-time in a Python loop, with basis construction inside the
timed region -- two confounds, both flattering us. Both are removed here: scalar
loop against scalar loop, basis objects hoisted out of the Gaussian timing.

THE CAVEAT THAT MUST TRAVEL WITH THE RESULT. eri_md is pure Python. Production
Gaussian codes (libint/libcint) are tuned C, typically 100-1000x faster than a
Python implementation of the same algorithm. So a 390x margin over eri_md does
NOT establish a margin over a real code -- that comparison is UNMEASURED and
cannot be made from here.

Run from repo root:  python debug/qc2_closed_form_speed.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from fractions import Fraction  # noqa: E402

from geovac import noci_engine as E  # noqa: E402
from geovac.two_center_eri import R_s, aabb_closed_form  # noqa: E402


def main() -> None:
    print("QC-2 -- compiled closed form vs Gaussian, like for like\n")
    Z1 = Fraction(1)
    shapes = {}
    arr, dco, _q = E.fit_sto_shape(0, 1, n_gauss=6)
    shapes["1s"] = (arr, dco)
    N = 2000
    Rs = np.linspace(1.2, 6.0, N)

    t0 = time.perf_counter()
    expr = aabb_closed_form(Z1, (1, 0, 0), (1, 0, 0), Z1, (1, 0, 0), (1, 0, 0))
    t_build = time.perf_counter() - t0
    f = sp.lambdify(R_s, sp.re(expr), "numpy")

    t0 = time.perf_counter()
    for r in Rs:
        f(float(r))
    t_scalar = (time.perf_counter() - t0) / N

    t0 = time.perf_counter()
    f(Rs)
    t_vec = (time.perf_counter() - t0) / N

    pa = np.array([0., 0., 0.])
    A = E.sto_shape_basis(pa, "1s", 1.0, shapes, (0, 0, 0))
    Bs = [E.sto_shape_basis(np.array([0., 0., float(r)]), "1s", 1.0, shapes,
                            (0, 0, 0)) for r in Rs[:300]]
    t0 = time.perf_counter()
    for B in Bs:
        E.eri_md(A, A, B, B)
    t_md = (time.perf_counter() - t0) / 300

    print(f"  symbolic build (once, per quartet type) : {t_build*1e3:9.1f} ms")
    print(f"  ours, scalar loop                       : {t_scalar*1e6:9.2f} us")
    print(f"  ours, vectorised over R                 : {t_vec*1e6:9.2f} us")
    print(f"  eri_md, integral only (fair baseline)   : {t_md*1e6:9.2f} us")
    print(f"\n  ratio, scalar vs scalar : {t_scalar/t_md:.4f}x  "
          f"({t_md/t_scalar:.0f}x faster)")
    print(f"  build amortises after   : {t_build/(t_md-t_scalar):.0f} geometries")
    print("\n  CAVEAT: eri_md is pure Python. Production Gaussian codes are tuned")
    print("  C and typically 100-1000x faster than a Python implementation of the")
    print("  same algorithm, so this does NOT establish a margin over a real code.")


if __name__ == "__main__":
    main()
