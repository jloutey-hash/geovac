"""Independently extend the c=3 box-rule ladder past the point I called a plateau.

DELTA #2's code reviewer sampled npts=400000 -- a point nobody had taken -- and
reports that it BREAKS the plateau I asserted on 2026-09-13, with the column
climbing 1.568e-7 -> 1.931e-7 -> 2.060e-7 toward a converged 2.163e-7.

It is also right that MY OWN printed table already contradicted me: it shows
2.0598e-07 at 600k against 1.5676e-07 at 250k -- a 31% RISE -- printed as a
"0.76x step" directly beneath the word PLATEAU, and I did not reconcile it.

This driver re-measures rather than accepting, because the last time a
measurement-based finding arrived I accepted it on the reviewer's sample and
was wrong.  Two routes:

  (A) the in-repo relative quantity, |T'|_1^off(c) vs the c=5 box at the SAME
      npts -- the thing the test measures;
  (B) the ABSOLUTE box truncation, each box against a far larger reference, so
      the moving reference cannot mask or manufacture a plateau.  Route (A)'s
      reference is itself npts-dependent, which is the reviewer's stated
      mechanism for the false flat spot.
"""
from __future__ import annotations

import numpy as np
import geovac.sturmian_secular as SS
import geovac.sturmian_variational as SV

Z_HE, NMAX, LMAX = 2.0, 8, 1


def tprime_off(box: float, npts: int) -> float:
    SV.set_grid(box, npts)
    cfgs = SS.build_configs(SV.family(NMAX, LMAX))
    M = SS.build_M(cfgs, Z=Z_HE)
    D = np.diag([Z_HE * c.Rnu for c in cfgs])
    Tp = M - D
    return float(np.abs(Tp).sum() - np.abs(np.diag(Tp)).sum())


def main() -> None:
    LAD = (100_000, 250_000, 400_000, 600_000)
    print("ROUTE A -- in-repo relative quantity: |c=3 box - c=5 box| / |c=5 box|")
    print(f"{'npts':>9} {'c3_abs':>18} {'c5_ref':>18} {'rel':>13} {'step':>8}")
    prev = None
    relA = {}
    for npts in LAD:
        a3 = tprime_off(3.0 * NMAX ** 2, npts)
        a5 = tprime_off(5.0 * NMAX ** 2, npts)
        rel = abs(a3 - a5) / abs(a5)
        relA[npts] = rel
        step = "--" if prev is None else f"{rel / prev:.3f}x"
        print(f"{npts:>9} {a3:>18.10f} {a5:>18.10f} {rel:>13.5e} {step:>8}")
        prev = rel

    print()
    print("ROUTE B -- ABSOLUTE truncation vs a c=12 box at high resolution")
    ref = tprime_off(12.0 * NMAX ** 2, 600_000)
    print(f"reference |T'|_1^off (c=12, 600k) = {ref:.10f}")
    print(f"{'c':>4} {'abs value':>18} {'rel vs ref':>13}")
    for c in (1.0, 2.0, 3.0, 5.0, 8.0):
        v = tprime_off(c * NMAX ** 2, 600_000)
        print(f"{c:>4.0f} {v:>18.10f} {abs(v - ref) / abs(ref):>13.5e}")

    print()
    print("VERDICT INPUTS")
    lo, hi = relA[100_000], relA[600_000]
    print(f"  route A 100k -> 600k: {lo:.5e} -> {hi:.5e}  ({hi / lo:.3f}x)")
    print(f"  is the 100k..250k pair the ONLY one under a 10% window?")
    ks = list(LAD)
    for a, b in zip(ks, ks[1:]):
        d = abs(relA[b] - relA[a]) / relA[a]
        print(f"    {a:>7} -> {b:>7}:  {d * 100:6.2f}%   {'PASSES <10%' if d < 0.10 else 'fails'}")


if __name__ == "__main__":
    main()
