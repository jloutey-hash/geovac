"""Confirm the CONVERGED c=3 box truncation on the graded mesh, independently.

My own uniform ladder stopped at 600k (2.060e-07) and was still climbing at
1.067x per step -- I read that as "no converged value" and recorded a 1.7x
disagreement between meshes.  That was wrong twice:

  * the ladder was still rising, not disagreeing -- extended to 1.5M it reaches
    2.14665e-07 with a step of 1.010x;
  * my "absolute" figure of 1.29e-07 used `ref = tprime_off(12 * n^2, 600000)`,
    a 768-bohr box at the SAME point count as the 192-bohr box it referenced --
    so the reference had 4x coarser spacing and carried ~9.2e-08 of its own
    grid error, comparable to the 2.16e-07 signal it was supposed to resolve.

The graded mesh `r = R t^2` resolves the near-origin region where the integrand
lives, so it converges at far fewer points.  This re-measures it here rather
than accepting the number, and cross-checks the c=5 absolute value against the
uniform route's converged figure -- if the two meshes really agree, those must
match to many digits.
"""
from __future__ import annotations

import numpy as np
import geovac.sturmian_secular as SS
import geovac.sturmian_variational as SV

Z, N = 2.0, 8


def tprime_off(box: float, npts: int, kind: str = "uni") -> float:
    SV.set_grid(box, npts, kind=kind)
    c = SS.build_configs(SV.family(N, 1))
    M = SS.build_M(c, Z=Z)
    D = np.diag([Z * x.Rnu for x in c])
    Tp = M - D
    return float(np.abs(Tp).sum() - np.abs(np.diag(Tp)).sum())


def main() -> None:
    print("GRADED mesh, c=3 vs c=5, relative -- does it converge?", flush=True)
    print(f"{'npts':>8} {'c3 abs':>18} {'c5 abs':>18} {'rel':>14} {'step':>8}",
          flush=True)
    prev = None
    for npts in (60_000, 120_000, 240_000):
        a3 = tprime_off(3.0 * N ** 2, npts, "grade")
        a5 = tprime_off(5.0 * N ** 2, npts, "grade")
        rel = abs(a3 - a5) / abs(a5)
        step = "--" if prev is None else f"{rel / prev:.3f}x"
        print(f"{npts:>8} {a3:>18.10f} {a5:>18.10f} {rel:>14.5e} {step:>8}",
              flush=True)
        prev = rel

    print(flush=True)
    print("CROSS-MESH CHECK -- the c=5 absolute value must agree if the meshes do:",
          flush=True)
    g5 = tprime_off(5.0 * N ** 2, 240_000, "grade")
    print(f"  graded  c=5 @240k : {g5:.10f}", flush=True)
    print(f"  uniform c=5 @1.5M : 17.2714838504  (reviewer's ladder)", flush=True)
    print(f"  relative gap      : {abs(g5 - 17.2714838504) / g5:.3e}", flush=True)

    print(flush=True)
    print("VERDICT", flush=True)
    print(f"  converged c=3 truncation (graded) = {prev:.5e}", flush=True)
    print(f"  uniform ladder at 1.5M            = 2.14665e-07", flush=True)
    print(f"  mesh-to-mesh spread               = "
          f"{abs(prev - 2.14665e-07) / prev * 100:.2f}%", flush=True)


if __name__ == "__main__":
    main()
