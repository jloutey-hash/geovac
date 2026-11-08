"""Pin the conical intersection of BeH2's configuration operator F=sum P_i and confirm the
Z2 Berry phase at multiple radii (Josh 2026-08-25).  Loads the cached overlaps (fast).

Expectation: on the symmetric line d1=d2 the swap symmetry sigma splits F into even/odd
sectors; an even and an odd band can CROSS (true degeneracy, gap->0) at a specific d* -- a
genuine conical intersection in the (d1,d2) plane.  A Berry loop around it should give PI at
EVERY radius that encloses it (a real topological invariant), 0 for loops that don't.
"""
from __future__ import annotations
import numpy as np
import beh2_bargmann_geometry as BG      # loads cached overlaps + projectors_geom


def gap_at(d1, d2):
    Ps = BG.projectors_geom(d1, d2)
    if Ps is None:
        return np.nan
    w = np.linalg.eigvalsh(Ps[0] + Ps[1] + Ps[2])
    return float(np.min(np.diff(w)))


def delta_at(d1, d2):
    Ps = BG.projectors_geom(d1, d2)
    return np.nan if Ps is None else float(np.trace(Ps[0] @ Ps[1] @ Ps[2]).real)


def berry(cx, cy, r, N=160):
    v0 = vp = None
    for i in range(N + 1):
        th = 2 * np.pi * i / N
        Ps = BG.projectors_geom(cx + r * np.cos(th), cy + r * np.sin(th))
        if Ps is None:
            return None
        w, V = np.linalg.eigh(Ps[0] + Ps[1] + Ps[2])
        if i == 0:
            k = int(np.argmin(np.diff(w))); v = V[:, k]; v0 = v.copy(); vp = v
        else:
            ov = V.T @ vp; j = int(np.argmax(np.abs(ov)))
            v = V[:, j] * np.sign(ov[j]); vp = v
    return float(np.sign(vp @ v0))


def main():
    # 1. locate the CI on the symmetric line d1=d2
    ds = np.linspace(1.5, 3.5, 500)
    g = np.array([gap_at(d, d) for d in ds])
    i = int(np.nanargmin(g)); dstar = ds[i]
    ds2 = np.linspace(dstar - 0.08, dstar + 0.08, 400)
    g2 = np.array([gap_at(d, d) for d in ds2])
    i2 = int(np.nanargmin(g2)); dstar = ds2[i2]; gmin = g2[i2]
    print(f"CI on symmetric line d1=d2: min F-gap = {gmin:.2e} at d* = {dstar:.5f}")
    print(f"  (gap 0.5 bohr away on the line: {gap_at(dstar - 0.5, dstar - 0.5):.4f} / "
          f"{gap_at(dstar + 0.5, dstar + 0.5):.4f})")
    print(f"  Delta at d*: {delta_at(dstar, dstar):+.5f}  "
          f"(Delta 0.3 away: {delta_at(dstar - 0.3, dstar - 0.3):+.4f} / "
          f"{delta_at(dstar + 0.3, dstar + 0.3):+.4f})  -> sign change across d* ?")

    # 2. Berry phase at multiple radii around the CI, plus far controls
    print("\nZ2 Berry phase around the CI (should be PI at EVERY enclosing radius):")
    for r in [0.12, 0.2, 0.3, 0.45, 0.6]:
        s = berry(dstar, dstar, r)
        print(f"  loop r={r:<4}: sign={s:+.0f}  => {'PI' if s and s < 0 else '0' if s else 'singular'}")
    print("far controls (should be 0):")
    for cx, cy, r in [(4.8, 4.8, 0.4), (dstar, dstar + 1.2, 0.3)]:
        s = berry(cx, cy, r)
        print(f"  loop @({cx:.2f},{cy:.2f}) r={r}: sign={s:+.0f}  => {'PI' if s and s < 0 else '0' if s else 'singular'}")


if __name__ == "__main__":
    main()
