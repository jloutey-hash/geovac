"""Close Paper 58's one uncorroborated census leg: the n_max = 3 g row.

Paper 58 states the permitted-density inflation as 13.8x at n_max=2 and 15.0x at
n_max=3, and labels the second "counted only" -- the MD corroboration ran at
n_max=2 because that basis needs only s and p functions.

n_max = 3 needs 3d, and geovac/noci_engine.py describes itself as an
"s/p" McMurchie-Davidson engine. Empirically d functions DO construct and
evaluate (the Hermite machinery is general), but d is undocumented and
untested there, so this script VALIDATES d before using it.

Part 1 -- d validation against an independent route.
  For a primitive Gaussian g_s = exp(-a|r-A|^2),
      d/dA_z   g_s = 2a (z-A_z) g_s            = 2a * g_pz
      d2/dA_z2 g_s = (4a^2 (z-A_z)^2 - 2a) g_s = 4a^2 * g_dzz - 2a * g_s
  so   g_dzz = ( d2/dA_z2 g_s + 2a g_s ) / (4 a^2).
  Every d integral is therefore reachable from SECOND CENTRE-DERIVATIVES of s
  integrals alone -- the same trick the engine's own suite uses to validate p
  against s. If MD's d agrees with that, d is trustworthy here.

Part 2 -- the n_max = 3 census, if and only if Part 1 passes.

--- RESULT 2026-08-11: Part 1 PASSED, Part 2 IS INVALID. DO NOT CITE ITS NUMBER.

Part 1 succeeded: MD's d path agrees with the independent FD-from-s route to
7.8e-9 (overlap) and 5.5e-8 (kinetic), near the 2nd-order FD floor. **d functions
in geovac/noci_engine.py are trustworthy despite the "s/p" docstring** -- that is
a genuine reusable finding.

Part 2 ran and reported 154176/614656 = 25.08% dense against the paper's counted
18.59%, a ratio of 1.349. That 35% gap is NOT a basis-convention artifact: the
n_max=2 comparison agreed to 1.4% (29.8 measured vs 29.4 counted). The gap is a
BUG IN THIS SCRIPT'S BASIS.

Diagnosis: SHELLS_NMAX3 uses five raw Cartesian d monomials -- (1,1,0), (1,0,1),
(0,1,1), (2,0,0), (0,0,2), i.e. dxy, dxz, dyz, dxx, dzz -- as a stand-in for the
five real l=2 harmonics. That is wrong twice over:
  * Cartesian d has SIX components; the pure l=2 set requires traceless
    combinations (d_z2 ~ 2zz-xx-yy, d_x2y2 ~ xx-yy). Raw xx and zz each carry
    an l=0 admixture, so this basis is not an l=2 set at all.
  * It includes dxx and dzz but omits dyy, so the basis is not even axially
    symmetric about z -- it cannot respect the m-structure the census counts.
Comparing its nonzero count against a complex-m symmetry census is therefore
meaningless.

Proper fix: build the five real harmonics as linear combinations of Cartesians
and expand every ERI over those combinations. That is correct but expensive here
-- a pure-d quartet becomes a sum over up to 3^4 = 81 Cartesian terms, so the
82,621 unique quartets would balloon to millions of eri_md calls (hours, not
minutes). The single-lmn interface of eri_md does not express contracted
harmonics directly.

STATUS: Paper 58's n_max=3 g-row figure remains COUNTED-ONLY, uncorroborated.
That labelling in the paper is correct and unchanged. Do not "fix" the paper on
the strength of the 25.08% figure above.

Run from repo root:  python debug/p58_eri_census_nmax3.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac import noci_engine as E  # noqa: E402

TOL = 1e-10
Z_A, Z_B, R = 3.0, 1.0, 3.0

# n_max = 3 per center: 1s 2s 2p(3) 3s 3p(3) 3d(5-ish via Cartesian 6, we use
# the 5 canonical Cartesian d used by the builder's real-harmonic count) = 14.
SHELLS_NMAX3 = [
    ("1s", "1s", (0, 0, 0)),
    ("2s", "2s", (0, 0, 0)),
    ("2px", "2p", (1, 0, 0)), ("2py", "2p", (0, 1, 0)), ("2pz", "2p", (0, 0, 1)),
    ("3s", "3s", (0, 0, 0)),
    ("3px", "3p", (1, 0, 0)), ("3py", "3p", (0, 1, 0)), ("3pz", "3p", (0, 0, 1)),
    ("3dxy", "3d", (1, 1, 0)), ("3dxz", "3d", (1, 0, 1)),
    ("3dyz", "3d", (0, 1, 1)), ("3dx2", "3d", (2, 0, 0)),
    ("3dz2", "3d", (0, 0, 2)),
]


# --------------------------------------------------------------------------
# Part 1: validate d against second center-derivatives of s
# --------------------------------------------------------------------------

def validate_d() -> bool:
    print("Part 1 -- validating d functions against 2nd centre-derivatives of s\n")
    a = 0.83
    al, co = np.array([a]), np.array([1.0])
    A = np.array([0.10, -0.20, 0.30])
    Bc = np.array([0.50, 0.40, -1.10])
    b_s = E.BasisFn(Bc, (0, 0, 0), np.array([1.27]), co)
    eps = 1e-4

    def s_at(Az):
        return E.BasisFn(np.array([A[0], A[1], Az]), (0, 0, 0), al, co)

    # normalization bookkeeping: BasisFn renormalizes, so work with the RATIO
    # of the analytic prefactors between the normalized d and normalized s.
    d_f = E.BasisFn(A, (0, 0, 2), al, co)
    s_f = E.BasisFn(A, (0, 0, 0), al, co)
    Nd, Ns = d_f.coeffs[0], s_f.coeffs[0]

    worst = 0.0
    for name, fn in (("overlap", E.overlap_md), ("kinetic", E.kinetic_md)):
        f_p, f_0, f_m = (fn(s_at(A[2] + eps), b_s), fn(s_at(A[2]), b_s),
                         fn(s_at(A[2] - eps), b_s))
        d2 = (f_p - 2 * f_0 + f_m) / eps ** 2
        # unnormalized-s second derivative -> unnormalized d_zz, then renormalize
        pred = (d2 / Ns + 2 * a * (f_0 / Ns)) / (4 * a ** 2) * Nd
        got = fn(d_f, b_s)
        err = abs(pred - got) / max(abs(got), 1e-30)
        worst = max(worst, err)
        print(f"  {name:<9} FD-from-s = {pred: .12e}   MD d = {got: .12e}"
              f"   rel err = {err:.2e}")

    ok = worst < 1e-5
    print(f"\n  worst relative error = {worst:.2e}  -> "
          f"{'d VALIDATED' if ok else 'd NOT TRUSTWORTHY'}")
    if not ok:
        print("  (FD is 2nd-order in eps=1e-4, so ~1e-6 is the floor; a failure")
        print("   here means the MD d path is wrong, not that FD is imprecise.)")
    return ok


# --------------------------------------------------------------------------
# Part 2: n_max = 3 census
# --------------------------------------------------------------------------

def build(shapes):
    pa, pb = np.array([0., 0., 0.]), np.array([0., 0., R])
    orbs, side = [], []
    for tag, pos, zeta in (("A", pa, Z_A), ("B", pb, 1.0)):
        for _lab, kind, lmn in SHELLS_NMAX3:
            orbs.append(E.sto_shape_basis(pos, kind, zeta, shapes, lmn))
            side.append(tag)
    return orbs, side


def census(orbs):
    m = len(orbs)
    nz = 0
    t0 = time.time()
    seen = {}
    for i in range(m):
        for j in range(i, m):
            for k in range(m):
                for l in range(k, m):
                    if (i, j) > (k, l):
                        continue
                    v = E.eri_md(orbs[i], orbs[j], orbs[k], orbs[l])
                    seen[(i, j, k, l)] = v
    # expand the 8-fold symmetry into a nonzero count over all m^4 slots
    for (i, j, k, l), v in seen.items():
        if abs(v) <= TOL:
            continue
        slots = {(i, j, k, l), (j, i, k, l), (i, j, l, k), (j, i, l, k),
                 (k, l, i, j), (l, k, i, j), (k, l, j, i), (l, k, j, i)}
        nz += len(slots)
    print(f"  unique quartets evaluated: {len(seen)}   in {time.time()-t0:.0f}s")
    return nz, m ** 4


def main() -> None:
    if not validate_d():
        print("\nSTOP: not running the n_max=3 census on an unvalidated d path.")
        return

    print("\nPart 2 -- n_max = 3 census (Z_A=3 / Z_B=1, R=3)\n")
    shapes = {}
    for kind, (l, n_r) in (("1s", (0, 1)), ("2s", (0, 2)), ("2p", (1, 2)),
                           ("3s", (0, 3)), ("3p", (1, 3)), ("3d", (2, 3))):
        arr, dco, q = E.fit_sto_shape(l, n_r)
        shapes[kind] = (arr, dco)
        print(f"  fit {kind}: <fit|STO> = {q:.7f}")

    orbs, _side = build(shapes)
    print(f"  M = {len(orbs)} orbitals\n")
    nz, dense = census(orbs)
    density = 100.0 * nz / dense
    print(f"\n  nonzero = {nz} / {dense}  ({density:.2f}% dense)")
    print(f"  paper's COUNTED value at n_max=3: 114280 / 614656 = 18.59%")
    print(f"  measured / counted ratio = {nz / 114280:.4f}")

    out = {"n_max": 3, "M": len(orbs), "nonzero": nz, "dense": dense,
           "density_pct": density, "counted_permitted": 114280,
           "counted_density_pct": 18.59, "tol": TOL,
           "caveat": "Gaussian-fitted floats, real Cartesian basis; the paper's "
                     "count is complex-m, and the two differ by the +/-m mixing."}
    p = REPO / "debug" / "data" / "p58_eri_census_nmax3.json"
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
