"""High-precision extension: chemistry vs the elliptic period class at weight 2.

The 46-digit run (debug/chemistry_elliptic_period_class.py) was DECISIVE-NEG at
weight 1 but UNDERPOWERED at weight 2 -- with |ring| = 10 the honest detectable
height at 46 digits is only ~10^4.6, and the decoy matched the real target
there, which is the signature of a search finding noise.

The 49 stored digits were a FORMATTING cap (`mp.nstr(..., 50)` inside
h2_pes_certify), not a computational limit. Recomputing the closed-form Newton
at higher dps lifts the ceiling: at ~130 digits a 10-term ring detects height
~10^13, so a null result becomes a genuine statement.

Digits are CERTIFIED by cross-precision agreement (two independent dps runs must
agree to the claimed digits), not merely printed.
"""

from __future__ import annotations

import os
import sys

import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.qfd_assemble import h2_closed_form_E
from routeC_T2_pslq_decisive import guarded_pslq
from chemistry_elliptic_period_class import ceiling, graded_ring

TARGET_DPS = 130


def constants(dps_list=(150, 170)):
    """Closed-form Newton on exact derivatives, at two precisions."""
    R = sp.Symbol("R", positive=True)
    E = h2_closed_form_E(R)
    fE = sp.lambdify(R, E, modules="mpmath")
    f = sp.lambdify(R, sp.diff(E, R), modules="mpmath")
    f2 = sp.lambdify(R, sp.diff(E, R, 2), modules="mpmath")

    runs = []
    for dps in dps_list:
        with mp.workdps(dps):
            x = mp.mpf("1.6")
            for _ in range(200):
                step = f(x) / f2(x)
                x -= step
                if abs(step) < mp.mpf(10) ** (-dps + 5):
                    break
            runs.append({"R_eq": +x, "D_e": mp.mpf(-1) - fE(x), "k": +f2(x)})
    return runs


def main():
    print("=== chemistry vs elliptic period class -- HIGH PRECISION ===")
    print(f"    recomputing H2 PES constants at dps 150/170 "
          f"(stored values were capped at 49 by nstr formatting)\n")

    runs = constants()
    lo, hi = runs

    targets = {}
    for key in ("R_eq", "D_e", "k"):
        with mp.workdps(TARGET_DPS + 40):
            diff = abs(hi[key] - lo[key])
            agree = mp.floor(-mp.log10(diff)) if diff > 0 else mp.inf
        print(f"  {key:5s} cross-precision agreement: ~{int(agree)} digits")
        if agree < TARGET_DPS:
            print(f"        !! only {int(agree)} < {TARGET_DPS} -- capping")
        targets[key] = hi[key]

    dps = min(TARGET_DPS, int(min(
        mp.floor(-mp.log10(abs(hi[k] - lo[k]))) for k in targets)) - 5)
    print(f"\n  certified working precision: {dps} dps")

    for wmax, names in ((2, ["pi", "vp"]),
                        (2, ["pi", "vp", "G"]),
                        (2, ["pi", "vp", "G", "P8"]),
                        (3, ["pi", "vp", "G"])):
        import chemistry_elliptic_period_class as base
        base.DPS = dps                      # ring built at the new precision
        ring = graded_ring(wmax, names)
        n = len(ring)
        hmax = 10.0 ** (dps / n)
        mc = int(max(10, min(10 ** 8, hmax)))
        print(f"\n{'='*66}")
        print(f"RING wt<={wmax} gens={names} |ring|={n}  "
              f"ceiling ~10^{dps/n:.2f}  searching to h<={mc}")

        # positive control from the ring's own elements
        keys = [k for k in ring if k != "1"][:3]
        planted = sum(c * ring[k] for c, k in zip([3, -5, 2], keys))
        rel = guarded_pslq(planted, ring, dps, mc, "CONTROL")
        if not (rel and max(abs(x) for x in rel) <= 40):
            print("   G2 FAIL -- skipping")
            continue
        print("   G2 PASS")

        for label, v in targets.items():
            with mp.workdps(dps + 20):
                decoy = v * (1 + mp.mpf(10) ** (-11)) + mp.euler / mp.mpf(10) ** 5
            print(f"\n   --- {label} ---")
            rr = guarded_pslq(+v, ring, dps, mc, f"{label} REAL")
            rd = guarded_pslq(+decoy, ring, dps, mc, f"{label} DECOY")
            hr = None if rr is None else max(abs(x) for x in rr)
            hd = None if rd is None else max(abs(x) for x in rd)
            if hr is None:
                print(f"       => DECISIVE-NEG: no relation with height <= {mc}")
            elif hd is not None and hd <= 4 * hr:
                print(f"       => INCONCLUSIVE: decoy matched "
                      f"(real {hr}, decoy {hd})")
            else:
                print(f"       => *** CANDIDATE *** real h={hr}, decoy h={hd}")


if __name__ == "__main__":
    main()
