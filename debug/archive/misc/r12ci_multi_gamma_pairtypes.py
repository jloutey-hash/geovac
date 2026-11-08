"""Does a SECOND correlation length buy more for a 3-electron atom than a 2-electron one?

Hypothesis under test (raised when Li's gamma_opt shifted 2.4x from Li+'s): a single global
gamma is the wrong ansatz for many-electron atoms because the distinct orbital-pair types
(1s1s tight, 1s2s / 2s2s diffuse) want different correlation ranges.  gamma cannot be tied to
a pair directly -- electrons are indistinguishable -- but giving the CI several geminals with
different gamma lets it BUILD a per-pair-type effective range by superposition.

  PREDICTION IF TRUE:  the 2nd geminal buys much more for Li (3 pair types) than for
                       Li+ (1 pair type, 1s^2).
  IF FALSE:            both gain similarly, and the gamma_opt shift needs another explanation.

CONTROL: two geminals with the SAME gamma is a redundant column and must buy ~nothing.
"""
import io
import itertools
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)

import r12ci_ne_engine as E  # noqa: E402

NG, LMAX, NX = 320, 10, 48
GAMS = [0.15, 0.35, 0.65, 1.10, 1.80]
CASES = [("Li+", 2, 0, 3.0, 4, 2.4), ("Li", 3, 1, 3.0, 4, 1.7)]

out = {"grid": dict(Ng=NG, Lmax=LMAX, nx=NX), "gammas": GAMS, "systems": {}}

for name, N, ms2, Z, ns, k in CASES:
    t0 = time.time()
    print("=" * 80)
    print(f"{name}:  N={N}, Z={Z}, ns={ns}, k={k}")
    print("=" * 80)
    S0, H0, dets = E.build(ns, N, Z, k, GAMS[0], ms2, Ng=NG, with_geminal=False,
                           Lmax=LMAX, nx=NX)
    e_plain, _ = E.solve(S0, H0)
    nd = len(dets)

    singles = {}
    for g in GAMS:
        S1, H1, _ = E.build(ns, N, Z, k, [g], ms2, Ng=NG, Lmax=LMAX, nx=NX)
        singles[g], _ = E.solve(S1, H1, nd=nd)
    best1_g = min(singles, key=singles.get)
    best1 = singles[best1_g]
    print(f"  plain           = {e_plain:.6f}")
    print("  1 geminal:      " + "  ".join(f"g={g}: {1000*(e_plain-singles[g]):.2f}"
                                           for g in GAMS) + "   (gain, mHa)")
    print(f"  best single     = {best1:.6f}  at gamma={best1_g}  "
          f"gain {1000*(e_plain-best1):.2f} mHa")

    # CONTROL: duplicate gamma must buy ~nothing
    Sc, Hc, _ = E.build(ns, N, Z, k, [best1_g, best1_g], ms2, Ng=NG, Lmax=LMAX, nx=NX)
    e_ctrl, _ = E.solve(Sc, Hc, nd=nd)
    print(f"  CONTROL (gamma duplicated) = {e_ctrl:.6f}   extra gain "
          f"{1000*(best1-e_ctrl):.4f} mHa  (must be ~0)")

    pairs_res = {}
    for ga, gb in itertools.combinations(GAMS, 2):
        S2, H2, _ = E.build(ns, N, Z, k, [ga, gb], ms2, Ng=NG, Lmax=LMAX, nx=NX)
        pairs_res[(ga, gb)], _ = E.solve(S2, H2, nd=nd)
    bp = min(pairs_res, key=pairs_res.get)
    best2 = pairs_res[bp]
    print(f"  best pair       = {best2:.6f}  at gamma={bp}  "
          f"gain {1000*(e_plain-best2):.2f} mHa")
    extra = 1000 * (best1 - best2)
    print(f"  EXTRA from the 2nd correlation length = {extra:.3f} mHa")
    print(f"  ({time.time()-t0:.0f} s)")
    print()
    out["systems"][name] = dict(
        N=N, k=k, plain=e_plain, singles={str(g): singles[g] for g in GAMS},
        best1=best1, best1_gamma=best1_g, control_dup=e_ctrl,
        control_extra_mHa=1000 * (best1 - e_ctrl),
        pairs={f"{a},{b}": v for (a, b), v in pairs_res.items()},
        best2=best2, best2_gammas=list(bp), extra_mHa=extra)

a, b = out["systems"]["Li+"], out["systems"]["Li"]
print("=" * 80)
print("VERDICT")
print("=" * 80)
print(f"  extra gain from a 2nd correlation length:")
print(f"    Li+ (2e, 1 pair type)   {a['extra_mHa']:.3f} mHa")
print(f"    Li  (3e, 3 pair types)  {b['extra_mHa']:.3f} mHa")
if a["extra_mHa"] > 1e-9:
    print(f"    ratio Li / Li+ = {b['extra_mHa'] / a['extra_mHa']:.2f}x")
print(f"  duplicate-gamma controls: Li+ {a['control_extra_mHa']:.4f}, "
      f"Li {b['control_extra_mHa']:.4f} mHa (must be ~0)")
out["verdict"] = dict(extra_Liplus=a["extra_mHa"], extra_Li=b["extra_mHa"])

os.makedirs("debug/data", exist_ok=True)
with io.open("debug/data/r12ci_multi_gamma_pairtypes.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/r12ci_multi_gamma_pairtypes.json")
