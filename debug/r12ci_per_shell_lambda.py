"""EXPERIMENT A: does a per-shell lambda (multi-zeta) close the Li accuracy gap?

Li's 71-86 mHa shortfall was traced to the SHARED exponent: one k must cover both the
tight 1s^2 core and the diffuse 2s valence, and moving that single knob swings the energy
by 542 mHa.  The corpus's multi-lambda machinery (geovac/shibuya_wulfman.py
`_hydrogenic_poly_coeffs_lam`) shows mixed exponents keep the closed-form integrals, so
the usual objection to multi-zeta does not obviously apply.

DISCRIMINATING CONTROL -- He.  Both He electrons live in ONE shell, so there is only one
length scale and per-shell lambda should buy nearly NOTHING.  Li has two shells and should
gain a lot.  If BOTH gain similarly, the "two length scales" explanation is wrong and the
improvement is just "more variational parameters".

Honest ceiling: this basis is s-only, so neither system can reach its exact energy --
angular correlation is absent.  For He the s-limit is known (-2.8790288), which calibrates
how much room is actually left.
"""
import io
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)

import r12ci_ne_engine as E  # noqa: E402

EXACT = {"He": -2.9037243770, "Li+": -7.2799134126, "Li": -7.4780603236}
S_LIMIT = {"He": -2.8790288}          # exact l=0-only limit; Li s-limit not quoted
NG, NX = 240, 48

# (name, N, ms2, Z, ns, shared-k scan grid)
CASES = [
    ("He", 2, 0, 2.0, 4, np.linspace(1.2, 2.6, 15)),
    ("Li+", 2, 0, 3.0, 4, np.linspace(1.8, 3.6, 15)),
    ("Li", 3, 1, 3.0, 4, np.linspace(1.0, 3.0, 15)),
]


def energy(ns, N, Z, ms2, lams, Ng=NG):
    S, H, _ = E.build(ns, N, Z, float(np.mean(lams)), 0.5, ms2, Ng=Ng,
                      with_geminal=False, nx=NX, lowdin=False, lams=list(lams))
    e, _ = E.solve(S, H)
    return e


out = {"grid": dict(Ng=NG, nx=NX), "exact": EXACT, "systems": {}}
for name, N, ms2, Z, ns, kgrid in CASES:
    t0 = time.time()
    print("=" * 78)
    print(f"{name}:  N={N}, Z={Z}, ns={ns}   (s-only basis)")
    print("=" * 78)

    # --- baseline: best SHARED exponent -------------------------------------
    es = [energy(ns, N, Z, ms2, [k] * ns) for k in kgrid]
    i0 = int(np.argmin(es))
    k_best, e_shared = float(kgrid[i0]), float(es[i0])
    print(f"  shared exponent : best k = {k_best:.3f}   E = {e_shared:.6f}")

    # --- free per-shell lambda ----------------------------------------------
    def obj(u):
        try:
            return energy(ns, N, Z, ms2, np.exp(u))
        except Exception:
            return 1e6

    best = None
    for start in ([k_best] * ns,
                  [k_best * 1.6] + [k_best * 0.6] * (ns - 1),
                  [Z * 0.9] + [Z / (i + 2) for i in range(ns - 1)]):
        res = minimize(obj, np.log(start), method="Nelder-Mead",
                       options=dict(maxiter=600, xatol=1e-4, fatol=1e-9))
        if best is None or res.fun < best.fun:
            best = res
    lam_opt = np.exp(best.x)
    e_free = float(best.fun)
    print(f"  per-shell lambda: {np.array2string(lam_opt, precision=3)}")
    print(f"                    E = {e_free:.6f}")
    gain = 1000 * (e_shared - e_free)
    print(f"  GAIN from per-shell lambda = {gain:.2f} mHa")
    print(f"  spread lambda_max/lambda_min = {lam_opt.max() / lam_opt.min():.2f}x")
    if name in S_LIMIT:
        print(f"  (s-only limit {S_LIMIT[name]:.6f}; remaining s-room "
              f"{1000 * (e_free - S_LIMIT[name]):.2f} mHa)")
    print(f"  vs exact: shared {1000*(e_shared-EXACT[name]):.1f} mHa, "
          f"free-lambda {1000*(e_free-EXACT[name]):.1f} mHa  ({time.time()-t0:.0f}s)")
    print()
    out["systems"][name] = dict(N=N, Z=Z, ns=ns, k_best=k_best, E_shared=e_shared,
                                lam_opt=lam_opt.tolist(), E_free=e_free,
                                gain_mHa=gain,
                                spread=float(lam_opt.max() / lam_opt.min()))

print("=" * 78)
print("VERDICT")
print("=" * 78)
for nm in ("He", "Li+", "Li"):
    d = out["systems"][nm]
    shells = "1 shell " if nm != "Li" else "2 shells"
    print(f"  {nm:>4} ({shells}):  gain {d['gain_mHa']:7.2f} mHa   "
          f"lambda spread {d['spread']:.2f}x")
gHe, gLi = out["systems"]["He"]["gain_mHa"], out["systems"]["Li"]["gain_mHa"]
print(f"\n  Li / He gain ratio = {gLi / gHe:.1f}x   "
      f"(prediction: >> 1 if the two-length-scale story is right)")
out["verdict"] = dict(gain_He=gHe, gain_Li=gLi, ratio=gLi / gHe)

os.makedirs("debug/data", exist_ok=True)
with io.open("debug/data/r12ci_per_shell_lambda.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/r12ci_per_shell_lambda.json")
