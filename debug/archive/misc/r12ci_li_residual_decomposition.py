"""Decompose Li's 33.2 mHa post-per-shell-lambda residual into angular / radial-degree /
cusp axes -- the analogue of debug/sprint_explicit_correlation_scoping_memo.md Sec. 4 (H2)
for the atomic engine, owed by CLAUDE.md Sec. 3 / the gamma-transfer memo Sec. 6+9.

Baseline (reproduced by this driver before anything else, per the current-state-check
rule): debug/r12ci_per_shell_lambda.py's Li result -- ns=4, s-only, free per-shell lambda,
NO geminal:  E = -7.444876 Ha,  33.2 mHa above exact.

This engine (debug/r12ci_ne_engine.py) has no l>0 orbital shapes -- adding them would be a
new method, out of scope (PI directive: measurement only).  So "angular" is bounded, not
built:  the radial-degree axis is pushed to the point of diminishing returns (more s
functions, same l=0 subspace), the cusp axis is the existing explicit-r12 geminal
machinery (which -- caveat stated up front, not discovered after the fact -- has
Legendre-moment content up to Lmax=20 via KernelBank.mom(), so it implicitly carries SOME
angular correlation through r_12 the way Hylleraas' 1929 He calculation did; a gamma-grid
scan cannot separate "cusp" from "the L>0 content of exp(-gamma r12)").  What survives
BOTH axes is the operational estimate of "angular correlation neither axis can reach.

One axis moves at a time from the ns=4 baseline, exactly as Sec. 4 did for H2.
"""
from __future__ import annotations

import io
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)

import r12ci_ne_engine as E  # noqa: E402

# Exact non-relativistic (infinite-nuclear-mass, clamped) Li ground state.  Cross-checked
# 2026-09-19 against an independent literature search returning -7.478060323452 Ha
# (Yan & Drake / King lineage) -- agrees with the corpus's registered EXACT["Li"] to every
# digit quoted.
EXACT_LI = -7.4780603236
Z, N, MS2 = 3.0, 3, 1
NG, NX, LMAX = 240, 48, 20
BASE_LAMS = [3.153, 2.208, 1.103, 1.505]   # per_shell_lambda.py's ns=4 optimum

out = {"exact": EXACT_LI, "grid": dict(Ng=NG, nx=NX, Lmax=LMAX), "baseline_lams": BASE_LAMS}


def energy_plain(lams, Ng=NG):
    ns = len(lams)
    S, H, _ = E.build(ns, N, Z, float(np.mean(lams)), 0.5, MS2, Ng=Ng,
                       with_geminal=False, nx=NX, lowdin=False, lams=list(lams))
    e, _ = E.solve(S, H)
    return e


def energy_gem(lams, gammas, Ng=NG, Lmax=LMAX):
    ns = len(lams)
    S, H, dets = E.build(ns, N, Z, float(np.mean(lams)), gammas, MS2, Ng=Ng,
                          with_geminal=True, nx=NX, lowdin=False, lams=list(lams), Lmax=Lmax)
    e, nk = E.solve(S, H, nd=len(dets))
    return e, nk


def mHa(e):
    return 1000.0 * (e - EXACT_LI)


print("=" * 78)
print("STEP 0 -- reproduce the per-shell-lambda Li baseline (current-state check)")
print("=" * 78)
t0 = time.time()
E0 = energy_plain(BASE_LAMS)
print(f"  E(ns=4, free per-shell lambda, no geminal) = {E0:.6f}  "
      f"({mHa(E0):.2f} mHa above exact)   [{time.time()-t0:.1f}s]")
assert abs(mHa(E0) - 33.2) < 0.5, "baseline did not reproduce -- STOP condition"
out["baseline"] = dict(E=E0, err_mHa=mHa(E0))
print()

# ---------------------------------------------------------------------------
# AXIS 1 -- radial degree: raise ns (more l=0 radial functions) at free per-shell lambda.
# ---------------------------------------------------------------------------
print("=" * 78)
print("AXIS 1 -- radial degree (ns), l_max fixed at 0, no geminal")
print("=" * 78)


def optimize_lams(start, maxiter):
    res = minimize(lambda u: energy_plain(np.exp(u)), np.log(start), method="Nelder-Mead",
                    options=dict(maxiter=maxiter, xatol=2e-4, fatol=2e-8))
    return np.exp(res.x), float(res.fun), res.nit


radial_rows = [dict(ns=4, lams=BASE_LAMS, E=E0, err_mHa=mHa(E0), gain_mHa=None)]
prev_lams, prev_e = list(BASE_LAMS), E0

for ns, polish_iter in ((5, 160), (6, 170), (7, 150)):
    t0 = time.time()
    # (a) freeze existing shells, 1-D optimize the NEW exponent (cheap, good warm start)
    def obj1d(x, _prev=tuple(prev_lams)):
        return energy_plain(list(_prev) + [float(x)])
    r1 = minimize_scalar(obj1d, bounds=(0.25, 7.0), method="bounded",
                          options=dict(xatol=2e-3))
    warm = prev_lams + [float(r1.x)]
    e_warm = energy_plain(warm)
    # (b) short joint polish (all ns exponents together) from that warm start
    lam_opt, e_opt, nit = optimize_lams(warm, maxiter=polish_iter)
    if e_opt > e_warm:   # Nelder-Mead can rarely fail to improve within the iter budget
        lam_opt, e_opt = np.array(warm), e_warm
    gain = 1000.0 * (prev_e - e_opt)
    dt = time.time() - t0
    print(f"  ns={ns}: warm(1D)={mHa(e_warm):7.2f} mHa -> polished={mHa(e_opt):7.2f} mHa   "
          f"gain={gain:6.2f} mHa   lams={np.array2string(lam_opt, precision=3)}   "
          f"[{dt:.0f}s, {nit} it]")
    radial_rows.append(dict(ns=ns, lams=lam_opt.tolist(), E=e_opt, err_mHa=mHa(e_opt),
                             gain_mHa=gain, wall_s=dt))
    prev_lams, prev_e = lam_opt.tolist(), e_opt

out["radial_ladder"] = radial_rows
radial_total_gain = 1000.0 * (E0 - prev_e)
print(f"\n  TOTAL radial-degree gain (ns 4->7) = {radial_total_gain:.2f} mHa   "
      f"residual now {mHa(prev_e):.2f} mHa")
ns7_lams, ns7_e = list(prev_lams), prev_e
print()

# ---------------------------------------------------------------------------
# AXIS 2 -- cusp: turn on the explicit r12 geminal, AT THE ns=4 BASELINE (apples-to-apples
# with the reported 33.2 mHa), lambda held fixed so only the geminal axis moves.
# ---------------------------------------------------------------------------
print("=" * 78)
print("AXIS 2 -- cusp (explicit r12 geminal), ns=4 baseline lambda held fixed")
print("=" * 78)

GAMMA_GRID = [0.10, 0.20, 0.30, 0.45, 0.65, 0.90, 1.10, 1.40, 1.70, 2.00, 2.50, 3.00]
rows1 = []
t0 = time.time()
for g in GAMMA_GRID:
    e, nk = energy_gem(BASE_LAMS, g)
    rows1.append(dict(gamma=g, E=e, err_mHa=mHa(e), kept=nk))
best1 = min(rows1, key=lambda r: r["E"])
print("  1-geminal scan:")
for r in rows1:
    flag = "  <-- best" if r is best1 else ""
    print(f"    gamma={r['gamma']:.2f}  E={r['E']:.6f}  err={r['err_mHa']:6.2f} mHa{flag}")
# refine around the grid optimum
g0 = best1["gamma"]
r2 = minimize_scalar(lambda g: energy_gem(BASE_LAMS, float(g))[0],
                      bounds=(max(0.02, g0 * 0.6), g0 * 1.6), method="bounded",
                      options=dict(xatol=1e-3))
e_ref, nk_ref = energy_gem(BASE_LAMS, float(r2.x))
gain_1gem = 1000.0 * (E0 - e_ref)
print(f"  refined optimum: gamma={float(r2.x):.4f}  E={e_ref:.6f}  err={mHa(e_ref):.2f} mHa  "
      f"GAIN={gain_1gem:.2f} mHa   [{time.time()-t0:.0f}s]")
out["cusp_1gem"] = dict(grid=rows1, gamma_opt=float(r2.x), E=e_ref, err_mHa=mHa(e_ref),
                         gain_mHa=gain_1gem)
print()

# 2-geminal (coarse grid): does a second correlation length buy more, as in Sec. 8 of the
# gamma-transfer memo (there: +0.394 mHa on a DIFFERENT lambda structure)?
print("  2-geminal coarse grid (does a second correlation length help here too?):")
t0 = time.time()
G2 = [0.15, 0.45, 0.90, 1.70, 2.50]
rows2 = []
best2 = None
for i, ga in enumerate(G2):
    for gb in G2[i:]:
        e, nk = energy_gem(BASE_LAMS, [ga, gb])
        rows2.append(dict(ga=ga, gb=gb, E=e, err_mHa=mHa(e)))
        if best2 is None or e < best2["E"]:
            best2 = rows2[-1]
gain_2gem_extra = 1000.0 * (e_ref - best2["E"])
print(f"    best pair (ga,gb)=({best2['ga']:.2f},{best2['gb']:.2f})  E={best2['E']:.6f}  "
      f"err={best2['err_mHa']:.2f} mHa   EXTRA over 1 geminal={gain_2gem_extra:.3f} mHa   "
      f"[{time.time()-t0:.0f}s]")
out["cusp_2gem"] = dict(grid=rows2, best=best2, extra_gain_mHa=gain_2gem_extra)
print()

# ---------------------------------------------------------------------------
# AXIS 1+2 combined -- geminal ON TOP OF the radial-degree-converged (ns=7) basis, to see
# whether the two axes are additive or whether they compete for the same correlation.
# ---------------------------------------------------------------------------
print("=" * 78)
print("COMBINED -- cusp geminal at the ns=7 radial-converged basis (interaction check)")
print("=" * 78)
t0 = time.time()
rows3 = []
for g in [0.20, 0.45, 0.65, 0.90, 1.10, 1.40, 1.70]:
    e, nk = energy_gem(ns7_lams, g)
    rows3.append(dict(gamma=g, E=e, err_mHa=mHa(e)))
best3 = min(rows3, key=lambda r: r["E"])
gain_combined = 1000.0 * (ns7_e - best3["E"])
print(f"  best gamma={best3['gamma']:.2f}  E={best3['E']:.6f}  err={best3['err_mHa']:.2f} mHa  "
      f"geminal gain AT ns=7 = {gain_combined:.2f} mHa (vs {gain_1gem:.2f} mHa at ns=4)   "
      f"[{time.time()-t0:.0f}s]")
out["combined_ns7_gem"] = dict(grid=rows3, best=best3, gain_mHa=gain_combined)
print()

# ---------------------------------------------------------------------------
# VERDICT -- decomposition table
# ---------------------------------------------------------------------------
print("=" * 78)
print("VERDICT -- decomposition of the 33.2 mHa residual")
print("=" * 78)
final_err_both_axes = best3["err_mHa"]
angular_estimate = final_err_both_axes
print(f"  baseline (ns=4, per-shell lambda, no geminal):        {mHa(E0):7.2f} mHa")
print(f"  radial-degree gain (ns 4->7):                         {radial_total_gain:7.2f} mHa")
print(f"  cusp/geminal gain (1 geminal, AT ns=4):                {gain_1gem:7.2f} mHa")
print(f"  cusp/geminal gain (1 geminal, AT ns=7 -- combined):    {gain_combined:7.2f} mHa")
print(f"  extra from a 2nd correlation length (AT ns=4):         {gain_2gem_extra:7.3f} mHa")
print(f"  --------------------------------------------------------------")
print(f"  residual surviving BOTH axes (ns=7 + best geminal):    {angular_estimate:7.2f} mHa")
print(f"  = {100*angular_estimate/33.2:.1f}% of the original 33.2 mHa")
out["verdict"] = dict(
    baseline_err_mHa=mHa(E0),
    radial_gain_mHa=radial_total_gain,
    cusp_gain_at_ns4_mHa=gain_1gem,
    cusp_gain_at_ns7_mHa=gain_combined,
    second_gamma_extra_mHa=gain_2gem_extra,
    residual_surviving_both_mHa=angular_estimate,
    pct_of_original_surviving=100 * angular_estimate / 33.2,
)

os.makedirs("debug/data", exist_ok=True)
with io.open("debug/data/r12ci_li_residual_decomposition.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/r12ci_li_residual_decomposition.json")
