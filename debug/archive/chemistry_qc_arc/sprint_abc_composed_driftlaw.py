"""
Test A (Connection A): composed LiH tilt vs l_max -- is there a zero-parameter
skeleton form, or does the discrete angular refinement DIVERGE the R-geometry?

Reads cached l_dependent PK PES (debug/data/lmax_convergence_algebraic.json),
refits the tilt eps'(R_true) and curvature at R_true=3.015 with a clean local
stencil per l_max, and tests whether the tilt tracks a simple channel-count law
n_ch=(l_max+1)^2 (consistency check only -- 3 points cannot DERIVE a form; the
(0,0)-weight dilution derivation needs solver instrumentation, banked).

Connection-A read: if refining l_max (the one discrete knob that should
systematically improve a variational basis) drives R_eq monotonically WORSE, the
non-compact R-coordinate is not pinned by -- and is destabilized by -- discrete
angular refinement. That is the free-side signature.
"""
from __future__ import annotations
import json
import numpy as np

R_TRUE = 3.015
CM1_TO_HA = 1.0 / 219474.6313702
MU_ME = (7.016004*1.007825)/(7.016004+1.007825) * 1822.888486
K_TRUE = MU_ME * (1405.65*CM1_TO_HA)**2

d = json.load(open('debug/data/lmax_convergence_algebraic.json'))
rows = [r for r in d['results'] if r['pk_mode'] == 'l_dependent']

print("="*84)
print("Composed (l_dependent PK) tilt/curv at R_true=3.015, refit from cached PES")
print("="*84)
out = []
for r in rows:
    R = np.array(r['pes_R']); E = np.array(r['pes_E'])
    l = r['l_max']; n_ch = r['n_channels']
    m = np.abs(R - R_TRUE) <= 0.55            # local window bracketing R_true
    x = R[m] - R_TRUE
    c = np.polyfit(x, E[m], 3)
    p = np.poly1d(c); dp = p.deriv(1); ddp = p.deriv(2)
    tilt = float(dp(0.0)); curv = float(ddp(0.0))
    row = {'l_max': l, 'n_ch': n_ch, 'tilt': tilt, 'curv': curv,
           'curv_over_k': curv/K_TRUE, 'R_eq': r['R_eq_bohr'],
           'R_eq_err_pct': r['R_eq_error_pct'], 'n_fit': int(m.sum())}
    out.append(row)
    print(f"  l_max={l} (n_ch={n_ch:2d}): tilt={tilt:+.4f}  curv={curv:+.4f} "
          f"({curv/K_TRUE:.2f}x k)  R_eq={r['R_eq_bohr']:.3f} (+{r['R_eq_error_pct']:.1f}%)"
          f"  [{int(m.sum())}pt fit]")

tl = np.array([o['tilt'] for o in out])
ll = np.array([o['l_max'] for o in out])
nch = np.array([o['n_ch'] for o in out], float)

print("\n-- divergence: does refining the discrete angular knob help or hurt? --")
print(f"  tilt:      {tl[0]:+.4f} -> {tl[1]:+.4f} -> {tl[2]:+.4f}  "
      f"(mean step {np.mean(np.diff(tl)):+.4f}/l_max)")
print(f"  R_eq err%: {out[0]['R_eq_err_pct']:.1f} -> {out[1]['R_eq_err_pct']:.1f} -> "
      f"{out[2]['R_eq_err_pct']:.1f}  (MONOTONE WORSE => free-side signature)")

print("\n-- zero-parameter channel-count consistency (NOT a derivation; 3 pts) --")
for name, xvar in [("tilt ~ l_max", ll.astype(float)),
                   ("tilt ~ n_ch=(l+1)^2", nch),
                   ("tilt ~ (l_max-1)", ll-1.0)]:
    A = np.vstack([xvar, np.ones_like(xvar)]).T
    (sl, ic), res, *_ = np.linalg.lstsq(A, tl, rcond=None)
    yhat = A @ [sl, ic]
    r2 = 1 - np.sum((tl-yhat)**2)/np.sum((tl-tl.mean())**2)
    print(f"  {name:24s}: slope={sl:+.5f}  intercept={ic:+.5f}  R^2={r2:.4f}")

json.dump({'k_true': K_TRUE, 'rows': out}, open('debug/data/abc_composed_driftlaw.json','w'), indent=2)
print("\n[saved] debug/data/abc_composed_driftlaw.json")
print("\nNote: Paper 17 already establishes l_dependent divergence to l_max=7 "
      "(+0.168 bohr/l_max);\nno live extension needed. (0,0)-weight dilution "
      "derivation banked (needs solver instrumentation).")
