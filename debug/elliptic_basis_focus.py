"""Pin beta*, test scale-robustness, quantify flatness. Is a CM ratio distinguished?"""
import os, sys, time, json
import numpy as np
from scipy.optimize import minimize_scalar
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from elliptic_basis_pilot import h2_energy
from geovac.sturmian_integrals import GoscinskianIntegrals

R = 1.4
gi = GoscinskianIntegrals(R, Lmax=14, nr=1600, nth=120, rmax=50.0)
SQRT2 = 2 ** 0.5

# 1) fine fixed-z_lo=1.0 grid to pin beta*
print("fine grid (z_lo=1.0):")
fine = []
for b in np.round(np.arange(1.25, 1.55, 0.025), 4):
    e = h2_energy([1.0, b], R, gi); fine.append((float(b), e))
    print(f"  beta={b:.3f}  E={e:.6f}")
bstar_fixed = min(fine, key=lambda t: t[1])
print(f"  -> beta*(z_lo=1.0) ~ {bstar_fixed[0]:.3f},  E={bstar_fixed[1]:.6f}")

# flatness: E span over beta in [1.30,1.45]
win = [e for b, e in fine if 1.30 <= b <= 1.45]
print(f"  flatness: E span over beta in [1.30,1.45] = {max(win)-min(win):.2e} Ha (grid noise ~2e-4)")

# 2) scale-robust: optimize z_lo at each of a few beta; is E*(beta) min still off sqrt2?
print("\nscale-relaxed E*(beta) = min over z_lo:")
star = {}
for b in [1.207, 1.30, 1.35, SQRT2, 1.45]:
    res = minimize_scalar(lambda z: h2_energy([z, z * b], R, gi),
                          bounds=(0.75, 1.6), method="bounded",
                          options={"xatol": 5e-3})
    star[round(b, 4)] = (float(res.x), float(res.fun))
    print(f"  beta={b:.4f}  z_lo*={res.x:.3f}  E*={res.fun:.6f}")
bstar_relaxed = min(star.items(), key=lambda kv: kv[1][1])
print(f"  -> best ratio (scale-relaxed) among tested: beta={bstar_relaxed[0]}  E*={bstar_relaxed[1][1]:.6f}")

# 3) high-grid decisive: E at beta* vs sqrt2 (disc-4) vs 1.207 (disc-8 p1)
print("\nhigh-grid decisive (z_lo=1.0):")
gih = GoscinskianIntegrals(R, Lmax=24, nr=3000, nth=200, rmax=60.0)
hd = {}
for label, b in [("beta*~1.35", 1.35), ("sqrt2 disc4", SQRT2), ("1.207 disc8", 1.207), ("1.099 disc8p2", 1.099)]:
    e = h2_energy([1.0, b], R, gih); hd[label] = (round(b, 4), e)
    print(f"  {label:16s} beta={b:.4f}  E={e:.6f}")
best_hd = min(hd.items(), key=lambda kv: kv[1][1])
print(f"  -> high-grid best among these: {best_hd[0]}  E={best_hd[1][1][1]:.6f}")
print(f"  E(sqrt2) - E(beta*=1.35) = {hd['sqrt2 disc4'][1]-hd['beta*~1.35'][1]:+.2e} Ha")

json.dump({"fine": fine, "scale_relaxed": {str(k): v for k, v in star.items()},
           "high_decisive": hd}, open(os.path.join(REPO, "debug/data/elliptic_basis_focus.json"), "w"), indent=1)
print("\nsaved debug/data/elliptic_basis_focus.json")
