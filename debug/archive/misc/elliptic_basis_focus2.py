"""True 2D exponent optimum + high-grid decisive points."""
import os, sys, time, json
import numpy as np
from scipy.optimize import minimize
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from elliptic_basis_pilot import h2_energy
from geovac.sturmian_integrals import GoscinskianIntegrals

R = 1.4; SQRT2 = 2 ** 0.5
gi = GoscinskianIntegrals(R, Lmax=14, nr=1600, nth=120, rmax=50.0)

# True 2-parameter optimum: free (z1,z2); the optimal RATIO is what CM would predict.
print("2D optimize E(z1,z2) [Nelder-Mead]:")
t0 = time.time()
res = minimize(lambda z: h2_energy([abs(z[0]), abs(z[1])], R, gi),
               x0=[1.0, 1.35], method="Nelder-Mead",
               options={"xatol": 3e-3, "fatol": 2e-6})
z1, z2 = sorted(abs(v) for v in res.x)
beta_opt = z2 / z1
print(f"  z* = ({z1:.4f}, {z2:.4f})   beta_opt = z2/z1 = {beta_opt:.4f}   E* = {res.fun:.6f}   [{time.time()-t0:.0f}s]")
print(f"  CM ratios: disc4/p2 sqrt2={SQRT2:.4f}  disc8/p2=1.0987  disc4/p1=2.0  disc8/p1=1.2071")
print(f"  beta_opt - sqrt2 = {beta_opt-SQRT2:+.4f}")

# high-grid decisive at fixed z_lo=1.0
print("\nhigh-grid decisive (z_lo=1.0):")
gih = GoscinskianIntegrals(R, Lmax=24, nr=3000, nth=200, rmax=60.0)
hd = {}
for label, b in [("beta*=1.35", 1.35), ("sqrt2(disc4)", SQRT2), ("1.207(disc8p1)", 1.207)]:
    e = h2_energy([1.0, b], R, gih); hd[label] = (round(b, 4), e)
    print(f"  {label:16s} E={e:.6f}")
print(f"  E(sqrt2)-E(1.35) = {hd['sqrt2(disc4)'][1]-hd['beta*=1.35'][1]:+.2e} Ha")

json.dump({"opt2d": {"z1": z1, "z2": z2, "beta_opt": beta_opt, "E": res.fun},
           "high_decisive": hd}, open(os.path.join(REPO, "debug/data/elliptic_basis_focus2.json"), "w"), indent=1)
print("saved.")
