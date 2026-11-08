"""Stage 3: with angular correlation IN the basis, is a CM exponent ratio variationally
special? Sweep beta=zeta_p/zeta_s in s+full-p H2; compare beta* to CM singular moduli.
The culminating 'does the bond's elliptic geometry help real (angular) accuracy' test."""
import os, sys, time, json
import numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from two_center_grid_lm import TwoCenterLM
from two_center_ci_lm import h2_energy_lm

R = 1.4; ZS = 1.2; EXACT = -1.17447
eng = TwoCenterLM(R, nr=1100, nu=24, nphi=24, rmax=50.0, Lmax=12, real=True)

def E(beta):
    zp = ZS * beta
    orbs = [(ZS,0,0,"A"),(ZS,0,0,"B"),
            (zp,1,0,"A"),(zp,1,0,"B"),(zp,1,1,"A"),(zp,1,1,"B"),(zp,1,-1,"A"),(zp,1,-1,"B")]
    return h2_energy_lm(orbs, eng, R)

# CM ratios beta = z_max/z_min via m = 1 - (z_min/z_max)^2 = m_CM   (c ~ 1/zeta^2)
m8 = 3 - 2 * 2 ** 0.5
CM = {"disc4 (m=1/2)": 2 ** 0.5, "disc8 (m=3-2rt2)": (1 - m8) ** -0.5,
      "disc4/p1": 2.0, "disc8/p1": 1 / (1 - m8)}

betas = [0.85, 1.0, 1.099, 1.2, 1.30, 1.414, 1.55, 1.75, 2.0, 2.3]
rows = []
print(f"s+full-p H2 (R={R}, zeta_s={ZS}); exact {EXACT}. Sweep beta=zeta_p/zeta_s:")
for b in betas:
    t0 = time.time(); e = E(b); rows.append((b, e))
    print(f"  beta={b:5.3f}  E_tot={e:.6f}  gap {e-EXACT:+.4f}  [{time.time()-t0:.0f}s]")

bstar, estar = min(rows, key=lambda t: t[1])
mstar = 1 - 1 / bstar ** 2 if bstar > 1 else 1 - bstar ** 2
print(f"\nvariational optimum: beta* = {bstar:.3f}  E={estar:.6f}  -> modulus m(beta*) = {mstar:.3f}")
print(f"CM singular moduli: m=0.5 (disc-4), m={m8:.3f} (disc-8)")
print("CM ratios:", {k: round(v, 3) for k, v in CM.items()})
near = [(k, v) for k, v in CM.items() if abs(v - bstar) < 0.06]
print(f"\nVERDICT: beta*={bstar:.3f} " + (f"NEAR CM {near}" if near else "is GENERIC (no CM within 0.06)")
      + f";  m(beta*)={mstar:.3f} vs CM {{0.5, {m8:.3f}}}")
json.dump({"rows": rows, "bstar": bstar, "mstar": mstar, "CM": CM},
          open(os.path.join(REPO, "debug/data/stage3_elliptic_angular.json"), "w"), indent=1)
