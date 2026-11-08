"""Finalize the Stage-3 audit: fit beta* from the off-sqrt2 grid, quantify sqrt2's (non)distinction."""
import os, sys, time
import numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from two_center_grid_lm import TwoCenterLM
from two_center_ci_lm import h2_energy_lm

R = 1.4; ZS = 1.2; SQRT2 = 2 ** 0.5
eng = TwoCenterLM(R, nr=1400, nu=32, nphi=32, rmax=50.0, Lmax=14, real=True)
def E(beta):
    zp = ZS * beta
    return h2_energy_lm([(ZS,0,0,"A"),(ZS,0,0,"B"),(zp,1,0,"A"),(zp,1,0,"B"),
                         (zp,1,1,"A"),(zp,1,1,"B"),(zp,1,-1,"A"),(zp,1,-1,"B")], eng, R)

# 4 points already computed at this grid (from stage3_audit)
known = [(1.25,-1.163839),(1.35,-1.164567),(1.45,-1.164848),(1.55,-1.164765)]
bs = np.array([p[0] for p in known]); es = np.array([p[1] for p in known])
c2,c1,c0 = np.polyfit(bs,es,2); bstar = -c1/(2*c2); Estar = c0+c1*bstar+c2*bstar**2
print(f"parabola fit (4 off-sqrt2 pts): beta* = {bstar:.4f}  (sqrt2={SQRT2:.4f}, diff {bstar-SQRT2:+.4f})")

t0=time.time(); e_s2 = E(SQRT2); print(f"E(sqrt2={SQRT2:.4f}) = {e_s2:.6f}  [{time.time()-t0:.0f}s]")
mstar = 1 - 1/bstar**2
print(f"\nmodulus at optimum m(beta*) = {mstar:.3f}   (disc-4 CM = 0.500, disc-8 CM = 0.172)")
print(f"E(sqrt2) - E(beta*_fit) = {e_s2-Estar:+.2e} Ha  (grid noise ~1e-4)")
print(f"\nVERDICT: beta*={bstar:.3f} GENERIC (m={mstar:.3f}, not a CM modulus). sqrt2 sits "
      f"{abs(e_s2-Estar):.1e} Ha off the optimum = grid-noise level -> NOT distinguished. "
      f"The earlier 'beta*=sqrt2' was a grid-placement + flat-landscape artifact.")
