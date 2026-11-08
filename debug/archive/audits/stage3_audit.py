"""AUDIT the beta*=sqrt2 result: is it a real minimum at the CM point, or a flat-landscape
grid coincidence? Fine sweep at validated grid, sqrt2 deliberately OFF the grid, parabola fit."""
import os, sys, time
import numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from two_center_grid_lm import TwoCenterLM
from two_center_ci_lm import h2_energy_lm

R = 1.4; ZS = 1.2; SQRT2 = 2 ** 0.5
eng = TwoCenterLM(R, nr=1400, nu=32, nphi=32, rmax=50.0, Lmax=14, real=True)

def E(beta, zs=ZS):
    zp = zs * beta
    orbs = [(zs,0,0,"A"),(zs,0,0,"B"),
            (zp,1,0,"A"),(zp,1,0,"B"),(zp,1,1,"A"),(zp,1,1,"B"),(zp,1,-1,"A"),(zp,1,-1,"B")]
    return h2_energy_lm(orbs, eng, R)

betas = [1.25, 1.35, 1.45, 1.55, 1.65]     # sqrt2=1.4142 deliberately NOT on the grid
pts = []
print("fine sweep (validated grid; sqrt2 off-grid):")
for b in betas:
    t0 = time.time(); e = E(b); pts.append((b, e))
    print(f"  beta={b:.3f}  E={e:.6f}  [{time.time()-t0:.0f}s]")

bs = np.array([p[0] for p in pts]); es = np.array([p[1] for p in pts])
c2, c1, c0 = np.polyfit(bs, es, 2)
bstar = -c1 / (2 * c2)
Estar = c0 + c1 * bstar + c2 * bstar ** 2
print(f"\nparabola fit: beta* = {bstar:.4f}   (sqrt2 = {SQRT2:.4f}, diff {bstar-SQRT2:+.4f})")

t0 = time.time(); e_s2 = E(SQRT2)
print(f"E(sqrt2) = {e_s2:.6f}   E(beta*_fit) ~ {Estar:.6f}   E(sqrt2)-E(beta*) = {e_s2-Estar:+.2e}")
# flatness: curvature -> how wide is the near-optimum region for a 1e-4 (grid-noise) energy change
half_width = (1e-4 / abs(c2)) ** 0.5
print(f"landscape curvature c2={c2:.4f}; beta-window for a 1e-4 energy change: +/-{half_width:.3f}")
print(f"\nREAD: if |beta*-sqrt2| << window and E(sqrt2)~E(beta*) within grid noise, sqrt2 is NOT")
print("distinguished (flat coincidence). A sharp min AT sqrt2 would need window << |grid points|.")
