"""Stage 1: does adding angular (p) functions drop H2 energy from the s-limit toward exact?
Turns the earlier INFERENCE ('the gap is angular') into a MEASUREMENT."""
import os, sys, time
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from two_center_grid_lm import TwoCenterLM
from two_center_ci_lm import h2_energy_lm

R = 1.4
EXACT = -1.17447
eng = TwoCenterLM(R, nr=1400, nu=32, nphi=32, rmax=50.0, Lmax=14, real=True)

print("real-harmonic isotropy check (same-center p, zeta=1.0): S=1, T & <1/rA> equal across m")
for m in (0, 1, -1):
    o = (1.0, 1, m, "A")
    print(f"  p(m={m:+d}): S={eng.overlap(o,o):.6f}  T={eng.kinetic(o,o):.6f}  <1/rA>={eng.coulomb_center(o,o,'A'):.6f}")

zs, zp = 1.20, 1.35
s   = [(zs,0,0,"A"),(zs,0,0,"B")]
spz = s   + [(zp,1,0,"A"),(zp,1,0,"B")]
sp  = spz + [(zp,1,1,"A"),(zp,1,1,"B"),(zp,1,-1,"A"),(zp,1,-1,"B")]
print(f"\nH2 R={R}, single-zeta (zeta_s={zs}, zeta_p={zp}); exact = {EXACT}")
prev = None
for label, orbs in [("s only", s), ("s + p_z (sigma)", spz), ("s + full p (sigma+pi)", sp)]:
    t0 = time.time()
    E = h2_energy_lm(orbs, eng, R)
    drop = f"  (drop {E-prev:+.4f})" if prev is not None else ""
    print(f"  {label:24s} E_tot = {E:.6f}   gap-to-exact {E-EXACT:+.4f}{drop}   [{time.time()-t0:.0f}s]")
    prev = E
