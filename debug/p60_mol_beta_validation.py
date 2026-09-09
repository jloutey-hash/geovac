"""DIAGNOSTIC: does the molecular Sturmian beta machinery validate?

For ONE electron the weighting potential IS the true potential, V_0 = V, so
    V C = V_0 B C   =>   B = 1.
Equivalently: solving [-1/2 grad^2 - E + beta v_0] phi = 0 at the EXACT H2+
electronic eigenvalue must return beta = 1 for the ground orbital.  That is a
sharp, reference-backed check on geovac/molecular_sturmian.py, and everything
downstream of it is worthless if it fails.

H2+ electronic energies (exact, Burrau/Bates; nuclear repulsion NOT included):
    R = 2.0 bohr : -1.1026342144 Ha   (1s sigma_g)
    R = 2.5      : -1.0021818
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from geovac.molecular_sturmian import compute_molecular_sturmian_betas

REF = {2.0: -1.1026342144, 2.5: -1.0021818}

print("beta at the EXACT H2+ electronic energy (target: beta = 1 for the ground orbital)")
print(" R     E_exact        p0        n  m  nsph nrad   beta      |beta-1|")
for R, E in REF.items():
    p0 = np.sqrt(-2.0 * E)
    res = compute_molecular_sturmian_betas(1.0, 1.0, R, p0, nmax=3)
    if not res:
        print(" %.1f  %.7f  %.6f   <no roots found>" % (R, E, p0))
        continue
    for (n, m, nsph, nrad, beta) in res[:3]:
        print(" %.1f  %.7f  %.6f   %d  %d   %d    %d   %.6f   %.2e"
              % (R, E, p0, n, m, nsph, nrad, beta, abs(beta - 1.0)))
    print()
