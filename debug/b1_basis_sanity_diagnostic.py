"""B.1 basis sanity diagnostic.

Two quick checks before committing to the explicit-core NaH HF sprint:

(4a) Hydrogenic eigenvalue spectrum at Z=11 for n=1,2,3 vs the [Ne] core's
     true ground-state energy after electron-electron repulsion.

(4b) Na 3s orbital radial structure at bare Z=11 vs Z_eff=2.2 (the screened
     Clementi-Raimondi value the existing frozen-core path uses).  If they're
     wildly different we have a real concern about HF on the hydrogenic basis.
"""

from __future__ import annotations

import numpy as np
import sys

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass


# (4a) Hydrogenic eigenvalues
print("=" * 72)
print("(4a) Hydrogenic eigenvalue spectrum at bare Z=11 (Na nucleus)")
print("=" * 72)
Z = 11
print(f"{'n':>3s}  {'E_hydrogenic':>15s}  {'occupancy':>10s}  {'sum':>15s}")
total_bare = 0.0
for n, occ in [(1, 2), (2, 8), (3, 1)]:  # 1s², 2s²+2p⁶, 3s¹
    E_per_e = -Z**2 / (2 * n**2)
    contrib = occ * E_per_e
    total_bare += contrib
    print(f"{n:3d}  {E_per_e:+15.4f}  {occ:>10d}  {contrib:+15.4f}")
print(f"\nTotal bare-Coulomb 11-electron energy (no e-e): {total_bare:+.4f} Ha")
print("Reference -- Na atom NIST experimental ground state: -162.255 Ha")
print(f"Difference (must equal -<V_ee>): {total_bare - (-162.255):+.4f} Ha")
print(f"  This is the e-e repulsion HF must reproduce.")

# (4b) Na 3s orbital extent
print()
print("=" * 72)
print("(4b) Na 3s orbital radial structure")
print("=" * 72)

# Hydrogenic 3s wavefunction radial part:
#   R_30(r; Z) = (1/(9*sqrt(3))) * Z^{3/2} * (6 - 6*Z*r + (Z*r)**2 / ... )
# We don't need the exact analytic form -- just the radial probability density
# peak location and the rough extent (~mean radius).
#
# For hydrogenic 3s, the mean radius is <r> = (3*a_0 / (2*Z)) * (3*n^2 - l*(l+1))
#                                          = (3 / (2*Z)) * (3*9 - 0)
#                                          = 40.5 / Z  bohr
# At Z=11 bare:    <r> = 3.68 bohr  (but this is for a 3s with full 11-charge
#                                    attraction; the bare 11-charge in fact
#                                    pulls the orbital inward.)

# Actually the correct formula for hydrogenic <r> for n=3, l=0:
#   <r>_{n,l} = (a_0 / Z) * [3*n^2 - l*(l+1)] / 2
#             = (1/Z) * (27 - 0) / 2 = 13.5 / Z
# At Z=11: <r> = 1.23 bohr
# At Z_eff=2.2: <r> = 6.14 bohr

import math
def mean_r_hydrogenic_3s(Z_eff):
    return 13.5 / Z_eff  # bohr

print(f"{'Z_eff':>8s}  {'mean_r_3s':>12s}  comment")
for Z_eff in [11.0, 5.5, 2.2, 1.0]:
    r = mean_r_hydrogenic_3s(Z_eff)
    if Z_eff == 11.0:
        c = "bare Na nucleus (no screening); orbital deeply inside core"
    elif Z_eff == 2.2:
        c = "Clementi-Raimondi screened (used by frozen-core spec)"
    elif Z_eff == 5.5:
        c = "partial screening (5 of 10 core e-)"
    else:
        c = "free 3s reference"
    print(f"  {Z_eff:6.2f}  {r:8.2f} bohr  {c}")

# Compare to [Ne] core extent
# [Ne] 2s and 2p mean r at Z=11:
print("\n[Ne] core orbital mean radii at Z=11 (bare):")
print(f"  <r>_{{1s, Z=11}} = {1.5 / 11:.3f} bohr (= 3/(2*Z) for n=1, l=0)")
print(f"  <r>_{{2s, Z=11}} = {(3*4 - 0)/2 / 11:.3f} bohr")
print(f"  <r>_{{2p, Z=11}} = {(3*4 - 2)/2 / 11:.3f} bohr")

# Where does the 3s peak relative to the [Ne] core?
print()
print("Spatial overlap question:")
print("  Hydrogenic 3s at bare Z=11 peaks around r ~ 1-2 bohr.")
print("  [Ne] core 2p extends to r ~ 0.5 bohr.")
print("  So the bare-Z hydrogenic 3s sits *just outside* the [Ne] core, but")
print("  with strong inner radial nodes inside the core region.")
print()
print("Implication for HF: the hydrogenic Na 3s at bare Z=11 will have")
print("strong inner-node spike inside the core volume; the [Ne] core electrons")
print("screen this to push the 3s outward to ~5-6 bohr in a real Na atom.")
print("HF self-consistency will move the 3s outward via the eri tensor.")

# Final verdict
print()
print("=" * 72)
print("Verdict on basis adequacy")
print("=" * 72)
print()
print("(4a) The hydrogenic Z=11 basis gives total bare-Coulomb energy")
print("     ~-1200 Ha for the 1s2 2s2 2p6 3s1 occupancy. The true Na atom")
print("     ground state is -162 Ha. HF on this basis needs to recover")
print("     ~+1000 Ha of electron-electron repulsion. This is enormous,")
print("     consistent with the standard observation that bare hydrogenic")
print("     basis at high Z is energetically far from physical.")
print()
print("     The bond block has the same issue: bare Z=11 valence orbital")
print("     is energetically far from the screened orbital that's physical.")
print()
print("     HOWEVER: bond energy differences (PES topology) depend on the")
print("     SCF self-consistency closing the gap correctly, not on the")
print("     absolute energy. HF can reach the right relative energies even")
print("     when the initial bare-Z guess is far off, provided the basis")
print("     has enough variational freedom.")
print()
print("(4b) Hydrogenic 3s at Z=11 (mean r ~ 1.2 bohr) is geometrically too")
print("     compact. The physical Na 3s is at ~6 bohr.  HF can pull this")
print("     outward via mixing with higher-n orbitals - but max_n=2 has NO")
print("     higher-n freedom, only 1s/2s/2p/(3s/3p/3d).")
print()
print("     WAIT - we need max_n=3 to have a 3s orbital at all.")
print("     max_n=2 only gives n=1 (1s) and n=2 (2s, 2p) - no n=3.")
print()
print("     ---> max_n MUST be at least 3 for explicit-core NaH at HF.")
print("          The 3s is the valence orbital; it cannot be omitted.")
print()
print("     At max_n=3 the orbital count is:")
print("       n=1: 1 (1s)")
print("       n=2: 4 (2s + 2p)")
print("       n=3: 9 (3s + 3p + 3d)")
print("       Total: 14 spatial = 28 spin-orbitals for 12 electrons.")
print()
print("     This is workable. Qubit count would be ~28 for HF (no FCI).")
print()
print("RECOMMENDATION: Proceed to B.1 implementation but with max_n=3 NaH spec")
print("                from the start. Plan agent's '4 days' estimate stays")
print("                valid; max_n=3 increases cross-block ERI cost by ~3^4")
print("                vs max_n=2 = 16x more eri terms, but with HF (not FCI)")
print("                the bottleneck is the SCF iteration not the integral build.")
