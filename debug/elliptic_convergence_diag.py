"""Diagnostic: does the two-center H2 radial convergence tail carry MODULAR structure
(resummable) or is it smooth/geometric (nothing to resum)?  (2026-08-23)

The elliptic period is a RADIAL-momentum object (angular part closed to Bessel j0), so the
s-only two-center engine is the right place to look for it.  Add K even-tempered 1s
exponents/center, K=1..4; read the increments dE(K)=E(K)-E(K-1).  Geometric (constant
ratio dE(K+1)/dE(K)) => smooth, modular idea has nothing to grab.  Non-geometric / stepped
=> possible modular structure.
"""
import os, sys, time, json
import numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from elliptic_basis_pilot import h2_energy
from geovac.sturmian_integrals import GoscinskianIntegrals

R = 1.4
gi = GoscinskianIntegrals(R, Lmax=14, nr=1600, nth=120, rmax=50.0)

def even_tempered(K, center=1.15, ratio=1.6):
    idx = np.arange(K) - (K - 1) / 2.0
    return list(center * ratio ** idx)

print(f"H2 R={R}, even-tempered 1s ladders (center=1.15, ratio=1.6):")
rows = []
for K in range(1, 5):
    exps = even_tempered(K)
    t0 = time.time()
    E = h2_energy(exps, R, gi)
    rows.append((K, E, exps))
    print(f"  K={K}  E_tot={E:.6f}   exps={[round(e,3) for e in exps]}   [{time.time()-t0:.0f}s]")

print("\nconvergence increments and ratios:")
prev = None; prevd = None
for K, E, _ in rows:
    if prev is not None:
        d = E - prev
        rat = (d / prevd) if prevd not in (None, 0) else float('nan')
        print(f"  dE({K}) = E({K})-E({K-1}) = {d:+.2e} Ha    ratio dE({K})/dE({K-1}) = {rat:.3f}")
        prevd = d
    prev = E

# geometric extrapolation to K->inf (if ratios ~ constant r): E_inf ~ E(K) + dE(K)*r/(1-r)
Es = [E for _, E, _ in rows]
if len(Es) >= 3:
    dlast = Es[-1] - Es[-2]; dprev = Es[-2] - Es[-3]
    r = dlast / dprev if dprev else float('nan')
    Einf = Es[-1] + dlast * r / (1 - r) if abs(r) < 1 else float('nan')
    print(f"\ngeometric-extrapolated s-limit E_inf ~ {Einf:.6f}  (last ratio r={r:.3f})")
    print(f"remaining tail |E(4)-E_inf| ~ {abs(Es[-1]-Einf):.2e} Ha")
print("\nread: ratios ~ constant => SMOOTH/geometric (modular idea inert). "
      "ratios jumping / non-monotone => structure worth resumming.")
json.dump({"R": R, "rows": [(K, E) for K, E, _ in rows]},
          open(os.path.join(REPO, "debug/data/elliptic_convergence_diag.json"), "w"), indent=1)
